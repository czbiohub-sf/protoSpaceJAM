#!/usr/bin/env python
"""
Build a variant-substituted genome (and its variant sidecar) from a VCF.

The output is consumed by protoSpaceJAM in two ways:

  1. materialized  - genome_files/fa_pickle/<out_genome_ver>/<chr>.pk
                     drop-in replacement for the reference pickles; get_seq needs no changes
  2. runtime       - genome_files/variant_sets/<out_genome_ver>/applied/<chr>.pk
                     the differences only; patched into each slice as it is fetched

(1) is the oracle that (2) is validated against, see 5_validate_variant_genome.py

Only length-preserving substitutions are applied.  Indels are recorded in the sidecar but never
written into the sequence: every other precomputed asset (parsed_gff3/loc2posType, ENST_codonPhases,
precomputed_gRNAs coordinates) is keyed to reference coordinates, so a single inserted or deleted
base would put all of them out of register.

No third-party dependencies beyond Biopython (already required): the VCF is read with the stdlib
gzip module, which handles BGZF fine for a sequential pass.

Example
-------
    python 4_make_variant_genome.py \
        --vcf   variant/czb_kolf.joint.GRCh38.small_variants.phased.norm.slivar.vcf.gz \
        --sample czML910_Parental_Ngn2-KOLF \
        --genome_ver GRCh38 \
        --out_genome_ver GRCh38_KOLF2.1J_hap1 \
        --validate_only
"""

import argparse
import array
import datetime
import gzip
import hashlib
import io
import json
import os
import pickle
import random
import sys
from collections import Counter, defaultdict

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


# ----------------------------------------------------------------------------------------------
# reasons a variant is recorded in the sidecar but not written into the sequence
# An indel is skipped either way, but *which haplotype it sits on* changes everything downstream:
#   - on the selected haplotype, the emitted sequence is knowingly wrong in that window, so any
#     homology arm or guide overlapping it must be flagged (Phase 1) or remapped (Phase 2)
#   - on the other haplotype, the selected haplotype really is reference-length there and the
#     emitted sequence is correct; no warning is warranted
SKIP_INDEL_SELECTED_HAP = "indel_selected_hap"
SKIP_INDEL_OTHER_HAP = "indel_other_hap"
SKIP_INDEL_UNPHASED = "indel_unphased"
SKIP_UNPHASED_HET = "unphased_het"
SKIP_OTHER_HAPLOTYPE = "other_haplotype"
SKIP_NO_CALL = "no_call"
SKIP_REF_HOM = "ref_hom"
SKIP_REF_MISMATCH = "ref_mismatch"
SKIP_SYMBOLIC = "symbolic_allele"
# Two substitutions whose REF footprints overlap cannot both be written; the second is demoted.
# Splitting a multiallelic (`bcftools norm -m -any`) turns a `1|2` into a `1|0` and an `0|1` at the
# same POS, so exactly one of the pair is on the selected haplotype and this never fires -- but if
# the input is ever built differently, the alternative is silently writing whichever came first.
SKIP_OVERLAPPING_SUB = "overlapping_substitution"
# A substitution we would have written, suppressed because it sits in a long homopolymer run.
# Unlike every other skip reason this one does NOT make the window unreliable: in runs >=7 the
# SNV Ti/Tv is 0.43-0.45, below the 0.5 random-error floor, so the reference base is the better
# bet and calling the window unreliable here would be crying wolf.  Recorded so the choice stays
# auditable.
SKIP_HOMOPOLYMER_SNV = "homopolymer_snv"

# zygosity codes stored in the sidecar
ZYG_HOM_ALT = "hom"
ZYG_HET = "het"


class Preflight(Exception):
    """A validation failure that must stop the build before anything is written."""


# ----------------------------------------------------------------------------------------------
# VCF reading


def open_maybe_gzip(path):
    """Open a plain or gzip/BGZF VCF for text reading."""
    with open(path, "rb") as probe:
        magic = probe.read(2)
    if magic == b"\x1f\x8b":
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def read_header(path):
    """Return (sample_names, contigs_declared) without consuming the record stream."""
    samples, contigs = [], []
    with open_maybe_gzip(path) as fh:
        for line in fh:
            if line.startswith("##contig="):
                field = line.split("ID=", 1)
                if len(field) == 2:
                    contigs.append(field[1].split(",", 1)[0].rstrip(">\n"))
            elif line.startswith("#CHROM"):
                cols = line.rstrip("\n").split("\t")
                if len(cols) > 9:
                    samples = cols[9:]
                break
            elif not line.startswith("#"):
                break
    return samples, contigs


def iter_records(path, sample_idx, chroms=None):
    """
    Stream (chrom, pos, ref, alt, info, fmt, sample_field) for one sample.

    The VCF is coordinate-sorted, so callers can rely on records arriving grouped by contig.
    """
    with open_maybe_gzip(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 10:
                continue
            chrom = cols[0]
            if chroms is not None and chrom not in chroms:
                continue
            yield chrom, int(cols[1]), cols[3], cols[4], cols[7], cols[8], cols[9 + sample_idx]


def info_int(info, key):
    """One integer INFO field, or None when absent or '.'."""
    for field in info.split(";"):
        if field.startswith(key) and field[len(key):len(key) + 1] == "=":
            value = field[len(key) + 1:]
            if value in (".", ""):
                return None
            try:
                return int(value)
            except ValueError:
                return None
    return None


def parse_genotype(fmt, sample_field):
    """
    Return (allele1, allele2, phased, phase_set).

    Alleles are ints, or None for a no-call.  phase_set is 0 when absent.
    """
    keys = fmt.split(":")
    vals = sample_field.split(":")
    rec = dict(zip(keys, vals))
    gt = rec.get("GT", "./.")
    phased = "|" in gt
    parts = gt.replace("|", "/").split("/")
    if len(parts) != 2:
        return None, None, phased, 0
    def _a(x):
        return None if x == "." else int(x)
    ps = rec.get("PS", "0")
    try:
        ps = int(ps)
    except ValueError:
        ps = 0
    return _a(parts[0]), _a(parts[1]), phased, ps


def classify(ref, alt, a1, a2, phased, haplotype, unphased_het):
    """
    Decide what to do with one record.

    Returns (action, zygosity_or_None, reason_or_None) where action is "apply" or "skip".
    """
    if alt in (".", "") or alt.startswith("<") or "," in alt or "*" in alt:
        # multiallelic or symbolic: the input is meant to be normalized+biallelic
        return "skip", None, SKIP_SYMBOLIC
    if a1 is None or a2 is None:
        return "skip", None, SKIP_NO_CALL
    if a1 == 0 and a2 == 0:
        return "skip", None, SKIP_REF_HOM
    if len(ref) != len(alt):
        if a1 == a2:                       # hom-alt indel: present on both haplotypes
            return "skip", ZYG_HOM_ALT, SKIP_INDEL_SELECTED_HAP
        if not phased:                     # cannot tell; assume it affects us
            return "skip", ZYG_HET, SKIP_INDEL_UNPHASED
        chosen = a1 if haplotype == 1 else a2
        return "skip", ZYG_HET, (
            SKIP_INDEL_SELECTED_HAP if chosen != 0 else SKIP_INDEL_OTHER_HAP
        )
    if a1 == a2:                                   # hom-alt: unambiguous
        return "apply", ZYG_HOM_ALT, None
    # heterozygous
    if not phased:
        if unphased_het == "apply_alt":
            return "apply", ZYG_HET, None
        return "skip", ZYG_HET, SKIP_UNPHASED_HET
    chosen = a1 if haplotype == 1 else a2
    if chosen == 0:
        return "skip", ZYG_HET, SKIP_OTHER_HAPLOTYPE
    return "apply", ZYG_HET, None


# ----------------------------------------------------------------------------------------------
# genome pickle IO


def genome_dir_for(genome_dir, genome_ver):
    return os.path.join(genome_dir, "fa_pickle", genome_ver)


def contig_path(genome_dir, genome_ver, chrom):
    return os.path.join(genome_dir_for(genome_dir, genome_ver), f"{chrom}.pk")


def list_contigs(genome_dir, genome_ver):
    d = genome_dir_for(genome_dir, genome_ver)
    if not os.path.isdir(d):
        raise Preflight(f"genome directory not found: {d}")
    return sorted(f[:-3] for f in os.listdir(d) if f.endswith(".pk"))


def load_contig(genome_dir, genome_ver, chrom):
    with open(contig_path(genome_dir, genome_ver, chrom), "rb") as fh:
        return pickle.load(fh)


def sha256_of(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


# ----------------------------------------------------------------------------------------------
# preflight


def preflight(args, samples, seq_cache):
    """
    Validate the VCF against the reference genome.  Raises Preflight on anything disqualifying.
    Returns a report dict.
    """
    report = {}

    if args.sample not in samples:
        raise Preflight(
            f"sample {args.sample!r} not in VCF.\nAvailable ({len(samples)}): "
            + ", ".join(samples)
        )
    sample_idx = samples.index(args.sample)
    report["sample"] = args.sample
    report["sample_index"] = sample_idx

    genome_contigs = set(list_contigs(args.genome_dir, args.genome_ver))
    report["genome_contigs"] = len(genome_contigs)

    # contig naming + REF concordance, sampled per contig
    seen = defaultdict(list)
    target = args.preflight_records_per_contig
    # The VCF is coordinate-sorted, so a contig is finished once the name changes.  We still
    # stream the whole file (contigs are only discoverable by reading), but stop *collecting*
    # from a contig at `target`, sampling evenly rather than taking the first N.
    stride = defaultdict(int)
    firstn = defaultdict(list)
    for chrom, pos, ref, alt, info, fmt, sfield in iter_records(args.vcf, sample_idx, args.chroms):
        stride[chrom] += 1
        if len(firstn[chrom]) < target:
            firstn[chrom].append((pos, ref))
        if len(seen[chrom]) < target and stride[chrom] % args.preflight_stride == 0:
            seen[chrom].append((pos, ref))
    # A contig with fewer than stride*target records samples nothing above (MT has 37, Y has 442
    # in the KOLF2.1J callset).  Backfill from the head of that contig, collected in the same
    # pass -- re-streaming a 4 M-record VCF once per small contig costs minutes.
    for chrom in list(stride):
        if not seen[chrom]:
            seen[chrom] = firstn[chrom]

    vcf_contigs = set(seen)
    report["vcf_contigs"] = sorted(vcf_contigs)
    missing = sorted(vcf_contigs - genome_contigs)
    report["vcf_contigs_absent_from_genome"] = missing
    if missing:
        hint = ""
        if any(c.startswith("chr") for c in missing) and not any(
            c.startswith("chr") for c in genome_contigs
        ):
            hint = (
                "\nThe VCF uses 'chr'-prefixed names and the genome does not. "
                "Rename with:  bcftools annotate --rename-chrs <map> "
            )
        raise Preflight(
            "VCF contigs absent from the genome: " + ", ".join(missing[:20]) + hint
        )

    checks = {}
    total_ok = total_n = 0
    for chrom, entries in sorted(seen.items()):
        seq = seq_cache(chrom)
        ok = 0
        examples = []
        for pos, ref in entries:
            got = seq[pos - 1 : pos - 1 + len(ref)].upper()
            if got == ref.upper():
                ok += 1
            elif len(examples) < 3:
                examples.append({"pos": pos, "vcf_ref": ref, "genome": got})
        checks[chrom] = {
            "checked": len(entries),
            "match": ok,
            "mismatch": len(entries) - ok,
            "examples": examples,
        }
        total_ok += ok
        total_n += len(entries)

    rate = (total_ok / total_n) if total_n else 0.0
    report["ref_concordance"] = {
        "checked": total_n,
        "match": total_ok,
        "rate": round(rate, 6),
        "per_contig": checks,
    }
    if total_n == 0:
        raise Preflight("no usable records found in the VCF")
    if rate < args.min_ref_concordance:
        bad = {c: v for c, v in checks.items() if v["mismatch"]}
        raise Preflight(
            f"REF-allele concordance {rate:.4f} is below --min_ref_concordance "
            f"{args.min_ref_concordance}.\nThis usually means the VCF was called against a "
            f"different assembly than {args.genome_ver}.\nOffending contigs: "
            + json.dumps(bad, indent=2)[:2000]
        )
    return report, sample_idx


# ----------------------------------------------------------------------------------------------
# build


def substitute_contig(seq_str, variants):
    """
    Apply length-preserving substitutions.

    variants: list of (pos1, ref, alt).  The ALT is written in the case of the REF base it
    replaces, so dna_sm soft-masking survives.
    Returns (new_seq_str, n_applied).
    """
    # a bytearray rather than list(seq_str): chr1 is 249 Mbp, and a list of that many single-char
    # objects costs ~2 GB of pointers where the bytearray costs 249 MB
    buf = bytearray(seq_str, "ascii")
    n = 0
    for pos, ref, alt in variants:
        i = pos - 1
        for k, base in enumerate(alt.encode("ascii")):
            # 0x20 is the ASCII case bit; carry the REF base's case so dna_sm masking survives
            buf[i + k] = (base | 0x20) if (buf[i + k] & 0x20) else (base & ~0x20)
        n += 1
    return buf.decode("ascii"), n


def write_contig(out_dir, chrom, record, seq_str, provenance):
    os.makedirs(out_dir, exist_ok=True)
    new_desc = record.description
    if provenance and provenance not in new_desc:
        new_desc = f"{new_desc} | {provenance}"
    out = SeqRecord(
        Seq(seq_str),
        id=record.id,
        name=record.name,
        description=new_desc,
    )
    tmp = os.path.join(out_dir, f".{chrom}.pk.tmp")
    with open(tmp, "wb") as fh:
        pickle.dump(out, fh, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(tmp, os.path.join(out_dir, f"{chrom}.pk"))


def link_unchanged_contig(genome_dir, base, out_ver, chrom):
    """
    A contig with no applied substitutions is byte-identical to the reference, so hardlink it
    rather than writing a second copy.  Falls back to a real copy across filesystems.
    Returns "link" or "copy".
    """
    src = contig_path(genome_dir, base, chrom)
    dst = contig_path(genome_dir, out_ver, chrom)
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    if os.path.exists(dst):
        os.remove(dst)
    try:
        os.link(src, dst)
        return "link"
    except OSError:
        import shutil

        shutil.copy2(src, dst)
        return "copy"


def read_nocall_bed(path, chroms=None):
    """
    Parse the no-call BED into {chrom: [[start1, end1], ...]}, merged and sorted.

    These are positions where the sample's genotype is *unknown*, not reference.  The single-sample
    VCF was made with `bcftools view --exclude-uncalled`, so they are simply absent from the record
    stream: without this file the pipeline emits reference there and reports the window reliable.
    """
    by_chrom = defaultdict(list)
    with open_maybe_gzip(path) as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 3:
                continue
            chrom = cols[0]
            if chroms is not None and chrom not in chroms:
                continue
            by_chrom[chrom].append((int(cols[1]) + 1, int(cols[2])))   # BED is 0-based half-open

    out = {}
    for chrom, intervals in by_chrom.items():
        intervals.sort()
        merged = []
        for start, end in intervals:
            # `+ 1` also merges book-ended runs, which is what a caller asking "is this window
            # covered" means by a region
            if merged and start <= merged[-1][1] + 1:
                merged[-1][1] = max(merged[-1][1], end)
            else:
                merged.append([start, end])
        out[chrom] = merged
    return out


def write_sidecar(path, payload):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "wb") as fh:
        pickle.dump(payload, fh, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(tmp, path)


def build(args, sample_idx, seq_cache, preflight_report, nocall=None):
    out_fa = genome_dir_for(args.genome_dir, args.out_genome_ver)
    out_vs = os.path.join(args.genome_dir, "variant_sets", args.out_genome_ver)
    provenance = (
        f"variant-substituted: {args.out_genome_ver} "
        f"[{args.sample} hap{args.haplotype}]"
    )

    per_contig = {}
    pending = []          # (pos, ref, alt, zyg, ps) to apply for the current contig
    skipped = []          # (pos, ref, alt, zyg, reason)
    current = None
    order = []

    def flush(chrom):
        if chrom is None:
            return
        record = load_contig(args.genome_dir, args.genome_ver, chrom)
        seq_str = str(record.seq)

        # re-verify every REF we are about to write; a mismatch demotes it to skipped
        verified = []
        for pos, ref, alt, zyg, ps in pending:
            if seq_str[pos - 1 : pos - 1 + len(ref)].upper() == ref.upper():
                verified.append((pos, ref, alt, zyg, ps))
            else:
                skipped.append((pos, ref, alt, zyg, SKIP_REF_MISMATCH))

        # Two substitutions cannot occupy the same base.  The sort is stable and the stream is
        # already coordinate-sorted, so the first record at a position wins and any overlapping
        # follower is recorded as unrepresented rather than silently overwriting it.
        verified.sort(key=lambda t: t[0])
        ap = []
        prev_end = 0
        for pos, ref, alt, zyg, ps in verified:
            if pos <= prev_end:
                skipped.append((pos, ref, alt, zyg, SKIP_OVERLAPPING_SUB))
                continue
            ap.append((pos, ref, alt, zyg, ps))
            prev_end = pos + len(ref) - 1

        # variant_annot bisects these arrays, so they have to be sorted; the REF-mismatch and
        # overlap demotions above were appended out of order
        skipped.sort(key=lambda t: t[0])

        new_seq, n_applied = substitute_contig(
            seq_str, [(p, r, a) for p, r, a, _, _ in ap]
        )
        assert len(new_seq) == len(seq_str), "substitution changed contig length"

        if not args.no_materialize and n_applied:
            write_contig(out_fa, chrom, record, new_seq, provenance)

        write_sidecar(
            os.path.join(out_vs, "applied", f"{chrom}.pk"),
            {
                "chrom": chrom,
                "pos": array.array("i", [p for p, _, _, _, _ in ap]),
                "ref": [r for _, r, _, _, _ in ap],
                "alt": [a for _, _, a, _, _ in ap],
                "zygosity": [z for _, _, _, z, _ in ap],
                "phase_set": array.array("i", [s for _, _, _, _, s in ap]),
            },
        )
        write_sidecar(
            os.path.join(out_vs, "skipped", f"{chrom}.pk"),
            {
                "chrom": chrom,
                "pos": array.array("i", [p for p, _, _, _, _ in skipped]),
                "ref": [r for _, r, _, _, _ in skipped],
                "alt": [a for _, _, a, _, _ in skipped],
                "zygosity": [z for _, _, _, z, _ in skipped],
                "reason": [x for _, _, _, _, x in skipped],
            },
        )
        reasons = Counter(x for _, _, _, _, x in skipped)
        affecting = sum(
            reasons.get(k, 0)
            for k in (SKIP_INDEL_SELECTED_HAP, SKIP_INDEL_UNPHASED, SKIP_UNPHASED_HET,
                      SKIP_REF_MISMATCH, SKIP_SYMBOLIC, SKIP_OVERLAPPING_SUB)
        )
        per_contig[chrom] = {
            "length": len(seq_str),
            "applied": n_applied,
            "skipped": len(skipped),
            "skipped_by_reason": dict(reasons),
            "sequence_known_wrong": affecting,
            # a contig with nothing applied is byte-identical to the reference and is
            # hardlinked rather than rewritten, so it carries no provenance string
            "linked_unchanged": n_applied == 0,
        }
        print(
            f"  {chrom:>12}  len={len(seq_str):>11,}  applied={n_applied:>8,}  "
            f"unrepresented={affecting:>7,}  "
            f"({', '.join(f'{k}={v}' for k, v in reasons.most_common())})",
            flush=True,
        )
        pending.clear()
        skipped.clear()

    hp_min = args.homopolymer_ref_min_len
    hp_tagged = 0
    for chrom, pos, ref, alt, info, fmt, sfield in iter_records(args.vcf, sample_idx, args.chroms):
        if chrom != current:
            flush(current)
            current = chrom
            order.append(chrom)
        a1, a2, phased, ps = parse_genotype(fmt, sfield)
        action, zyg, reason = classify(
            ref, alt, a1, a2, phased, args.haplotype, args.unphased_het
        )
        if hp_min:
            hp_len = info_int(info, "HPLEN")
            if hp_len is not None:
                hp_tagged += 1
                # Only substitutions ever reach "apply", so this cannot change an indel's
                # classification -- deliberately.  Ti/Tv is an SNV statistic; it says nothing
                # about whether an indel inside a homopolymer is real, and homopolymers are
                # exactly where real germline indels concentrate.
                if hp_len >= hp_min and action == "apply":
                    action, reason = "skip", SKIP_HOMOPOLYMER_SNV
        if action == "apply":
            pending.append((pos, ref, alt, zyg, ps))
        elif reason not in (SKIP_REF_HOM, SKIP_NO_CALL):
            # ref-hom and no-call are the bulk and carry no design information; drop them
            skipped.append((pos, ref, alt, zyg, reason))
    flush(current)

    if hp_min and not hp_tagged:
        print(
            f"  WARNING: --homopolymer_ref_min_len {hp_min} was given but not one record carries "
            f"an HPLEN INFO tag, so the policy did nothing.  Annotate the VCF first:\n"
            f"           bcftools annotate -a <homopolymer.bed.gz> "
            f"-c CHROM,FROM,TO,INFO/HPLEN,INFO/HPBASE -h hp.hdr",
            flush=True,
        )

    # No-call regions: unknown genotype, emitted as reference.  Written as a third shard kind so
    # variant_annot can flag any window overlapping one without needing to read a BED.
    nocall_regions = nocall_bp = 0
    if nocall:
        for chrom, intervals in sorted(nocall.items()):
            write_sidecar(
                os.path.join(out_vs, "nocall", f"{chrom}.pk"),
                {
                    "chrom": chrom,
                    "start": array.array("i", [s for s, _ in intervals]),
                    "end": array.array("i", [e for _, e in intervals]),
                },
            )
            nocall_regions += len(intervals)
            nocall_bp += sum(e - s + 1 for s, e in intervals)
        print(f"  no-call regions: {nocall_regions:,} spanning {nocall_bp:,} bp")

    # Contigs absent from the VCF (or present but with nothing applied) are identical to the
    # reference.  The materialized genome must still contain all of them, or get_seq() will fail
    # on any scaffold.  Hardlink rather than duplicate ~2.9 GB.
    linked = []
    if not args.no_materialize:
        all_contigs = list_contigs(args.genome_dir, args.genome_ver)
        wanted = [c for c in all_contigs if args.chroms is None or c in args.chroms]
        for chrom in wanted:
            if per_contig.get(chrom, {}).get("applied"):
                continue
            link_unchanged_contig(args.genome_dir, args.genome_ver, args.out_genome_ver, chrom)
            linked.append(chrom)
        if linked:
            print(f"  linked {len(linked)} unchanged contig(s) not present in the VCF")

    manifest = {
        "out_genome_ver": args.out_genome_ver,
        "base_genome_ver": args.genome_ver,
        "source_vcf": os.path.abspath(args.vcf),
        "source_vcf_sha256": sha256_of(args.vcf) if args.checksum else None,
        "sample": args.sample,
        "haplotype": args.haplotype,
        "unphased_het_policy": args.unphased_het,
        "homopolymer_ref_min_len": args.homopolymer_ref_min_len,
        "substitutions": "length-preserving only (SNV/MNV); indels recorded, never applied",
        "materialized": not args.no_materialize,
        "built": datetime.datetime.now().isoformat(timespec="seconds"),
        "builder": os.path.basename(__file__),
        "contigs": order,
        "contigs_linked_unchanged": linked,
        "nocall_bed": os.path.abspath(args.nocall_bed) if args.nocall_bed else None,
        "nocall_bed_sha256": (
            sha256_of(args.nocall_bed) if (args.nocall_bed and args.checksum) else None
        ),
        "nocall_regions": nocall_regions,
        "nocall_bp": nocall_bp,
        "per_contig": per_contig,
        "preflight": preflight_report,
    }
    os.makedirs(out_vs, exist_ok=True)
    with open(os.path.join(out_vs, "manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=2)

    with open(os.path.join(out_vs, "summary.tsv"), "w") as fh:
        fh.write("contig\tlength\tapplied\tskipped\tsequence_known_wrong\tskipped_detail\n")
        for c in order:
            v = per_contig[c]
            detail = ";".join(f"{k}={n}" for k, n in sorted(v["skipped_by_reason"].items()))
            fh.write(
                f"{c}\t{v['length']}\t{v['applied']}\t{v['skipped']}"
                f"\t{v['sequence_known_wrong']}\t{detail}\n"
            )

    return manifest


# ----------------------------------------------------------------------------------------------


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--vcf", required=True, help="bgzipped or plain VCF")
    p.add_argument("--sample", required=True, help="sample column to use")
    p.add_argument("--genome_ver", default="GRCh38", help="source genome under --genome_dir")
    p.add_argument("--out_genome_ver", help="name of the substituted genome, e.g. GRCh38_KOLF2.1J_hap1")
    p.add_argument("--genome_dir", default="genome_files")
    p.add_argument("--haplotype", type=int, choices=(1, 2), default=1)
    p.add_argument(
        "--unphased_het",
        choices=("skip", "apply_alt"),
        default="skip",
        help="what to do with unphased heterozygous sites (default: skip and record)",
    )
    p.add_argument("--chroms", help="comma-separated subset, e.g. 21,22")
    p.add_argument(
        "--homopolymer_ref_min_len",
        type=int,
        default=7,
        help="Do not apply a substitution whose HPLEN INFO tag is >= this; emit the reference base "
        "instead (default: 7, the run length at which SNV Ti/Tv drops below the 0.5 random-error "
        "floor).  0 disables.  Inert on a VCF with no HPLEN tag.  Indels are unaffected -- they "
        "are never applied anyway, and Ti/Tv says nothing about indel accuracy.",
    )
    p.add_argument(
        "--nocall_bed",
        help="BED of positions where this sample has no genotype call.  Single-sample VCFs are "
        "usually made with --exclude-uncalled, which deletes them; without this file those "
        "windows are emitted as reference and reported reliable.",
    )
    p.add_argument("--validate_only", action="store_true", help="run preflight, write nothing")
    p.add_argument(
        "--no_materialize",
        action="store_true",
        help="write only the sidecar, skip the ~2.9 GB substituted genome",
    )
    p.add_argument("--min_ref_concordance", type=float, default=0.99)
    p.add_argument("--preflight_records_per_contig", type=int, default=400)
    p.add_argument("--preflight_stride", type=int, default=7,
                   help="sample every Nth record per contig during preflight (spreads the check)")
    p.add_argument("--checksum", action="store_true", help="sha256 the source VCF for the manifest")
    p.add_argument("--force", action="store_true", help="overwrite existing output")
    a = p.parse_args(argv)
    a.chroms = set(a.chroms.split(",")) if a.chroms else None
    if not a.validate_only and not a.out_genome_ver:
        p.error("--out_genome_ver is required unless --validate_only")
    return a


def main(argv=None):
    args = parse_args(argv)

    if not os.path.isfile(args.vcf):
        sys.exit(f"ERROR: VCF not found: {args.vcf}")

    if args.out_genome_ver and not args.validate_only and not args.force:
        existing = genome_dir_for(args.genome_dir, args.out_genome_ver)
        if os.path.isdir(existing) and os.listdir(existing):
            sys.exit(f"ERROR: {existing} already exists and is not empty; use --force")

    # one-chromosome-at-a-time cache: preflight walks contigs in sorted order, build in VCF order
    _cache = {}

    def seq_cache(chrom):
        if chrom not in _cache:
            _cache.clear()
            _cache[chrom] = str(load_contig(args.genome_dir, args.genome_ver, chrom).seq)
        return _cache[chrom]

    if args.nocall_bed and not os.path.isfile(args.nocall_bed):
        sys.exit(f"ERROR: --nocall_bed not found: {args.nocall_bed}")

    samples, _declared = read_header(args.vcf)
    print(f"VCF: {args.vcf}")
    print(f"  samples in file: {len(samples)}")

    try:
        report, sample_idx = preflight(args, samples, seq_cache)
    except Preflight as e:
        sys.exit(f"\nPREFLIGHT FAILED\n{e}\n")

    rc = report["ref_concordance"]
    print("\npreflight")
    print(f"  sample index          : {report['sample_index']}")
    print(f"  contigs in VCF        : {len(report['vcf_contigs'])}")
    print(f"  contigs in genome     : {report['genome_contigs']}")
    print(f"  REF concordance       : {rc['match']}/{rc['checked']} = {rc['rate']:.4f}")
    for c, v in sorted(rc["per_contig"].items()):
        flag = "" if v["mismatch"] == 0 else f"   <-- {v['mismatch']} MISMATCH {v['examples']}"
        print(f"      {c:>12}  {v['match']}/{v['checked']}{flag}")

    if args.validate_only:
        print("\n--validate_only: nothing written.")
        return 0

    nocall = None
    if args.nocall_bed:
        nocall = read_nocall_bed(args.nocall_bed, args.chroms)
        print(f"\nno-call BED: {args.nocall_bed}")
        print(f"  contigs {len(nocall)}, merged regions "
              f"{sum(len(v) for v in nocall.values()):,}")
        stray = sorted(set(nocall) - set(list_contigs(args.genome_dir, args.genome_ver)))
        if stray:
            sys.exit(f"ERROR: --nocall_bed contigs absent from the genome: {', '.join(stray[:20])}")

    _cache.clear()
    print(f"\nbuilding {args.out_genome_ver}  (haplotype {args.haplotype}, "
          f"unphased_het={args.unphased_het}, materialize={not args.no_materialize})")
    print("  " + (
        f"substitutions inside homopolymer runs >= {args.homopolymer_ref_min_len} bp "
        f"emit the reference base"
        if args.homopolymer_ref_min_len else "homopolymer policy off"
    ))
    manifest = build(args, sample_idx, seq_cache, report, nocall=nocall)

    tot_a = sum(v["applied"] for v in manifest["per_contig"].values())
    tot_s = sum(v["skipped"] for v in manifest["per_contig"].values())
    print(f"\ndone.  applied={tot_a:,}  recorded-but-not-applied={tot_s:,}")
    print(f"  genome  : {genome_dir_for(args.genome_dir, args.out_genome_ver)}")
    print(f"  sidecar : {os.path.join(args.genome_dir, 'variant_sets', args.out_genome_ver)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
