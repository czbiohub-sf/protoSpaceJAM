"""
Runtime variant layer for protoSpaceJAM.

Two jobs:

  1. patch()          - substitute variant alleles into a reference slice as it is fetched,
                        so donors and guide sequences are personalized without materializing
                        a whole genome.  This is the counterpart to the substituted pickles
                        written by precompute/variant_sets/4_make_variant_genome.py, and must produce
                        byte-identical results (see 6_validate_runtime_patch.py).

  2. variants_in()    - report which variants intersect a window, and whether the emitted
                        sequence is trustworthy there.  Used for guide penalties, result.csv
                        columns and GenBank features.

Only length-preserving substitutions are ever applied; see 4_make_variant_genome.py for why.
An indel on the *selected* haplotype means the emitted sequence is knowingly wrong in that
window - callers must surface that rather than silently designing against it.
"""

import bisect
import json
import os
import pickle

# reasons, mirrored from precompute/variant_sets/4_make_variant_genome.py
SKIP_INDEL_SELECTED_HAP = "indel_selected_hap"
SKIP_INDEL_OTHER_HAP = "indel_other_hap"
SKIP_INDEL_UNPHASED = "indel_unphased"
SKIP_UNPHASED_HET = "unphased_het"
SKIP_OTHER_HAPLOTYPE = "other_haplotype"
SKIP_REF_MISMATCH = "ref_mismatch"
SKIP_SYMBOLIC = "symbolic_allele"
SKIP_OVERLAPPING_SUB = "overlapping_substitution"

#: A substitution the builder declined to write because it sits inside a long homopolymer run.
#: Deliberately **not** in SEQUENCE_UNRELIABLE: in runs >=7 bp the SNV Ti/Tv is 0.43-0.45, below
#: the 0.5 random-error floor, so the reference base is the more likely truth and the emitted
#: sequence is the best guess available.  Kept in the sidecar so the decision is auditable, but
#: not surfaced as a warning -- see 4_make_variant_genome.py --homopolymer_ref_min_len.
SKIP_HOMOPOLYMER_SNV = "homopolymer_snv"

#: A region where the sample has no genotype call at all.  Single-sample VCFs are extracted with
#: `bcftools view --exclude-uncalled`, so these positions are simply absent from the record stream
#: and the pipeline would otherwise emit reference and call the window reliable.  Unknown is not
#: reference: the builder writes them into a `nocall/` shard from the caller's BED.
SKIP_NO_CALL_REGION = "no_call_region"

#: skipped-variant reasons that mean the emitted sequence does not represent the
#: selected haplotype in that window
SEQUENCE_UNRELIABLE = frozenset(
    {SKIP_INDEL_SELECTED_HAP, SKIP_INDEL_UNPHASED, SKIP_UNPHASED_HET,
     SKIP_REF_MISMATCH, SKIP_SYMBOLIC, SKIP_OVERLAPPING_SUB, SKIP_NO_CALL_REGION}
)


class Variant:
    """One variant intersecting a queried window."""

    __slots__ = (
        "chrom", "pos", "ref", "alt", "zygosity", "phase_set", "applied", "reason", "_end",
    )

    def __init__(self, chrom, pos, ref, alt, zygosity, phase_set=0, applied=True, reason=None,
                 end=None):
        self.chrom = chrom
        self.pos = pos                  # 1-based
        self.ref = ref
        self.alt = alt
        self.zygosity = zygosity        # "hom" | "het" | None
        self.phase_set = phase_set
        self.applied = applied          # True = written into the sequence
        self.reason = reason            # why not, when applied is False
        # a no-call region has no REF allele to measure, so its span is carried explicitly
        self._end = end

    @property
    def is_indel(self):
        return len(self.ref) != len(self.alt)

    @property
    def is_region(self):
        """True for a span (no-call) rather than a called allele."""
        return self.reason == SKIP_NO_CALL_REGION

    @property
    def end(self):
        """1-based inclusive end of the REF footprint."""
        if self._end is not None:
            return self._end
        return self.pos + len(self.ref) - 1

    @property
    def makes_sequence_unreliable(self):
        return (not self.applied) and self.reason in SEQUENCE_UNRELIABLE

    def __repr__(self):
        state = "applied" if self.applied else f"skipped:{self.reason}"
        return f"{self.chrom}:{self.pos} {self.ref}>{self.alt} {self.zygosity} [{state}]"

    def short(self):
        if self.is_region:
            span = self.pos if self.end == self.pos else f"{self.pos}-{self.end}"
            return f"{self.chrom}:{span}(no-call)"
        z = self.zygosity or "?"
        return f"{self.chrom}:{self.pos}{self.ref}>{self.alt}({z})"


class VariantSet:
    """
    Lazily-loaded per-contig variant shards for one sample/haplotype.

    Instantiate once per run and pass down to get_seq() as `variant_ctx`.
    """

    def __init__(self, set_name, genome_dir="genome_files", cache_contigs=2):
        self.set_name = set_name
        self.root = os.path.join(genome_dir, "variant_sets", set_name)
        manifest_path = os.path.join(self.root, "manifest.json")
        if not os.path.isfile(manifest_path):
            raise FileNotFoundError(
                f"variant set {set_name!r} not found (expected {manifest_path})"
            )
        with open(manifest_path) as fh:
            self.manifest = json.load(fh)
        self.base_genome_ver = self.manifest["base_genome_ver"]
        self.sample = self.manifest["sample"]
        self.haplotype = self.manifest["haplotype"]
        self._cache = {}
        self._cache_contigs = max(1, cache_contigs)

    # -- shard loading ---------------------------------------------------------------------

    def _load(self, chrom):
        if chrom in self._cache:
            return self._cache[chrom]
        shard = {"applied": None, "skipped": None, "nocall": None}
        for kind in ("applied", "skipped", "nocall"):
            path = os.path.join(self.root, kind, f"{chrom}.pk")
            if os.path.isfile(path):
                with open(path, "rb") as fh:
                    shard[kind] = pickle.load(fh)
        # An MNV can begin *before* the slice and still overlap it, so queries must widen the
        # left edge by the longest REF present.  Precompute it once per contig.
        for kind in ("applied", "skipped"):
            if shard[kind]:
                shard[kind]["_max_ref"] = max(
                    (len(r) for r in shard[kind]["ref"]), default=1
                )
        # contigs with no variants have no shard at all; represent as empty
        if len(self._cache) >= self._cache_contigs:
            self._cache.clear()
        self._cache[chrom] = shard
        return shard

    @staticmethod
    def _range(positions, start, end):
        """Index range of positions within [start, end], both 1-based inclusive."""
        lo = bisect.bisect_left(positions, start)
        hi = bisect.bisect_right(positions, end)
        return lo, hi

    # -- the hot path ----------------------------------------------------------------------

    def patch(self, chrom, start, seq):
        """
        Apply substitutions to `seq`, which is the reference slice beginning at 1-based `start`.

        Case is preserved: the ALT base is written in the case of the REF base it replaces, so
        dna_sm soft-masking survives, exactly as 4_make_variant_genome.py does it.
        """
        shard = self._load(chrom)["applied"]
        if not shard or not len(shard["pos"]):
            return seq
        end = start + len(seq) - 1
        widen = shard.get("_max_ref", 1) - 1
        lo, hi = self._range(shard["pos"], start - widen, end)
        if lo == hi:
            return seq

        chars = None
        positions, alts = shard["pos"], shard["alt"]
        for i in range(lo, hi):
            off = positions[i] - start
            alt = alts[i]
            if off + len(alt) <= 0:
                continue                      # widened past it; ends before the slice begins
            if chars is None:
                chars = list(seq)
            for k, base in enumerate(alt):
                j = off + k
                if 0 <= j < len(chars):       # a substitution may straddle either slice edge
                    chars[j] = base.lower() if chars[j].islower() else base.upper()
        return seq if chars is None else "".join(chars)

    # -- reporting -------------------------------------------------------------------------

    def variants_in(self, chrom, start, end, include_skipped=True):
        """
        All variants whose REF footprint intersects [start, end], 1-based inclusive.

        Skipped variants are included by default because an indel on the selected haplotype is
        the single most important thing to report about a window.
        """
        shard = self._load(chrom)
        out = []

        applied = shard["applied"]
        if applied and len(applied["pos"]):
            # an applied MNV can start before `start` and still overlap it, same as patch()
            widen = applied.get("_max_ref", 1) - 1
            lo, hi = self._range(applied["pos"], start - widen, end)
            for i in range(lo, hi):
                if applied["pos"][i] + len(applied["ref"][i]) - 1 < start:
                    continue
                out.append(
                    Variant(
                        chrom,
                        applied["pos"][i],
                        applied["ref"][i],
                        applied["alt"][i],
                        applied["zygosity"][i],
                        applied["phase_set"][i],
                        applied=True,
                    )
                )

        if include_skipped:
            skipped = shard["skipped"]
            if skipped and len(skipped["pos"]):
                # a skipped indel can start before `start` and still overlap it
                widen = skipped.get("_max_ref", 1) - 1
                lo, hi = self._range(skipped["pos"], start - widen, end)
                for i in range(lo, hi):
                    pos, ref = skipped["pos"][i], skipped["ref"][i]
                    if pos + len(ref) - 1 < start:
                        continue
                    out.append(
                        Variant(
                            chrom,
                            pos,
                            ref,
                            skipped["alt"][i],
                            skipped["zygosity"][i],
                            0,
                            applied=False,
                            reason=skipped["reason"][i],
                        )
                    )

        if include_skipped:
            for lo_r, hi_r in self._nocall_overlaps(shard["nocall"], start, end):
                out.append(
                    Variant(
                        chrom, lo_r, ".", ".", None, 0,
                        applied=False, reason=SKIP_NO_CALL_REGION, end=hi_r,
                    )
                )

        out.sort(key=lambda v: (v.pos, not v.applied))
        return out

    @staticmethod
    def _nocall_overlaps(shard, start, end):
        """(lo, hi) pairs from the merged no-call shard intersecting [start, end]."""
        if not shard or not len(shard["start"]):
            return []
        starts, ends = shard["start"], shard["end"]
        # the regions are merged and non-overlapping, so both arrays are sorted
        lo = bisect.bisect_left(ends, start)        # first region ending at or after `start`
        hi = bisect.bisect_right(starts, end)       # first region starting after `end`
        return [(starts[i], ends[i]) for i in range(lo, hi)]

    def describe_window(self, chrom, start, end):
        """Summary used for guide penalties and result.csv columns."""
        vs = self.variants_in(chrom, start, end)
        applied = [v for v in vs if v.applied]
        unreliable = [v for v in vs if v.makes_sequence_unreliable]
        return {
            "n_applied": len(applied),
            "n_unreliable": len(unreliable),
            "n_nocall": sum(1 for v in unreliable if v.is_region),
            "has_hom": any(v.zygosity == "hom" for v in applied),
            "has_het": any(v.zygosity == "het" for v in applied),
            "sequence_reliable": not unreliable,
            "applied": applied,
            "unreliable": unreliable,
            "summary": "; ".join(v.short() for v in vs) if vs else "",
        }


def load_variant_set(set_name, genome_dir="genome_files", expected_genome_ver=None):
    """Convenience loader with the base-genome consistency check."""
    vset = VariantSet(set_name, genome_dir=genome_dir)
    if expected_genome_ver and vset.base_genome_ver != expected_genome_ver:
        raise ValueError(
            f"variant set {set_name!r} was built against {vset.base_genome_ver}, "
            f"but --genome_ver is {expected_genome_ver}"
        )
    return vset


# ==================================================================================================
# guide-level assessment
#
# The precomputed gRNA tables are keyed to reference coordinates and carry reference sequence and
# reference off-target scores.  A guide whose PAM is destroyed in this cell line is still in the
# table, still scored, and still looks fine.  Everything below exists to catch that.
# ==================================================================================================

#: PAM-proximal bases of the protospacer.  Mismatches here are far more disruptive than distal ones.
SEED_LEN = 12

#: PAMs that sit 5' of the protospacer (Cas12a-family) rather than 3' (Cas9-family)
FIVE_PRIME_PAMS = frozenset({"TTTV", "TTTN"})

#: weight applied to a guide, by where the variant falls and whether both alleles carry it.
#: het is penalized *harder* than hom-alt in the protospacer: with hom-alt the personalized guide
#: matches both alleles, with het no single guide can.
GUIDE_WEIGHTS = {
    ("pam", "hom"): 0.0,       # PAM gone in these cells -- guide is dead
    ("pam", "het"): 0.0,       # cuts one allele only
    ("seed", "hom"): 0.5,
    ("seed", "het"): 0.25,
    ("distal", "hom"): 0.9,
    ("distal", "het"): 0.75,
}

#: a variant we could not represent *and* that changes length invalidates the coordinate model the
#: whole design rests on, so the guide is unusable regardless of where in the window it falls
WEIGHT_COORDINATE_MODEL_INVALID = 0.0


def pam_is_5prime(pam):
    return str(pam).upper() in FIVE_PRIME_PAMS


def guide_subwindows(start, end, strand, pam_len=3, five_prime_pam=False):
    """
    Split a precomputed-table guide row into PAM / seed / distal genomic windows.

    `start` and `end` are the table's own columns: 1-based inclusive, and start > end on the
    minus strand.  Their span is the whole protospacer+PAM (23 nt for NGG/NGA, 24 for TTTV), so
    the protospacer length is derived rather than assumed.

    Returns 1-based inclusive (lo, hi) tuples keyed "guide", "pam", "seed"; "seed" is None when
    the protospacer is shorter than SEED_LEN.  Windows are genomic, so lo <= hi on both strands.
    """
    lo, hi = (start, end) if start <= end else (end, start)
    span = hi - lo + 1
    plus = str(strand) in ("+", "1")

    def to_genomic(off_start, off_end):
        """offsets are 0-based from the 5' end of the guide, in guide orientation"""
        if plus:
            return (lo + off_start, lo + off_end)
        return (hi - off_end, hi - off_start)

    proto_len = span - pam_len
    if five_prime_pam:
        pam = to_genomic(0, pam_len - 1)
        seed = to_genomic(pam_len, pam_len + SEED_LEN - 1) if proto_len >= SEED_LEN else None
    else:
        pam = to_genomic(span - pam_len, span - 1)
        seed = to_genomic(proto_len - SEED_LEN, proto_len - 1) if proto_len >= SEED_LEN else None

    return {"guide": (lo, hi), "pam": pam, "seed": seed}


def region_of(variant, windows):
    """Which part of the guide a variant's REF footprint touches; PAM wins, then seed."""
    v_lo, v_hi = variant.pos, variant.end

    def hits(win):
        return win is not None and v_lo <= win[1] and v_hi >= win[0]

    if hits(windows["pam"]):
        return "pam"
    if hits(windows["seed"]):
        return "seed"
    return "distal"


def assess_guide(vset, chrom, start, end, strand, pam="NGG"):
    """
    Score one precomputed-table guide row against the variant set.

    Returns a dict with:
        weight       multiplier for final_weight; 1.0 when the guide is untouched
        stale        True when the reference MIT/CFD scores no longer describe this guide
        flags        short ";"-joined reasons, for the result.csv `variant_warnings` column
        summary      ";"-joined variant shorthand, for `variants_in_protospacer_PAM`
        variants     the Variant objects, for variants_report.csv

    `stale` is deliberately independent of `weight`: a distal hom-alt substitution barely moves
    the weight but still changes the off-target profile, and we cannot recompute it without
    re-running bwa.  Never present a reference score as if it still applied.
    """
    windows = guide_subwindows(
        start, end, strand, pam_len=len(str(pam)), five_prime_pam=pam_is_5prime(pam)
    )
    lo, hi = windows["guide"]
    variants = vset.variants_in(chrom, lo, hi)
    if not variants:
        return {
            "weight": 1.0,
            "stale": False,
            "flags": "",
            "summary": "",
            "variants": [],
        }

    weight = 1.0
    flags = []
    for v in variants:
        region = region_of(v, windows)
        if v.applied:
            zyg = v.zygosity if v.zygosity in ("hom", "het") else "het"  # unknown -> conservative
            w = GUIDE_WEIGHTS[(region, zyg)]
            flags.append(f"{region}_{zyg}")
        elif v.makes_sequence_unreliable:
            if v.is_indel:
                # the emitted sequence is the wrong length here; start/end/cut positions, the
                # phase list and every `pos - Lstart` index downstream are out of register
                w = WEIGHT_COORDINATE_MODEL_INVALID
            else:
                # a substitution we could not phase: same footprint, unknown allele -> treat as het
                w = GUIDE_WEIGHTS[(region, "het")]
            flags.append(f"{region}_{v.reason}")
        else:
            continue  # other-haplotype indels: hap1 really is reference here
        weight = min(weight, w)

    reportable = [v for v in variants if v.applied or v.makes_sequence_unreliable]
    return {
        "weight": weight,
        "stale": bool(reportable),
        "flags": ";".join(sorted(set(flags))),
        "summary": ";".join(v.short() for v in reportable),
        "variants": reportable,
    }


def list_variant_sets(genome_dir="genome_files"):
    """
    Every built variant set, as plain dicts for the portal's /api/variant_sets/ endpoint.

    Reads the manifests only; never touches the shards.
    """
    root = os.path.join(genome_dir, "variant_sets")
    out = []
    if not os.path.isdir(root):
        return out
    for name in sorted(os.listdir(root)):
        manifest_path = os.path.join(root, name, "manifest.json")
        if not os.path.isfile(manifest_path):
            continue
        try:
            with open(manifest_path) as fh:
                m = json.load(fh)
        except (ValueError, OSError):
            continue
        per_contig = m.get("per_contig", {}) or {}
        applied = sum(v.get("applied", 0) for v in per_contig.values())
        recorded = sum(v.get("skipped", 0) for v in per_contig.values())
        by_reason = {}
        for v in per_contig.values():
            for reason, n in (v.get("skipped_by_reason") or {}).items():
                by_reason[reason] = by_reason.get(reason, 0) + n
        indels = sum(
            by_reason.get(r, 0)
            for r in (SKIP_INDEL_SELECTED_HAP, SKIP_INDEL_OTHER_HAP, SKIP_INDEL_UNPHASED)
        )
        span = sum(v.get("length", 0) for v in per_contig.values())
        out.append(
            {
                "name": name,
                "base_genome_ver": m.get("base_genome_ver"),
                "sample": m.get("sample"),
                "haplotype": m.get("haplotype"),
                "materialized": m.get("materialized", False),
                "built": m.get("built"),
                "applied": applied,
                "recorded_not_applied": recorded,
                "indels_recorded": indels,
                "sequence_known_wrong": sum(
                    v.get("sequence_known_wrong", 0) for v in per_contig.values()
                ),
                "nocall_regions": m.get("nocall_regions", 0),
                # substituted bases per kb of the contigs the VCF covers; the headline number for
                # "how personalized is this genome"
                "applied_per_kb": round(applied / (span / 1000.0), 3) if span else None,
                "source_vcf": os.path.basename(m.get("source_vcf") or ""),
            }
        )
    return out
