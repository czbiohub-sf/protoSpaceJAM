#!/usr/bin/env python
"""
Validate a variant-substituted genome built by 4_make_variant_genome.py.

Checks (per contig):

  1. lossless  - length, id, name unchanged; provenance recorded in the description
  2. applied   - every sidecar 'applied' position read REF before and reads ALT after,
                 with dna_sm soft-mask case preserved
  3. minimal   - the ONLY differences between reference and substituted sequence are the
                 applied substitutions (nothing else was touched)
  4. skipped   - every sidecar 'skipped' variant is untouched, except where its REF span
                 legitimately overlaps a different applied variant (common: a het indel on
                 one haplotype spanning a het SNV on the other)
  5. policy    - a policy the manifest claims was in force actually left a mark.  Check 3 already
                 proves that a homopolymer-suppressed substitution emits the reference base; this
                 catches the opposite failure, a rebuild where the HPLEN annotation went missing
                 and the policy silently did nothing.

Exit code is non-zero if any check fails.

    python 5_validate_variant_genome.py --genome_ver GRCh38 \
        --out_genome_ver GRCh38_KOLF2.1J_hap1 --chroms 21,22
"""

import argparse
import json
import os
import pickle
import sys


def load(path):
    with open(path, "rb") as fh:
        return pickle.load(fh)


def check_contig(genome_dir, base, new, chrom, linked=False):
    ref_rec = load(os.path.join(genome_dir, "fa_pickle", base, f"{chrom}.pk"))
    new_rec = load(os.path.join(genome_dir, "fa_pickle", new, f"{chrom}.pk"))
    vs = os.path.join(genome_dir, "variant_sets", new)
    applied = load(os.path.join(vs, "applied", f"{chrom}.pk"))
    skipped = load(os.path.join(vs, "skipped", f"{chrom}.pk"))

    rs, ns = str(ref_rec.seq), str(new_rec.seq)
    problems = []

    # 1. lossless
    if len(rs) != len(ns):
        problems.append(f"length changed: {len(rs)} -> {len(ns)}")
    if new_rec.id != ref_rec.id or new_rec.name != ref_rec.name:
        problems.append(f"id/name changed: {ref_rec.id}/{ref_rec.name} -> {new_rec.id}/{new_rec.name}")
    if linked:
        # nothing was applied, so the contig is hardlinked from the reference; assert it really
        # is identical rather than expecting a provenance string it cannot have
        if str(new_rec.seq) != rs:
            problems.append("contig marked unchanged but differs from the reference")
    elif "variant-substituted" not in new_rec.description:
        problems.append("provenance missing from description")

    # 2. applied
    n_app = len(applied["pos"])
    bad_ref = bad_alt = bad_case = 0
    changed = set()
    for pos, ref, alt in zip(applied["pos"], applied["ref"], applied["alt"]):
        i = pos - 1
        if rs[i : i + len(ref)].upper() != ref.upper():
            bad_ref += 1
        if ns[i : i + len(alt)].upper() != alt.upper():
            bad_alt += 1
        for k in range(len(alt)):
            if rs[i + k].islower() != ns[i + k].islower():
                bad_case += 1
            if rs[i + k] != ns[i + k]:
                changed.add(i + k)
    if bad_ref:
        problems.append(f"{bad_ref}/{n_app} applied variants did not match REF beforehand")
    if bad_alt:
        problems.append(f"{bad_alt}/{n_app} applied variants do not read ALT afterwards")
    if bad_case:
        problems.append(f"{bad_case}/{n_app} applied variants lost soft-mask case")

    # 3. minimal - instead of scanning 3 Gb char-by-char in Python, revert the known applied
    # positions in the substituted sequence and assert the result is byte-identical to the
    # reference.  One C-level comparison rather than a per-base loop.
    n_diffs = len(changed)
    if len(rs) == len(ns):
        rb = rs.encode()
        reverted = bytearray(ns.encode())
        for pos, ref, alt in zip(applied["pos"], applied["ref"], applied["alt"]):
            i = pos - 1
            reverted[i : i + len(alt)] = rb[i : i + len(alt)]
        if bytes(reverted) != rb:
            # locate the first offending block for a useful message
            first = next(
                (i for i in range(len(rb)) if reverted[i] != rb[i]), None
            )
            problems.append(
                f"differences outside the applied variants, first at offset {first}"
            )

    # 4. skipped, allowing legitimate overlap with a different applied variant
    applied_offsets = changed
    touched_illegally = 0
    for pos, ref in zip(skipped["pos"], skipped["ref"]):
        span = range(pos - 1, pos - 1 + len(ref))
        for i in span:
            if i < len(rs) and rs[i] != ns[i] and i not in applied_offsets:
                touched_illegally += 1
                break
    if touched_illegally:
        problems.append(f"{touched_illegally} skipped variants altered outside any applied variant")

    overlapped = sum(
        1
        for pos, ref in zip(skipped["pos"], skipped["ref"])
        if any(i in applied_offsets for i in range(pos - 1, pos - 1 + len(ref)))
    )

    return {
        "chrom": chrom,
        "length": len(ns),
        "applied": n_app,
        "skipped": len(skipped["pos"]),
        "skipped_overlapping_applied": overlapped,
        "diff_positions": n_diffs,
        "problems": problems,
    }


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--genome_dir", default="genome_files")
    p.add_argument("--genome_ver", default="GRCh38", help="base/reference genome")
    p.add_argument("--out_genome_ver", required=True, help="substituted genome to validate")
    p.add_argument("--chroms", help="comma-separated subset; default = all in the manifest")
    a = p.parse_args(argv)

    vs = os.path.join(a.genome_dir, "variant_sets", a.out_genome_ver)
    with open(os.path.join(vs, "manifest.json")) as fh:
        manifest = json.load(fh)
    chroms = a.chroms.split(",") if a.chroms else manifest["contigs"]

    print(f"validating {a.out_genome_ver} against {a.genome_ver}")
    print(f"  sample {manifest['sample']}  haplotype {manifest['haplotype']}\n")
    print(f"{'contig':>8} {'length':>13} {'applied':>9} {'skipped':>9} {'overlap':>8} {'diffs':>8}  status")

    failed = 0
    for chrom in chroms:
        linked = manifest["per_contig"].get(chrom, {}).get("linked_unchanged", False) or (
            chrom in manifest.get("contigs_linked_unchanged", [])
        )
        r = check_contig(a.genome_dir, a.genome_ver, a.out_genome_ver, chrom, linked=linked)
        status = ("ok" if not r["problems"] else "FAIL") + (" (linked)" if linked else "")
        failed += bool(r["problems"])
        print(
            f"{r['chrom']:>8} {r['length']:>13,} {r['applied']:>9,} {r['skipped']:>9,} "
            f"{r['skipped_overlapping_applied']:>8,} {r['diff_positions']:>8,}  {status}"
        )
        for msg in r["problems"]:
            print(f"           - {msg}")

    # 5. the policies the manifest claims were in force must have left a trace
    policy_problems = []
    if not a.chroms:
        by_reason = {}
        for v in manifest["per_contig"].values():
            for reason, n in (v.get("skipped_by_reason") or {}).items():
                by_reason[reason] = by_reason.get(reason, 0) + n
        hp_min = manifest.get("homopolymer_ref_min_len", 0)
        if hp_min:
            n_hp = by_reason.get("homopolymer_snv", 0)
            print(f"\nhomopolymer policy: runs >= {hp_min} bp, "
                  f"{n_hp:,} substitutions suppressed to reference")
            if not n_hp:
                policy_problems.append(
                    f"manifest says --homopolymer_ref_min_len {hp_min} but not one substitution "
                    "was suppressed; the source VCF was probably rebuilt without its HPLEN tags"
                )
        if manifest.get("nocall_bed") and not manifest.get("nocall_regions"):
            policy_problems.append(
                "manifest names a --nocall_bed but recorded zero no-call regions"
            )
        for msg in policy_problems:
            print(f"  - {msg}")

    detail = []
    if failed:
        detail.append(f"{failed} contig(s)")
    if policy_problems:
        detail.append(f"{len(policy_problems)} policy check(s)")
    print("\nRESULT:", "PASS" if not detail else "FAIL (" + ", ".join(detail) + ")")
    return 1 if detail else 0


if __name__ == "__main__":
    sys.exit(main())
