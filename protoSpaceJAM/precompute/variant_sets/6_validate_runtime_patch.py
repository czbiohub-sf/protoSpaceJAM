#!/usr/bin/env python
"""
Oracle-equivalence test for the runtime variant layer.

The materialized genome written by 4_make_variant_genome.py is the oracle: substitution happened
once, offline, and was verified exhaustively by 5_validate_variant_genome.py.  This script asserts
that the runtime path produces byte-identical sequence:

    get_seq(reference pickle) |> VariantSet.patch()      ==      get_seq(substituted pickle)

for randomly drawn windows, on both strands, deliberately including windows that straddle applied
variants and contig edges.  Any off-by-one in patch()'s `pos - start` arithmetic fails here rather
than silently in a donor.

    python 6_validate_runtime_patch.py --variant_set GRCh38_KOLF2.1J_hap1 --windows 200000
"""

import argparse
import os
import pickle
import random
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..")))

from protoSpaceJAM.util.variant_annot import VariantSet  # noqa: E402

COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def revcom(s):
    return s.translate(COMPLEMENT)[::-1]


def load_seq(genome_dir, genome_ver, chrom):
    with open(os.path.join(genome_dir, "fa_pickle", genome_ver, f"{chrom}.pk"), "rb") as fh:
        return str(pickle.load(fh).seq)


def draw_windows(rng, seq_len, applied_pos, n):
    """
    Mix of window kinds:
      - uniformly random
      - centred on an applied variant (the interesting case)
      - straddling a variant at every offset within the window
      - clamped to both contig edges
    """
    out = []
    edge = [1, 2, max(1, seq_len - 3000)]
    for e in edge:
        for length in (20, 23, 500, 2000):
            out.append((e, length))
    for _ in range(n // 2):
        length = rng.choice([20, 23, 50, 200, 500, 1000, 2000])
        start = rng.randint(1, max(1, seq_len - length))
        out.append((start, length))
    if applied_pos:
        for _ in range(n - len(out)):
            v = applied_pos[rng.randrange(len(applied_pos))]
            length = rng.choice([20, 23, 50, 200, 500, 1000, 2000])
            # offset so the variant lands anywhere in the window, including both edges
            off = rng.randint(0, length - 1)
            start = max(1, min(v - off, max(1, seq_len - length)))
            out.append((start, length))
    return out[:max(n, len(edge) * 4)]


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--genome_dir", default="genome_files")
    p.add_argument("--variant_set", required=True)
    p.add_argument("--chroms", help="comma-separated; default = contigs with applied variants")
    p.add_argument("--windows", type=int, default=100000, help="windows per contig")
    p.add_argument("--seed", type=int, default=0)
    a = p.parse_args(argv)

    vset = VariantSet(a.variant_set, genome_dir=a.genome_dir)
    base = vset.base_genome_ver
    chroms = a.chroms.split(",") if a.chroms else [
        c for c, v in vset.manifest["per_contig"].items() if v["applied"]
    ]

    print(f"oracle : fa_pickle/{a.variant_set}")
    print(f"runtime: fa_pickle/{base} + variant_sets/{a.variant_set}")
    print(f"{'contig':>8} {'windows':>9} {'w/variant':>10} {'mismatch':>9}  status")

    rng = random.Random(a.seed)
    total_w = total_v = total_bad = 0
    failures = []

    for chrom in chroms:
        ref = load_seq(a.genome_dir, base, chrom)
        oracle = load_seq(a.genome_dir, a.variant_set, chrom)
        shard = vset._load(chrom)["applied"]
        applied_pos = list(shard["pos"]) if shard else []

        windows = draw_windows(rng, len(ref), applied_pos, a.windows)
        bad = with_var = 0
        for start, length in windows:
            end = start + length
            if end - 1 > len(ref):
                continue
            plain = ref[start - 1 : end - 1]
            expect = oracle[start - 1 : end - 1]
            got = vset.patch(chrom, start, plain)
            if plain != expect:
                with_var += 1
            for strand_got, strand_expect in (
                (got, expect),
                (revcom(got), revcom(expect)),
            ):
                if strand_got != strand_expect:
                    bad += 1
                    if len(failures) < 5:
                        failures.append((chrom, start, length, strand_expect[:60], strand_got[:60]))
                    break
        status = "ok" if bad == 0 else "FAIL"
        print(f"{chrom:>8} {len(windows):>9,} {with_var:>10,} {bad:>9,}  {status}")
        total_w += len(windows)
        total_v += with_var
        total_bad += bad

    print(f"\n{'TOTAL':>8} {total_w:>9,} {total_v:>10,} {total_bad:>9,}")
    for chrom, start, length, exp, got in failures:
        print(f"  {chrom}:{start}+{length}\n    expect {exp}\n    got    {got}")
    print("\nRESULT:", "PASS" if total_bad == 0 else f"FAIL ({total_bad} mismatches)")
    return 1 if total_bad else 0


if __name__ == "__main__":
    sys.exit(main())
