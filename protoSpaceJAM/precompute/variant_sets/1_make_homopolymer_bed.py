#!/usr/bin/env python
"""
Homopolymer runs in a reference FASTA -> BED on stdout.

Step 1 of the variant-set pipeline (see README.md).  The output feeds `bcftools annotate` in step 2,
which stamps an HPLEN/HPBASE INFO tag onto every variant sitting in a run.  The builder in step 4
reads that tag off the record; it has no BED-intersection capability and a 148 MB BED would be
unusable from a stdlib-only code path.

Columns: chrom, start (0-based), end (exclusive), run length, run base.

Reference-only, so it is built once and reused for every sample called against that genome.

    python 1_make_homopolymer_bed.py --fasta GRCh38.fa --min_len 5 | bgzip > hp_ge5.bed.gz
    tabix -p bed hp_ge5.bed.gz

~2 min and ~500 MB peak RSS on GRCh38 (chr1 dominates both).
"""

import argparse
import re
import sys


def iter_contigs(path):
    """Stream (name, sequence) out of a FASTA without holding more than one contig."""
    name, chunks = None, []
    with open(path, "r") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                # the id is everything up to the first whitespace, matching samtools faidx
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def main(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--fasta", required=True, help="reference FASTA (may be soft-masked)")
    p.add_argument("--min_len", type=int, default=5, help="shortest run to emit (default: 5)")
    p.add_argument("--out", help="output file; default stdout, so you can pipe to bgzip")
    a = p.parse_args(argv)

    # The FASTA is dna_sm: soft-masked bases are lower case and a run can straddle the boundary
    # between masked and unmasked sequence, so match on the upper-cased sequence.  N runs are not
    # homopolymers in the sense that matters here (they are assembly gaps), so they are excluded.
    pattern = re.compile("|".join(f"{b}{{{a.min_len},}}" for b in "ACGT"))

    out = open(a.out, "w") if a.out else sys.stdout
    runs = total_bp = 0
    try:
        for name, seq in iter_contigs(a.fasta):
            seq = seq.upper()
            n_here = 0
            for m in pattern.finditer(seq):
                start, end = m.start(), m.end()
                out.write(f"{name}\t{start}\t{end}\t{end - start}\t{seq[start]}\n")
                n_here += 1
                total_bp += end - start
            runs += n_here
            print(f"  {name:>14}  {len(seq):>12,} bp  {n_here:>9,} runs", file=sys.stderr)
    finally:
        if a.out:
            out.close()

    print(f"\n{runs:,} runs >= {a.min_len} bp, {total_bp:,} bp total", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
