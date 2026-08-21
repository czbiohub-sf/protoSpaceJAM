#!/usr/bin/env bash
# Step 3 of the variant-set pipeline (see README.md).
#
# BED of the positions where one sample has NO genotype call.
#
# This exists because of an asymmetry that is easy to miss.  Step 2 extracts the sample with
# `bcftools view --exclude-uncalled`, which deletes every `./.` record -- so those sites are simply
# absent from the filtered VCF, indistinguishable from sites where the sample matches the
# reference.  They are not the same thing: one is "we looked and it is reference", the other is
# "we do not know".  Without this file, step 4 emits reference there and reports the window
# reliable, which is a confident statement about data that does not exist.
#
# Must be run against the ORIGINAL (pre-extraction) VCF -- the information is gone from the
# filtered one by construction.
#
#   ./3_make_nocall_bed.sh --vcf joint.vcf.gz --sample czML1175_KOLF2.1J-Parental \
#       --out kolf.nocall.bed.gz
#
# Requires bcftools, bgzip and tabix on PATH.

set -euo pipefail

VCF= SAMPLE= OUT=
while [[ $# -gt 0 ]]; do
    case "$1" in
        --vcf)    VCF=$2;    shift 2 ;;
        --sample) SAMPLE=$2; shift 2 ;;
        --out)    OUT=$2;    shift 2 ;;
        -h|--help) sed -n '2,20p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "unknown argument: $1" >&2; exit 1 ;;
    esac
done
[[ -n $VCF && -n $SAMPLE && -n $OUT ]] || {
    echo "usage: $0 --vcf FILE --sample NAME --out FILE.bed.gz" >&2; exit 1; }

# The span is length(REF), so a no-call at a deletion record covers the whole deleted stretch.
# Contig names are normalized the same way step 2 does it, so the two agree.
bcftools query -s "$SAMPLE" -f '%CHROM\t%POS\t%REF\t[%GT]\n' "$VCF" \
| awk -F'\t' '$4=="./." || $4==".|." || $4=="." {
    c=$1; sub(/^chr/,"",c); if (c=="M") c="MT";
    print c "\t" ($2-1) "\t" ($2-1+length($3))
  }' \
| sort -k1,1 -k2,2n \
| bgzip > "$OUT"
tabix -p bed -f "$OUT"

n=$(bgzip -dc "$OUT" | wc -l | tr -d ' ')
bp=$(bgzip -dc "$OUT" | awk '{s+=$3-$2} END{print s+0}')
echo "wrote $OUT"
echo "  $n intervals, $bp bp"
echo "  (step 4 merges book-ended and overlapping intervals when it builds the nocall/ shard)"

if [[ $n -eq 0 ]]; then
    cat >&2 <<'WARN'

WARNING: zero no-call intervals.  A real WGS callset always has some, so the likely cause is that
         --vcf already had its uncalled sites stripped -- i.e. you pointed this at the extracted
         single-sample file rather than at the joint callset it came from.  Point it at the
         original.  An empty BED here is worse than no BED: step 4 will record a --nocall_bed in
         the manifest and 5_validate_variant_genome.py will fail the asset for having none.
WARN
    exit 1
fi
