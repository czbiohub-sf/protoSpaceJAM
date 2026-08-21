#!/usr/bin/env bash
# Step 2 of the variant-set pipeline (see README.md).
#
# Joint or single-sample VCF  ->  one sample, Ensembl contig names, normalized, quality-filtered,
# multiallelics split, homopolymer context annotated.  The output is what step 4 consumes.
#
# Every threshold here is argued for in README.md "The filter, and why each threshold".  Change
# them by flag, not by editing, so the command line stays a record of what was built.
#
#   ./2_filter_callset.sh --vcf joint.vcf.gz --sample czML1175_KOLF2.1J-Parental \
#       --fasta GRCh38.fa --hp_bed hp_ge5.bed.gz --out kolf.filtered.vcf.gz
#
# Requires bcftools >= 1.15 and tabix on PATH.  ~6 min on a 5.4 M-record callset.

set -euo pipefail

GQ=20
DP_MIN=25
DP_MAX=120
# GRCh38 X pseudoautosomal boundaries.  Outside these, a male sample cannot be heterozygous on X.
PAR1_END=2781479
PAR2_START=155701383
VCF= SAMPLE= FASTA= HP_BED= OUT=

usage() {
    sed -n '2,15p' "$0" | sed 's/^# \{0,1\}//'
    cat <<USAGE

required:
  --vcf FILE       source VCF (bgzipped, indexed).  May be multi-sample.
  --sample NAME    sample column to extract
  --fasta FILE     reference FASTA the VCF was called against, for left-alignment
  --out FILE       output .vcf.gz

optional:
  --hp_bed FILE    homopolymer BED from step 1.  Strongly recommended: without it the built
                   variant set carries no HPLEN tags and step 4's homopolymer policy silently
                   does nothing (step 5 will fail the asset for it).
  --gq INT         minimum genotype quality (default: $GQ)
  --dp_min INT     minimum depth (default: $DP_MIN)
  --dp_max INT     maximum depth, MT exempt (default: $DP_MAX)
USAGE
    exit "${1:-1}"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --vcf)     VCF=$2;     shift 2 ;;
        --sample)  SAMPLE=$2;  shift 2 ;;
        --fasta)   FASTA=$2;   shift 2 ;;
        --hp_bed)  HP_BED=$2;  shift 2 ;;
        --out)     OUT=$2;     shift 2 ;;
        --gq)      GQ=$2;      shift 2 ;;
        --dp_min)  DP_MIN=$2;  shift 2 ;;
        --dp_max)  DP_MAX=$2;  shift 2 ;;
        -h|--help) usage 0 ;;
        *) echo "unknown argument: $1" >&2; usage ;;
    esac
done
[[ -n $VCF && -n $SAMPLE && -n $FASTA && -n $OUT ]] || usage

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

# ------------------------------------------------------------------ contig name map
# The joint callset is chr-prefixed; protoSpaceJAM's genome pickles use Ensembl names.  Derive the
# map from the two files rather than shipping a static one: a hand-written map silently maps
# contigs the reference does not have (chrEBV -> EBV is the classic), which breaks `norm`.
[[ -f "$FASTA.fai" ]] || samtools faidx "$FASTA"
cut -f1 "$FASTA.fai" | sort -u > "$work/ref_contigs"
bcftools view -h "$VCF" | sed -n 's/^##contig=<ID=\([^,>]*\).*/\1/p' | sort -u > "$work/vcf_contigs"

: > "$work/rename_chrs.txt"
while read -r c; do
    case "$c" in
        chrM)  t=MT ;;
        chr*)  t=${c#chr} ;;
        *)     continue ;;          # already Ensembl-style, nothing to do
    esac
    # only rename onto a name the reference actually has
    if grep -qxF "$t" "$work/ref_contigs"; then printf '%s\t%s\n' "$c" "$t" >> "$work/rename_chrs.txt"; fi
done < "$work/vcf_contigs"
echo "contig renames: $(wc -l < "$work/rename_chrs.txt")"

# ------------------------------------------------------------------ extract, normalize, filter
# Order matters twice over:
#   - rename first, or the CHROM=="MT"/"X"/"Y" clauses below match nothing and look like a filter
#     that works but removes everything
#   - norm before view, so the GQ/DP tests see one ALT per record
# De Morgan the haploid-contig rules: `!(...)` is not valid in a bcftools filter expression.
echo "[1/2] extract $SAMPLE, rename, normalize, filter (GQ>=$GQ, DP $DP_MIN-$DP_MAX)"
bcftools view -s "$SAMPLE" --exclude-uncalled --trim-alt-alleles "$VCF" -Ou \
| bcftools annotate --rename-chrs "$work/rename_chrs.txt" -Ou \
| bcftools norm -f "$FASTA" -m -any -a -Ou \
| bcftools view -i "GT=\"alt\" && GT!=\"mis\" && FILTER!=\"MONOALLELIC\"
    && FMT/GQ>=$GQ && FMT/DP>=$DP_MIN
    && (CHROM==\"MT\" || FMT/DP<=$DP_MAX)
    && (CHROM!=\"Y\" || GT!=\"het\")
    && (CHROM!=\"X\" || POS<=$PAR1_END || POS>=$PAR2_START || GT!=\"het\")" \
  -Oz -o "$work/filtered.vcf.gz"
bcftools index -t "$work/filtered.vcf.gz"

# ------------------------------------------------------------------ homopolymer context
if [[ -n $HP_BED ]]; then
    echo "[2/2] annotate HPLEN/HPBASE from $HP_BED"
    cat > "$work/hp.hdr" <<'HDR'
##INFO=<ID=HPLEN,Number=1,Type=Integer,Description="Length of the homopolymer run (>=5bp) overlapping this position; absent if none">
##INFO=<ID=HPBASE,Number=1,Type=String,Description="Base composing the overlapping homopolymer run">
HDR
    bcftools annotate -a "$HP_BED" -c CHROM,FROM,TO,INFO/HPLEN,INFO/HPBASE \
        -h "$work/hp.hdr" "$work/filtered.vcf.gz" -Oz -o "$OUT"
else
    echo "[2/2] no --hp_bed given; SKIPPING homopolymer annotation" >&2
    cp "$work/filtered.vcf.gz" "$OUT"
fi
bcftools index -t -f "$OUT"

# ------------------------------------------------------------------ report
echo
echo "wrote $OUT"
echo "  records : $(bcftools index -n "$OUT")"
echo "  contigs : $(tabix -l "$OUT" | tr '\n' ' ')"
[[ -n $HP_BED ]] && echo "  HPLEN>=5: $(bcftools view -H -i 'INFO/HPLEN>=5' "$OUT" | wc -l | tr -d ' ')"
bcftools stats "$OUT" | grep -E '^SN|^TSTV' | cut -f3-
echo
echo "Ti/Tv near 2.0 and MT present in the contig list are the two things worth eyeballing."
