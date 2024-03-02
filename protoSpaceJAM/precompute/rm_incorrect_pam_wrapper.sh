#!/bin/bash

# Check if an argument is provided
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <genome-identifier> <genome-fasta-name> <pam> <size-of-array-job> <throttle>"
    exit 1
fi

# prerequisite commands
#rm -rf rm_incorrect_pam.slurm.out && mkdir rm_incorrect_pam.slurm.out
#job_array_size_h38=$(cat gRNAs/gRNA_GRCh38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz.out.split.tab | wc -l)
#job_array_size_m39=$(cat gRNAs/gRNA_GRCm39/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz.out.split.tab | wc -l)
#job_array_size_z11=$(cat gRNAs/gRNA_GRCz11/Danio_rerio.GRCz11.dna_sm.primary_assembly.fa.gz.out.split.tab | wc -l)

# three possible commands
# bash rm_incorrect_pam_wrapper.sh GRCh38 Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa "${pam}" "${job_array_size_h38}" 10
# bash rm_incorrect_pam_wrapper.sh GRCm39 Mus_musculus.GRCm39.dna_sm.primary_assembly.fa "${pam}" "${job_array_size_m39}" 10
# bash rm_incorrect_pam_wrapper.sh GRCz11 Danio_rerio.GRCz11.dna_sm.primary_assembly.fa "${pam}" "${job_array_size_z11}" 10



genome_identifier=$1
genome_fa_name=$2
pam=$3
size_of_array_job=$4
throttle=$5

# Path to the original sbatch script
sbatch_script="scripts/sbatch_templates/rm_incorrect_pam.template.sh"
tmp_sbatch_script="${sbatch_script}.${genome_identifier}.sh"

# make a copy of the original sbatch script
cp "$sbatch_script" "$tmp_sbatch_script"

# Use 'sed' to modify the array range in the sbatch script
sed -i "s/--array=1-placeholder1%placeholder2/--array=1-${size_of_array_job}%${throttle}/" "${tmp_sbatch_script}"
sed -i "s/--job-name=gRNA_score/--job-name=gRNA_score_${genome_identifier}/" "${tmp_sbatch_script}"

# create output dir
outdir="gRNAs/gRNA_${genome_identifier}/gRNA.tab.gz.split.BwaMapped_rm_incorrect_PAM"
mkdir -p $outdir

# submit the sbatch job
sbatch ${tmp_sbatch_script} ${genome_identifier} ${genome_fa_name} . ${pam}

# sleep 1 seconds
sleep 1

# delete the temporary sbatch script
rm "${tmp_sbatch_script}"