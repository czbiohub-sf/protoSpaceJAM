#!/bin/bash

# Check if less than 5 arguments are provided
if [ "$#" -lt 5 ]; then
    echo "Usage: $0 <genome-identifier> <genome-fasta-name> <pam> <size-of-array-job> <throttle> <cpu>"
    exit 1
fi

genome_identifier=$1
genome_fa_name=$2
pam=$3
size_of_array_job=$4
throttle=$5

# Path to the original sbatch script
sbatch_script="scripts/sbatch_templates/score_gRNA.template.sh"
tmp_sbatch_script="${sbatch_script}.${genome_identifier}.sh"

# make a copy of the original sbatch script
cp "$sbatch_script" "$tmp_sbatch_script"

# Use 'sed' to modify the array range in the sbatch script
sed -i "s/--array=1-placeholder1%placeholder2/--array=1-${size_of_array_job}%${throttle}/" "${tmp_sbatch_script}"
sed -i "s/--job-name=gRNA_score/--job-name=gRNA_score_${genome_identifier}/" "${tmp_sbatch_script}"

# if the 6th argument is cpu then replace "preempted" with "cpu"
if [ "$6" == "cpu" ]; then
    sed -i "s/preempted/cpu/" "${tmp_sbatch_script}"
fi

# submit the sbatch job
sbatch ${tmp_sbatch_script} ${genome_identifier} ${genome_fa_name} . ${pam}

# sleep 1 seconds
sleep 1

# delete the temporary sbatch script
rm "${tmp_sbatch_script}"