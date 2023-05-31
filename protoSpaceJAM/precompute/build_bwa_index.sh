#!/bin/bash
genome_file_dir="genome_files"
outdir="${genome_file_dir}/indexes_bwa"
# Check if the output directory exists
if [ ! -d ${outdir} ]; then
    mkdir ${outdir}
fi
#make a copy of the genome file
genome_file="${2}.${1}.dna_sm.primary_assembly.fa.gz"
genome_file_path="${genome_file_dir}/${genome_file}"
cp ${genome_file_path} ${outdir}/${genome_file}
cd ${outdir}
# get chromosome sizes
python ../../scripts/get_fagz_sizes.py --fastagz ${genome_file}
mv ${genome_file}.sizes ${genome_file%.gz}.sizes
#bwa
bwa_bin="../../utils/FindOfftargetBwa/bin/Linux/bwa"
chmod a+xX ${bwa_bin}
gunzip ${genome_file}
${bwa_bin} index ${genome_file%.gz}
mv ${genome_file%.gz} ../${genome_file%.gz} # move the genome fasta file to genome_files/, this is needed for the post-bwa processing in the bwa scripts
