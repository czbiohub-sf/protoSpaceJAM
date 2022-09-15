#!/bin/bash

#SBATCH --job-name=gRNA_split
#SBATCH --time=14-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --cpus-per-task=2
#SBATCH -e slurm_gRNA_split_%j.err
#SBATCH -o slurm_gRNA_split-%j.out

module load anaconda
conda activate protospaceX

working_dir=/hpc/projects/data_lg/duo.peng/protospaceXS/precomuted_gRNAs
cd ${working_dir}

########
#GRCh38#
########
dir=gRNA_GRCh38
genome_file="Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz"
script_folder="../scripts"
python ${script_folder}/gRNA_scan_split4mapping.py --dir ${dir}/gRNA.tab.gz --tab ${dir}/${genome_file}.out.tab --part_size 200000

########
#GRCm39#
########
dir=gRNA_GRCm39
genome_file="Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz"
script_folder="../scripts"
python ${script_folder}/gRNA_scan_split4mapping.py --dir ${dir}/gRNA.tab.gz --tab ${dir}/${genome_file}.out.tab --part_size 200000

########
#GRCz11#
########
dir=gRNA_GRCz11
genome_file="Danio_rerio.GRCz11.dna_sm.primary_assembly.fa.gz"
script_folder="../scripts"
python ${script_folder}/gRNA_scan_split4mapping.py --dir ${dir}/gRNA.tab.gz --tab ${dir}/${genome_file}.out.tab --part_size 200000