#!/bin/bash

#SBATCH --job-name=gRNA_scan
#SBATCH --time=14-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=128G
#SBATCH --cpus-per-task=2
#SBATCH -e slurm_gRNA_scan_%j.err
#SBATCH -o slurm_gRNA_scan-%j.out

module load anaconda
conda activate protospaceX

working_dir=/hpc/projects/data_lg/duo.peng/protospaceXS/precomputed_gRNAs
cd ${working_dir}

########
#GRCh38#
########
#prepare output dir
outdir=gRNA_GRCh38
rm -rf ${outdir}
mkdir ${outdir}
#get all gRNAs
script_folder="../scripts"
genome_file_dir="../genome_files"
genome_file="Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz"
python ${script_folder}/gRNA_scan_get_all_gRNA.py --fastagz ${genome_file_dir}/${genome_file} --outdir ${outdir}

########
#GRCm39#
########
#prepare output dir
outdir=gRNA_GRCm39
rm -rf ${outdir}
mkdir ${outdir}
#get all gRNAs
script_folder="../scripts"
genome_file_dir="../genome_files"
genome_file="Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz"
python ${script_folder}/gRNA_scan_get_all_gRNA.py --fastagz ${genome_file_dir}/${genome_file} --outdir ${outdir}

########
#GRCz11#
########
#prepare output dir
outdir=gRNA_GRCz11
rm -rf ${outdir}
mkdir ${outdir}
#get all gRNAs
script_folder="../scripts"
genome_file_dir="../genome_files"
genome_file="Danio_rerio.GRCz11.dna_sm.primary_assembly.fa.gz"
python ${script_folder}/gRNA_scan_get_all_gRNA.py --fastagz ${genome_file_dir}/${genome_file} --outdir ${outdir}