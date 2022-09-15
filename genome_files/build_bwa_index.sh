#!/bin/bash

#SBATCH --job-name=Bwa_idx
#SBATCH --time=14-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH -e bwa_idx-%j.err
#SBATCH -o bwa_idx-%j.out

release=release-107
release_num=107

module load anaconda
conda activate protospaceX

########
#GRCh38#
########
genome_file="Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz"
cp ${genome_file} indexes_bwa/${genome_file_fagz}
cd indexes_bwa
gunzip ${genome_file}
chmod a+xX ../../utils/FindOfftargetBwa/bin/Linux/bwa
../../utils/FindOfftargetBwa/bin/Linux/bwa index ${genome_file%.gz}
cd ..
python get_fagz_sizes.py --fastagz ${genome_file}
mv ${genome_file}.sizes indexes_bwa/${genome_file%.gz}.sizes

########
#GRCm39#
########
genome_file="Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz"
cp ${genome_file} indexes_bwa/${genome_file_fagz}
cd indexes_bwa
gunzip ${genome_file}
chmod a+xX ../../utils/FindOfftargetBwa/bin/Linux/bwa
../../utils/FindOfftargetBwa/bin/Linux/bwa index ${genome_file%.gz}
cd ..
python get_fagz_sizes.py --fastagz ${genome_file}
mv ${genome_file}.sizes indexes_bwa/${genome_file%.gz}.sizes

########
#GRCz11#
########
genome_file="Danio_rerio.GRCz11.dna_sm.primary_assembly.fa.gz"
cp ${genome_file} indexes_bwa/${genome_file_fagz}
cd indexes_bwa
gunzip ${genome_file}
chmod a+xX ../../utils/FindOfftargetBwa/bin/Linux/bwa
../../utils/FindOfftargetBwa/bin/Linux/bwa index ${genome_file%.gz}
cd ..
python get_fagz_sizes.py --fastagz ${genome_file}
mv ${genome_file}.sizes indexes_bwa/${genome_file%.gz}.sizes