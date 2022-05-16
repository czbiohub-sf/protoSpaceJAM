#!/bin/bash

#SBATCH --job-name=CFD
#SBATCH --time=14-00:00:00
#SBATCH --array=1-76%1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -e slurm.out/slurm-%j.err
#SBATCH -o slurm.out/slurm-%j.out

# declare arrays
declare -a filename=( 
GL000008.2.tab.mapped.gz \
GL000213.1.tab.mapped.gz \
GL000218.1.tab.mapped.gz \
GL000219.1.tab.mapped.gz \
GL000220.1.tab.mapped.gz \
GL000225.1.tab.mapped.gz \
KI270305.1.tab.mapped.gz \
KI270311.1.tab.mapped.gz \
KI270316.1.tab.mapped.gz \
KI270330.1.tab.mapped.gz \
KI270334.1.tab.mapped.gz \
KI270362.1.tab.mapped.gz \
KI270375.1.tab.mapped.gz \
KI270378.1.tab.mapped.gz \
KI270379.1.tab.mapped.gz \
KI270381.1.tab.mapped.gz \
KI270382.1.tab.mapped.gz \
KI270384.1.tab.mapped.gz \
KI270391.1.tab.mapped.gz \
KI270393.1.tab.mapped.gz \
KI270394.1.tab.mapped.gz \
KI270395.1.tab.mapped.gz \
KI270419.1.tab.mapped.gz \
KI270420.1.tab.mapped.gz \
KI270422.1.tab.mapped.gz \
KI270425.1.tab.mapped.gz \
KI270435.1.tab.mapped.gz \
KI270438.1.tab.mapped.gz \
KI270442.1.tab.mapped.gz \
KI270467.1.tab.mapped.gz \
KI270468.1.tab.mapped.gz \
KI270508.1.tab.mapped.gz \
KI270512.1.tab.mapped.gz \
KI270519.1.tab.mapped.gz \
KI270528.1.tab.mapped.gz \
KI270539.1.tab.mapped.gz \
KI270548.1.tab.mapped.gz \
KI270579.1.tab.mapped.gz \
KI270580.1.tab.mapped.gz \
KI270582.1.tab.mapped.gz \
KI270583.1.tab.mapped.gz \
KI270584.1.tab.mapped.gz \
KI270591.1.tab.mapped.gz \
KI270708.1.tab.mapped.gz \
KI270712.1.tab.mapped.gz \
KI270714.1.tab.mapped.gz \
KI270715.1.tab.mapped.gz \
KI270716.1.tab.mapped.gz \
KI270719.1.tab.mapped.gz \
KI270720.1.tab.mapped.gz \
KI270722.1.tab.mapped.gz \
KI270723.1.tab.mapped.gz \
KI270724.1.tab.mapped.gz \
KI270725.1.tab.mapped.gz \
KI270726.1.tab.mapped.gz \
KI270727.1.tab.mapped.gz \
KI270728.1.tab.mapped.gz \
KI270729.1.tab.mapped.gz \
KI270730.1.tab.mapped.gz \
KI270731.1.tab.mapped.gz \
KI270732.1.tab.mapped.gz \
KI270733.1.tab.mapped.gz \
KI270734.1.tab.mapped.gz \
KI270735.1.tab.mapped.gz \
KI270736.1.tab.mapped.gz \
KI270737.1.tab.mapped.gz \
KI270739.1.tab.mapped.gz \
KI270740.1.tab.mapped.gz \
KI270742.1.tab.mapped.gz \
KI270745.1.tab.mapped.gz \
KI270749.1.tab.mapped.gz \
KI270751.1.tab.mapped.gz \
KI270755.1.tab.mapped.gz \
KI270756.1.tab.mapped.gz \
MT.tab.mapped.gz \
Y.part9.tab.mapped.gz \
)


declare -x idx=$(( ${SLURM_ARRAY_TASK_ID} -1))

module load anaconda
conda activate sklearn0

#check index and filename
echo "index -> filename/gRNA_count"
echo "$idx -> ${filename[$idx]}/${gRNA_count[$idx]}"

#setting directories
data_lg="/hpc/projects/data_lg"
working_dir="${data_lg}/duo.peng/protospaceXS/gRNA_hg38"
script_folder="../scripts"
cd $working_dir

#main command
python ${script_folder}/1.3.calcScores.py --gzdir gRNA.tab.gz.split.BwaMapped --gzfile ${filename[$idx]}

