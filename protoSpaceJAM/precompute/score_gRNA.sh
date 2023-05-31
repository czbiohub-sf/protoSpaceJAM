#!/bin/bash

#SBATCH --job-name=Score
#SBATCH --time=14-00:00:00
#SBATCH --array=1-1700%500
#SBATCH --partition preempted
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=12G
#SBATCH --cpus-per-task=1
#SBATCH -e score.slurm.out/slurm-%A_%a.err
#SBATCH -o score.slurm.out/slurm-%A_%a.out

module load anaconda
conda activate protospacejam

#setting variables
declare -x idx=$(( ${SLURM_ARRAY_TASK_ID} -1))


#setting directories
working_dir=${3}
gzdir="gRNAs/gRNA_${1}/gRNA.tab.gz.split.BwaMapped"
tabfile="gRNAs/gRNA_${1}/${2}.gz.out.split.tab"
script_folder="./scripts"
cd $working_dir


# read file array from tab file
readarray -t filenames < <(cat $tabfile | cut -f6)
filenames=("${filenames[@]:1}") #skip header
echo "number of files: ${#filenames[@]}"

# read gRNA num from tab file
readarray -t gRNA_counts < <(cat $tabfile | cut -f4)
gRNA_counts=("${gRNA_counts[@]:1}") #skip header
echo "number of gRNA counts: ${#gRNA_counts[@]}"

#derive the name of the mapped gzfile 
gzfile=${filenames[$idx]} 
mappedgzfile="${gzfile::-7}.tab.mapped.gz"

#main command
python ${script_folder}/gRNA_calcScores.py --gzdir ${gzdir} --gzfile ${mappedgzfile} --skip_eff_score # skip efficiency score for now, just calculate the off-target score
