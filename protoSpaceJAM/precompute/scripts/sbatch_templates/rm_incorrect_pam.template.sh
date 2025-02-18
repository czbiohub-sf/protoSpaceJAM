#!/bin/bash

#SBATCH --job-name=rm_incorrect_pam
#SBATCH --time=2-00:00:00
#SBATCH --array=1-placeholder1%placeholder2
#SBATCH --nodes=1
#SBATCH --partition preempted
#SBATCH --ntasks=1
#SBATCH --mem=12G
#SBATCH --cpus-per-task=1
#SBATCH -e rm_incorrect_pam.slurm.out/slurm-%A_%a.err
#SBATCH -o rm_incorrect_pam.slurm.out/slurm-%A_%a.out

module load anaconda
conda activate sklearn0 # Python 3.5 env used to run score predictions made with older version of sklearn


#setting variables
declare -x idx=$(( ${SLURM_ARRAY_TASK_ID} -1))


#setting directories
working_dir=${3}
gzdir="gRNAs/gRNA_${1}/gRNA.tab.gz.split.BwaMapped"
outdir="gRNAs/gRNA_${1}/gRNA.tab.gz.split.BwaMapped_rm_incorrect_PAM"
tabfile="gRNAs/gRNA_${1}/${2}.gz.out.split.tab"
script_folder="./scripts"
pam="${4}"
cd $working_dir


# read file array from tab file
readarray -t filenames < <(cat $tabfile | cut -f6)
filenames=("${filenames[@]:1}") #skip header
echo "number of files: ${#filenames[@]}"

echo "current index is: ${idx}"
echo "current file name is: ${filenames[${idx}]}"

#derive the name of the mapped gzfile 
gzfile=${filenames[$idx]} 
mappedgzfile="${gzfile::-7}.tab.mapped.gz"


# main command
python ${script_folder}/rm_incorrect_PAM.py --pam $pam --infile "${gzdir}/${mappedgzfile}" --outfile "${outdir}/${mappedgzfile}"







