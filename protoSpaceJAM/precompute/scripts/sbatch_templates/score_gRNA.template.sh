#!/bin/bash

#SBATCH --job-name=gRNA_score
#SBATCH --time=14-00:00:00
#SBATCH --array=1-placeholder1%placeholder2
#SBATCH --nodes=1
#SBATCH --partition preempted
#SBATCH --ntasks=1
#SBATCH --mem=12G
#SBATCH --cpus-per-task=1
#SBATCH -e score.slurm.out/slurm-%A_%a.err
#SBATCH -o score.slurm.out/slurm-%A_%a.out

module load anaconda
conda activate sklearn0 # Python 3.5 env used to run score predictions made with older version of sklearn


#setting variables
declare -x idx=$(( ${SLURM_ARRAY_TASK_ID} -1))


#setting directories
working_dir=${3}
gzdir="gRNAs/gRNA_${1}/gRNA.tab.gz.split.BwaMapped"
scoreddir="gRNAs/gRNA_${1}/gRNA.tab.gz.split.BwaMapped.scored"
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

# read gRNA num from tab file
readarray -t gRNA_counts < <(cat $tabfile | cut -f4)
gRNA_counts=("${gRNA_counts[@]:1}") #skip header
echo "number of gRNA in current file: ${gRNA_counts[${idx}]}"

#derive the name of the mapped gzfile 
gzfile=${filenames[$idx]} 
mappedgzfile="${gzfile::-7}.tab.mapped.gz"

#derive the name of the scored gzfile
scoredgzfile="${gzfile::-7}.tab.mapped.scored.gz"

#main command
if [ -e ${scoreddir}/${scoredgzfile} ] && $(gzip -t ${scoreddir}/${scoredgzfile}); then # check if scored gzfile exists
    echo "${scoreddir}/${scoredgzfile} exists and is intact, skipping" >> gRNAs/gRNA_${1}/score_log.skipped.txt
    exit 0
fi
echo "${scoreddir}/${scoredgzfile} does not exist or is not intact, reprocessing" >> gRNAs/gRNA_${1}/score_log.reprocess.txt

#remove previous (incomplete) scoring gzfile
#rm -rf ${scoreddir}/${scoredgzfile}

python ${script_folder}/gRNA_calcScores.py --gzdir ${gzdir} --gzfile ${mappedgzfile} --pam ${pam} --skip_eff_score # skip efficiency score for now, just calculate the off-target score