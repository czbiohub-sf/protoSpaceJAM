#!/bin/bash

#SBATCH --job-name=Bwa
#SBATCH --time=14-00:00:00
#SBATCH --array=1-1700%1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH -e map.slurm.out/slurm-%A_%a.err
#SBATCH -o map.slurm.out/slurm-%A_%a.out

module load anaconda
conda activate protospaceX

#setting variables
declare -x idx=$(( ${SLURM_ARRAY_TASK_ID} -1))
working_dir=/hpc/projects/data_lg/duo.peng/protospaceXS/precomuted_gRNAs
script_folder="../scripts"
gzdir="${1}/gRNA.tab.gz.split" # need to change this for other genomes
dirmapped="${1}/gRNA.tab.gz.split.BwaMapped"
tabfile="${1}/${2}.gz.out.split.tab"
genome_fa=$2
MINWAIT=0
MAXWAIT=10

cd $working_dir

# read file array from tab file
readarray -t filenames < <(cat $tabfile | cut -f6)
filenames=("${filenames[@]:1}") #skip header
echo "number of files: ${#filenames[@]}"

# read gRNA num from tab file
readarray -t gRNA_counts < <(cat $tabfile | cut -f4)
gRNA_counts=("${gRNA_counts[@]:1}") #skip header
echo "number of gRNA counts: ${#gRNA_counts[@]}"

#check index and filename
echo "index -> filename/gRNA_count"
echo "$idx -> ${filenames[$idx]}/${gRNA_counts[$idx]}"

#predict the name of the mapped gzfile to check it existence, if found, will skip mapping
gzfile=${filenames[$idx]} 
mappedgzfile="${gzfile::-7}.tab.mapped.gz"

#main command
read linecount < <(gunzip -c ${dirmapped}/${mappedgzfile} | wc -l)
if [ ${linecount} == ${gRNA_counts[$idx]} ] ; then  # mapped gzfile contains all the gRNAs, skip mapping
    sleep $((MINWAIT+RANDOM % (MAXWAIT-MINWAIT)))
    echo "${dirmapped}/${mappedgzfile} has ${linecount} lines, expecting ${gRNA_counts[$idx]}" >> ${1}.mapping_log.skipped.txt
else #do mapping
    sleep $((MINWAIT+RANDOM % (MAXWAIT-MINWAIT)))
    echo "${dirmapped}/${mappedgzfile} has ${linecount} lines, expecting ${gRNA_counts[$idx]}" >> ${1}.mapping_log.mapped.txt
    python ${script_folder}/gRNA_scan_bwa.py --gzdir ${gzdir} --gzfile ${filenames[$idx]} --gRNA_count ${gRNA_counts[$idx]} --thread ${SLURM_CPUS_PER_TASK} --genome_fa ${genome_fa}
fi
