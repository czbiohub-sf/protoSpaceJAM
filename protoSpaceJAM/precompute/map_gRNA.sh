#!/bin/bash

#SBATCH --job-name=Bwa
#SBATCH --time=14-00:00:00
#SBATCH --array=1-1700%1700
#SBATCH --partition preempted
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH -e bwa.slurm.out/slurm-%A_%a.err
#SBATCH -o bwa.slurm.out/slurm-%A_%a.out

module load anaconda
conda activate protospacejam

#setting variables
declare -x idx=$(( ${SLURM_ARRAY_TASK_ID} -1))
working_dir=${3}
script_folder="./scripts"
gzdir="gRNAs/gRNA_${1}/gRNA.tab.gz.split" 
dirmapped="gRNAs/gRNA_${1}/gRNA.tab.gz.split.BwaMapped"
tabfile="gRNAs/gRNA_${1}/${2}.gz.out.split.tab"
bwa_idx="genome_files/indexes_bwa/$2"
genome_fa="genome_files/$2"
MINWAIT=0
MAXWAIT=5

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
echo "path to mapped gzfile: ${dirmapped}/${mappedgzfile}"

#main command
sleep $((MINWAIT+RANDOM % (MAXWAIT-MINWAIT)))
if [ -e ${dirmapped}/${mappedgzfile} ] && $(gzip -t ${dirmapped}/${mappedgzfile}); then # check if mapped gzfile exists
    read linecount < <(gunzip -c ${dirmapped}/${mappedgzfile} | wc -l)
    if [ ${linecount} == ${gRNA_counts[$idx]} ] ; then  # check if the mapped gzfile contains all the gRNAs, skip mapping if it does
    echo "${dirmapped}/${mappedgzfile} has ${linecount} lines, expecting ${gRNA_counts[$idx]}" >> gRNAs/gRNA_${1}/mapping_log.skipped.txt
    exit 0
    fi
    echo "${dirmapped}/${mappedgzfile} has ${linecount} lines, expecting ${gRNA_counts[$idx]}" >> gRNAs/gRNA_${1}/mapping_log.reprocess.txt
fi
#remove previous (incomplete) mapped gzfile
rm -rf ${dirmapped}/${gzfile::-7}.tab.*
#execute bwa mapping
echo "processing ${mappedgzfile}" >> gRNAs/gRNA_${1}/mapping_log.process.txt
python ${script_folder}/bwa_map.py --gzdir ${gzdir} --gzfile ${filenames[$idx]} --gRNA_count ${gRNA_counts[$idx]}  --genome_fa ${genome_fa} --bwa_idx ${bwa_idx} --thread ${SLURM_CPUS_PER_TASK}

