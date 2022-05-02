data_lg="/hpc/projects/data_lg"
working_dir="${data_lg}/duo.peng/protospaceXS/gRNA_hg38"
script_folder="../scripts"
genome_file_dir_remote="http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/"
genome_file="Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz"
genome_file_dir="../genome_files"
offtarget_script_dir="/home/duo.peng/github_repos/protospaceXS/FindOfftargetBwa"
cd $working_dir

wget ${genome_file_dir_remote}${genome_file} --directory-prefix ${genome_file_dir}
python ${script_folder}/1.0.get_all_human_gRNA.py --fastagz ${genome_file_dir}/${genome_file}
python ${script_folder}/circos.get.gRNAnum.py --dir gRNA.tab.gz --tab ${genome_file_dir}/${genome_file}.out.tab
python ${script_folder}/1.1.split4mapping.py --dir gRNA.tab.gz --tab ${genome_file_dir}/${genome_file}.out.tab --part_size 200000

#parallelize: (add thread information)
python ${script_folder}/1.2.getOfftarget.py --gzfile MT.tab.gz --gzdir gRNA.tab.gz.split --gRNA_count 2193 --genome_fa Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa --thread 2
python ${script_folder}/1.2.getOfftarget.py --gzfile 1.part0.tab.gz --gzdir gRNA.tab.gz.split --gRNA_count 200000 --genome_fa Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa --thread 2
#python ${script_folder}/1.2.getOfftarget.py --gzdir gRNA.tab.gz.split --gzfile ${filename[$idx]}  --gRNA_count ${gRNA_count[$idx]} --thread ${SLURM_CPUS_PER_TASK} --genome_fa Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa

#batch job,  need to first throttle to 1 (then increase throttle), this avoids simultaneous creating of output dir and cause some jobs to fail
sbatch array.job.bwa.sh

#check if all gRNAs are mapped by bwa
sbatch array.job.bwa_check.sh

#calculate scores. TODO: make it cluster array job
#conda activate protospacerX
#python ${script_folder}/1.3.calcSpecificityScores.py --gzdir gRNA.tab.gz.split.BwaMapped --gzfile 1.part0.tab.mapped.gz

#calculate scores. 
module load anaconda
conda activate sklearn0
python ${script_folder}/1.3.calcScores.py --gzdir gRNA.tab.gz.split.BwaMapped --gzfile MT.part0.tab.mapped.gz
