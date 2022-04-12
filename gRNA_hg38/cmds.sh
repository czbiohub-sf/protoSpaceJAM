working_dir="/home/duo.peng/github_repos/protospaceXS/gRNA_hg38"
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

#parallelize:
python ${script_folder}/1.2.getOfftarget.py --gzfile 1.part0.tab.gz --gzdir gRNA.tab.gz.split \
--gRNA_count 2000000 --genome_fa Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa


