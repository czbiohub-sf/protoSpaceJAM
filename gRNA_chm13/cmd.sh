script_folder="../scripts"
genome_file_dir_remote="https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/"
genome_file="chm13v2.0.fa.gz"
genome_file_dir="../genome_files"

wget ${genome_file_dir_remote}${genome_file} --directory-prefix ${genome_file_dir}
python ${script_folder}/1.0.get_all_human_gRNA.py --fastagz ${genome_file_dir}/${genome_file}


python circos.get.gRNAnum.py --dir gRNA.tab.gz --tab chm13v2.0.fa.gz.out.tab
python 1.1.split4mapping.py --dir gRNA.tab.gz --tab chm13v2.0.fa.gz.out.tab