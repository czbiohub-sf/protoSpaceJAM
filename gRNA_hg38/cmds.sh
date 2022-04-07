script_folder="../scripts"
genome_file_dir_remote="http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/"
genome_file="Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz"
genome_file_dir="../genome_files"

wget ${genome_file_dir_remote}${genome_file} --directory-prefix ${genome_file_dir}
python ${script_folder}/1.0.get_all_human_gRNA.py --fastagz ${genome_file_dir}/${genome_file}


python circos.get.gRNAnum.py --dir gRNA.tab.gz --tab chm13v2.0.fa.gz.out.tab
python 1.1.split4mapping.py --dir gRNA.tab.gz --tab chm13v2.0.fa.gz.out.tab




#python circos.get.gRNAnum.py --dir gRNA.tab.gz --tab Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz.out.tab
#python 1.1.split4mapping.py --dir gRNA.tab.gz --tab Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz.out.tab