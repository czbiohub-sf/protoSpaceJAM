release=release-107
release_num=107

########
#GRCh38#
########
#fasta
genome_file_dir_remote="http://ftp.ensembl.org/pub/${release}/fasta/homo_sapiens/dna/"
genome_file="Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz"
rm ${genome_file} 
wget ${genome_file_dir_remote}${genome_file} 
#gff3
genome_file_dir_remote="http://ftp.ensembl.org/pub/${release}/gff3/homo_sapiens/"
genome_file="Homo_sapiens.GRCh38.${release_num}.gff3.gz"
rm ${genome_file} 
wget ${genome_file_dir_remote}${genome_file} 

########
#GRCm39#
########
genome_file_dir_remote="http://ftp.ensembl.org/pub/${release}/fasta/mus_musculus/dna/"
genome_file="Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz"
rm ${genome_file} 
wget ${genome_file_dir_remote}${genome_file} 
#gff3
genome_file_dir_remote="http://ftp.ensembl.org/pub/${release}/gff3/mus_musculus/"
genome_file="Mus_musculus.GRCm39.${release_num}.gff3.gz"
rm ${genome_file} 
wget ${genome_file_dir_remote}${genome_file} 

########
#GRCz11#
########
genome_file_dir_remote="http://ftp.ensembl.org/pub/${release}/fasta/danio_rerio/dna/"
genome_file="Danio_rerio.GRCz11.dna_sm.primary_assembly.fa.gz"
rm ${genome_file} 
wget ${genome_file_dir_remote}${genome_file} 
#gff3
genome_file_dir_remote="http://ftp.ensembl.org/pub/${release}/gff3/danio_rerio/"
genome_file="Danio_rerio.GRCz11.${release_num}.gff3.gz"
rm ${genome_file} 
wget ${genome_file_dir_remote}${genome_file} 

echo "Download done"