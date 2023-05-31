module load anaconda # hpc
conda activate protospacejam

##################
#download genomes#
##################
bash download_genomes.sh # download genomes from ensembl, output directory: genomes
bash serialize_fa.sh # serialize genome fasta for fast access

####################################
#extract gene models and codon info#
####################################
python scripts/extract_gene_models_info.py --gff3_gz genome_files/Homo_sapiens.GRCh38.109.gff3.gz --out_dir genome_files/parsed_gff3/GRCh38
python scripts/extract_gene_models_info.py --gff3_gz genome_files/Mus_musculus.GRCm39.109.gff3.gz --out_dir genome_files/parsed_gff3/GRCm39
python scripts/extract_gene_models_info.py --gff3_gz genome_files/Danio_rerio.GRCz11.109.gff3.gz --out_dir genome_files/parsed_gff3/GRCz11

###############
#scan for gRNA#
###############
mkdir gRNAs
bash search_gRNAs.sh GRCh38 Homo_sapiens
bash search_gRNAs.sh GRCm39 Mus_musculus
bash search_gRNAs.sh GRCz11 Danio_rerio

###############
#build indexes#
###############
bash build_bwa_index.sh GRCh38 Homo_sapiens
bash build_bwa_index.sh GRCm39 Mus_musculus
bash build_bwa_index.sh GRCz11 Danio_rerio

#######
# map #
#######
#split 
bash split4mapping.sh GRCh38 Homo_sapiens
bash split4mapping.sh GRCm39 Mus_musculus
bash split4mapping.sh GRCz11 Danio_rerio
#bwa
rm -rf bwa.slurm.out && mkdir bwa.slurm.out
chmod -R a+xX /hpc/projects/data_lg/duo.peng/pJAM_precompute/utils/FindOfftargetBwa/bin
sbatch map_gRNA.sh GRCh38 Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa .
sbatch map_gRNA.sh GRCm39 Mus_musculus.GRCm39.dna_sm.primary_assembly.fa .
sbatch map_gRNA.sh GRCz11 Danio_rerio.GRCz11.dna_sm.primary_assembly.fa .


#########
# score #
#########
sbatch score_gRNA.sh GRCh38 Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa .
sbatch score_gRNA.sh GRCm39 Mus_musculus.GRCm39.dna_sm.primary_assembly.fa .
sbatch score_gRNA.sh GRCz11 Danio_rerio.GRCz11.dna_sm.primary_assembly.fa .
#index
python scripts/index_gRNA_results.py --gzdir gRNAs/gRNA_GRCh38/gRNA.tab.gz.split.BwaMapped.scored


#####################################################
#clean up and move precomputed files to destination #
#####################################################
bash post_processing.sh
