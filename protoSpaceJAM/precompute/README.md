# Precompute gRNA and genome information
This page provides step-by-step instructions for precomputing gRNAs for a user-specified PAM  
If you prefer downloading the precomputed gRNA, following the instructions [here](https://github.com/czbiohub/protoSpaceJAM#download-and-unzip-pre-computed-data)

## Prerequisite
All scripts and binary files in `protoSpaceJAM/utils` should be excutable

</br>
## Prepare genomes
Download genome sequence and annotation files 
```shell
cd precompute
bash download_genomes.sh
```
Extract gene models and codon info
```shell
module load anaconda
conda activate protospacejam
python scripts/extract_gene_models_info.py --gff3_gz genome_files/Homo_sapiens.GRCh38.109.gff3.gz --out_dir genome_files/parsed_gff3/GRCh38
python scripts/extract_gene_models_info.py --gff3_gz genome_files/Mus_musculus.GRCm39.109.gff3.gz --out_dir genome_files/parsed_gff3/GRCm39
python scripts/extract_gene_models_info.py --gff3_gz genome_files/Danio_rerio.GRCz11.109.gff3.gz --out_dir genome_files/parsed_gff3/GRCz11
bash serialize_fa.sh # serialize genome fasta for fast access
```

</br>
## Specify PAM  
PAM sequences should be specified in this format: a single on-target PAM followed by "|" and one or more off-target PAM (off-target PAMs should be separated by comma)  
for 5' PAMs, set `pamloc` to "5", and for 3' PAMs, set `pamloc` = "3"    
#### for SpCas9
```shell
pam="NGG|NGA,NAG" && pamloc="3"
```
#### for SpCas9-VQR
```shell
pam="NGA|NGG" && pamloc="3"
```
#### for enAsCas12a
```shell
pam="TTTV|TTTN" && pamloc="5"
```
</br>
## search gRNA  
This step takes ~1.5 hours for each genome  
```shell
mkdir gRNAs
bash search_gRNAs.sh GRCh38 Homo_sapiens "${pam}" "${pamloc}"
bash search_gRNAs.sh GRCm39 Mus_musculus "${pam}" "${pamloc}"
bash search_gRNAs.sh GRCz11 Danio_rerio "${pam}" "${pamloc}"
```
</br>

## Align gRNAs to the genome
### Build bwa indexes
```shell
bash build_bwa_index.sh GRCh38 Homo_sapiens
bash build_bwa_index.sh GRCm39 Mus_musculus
bash build_bwa_index.sh GRCz11 Danio_rerio
```
### Split the gRNA files into small chunks, for parallele computing on the hpc
This step takes ~50 minutes for each genome.
```shell
bash split4mapping.sh GRCh38 Homo_sapiens
bash split4mapping.sh GRCm39 Mus_musculus
bash split4mapping.sh GRCz11 Danio_rerio
```
### Run bwa on the hpc

prep work which includes determining the size of the array jobs for each genome
```shell
rm -rf bwa.slurm.out && mkdir bwa.slurm.out
chmod -R a+xX ./utils/FindOfftargetBwa/bin
job_array_size_h38=$(cat gRNAs/gRNA_GRCh38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz.out.split.tab | wc -l)
job_array_size_m39=$(cat gRNAs/gRNA_GRCm39/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz.out.split.tab | wc -l)
job_array_size_z11=$(cat gRNAs/gRNA_GRCz11/Danio_rerio.GRCz11.dna_sm.primary_assembly.fa.gz.out.split.tab | wc -l)
```
The following commands submit three batch jobs to the hpc, each takes ~20 hours to finish with 1000 cores  

```shell
bash map_gRNA_wrapper.sh GRCh38 Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa "${pam}" "${job_array_size_h38}" 100
bash map_gRNA_wrapper.sh GRCm39 Mus_musculus.GRCm39.dna_sm.primary_assembly.fa "${pam}" "${job_array_size_m39}" 100
bash map_gRNA_wrapper.sh GRCz11 Danio_rerio.GRCz11.dna_sm.primary_assembly.fa "${pam}" "${job_array_size_z11}" 100
```
## Compute off-target score
submit three batch jobs to the hpc, takes ~10 hours to finish with 1000 cores
```shell
rm -rf score.slurm.out && mkdir score.slurm.out
bash score_gRNA_wrapper.sh GRCh38 Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa "${pam}" "${job_array_size_h38}" 100
bash score_gRNA_wrapper.sh GRCm39 Mus_musculus.GRCm39.dna_sm.primary_assembly.fa "${pam}" "${job_array_size_m39}" 100
bash score_gRNA_wrapper.sh GRCz11 Danio_rerio.GRCz11.dna_sm.primary_assembly.fa "${pam}" "${job_array_size_z11}" 100
```
index the off-target scores
```shell
python scripts/index_gRNA_results.py --gzdir gRNAs/gRNA_GRCh38/gRNA.tab.gz.split.BwaMapped.scored
python scripts/index_gRNA_results.py --gzdir gRNAs/gRNA_GRCm39/gRNA.tab.gz.split.BwaMapped.scored
python scripts/index_gRNA_results.py --gzdir gRNAs/gRNA_GRCz11/gRNA.tab.gz.split.BwaMapped.scored
```
## Clean up and move precomputed files to destination
```shell
bash post_processing.sh "${pam}"
```
