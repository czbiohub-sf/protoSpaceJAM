# Precompute gRNA and genome information
If you prefer downloading the precomputed gRNA, following the instructions [here](https://github.com/czbiohub/protoSpaceJAM#download-and-unzip-pre-computed-data)

## Prerequisite
All scripts and binary files in protoSpaceJAM/utils should be excutable

## Prepare genomes
Download genome sequence and annotation files 
```shell
cd precompute
bash download_genomes.sh
```
Extract gene models and codon info (for protoSpaceJAM use)
```shell
module load anaconda
conda activate protospacejam
python scripts/extract_gene_models_info.py --gff3_gz genome_files/Homo_sapiens.GRCh38.109.gff3.gz --out_dir genome_files/parsed_gff3/GRCh38
python scripts/extract_gene_models_info.py --gff3_gz genome_files/Mus_musculus.GRCm39.109.gff3.gz --out_dir genome_files/parsed_gff3/GRCm39
python scripts/extract_gene_models_info.py --gff3_gz genome_files/Danio_rerio.GRCz11.109.gff3.gz --out_dir genome_files/parsed_gff3/GRCz11
bash serialize_fa.sh # serialize genome fasta for fast access
```
</br>

## Search for gRNAs in the genomes
This step takes ~1.5 hours for each genome  
PAM (if other than NGG) can be specificed by using passing the --pam argument to gRNA_scan_get_all_gRNA.py in `search_gRNAs.sh`  
```shell
mkdir gRNAs
bash search_gRNAs.sh GRCh38 Homo_sapiens
bash search_gRNAs.sh GRCm39 Mus_musculus
bash search_gRNAs.sh GRCz11 Danio_rerio
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

Prep work before running bwa to map the gRNAs to the genome.
```shell
rm -rf bwa.slurm.out && mkdir bwa.slurm.out
chmod -R a+xX ./utils/FindOfftargetBwa/bin
job_array_size_h38=$(ls -A gRNAs/gRNA_GRCh38/gRNA.tab.gz.split | wc -l)
job_array_size_m39=$(ls -A gRNAs/gRNA_GRCm39/gRNA.tab.gz.split | wc -l)
job_array_size_z11=$(ls -A gRNAs/gRNA_GRCz11/gRNA.tab.gz.split | wc -l)
echo "Please set the human gRNA job array size to: $job_array_size_h38"
echo "Please set the mouse gRNA job array size to: $job_array_size_m39"
echo "Please set the zebrafish gRNA job array size to: $job_array_size_z11"
```
The following commands submit three batch jobs to the hpc, each takes ~20 hours to finish with 1000 cores  
please set the correct array job size in `map_gRNA.sh` before executing the following commands.  
*For example, to create an array job of 3217 (and run 1000 simultaneously): `--array=1-3217s%1000`*
```shell
sbatch map_gRNA.sh GRCh38 Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa .
sbatch map_gRNA_mm39.sh GRCm39 Mus_musculus.GRCm39.dna_sm.primary_assembly.fa .
sbatch map_gRNA_z11.sh GRCz11 Danio_rerio.GRCz11.dna_sm.primary_assembly.fa .
```
## Compute off-target score
submit three batch jobs to the hpc, takes ~10 hours to finish with 1000 cores
```shell
sbatch score_gRNA.sh GRCh38 Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa .
sbatch score_gRNA.sh GRCm39 Mus_musculus.GRCm39.dna_sm.primary_assembly.fa .
sbatch score_gRNA.sh GRCz11 Danio_rerio.GRCz11.dna_sm.primary_assembly.fa .
#index the off-target scores
python scripts/index_gRNA_results.py --gzdir gRNAs/gRNA_GRCh38/gRNA.tab.gz.split.BwaMapped.scored
```
## Clean up and move precomputed files to destination
```shell
bash post_processing.sh
```