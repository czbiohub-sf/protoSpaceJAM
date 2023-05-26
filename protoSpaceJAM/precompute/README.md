# Precompute gRNA and genome information
If you prefer downloading the precomputed gRNA, following the instructions [here](https://github.com/czbiohub/protoSpaceJAM#download-and-unzip-pre-computed-data)

## Prerequisite
All scripts and binary files in protoSpaceJAM/utils should be excutable

## Prepare genomes
Download genome sequence and annotation files 
```
cd precompute
bash download_genomes.sh
```
Extract gene models and codon info (for protoSpaceJAM use)
```
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
```
mkdir gRNAs
bash search_gRNAs.sh GRCh38 Homo_sapiens
bash search_gRNAs.sh GRCm39 Mus_musculus
bash search_gRNAs.sh GRCz11 Danio_rerio
```
</br>

## Align gRNAs to the genome
### Build bwa indexes
```
bash build_bwa_index.sh GRCh38 Homo_sapiens
bash build_bwa_index.sh GRCm39 Mus_musculus
bash build_bwa_index.sh GRCz11 Danio_rerio
```
### Split the gRNA files into small chunks, for parallele computing on the hpc
This step takes ~50 minutes for each genome.
```
bash split4mapping.sh GRCh38 Homo_sapiens
bash split4mapping.sh GRCm39 Mus_musculus
bash split4mapping.sh GRCz11 Danio_rerio
```
### Run bwa on the hpc
submit three batch jobs to the hpc, takes ~20 hours to finish with 1700 cores
```
rm -rf bwa.slurm.out && mkdir bwa.slurm.out
chmod -R a+xX /hpc/projects/data_lg/duo.peng/pJAM_precompute/utils/FindOfftargetBwa/bin
sbatch map_gRNA.sh GRCh38 Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa .
sbatch map_gRNA.sh GRCh38 Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa .
sbatch map_gRNA.sh GRCh38 Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa .
```
## Compute off-target score
submit three batch jobs to the hpc, takes ~20 hours to finish with 1700 cores
```
rm -rf score.slurm.out && mkdir score.slurm.out
sbatch score_gRNA.sh GRCh38 Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa .
sbatch score_gRNA.sh GRCm39 Mus_musculus.GRCm39.dna_sm.primary_assembly.fa .
sbatch score_gRNA.sh GRCz11 Danio_rerio.GRCz11.dna_sm.primary_assembly.fa .
```








</br></br></br></br></br></br></br></br></br></br>
Serialize fasta file for fast I/O
```
#in genomefiles
bash serialize_fa.sh
```
Get all gRNAs
```
cd precomputed_gRNAs
bash 1.get_gRNAs.sh  #Note: comment out line #12 if not using an hpc cluster
```
Split into chunks
```
# in precomputed_gRNAs
bash 2.split4mapping.sh
```
Map to the genome  
This step is computationally intensive, it's recommended to run on an hpc cluster (with the slurm scheduler in this example)

```
working_dir="/PATH/TO/REPO/precomputed_gRNAs"
cd $working_dir
rm -rf map.slurm.out
mkdir map.slurm.out
sbatch 3.map_gRNA.sh gRNA_GRCh38 Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa $working_dir
sbatch 3.map_gRNA.sh gRNA_GRCm39 Mus_musculus.GRCm39.dna_sm.primary_assembly.fa $working_dir
sbatch 3.map_gRNA.sh gRNA_GRCz11 Danio_rerio.GRCz11.dna_sm.primary_assembly.fa $working_dir
```
Get off-target score and efficiency score predictions  
This step is computationally intensive, it's recommended to run on an hpc cluster (with the slurm scheduler in this example)
```
working_dir="/PATH/TO/REPO/precomputed_gRNAs"
cd $working_dir
rm -rf score.slurm.out
mkdir score.slurm.out
sbatch 4.score_gRNA.sh gRNA_GRCh38 Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa $working_dir
sbatch 4.score_gRNA.sh gRNA_GRCm39 Mus_musculus.GRCm39.dna_sm.primary_assembly.fa $working_dir
sbatch 4.score_gRNA.sh gRNA_GRCz11 Danio_rerio.GRCz11.dna_sm.primary_assembly.fa $working_dir
```
bash index the precompute gRNAs for fast access
```
# in precomputed_gRNAs
5.index_result_files.sh
```
