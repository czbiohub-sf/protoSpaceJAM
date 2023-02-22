### ** Optional **: Precompute gRNA  
If you prefer downloading the precomputed gRNA, following the instructions [here](https://github.com/czbiohub/protoSpaceJAM)

### Prerequisite
All scripts and binary files in protoSpaceJAM/utils should be excutable

### Genome download and preprocessing
Download default genomes:   
GRCh38, GRCm39, GRCz11 (fasta + gff3)
```
cd ../genomefiles
bash download_genomes.sh
```
Preprocess GFF3 annotation file
```
bash preprocess_GFF3.sh
```
build bwa indexes Optional: for precomputing gRNA
```
bash build_bwa_index.sh
```
Serialize fasta file for fast I/O
```
#in genomefiles
bash serialize_fa.sh
```
Get all gRNAs
```
cd precomuted_gRNAs
bash 1.get_gRNAs.sh  #Note: comment out line #12 if not using an hpc cluster
```
Split into chunks
```
# in precomuted_gRNAs
bash 2.split4mapping.sh
```
Map to the genome  
This step is computationally intensive, it's recommended to run on an hpc cluster (with the slurm scheduler in this example)

```
working_dir="/hpc/projects/data_lg/duo.peng/protospaceXS/precomuted_gRNAs"
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
working_dir="/hpc/projects/data_lg/duo.peng/protospaceXS/precomuted_gRNAs"
cd $working_dir
rm -rf score.slurm.out
mkdir score.slurm.out
sbatch 4.score_gRNA.sh gRNA_GRCh38 Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa $working_dir
sbatch 4.score_gRNA.sh gRNA_GRCm39 Mus_musculus.GRCm39.dna_sm.primary_assembly.fa $working_dir
sbatch 4.score_gRNA.sh gRNA_GRCz11 Danio_rerio.GRCz11.dna_sm.primary_assembly.fa $working_dir
```
bash index the precompute gRNAs for fast access
```
# in precomuted_gRNAs
5.index_result_files.sh
```
