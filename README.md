# ProtospaceXS
- Guide RNA and repair donor design
- ssODN recoding to facilitate payload insertion and prevent recut
- Pre-compute all CRISPR guide RNAs and their properties for fast runtime

## Usage

make sure all script/binary files in folder utils are excutable
### Genome download and preprocessing
Download default genomes:   
GRCh38, GRCm39, GRCz11 (fasta + gff3)
```
cd genomefiles
bash download_genomes.sh
```
Preprocess GFF3 annotation file
```
#in directory genomefiles
bash preprocess_GFF3.sh
```
build bwa indexes
```
bash build_bwa_index.sh
```
Serialize fasta file for fast I/O
```
#in directory genomefiles
bash serialize_fa.sh
```
### Precompute gRNA
get all gRNAs
```
cd ..
cd precomuted_gRNAs
bash 1.get_gRNAs.sh  #Note: comment out line #12 if not using an hpc cluster
```
split into chunks
```
bash 2.split4mapping.sh
```
map to the genome
```
mkdir slurm.out
sbatch 3.map_gRNA.GRCh38.sh
```
get off-target score and efficiency score predictions
```
bash 4.score.sh
```
