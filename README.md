# ProtospaceX  
A standalone program to design Guide RNA and repair donors for CRISPR knock-in experiments  

## Key features:  
- Fully standalone, no calling to other bioinformatics servers
- Sophisticated guide RNA ranking system
  - specificity weight
  - Penalize cuts near exon-intron junctions etc.
  - Penalize cuts far away from the payload insertion site
- Sophisticated DNA donor design
  - recode to prevent recut
  - recode to facilitate payload insertion and prevent recut
  - center the DNA donor around the region containing the payload and recoded bases. 
  - enforce maximum length of DNA donor
  - scan and trim hard-to-synthesis motifs (coming soon)
- Pre-compute all CRISPR guide RNAs and their properties for fast runtime


## Usage

### Clone the repository
```
git clone https://github.com/czbiohub/protospaceX
```
### Go the repository directory, switch the branch if running branch other than master
```
cd protospaceXS
git checkout <branch you'd like to run>
```
### Create conda environment and activate it
```
conda env create -f environment1.yml
conda env create -f environment2.yml
conda activate protospaceX
```
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
### Precompute gRNA (optional)
Get all gRNAs
```
cd ..
cd precomuted_gRNAs
bash 1.get_gRNAs.sh  #Note: comment out line #12 if not using an hpc cluster
```
Split into chunks
```
bash 2.split4mapping.sh
```
Map to the genome  
This step is computationally intensive, it's recommended to run on an hpc cluster

```
mkdir slurm.out
sbatch 3.map_gRNA.GRCh38.sh
```
Get off-target score and efficiency score predictions  
This step is computationally intensive, it's recommended to run on an hpc cluster
```
conda deactivate
conda activate sklearn0
bash 4.score.GRCh38.sh
```

### Predict gRNA, DNA donor
```
conda deactivate
conda activate protospaceX
python main.py --path2csv input/test_protospaceX.csv --ssODN_max_size 200 --recoding_all
```
