# ProtospaceX  
A standalone program to design guide RNA and repair donors for CRISPR knock-in experiments  

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
  - automated selection of of ssDNA donor for maximum chance of payload insertion and prevent recut
  - scan and trim hard-to-synthesis motifs (coming soon)
- Use pre-computed genome-wide guide RNAs and their properties for fast runtime


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
#in genomefiles
bash preprocess_GFF3.sh
```
build bwa indexes
```
bash build_bwa_index.sh
```
Serialize fasta file for fast I/O
```
#in genomefiles
bash serialize_fa.sh
```
### Precompute gRNA (optional)
Get all gRNAs
```
cd ../precomuted_gRNAs
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
# in precomuted_gRNAs
mkdir map.slurm.out
#human
sbatch 3.map_gRNA.sh gRNA_GRCh38 Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
#mouse
sbatch 3.map_gRNA.sh gRNA_GRCm39 Mus_musculus.GRCm39.dna_sm.primary_assembly.fa
#zebrafish
sbatch 3.map_gRNA.sh gRNA_GRCz11 Danio_rerio.GRCz11.dna_sm.primary_assembly.fa
```
Get off-target score and efficiency score predictions  
This step is computationally intensive, it's recommended to run on an hpc cluster
```
# in precomuted_gRNAs
conda deactivate
conda activate sklearn0
bash 4.score.GRCh38.sh
```
index the precompute gRNAs for fast access
```
# in precomuted_gRNAs
5.index_result_files.sh
```

### Predict gRNA, DNA donor
```
conda deactivate
conda activate protospaceX
python main.py --path2csv input/test_protospaceX.csv --ssODN_max_size 200 --recoding_all
```
