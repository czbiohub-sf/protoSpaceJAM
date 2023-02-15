![Webserver dependencies-0-brightgreen](https://user-images.githubusercontent.com/4129442/198696112-92ecc372-f3b5-4498-8cd9-4a01de0f851b.svg)
![status-active development-blueviolet](https://user-images.githubusercontent.com/4129442/198695999-a70bcd5f-c52e-4895-a1e7-d6b0da132812.svg)

# protoSpaceJAM 
A standalone program to design guide RNA and repair donors for CRISPR knock-in experiments  

## Key features:  
- Fully standalone, no calling to other bioinformatics servers
- Sophisticated guide RNA ranking system
  - Specificity weight
  - Penalize cuts near exon-intron junctions etc.
  - Penalize cuts far away from the payload insertion site
- Sophisticated DNA donor design
  - Recode to prevent recut
  - Recode to facilitate payload insertion and prevent recut
  - Center the DNA donor around the region containing the payload and recoded bases. 
  - Enforce maximum length of DNA donor
  - Automated selection of of ssDNA donor for maximum chance of payload insertion and prevent recut
  - Scan and trim hard-to-synthesis motifs (coming soon)
- Use/generate pre-computed genome-wide guide RNAs and their properties for fast runtime


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

### Predict gRNA, DNA donor
```
conda deactivate
conda activate protospaceX
python main.py --path2csv input/test_protospaceX.csv --ssODN_max_size 200 --recoding_all
```

### ** Optional **: Precompute gRNA (precomputed results available for download)
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

