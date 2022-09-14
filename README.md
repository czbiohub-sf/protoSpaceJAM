# ProtospaceXS
### Pre-compute all CRISPR guide RNAs and their properties in the human genome
### Rationale: The pre-computed guide RNAs can be later selected/filtered for various applications
##### "S" in the name stands for standalone (as opposed to protospaceX that relied on other web servers)


### Precompute gRNA
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

### Serialize fasta file for fast I/O
GRCh38  
```
scripts\serialize_fa.py --fastagz genome_files/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
```
GRCz11  
```
scripts\serialize_fa.py --fastagz genome_files/Danio_rerio.GRCz11.dna_sm.primary_assembly.fa.gz
```
