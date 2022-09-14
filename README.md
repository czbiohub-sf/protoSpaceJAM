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

Serialize fasta file for fast I/O
```
#in directory genomefiles
bash serialize_fa.sh
```
