# ProtospaceXS
- Guide RNA and repair donor design
- Rank guide RNAs for various applications
- ssODN design, and recoding to facilitate payload insertion and prevent recut
- Pre-compute all CRISPR guide RNAs and their properties for fast runtime

## Usage
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
