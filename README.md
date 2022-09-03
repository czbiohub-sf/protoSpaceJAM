# ProtospaceXS
### Pre-compute all CRISPR guide RNAs and their properties in the human genome
### Rationale: The pre-computed guide RNAs can be later selected/filtered for various applications
##### "S" in the name stands for standalone (as opposed to protospaceX that relied on other web servers)

### TO DO
- [x] Ignore (non-chimeric) homology arms in sliding windowing analysis
- [ ] Fix 888999 ignored by synoymoous mutations

### Precompute gRNA
GRCh38  
```

```
GRCz11  
```

```

### Preprocess GFF3 annotation file
GRCh38  
```
scripts\process_ENST_info.py --release 107 --genome_ver GRCh38
```
GRCz11  
```
scripts\process_ENST_info.py --release 107 --genome_ver GRCz11
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
