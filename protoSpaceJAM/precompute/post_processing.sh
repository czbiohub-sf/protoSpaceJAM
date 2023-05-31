rm genome_files/*.gz
rm genome_files/*.fa

rm genome_files/parsed_gff3/GRCh38/debug.txt 
rm genome_files/parsed_gff3/GRCm39/debug.txt 
rm genome_files/parsed_gff3/GRCz11/debug.txt

rm -r gRNAs/gRNA_GRCh38/gRNA.tab.gz
rm -r gRNAs/gRNA_GRCh38/gRNA.tab.gz.split
rm -r gRNAs/gRNA_GRCh38/gRNA.tab.gz.split.BwaMapped
rm gRNAs/gRNA_GRCh38/mapping_log.process.txt
rm gRNAs/gRNA_GRCh38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz.out.tab
rm gRNAs/gRNA_GRCh38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz.out.split.tab

rm -r gRNAs/gRNA_GRCm39/gRNA.tab.gz
rm -r gRNAs/gRNA_GRCm39/gRNA.tab.gz.split
rm -r gRNAs/gRNA_GRCm39/gRNA.tab.gz.split.BwaMapped
rm gRNAs/gRNA_GRCm39/mapping_log.process.txt
rm gRNAs/gRNA_GRCm39/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz.out.tab
rm gRNAs/gRNA_GRCm39/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz.out.split.tab

rm -r gRNAs/gRNA_GRCz11/gRNA.tab.gz
rm -r gRNAs/gRNA_GRCz11/gRNA.tab.gz.split
rm -r gRNAs/gRNA_GRCz11/gRNA.tab.gz.split.BwaMapped
rm gRNAs/gRNA_GRCz11/mapping_log.process.txt
rm gRNAs/gRNA_GRCz11/Danio_rerio.GRCz11.dna_sm.primary_assembly.fa.gz.out.tab
rm gRNAs/gRNA_GRCz11/Danio_rerio.GRCz11.dna_sm.primary_assembly.fa.gz.out.split.tab

mv genome_files/ ../genome_files
mv gRNAs/ ../precomputed_gRNAs
