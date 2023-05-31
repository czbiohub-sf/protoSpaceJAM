#GRCh38
rm -rf genome_files/fa_pickle/GRCh38
mkdir genome_files/fa_pickle/GRCh38
python ./scripts/serialize_fa.py --fastagz genome_files/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz --outdir genome_files/fa_pickle/GRCh38

#GRCm39
rm -rf genome_files/fa_pickle/GRCm39
mkdir genome_files/fa_pickle/GRCm39
python ./scripts/serialize_fa.py --fastagz genome_files/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz --outdir genome_files/fa_pickle/GRCm39

#GRCz11
rm -rf genome_files/fa_pickle/GRCz11
mkdir genome_files/fa_pickle/GRCz11
python ./scripts/serialize_fa.py --fastagz genome_files/Danio_rerio.GRCz11.dna_sm.primary_assembly.fa.gz --outdir genome_files/fa_pickle/GRCz11

