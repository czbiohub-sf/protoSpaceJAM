#GRCh38
rm -rf fa_pickle/GRCh38
mkdir fa_pickle/GRCh38
python ../scripts/serialize_fa.py --fastagz Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz --outdir fa_pickle/GRCh38

#GRCm39
rm -rf fa_pickle/GRCm39
mkdir fa_pickle/GRCm39
python ../scripts/serialize_fa.py --fastagz Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz --outdir fa_pickle/GRCm39

#GRCz11
rm -rf fa_pickle/GRCz11
mkdir fa_pickle/GRCz11
python ../scripts/serialize_fa.py --fastagz Danio_rerio.GRCz11.dna_sm.primary_assembly.fa.gz --outdir fa_pickle/GRCz11

