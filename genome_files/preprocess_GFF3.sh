#GRCh38
rm -rf parsed_gff3/GRCh38
mkdir parsed_gff3/GRCh38
python ../scripts/process_ENST_info.py --gff3_gz Homo_sapiens.GRCh38.107.gff3.gz --out_dir parsed_gff3/GRCh38

#GRCm39
rm -rf parsed_gff3/GRCm39
mkdir parsed_gff3/GRCm39
python ../scripts/process_ENST_info.py --gff3_gz Mus_musculus.GRCm39.107.gff3.gz --out_dir parsed_gff3/GRCm39

#GRCz11
rm -rf parsed_gff3/GRCz11
mkdir parsed_gff3/GRCz11
python ../scripts/process_ENST_info.py --gff3_gz Danio_rerio.GRCz11.107.gff3.gz --out_dir parsed_gff3/GRCz11

