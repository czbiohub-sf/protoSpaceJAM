#download fa.gz
wget http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz
wget http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.2.fa.gz
wget http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.3.fa.gz
wget http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.4.fa.gz
wget http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.5.fa.gz
wget http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.6.fa.gz
wget http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.7.fa.gz
wget http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.8.fa.gz
wget http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.9.fa.gz
wget http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.10.fa.gz
wget http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.11.fa.gz
wget http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.12.fa.gz
wget http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.13.fa.gz
wget http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.14.fa.gz
wget http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz
wget http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.16.fa.gz
wget http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.17.fa.gz
wget http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.18.fa.gz
wget http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.19.fa.gz
wget http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.20.fa.gz
wget http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.21.fa.gz
wget http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz
wget http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz
wget http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.X.fa.gz
wget http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.Y.fa.gz
wget http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.nonchromosomal.fa.gz

#serilize a fa.gz
python serialize_fa.py  --fastagz Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz
python serialize_fa.py  --fastagz Homo_sapiens.GRCh38.dna.chromosome.2.fa.gz
python serialize_fa.py  --fastagz Homo_sapiens.GRCh38.dna.chromosome.3.fa.gz
python serialize_fa.py  --fastagz Homo_sapiens.GRCh38.dna.chromosome.4.fa.gz
python serialize_fa.py  --fastagz Homo_sapiens.GRCh38.dna.chromosome.5.fa.gz
python serialize_fa.py  --fastagz Homo_sapiens.GRCh38.dna.chromosome.6.fa.gz
python serialize_fa.py  --fastagz Homo_sapiens.GRCh38.dna.chromosome.7.fa.gz
python serialize_fa.py  --fastagz Homo_sapiens.GRCh38.dna.chromosome.8.fa.gz
python serialize_fa.py  --fastagz Homo_sapiens.GRCh38.dna.chromosome.9.fa.gz
python serialize_fa.py  --fastagz Homo_sapiens.GRCh38.dna.chromosome.10.fa.gz
python serialize_fa.py  --fastagz Homo_sapiens.GRCh38.dna.chromosome.11.fa.gz
python serialize_fa.py  --fastagz Homo_sapiens.GRCh38.dna.chromosome.12.fa.gz
python serialize_fa.py  --fastagz Homo_sapiens.GRCh38.dna.chromosome.13.fa.gz
python serialize_fa.py  --fastagz Homo_sapiens.GRCh38.dna.chromosome.14.fa.gz
python serialize_fa.py  --fastagz Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz
python serialize_fa.py  --fastagz Homo_sapiens.GRCh38.dna.chromosome.16.fa.gz
python serialize_fa.py  --fastagz Homo_sapiens.GRCh38.dna.chromosome.17.fa.gz
python serialize_fa.py  --fastagz Homo_sapiens.GRCh38.dna.chromosome.18.fa.gz
python serialize_fa.py  --fastagz Homo_sapiens.GRCh38.dna.chromosome.19.fa.gz
python serialize_fa.py  --fastagz Homo_sapiens.GRCh38.dna.chromosome.20.fa.gz
python serialize_fa.py  --fastagz Homo_sapiens.GRCh38.dna.chromosome.21.fa.gz
python serialize_fa.py  --fastagz Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz
python serialize_fa.py  --fastagz Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz
python serialize_fa.py  --fastagz Homo_sapiens.GRCh38.dna.chromosome.X.fa.gz
python serialize_fa.py  --fastagz Homo_sapiens.GRCh38.dna.chromosome.Y.fa.gz
python serialize_fa.py  --fastagz Homo_sapiens.GRCh38.dna.nonchromosomal.fa.gz

#compare loading speed
#python load_fa.py --fastagz hg38_byChr/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz
#python load_pk.py --fastagzpk hg38_byChr/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz.pk
