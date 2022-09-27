rm -rf map.slurm.out
mkdir map.slurm.out
sbatch 3.map_gRNA.sh gRNA_GRCh38 Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa

rm -rf map.slurm.out
mkdir map.slurm.out
sbatch 3.map_gRNA.sh gRNA_GRCm39 Mus_musculus.GRCm39.dna_sm.primary_assembly.fa

rm -rf map.slurm.out
mkdir map.slurm.out
sbatch 3.map_gRNA.sh gRNA_GRCz11 Danio_rerio.GRCz11.dna_sm.primary_assembly.fa