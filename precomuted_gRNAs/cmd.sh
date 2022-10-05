#######
# map #
#######
working_dir="/hpc/projects/data_lg/duo.peng/protospaceXS/precomuted_gRNAs"
rm -rf map.slurm.out
mkdir map.slurm.out
sbatch 3.map_gRNA.sh gRNA_GRCh38 Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa $working_dir

rm -rf map.slurm.out
mkdir map.slurm.out
sbatch 3.map_gRNA.sh gRNA_GRCm39 Mus_musculus.GRCm39.dna_sm.primary_assembly.fa $working_dir

rm -rf map.slurm.out
mkdir map.slurm.out
sbatch 3.map_gRNA.sh gRNA_GRCz11 Danio_rerio.GRCz11.dna_sm.primary_assembly.fa $working_dir

#########
# score #
#########
working_dir="/hpc/projects/data_lg/duo.peng/protospaceXS/precomuted_gRNAs"

rm -rf score.slurm.out
mkdir score.slurm.out
sbatch 4.score_gRNA.sh gRNA_GRCh38 Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa $working_dir

rm -rf score.slurm.out
mkdir score.slurm.out
sbatch 4.score_gRNA.sh gRNA_GRCm39 Mus_musculus.GRCm39.dna_sm.primary_assembly.fa $working_dir

rm -rf score.slurm.out
mkdir score.slurm.out
sbatch 4.score_gRNA.sh gRNA_GRCz11 Danio_rerio.GRCz11.dna_sm.primary_assembly.fa $working_dir


