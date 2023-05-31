#prepare output dir
outdir=gRNAs/gRNA_${1}
rm -rf ${outdir}
#scan gRNAs
script_folder="./scripts"
genome_file_dir="genome_files"
genome_file="${2}.${1}.dna_sm.primary_assembly.fa.gz"
python ${script_folder}/gRNA_scan_get_all_gRNA.py --fastagz ${genome_file_dir}/${genome_file} --outdir ${outdir}
