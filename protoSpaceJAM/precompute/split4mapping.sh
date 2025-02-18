#split4mapping.sh <chr> <genome_file>
dir=gRNAs/gRNA_${1}
genome_file="${2}.${1}.dna_sm.toplevel.fa.gz"
genome_file_path="${dir}/${genome_file}"
# check if the genome file exists, if not try primary assembly
if [ ! -f ${genome_file_path} ]; then
    genome_file="${2}.${1}.dna_sm.primary_assembly.fa.gz"
fi

script_folder="./scripts"
python ${script_folder}/split4mapping.py --dir ${dir}/gRNA.tab.gz --tab ${dir}/${genome_file}.out.tab --part_size 200000