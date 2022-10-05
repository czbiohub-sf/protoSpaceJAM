########
#GRCh38#
########
dir=gRNA_GRCh38
script_folder="../scripts"
python ${script_folder}/gRNA_index_result_files.py --gzdir ${dir}/gRNA.tab.gz.split.BwaMapped.scored

########
#GRCm39#
########
dir=gRNA_GRCm39
script_folder="../scripts"
python ${script_folder}/gRNA_index_result_files.py --gzdir ${dir}/gRNA.tab.gz.split.BwaMapped.scored


########
#GRCz11#
########
dir=gRNA_GRCz11
script_folder="../scripts"
python ${script_folder}/gRNA_index_result_files.py --gzdir ${dir}/gRNA.tab.gz.split.BwaMapped.scored