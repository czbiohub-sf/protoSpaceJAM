script_folder="./scripts"
python ${script_folder}/index_gRNA_results.py --gzdir gRNA/gRNA_GRCh38/gRNA.tab.gz.split.BwaMapped.scored
python ${script_folder}/index_gRNA_results.py --gzdir gRNA/gRNA_GRCm39/gRNA.tab.gz.split.BwaMapped.scored
python ${script_folder}/index_gRNA_results.py --gzdir gRNA/gRNA_GRCz11/gRNA.tab.gz.split.BwaMapped.scored