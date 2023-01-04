#This script retrieves all N-,C-terminus tagging sites (N = end of ATG, C = start of stop codon)
#Default parameters will run the GRCh38 genome

import sys
sys.path.insert(1, '../../')
import os.path
from scripts.utils import *
import traceback

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

def parse_args():
    parser= MyParser(description='ProtospaceX')
    parser.add_argument('--genome_ver', default="GRCh38", type=str, help='currently supports three genomes: GRCh38, GRCm39, GRCz11, ', metavar='')
    parser.add_argument('--path2csv',   default="../../input/mart_export_canonical_proteincoding.csv", type=str,help='path to a csv file containing ENST information\n *required columns*: Ensemble_ID',metavar='')
    parser.add_argument('--outdir',   default="logs", type=str, help='output directory')

    config = parser.parse_args()
    return config

#configs
config = vars(parse_args())
genome = config['genome_ver']
#####################
##      main       ##
#####################
def main():
    try:
        #load gene model info
        print("loading gene model info")
        ENST_info = read_pickle_files(os.path.join("..","..","genome_files","parsed_gff3", config['genome_ver'],"ENST_info.pickle"))

        #load ENST list (use default: input/mart_export_canonical_proteincoding.csv)
        if os.path.isfile(config['path2csv']):
            log.info(f"begin processing user-supplied list of gene IDs in file {config['path2csv']}")
            df = pd.read_csv(os.path.join(config['path2csv']))
            # check csv columns
            keys2check = set(["Ensemble_ID"])
            if not keys2check.issubset(df.columns):
                log.error(f"Missing columns in the input csv file\n Required columns:\"Ensemble_ID\"")
                log.info(f"Please fix the input csv file and try again")
                sys.exit()
        else:
            file=config['path2csv']
            sys.exit("cannot open csv file {file}")

        #loop through each ENST, and get the coordinates (end of ATG, start of stop codon)
        o = open(f"canonical_N_C_sites_list_{genome}.csv", 'w')
        o.write("gene_ID_or_name,ref,chr,Ensemble_transcript_strand,coordinate,type")
        transcript_count = 0
        for index, row in df.iterrows():
            ENST_ID = row["Ensemble_ID"]
            try:
                ATG_loc, stop_loc = get_start_stop_loc(ENST_ID, ENST_info)
                end_of_ATG_loc = get_end_pos_of_ATG(ATG_loc)  # [chr, pos,strand]
                start_of_stop_loc = get_start_pos_of_stop(stop_loc)
                #print(f"{ENST_ID} {end_of_ATG_loc} {start_of_stop_loc}")
                o.write(f"\n{ENST_ID},ensembl_{genome}_latest,{end_of_ATG_loc[0]},{end_of_ATG_loc[2]},{end_of_ATG_loc[1]},end_of_ATG")
                o.write(f"\n{ENST_ID},ensembl_{genome}_latest,{start_of_stop_loc[0]},{end_of_ATG_loc[2]},{start_of_stop_loc[1]},start_of_STOP")
            except Exception as e:
                pass

            transcript_count += 1
            if transcript_count % 1000 == 0:
                print(f"processed {transcript_count} ENST IDs")

        print(f"processed {transcript_count} ENST IDs")
        o.close()

    except Exception as e:
        print("Unexpected error:", str(sys.exc_info()))
        traceback.print_exc()
        print("additional information:", e)
        PrintException()

##########################
## function definitions ##
##########################
def mkdir(mypath):
    if not os.path.exists(mypath):
        os.makedirs(mypath)

def get_end_pos_of_ATG(ATG_loc):
    """
    input: ATG_loc              [chr,start,end,strand] #start < end
    output: the pos of G in ATG [chr,pos,strand]       #start < end
    """
    strand = ATG_loc[3]
    if str(strand) == "+" or str(strand) == "1":
        return [ATG_loc[0],ATG_loc[2],ATG_loc[3]]
    else:
        return [ATG_loc[0],ATG_loc[1],ATG_loc[3]]

def PrintException():
    exc_type, exc_obj, tb = sys.exc_info()
    f = tb.tb_frame
    lineno = tb.tb_lineno
    filename = f.f_code.co_filename
    linecache.checkcache(filename)
    line = linecache.getline(filename, lineno, f.f_globals)
    print('EXCEPTION IN ({}, LINE {} "{}"): {}'.format(filename, lineno, line.strip(), exc_obj))

if __name__ == "__main__": main()

