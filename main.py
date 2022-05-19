import os.path
import pandas as pd
from subprocess import Popen
from Bio import SeqIO
from Bio.Seq import reverse_complement
import csv
import argparse
import sys
import linecache
import datetime
import gzip
import shutil
import gc
import math
import re
import pickle
from scripts.utils import *

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

def parse_args():
    parser= MyParser(description='This scripts creates a mapping of ENST to chr from gff3')
    parser.add_argument('--genome_ver', default="hg38", type=str, help='pickle file containing the ENST_info dict', metavar='')
    parser.add_argument('--path2csv', default="input.csv", type=str,
                        help='path to a csv file containing ENST information\n *required columns*: Ensemble_ID',
                        metavar='')
    config = parser.parse_args()
    return config

logging.setLoggerClass(ColoredLogger)
#logging.basicConfig()
log = logging.getLogger("get_scores")
log.propagate = False
log.setLevel(logging.INFO) #set the level of warning displayed

config = vars(parse_args())

#####################
##      main       ##
#####################
def main():
    try:
        starttime = datetime.datetime.now()

        #load gene model info
        ENST_info = read_pickle_files(os.path.join("genome_files","gff3",config['genome_ver'],"ENST_info.pickle"))

        #load gRNA info index (loc-to-filename)
        loc2file_index = read_pickle_files(os.path.join(f"gRNA_{config['genome_ver']}","loc2file_index.pickle"))

        #load ENST list
        df = pd.read_csv(os.path.join(config['path2csv']))
        # check csv columns
        keys2check = set(["Ensemble_ID"])
        if not keys2check.issubset(df.columns):
            log.error(f"Missing columns in the input csv file\n Required columns:\"Ensemble_ID\"")
            log.info(f"Please fix the input csv file and try again")
            sys.exit()

        #loop through each ENST
        for index, row in df.iterrows():
            ENST_ID = row["Ensemble_ID"]
            print(f"processing {ENST_ID}", flush=True)
            gRNAs = get_gRNAs(ENST_ID,ENST_info)

        #write csv out

        endtime = datetime.datetime.now()
        elapsed_sec = endtime - starttime
        elapsed_min = elapsed_sec.seconds / 60

        print(f"finished in {elapsed_min:.2f} min ({elapsed_sec} sec), pickle file(s) loaded", flush=True)
        log.info(f"finished in {elapsed_min:.2f} min ({elapsed_sec} sec) , pickle file(s) loaded")

    except Exception as e:
        print("Unexpected error:", str(sys.exc_info()))
        print("additional information:", e)
        PrintException()

##########################
## function definitions ##
##########################
def read_pickle_files(file):
    with open(file, 'rb') as handle:
        mydict = pickle.load(handle)
    return mydict

def get_gRNAs(ENST_ID,ENST_info):
    ATG_loc, stop_loc = get_start_stop_loc(ENST_ID,ENST_info)

    seq = get_seq(ATG_loc[0],ATG_loc[1],ATG_loc[2],ATG_loc[3])
    print (ATG_loc)
    print (f"start codon {seq}")

    seq = get_seq(stop_loc[0],stop_loc[1],stop_loc[2],stop_loc[3])
    print (stop_loc)
    print(f"stop codon {seq}")

def get_start_stop_loc(ENST_ID,ENST_info):
    my_transcript = ENST_info[ENST_ID]  # get the seq record
    # constructing the list of cds
    cdsList = [feat for feat in my_transcript.features if feat.type == 'CDS']
    CDS_first = cdsList[0]
    CDS_last = cdsList[len(cdsList) - 1]
    #get start codon location
    if CDS_first.strand==1:
        ATG_loc = [CDS_first.location.ref, CDS_first.location.start+0,CDS_first.location.start+2, 1] # format [start, end, strand]
    else:
        stop_loc = [CDS_first.location.ref, CDS_first.location.start+0,CDS_first.location.start+2, -1]
    #get stop codon location
    if CDS_last.strand==1:
        stop_loc = [CDS_last.location.ref, CDS_last.location.end-2,CDS_last.location.end+0, 1]
    else:
        ATG_loc = [CDS_last.location.ref, CDS_last.location.end-2,CDS_last.location.end+0, -1]

    return([ATG_loc, stop_loc])

def get_seq(chr,start,end,strand):
    '''
    chr
    start
    end
    strand 1 or -1 (str)
    '''
    chr_file_path = os.path.join("genome_files",f"{config['genome_ver']}_byChr",f"{chr}.pk")
    print(f"opening file {chr_file_path}")
    if os.path.isfile(chr_file_path):
        #read file
        chr_seqrecord = read_pickle_files(chr_file_path)
        subseq = chr_seqrecord.seq._data.decode()[(start-1):end]
        if strand == "-1" or strand == -1:
            return(reverse_complement(subseq))
        else:
            return(subseq)
    else:
        sys.exit(f"ERROR: file not found: {chr_file_path}")




def PrintException():
    exc_type, exc_obj, tb = sys.exc_info()
    f = tb.tb_frame
    lineno = tb.tb_lineno
    filename = f.f_code.co_filename
    linecache.checkcache(filename)
    line = linecache.getline(filename, lineno, f.f_globals)
    print('EXCEPTION IN ({}, LINE {} "{}"): {}'.format(filename, lineno, line.strip(), exc_obj))

if __name__ == "__main__": main()




