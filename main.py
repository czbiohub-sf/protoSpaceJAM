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
        freq_dict = dict()

        #load gene model info
        log.info("loading gene model info")
        ENST_info = read_pickle_files(os.path.join("genome_files","gff3",config['genome_ver'],"ENST_info.pickle"))

        #load gRNA info index (mapping of chromosomal location to file parts)
        log.info("loading the mapping of chromosomal location to file parts")
        loc2file_index = read_pickle_files(os.path.join(f"gRNA_{config['genome_ver']}","loc2file_index.pickle"))

        #load ENST list (the user input list)
        log.info("begin processing user-supplied list of gene IDs")
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
            log.info(f"processing {ENST_ID}", flush=True)
            gRNAs = get_gRNAs(ENST_ID = ENST_ID, ENST_info= ENST_info, freq_dict = freq_dict, loc2file_index= loc2file_index, dist = 100)

        #write csv out

        endtime = datetime.datetime.now()
        elapsed_sec = endtime - starttime
        elapsed_min = elapsed_sec.seconds / 60

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

def get_gRNAs(ENST_ID, ENST_info, freq_dict, loc2file_index, dist = 100):
    """
    input
        ENST_ID:
        ENST_info:
    return:
        a dictionary of guide RNAs which cuts <[dist] to end of start codon
        a dictionary of guide RNAs which cuts <[dist] to start of stop codon
    """
    #get location of start and stop location
    ATG_loc, stop_loc = get_start_stop_loc(ENST_ID,ENST_info)

    # #get start codon seq (for debug)
    # seq = get_seq(ATG_loc[0],ATG_loc[1],ATG_loc[2],ATG_loc[3])
    # print (ATG_loc)
    # print (f"start codon {seq}")
    # update_dict_count(seq,freq_dict)
    #
    # #get the exact stop codon (for debug)
    # seq = get_seq(stop_loc[0],stop_loc[1],stop_loc[2],stop_loc[3])
    # print (stop_loc)
    # print(f"stop codon {seq}")
    # update_dict_count(seq, freq_dict)

    #get gRNAs around the start codon
    log.debug(f"ATG_loc: {ATG_loc}")
    end_of_ATG_loc = get_end_pos_of_ATG(ATG_loc) # [chr, pos,strand]
    log.debug(f"end of the ATG: {end_of_ATG_loc}")
    df_gRNAs_ATG = get_gRNAs_near_loc(loc=end_of_ATG_loc,dist=100, loc2file_index=loc2file_index)

    #get gRNAs around the stop codon
    log.debug(f"stop_loc: {stop_loc}")
    start_of_stop_loc = get_start_pos_of_stop(stop_loc)
    log.debug(f"start of stop: {start_of_stop_loc}")# [chr, pos,strand]
    df_gRNAs_stop = get_gRNAs_near_loc(loc=start_of_stop_loc,dist=100, loc2file_index=loc2file_index)
    log.debug("a")

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

def get_start_pos_of_stop(stop_loc):
    """
    input: stop_loc                                 [chr,start,end,strand]  #start < end
    output: the pos of first base in the stop codon [chr,pos,strand]        #start < end
    """
    strand = stop_loc[3]
    if str(strand) == "+" or str(strand) == "1":
        return [stop_loc[0],stop_loc[1],stop_loc[3]]
    else:
        return [stop_loc[0],stop_loc[2],stop_loc[3]]

def get_start_stop_loc(ENST_ID,ENST_info):
    """
    input
    output: a list of two items
            ATG_loc: [chr,start,end,strand]  #start < end
            stop_loc: [chr,start,end,strand] #start < end
    """
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

def get_gRNAs_near_loc(loc,dist, loc2file_index):
    """
    input
        loc: [chr,pos,strand]
        dist: distance
    return:
        a dataframe of guide RNAs which cuts <[dist] to the loc
    """
    chr = loc[0]
    pos = loc[1]
    chr_dict = loc2file_index[chr]
    target_files =[] # a list of file names containing gRNAs near loc
    #lookup the file
    for key in chr_dict.keys():
        file_start = key.split("-")[0]
        file_end =  key.split("-")[1]
        interval = [file_start, file_end]
        if in_interval(pos, interval) or in_interval(pos-1000, interval) or in_interval(pos+1000, interval):
            target_files.append(chr_dict[key])
    #print(target_files)
    #load the gRNAs in the file
    dfs =[]
    for file in target_files:
        file_path = os.path.join(f"gRNA_{config['genome_ver']}","gRNA.tab.gz.split.BwaMapped.scored",file)
        df_tmp = pd.read_csv(file_path, sep="\t", compression='infer', header=None, names = ["seq","pam","start","end", "strand", "CSS", "ES"])
        dfs.append(df_tmp)
    df_gRNA = pd.concat(dfs)

    #subset gRNA based on strand  !ATTN: start > end when strand is '-'
    df_gRNA_on_sense = df_gRNA[(df_gRNA['strand'] == '+')]
    df_gRNA_on_antisense = df_gRNA[(df_gRNA['strand'] == '-')]
    #subset gRNAs and retain those cuts <[dist] to the loc
    df_gRNA_on_sense = df_gRNA_on_sense[(df_gRNA_on_sense['start'] > (pos-17-dist)) & (df_gRNA_on_sense['start'] < (pos-17+dist))]
    df_gRNA_on_antisense = df_gRNA_on_antisense[(df_gRNA_on_antisense['start'] > (pos+17-dist)) & (df_gRNA_on_antisense['start'] < (pos+17+dist))]

    return(pd.concat([df_gRNA_on_sense,df_gRNA_on_antisense]))

def in_interval(pos,interval):
    """
    check if pos in is interval
    input
        pos
        interval: [start,end]
    return: boolean
    """
    pos = int(pos)
    interval[0] = int(interval[0])
    interval[1] = int(interval[1])
    if pos >= interval[0] and pos <= interval[1]:
        log.debug(f"{pos} is in {interval}")
        return True
    else:
        return False

def update_dict_count(key,dict): # update the dictionary that keeps the count of each string (key)
    if key in dict.keys():
        dict[key]+=1
    else:
        dict[key]=1
    return dict

def get_seq(chr,start,end,strand):
    '''
    chr
    start
    end
    strand 1 or -1 (str)
    '''
    chr_file_path = os.path.join("genome_files",f"{config['genome_ver']}_byChr",f"{chr}.pk")
    log.debug(f"opening file {chr_file_path}")
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




