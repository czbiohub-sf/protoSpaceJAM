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
log = logging.getLogger("ProtospaceX")
log.propagate = False
log.setLevel(logging.DEBUG) #set the level of warning displayed

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

        elapsed = cal_elapsed_time(starttime,datetime.datetime.now())
        log.info(f"finished loading in {elapsed[0]:.2f} min ({elapsed[1]} sec)")

        #loop through each ENST
        for index, row in df.iterrows():
            ENST_ID = row["Ensemble_ID"]
            log.info(f"processing {ENST_ID}")
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
def cal_elapsed_time(starttime,endtime):
    """
    output: [elapsed_min,elapsed_sec]
    """
    elapsed_sec = endtime - starttime
    elapsed_min = elapsed_sec.seconds / 60
    return[elapsed_min,elapsed_sec]

def read_pickle_files(file):
    with open(file, 'rb') as handle:
        mydict = pickle.load(handle)
    return mydict

def get_gRNAs(ENST_ID, ENST_info, freq_dict, loc2file_index, dist = 100):
    """
    input
        ENST_ID: ENST ID
        ENST_info: gene model info loaded from pickle file
    return:
        a dictionary of guide RNAs which cuts <[dist] to end of start codon
        a dictionary of guide RNAs which cuts <[dist] to start of stop codon
    """
    #get location of start and stop location
    ATG_loc, stop_loc = get_start_stop_loc(ENST_ID,ENST_info)

    #(for debug)
    # #get start codon seq
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

    ##################################
    #get gRNAs around the start codon#
    ##################################
    #get the start codon chromosomal location
    log.debug(f"ATG_loc: {ATG_loc}")
    end_of_ATG_loc = get_end_pos_of_ATG(ATG_loc) # [chr, pos,strand]
    log.debug(f"end of the ATG: {end_of_ATG_loc}")
    # get gRNA around the chromosomeal location (near ATG)
    df_gRNAs_ATG = get_gRNAs_near_loc(loc=end_of_ATG_loc,dist=100, loc2file_index=loc2file_index)
    # rank gRNAs
    ranked_df_gRNAs_ATG = rank_gRNAs_for_tagging(end_of_ATG_loc, df_gRNAs_ATG)


    ##################################
    #get gRNAs around the stop  codon#
    ##################################
    #get gRNAs around the stop codon
    log.debug(f"stop_loc: {stop_loc}")
    start_of_stop_loc = get_start_pos_of_stop(stop_loc)
    log.debug(f"start of stop: {start_of_stop_loc}")# [chr, pos,strand]
    # get gRNA around the chromosomeal location (near stop location)
    df_gRNAs_stop = get_gRNAs_near_loc(loc=start_of_stop_loc,dist=100, loc2file_index=loc2file_index)
    # rank gRNAs
    ranked_df_gRNAs_stop = rank_gRNAs_for_tagging(start_of_stop_loc, df_gRNAs_stop)


def rank_gRNAs_for_tagging(loc,gRNA_df):
    """
    input:  loc         [chr,pos,strand]  #start < end
            gRNA_df     pandas dataframe, *unranked*   columns: "seq","pam","start","end", "strand", "CSS", "ES"  !! neg strand: start > end
    output: gRNA_df     pandas dataframe *ranked*      columns: "seq","pam","start","end", "strand", "CSS", "ES"  !! neg strand: start > end
    """
    insPos = loc[1]
    #assign a score to each gRNA
    for index, row in gRNA_df.iterrows():
        start = row[2]
        end = row[3]
        strand = row[4]
        #Get cut to insert distance
        cutPos = None
        if strand == "+" or strand == "1" or strand == 1:
            cutPos = start + 16
        else:
            cutPos = start - 17
        cut2insDist = cutPos - insPos



        #calc. specificity_weight
        CSS = row[5]
        specificity_weight = _specificity_weight(CSS)

        #calc. distance_weight
        distance_weight = _dist_weight(hdr_dist = cut2insDist)

        #get position_weight


        log.debug(f"strand {strand} {start}-{end} cutPos {cutPos} loc {loc} cut2insDist {cut2insDist} CSS {CSS} specificity_weight {specificity_weight} distance_weight {distance_weight}")
        # final score = (specificity_weight)^alpha * distance_weight * position_weight


        #rank gRNAs based on the score

    pass

def _dist_weight(hdr_dist: int, _dist_weight_variance = 55) -> float:
    """
    taken from https://github.com/czbiohub/crispycrunch
    >>> _dist_weight(0)
    1.0
    >>> _dist_weight(5)
    0.7967034698934616
    >>> _dist_weight(10)
    0.402890321529133
    >>> _dist_weight(-20)
    0.026347980814448734
    """
    variance = _dist_weight_variance

    hdr_dist = abs(hdr_dist)  # make symmetric
    assert hdr_dist >= 0 and hdr_dist <= 100  # 100 is resonable upper bound

    # Returns a gaussian
    weight = math.exp((-1 * hdr_dist**2) / (2 * variance))
    assert weight >= 0 and weight <= 1
    return weight

def _specificity_weight(specificity_score: float, _specificity_weight_low = 45, _specificity_weight_high = 65 ):
    """
    taken from https://github.com/czbiohub/crispycrunch
    >>> _specificity_weight(20)
    0
    >>> _specificity_weight(60)
    0.75
    >>> _specificity_weight(80)
    1
    """
    low = _specificity_weight_low
    high = _specificity_weight_high

    assert specificity_score >= 0 and specificity_score <= 100
    if specificity_score <= low:
        return 0
    elif specificity_score >= high:
        return 1
    else:
        return 1 / (high - low) * (specificity_score - low)

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
    Get the chromosomal location of start and stop codons
    input: ENST_ID, ENST_info
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
        a dataframe of guide RNAs which cuts <[dist] to the loc, the columns are "seq","pam","start","end", "strand", "CSS", "ES"
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

    #add the cut position


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
        #log.debug(f"{pos} is in {interval}")
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




