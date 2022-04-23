# -*- coding: future_fstrings -*-
import os.path
import pandas as pd
from subprocess import Popen
from Bio import SeqIO
import csv
import argparse
import sys
import linecache
import datetime
import gzip
import shutil
import gc
import math
sys.path.insert(1, '..')
#from gRNA_search import *
#from utils import *
from crisporEffScores.crisporEffScores import *

#################
#custom logging #
#################
import logging

BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE = range(8)
#The background is set with 40 plus the number of the color, and the foreground with 30

#These are the sequences need to get colored ouput
RESET_SEQ = "\033[0m"
COLOR_SEQ = "\033[1;%dm"
BOLD_SEQ = "\033[1m"

def formatter_message(message, use_color = True):
    if use_color:
        message = message.replace("$RESET", RESET_SEQ).replace("$BOLD", BOLD_SEQ)
    else:
        message = message.replace("$RESET", "").replace("$BOLD", "")
    return message

COLORS = {
    'WARNING': YELLOW,
    'INFO': WHITE,
    'DEBUG': BLUE,
    'CRITICAL': YELLOW,
    'ERROR': RED
}

class ColoredFormatter(logging.Formatter):
    def __init__(self, msg, use_color = True):
        logging.Formatter.__init__(self, msg)
        self.use_color = use_color

    def format(self, record):
        levelname = record.levelname
        if self.use_color and levelname in COLORS:
            levelname_color = COLOR_SEQ % (30 + COLORS[levelname]) + levelname + RESET_SEQ
            record.levelname = levelname_color
        return logging.Formatter.format(self, record)

# Custom logger class with multiple destinations
class ColoredLogger(logging.Logger):
    FORMAT = "[$BOLD%(name)-1s$RESET][%(levelname)-1s]  %(message)s " #($BOLD%(filename)s$RESET:%(lineno)d)
    COLOR_FORMAT = formatter_message(FORMAT, True)
    def __init__(self, name):
        logging.Logger.__init__(self, name, logging.DEBUG)

        color_formatter = ColoredFormatter(self.COLOR_FORMAT)

        console = logging.StreamHandler()
        console.setFormatter(color_formatter)

        self.addHandler(console)
        return


#add to sys path, so the scoring modules can load
scoring_script_dir = os.path.join(sys.path[0],"..","crisporEffScores")
sys.path.insert(1, scoring_script_dir)
from cfd import *
from crisporEffScores import *

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

def parse_args():
    parser= MyParser(description='')
    parser.add_argument('--gzfile', default="", type=str, help='path to the gzfile', metavar='')
    parser.add_argument('--gzdir', default="", type=str, help='path to the dir containing gzfile', metavar='')
    parser.add_argument('--genome_fa', default="", type=str, help='name of genome fasta file', metavar='')
    config = parser.parse_args()
    if len(sys.argv)==1: # print help message if arguments are not valid
        parser.print_help()
        sys.exit(1)
    return config

logging.setLoggerClass(ColoredLogger)
#logging.basicConfig()
log = logging.getLogger("get_scores")
log.propagate = False
log.setLevel(logging.INFO) #set the level of warning displayed

config = vars(parse_args())
pam_len = 3
protospacer_len = 20
n_gRNA_per_chunk = 200
#####################
##      main       ##
#####################
def main():
    try:
        #check input
        if (config["gzfile"] is None or config["gzfile"] == ""):
            log.error(f"please specify --gzfile")
            sys.exit("Please fix the error(s) above and rerun the script")

        starttime = datetime.datetime.now()
        gRNA_count = 0
        #keep a dictionary of multi-mapping guides for current file
        gRNA_match = dict()

        #make output dir
        outdir_path =config["gzdir"] + ".scored"
        if os.path.isdir(outdir_path):
            #shutil.rmtree(outdir_path)  # remove existing dir
            pass
        else:
            os.makedirs(outdir_path)

        #process gRNA ()

        file = os.path.join(config["gzdir"],config["gzfile"])      
        Chr = file.rstrip('.tab.gz').split("/")[1].split(r'.')[0]
        #print(Chr)
        file_gRNA_count = 0

        with gzip.open(file, "rt") as fh, gzip.open(os.path.join(outdir_path,f"{config['gzfile'].rstrip('.gz')}.scored.gz"), "wt") as wgfh:

            for line in fh:
                #print(line)
                #read gRNA from gzip file
                seq = line.split("\t")[0]
                pam = line.split("\t")[1]
                st = line.split("\t")[2]
                en = line.split("\t")[3]
                strand = line.split("\t")[4].rstrip()
                leftFlank = line.split("\t")[5]
                rightFlank = line.split("\t")[6]
                gRNA_name = f"{Chr}_{st}_{en}_{strand}_{pam}"
                gRNA_name2 = f"{Chr}:{st}-{en}_{strand}"

                #####################
                #calculate cfd score#
                #####################
                if len(line.split("\t")[7])<=2:       # somes gRNAs does not have a bwa mapping, for example:  nnnnnnnnnnnnnngattca      tgg     257653  257675  + 
                    bwa_mapping = []
                else:
                    bwa_mapping = eval(line.split("\t")[7])

                #get bwa_mapping seqs
                bwa_mapping_seqs = [i.split('|')[2] for i in bwa_mapping]
                bwa_mapping_coords = [i.split('|')[0] for i in bwa_mapping]
                
                #remove on-target seq from bwa_mapping
                if (gRNA_name2 in bwa_mapping_coords):
                    idx = bwa_mapping_coords.index(gRNA_name2)
                    del(bwa_mapping_coords[idx])
                    del(bwa_mapping_seqs[idx])
                else:
                    print(f"warning: bwa hits does not contain the ontarget of {gRNA_name}")

                cfd_scores = [calcCfdScore(f"{seq}{pam}", otseq) for otseq in bwa_mapping_seqs]

                #print out the cfd scores together with alignments
                # print(gRNA_name2)
                # for idx, otseq in enumerate(bwa_mapping_seqs):
                #     score = calcCfdScore(f"{seq}{pam}", otseq)
                #     coord = bwa_mapping_coords[idx]
                #     print(f"{seq}{pam}\t{gRNA_name2}\n{otseq}\t{coord} \tCFD: {score}\n")

                #calc specificity score for the current guide
                cfd_scores = [s for s in cfd_scores if s]  #remove entries where cfd score failed to calcuate. For example, there is an offending "N" in NGTGGTGGTGGTAGAGGTGGTGG at position 6:167591372-167591394_-
                guideCfdScore = calcMitGuideScore(sum(cfd_scores))
                #print(guideCfdScore)

                ######################
                #calculate eff scores#
                ######################
                
                #input for calcAllScores() is  100bp sequences (50bp 5' of PAM, 50bp 3' of PAM) calculate all efficiency scores
                len_L = 50 - len(seq)
                len_R = 50 - len(pam)
                gRNA_with_flank = leftFlank[-len_L:] + seq + pam + rightFlank[:len_R]
                #print(len(gRNA_with_flank))
                #test 
                #print(sorted(calcAllScores(["CCACGTCTCCACACATCAGCACAACTACGCAGCGCCTCCCTCCACTCGGAAGGACTATCCTGCTGCCAAGAGGGTCAAGTTGGACAGTGTCAGAGTCCTG"]).items()))
                scores = calcAllScores([gRNA_with_flank])
                scores_flattened = flatten_score_dict(scores)
                #print(scores_flattened)

                ################################################
                #write to scores (with gRNA info) to a new file#
                ################################################
                wgfh.write(f"{seq}\t{pam}\t{st}\t{en}\t{strand}\t{guideCfdScore}\t{scores_flattened}\n")

                gRNA_count += 1

        endtime = datetime.datetime.now()
        elapsed_sec = endtime - starttime
        elapsed_min = elapsed_sec.seconds / 60
        log.info(f"finished in {elapsed_min:.2f} min, scored {gRNA_count} gRNAs in file {file}")
        print(f"finished in {elapsed_min:.2f} min, scored {gRNA_count} gRNAs in file {file}", flush=True)

    except Exception as e:
        print("Unexpected error:", str(sys.exc_info()))
        print("additional information:", e)
        PrintException()

##########################
## function definitions ##
##########################
def PrintException():
    exc_type, exc_obj, tb = sys.exc_info()
    f = tb.tb_frame
    lineno = tb.tb_lineno
    filename = f.f_code.co_filename
    linecache.checkcache(filename)
    line = linecache.getline(filename, lineno, f.f_globals)
    print('EXCEPTION IN ({}, LINE {} "{}"): {}'.format(filename, lineno, line.strip(), exc_obj))

def flatten_score_dict(score_dict):
    score_string = ""
    if 'fusi' in score_dict.keys():
        score_string = score_string + "|" + str(score_dict['fusi'][0])
    else:
        score_string = score_string + "|" + "NA"
    if 'fusiOld' in score_dict.keys():
        score_string = score_string + "|" + str(score_dict['fusiOld'][0])
    else:
        score_string = score_string + "|" + "NA"
    if 'chariRank' in score_dict.keys():
        score_string = score_string + "|" + str(score_dict['chariRank'][0])
    else:
        score_string = score_string + "|" + "NA"
    if 'ssc' in score_dict.keys():
        score_string = score_string + "|" + str(score_dict['ssc'][0])
    else:
        score_string = score_string + "|" + "NA"
    if 'wuCrispr' in score_dict.keys():
        score_string = score_string + "|" + str(score_dict['wuCrispr'][0])
    else:
        score_string = score_string + "|" + "NA"
    if 'doench' in score_dict.keys():
        score_string = score_string + "|" + str(score_dict['doench'][0])
    else:
        score_string = score_string + "|" + "NA"
    if 'wang' in score_dict.keys():
        score_string = score_string + "|" + str(score_dict['wang'][0])
    else:
        score_string = score_string + "|" + "NA"
    if 'crisprScan' in score_dict.keys(): #Moreno-Mateos
        score_string = score_string + "|" + str(score_dict['crisprScan'][0])
    else:
        score_string = score_string + "|" + "NA"
    if 'aziInVitro' in score_dict.keys():
        score_string = score_string + "|" + str(score_dict['aziInVitro'][0])
    else:
        score_string = score_string + "|" + "NA"
    if 'ccTop' in score_dict.keys():
        score_string = score_string + "|" + str(score_dict['ccTop'][0])
    else:
        score_string = score_string + "|" + "NA"
    if 'chariRaw' in score_dict.keys():
        score_string = score_string + "|" + str(score_dict['chariRaw'][0])
    else:
        score_string = score_string + "|" + "NA"
    if 'housden' in score_dict.keys():
        score_string = score_string + "|" + str(score_dict['housden'][0])
    else:
        score_string = score_string + "|" + "NA"                

    return score_string

if __name__ == "__main__": main()

