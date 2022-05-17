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

        #load gene model
        ENST_info = read_pickle_files(os.path.join("genome_files","gff3",config['genome_ver'],"ENST_info.pickle"))

        #load gRNA file
        loc2file_index = read_pickle_files(os.path.join(f"gRNA_{config['genome_ver']}","loc2file_index.pickle"))

        #load ENST list

        #loop through each ENST

        #write csv out

        endtime = datetime.datetime.now()
        elapsed_sec = endtime - starttime
        elapsed_min = elapsed_sec.seconds / 60
        log.info(f"finished in {elapsed_min:.2f} min ({elapsed_sec} sec) , pickle file(s) loaded")
        print(f"finished in {elapsed_min:.2f} min ({elapsed_sec} sec), pickle file(s) loaded", flush=True)

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

#def read_csv(file):
#    with open(file,"r") as handle:


def PrintException():
    exc_type, exc_obj, tb = sys.exc_info()
    f = tb.tb_frame
    lineno = tb.tb_lineno
    filename = f.f_code.co_filename
    linecache.checkcache(filename)
    line = linecache.getline(filename, lineno, f.f_globals)
    print('EXCEPTION IN ({}, LINE {} "{}"): {}'.format(filename, lineno, line.strip(), exc_obj))

if __name__ == "__main__": main()




