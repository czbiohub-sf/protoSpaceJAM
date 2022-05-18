import os.path
import pandas as pd
from Bio import SeqIO
import csv
import argparse
import sys
import linecache
import datetime
import gzip
import shutil
import gc
import re
import pickle

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

def parse_args():
    parser= MyParser(description='This script get all the gRNA in the human genome, please specify either --fastagz or --fasta')
    parser.add_argument('--fastagzpk', default="", type=str, help='path to the human genome fasta file in pk format', metavar='')
    config = parser.parse_args()
    if len(sys.argv)==1: # print help message if arguments are not valid
        parser.print_help()
        sys.exit(1)
    return config

config = vars(parse_args())

#####################
##      main       ##
#####################
def main():
    try:
        #check input
        if (config["fastagzpk"] is None or config["fastagzpk"] == "") and (config["fasta"] is None or config["fasta"] == ""):
            log.error(f"please specify either --fastagzpk or --fasta")
            sys.exit("Please fix the error(s) above and rerun the script")

        starttime = datetime.datetime.now()
        seq_count = 0


        infile = config["fastagzpk"]
        if infile is None or infile == "":
            infile = config["fasta"]

        fa_pk = read_pickle_files(config["fastagzpk"])


        endtime = datetime.datetime.now()
        elapsed_sec = endtime - starttime
        elapsed_min = elapsed_sec.seconds / 60
        #log.info(f"finished in {elapsed_min:.2f} min, processed {seq_count} seqs")
        print(f"finished in {elapsed_min:.2f} min ({elapsed_sec} sec)")

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

def PrintException():
    exc_type, exc_obj, tb = sys.exc_info()
    f = tb.tb_frame
    lineno = tb.tb_lineno
    filename = f.f_code.co_filename
    linecache.checkcache(filename)
    line = linecache.getline(filename, lineno, f.f_globals)
    print('EXCEPTION IN ({}, LINE {} "{}"): {}'.format(filename, lineno, line.strip(), exc_obj))


if __name__ == "__main__": main()

