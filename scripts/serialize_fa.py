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
    parser= MyParser(description='This script converts fasta sequences into pickle files for faster IO ')
    parser.add_argument('--fastagz', required = False, default="genome_files/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz", type=str, help='path to the the genome fa.gz file', metavar='')
    parser.add_argument('--fasta', required = False, default="", type=str, help='path to the the genome fa file', metavar='')
    config = parser.parse_args()
    return config

config = vars(parse_args())

#####################
##      main       ##
#####################
def main():
    try:
        #check input
        if (config["fastagz"] is None or config["fastagz"] == "") and (config["fasta"] is None or config["fasta"] == ""):
            log.error(f"please specify a valid fa.gz file")
            sys.exit("Please fix the error(s) above and rerun the script")

        starttime = datetime.datetime.now()
        seq_count = 0

        infile = config["fastagz"]
        if infile is None or infile == "":
            infile = config["fasta"]

        outdir = os.path.join("genome_files","fa_pickle",infile.split(".")[1])
        make_output_dir_if_nonexistent(outdir)

        with gzip.open(infile, "rt") as f:
            #wfh.write(f"seq_id\tlength\n")#header
            for record in SeqIO.parse(f, "fasta"):
                name = record.id
                # write to file
                pickle_filename = f"{name}.pk"
                with open(os.path.join(outdir,pickle_filename), 'wb') as handle:
                    pickle.dump(record, handle, protocol=pickle.HIGHEST_PROTOCOL)
                #print(f"chr:\t{record.id}\tlength:\t{len(record.seq)}\n}", flush=True)
                seq_count += 1

        endtime = datetime.datetime.now()
        elapsed_sec = endtime - starttime
        elapsed_min = elapsed_sec.seconds / 60
        #log.info(f"finished in {elapsed_min:.2f} min, processed {seq_count} seqs")
        print(f"finished in {elapsed_min:.2f} min ({elapsed_sec} sec), processed {seq_count} seqs")

    except Exception as e:
        print("Unexpected error:", str(sys.exc_info()))
        print("additional information:", e)
        PrintException()

##########################
## function definitions ##
##########################
def make_output_dir_if_nonexistent(path):
    if (not os.path.exists(path)): #output doesn't exist
      os.makedirs(path)

def PrintException():
    exc_type, exc_obj, tb = sys.exc_info()
    f = tb.tb_frame
    lineno = tb.tb_lineno
    filename = f.f_code.co_filename
    linecache.checkcache(filename)
    line = linecache.getline(filename, lineno, f.f_globals)
    print('EXCEPTION IN ({}, LINE {} "{}"): {}'.format(filename, lineno, line.strip(), exc_obj))


if __name__ == "__main__": main()
