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
#from BLAST_utils import check_blastDB_human

sys.path.insert(1, '..')
from gRNA_search import *
from utils import *

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

def parse_args():
    parser= MyParser(description='')
    parser.add_argument('--gzfile', default="", type=str, help='path to the gzfile', metavar='')
    parser.add_argument('--gzdir', default="", type=str, help='path to the dir containing gzfile', metavar='')
    parser.add_argument('--gRNA_count', default="", type=int, help='path to the gzfile', metavar='')
    parser.add_argument('--genome_fa', default="", type=str, help='name of genome fasta file', metavar='')
    parser.add_argument('--thread', default="", type=str, help='num of threads to use for bwa', metavar='')
    config = parser.parse_args()
    if len(sys.argv)==1: # print help message if arguments are not valid
        parser.print_help()
        sys.exit(1)
    return config

logging.setLoggerClass(ColoredLogger)
#logging.basicConfig()
log = logging.getLogger("map_gRNA")
log.propagate = False
log.setLevel(logging.INFO) #set the level of warning displayed

config = vars(parse_args())
genome_fa = config["genome_fa"]
thread2use = config["thread"]
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
        outdir_path =config["gzdir"] + ".BwaMapped_failedCheck"
        if os.path.isdir(outdir_path):
            #shutil.rmtree(outdir_path)  # remove existing dir
            pass
        else:
            os.makedirs(outdir_path)

        #process gRNA (bwa and parse in chunks)
        total_chunk_count = str(math.ceil(int(str(config["gRNA_count"]).strip())/n_gRNA_per_chunk))
        file = os.path.join(config["gzdir"],config["gzfile"])      
        Chr = file.rstrip('.tab.gz').split("/")[1].split(r'.')[0]

        gRNA = dict()
        mapped_gRNA = dict()

        with  gzip.open(file, "rt") as infh:
            for line in infh:
                #print(line)
                #read gRNA from gzip file
                seq = line.split("\t")[0]
                pam = line.split("\t")[1]
                st = line.split("\t")[2]
                en = line.split("\t")[3]
                strand = line.split("\t")[4].rstrip()
                gRNA_name = f"{Chr}_{st}_{en}_{strand}_{pam}"
                gRNA_name2 = f"{Chr}:{st}-{en}_{strand}"
                gRNA[gRNA_name] = 1
                gRNA_count += 1


        with gzip.open(os.path.join(config["gzdir"] + ".BwaMapped",f"{config['gzfile'].rstrip('.gz')}.mapped.gz"), "rt") as outfh:
            for line in outfh:
                #print(line)
                #read gRNA from gzip file
                seq = line.split("\t")[0]
                pam = line.split("\t")[1]
                st = line.split("\t")[2]
                en = line.split("\t")[3]
                strand = line.split("\t")[4].rstrip()
                gRNA_name = f"{Chr}_{st}_{en}_{strand}_{pam}"
                gRNA_name2 = f"{Chr}:{st}-{en}_{strand}"
                mapped_gRNA[gRNA_name] = 1


        if len(mapped_gRNA) < len(gRNA):
            with gzip.open(os.path.join(outdir_path,f"{config['gzfile'].rstrip('.gz')}.mapped.gz"), "wt") as gzfh:
                gzfh.write(f"scored_gRNA_count < gRNA_count")    


        endtime = datetime.datetime.now()
        elapsed_sec = endtime - starttime
        elapsed_min = elapsed_sec.seconds / 60
        log.info(f"finished in {elapsed_min:.2f} min, checked {gRNA_count} gRNAs in file {file}")
        print(f"finished in {elapsed_min:.2f} min, checked {gRNA_count} gRNAs in file {file}", flush=True)

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

def BLAST(fastfile):
    numThread2use = max([1, os.cpu_count() - 2])  # for BLAST, use all CPUs except 2
    #numThread2use = 1
    BLAST_bin, exe_suffix, BLAST_db_path = check_blastDB_human()
    query = fastfile
    # start BLAST
    cmd = [f"{BLAST_bin}blastn{exe_suffix}", "-task", "blastn-short", "-query", f"{query}", "-db", f"{BLAST_db_path}",
           "-num_threads", f"{numThread2use}", "-perc_identity", "100", "-outfmt",
           "6 qseqid sseqid qstart qend sstart send pident", "-out", f"{query}.out", "-ungapped"]
    p = Popen(cmd, universal_newlines=True)
    p.communicate()  # now wait for the process to finish
    return str(p.returncode)


if __name__ == "__main__": main()

