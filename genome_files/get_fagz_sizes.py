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

from protoSpaceJAM.util.utils import MyParser


def parse_args():
    parser = MyParser(
        description="This script get all the gRNA in the human genome, please specify either --fastagz or --fasta"
    )
    parser.add_argument(
        "--fastagz",
        default="",
        type=str,
        help="path to the human genome fasta file",
        metavar="",
    )
    config = parser.parse_args()
    if len(sys.argv) == 1:  # print help message if arguments are not valid
        parser.print_help()
        sys.exit(1)
    return config


protosp_len = 20
PAM = "NGG"

config = vars(parse_args())

#####################
##      main       ##
#####################
def main():
    try:
        # check input
        if (config["fastagz"] is None or config["fastagz"] == "") and (
            config["fasta"] is None or config["fasta"] == ""
        ):
            log.error(f"please specify either --fastagz or --fasta")
            sys.exit("Please fix the error(s) above and rerun the script")

        starttime = datetime.datetime.now()
        seq_count = 0

        infile = config["fastagz"]
        if infile is None or infile == "":
            infile = config["fasta"]

        with gzip.open(infile, "rt") as f, open(f"{infile}.sizes", "w") as wfh:
            # wfh.write(f"seq_id\tlength\n")#header
            for record in SeqIO.parse(f, "fasta"):

                seq_count += 1
                # print(f"chr:\t{record.id}\tlength:\t{len(record.seq)}\n}", flush=True)

                # write tab file
                wfh.write(f"{record.id}\t{len(record.seq)}\n")
                wfh.flush()

        endtime = datetime.datetime.now()
        elapsed_sec = endtime - starttime
        elapsed_min = elapsed_sec.seconds / 60
        # log.info(f"finished in {elapsed_min:.2f} min, processed {seq_count} seqs")
        # print(f"finished in {elapsed_min:.2f} min, processed {gene_count} genes, designed {gRNA_count} gRNAs")

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
    print(
        'EXCEPTION IN ({}, LINE {} "{}"): {}'.format(
            filename, lineno, line.strip(), exc_obj
        )
    )


if __name__ == "__main__":
    main()
