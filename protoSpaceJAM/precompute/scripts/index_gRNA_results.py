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
import pickle

from utils import *


def parse_args():
    parser = MyParser(description="")
    parser.add_argument(
        "--gzdir",
        default="",
        type=str,
        help="path to the dir containing gzfiles",
        metavar="",
    )
    config = parser.parse_args()
    return config


logging.setLoggerClass(ColoredLogger)
# logging.basicConfig()
log = logging.getLogger("index_gRNA_results")
log.propagate = False
log.setLevel(logging.INFO)  # set the level of warning displayed

config = vars(parse_args())
pam_len = 3
protospacer_len = 20
n_gRNA_per_chunk = 200
#####################
##      main       ##
#####################
def main():
    try:
        # check input
        if config["gzdir"] is None or config["gzdir"] == "":
            log.error(f"please specify --gzdir")
            sys.exit("Please fix the error(s) above and rerun the script")

        starttime = datetime.datetime.now()
        file_count = 0
        # make a dict for the index information will be collected
        file_index = dict()

        for filename in os.listdir(config["gzdir"]):
            file_path = os.path.join(config["gzdir"], filename)

            Chr = filename.split(".tab")[0].split(".part")[0]
            print(f"file: {filename}\tchr: {Chr}")

            if Chr not in file_index.keys():
                file_index[Chr] = dict()

            min_pos = 999999999
            max_pos = 0
            with gzip.open(file_path, "rt") as infh:
                for line in infh:
                    # print(line)
                    # read gRNA from gzip file
                    # seq = line.split("\t")[0]
                    # pam = line.split("\t")[1]
                    st = int(line.split("\t")[2])
                    en = int(line.split("\t")[3])

                    min_pos = update_min(min_pos, st)
                    min_pos = update_min(min_pos, en)
                    max_pos = update_max(max_pos, st)
                    max_pos = update_max(max_pos, en)

            file_index[Chr][f"{min_pos}-{max_pos}"] = filename
            file_count += 1

        # write dict to file
        with open(
            os.path.join(config["gzdir"], "loc2file_index.pickle"), "wb"
        ) as handle:
            pickle.dump(file_index, handle, protocol=pickle.HIGHEST_PROTOCOL)

        endtime = datetime.datetime.now()
        elapsed_sec = endtime - starttime
        elapsed_min = elapsed_sec.seconds / 60
        log.info(f"finished in {elapsed_min:.2f} min, processed {file_count} files")
        print(
            f"finished in {elapsed_min:.2f} min, processed {file_count} files",
            flush=True,
        )

    except Exception as e:
        print("Unexpected error:", str(sys.exc_info()))
        print("additional information:", e)
        PrintException()


##########################
## function definitions ##
##########################
def update_min(min, number):
    if number <= min:
        return number
    else:
        return min


def update_max(max, number):
    if number >= max:
        return number
    else:
        return max


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
