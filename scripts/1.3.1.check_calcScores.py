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

sys.path.insert(1, "..")
# from gRNA_search import *
# from utils import *
from crisporEffScores.crisporEffScores import *

#################
# custom logging #
#################
import logging

BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE = range(8)
# The background is set with 40 plus the number of the color, and the foreground with 30

# These are the sequences need to get colored ouput
RESET_SEQ = "\033[0m"
COLOR_SEQ = "\033[1;%dm"
BOLD_SEQ = "\033[1m"


def formatter_message(message, use_color=True):
    if use_color:
        message = message.replace("$RESET", RESET_SEQ).replace("$BOLD", BOLD_SEQ)
    else:
        message = message.replace("$RESET", "").replace("$BOLD", "")
    return message


COLORS = {
    "WARNING": YELLOW,
    "INFO": WHITE,
    "DEBUG": BLUE,
    "CRITICAL": YELLOW,
    "ERROR": RED,
}


class ColoredFormatter(logging.Formatter):
    def __init__(self, msg, use_color=True):
        logging.Formatter.__init__(self, msg)
        self.use_color = use_color

    def format(self, record):
        levelname = record.levelname
        if self.use_color and levelname in COLORS:
            levelname_color = (
                COLOR_SEQ % (30 + COLORS[levelname]) + levelname + RESET_SEQ
            )
            record.levelname = levelname_color
        return logging.Formatter.format(self, record)


# Custom logger class with multiple destinations
class ColoredLogger(logging.Logger):
    FORMAT = "[$BOLD%(name)-1s$RESET][%(levelname)-1s]  %(message)s "  # ($BOLD%(filename)s$RESET:%(lineno)d)
    COLOR_FORMAT = formatter_message(FORMAT, True)

    def __init__(self, name):
        logging.Logger.__init__(self, name, logging.DEBUG)

        color_formatter = ColoredFormatter(self.COLOR_FORMAT)

        console = logging.StreamHandler()
        console.setFormatter(color_formatter)

        self.addHandler(console)
        return


# add to sys path, so the scoring modules can load
scoring_script_dir = os.path.join(sys.path[0], "..", "crisporEffScores")
sys.path.insert(1, scoring_script_dir)
from cfd import *
from crisporEffScores import *


class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write("error: %s\n" % message)
        self.print_help()
        sys.exit(2)


def parse_args():
    parser = MyParser(description="")
    parser.add_argument(
        "--gzfile", default="", type=str, help="path to the gzfile", metavar=""
    )
    parser.add_argument(
        "--gzdir",
        default="",
        type=str,
        help="path to the dir containing gzfile",
        metavar="",
    )
    parser.add_argument(
        "--genome_fa",
        default="",
        type=str,
        help="name of genome fasta file",
        metavar="",
    )
    config = parser.parse_args()
    if len(sys.argv) == 1:  # print help message if arguments are not valid
        parser.print_help()
        sys.exit(1)
    return config


logging.setLoggerClass(ColoredLogger)
# logging.basicConfig()
log = logging.getLogger("get_scores")
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
        if config["gzfile"] is None or config["gzfile"] == "":
            log.error(f"please specify --gzfile")
            sys.exit("Please fix the error(s) above and rerun the script")

        starttime = datetime.datetime.now()
        gRNA_count = 0
        scored_gRNA_count = 0
        # keep a dictionary of multi-mapping guides for current file
        gRNA = dict()
        scored_gRNA = dict()

        # make output dir
        outdir_path = config["gzdir"] + ".scored_failedCheck"
        if os.path.isdir(outdir_path):
            # shutil.rmtree(outdir_path)  # remove existing dir
            pass
        else:
            os.makedirs(outdir_path)

        file = os.path.join(config["gzdir"], config["gzfile"])
        Chr = file.rstrip(".tab.gz").split("/")[1].split(r".")[0]
        # print(Chr)

        with gzip.open(file, "rt") as infh:
            for line in infh:
                # print(line)
                # read gRNA from gzip file
                seq = line.split("\t")[0]
                pam = line.split("\t")[1]
                st = line.split("\t")[2]
                en = line.split("\t")[3]
                strand = line.split("\t")[4].rstrip()
                gRNA_name = f"{Chr}_{st}_{en}_{strand}_{pam}"
                gRNA_name2 = f"{Chr}:{st}-{en}_{strand}"
                gRNA[gRNA_name] = 1
                gRNA_count += 1

        result_file = os.path.join(
            config["gzdir"] + ".scored", f"{config['gzfile'].rstrip('.gz')}.scored.gz"
        )
        if os.path.isfile(result_file):  # score result file exists
            with gzip.open(result_file, "rt") as outfh:
                for line in outfh:
                    # print(line)
                    # read gRNA from gzip file
                    seq = line.split("\t")[0]
                    pam = line.split("\t")[1]
                    st = line.split("\t")[2]
                    en = line.split("\t")[3]
                    strand = line.split("\t")[4].rstrip()
                    gRNA_name = f"{Chr}_{st}_{en}_{strand}_{pam}"
                    gRNA_name2 = f"{Chr}:{st}-{en}_{strand}"
                    scored_gRNA[gRNA_name] = 1
            if len(scored_gRNA) < len(gRNA):  # check the number of scored gRNAs
                with gzip.open(
                    os.path.join(outdir_path, config["gzfile"]), "wt"
                ) as gzfh:
                    gzfh.write(f"scored_gRNA_count < gRNA_count")
        else:  # score result file doesn't exist
            with gzip.open(os.path.join(outdir_path, config["gzfile"]), "wt") as gzfh:
                gzfh.write(f"scored gRNA file doesn't exist")

        endtime = datetime.datetime.now()
        elapsed_sec = endtime - starttime
        elapsed_min = elapsed_sec.seconds / 60
        log.info(
            f"finished in {elapsed_min:.2f} min, checked {gRNA_count} gRNAs in file {file}"
        )
        print(
            f"finished in {elapsed_min:.2f} min, checked {gRNA_count} gRNAs in file {file}",
            flush=True,
        )

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
