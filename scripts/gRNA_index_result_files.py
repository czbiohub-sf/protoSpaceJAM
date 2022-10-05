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
sys.path.insert(1, '..')
#from gRNA_search import *
#from utils import *

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
# scoring_script_dir = os.path.join(sys.path[0],"..","crisporEffScores")
# sys.path.insert(1, scoring_script_dir)
# from cfd import *
# from crisporEffScores import *

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

def parse_args():
    parser= MyParser(description='')
    parser.add_argument('--gzdir', default="", type=str, help='path to the dir containing gzfiles', metavar='')
    config = parser.parse_args()
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
        if (config["gzdir"] is None or config["gzdir"] == ""):
            log.error(f"please specify --gzdir")
            sys.exit("Please fix the error(s) above and rerun the script")

        starttime = datetime.datetime.now()
        file_count = 0
        #make a dict for the index information will be collected
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
                    #print(line)
                    #read gRNA from gzip file
                    #seq = line.split("\t")[0]
                    #pam = line.split("\t")[1]
                    st = int(line.split("\t")[2])
                    en = int(line.split("\t")[3])

                    min_pos = update_min(min_pos, st)
                    min_pos = update_min(min_pos, en)
                    max_pos = update_max(max_pos, st)
                    max_pos = update_max(max_pos, en)

            file_index[Chr][f"{min_pos}-{max_pos}"] = filename
            file_count+=1

        # write dict to file
        with open(os.path.join(config["gzdir"],'loc2file_index.pickle'), 'wb') as handle:
            pickle.dump(file_index, handle, protocol=pickle.HIGHEST_PROTOCOL)

        endtime = datetime.datetime.now()
        elapsed_sec = endtime - starttime
        elapsed_min = elapsed_sec.seconds / 60
        log.info(f"finished in {elapsed_min:.2f} min, processed {file_count} files")
        print(f"finished in {elapsed_min:.2f} min, processed {file_count} files", flush=True)

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
    print('EXCEPTION IN ({}, LINE {} "{}"): {}'.format(filename, lineno, line.strip(), exc_obj))

if __name__ == "__main__": main()

