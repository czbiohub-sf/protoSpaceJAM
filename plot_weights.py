import os.path
import pandas as pd
import argparse
import sys
import linecache
import datetime
import matplotlib.pyplot as plt
import numpy as np
import pickle
from scripts.utils import *

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

def parse_args():
    parser= MyParser(description='This scripts creates a mapping of ENST to chr from gff3')
    parser.add_argument('--num_to_process', default="300", type=str, help='this parameter decides which file to load, the files have name start/stop_gRNAs_of_{num_to_process}_genes.pickle', metavar='')
    config = parser.parse_args()
    return config

logging.setLoggerClass(ColoredLogger)
#logging.basicConfig()
log = logging.getLogger("plot weights")
log.propagate = False
log.setLevel(logging.INFO) #set the level of warning displayed
#log.setLevel(logging.DEBUG) #set the level of warning displayed

config = vars(parse_args())

#####################
##      main       ##
#####################
def main():
    try:
        starttime = datetime.datetime.now()


        #load gRNA dataframes
        log.info("loading pickle files")
        num_to_process = config["num_to_process"]
        start_gRNA_df = read_pickle_files(os.path.join(f"start_gRNAs_of_{num_to_process}_genes.pickle"))
        stop_gRNA_df = read_pickle_files(os.path.join(f"stop_gRNAs_of_{num_to_process}_genes.pickle"))

        #plot dist_weight
        plot_hist(df = start_gRNA_df,  col = "dist_weight", bin_num=40, num_to_process = num_to_process)
        plot_hist(df = start_gRNA_df,  col = "pos_weight", bin_num=40, num_to_process = num_to_process)
        plot_hist(df = start_gRNA_df,  col = "spec_weight", bin_num=40, num_to_process = num_to_process)
        plot_hist(df = start_gRNA_df,  col = "final_weight", bin_num=40, num_to_process = num_to_process)

        #plot dist_weight (with overflow bin)
        # lower = 0
        # upper = 0.005
        # bin_size = 0.0005
        # bins = np.arange(lower, upper, bin_size)
        # fig, ax = plt.subplots(figsize=(9, 5))
        # _, bins, patches = plt.hist([np.clip(d, bins[0], bins[-1])],
        #                             # normed=1,  # normed is deprecated; replace with density
        #                             density=True,
        #                             bins=bins, color=['#607c8e'])
        # xlabels = bins[1:].astype(str)
        # xlabels[-1] += '+'
        # N_labels = len(xlabels)
        # plt.xlim([lower, upper])
        # plt.xticks(bin_size * np.arange(N_labels) )
        # ax.set_xticklabels(xlabels)
        # plt.yticks([])
        # plt.title('')
        # plt.setp(patches, linewidth=0)
        # plt.legend(loc='upper left')
        # fig.tight_layout()
        #
        # plt.savefig(f"{name}.overflow.png")



        #report time used
        elapsed = cal_elapsed_time(starttime,datetime.datetime.now())
        log.info(f"finished loading in {elapsed[0]:.2f} min ({elapsed[1]} sec)")



        #
        endtime = datetime.datetime.now()
        elapsed_sec = endtime - starttime
        elapsed_min = elapsed_sec.seconds / 60

        log.info(f"finished in {elapsed_min:.2f} min ({elapsed_sec} sec)")


    except Exception as e:
        print("Unexpected error:", str(sys.exc_info()))
        print("additional information:", e)
        PrintException()

##########################
## function definitions ##
##########################
def plot_hist(df, col, bin_num, num_to_process):
    name = col
    d = df[name]
    data_volume = len(d)
    d.plot.hist(grid=True, bins=bin_num, rwidth=0.9, color='#607c8e')
    plt.title(f"Distribution of {name} for {data_volume} gRNAs in {num_to_process} genes")
    plt.xlabel(name)
    plt.ylabel('Counts')
    plt.grid(axis='y', alpha=0.75)
    plt.savefig(f"{name}.{bin_num}bins.png")
    plt.close()

def cal_elapsed_time(starttime,endtime):
    """
    output: [elapsed_min,elapsed_sec]
    """
    elapsed_sec = endtime - starttime
    elapsed_min = elapsed_sec.seconds / 60
    return[elapsed_min,elapsed_sec]

def read_pickle_files(file):
    if os.path.isfile(file):
        with open(file, 'rb') as handle:
            mydict = pickle.load(handle)
        return mydict
    else:
        sys.exit(f"Cannot open file: {file}")

def PrintException():
    exc_type, exc_obj, tb = sys.exc_info()
    f = tb.tb_frame
    lineno = tb.tb_lineno
    filename = f.f_code.co_filename
    linecache.checkcache(filename)
    line = linecache.getline(filename, lineno, f.f_globals)
    print('EXCEPTION IN ({}, LINE {} "{}"): {}'.format(filename, lineno, line.strip(), exc_obj))

if __name__ == "__main__": main()



