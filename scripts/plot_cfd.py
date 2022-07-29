import os.path
import pandas as pd
import argparse
import sys
import linecache
import datetime
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pickle
from scripts.utils import *
from scipy.stats import gaussian_kde

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

def parse_args():
    parser= MyParser(description='This scripts creates a mapping of ENST to chr from gff3')
    parser.add_argument('--num_to_process', default="all", type=str, help='this parameter decides which file to load, the files have name start_cfd1_[num_to_process].pickle', metavar='')
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
        start_cfd1 = read_pickle_files(os.path.join(f"../pickles/start_cfd1_{num_to_process}_genes.pickle"))
        start_cfd2 = read_pickle_files(os.path.join(f"../pickles/start_cfd2_{num_to_process}_genes.pickle"))
        start_cfd3 = read_pickle_files(os.path.join(f"../pickles/start_cfd3_{num_to_process}_genes.pickle"))
        start_cfd4 = read_pickle_files(os.path.join(f"../pickles/start_cfd4_{num_to_process}_genes.pickle"))
        start_cfdfinal = read_pickle_files(os.path.join(f"../pickles/start_cfdfinal_{num_to_process}_genes.pickle"))

        stop_cfd1 = read_pickle_files(os.path.join(f"../pickles/start_cfd1_{num_to_process}_genes.pickle"))
        stop_cfd2 = read_pickle_files(os.path.join(f"../pickles/start_cfd2_{num_to_process}_genes.pickle"))
        stop_cfd3 = read_pickle_files(os.path.join(f"../pickles/start_cfd3_{num_to_process}_genes.pickle"))
        stop_cfd4 = read_pickle_files(os.path.join(f"../pickles/start_cfd4_{num_to_process}_genes.pickle"))
        #stop_cfdfinal = read_pickle_files(os.path.join(f"../pickles/start_cfdfinal_{num_to_process}_genes.pickle"))

        start_cfd1_df = pd.DataFrame(start_cfd1, columns=['start_cfd1'])
        start_cfd2_df = pd.DataFrame(start_cfd2, columns=['start_cfd2'])
        start_cfd3_df = pd.DataFrame(start_cfd3, columns=['start_cfd3'])
        start_cfd4_df = pd.DataFrame(start_cfd4, columns=['start_cfd4'])
        start_cfdfinal_df = pd.DataFrame(start_cfdfinal, columns=['start_cfd-final'])

        stop_cfd1_df = pd.DataFrame(stop_cfd1, columns=['stop_cfd1'])
        stop_cfd2_df = pd.DataFrame(stop_cfd2, columns=['stop_cfd2'])
        stop_cfd3_df = pd.DataFrame(stop_cfd3, columns=['stop_cfd3'])
        stop_cfd4_df = pd.DataFrame(stop_cfd4, columns=['stop_cfd4'])
        #stop_cfdfinal_df = pd.DataFrame(stop_cfdfinal, columns=['cfd'])

        #prepare df to plot
        #df2plot = best_stop_gRNA_df
        #df2plot = subset_df(df = df2plot, col ="final_weight", max =1, min=0)

        #plot params
        bin_num = 30

        #plot histogram of weight
        plot_hist(df = start_cfd1_df,  col = "start_cfd1", bin_num=bin_num, num_to_process = num_to_process)
        plot_hist(df = start_cfd2_df,  col = "start_cfd2", bin_num=bin_num, num_to_process = num_to_process)
        plot_hist(df = start_cfd3_df,  col = "start_cfd3", bin_num=bin_num, num_to_process = num_to_process)
        plot_hist(df = start_cfd4_df,  col = "start_cfd4", bin_num=bin_num, num_to_process = num_to_process)
        plot_hist(df = start_cfdfinal_df, col="start_cfd-final", bin_num=bin_num, num_to_process=num_to_process)

        #plot scatterplot of weight
        #plot_scatter(df = df2plot, col1 = "dist_weight", col2 = "pos_weight",  num_to_process = num_to_process)
        #plot_scatter(df = df2plot, col1 = "dist_weight", col2 = "spec_weight", num_to_process=num_to_process)
        #plot_scatter(df = df2plot, col1 = "pos_weight", col2 = "spec_weight",  num_to_process=num_to_process)


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

def subset_df(df, col, max, min):
    df = df[df[col]<=max]
    df = df[df[col]>=min]
    return(df)

def plot_hist(df, col, bin_num, num_to_process):
    name = col
    d = df[name]
    data_volume = len(d)
    d.plot.hist(grid=True, bins=bin_num, rwidth=0.9, color='#607c8e', range=[0, 1])
    mpl.rcParams.update({'figure.autolayout': True})
    if name == "final_weight":
        name = "final_score"
    plt.title(f"Distribution of {name} for {data_volume} gRNAs in {num_to_process} transcripts", fontsize = 13)
    plt.xlabel(name, fontsize = 14)
    plt.ylabel('Counts', fontsize = 14)
    mpl.rcParams["axes.labelsize"] = 14
    mpl.rcParams["axes.titlesize"] = 14
    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)
    plt.grid(axis='y', alpha=0.75)
    plt.savefig(f"../plots/{name}.{bin_num}bins.png",bbox_inches = "tight")
    plt.close()

def plot_scatter(df, col1,col2, num_to_process):
    name1 = col1
    name2 = col2
    d1 = df[name1]
    d2 = df[name2]
    data_volume = len(d1)

    mpl.rcParams.update({'figure.autolayout': True})
    mpl.rcParams["axes.labelsize"] = 14
    mpl.rcParams["axes.titlesize"] = 14

    # Calculate the point density
    xy = np.vstack([d1, d2])
    z = gaussian_kde(xy)(xy)
    plt.scatter(d1,d2, c=z, alpha = 1, s=12, marker='o', cmap = 'coolwarm')
    plt.title(f"{name1} vs {name2} for {data_volume} gRNAs in {num_to_process} transcripts", fontsize = 13)
    plt.xlabel(name1)
    plt.ylabel(name2)
    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)
    #plt.grid(axis='y', alpha=0.75)
    cb = plt.colorbar(shrink = 0.5)
    plt.savefig(f"../plots/{name1}.vs.{name2}.png")
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




