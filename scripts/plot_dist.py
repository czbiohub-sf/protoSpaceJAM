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
import seaborn as sns
import traceback

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

        #read data
        df = pd.read_csv("../logs/result.csv", low_memory=False)
        df["gstrand"] = [-1 if i>0 else 1 for i in (df["gRNA_start"] - df["gRNA_end"])]

        df_N = df[df["terminus"]=="N"]
        df_C = df[df["terminus"]=="C"]

        #plot params
        bin_num = 33

        #plot histogram
        plot_hist(lst = df_N["distance_between_cut_and_edit"].values,  label = "distance_between_cut_and_edit (N term)", bin_num=bin_num)
        plot_hist(lst = df_C["distance_between_cut_and_edit"].values,  label = "distance_between_cut_and_edit (C term)", bin_num=bin_num)



        #report time used
        endtime = datetime.datetime.now()
        elapsed_sec = endtime - starttime
        elapsed_min = elapsed_sec.seconds / 60

        log.info(f"finished in {elapsed_min:.2f} min ({elapsed_sec} sec)")


    except Exception as e:
        print("Unexpected error:", str(sys.exc_info()))
        print("additional information:", e)
        traceback.print_exc()
        PrintException()

##########################
## function definitions ##
##########################

def subset_df(df, col, max, min):
    df = df[df[col]<=max]
    df = df[df[col]>=min]
    return(df)

def plot_violin(df):
    fig, ax = plt.subplots()
    ax = sns.violinplot(x="phase", y="cfd", hue="terminus",data=df, palette="muted", split=True, cut = 0)
    mpl.rcParams.update({'figure.autolayout': True})
    ax.set_title(f"Distribution of recut cfd scores", fontsize = 13)
    ax.set_xlabel("phase", fontsize = 14)
    ax.set_ylabel('recut cfd', fontsize = 14)
    ax.set_xticklabels(["1","2","3","4"])
    mpl.rcParams["axes.labelsize"] = 14
    mpl.rcParams["axes.titlesize"] = 14
    plt.xticks(fontsize = 14)
    plt.xticks(fontsize = 14)
    plt.show()
    plt.savefig(f"../plots/cfd.violin.png",bbox_inches = "tight")
    plt.close()

def plot_lines(df):
    fig, ax = plt.subplots()
    ax = sns.lineplot(x="phase", y="cfd", hue="terminus", data=df, ci=99)
    mpl.rcParams.update({'figure.autolayout': True})
    ax.set_title(f"recut cfd scores", fontsize = 13)
    ax.set_xlabel("phase", fontsize = 14)
    ax.set_ylabel('recut cfd', fontsize = 14)
    ax.set_xticklabels(["1","2","3","4"])
    mpl.rcParams["axes.labelsize"] = 14
    mpl.rcParams["axes.titlesize"] = 14
    plt.xticks(fontsize = 14)
    plt.xticks(fontsize = 14)
    plt.show()
    plt.savefig(f"../plots/cfd.violin.png",bbox_inches = "tight")
    plt.close()

def plot_lines2(df):
    fig, ax = plt.subplots()
    ax = sns.lineplot(x="phase", y="cfd", hue="terminus", data=df, units="ID",estimator=None,)
    mpl.rcParams.update({'figure.autolayout': True})
    ax.set_title(f"recut cfd scores", fontsize = 13)
    ax.set_xlabel("phase", fontsize = 14)
    ax.set_ylabel('recut cfd', fontsize = 14)
    ax.set_xticklabels(["1","2","3","4"])
    mpl.rcParams["axes.labelsize"] = 14
    mpl.rcParams["axes.titlesize"] = 14
    plt.xticks(fontsize = 14)
    plt.xticks(fontsize = 14)
    plt.show()
    plt.savefig(f"../plots/cfd.violin.png",bbox_inches = "tight")
    plt.close()

def plot_hist(lst, label, bin_num):
    name = label
    d = lst
    d = pd.DataFrame(d)
    data_volume = len(d)
    d.plot.hist(grid=True, bins=bin_num, rwidth=0.9, color='#607c8e', range = (-40,40), legend=None)
    mpl.rcParams.update({'figure.autolayout': True})
    plt.title(f"Distribution of {name}\n for {data_volume} gRNAs", fontsize = 18)
    plt.xlabel(name, fontsize = 16)
    plt.ylabel('Counts', fontsize = 16)
    mpl.rcParams["axes.labelsize"] = 16
    mpl.rcParams["axes.titlesize"] = 16
    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)
    plt.grid(axis='y', alpha=0.75)
    plt.ylim((0,9000))
    plt.xlim((-40,40))
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



