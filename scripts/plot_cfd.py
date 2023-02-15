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
        sys.stderr.write("error: %s\n" % message)
        self.print_help()
        sys.exit(2)


def parse_args():
    parser = MyParser(
        description="This scripts creates a mapping of ENST to chr from gff3"
    )
    parser.add_argument(
        "--num_to_process",
        default="all",
        type=str,
        help="this parameter decides which file to load, the files have name start_cfd1_[num_to_process].pickle",
        metavar="",
    )
    config = parser.parse_args()
    return config


logging.setLoggerClass(ColoredLogger)
# logging.basicConfig()
log = logging.getLogger("plot weights")
log.propagate = False
log.setLevel(logging.INFO)  # set the level of warning displayed
# log.setLevel(logging.DEBUG) #set the level of warning displayed

config = vars(parse_args())

#####################
##      main       ##
#####################
def main():
    try:
        starttime = datetime.datetime.now()

        # read data
        df_N = pd.read_csv("../logs/out_Nterm_recut_cfd.csv", low_memory=False)
        df_C = pd.read_csv("../logs/out_Cterm_recut_cfd.csv", low_memory=False)
        df_N["terminus"] = "N"

        df_N_melt = pd.melt(
            df_N,
            id_vars=["ID"],
            value_vars=["cfd1", "cfd2", "cfd3", "cfd4"],
            value_name="cfd",
            var_name="phase",
        )
        df_N_melt["terminus"] = "N"
        df_C_melt = pd.melt(
            df_C,
            id_vars=["ID"],
            value_vars=["cfd1", "cfd2", "cfd3", "cfd4"],
            value_name="cfd",
            var_name="phase",
        )
        df_C_melt["terminus"] = "C"

        df = pd.concat([df_N_melt, df_C_melt])
        plot_violin(df)

        df_N_melt["ID"] = df_N_melt["ID"] + "_N"
        df_C_melt["ID"] = df_C_melt["ID"] + "_C"
        df = pd.concat([df_N_melt, df_C_melt], ignore_index=True)
        plot_lines(df)
        # plot_lines2(df)

        cfd1_N = df_N["cfd1"].values
        cfd2_N = df_N["cfd2"].values
        cfd3_N = df_N["cfd3"].values
        cfd4_N = df_N["cfd4"].values
        #
        cfd1_C = df_C["cfd1"].values
        cfd2_C = df_C["cfd2"].values
        cfd3_C = df_C["cfd3"].values
        cfd4_C = df_C["cfd4"].values
        #
        cfd1 = np.concatenate((cfd1_N, cfd1_C), axis=None)
        cfd2 = np.concatenate((cfd2_N, cfd2_C), axis=None)
        cfd3 = np.concatenate((cfd3_N, cfd3_C), axis=None)
        cfd4 = np.concatenate((cfd4_N, cfd4_C), axis=None)

        # plot params
        bin_num = 33
        num_to_process = len(cfd1)

        # plot histogram
        plot_hist(lst=cfd1, label="cfd1", bin_num=bin_num)
        plot_hist(lst=cfd2, label="cfd2", bin_num=bin_num)
        plot_hist(lst=cfd3, label="cfd3", bin_num=bin_num)
        plot_hist(lst=cfd4, label="cfd4", bin_num=bin_num)

        plot_hist(lst=cfd1_N, label="cfd1_N", bin_num=bin_num)
        plot_hist(lst=cfd2_N, label="cfd2_N", bin_num=bin_num)
        plot_hist(lst=cfd3_N, label="cfd3_N", bin_num=bin_num)
        plot_hist(lst=cfd4_N, label="cfd4_N", bin_num=bin_num)
        plot_hist(lst=cfd1_C, label="cfd1_C", bin_num=bin_num)
        plot_hist(lst=cfd2_C, label="cfd2_C", bin_num=bin_num)
        plot_hist(lst=cfd3_C, label="cfd3_C", bin_num=bin_num)
        plot_hist(lst=cfd4_C, label="cfd4_C", bin_num=bin_num)
        # print stats
        print(f"number of cfd<0.03 in cfd1: {sum(i < 0.03 for i in cfd1)}")
        print(f"number of cfd<0.03 in cfd2: {sum(i < 0.03 for i in cfd2)}")
        print(f"number of cfd<0.03 in cfd3: {sum(i < 0.03 for i in cfd3)}")
        print(f"number of cfd<0.03 in cfd4: {sum(i < 0.03 for i in cfd4)}")

        print(f"number of cfd<0.03 in cfd1_N: {sum(i < 0.03 for i in cfd1_N)}")
        print(f"number of cfd<0.03 in cfd2_N: {sum(i < 0.03 for i in cfd2_N)}")
        print(f"number of cfd<0.03 in cfd3_N: {sum(i < 0.03 for i in cfd3_N)}")
        print(f"number of cfd<0.03 in cfd4_N: {sum(i < 0.03 for i in cfd4_N)}")

        print(f"number of cfd<0.03 in cfd1_C: {sum(i < 0.03 for i in cfd1_C)}")
        print(f"number of cfd<0.03 in cfd2_C: {sum(i < 0.03 for i in cfd2_C)}")
        print(f"number of cfd<0.03 in cfd3_C: {sum(i < 0.03 for i in cfd3_C)}")
        print(f"number of cfd<0.03 in cfd4_C: {sum(i < 0.03 for i in cfd4_C)}")
        # plot scatterplot of weight
        # plot_scatter(df = df2plot, col1 = "dist_weight", col2 = "pos_weight",  num_to_process = num_to_process)
        # plot_scatter(df = df2plot, col1 = "dist_weight", col2 = "spec_weight", num_to_process=num_to_process)
        # plot_scatter(df = df2plot, col1 = "pos_weight", col2 = "spec_weight",  num_to_process=num_to_process)

        # plot dist_weight (with overflow bin)
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

        # report time used
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
    df = df[df[col] <= max]
    df = df[df[col] >= min]
    return df


def plot_violin(df):
    fig, ax = plt.subplots()
    ax = sns.violinplot(
        x="phase", y="cfd", hue="terminus", data=df, palette="muted", split=True, cut=0
    )
    mpl.rcParams.update({"figure.autolayout": True})
    ax.set_title(f"Distribution of recut cfd scores", fontsize=13)
    ax.set_xlabel("phase", fontsize=14)
    ax.set_ylabel("recut cfd", fontsize=14)
    ax.set_xticklabels(["1", "2", "3", "4"])
    mpl.rcParams["axes.labelsize"] = 14
    mpl.rcParams["axes.titlesize"] = 14
    plt.xticks(fontsize=14)
    plt.xticks(fontsize=14)
    plt.show()
    plt.savefig(f"../plots/cfd.violin.png", bbox_inches="tight")
    plt.close()


def plot_lines(df):
    fig, ax = plt.subplots()
    ax = sns.lineplot(x="phase", y="cfd", hue="terminus", data=df, ci=99)
    mpl.rcParams.update({"figure.autolayout": True})
    ax.set_title(f"recut cfd scores", fontsize=13)
    ax.set_xlabel("phase", fontsize=14)
    ax.set_ylabel("recut cfd", fontsize=14)
    ax.set_xticklabels(["1", "2", "3", "4"])
    mpl.rcParams["axes.labelsize"] = 14
    mpl.rcParams["axes.titlesize"] = 14
    plt.xticks(fontsize=14)
    plt.xticks(fontsize=14)
    plt.show()
    plt.savefig(f"../plots/cfd.violin.png", bbox_inches="tight")
    plt.close()


def plot_lines2(df):
    fig, ax = plt.subplots()
    ax = sns.lineplot(
        x="phase", y="cfd", hue="terminus", data=df, units="ID", estimator=None,
    )
    mpl.rcParams.update({"figure.autolayout": True})
    ax.set_title(f"recut cfd scores", fontsize=13)
    ax.set_xlabel("phase", fontsize=14)
    ax.set_ylabel("recut cfd", fontsize=14)
    ax.set_xticklabels(["1", "2", "3", "4"])
    mpl.rcParams["axes.labelsize"] = 14
    mpl.rcParams["axes.titlesize"] = 14
    plt.xticks(fontsize=14)
    plt.xticks(fontsize=14)
    plt.show()
    plt.savefig(f"../plots/cfd.violin.png", bbox_inches="tight")
    plt.close()


def plot_hist(lst, label, bin_num):
    name = label
    d = lst
    d = pd.DataFrame(d)
    data_volume = len(d)
    d.plot.hist(
        grid=True, bins=bin_num, rwidth=0.9, color="#607c8e", range=[0, 1], legend=None
    )
    mpl.rcParams.update({"figure.autolayout": True})
    plt.title(f"Distribution of {name} for {data_volume} gRNAs", fontsize=18)
    plt.xlabel(name, fontsize=16)
    plt.ylabel("Counts", fontsize=16)
    mpl.rcParams["axes.labelsize"] = 16
    mpl.rcParams["axes.titlesize"] = 16
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.grid(axis="y", alpha=0.75)
    # plt.ylim((0,2000))
    plt.savefig(f"../plots/{name}.{bin_num}bins.png", bbox_inches="tight")
    plt.close()


def plot_scatter(df, col1, col2, num_to_process):
    name1 = col1
    name2 = col2
    d1 = df[name1]
    d2 = df[name2]
    data_volume = len(d1)

    mpl.rcParams.update({"figure.autolayout": True})
    mpl.rcParams["axes.labelsize"] = 14
    mpl.rcParams["axes.titlesize"] = 14

    # Calculate the point density
    xy = np.vstack([d1, d2])
    z = gaussian_kde(xy)(xy)
    plt.scatter(d1, d2, c=z, alpha=1, s=12, marker="o", cmap="coolwarm")
    plt.title(
        f"{name1} vs {name2} for {data_volume} gRNAs in {num_to_process} transcripts",
        fontsize=13,
    )
    plt.xlabel(name1)
    plt.ylabel(name2)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    # plt.grid(axis='y', alpha=0.75)
    cb = plt.colorbar(shrink=0.5)
    plt.savefig(f"../plots/{name1}.vs.{name2}.png")
    plt.close()


def cal_elapsed_time(starttime, endtime):
    """
    output: [elapsed_min,elapsed_sec]
    """
    elapsed_sec = endtime - starttime
    elapsed_min = elapsed_sec.seconds / 60
    return [elapsed_min, elapsed_sec]


def read_pickle_files(file):
    if os.path.isfile(file):
        with open(file, "rb") as handle:
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
    print(
        'EXCEPTION IN ({}, LINE {} "{}"): {}'.format(
            filename, lineno, line.strip(), exc_obj
        )
    )


if __name__ == "__main__":
    main()
