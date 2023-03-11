import os.path
import sys
import numpy as np
import errno, os, stat

sys.path.insert(1, "../..")
from utils import *


def parse_args():
    parser = MyParser(
        description="This script splits the gz files into parts suitable for parallel computing, please specify either --fastagz or --fasta"
    )
    parser.add_argument(
        "--dir",
        default="",
        type=str,
        help="path to the directory containing tabulated gRNA outputs",
        metavar="",
    )
    parser.add_argument(
        "--tab",
        default="",
        type=str,
        help="path to the tab file (containing gRNA info)",
        metavar="",
    )
    config = parser.parse_args()
    if len(sys.argv) == 1:  # print help message if arguments are not valid
        parser.print_help()
        sys.exit(1)
    return config


logging.setLoggerClass(ColoredLogger)
# logging.basicConfig()
log = logging.getLogger("map_gRNA")
log.propagate = False
log.setLevel(logging.INFO)  # set the level of warning displayed

config = vars(parse_args())
part_size = 200000

#####################
##      main       ##
#####################
def main():
    try:
        # check input
        if (config["dir"] is None or config["dir"] == "") or (
            config["tab"] is None or config["tab"] == ""
        ):
            log.error(f"please specify both --dir or --tab")
            sys.exit("Please fix the error(s) above and rerun the script")

        # make output directory
        path2dir = config["dir"] + ".circos"
        if os.path.isdir(path2dir):
            shutil.rmtree(path2dir, onerror=handleRemoveReadonly)  # remove existing dir
        # if not os.path.isdir(path2dir):
        os.makedirs(path2dir)

        # read chr information
        df = pd.read_csv(os.path.join(config["tab"]), sep="\t")
        if (not "Chr" in df.columns) or (not "gRNA_count" in df.columns):
            log.error(
                f"The csv file does not contain either a column named: Chr or gRNA_count"
            )
            sys.exit("Please fix the error(s) above and rerun the script")

        for bin_size in [25000, 50000, 100000, 250000, 500000, 1000000]:
            print(f"bin_size={str(int(bin_size / 1000))}kb")
            bin_gRNA_num(bin_size=bin_size, df=df, outdir=path2dir)
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


def bin_gRNA_num(bin_size, df, outdir):
    starttime = datetime.datetime.now()
    chr_count = 0
    gRNA_count = 0
    out = f"circos.{str(int(bin_size / 1000))}kb.tab"
    out = os.path.join(outdir, out)
    with open(out, "w") as wfh:
        for index, row in df.iterrows():
            current_chr = str(row["Chr"]).strip()
            if len(current_chr) > 3:  # filter out contigs
                continue
            chr_count += 1
            current_chr = str(row["Chr"]).strip()
            current_chr_len = str(row["length"]).strip()
            tab_gRNA_count = int(str(row["gRNA_count"]).strip())
            # print(f"processing chr {current_chr} (length: {current_chr_len}) with {tab_gRNA_count} gRNAs")
            gRNA_count += tab_gRNA_count

            # get the bins
            bins = list(range(0, int(current_chr_len), bin_size))

            # start going through gRNAs
            list_gRNA_pos = []
            with gzip.open(
                os.path.join(config["dir"], f"{current_chr}.tab.gz"), "rt"
            ) as gzfh:
                for line in gzfh:
                    st = int(line.split("\t")[2])
                    en = int(line.split("\t")[3])
                    gRNA_pos = (st + en) / 2
                    list_gRNA_pos.append(gRNA_pos)

            # put count num gRNA in bins
            n, bins = np.histogram(list_gRNA_pos, bins=bins)
            # print(bins[0:10])
            # print(n[0:10])
            # print(f"{len(bins)} {len(n)}")
            # write to file
            for i in range(len(bins)):
                if not i == (len(bins) - 1):
                    wfh.write(f"hs{current_chr}\t{bins[i]}\t{bins[i + 1]}\t{n[i]}\n")

    endtime = datetime.datetime.now()
    elapsed_sec = endtime - starttime
    elapsed_min = elapsed_sec.seconds / 60
    log.info(
        f"bin_size={str(int(bin_size / 1000))}kb, finished in {elapsed_min:.2f} min, processed {chr_count} chromosomes, {gRNA_count} gRNAs"
    )
    print(
        f"bin_size={str(int(bin_size / 1000))}kb, finished in {elapsed_min:.2f} min, processed {chr_count} chromosomes, {gRNA_count} gRNAs"
    )


def handleRemoveReadonly(func, path, exc):
    excvalue = exc[1]
    if func in (os.rmdir, os.remove) and excvalue.errno == errno.EACCES:
        os.chmod(path, stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO)  # 0777
        func(path)
    else:
        raise


if __name__ == "__main__":
    main()
