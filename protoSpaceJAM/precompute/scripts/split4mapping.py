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
import gc
import math
import errno, os, stat, shutil
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
    parser.add_argument(
        "--part_size",
        default="",
        type=int,
        help="number of gRNAs per part, a good start is 200000",
        metavar="",
    )
    config = parser.parse_args()
    if len(sys.argv) == 1:  # print help message if arguments are not valid
        parser.print_help()
        sys.exit(1)
    return config


logging.setLoggerClass(ColoredLogger)
# logging.basicConfig()
log = logging.getLogger("split files")
log.propagate = False
log.setLevel(logging.INFO)  # set the level of warning displayed

config = vars(parse_args())
part_size = config["part_size"]

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

        starttime = datetime.datetime.now()
        chr_count = 0
        gRNA_count = 0

        # make output directory
        path2dir = config["dir"] + ".split"
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

        # process gRNAs by chr(file)
        new_tab = config["tab"].rstrip("tab") + "split.tab"
        with open(os.path.join(new_tab), "w") as wfh:
            wfh.write(
                "chr_part\tlength\trepeat_length\tgRNA_count\tdirectory\tfilename\n"
            )
            wfh.flush()
            for index, row in df.iterrows():
                chr_count += 1
                current_chr = str(row["Chr"]).strip()
                tab_gRNA_count = int(str(row["gRNA_count"]).strip())
                print(f"processing chr {current_chr} with {tab_gRNA_count} gRNAs")

                # chrs that don't need to be split
                if tab_gRNA_count <= part_size:
                    shutil.copy2(
                        os.path.join(config["dir"], f"{current_chr}.tab.gz"),
                        os.path.join(path2dir, f"{current_chr}.tab.gz"),
                    )
                    row_content = [str(i) for i in row.tolist()]
                    wfh.write("\t".join(row_content))
                    wfh.write(f"\t{path2dir}\t{current_chr}.tab.gz\n")
                    wfh.flush()
                    gRNA_count += tab_gRNA_count
                # chrs that need to be split
                else:
                    with gzip.open(
                        os.path.join(config["dir"], f"{current_chr}.tab.gz"), "rt"
                    ) as gzfh:
                        part = 0
                        file_gRNA_count = 0
                        part_gRNA_count = 0
                        gzwfh = gzip.open(
                            os.path.join(path2dir, f"{current_chr}.part{part}.tab.gz"),
                            "wt",
                        )
                        for line in gzfh:
                            # determine if this is the end of current chunk
                            if (
                                file_gRNA_count % part_size == 0
                                and file_gRNA_count != 0
                            ):
                                gzwfh.close()
                                wfh.write(
                                    f"{current_chr}.part{part}\tNA\tNA\t{part_gRNA_count}\t{path2dir}\t{current_chr}.part{part}.tab.gz\n"
                                )
                                wfh.flush()
                                print(
                                    f"--processed chr {current_chr} part {part} with {part_gRNA_count} gRNAs"
                                )
                                part += 1
                                gzwfh = gzip.open(
                                    os.path.join(
                                        path2dir, f"{current_chr}.part{part}.tab.gz"
                                    ),
                                    "wt",
                                )
                                part_gRNA_count = 0
                            gzwfh.write(line)
                            file_gRNA_count += 1
                            part_gRNA_count += 1
                        # write the last part
                        if file_gRNA_count % part_size != 0 and file_gRNA_count != 0:
                            print(
                                f"--processed chr {current_chr} part {part} with {part_gRNA_count} gRNAs"
                            )
                            wfh.write(
                                f"{current_chr}.part{part}\tNA\tNA\t{part_gRNA_count}\t{path2dir}\t{current_chr}.part{part}.tab.gz\n"
                            )
                            wfh.flush()
                    # check if the gRNA count in split files add up to that in the tab file
                    if file_gRNA_count != tab_gRNA_count:
                        sys.exit(
                            f"error: got {file_gRNA_count} from gz file, and {tab_gRNA_count} from tab file"
                        )
                    gRNA_count += file_gRNA_count
                wfh.flush()

        endtime = datetime.datetime.now()
        elapsed_sec = endtime - starttime
        elapsed_min = elapsed_sec.seconds / 60
        log.info(
            f"finished in {elapsed_min:.2f} min, processed {chr_count} chromosomes, {gRNA_count} gRNAs"
        )
        print(
            f"finished in {elapsed_min:.2f} min, processed {chr_count} chromosomes, {gRNA_count} gRNAs"
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


def BLAST(fastfile):
    numThread2use = max([1, os.cpu_count() - 2])  # for BLAST, use all CPUs except 2
    BLAST_bin, exe_suffix, BLAST_db_path = check_blastDB_human()
    query = fastfile
    # start BLAST
    cmd = [
        f"{BLAST_bin}blastn{exe_suffix}",
        "-task",
        "blastn-short",
        "-query",
        f"{query}",
        "-db",
        f"{BLAST_db_path}",
        "-num_threads",
        f"{numThread2use}",
        "-perc_identity",
        "100",
        "-outfmt",
        "6 qseqid sseqid qstart qend sstart send pident",
        "-out",
        f"{query}.out",
        "-ungapped",
    ]
    p = Popen(cmd, universal_newlines=True)
    p.communicate()  # now wait for the process to finish
    return str(p.returncode)


def handleRemoveReadonly(func, path, exc):
    excvalue = exc[1]
    if func in (os.rmdir, os.remove) and excvalue.errno == errno.EACCES:
        os.chmod(path, stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO)  # 0777
        func(path)
    else:
        raise


if __name__ == "__main__":
    main()
