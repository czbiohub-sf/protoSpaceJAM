import os.path
import pandas as pd
from Bio import SeqIO
import csv
import sys
import linecache
import datetime
import gzip
import shutil
import gc
import re

from utils import *
from gRNA_search import *

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
    parser.add_argument(
        "--fasta",
        default="",
        type=str,
        help="path to the human genome fasta.gz file",
        metavar="",
    )
    parser.add_argument(
        "--outdir",
        default="",
        type=str,
        help="path to the human genome fasta.gz file",
        metavar="",
    )
    config = parser.parse_args()
    if len(sys.argv) == 1:  # print help message if arguments are not valid
        parser.print_help()
        sys.exit(1)
    return config


protosp_len = 20
PAM = "NGG"

logging.setLoggerClass(ColoredLogger)
# logging.basicConfig()
log = logging.getLogger("Scan gRNA")
log.propagate = False
log.setLevel(logging.INFO)  # set the level of warning displayed

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
        chr_count = 0
        gRNA_count = 0

        infile = config["fastagz"]
        if infile is None or infile == "":
            infile = config["fasta"]

        N = re.compile("N")

        # make output dir
        outdir = config["outdir"]
        outdir_path = os.path.join(outdir, f"gRNA.tab.gz")
        if os.path.isdir(outdir_path):
            shutil.rmtree(outdir_path)  # remove existing dir
        os.makedirs(outdir_path)

        infile_basename = os.path.basename(infile)
        with gzip.open(infile, "rt") as f, open(
            os.path.join(outdir, f"{infile_basename}.out.tab"), "w"
        ) as wfh:
            wfh.write(f"Chr\tlength\trepeat_length\tgRNA_count\n")  # header
            for record in SeqIO.parse(f, "fasta"):
                seq = str(record.seq)

                # get gRNAs
                res_gRNA_list = search_gRNA(
                    protosp_len=protosp_len, PAM=PAM, search_in=seq
                )

                # process gRNAs
                chr_count += 1
                print(
                    f"chr:\t{record.id}\tlength:\t{len(record.seq)}\t{len(res_gRNA_list)}",
                    flush=True,
                )

                # write gRNA to file
                chr_gRNA_count = 0
                gzip_path = os.path.join(outdir_path, f"{record.id}.tab.gz")
                with gzip.open(gzip_path, "wt") as gwfh:
                    for res_gRNA in res_gRNA_list:
                        if re.search(N, res_gRNA.protospacer) is not None:
                            continue
                        gwfh.write(
                            f"{res_gRNA.protospacer}\t{res_gRNA.pam}\t{res_gRNA.g_st}\t{res_gRNA.g_en}\t{res_gRNA.g_strand}\t{res_gRNA.protospacer5p_flank}\t{res_gRNA.pam3p_flank}\n"
                        )
                        gRNA_count += 1
                        chr_gRNA_count += 1

                # lowercase count
                lowercase_characters = re.findall(r"[a-z]", seq)
                lc_count = len(lowercase_characters)

                # write tab file
                wfh.write(
                    f"{record.id}\t{len(record.seq)}\t{lc_count}\t{chr_gRNA_count}\n"
                )
                wfh.flush()

        endtime = datetime.datetime.now()
        elapsed_sec = endtime - starttime
        elapsed_min = elapsed_sec.seconds / 60
        log.info(
            f"finished in {elapsed_min:.2f} min, processed {chr_count} chromosomes, found {gRNA_count} gRNAs"
        )
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
