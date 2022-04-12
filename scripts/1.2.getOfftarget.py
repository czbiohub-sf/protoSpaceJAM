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
#from BLAST_utils import check_blastDB_human

sys.path.insert(1, '..')
from gRNA_search import *
from utils import *

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

def parse_args():
    parser= MyParser(description='')
    parser.add_argument('--gzfile', default="", type=str, help='path to the gzfile', metavar='')
    parser.add_argument('--gzdir', default="", type=str, help='path to the dir containing gzfile', metavar='')
    parser.add_argument('--gRNA_count', default="", type=int, help='path to the gzfile', metavar='')
    parser.add_argument('--genome_fa', default="", type=str, help='name of genome fasta file', metavar='')
    config = parser.parse_args()
    if len(sys.argv)==1: # print help message if arguments are not valid
        parser.print_help()
        sys.exit(1)
    return config

logging.setLoggerClass(ColoredLogger)
#logging.basicConfig()
log = logging.getLogger("map_gRNA")
log.propagate = False
log.setLevel(logging.INFO) #set the level of warning displayed

config = vars(parse_args())
genome_fa = config["genome_fa"]
pam_len = 3
protospacer_len = 20
n_gRNA_per_chunk = 200
#####################
##      main       ##
#####################
def main():
    try:
        #check input
        if (config["gzfile"] is None or config["gzfile"] == ""):
            log.error(f"please specify --gzfile")
            sys.exit("Please fix the error(s) above and rerun the script")

        starttime = datetime.datetime.now()
        gRNA_count = 0
        #keep a dictionary of multi-mapping guides for current file
        gRNA_match = dict()

        #make output dir
        outdir_path =config["gzdir"] + ".BwaMapped"
        if os.path.isdir(outdir_path):
            #shutil.rmtree(outdir_path)  # remove existing dir
            pass
        else:
            os.makedirs(outdir_path)

        #process gRNA (bwa and parse in chunks)
        total_chunk_count = str(math.ceil(int(str(config["gRNA_count"]).strip())/n_gRNA_per_chunk))
        file = os.path.join(config["gzdir"],config["gzfile"])      
        Chr = file.rstrip('.tab.gz').split("/")[1].split(r'.')[0]
        #print(Chr)
        file_gRNA_count = 0
        n_chunk = str(math.floor(file_gRNA_count / n_gRNA_per_chunk))
        with gzip.open(file, "rt") as fh, gzip.open(os.path.join(outdir_path,f"{config['gzfile'].rstrip('.gz')}.mapped.gz"), "wt") as wgfh:
            wfh = open(f"{file}.{n_chunk}.fa", "w")
            wgfh_tmp = gzip.open(os.path.join(outdir_path,f"{config['gzfile'].rstrip('.gz')}.tmp.gz"), "wt")
            for line in fh:
                #print(line)
                file_gRNA_count +=1

                #read gRNA from gzip file
                seq = line.split("\t")[0]
                pam = line.split("\t")[1]
                st = line.split("\t")[2]
                en = line.split("\t")[3]
                strand = line.split("\t")[4].rstrip()
                gRNA_name = f"{Chr}_{st}_{en}_{strand}_{pam}"

                #write gRNA to fasta file, WITHOUT pam 
                wfh.write(f">{gRNA_name}\n{seq}\n")
                #write gz lines to temp file
                wgfh_tmp.write(line)


                # determine if this is the end of current chunk
                if file_gRNA_count%n_gRNA_per_chunk==0:
                    starttime_chunk = datetime.datetime.now()

                    #find off targets
                    wfh.close() #close the fasta file handle
                    wgfh_tmp.close()
                    fastafile = f"{file}.{n_chunk}.fa"
                    fastafile_dir = os.path.split(fastafile)[0]
                    fastafile_name = os.path.split(fastafile)[1]
                    command = ["python", "../FindOfftargetBwa/findOfftargetsBwa.py", 
                    f"--fa", f"{fastafile_name}",
                    f"--fa_dir", f"{fastafile_dir}", 
                    f"--bin_dir", f"../FindOfftargetBwa/bin/Linux", 
                    f"--script_dir", f"../FindOfftargetBwa/bin/", 
                    f"--bwa_idx", f"../genome_files/indexes_bwa/{genome_fa}", 
                    f"--genome_fa", f"{genome_fa}", 
                    f"--guideLen", f"20"]

                    p = Popen(command)
                    p.communicate()  # wait for the commands to process

                    os.remove(f"{fastafile}")

                    #parse bwa out
                    Bwa_mapping = dict()
                    with open(f"{fastafile}_bwa/bwa.out.bed", "r") as BEDOUT:
                        for line in BEDOUT:
                            #sample line: 1	12629	12652	1_12630_12652_+_AGG|+|0|GCATGCCCTTCCCCAGCATCAGG|0|0|0	1 12.63 Kbp
                            fields = line.split("\t")
                            matchChr=fields[0]
                            matchStart= str(int(fields[1]) + 1)
                            matchEnd=fields[2]
                            matchInfo="|".join(fields[3].split("|")[1:])
                            gRNAname=fields[3].split("|")[0]
                            if gRNAname in Bwa_mapping.keys():
                                Bwa_mapping[gRNAname].append(f"{matchChr}:{matchStart}-{matchEnd}_{matchInfo}")
                            else:
                                Bwa_mapping[gRNAname] = [f"{matchChr}:{matchStart}-{matchEnd}_{matchInfo}"]

                    #write mapping results to a set of new gzip files
                    with gzip.open(os.path.join(outdir_path,f"{config['gzfile'].rstrip('.gz')}.tmp.gz"), "rt") as tmp_gz_fh:
                        for tmp_lines in tmp_gz_fh:
                            #read gRNA from gzip file
                            seq = tmp_lines.split("\t")[0]
                            pam = tmp_lines.split("\t")[1]
                            st = tmp_lines.split("\t")[2]
                            en = tmp_lines.split("\t")[3]
                            strand = tmp_lines.split("\t")[4].rstrip()
                            gRNA_name = f"{Chr}_{st}_{en}_{strand}_{pam}"
                            match_list = ""
                            if gRNA_name in Bwa_mapping.keys():
                                match_list = Bwa_mapping[gRNA_name]
                            wgfh.write(f"{tmp_lines.rstrip()}\t{match_list}\n")
                    
                    #remove bwa working folder
                    shutil.rmtree(f"{fastafile}_bwa")

                    endtime_chunk = datetime.datetime.now()
                    elapsed_sec_chunk = endtime_chunk - starttime_chunk
                    elapsed_min_chunk = elapsed_sec_chunk.seconds / 60

                    print(f"Finished processing chunk {n_chunk} (total {total_chunk_count}), chunk size:{n_gRNA_per_chunk}, took {elapsed_min_chunk:.2f} min")

                    #start new chunk
                    n_chunk = str(math.floor(file_gRNA_count/n_gRNA_per_chunk))
                    wfh = open(f"{file}.{n_chunk}.fa", "w") #open new file handle
                    wgfh_tmp = gzip.open(os.path.join(outdir_path,f"{config['gzfile'].rstrip('.gz')}.tmp.gz"), "wt")
            
            #process last chunk
            if file_gRNA_count%n_gRNA_per_chunk!=0:
                starttime_chunk = datetime.datetime.now()

                #find off targets
                wfh.close() #close the fasta file handle
                wgfh_tmp.close()
                fastafile = f"{file}.{n_chunk}.fa"
                fastafile_dir = os.path.split(fastafile)[0]
                fastafile_name = os.path.split(fastafile)[1]
                command = ["python", "../FindOfftargetBwa/findOfftargetsBwa.py", 
                f"--fa", f"{fastafile_name}",
                f"--fa_dir", f"{fastafile_dir}", 
                f"--bin_dir", f"../FindOfftargetBwa/bin/Linux", 
                f"--script_dir", f"../FindOfftargetBwa/bin/", 
                f"--bwa_idx", f"../genome_files/indexes_bwa/{genome_fa}", 
                f"--genome_fa", f"{genome_fa}", 
                f"--guideLen", f"20"]

                p = Popen(command)
                p.communicate()  # wait for the commands to process

                os.remove(f"{fastafile}")
                #parse bwa out
                Bwa_mapping = dict()
                with open(f"{fastafile}_bwa/bwa.out.bed", "r") as BEDOUT:
                    for line in BEDOUT:
                        #sample line: 1	12629	12652	1_12630_12652_+_AGG|+|0|GCATGCCCTTCCCCAGCATCAGG|0|0|0	1 12.63 Kbp
                        fields = line.split("\t")
                        matchChr=fields[0]
                        matchStart= str(int(fields[1]) + 1)
                        matchEnd=fields[2]
                        matchInfo="|".join(fields[3].split("|")[1:])
                        gRNAname=fields[3].split("|")[0]
                        if gRNAname in Bwa_mapping.keys():
                            Bwa_mapping[gRNAname].append(f"{matchChr}:{matchStart}-{matchEnd}_{matchInfo}")
                        else:
                            Bwa_mapping[gRNAname] = [f"{matchChr}:{matchStart}-{matchEnd}_{matchInfo}"]
                #write mapping results to a set of new gzip files
                with gzip.open(os.path.join(outdir_path,f"{config['gzfile'].rstrip('.gz')}.tmp.gz"), "rt") as tmp_gz_fh:
                    for tmp_lines in tmp_gz_fh:
                        #read gRNA from gzip file
                        seq = tmp_lines.split("\t")[0]
                        pam = tmp_lines.split("\t")[1]
                        st = tmp_lines.split("\t")[2]
                        en = tmp_lines.split("\t")[3]
                        strand = tmp_lines.split("\t")[4].rstrip()
                        gRNA_name = f"{Chr}_{st}_{en}_{strand}_{pam}"
                        match_list = ""
                        if gRNA_name in Bwa_mapping.keys():
                            match_list = Bwa_mapping[gRNA_name]
                        wgfh.write(f"{tmp_lines.rstrip()}\t{match_list}\n")

                #remove bwa working folder
                shutil.rmtree(f"{fastafile}_bwa")

                endtime_chunk = datetime.datetime.now()
                elapsed_sec_chunk = endtime_chunk - starttime_chunk
                elapsed_min_chunk = elapsed_sec_chunk.seconds / 60

                print(f"Finished processing chunk {n_chunk} (total {total_chunk_count}), chunk size:{n_gRNA_per_chunk}, took {elapsed_min_chunk:.2f} min")

                #remove tmp gz file
                os.remove(os.path.join(outdir_path,f"{config['gzfile'].rstrip('.gz')}.tmp.gz"))

            gRNA_count += file_gRNA_count

        endtime = datetime.datetime.now()
        elapsed_sec = endtime - starttime
        elapsed_min = elapsed_sec.seconds / 60
        log.info(f"finished in {elapsed_min:.2f} min, bwa-mapped {gRNA_count} gRNAs in file {file}")
        print(f"finished in {elapsed_min:.2f} min, bwa-mapped {gRNA_count} gRNAs in file {file}")

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
    print('EXCEPTION IN ({}, LINE {} "{}"): {}'.format(filename, lineno, line.strip(), exc_obj))

def BLAST(fastfile):
    numThread2use = max([1, os.cpu_count() - 2])  # for BLAST, use all CPUs except 2
    #numThread2use = 1
    BLAST_bin, exe_suffix, BLAST_db_path = check_blastDB_human()
    query = fastfile
    # start BLAST
    cmd = [f"{BLAST_bin}blastn{exe_suffix}", "-task", "blastn-short", "-query", f"{query}", "-db", f"{BLAST_db_path}",
           "-num_threads", f"{numThread2use}", "-perc_identity", "100", "-outfmt",
           "6 qseqid sseqid qstart qend sstart send pident", "-out", f"{query}.out", "-ungapped"]
    p = Popen(cmd, universal_newlines=True)
    p.communicate()  # now wait for the process to finish
    return str(p.returncode)


if __name__ == "__main__": main()

