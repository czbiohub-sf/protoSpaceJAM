import os.path
import pandas as pd
from subprocess import Popen
from Bio import SeqIO
from Bio.Seq import reverse_complement
import csv
import argparse
import sys
import linecache
import datetime
import gzip
import shutil
import gc
import math
import re
import pickle
from scripts.utils import *
import traceback

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

def parse_args():
    parser= MyParser(description='This scripts creates a mapping of ENST to chr from gff3')
    parser.add_argument('--genome_ver', default="hg38", type=str, help='pickle file containing the ENST_info dict', metavar='')
    parser.add_argument('--path2csv', default="input/mart_export_canonical_proteincoding.csv", type=str,
                        help='path to a csv file containing ENST information\n *required columns*: Ensemble_ID',
                        metavar='')
    config = parser.parse_args()
    return config

logging.setLoggerClass(ColoredLogger)
#logging.basicConfig()
log = logging.getLogger("ProtospaceX")
log.propagate = False
log.setLevel(logging.INFO) #set the level of warning displayed
#log.setLevel(logging.DEBUG) #set the level of warning displayed

max_cut2ins_dist = 50
HDR_arm_len = 100

mNG2_11 = "ACCGAGCTCAACTTCAAGGAGTGGCAAAAGGCCTTTACCGATATGATG"
tag = mNG2_11

config = vars(parse_args())

#####################
##      main       ##
#####################
def main():
    try:
        starttime = datetime.datetime.now()
        freq_dict = dict()

        #load gRNA info index (mapping of chromosomal location to file parts)
        log.info("loading the mapping of chromosomal location to (gRNA) file parts")
        loc2file_index = read_pickle_files(os.path.join(f"gRNA_{config['genome_ver']}","loc2file_index.pickle"))

        #load chr location to type (e.g. UTR, cds, exon/intron junction) mappings
        log.info("loading the mapping of chromosomal location to type (e.g. UTR, cds, exon/intron junction)")
        loc2posType = read_pickle_files(os.path.join("genome_files","gff3",config['genome_ver'],"loc2posType.pickle"))

        #load gene model info
        log.info("loading gene model info")
        ENST_info = read_pickle_files(os.path.join("genome_files","gff3",config['genome_ver'],"ENST_info.pickle"))

        #load codon phase info
        log.info("loading codon phase info")
        ENST_PhaseInCodon = read_pickle_files(os.path.join("genome_files","gff3",config['genome_ver'],"ENST_codonPhase.pickle"))

        #report time used
        elapsed = cal_elapsed_time(starttime,datetime.datetime.now())
        log.info(f"finished loading in {elapsed[0]:.2f} min ({elapsed[1]} sec)")

        #report the number of ENSTs which has ATG at the end of the exon
        ExonEnd_ATG_count,ExonEnd_ATG_list = count_ATG_at_exonEnd(ENST_info)

        #open log files
        recut_CFD_pass = open("logs/recut_CFD_pass.txt", "w")
        recut_CFD_fail = open("logs/recut_CFD_fail.txt", "w")


        #load ENST list (the user input list or the whole transcriptome)
        if os.path.isfile(config['path2csv']):
            log.info(f"begin processing user-supplied list of gene IDs in file {config['path2csv']}")
            df = pd.read_csv(os.path.join(config['path2csv']))
            # check csv columns
            keys2check = set(["Ensemble_ID"])
            if not keys2check.issubset(df.columns):
                log.error(f"Missing columns in the input csv file\n Required columns:\"Ensemble_ID\"")
                log.info(f"Please fix the input csv file and try again")
                sys.exit()
        else:
            log.warning(f"The input file {config['path2csv']} is not found, using the whole human transcriptome")
            input("Press Enter to continue...")
            df = pd.DataFrame(ENST_info.keys(), columns = ["Ensemble_ID"]) # create data frame from ENST_info

        #dataframes to store best gRNAs
        best_start_gRNAs = pd.DataFrame()
        best_stop_gRNAs = pd.DataFrame()

        #list to store IDs failed to yield gRNAs
        start_failed = []
        stop_failed = []

        #list to store IDs whose gRNA is outside of the HDR window
        #gRNA_out_of_arms = dict()
        #gRNA_out_of_arms["start"]=dict()
        #gRNA_out_of_arms["stop"]=dict()

        #list to store CFD scores
        start_cfd1 = []
        start_cfd2 = []
        start_cfd3 = []
        start_cfd4 = []
        start_cfdfinal = []
        stop_cfd1 = []
        stop_cfd2 = []
        stop_cfd3 = []
        stop_cfd4 = []
        stop_cfdfinal = []

        #loop through each ENST
        transcript_count = 0
        protein_coding_transcripts_count = 0
        df_start_gRNAs = pd.DataFrame()
        df_stop_gRNAs = pd.DataFrame()
        for index, row in df.iterrows():
            ENST_ID = row["Ensemble_ID"]
            if not ENST_ID in ENST_info.keys():
                log.warning(f"skipping {ENST_ID} b/c transcript is not in the annotated ENST collection (excluding those on chr_patch_hapl_scaff)")
                continue
            transcript_type = ENST_info[ENST_ID].description.split("|")[1]
            if transcript_type == "protein_coding" and ENST_ID == "ENST00000372781":
                # if not ENST_ID in ExonEnd_ATG_list: # only process edge cases in which genes with ATG are at the end of exons
                #     continue
                log.info(f"processing {ENST_ID}\ttranscript type: {transcript_type}")
                ranked_df_gRNAs_ATG, ranked_df_gRNAs_stop = get_gRNAs(ENST_ID = ENST_ID, ENST_info= ENST_info, freq_dict = freq_dict, loc2file_index= loc2file_index, loc2posType = loc2posType, dist = max_cut2ins_dist, genome_ver=config["genome_ver"])
                if ranked_df_gRNAs_ATG.empty == False:
                    df_start_gRNAs = pd.concat([df_start_gRNAs,ranked_df_gRNAs_ATG])
                else:
                    start_failed.append(ENST_ID)
                if ranked_df_gRNAs_stop.empty == False:
                    df_stop_gRNAs = pd.concat([df_stop_gRNAs,ranked_df_gRNAs_stop])
                else:
                    stop_failed.append(ENST_ID)
                #get best gRNA
                best_start_gRNA = ranked_df_gRNAs_ATG[ranked_df_gRNAs_ATG["final_weight"] == ranked_df_gRNAs_ATG['final_weight'].max()]
                best_stop_gRNA = ranked_df_gRNAs_stop[ranked_df_gRNAs_stop["final_weight"] == ranked_df_gRNAs_stop['final_weight'].max()]
                if best_start_gRNA.empty == False:
                    if best_start_gRNA.shape[0] > 1: # multiple best scoring gRNA
                        best_start_gRNA = best_start_gRNA[best_start_gRNA["CSS"] == best_start_gRNA["CSS"].max()] # break the tie by CSS score
                        best_start_gRNA = best_start_gRNA.head(1) #get the first row in case of ties

                    #get HDR template
                    HDR_template = get_HDR_template(df = best_start_gRNA, ENST_info = ENST_info, type = "start", ENST_PhaseInCodon = ENST_PhaseInCodon, HDR_arm_len=HDR_arm_len, genome_ver=config["genome_ver"], tag = tag)
                    #append cfd score to list for plotting
                    if hasattr(HDR_template,"cdf_score_post_mut_ins"):
                        start_cfd1.append(HDR_template.cdf_score_post_mut_ins)
                    if hasattr(HDR_template,"cdf_score_post_mut2"):
                        start_cfd2.append(HDR_template.cdf_score_post_mut2)
                    if hasattr(HDR_template,"cdf_score_post_mut3"):
                        start_cfd3.append(HDR_template.cdf_score_post_mut3)
                    if hasattr(HDR_template,"cdf_score_post_mut4"):
                        start_cfd4.append(HDR_template.cdf_score_post_mut4)
                    if hasattr(HDR_template,"final_cfd"):
                        start_cfdfinal.append(HDR_template.final_cfd)
                    #write log
                    this_log = f"{HDR_template.info}{HDR_template.info_arm}{HDR_template.info_p1}{HDR_template.info_p2}{HDR_template.info_p3}{HDR_template.info_p4}final CFD:{HDR_template.final_cfd:.4f}\nbefore mutation: {HDR_template.ODN_vanillia}\n  after mutation:{HDR_template.ODN_postMut}\n     final ssODN:{HDR_template.ODN_postMut_ss}\n"
                    if HDR_template.final_cfd < 0.03:
                        recut_CFD_pass.write(this_log)
                    else:
                        recut_CFD_fail.write(this_log)

                    best_start_gRNAs = pd.concat([best_start_gRNAs, best_start_gRNA]) #append the best gRNA to the final df
                if best_stop_gRNA.empty == False:
                    if best_stop_gRNA.shape[0] > 1: # multiple best scoring gRNA
                        best_stop_gRNA = best_stop_gRNA[best_stop_gRNA["CSS"] == best_stop_gRNA["CSS"].max()] # break the tie by CSS score
                        best_stop_gRNA = best_stop_gRNA.head(1) #get the first row in case of ties

                    #get HDR template
                    HDR_template = get_HDR_template(df=best_stop_gRNA, ENST_info=ENST_info, type="stop", ENST_PhaseInCodon = ENST_PhaseInCodon, HDR_arm_len = HDR_arm_len, genome_ver=config["genome_ver"], tag = tag)
                    # append cfd score to list for plotting
                    if hasattr(HDR_template,"cdf_score_post_mut_ins"):
                        stop_cfd1.append(HDR_template.cdf_score_post_mut_ins)
                    if hasattr(HDR_template,"cdf_score_post_mut2"):
                        stop_cfd2.append(HDR_template.cdf_score_post_mut2)
                    if hasattr(HDR_template,"cdf_score_post_mut3"):
                        stop_cfd3.append(HDR_template.cdf_score_post_mut3)
                    if hasattr(HDR_template,"cdf_score_post_mut4"):
                        stop_cfd4.append(HDR_template.cdf_score_post_mut4)
                    if hasattr(HDR_template,"final_cfd"):
                        stop_cfdfinal.append(HDR_template.final_cfd)
                    #write log
                    this_log = f"{HDR_template.info}{HDR_template.info_arm}{HDR_template.info_p1}{HDR_template.info_p2}{HDR_template.info_p3}{HDR_template.info_p4}final CFD:{HDR_template.final_cfd:.4f}\nbefore mutation: {HDR_template.ODN_vanillia}\n  after mutation:{HDR_template.ODN_postMut}\n     final ssODN:{HDR_template.ODN_postMut_ss}\n"
                    if HDR_template.final_cfd < 0.03:
                        recut_CFD_pass.write(this_log)
                    else:
                        recut_CFD_fail.write(this_log)

                    best_stop_gRNAs = pd.concat([best_stop_gRNAs, best_stop_gRNA])
                protein_coding_transcripts_count +=1
            else:
                log.info(f"skipping {ENST_ID} transcript type: {transcript_type} b/c transcript is not protein_coding")
            transcript_count +=1
            #report progress
            if protein_coding_transcripts_count%10 == 0 and protein_coding_transcripts_count != 0:
                endtime = datetime.datetime.now()
                elapsed_sec = endtime - starttime
                elapsed_min = elapsed_sec.seconds / 60
                log.info(f"processed {protein_coding_transcripts_count}/{transcript_count} transcripts, elapsed time {elapsed_min:.2f} min ({elapsed_sec} sec)")
                #gnum = len(gRNA_out_of_arms["start"]) + len(gRNA_out_of_arms["stop"])
                #log.info(f"number of gRNAs outside the HDR arm: {gnum}")

            ################
            #early stopping#
            ################
            # if ENST_ID == "ENST00000360426":
            #    sys.exit()
            #num_to_process = 2000
            # if protein_coding_transcripts_count >=num_to_process:
            #     break

        #write csv out
        endtime = datetime.datetime.now()
        elapsed_sec = endtime - starttime
        elapsed_min = elapsed_sec.seconds / 60

        if 'num_to_process' in locals():
            pass
        else:
            num_to_process = "all"

        # write gRNA dfs to file
        with open(f"pickles/start_gRNAs_of_{num_to_process}_genes.pickle", 'wb') as handle:
            pickle.dump(df_start_gRNAs, handle, protocol=pickle.HIGHEST_PROTOCOL)
        with open(f"pickles/stop_gRNAs_of_{num_to_process}_genes.pickle", 'wb') as handle:
            pickle.dump(df_stop_gRNAs, handle, protocol=pickle.HIGHEST_PROTOCOL)

        # write best gRNA dfs to file
        with open(f"pickles/best_start_gRNAs_of_{num_to_process}_genes.pickle", 'wb') as handle:
            pickle.dump(best_start_gRNAs, handle, protocol=pickle.HIGHEST_PROTOCOL)
        with open(f"pickles/best_stop_gRNAs_of_{num_to_process}_genes.pickle", 'wb') as handle:
            pickle.dump(best_stop_gRNAs, handle, protocol=pickle.HIGHEST_PROTOCOL)

        # write failed ENSTs to file
        with open(f"pickles/start_failed_IDs_of_{num_to_process}_genes.pickle", 'wb') as handle:
            pickle.dump(start_failed, handle, protocol=pickle.HIGHEST_PROTOCOL)
        with open(f"pickles/stop_failed_IDs_of_{num_to_process}_genes.pickle", 'wb') as handle:
            pickle.dump(stop_failed, handle, protocol=pickle.HIGHEST_PROTOCOL)

        # write ENSTs (whose gRNA is outside of the default HDR arm) to file
        #with open(f"pickles/gRNA_out_of_arms_{num_to_process}_genes.pickle", 'wb') as handle:
        #    pickle.dump(gRNA_out_of_arms, handle, protocol=pickle.HIGHEST_PROTOCOL)

        # write lists of cfd to file
        with open(f"pickles/start_cfd1_{num_to_process}_genes.pickle", 'wb') as handle:
            pickle.dump(start_cfd1, handle, protocol=pickle.HIGHEST_PROTOCOL)
        with open(f"pickles/start_cfd2_{num_to_process}_genes.pickle", 'wb') as handle:
            pickle.dump(start_cfd2, handle, protocol=pickle.HIGHEST_PROTOCOL)
        with open(f"pickles/start_cfd3_{num_to_process}_genes.pickle", 'wb') as handle:
            pickle.dump(start_cfd3, handle, protocol=pickle.HIGHEST_PROTOCOL)
        with open(f"pickles/start_cfd4_{num_to_process}_genes.pickle", 'wb') as handle:
            pickle.dump(start_cfd4, handle, protocol=pickle.HIGHEST_PROTOCOL)
        with open(f"pickles/start_cfdfinal_{num_to_process}_genes.pickle", 'wb') as handle:
            pickle.dump(start_cfdfinal, handle, protocol=pickle.HIGHEST_PROTOCOL)

        with open(f"pickles/stop_cfd1_{num_to_process}_genes.pickle", 'wb') as handle:
            pickle.dump(stop_cfd1, handle, protocol=pickle.HIGHEST_PROTOCOL)
        with open(f"pickles/stop_cfd2_{num_to_process}_genes.pickle", 'wb') as handle:
            pickle.dump(stop_cfd2, handle, protocol=pickle.HIGHEST_PROTOCOL)
        with open(f"pickles/stop_cfd3_{num_to_process}_genes.pickle", 'wb') as handle:
            pickle.dump(stop_cfd3, handle, protocol=pickle.HIGHEST_PROTOCOL)
        with open(f"pickles/stop_cfd4_{num_to_process}_genes.pickle", 'wb') as handle:
            pickle.dump(stop_cfd4, handle, protocol=pickle.HIGHEST_PROTOCOL)
        with open(f"pickles/stop_cfdfinal_{num_to_process}_genes.pickle", 'wb') as handle:
            pickle.dump(stop_cfdfinal, handle, protocol=pickle.HIGHEST_PROTOCOL)

        log.info(f"finished in {elapsed_min:.2f} min ({elapsed_sec} sec) , processed {protein_coding_transcripts_count}/{transcript_count} transcripts\nnonprotein-coding transcripts were skipped")

        recut_CFD_pass.close()
        recut_CFD_fail.close()

    except Exception as e:
        print("Unexpected error:", str(sys.exc_info()))
        traceback.print_exc()
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

if __name__ == "__main__": main()




