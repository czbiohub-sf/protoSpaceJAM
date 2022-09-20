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
    parser= MyParser(description='ProtospaceX')
    parser.add_argument('--genome_ver', default="GRCh38", type=str, help='pickle file containing the ENST_info dict', metavar='')
    parser.add_argument('--path2csv', default="input/mart_export_canonical_proteincoding.csv", type=str,help='path to a csv file containing ENST information\n *required columns*: Ensemble_ID',metavar='')

    #gRNA
    parser.add_argument('--num_gRNA_per_term',  default=1, type=int, help='payload at the N terminus', metavar='')

    #donor
    parser.add_argument('--HA_len',  default=500, help='length of the homology arm on one side', type=int, metavar='')
    parser.add_argument('--ssODN_max_size', type=int, help='length restraint of the ssODN (both arms + payload), setting this option will center the ssODN with respect to the payload and the recoded region', metavar='')

    #payload
    parser.add_argument('--payload', default="", type=str, help='payload, overrides --Npayloadf and --Cpayload', metavar='')
    parser.add_argument('--Npayload', default="ACCGAGCTCAACTTCAAGGAGTGGCAAAAGGCCTTTACCGATATGATGGGTGGCGGATTGGAAGTTTTGTTTCAAGGTCCAGGAAGTGGT", type=str, help='payload at the N terminus', metavar='')
    parser.add_argument('--Cpayload',  default="GGTGGCGGATTGGAAGTTTTGTTTCAAGGTCCAGGAAGTGGTACCGAGCTCAACTTCAAGGAGTGGCAAAAGGCCTTTACCGATATGATG", type=str, help='payload at the N terminus', metavar='')

    config = parser.parse_args()
    return config

logging.setLoggerClass(ColoredLogger)
#logging.basicConfig()
log = logging.getLogger("ProtospaceX")
log.propagate = False
log.setLevel(logging.INFO) #set the level of warning displayed
#log.setLevel(logging.DEBUG) #set the level of warning displayed

config = vars(parse_args())
gRNA_num_out = config['num_gRNA_per_term']
max_cut2ins_dist = 50 #deprecated?
HDR_arm_len = config['HA_len']
ssODN_max_size = config["ssODN_max_size"]

#check if HA_len is too short to satisfy ssODN_max_size
if ssODN_max_size is not None:
    max_payload_size = max([len(config["Npayload"]),len(config["Cpayload"])])
    derived_HDR_arm_len = ssODN_max_size- max_payload_size / 2
    if derived_HDR_arm_len >= HDR_arm_len:
        print(f"HA_len={HDR_arm_len} is to short to meet the requirement of ssODN_max_size={ssODN_max_size}, payload size={max_payload_size}\n ssODN_max_size={ssODN_max_size} requires HA_len = ssODN_max_size- max_payload_size / 2 = {derived_HDR_arm_len}")
        HDR_arm_len = derived_HDR_arm_len + 100
        print(f"HA_len is adjusted to {HDR_arm_len}")

#####################
##      main       ##
#####################
def main():
    try:
        starttime = datetime.datetime.now()
        freq_dict = dict()

        #load gRNA info index (mapping of chromosomal location to file parts)
        log.info("loading the mapping of chromosomal location to (gRNA) file parts")
        loc2file_index = read_pickle_files(os.path.join("genome_files", "fa_pickle", config['genome_ver'],"loc2file_index.pickle"))

        #load chr location to type (e.g. UTR, cds, exon/intron junction) mappings
        log.info("loading the mapping of chromosomal location to type (e.g. UTR, cds, exon/intron junction)")
        loc2posType = read_pickle_files(os.path.join("genome_files","parsed_gff3", config['genome_ver'],"loc2posType.pickle"))

        #load gene model info
        log.info("loading gene model info")
        ENST_info = read_pickle_files(os.path.join("genome_files","parsed_gff3", config['genome_ver'],"ENST_info.pickle"))

        #load codon phase info
        log.info("loading codon phase info")
        ENST_PhaseInCodon = read_pickle_files(os.path.join("genome_files","parsed_gff3", config['genome_ver'],"ENST_codonPhase.pickle"))

        #report time used
        elapsed = cal_elapsed_time(starttime,datetime.datetime.now())
        log.info(f"finished loading in {elapsed[0]:.2f} min ({elapsed[1]} sec)")

        #report the number of ENSTs which has ATG at the end of the exon
        ExonEnd_ATG_count,ExonEnd_ATG_list = count_ATG_at_exonEnd(ENST_info)

        #open log files
        recut_CFD_all = open("logs/recut_CFD_all.txt", "w")
        recut_CFD_fail = open("logs/recut_CFD_fail.txt", "w")
        csvout_N = open("logs/out_Nterm_recut_cfd.csv", "w")
        csvout_C = open("logs/out_Cterm_recut_cfd.csv", "w")
        csvout_header = "ID,cfd1,cfd2,cfd3,cfd4,max_of_cfd4_cfdScan,cfd_max\n"
        csvout_N.write(csvout_header)
        csvout_C.write(csvout_header)
        fiveUTR_log = open("logs/fiveUTR.txt", "w")

        #open result file and write header
        csvout_res = open("logs/result.csv", "w")
        csvout_res.write(f"ID,chr,transcript_type,name,terminus,gRNA_seq,PAM,gRNA_start,gRNA_end,distance_between_cut_and_edit,CFD_specificity_score,specificity_weight,distance_weight,position_weight,final_weight,cfd_after_recoding,cfd_after_windowScan_and_recoding,max_recut_cfd,ssODN,effective_HA_len\n")

        #dataframes to store best gRNAs
        best_start_gRNAs = pd.DataFrame()
        best_stop_gRNAs = pd.DataFrame()

        #logging cfd score, failed gRNAs etc
        start_info = info()
        stop_info = info()

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

        #loop through each ENST
        transcript_count = 0
        protein_coding_transcripts_count = 0
        target_terminus = "all"
        for index, row in df.iterrows():
            ENST_ID = row["Ensemble_ID"]
            if "Target Terminus" in df.columns:
                target_terminus = row["Target Terminus"]
                if target_terminus!="N" and target_terminus!="C" and target_terminus!="all":
                    sys.exit(f"invalid target terminus: {target_terminus}")
            if not ENST_ID in ENST_info.keys():
                log.warning(f"skipping {ENST_ID} b/c transcript is not in the annotated ENST collection (excluding those on chr_patch_hapl_scaff)")
                continue
            transcript_type = ENST_info[ENST_ID].description.split("|")[1]
            if transcript_type == "protein_coding": # and ENST_ID == "ENST00000329276":
                # if not ENST_ID in ExonEnd_ATG_list: # only process edge cases in which genes with ATG are at the end of exons
                #     continue
                log.info(f"processing {ENST_ID}\ttranscript type: {transcript_type}")
                csvout_N.write(ENST_ID)
                csvout_C.write(ENST_ID)
                if hasattr(ENST_info[ENST_ID],"name"):
                    name = ENST_info[ENST_ID].name
                else:
                    name = ""
                row_prefix = f"{ENST_ID},{ENST_info[ENST_ID].chr},{transcript_type},{name}"

                #get gRNAs
                ranked_df_gRNAs_ATG, ranked_df_gRNAs_stop = get_gRNAs(ENST_ID = ENST_ID, ENST_info= ENST_info, freq_dict = freq_dict, loc2file_index= loc2file_index, loc2posType = loc2posType, dist = max_cut2ins_dist, genome_ver=config["genome_ver"])
                if ranked_df_gRNAs_ATG.empty == True:
                    start_info.failed.append(ENST_ID)
                    csvout_N.write(",,,,,\n")
                if ranked_df_gRNAs_stop.empty == True:
                    stop_info.failed.append(ENST_ID)
                    csvout_C.write(",,,,,\n")

                ##################################
                #best start gRNA and HDR template#
                ##################################
                if target_terminus=="all" or target_terminus=="N":
                    for i in range(0,min([gRNA_num_out, ranked_df_gRNAs_ATG.shape[0]])):
                        # if best_start_gRNA.shape[0] > 1: # multiple best scoring gRNA
                        #     best_start_gRNA = best_start_gRNA[best_start_gRNA["CSS"] == best_start_gRNA["CSS"].max()] # break the tie by CSS score
                        #     best_start_gRNA = best_start_gRNA.head(1) #get the first row in case of ties
                        current_gRNA = ranked_df_gRNAs_ATG.iloc[[i]]

                        #get HDR template
                        HDR_template = get_HDR_template(df = current_gRNA, ENST_info = ENST_info, type = "start", ENST_PhaseInCodon = ENST_PhaseInCodon, HDR_arm_len=HDR_arm_len, genome_ver=config["genome_ver"], tag = config["Npayload"], loc2posType = loc2posType, ssODN_max_size = ssODN_max_size)

                        # append the best gRNA to the final df
                        if i==0:
                            best_start_gRNAs = pd.concat([best_start_gRNAs, current_gRNA])

                        #append cfd score to list for plotting
                        cfd1 = HDR_template.cdf_score_post_mut_ins
                        if not hasattr(HDR_template,"cdf_score_post_mut2"):
                            cfd2 = cfd1
                        else:
                            cfd2 = HDR_template.cdf_score_post_mut2
                        if not hasattr(HDR_template,"cdf_score_post_mut3"):
                            cfd3 = cfd2
                        else:
                            cfd3 = HDR_template.cdf_score_post_mut3
                        if not hasattr(HDR_template,"cdf_score_post_mut4"):
                            cfd4 = cfd3
                        else:
                            cfd4 = HDR_template.cdf_score_post_mut4
                        start_info.cfd4.append(cfd4)
                        cfd_scan = 0
                        if hasattr(HDR_template,"cdf_score_highest_in_win_scan"):
                            cfd_scan = HDR_template.cdf_score_highest_in_win_scan

                        cfdfinal = HDR_template.final_cfd

                        #write csv
                        csvout_N.write(f",{cfd1:.6f},{cfd2:.6f},{cfd3:.6f},{cfd4:.6f},{cfd_scan:.6f},{cfdfinal:.6f}\n")
                        CSS, seq, pam, s, e, cut2ins_dist, spec_weight, dist_weight, pos_weight, final_weight = get_res(current_gRNA)
                        ssODN = HDR_template.ODN_final_ss
                        csvout_res.write(f"{row_prefix},N,{seq},{pam},{s},{e},{cut2ins_dist},{CSS},{spec_weight:.6f},{dist_weight:.6f},{pos_weight:.6f},{final_weight:.6f},{cfd4:.6f},{cfd_scan:.6f},{cfdfinal:.6f},{ssODN},{HDR_template.effective_HA_len}\n")


                        #write log
                        this_log = f"{HDR_template.info}{HDR_template.info_arm}{HDR_template.info_p1}{HDR_template.info_p2}{HDR_template.info_p3}{HDR_template.info_p4}{HDR_template.info_p5}--------------------final CFD:{HDR_template.final_cfd:.6f}\n    ssODN before any recoding:{HDR_template.ODN_vanillia}\n     ssODN after all recoding:{HDR_template.ODN_postMut}\nssODN centered(if applicable):{HDR_template.ODN_postMut_centered}\n          ssODN (best strand):{HDR_template.ODN_final_ss}\n\n"
                        recut_CFD_all.write(this_log)
                        if HDR_template.final_cfd > 0.03:
                            recut_CFD_fail.write(this_log)

                        if hasattr(HDR_template,"info_phase4_5UTR"):
                            fiveUTR_log.write(f"phase4_UTR\t{HDR_template.info_phase4_5UTR[0]}\t{HDR_template.info_phase4_5UTR[1]}\n")
                        if hasattr(HDR_template,"info_phase5_5UTR"):
                            fiveUTR_log.write(f"phase5_UTR\t{HDR_template.info_phase5_5UTR[0]}\t{HDR_template.info_phase5_5UTR[1]}\n")

                #################################
                #best stop gRNA and HDR template#
                #################################
                if target_terminus=="all" or target_terminus=="C":
                    for i in range(0,min([gRNA_num_out, ranked_df_gRNAs_stop.shape[0]])):
                        # if best_stop_gRNA.shape[0] > 1: # multiple best scoring gRNA
                        #     best_stop_gRNA = best_stop_gRNA[best_stop_gRNA["CSS"] == best_stop_gRNA["CSS"].max()] # break the tie by CSS score
                        #     best_stop_gRNA = best_stop_gRNA.head(1) #get the first row in case of ties
                        current_gRNA = ranked_df_gRNAs_stop.iloc[[i]]

                        #get HDR template
                        HDR_template = get_HDR_template(df=current_gRNA, ENST_info=ENST_info, type="stop", ENST_PhaseInCodon = ENST_PhaseInCodon, HDR_arm_len = HDR_arm_len, genome_ver=config["genome_ver"], tag = config["Cpayload"], loc2posType = loc2posType, ssODN_max_size = ssODN_max_size)

                        # append the best gRNA to the final df
                        best_stop_gRNAs = pd.concat([best_stop_gRNAs, current_gRNA])

                        #append cfd score to list for plotting
                        cfd1 = HDR_template.cdf_score_post_mut_ins
                        if not hasattr(HDR_template,"cdf_score_post_mut2"):
                            cfd2 = cfd1
                        else:
                            cfd2 = HDR_template.cdf_score_post_mut2
                        if not hasattr(HDR_template,"cdf_score_post_mut3"):
                            cfd3 = cfd2
                        else:
                            cfd3 = HDR_template.cdf_score_post_mut3
                        if not hasattr(HDR_template,"cdf_score_post_mut4"):
                            cfd4 = cfd3
                        else:
                            cfd4 = HDR_template.cdf_score_post_mut4

                        cfd_scan = 0
                        if hasattr(HDR_template,"cdf_score_highest_in_win_scan"):
                            cfd_scan = HDR_template.cdf_score_highest_in_win_scan

                        cfdfinal = HDR_template.final_cfd

                        #write csv
                        csvout_C.write(f",{cfd1:.6f},{cfd2:.6f},{cfd3:.6f},{cfd4:.6f},{cfd_scan:.6f},{cfdfinal:.6f}\n")
                        CSS, seq, pam, s, e, cut2ins_dist, spec_weight, dist_weight, pos_weight, final_weight = get_res(current_gRNA)
                        ssODN = HDR_template.ODN_final_ss
                        #csvout_res.write(f"{row_prefix},C,{seq},{pam},{s},{e},{cut2ins_dist},{CSS},{spec_weight:.6f},{dist_weight:.6f},{pos_weight:.6f},{final_weight:.6f},{cfd1:.6f},{cfd2:.6f},{cfd3:.6f},{cfd4:.6f},{cfd_scan:.6f},{cfdfinal:.6f},{ssODN}\n")
                        csvout_res.write(f"{row_prefix},N,{seq},{pam},{s},{e},{cut2ins_dist},{CSS},{spec_weight:.6f},{dist_weight:.6f},{pos_weight:.6f},{final_weight:.6f},{cfd4:.6f},{cfd_scan:.6f},{cfdfinal:.6f},{ssODN},{HDR_template.effective_HA_len}\n")

                        #write log
                        this_log = f"{HDR_template.info}{HDR_template.info_arm}{HDR_template.info_p1}{HDR_template.info_p2}{HDR_template.info_p3}{HDR_template.info_p4}{HDR_template.info_p5}--------------------final CFD:{HDR_template.final_cfd:.6f}\n   ssODN before any recoding:{HDR_template.ODN_vanillia}\n    ssODN after all recoding:{HDR_template.ODN_postMut}\n             ssODN centered:{HDR_template.ODN_postMut_centered}\nssODN centered (best strand):{HDR_template.ODN_final_ss}\n\n"
                        recut_CFD_all.write(this_log)
                        if HDR_template.final_cfd > 0.03:
                            recut_CFD_fail.write(this_log)

                        if hasattr(HDR_template,"info_phase4_5UTR"):
                            fiveUTR_log.write(f"phase4_UTR\t{HDR_template.info_phase4_5UTR[0]}\t{HDR_template.info_phase4_5UTR[1]}\n")
                        if hasattr(HDR_template,"info_phase5_5UTR"):
                            fiveUTR_log.write(f"phase5_UTR\t{HDR_template.info_phase5_5UTR[0]}\t{HDR_template.info_phase5_5UTR[1]}\n")


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
            # if ENST_ID == "ENST00000562221":
            #     sys.exit()
            # num_to_process = 100
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

        # write best gRNA dfs to file
        with open(f"pickles/best_start_gRNAs_of_{num_to_process}_genes.pickle", 'wb') as handle:
            pickle.dump(best_start_gRNAs, handle, protocol=pickle.HIGHEST_PROTOCOL)
        with open(f"pickles/best_stop_gRNAs_of_{num_to_process}_genes.pickle", 'wb') as handle:
            pickle.dump(best_stop_gRNAs, handle, protocol=pickle.HIGHEST_PROTOCOL)

        # write failed ENSTs to file
        with open(f"pickles/start_failed_IDs_of_{num_to_process}_genes.pickle", 'wb') as handle:
            pickle.dump(start_info.failed, handle, protocol=pickle.HIGHEST_PROTOCOL)
        with open(f"pickles/stop_failed_IDs_of_{num_to_process}_genes.pickle", 'wb') as handle:
            pickle.dump(stop_info.failed, handle, protocol=pickle.HIGHEST_PROTOCOL)

        # write ENSTs (whose gRNA is outside of the default HDR arm) to file
        #with open(f"pickles/gRNA_out_of_arms_{num_to_process}_genes.pickle", 'wb') as handle:
        #    pickle.dump(gRNA_out_of_arms, handle, protocol=pickle.HIGHEST_PROTOCOL)

        log.info(f"finished in {elapsed_min:.2f} min ({elapsed_sec} sec) , processed {protein_coding_transcripts_count}/{transcript_count} transcripts\nnonprotein-coding transcripts were skipped")

        recut_CFD_all.close()
        recut_CFD_fail.close()
        csvout_N.close()
        csvout_C.close()

    except Exception as e:
        print("Unexpected error:", str(sys.exc_info()))
        traceback.print_exc()
        print("additional information:", e)
        PrintException()

##########################
## function definitions ##
##########################
class info:
    '''
    info log class
    '''
    def __init__(self)->None:
        self.cfd1 = []
        self.cfd2 = []
        self.cfd3 = []
        self.cfd4 = []
        self.cfdfinal = []
        self.failed = []

def get_res(best_start_gRNA):
    CSS = best_start_gRNA["CSS"].values[0]
    seq = best_start_gRNA["seq"].values[0]
    pam = best_start_gRNA["pam"].values[0]
    s = best_start_gRNA["start"].values[0]
    e = best_start_gRNA["end"].values[0]
    cut2ins_dist = best_start_gRNA["Cut2Ins_dist"].values[0]
    spec_weight = best_start_gRNA["spec_weight"].values[0]
    dist_weight = best_start_gRNA["dist_weight"].values[0]
    pos_weight = best_start_gRNA["pos_weight"].values[0]
    final_weight = best_start_gRNA["final_weight"].values[0]
    return([CSS,seq,pam,s,e,cut2ins_dist,spec_weight,dist_weight,pos_weight,final_weight])

def PrintException():
    exc_type, exc_obj, tb = sys.exc_info()
    f = tb.tb_frame
    lineno = tb.tb_lineno
    filename = f.f_code.co_filename
    linecache.checkcache(filename)
    line = linecache.getline(filename, lineno, f.f_globals)
    print('EXCEPTION IN ({}, LINE {} "{}"): {}'.format(filename, lineno, line.strip(), exc_obj))

if __name__ == "__main__": main()




