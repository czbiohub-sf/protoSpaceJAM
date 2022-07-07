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
from scripts.hdr import *
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
        gRNA_out_of_arms = dict()
        gRNA_out_of_arms["start"]=dict()
        gRNA_out_of_arms["stop"]=dict()

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
            if transcript_type == "protein_coding":
                log.info(f"processing {ENST_ID}\ttranscript type: {transcript_type}")
                ranked_df_gRNAs_ATG, ranked_df_gRNAs_stop = get_gRNAs(ENST_ID = ENST_ID, ENST_info= ENST_info, freq_dict = freq_dict, loc2file_index= loc2file_index, loc2posType = loc2posType, dist = max_cut2ins_dist)
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
                    HDR_template = get_HDR_template(df = best_start_gRNA, ENST_info = ENST_info, type = "start", ENST_PhaseInCodon = ENST_PhaseInCodon)

                    best_start_gRNAs = pd.concat([best_start_gRNAs, best_start_gRNA]) #append the best gRNA to the final df
                if best_stop_gRNA.empty == False:
                    if best_stop_gRNA.shape[0] > 1: # multiple best scoring gRNA
                        best_stop_gRNA = best_stop_gRNA[best_stop_gRNA["CSS"] == best_stop_gRNA["CSS"].max()] # break the tie by CSS score
                        best_stop_gRNA = best_stop_gRNA.head(1) #get the first row in case of ties

                    #get HDR template
                    HDR_template = get_HDR_template(df=best_stop_gRNA, ENST_info=ENST_info, type="stop", ENST_PhaseInCodon = ENST_PhaseInCodon)

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
                gnum = len(gRNA_out_of_arms["start"]) + len(gRNA_out_of_arms["stop"])
                log.info(f"number of gRNAs outside the HDR arm: {gnum}")

            ################
            #early stopping#
            ################
            #if ENST_ID == "ENST00000360426":
            #    sys.exit()
            # num_to_process = 18000
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
        with open(f"pickles/gRNA_out_of_arms_{num_to_process}_genes.pickle", 'wb') as handle:
            pickle.dump(gRNA_out_of_arms, handle, protocol=pickle.HIGHEST_PROTOCOL)

        log.info(f"finished in {elapsed_min:.2f} min ({elapsed_sec} sec) , processed {protein_coding_transcripts_count}/{transcript_count} transcripts\nnonprotein-coding transcripts were skipped")

    except Exception as e:
        print("Unexpected error:", str(sys.exc_info()))
        traceback.print_exc()
        print("additional information:", e)
        PrintException()

##########################
## function definitions ##
##########################
def get_phase_in_codon(Chr,Pos,ENST_ID, ENST_PhaseInCodon):
    """
    get the phase in codon for the current position and flanking 2 bp
    returns a dictionary:
    {0:[-1/-2/-3/1/2/3/0] 0 stands for not coding sequence (current position is not in a codon)
    -1:[-1/-2/-3/1/2/3/0]
    -2:[-1/-2/-3/1/2/3/0]
    +1:[-1/-2/-3/1/2/3/0]
    +2:[-1/-2/-3/1/2/3/0]}
    """
    mydict = {-2:0, -1:0, 0:0, 1:0, 2:0}
    if Chr in ENST_PhaseInCodon.keys():
        Chr_dict = ENST_PhaseInCodon[Chr]
        if Pos in Chr_dict.keys():
            if ENST_ID in Chr_dict[Pos].keys():
                mydict[0] = Chr_dict[Pos][ENST_ID]
        if Pos+1 in Chr_dict.keys():
            if ENST_ID in Chr_dict[Pos+1].keys():
                mydict[+1] = Chr_dict[Pos+1][ENST_ID]
        if Pos+2 in Chr_dict.keys():
            if ENST_ID in Chr_dict[Pos+2].keys():
                mydict[+2] = Chr_dict[Pos+2][ENST_ID]
        if Pos-1 in Chr_dict.keys():
            if ENST_ID in Chr_dict[Pos-1].keys():
                mydict[-1] = Chr_dict[Pos-1][ENST_ID]
        if Pos-2 in Chr_dict.keys():
            if ENST_ID in Chr_dict[Pos-2].keys():
                mydict[-2] = Chr_dict[Pos-2][ENST_ID]
    return mydict

def get_phase_in_codon0(Chr,Pos,ENST_ID, ENST_PhaseInCodon):
    """
    get the phase in codon for the current position and flanking 2 bp
    returns an int:
    -1/-2/-3/1/2/3/0 0 stands for not coding sequence (current position is not in a codon)
    """
    myInt = 0
    if Chr in ENST_PhaseInCodon.keys():
        Chr_dict = ENST_PhaseInCodon[Chr]
        if Pos in Chr_dict.keys():
            if ENST_ID in Chr_dict[Pos].keys():
                myInt = Chr_dict[Pos][ENST_ID]
    return myInt

def get_range(start,end): #TODO phase should be ENST specific
    """
    return a list of number, start to end, step size = 1, -1 (if start > end)
    """
    if start<=end:
        return list(range(start,end+1,1))
    else:
        return list(range(start,end-1,-1))



def get_HDR_template(df, ENST_info,type,ENST_PhaseInCodon):
    for index, row in df.iterrows():
        ENST_ID = row["ID"]
        ENST_strand = ENST_info[ENST_ID].features[0].strand
        Chr = row["chr"]
        InsPos = row["Insert_pos"] # InsPos is the first letter of stop codon "T"AA or the last letter of the start codon AT"G"
        gStart = row["start"]
        gStrand = convert_strand(row["strand"]) #gRNA strand
        CutPos = get_cut_pos(gStart,gStrand)

        ##########################
        # important debug info
        # print(f"ENST: {ENST_ID} ENST_strand: {ENST_strand} chr {Chr} type {type}-tagging InsPos {InsPos} gStrand {gStrand} CutPos {CutPos} gStart {gStart} ")
        ##########################
        #CutPos_phase = get_phase_in_codon(Chr=Chr, Pos=InsPos, ENST_PhaseInCodon=ENST_PhaseInCodon)
        #print(f"CutPos phase: {CutPos_phase}")

        #get target_seq
        #target_seq = get_target_seq(Chr= Chr, InsPos = InsPos, gRNAstrand = gStrand, CutPos = CutPos , type = type, ENST_ID = ENST_ID, ENST_info = ENST_info)

        leftArm, rightArm, left_start, left_end, right_start, right_end = get_HDR_arms(loc = [Chr,InsPos,ENST_strand], half_len = HDR_arm_len, type = type) # start>end is possible
        left_start_phase = get_phase_in_codon(Chr=Chr, Pos=left_start, ENST_ID=ENST_ID, ENST_PhaseInCodon=ENST_PhaseInCodon)
        left_end_phase = get_phase_in_codon(Chr=Chr, Pos=left_end, ENST_ID=ENST_ID, ENST_PhaseInCodon=ENST_PhaseInCodon)
        right_start_phase = get_phase_in_codon(Chr=Chr, Pos=right_start, ENST_ID=ENST_ID, ENST_PhaseInCodon=ENST_PhaseInCodon)
        right_end_phase = get_phase_in_codon(Chr=Chr, Pos=right_end, ENST_ID=ENST_ID, ENST_PhaseInCodon=ENST_PhaseInCodon)
        # print(f"{leftArm}\t"
        #       f"{left_start}({left_start_phase})\t"
        #       f"{left_end}({left_end_phase})\t"
        #       f"{rightArm}\t"
        #       f"{right_start}({right_start_phase})\t"
        #       f"{right_end}({right_end_phase})")

        left_Arm_Phases = [get_phase_in_codon0(Chr=Chr, Pos=i,ENST_ID=ENST_ID, ENST_PhaseInCodon=ENST_PhaseInCodon) for i in get_range(left_start,left_end)]
        right_Arm_Phases = [get_phase_in_codon0(Chr=Chr, Pos=i,ENST_ID=ENST_ID, ENST_PhaseInCodon=ENST_PhaseInCodon) for i in get_range(right_start,right_end)]

        myflank = HDR_flank(left_flk_seq = leftArm , right_flk_seq = rightArm,
                            left_flk_coord_lst = [left_start, left_end], right_flk_coord_lst = [right_start, right_end],
                            left_flk_phases = left_Arm_Phases, right_flk_phases = right_Arm_Phases,
                            type= type, ENST_ID= ENST_ID, ENST_strand=ENST_strand, gStart= gStart, gStrand= gStrand)

        #log IDs whose gRNA is not in the default-size HDR arms
        if myflank.entire_gRNA_in_HDR_arms == False:
            gRNA_out_of_arms[type][ENST_ID] = False


# def get_target_seq(Chr, InsPos, gRNAstrand, CutPos, type, ENST_ID, ENST_info):
#     """
#     The target sequence is used to instantiate the HDR class (gdingle)
#     The target sequence should be in the direction of the gene. Reading from
#     left to right, it should have either a ATG or one of TAG, TGA, or TAA.
#     """
#     ATG_loc, stop_loc = get_start_stop_loc(ENST_ID, ENST_info)
#     if type == "start":
#         target_codon_loc = ATG_loc
#         target_codon_strand = ATG_loc[3]
#         if target_codon_strand == 1:
#             ATG_loc[2]
#     elif type == "stop":
#         target_codon_loc = stop_loc
#         target_codon_strand = stop_loc[3]
#     else:
#         sys.exit(f"unknown type: {type}")

def convert_strand(strand):
    #convert strand from +/- to 1/-1
    if strand == "+":
        return(1)
    elif strand == "-":
        return(-1)
    else:
        return(f"input strand:{strand} needs to be +/-")

def get_HDR_arms(loc, half_len, type):
    """
    input:  loc         [chr,pos,strand]  #start < end , strand is the coding strand
            half_len      length of the HDR arm (one sided)
    return: HDR arms -> in coding strand <-
            [5'arm, 3'arm]
    """
    Chr,Pos,Strand = loc
    #get arms
    if type == "start":
        if Strand == 1:
            vanilla_left_arm = get_seq(chr = Chr, start = Pos-half_len+1, end = Pos+1, strand = 1)
            vanilla_right_arm = get_seq(chr = Chr, start = Pos+1, end = Pos+half_len+1, strand = 1)
            return([vanilla_left_arm,vanilla_right_arm, Pos-half_len+1,Pos,Pos+1,Pos+half_len])
        elif Strand == -1:
            vanilla_left_arm = get_seq(chr = Chr, start = Pos-half_len, end = Pos, strand = 1)
            vanilla_right_arm = get_seq(chr = Chr, start = Pos, end = Pos+half_len, strand = 1)
            return([reverse_complement(vanilla_right_arm),reverse_complement(vanilla_left_arm),Pos+half_len-1,Pos,Pos-1,Pos-half_len])
        else:
            sys.exit(f"unknown strand: {Strand}, acceptable values are -1 and 1")
    elif type == "stop":
        if Strand == 1:
            vanilla_left_arm = get_seq(chr = Chr, start = Pos-half_len, end = Pos, strand = 1)
            vanilla_right_arm = get_seq(chr = Chr, start = Pos, end = Pos+half_len, strand = 1)
            return([vanilla_left_arm,vanilla_right_arm,Pos-half_len,Pos-1,Pos,Pos+half_len-1])
        elif Strand == -1:
            vanilla_left_arm = get_seq(chr = Chr, start = Pos-half_len+1, end = Pos+1, strand = 1)
            vanilla_right_arm = get_seq(chr = Chr, start = Pos+1, end = Pos+half_len+1, strand = 1)
            return([reverse_complement(vanilla_right_arm),reverse_complement(vanilla_left_arm),Pos+half_len,Pos+1,Pos,Pos-half_len+1])
    else:
        sys.exit("unknown type {type}, acceptable values: start, stop")

    #rev_com ajustment
    if Strand == 1:
        return([vanilla_left_arm,vanilla_right_arm])
    elif Strand == -1:
        return([reverse_complement(vanilla_right_arm),reverse_complement(vanilla_left_arm)])
    else:
        sys.exit(f"unknown strand: {Strand}, acceptable values are -1 and 1")

def get_gRNAs(ENST_ID, ENST_info, freq_dict, loc2file_index, loc2posType, dist=50):
    """
    input
        ENST_ID: ENST ID
        ENST_info: gene model info loaded from pickle file
        loc2file_index: mapping of location to the file part that stores gRNAs in that region
        loc2posType: mapping of location to location type (e.g. 5UTR etc)
        dist: max cut to insert distance (default 50)
    return:
        a dictionary of guide RNAs which cuts <[dist] to end of start codon
        a dictionary of guide RNAs which cuts <[dist] to start of stop codon
    """
    # get location of start and stop location
    ATG_loc, stop_loc = get_start_stop_loc(ENST_ID, ENST_info)



    ##################################
    # get gRNAs around the start codon#
    ##################################
    # get the start codon chromosomal location
    log.debug(f"ATG_loc: {ATG_loc}")
    end_of_ATG_loc = get_end_pos_of_ATG(ATG_loc)  # [chr, pos,strand]
    log.debug(f"end of the ATG: {end_of_ATG_loc}")
    # get gRNA around the chromosomeal location (near ATG)
    df_gRNAs_ATG = get_gRNAs_near_loc(loc=end_of_ATG_loc, dist=dist, loc2file_index=loc2file_index)
    # rank gRNAs
    ranked_df_gRNAs_ATG = rank_gRNAs_for_tagging(loc=end_of_ATG_loc, gRNA_df=df_gRNAs_ATG, loc2posType=loc2posType, ENST_ID = ENST_ID)

    ##################################
    # get gRNAs around the stop  codon#
    ##################################
    # get gRNAs around the stop codon
    log.debug(f"stop_loc: {stop_loc}")
    start_of_stop_loc = get_start_pos_of_stop(stop_loc)
    log.debug(f"start of stop: {start_of_stop_loc}")  # [chr, pos,strand]
    # get gRNA around the chromosomeal location (near stop location)
    df_gRNAs_stop = get_gRNAs_near_loc(loc=start_of_stop_loc, dist=dist, loc2file_index=loc2file_index)
    # rank gRNAs
    ranked_df_gRNAs_stop = rank_gRNAs_for_tagging(loc=start_of_stop_loc, gRNA_df=df_gRNAs_stop, loc2posType=loc2posType, ENST_ID = ENST_ID)

    return ([ranked_df_gRNAs_ATG, ranked_df_gRNAs_stop])

def rank_gRNAs_for_tagging(loc,gRNA_df, loc2posType, ENST_ID, alpha = 1):
    """
    input:  loc         [chr,pos,strand]  #start < end
            gRNA_df     pandas dataframe, *unranked*   columns: "seq","pam","start","end", "strand", "CSS", "ES"  !! neg strand: start > end
            alpha       specificity weight is raised to the power of alpha
    output: gRNA_df     pandas dataframe *ranked*      columns: "seq","pam","start","end", "strand", "CSS", "ES"  !! neg strand: start > end
    """
    insPos = loc[1]
    Chr = loc[0]

    col_spec_weight = []
    col_dist_weight = []
    col_pos_weight = []
    col_final_weight = []
    Chrs = []
    ENSTs = []
    InsertPos = []
    #assign a score to each gRNA
    for index, row in gRNA_df.iterrows():
        start = row[2]
        end = row[3]
        strand = row[4]
        #Get cut to insert distance
        cutPos = get_cut_pos(start,strand)
        cut2insDist = cutPos - insPos

        #calc. specificity_weight
        CSS = row[5]
        specificity_weight = _specificity_weight(CSS)
        col_spec_weight.append(specificity_weight)

        #calc. distance_weight
        distance_weight = _dist_weight(hdr_dist = cut2insDist)
        col_dist_weight.append(distance_weight)

        #get position_weight
        position_type = _get_position_type(chr = Chr, ID = ENST_ID, pos = cutPos, loc2posType = loc2posType)
        position_weight = _position_weight(position_type)
        col_pos_weight.append(position_weight)

        #add info to the df
        Chrs.append(Chr)
        ENSTs.append(ENST_ID)
        InsertPos.append(insPos)

        log.debug(f"strand {strand} {start}-{end} cutPos {cutPos} insert_loc {loc} cut2insDist {cut2insDist} distance_weight {distance_weight:.2f} CFD_score {CSS} specificity_weight {specificity_weight} pos_type {position_type} position_weight {position_weight}")

        final_score = float(pow(specificity_weight,alpha)) * float(distance_weight) * float(position_weight)
        col_final_weight.append(final_score)

    #add info and weight columns to the df
    gRNA_df["chr"] = Chrs
    gRNA_df["ID"] = ENSTs
    gRNA_df["Insert_pos"] = InsertPos
    gRNA_df["spec_weight"] = col_spec_weight
    gRNA_df["dist_weight"] = col_dist_weight
    gRNA_df["pos_weight"] = col_pos_weight
    gRNA_df["final_weight"] = col_final_weight


    #rank gRNAs based on the score
    gRNA_df['final_pct_rank'] = gRNA_df['final_weight'].rank(pct=True)

    return gRNA_df

def get_cut_pos(start,strand):
    """
    start:gRNA start
    strand:gRNA strand
    """
    if strand == "+" or strand == "1" or strand == 1:
        cutPos = start + 16
    else:
        cutPos = start - 16
    return(cutPos)

def _get_position_type(chr, ID, pos, loc2posType):
    """
    input: mostly self-explanatory, loc2posType is a dictionary that translates location into types (e.g. exon/intron junctions etc)
    return a list of types for the input position/ID combination
    """
    chr_dict = loc2posType[chr]
    if not ID in chr_dict.keys():
        return []
    else:
        types = []
        mapping_dict = chr_dict[ID]
        for key in mapping_dict.keys():
            if in_interval_leftrightInclusive(pos,key):
                types.append(mapping_dict[key])
        return types

def _position_weight(types):
    """
    input: a list of types
    output: the lowest weight among all the types
    """
    mapping_dict = {"5UTR":0.4,
                    "3UTR":1,
                    "cds":1,
                    "within_2bp_of_exon_intron_junction":0.01,
                    "within_2bp_of_intron_exon_junction":0.01,
                    "3N4bp_up_of_exon_intron_junction":0.1,
                    "3_to_6bp_down_of_exon_intron_junction":0.1,
                    "3N4bp_up_of_intron_exon_junction":0.1,
                    "3N4bp_down_of_intron_exon_junction":0.5}
    lowest_weight = 1
    for type in types:
        if type in mapping_dict.keys():
            weight = mapping_dict[type]
            if weight < lowest_weight:
                lowest_weight = weight
        else:
            sys.exit(f"unexpected position type: {type}")
    return(lowest_weight)

def _dist_weight(hdr_dist: int, _dist_weight_variance = 55) -> float:
    """
    taken from https://github.com/czbiohub/crispycrunch
    >>> _dist_weight(0)
    1.0
    >>> _dist_weight(5)
    0.7967034698934616
    >>> _dist_weight(10)
    0.402890321529133
    >>> _dist_weight(-20)
    0.026347980814448734
    """
    variance = _dist_weight_variance

    hdr_dist = abs(hdr_dist)  # make symmetric
    assert hdr_dist >= 0 and hdr_dist <= 100  # 100 is resonable upper bound

    # Returns a gaussian
    weight = math.exp((-1 * hdr_dist**2) / (2 * variance))
    assert weight >= 0 and weight <= 1
    return weight

def _specificity_weight(specificity_score: float, _specificity_weight_low = 45, _specificity_weight_high = 65 ):
    """
    taken from https://github.com/czbiohub/crispycrunch
    >>> _specificity_weight(20)
    0
    >>> _specificity_weight(60)
    0.75
    >>> _specificity_weight(80)
    1
    """
    low = _specificity_weight_low
    high = _specificity_weight_high

    assert specificity_score >= 0 and specificity_score <= 100
    if specificity_score <= low:
        return 0
    elif specificity_score >= high:
        return 1
    else:
        return 1 / (high - low) * (specificity_score - low)

def get_end_pos_of_ATG(ATG_loc):
    """
    input: ATG_loc              [chr,start,end,strand] #start < end
    output: the pos of G in ATG [chr,pos,strand]       #start < end
    """
    strand = ATG_loc[3]
    if str(strand) == "+" or str(strand) == "1":
        return [ATG_loc[0],ATG_loc[2],ATG_loc[3]]
    else:
        return [ATG_loc[0],ATG_loc[1],ATG_loc[3]]

def get_start_pos_of_stop(stop_loc):
    """
    input: stop_loc                                 [chr,start,end,strand]  #start < end
    output: the pos of first base in the stop codon [chr,pos,strand]        #start < end
    """
    strand = stop_loc[3]
    if str(strand) == "+" or str(strand) == "1":
        return [stop_loc[0],stop_loc[1],stop_loc[3]]
    else:
        return [stop_loc[0],stop_loc[2],stop_loc[3]]

def get_start_stop_loc(ENST_ID,ENST_info):
    """
    Get the chromosomal location of start and stop codons
    input: ENST_ID, ENST_info
    output: a list of two items
            ATG_loc: [chr,start,end,strand]  #start < end
            stop_loc: [chr,start,end,strand] #start < end
    """
    my_transcript = ENST_info[ENST_ID]  # get the seq record
    # constructing the list of cds
    cdsList = [feat for feat in my_transcript.features if feat.type == 'CDS']
    CDS_first = cdsList[0]
    CDS_last = cdsList[len(cdsList) - 1]
    #get start codon location
    if CDS_first.strand==1:
        ATG_loc = [CDS_first.location.ref, CDS_first.location.start+0,CDS_first.location.start+2, 1] # format [start, end, strand]
    else:
        stop_loc = [CDS_first.location.ref, CDS_first.location.start+0,CDS_first.location.start+2, -1]
    #get stop codon location
    if CDS_last.strand==1:
        stop_loc = [CDS_last.location.ref, CDS_last.location.end-2,CDS_last.location.end+0, 1]
    else:
        ATG_loc = [CDS_last.location.ref, CDS_last.location.end-2,CDS_last.location.end+0, -1]

    return([ATG_loc, stop_loc])

def get_gRNAs_near_loc(loc,dist, loc2file_index):
    """
    input
        loc: [chr,pos,strand]
        dist: max cut to loc distance
    return:
        a dataframe of guide RNAs which cuts <[dist] to the loc, the columns are "seq","pam","start","end", "strand", "CSS", "ES"
    """
    chr = loc[0]
    pos = loc[1]
    chr_dict = loc2file_index[chr]
    target_files =[] # a list of file names containing gRNAs near loc
    #lookup the file
    for key in chr_dict.keys():
        file_start = key.split("-")[0]
        file_end =  key.split("-")[1]
        interval = [file_start, file_end]
        if in_interval(pos, interval) or in_interval(pos-1000, interval) or in_interval(pos+1000, interval):
            target_files.append(chr_dict[key])
    #print(target_files)
    #load the gRNAs in the file
    dfs =[]
    for file in target_files:
        file_path = os.path.join(f"gRNA_{config['genome_ver']}","gRNA.tab.gz.split.BwaMapped.scored",file)
        df_tmp = pd.read_csv(file_path, sep="\t", compression='infer', header=None, names = ["seq","pam","start","end", "strand", "CSS", "ES"])
        dfs.append(df_tmp)
    df_gRNA = pd.concat(dfs)

    #subset gRNA based on strand  !ATTN: start > end when strand is '-'
    df_gRNA_on_sense = df_gRNA[(df_gRNA['strand'] == '+')]
    df_gRNA_on_antisense = df_gRNA[(df_gRNA['strand'] == '-')]

    #subset gRNAs and retain those cuts <[dist] to the loc
    df_gRNA_on_sense = df_gRNA_on_sense[(df_gRNA_on_sense['start'] > (pos-17-dist)) & (df_gRNA_on_sense['start'] < (pos-17+dist))]
    df_gRNA_on_antisense = df_gRNA_on_antisense[(df_gRNA_on_antisense['start'] > (pos+17-dist)) & (df_gRNA_on_antisense['start'] < (pos+17+dist))]

    return(pd.concat([df_gRNA_on_sense,df_gRNA_on_antisense]))

def in_interval(pos,interval):
    """
    check if pos in is interval
    input
        pos
        interval: [start,end]
    return: boolean
    """
    pos = int(pos)
    interval = list(interval)
    interval[0] = int(interval[0])
    interval[1] = int(interval[1])
    if pos >= interval[0] and pos <= interval[1]:
        #log.debug(f"{pos} is in {interval}")
        return True
    else:
        return False

def in_interval_leftrightInclusive(pos,interval):
    """
    check if pos in is interval
    input
        pos
        interval: [start,end]
    return: boolean
    """
    pos = int(pos)
    interval = list(interval)
    interval[0] = int(interval[0])
    interval[1] = int(interval[1])
    if pos >= interval[0] and pos <= interval[1]:
        #log.debug(f"{pos} is in {interval}")
        return True
    else:
        return False

def update_dict_count(key,dict): # update the dictionary that keeps the count of each string (key)
    if key in dict.keys():
        dict[key]+=1
    else:
        dict[key]=1
    return dict

def get_seq(chr,start,end,strand):
    '''
    chr
    start (1-indexed)
    end (the end position is not included
    strand 1 or -1 (str)
    '''
    chr_file_path = os.path.join("genome_files",f"{config['genome_ver']}_byChr",f"{chr}.pk")
    log.debug(f"opening file {chr_file_path}")
    if os.path.isfile(chr_file_path):
        #read file
        chr_seqrecord = read_pickle_files(chr_file_path)
        subseq = chr_seqrecord.seq._data.decode()[(start-1):(end-1)] # use -1 to convert 1-index to 0-index
        if strand == "-1" or strand == -1:
            return(reverse_complement(subseq))
        else:
            return(subseq)
    else:
        sys.exit(f"ERROR: file not found: {chr_file_path}")

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




