# -*- coding: future_fstrings -*-
import os.path
import pickle
import pandas as pd
from Bio import SeqIO
import csv
import argparse
import sys
import linecache
import datetime
import gzip
import shutil
import re
import math
import ast
from pathlib import Path
from utils import *


scoring_script_dir = os.path.join(sys.path[0], "..", "utils", "crisporScores")
sys.path.insert(1, scoring_script_dir)

from crisporEffScores import *

def parse_args():
    parser = MyParser(description='')
    parser.add_argument('--gzfile', default="", type=str, help='path to the gzfile', metavar='')
    parser.add_argument('--gzdir', default="", type=str, help='path to the dir containing gzfile', metavar='')
    parser.add_argument('--skip_eff_score', default = False, action='store_true', help='skip calculation of the efficieny score')
    parser.add_argument("--pam", default="NGG", type=str,  help="PAM sequence, default: NGG", metavar="",
    )

    config = parser.parse_args()
    if len(sys.argv)==1: # print help message if arguments are not valid
        parser.print_help()
        sys.exit(1)
    return config

logging.setLoggerClass(ColoredLogger)
#logging.basicConfig()
log = logging.getLogger("calcScores")
log.propagate = False
log.setLevel(logging.INFO) #set the level of warning displayed

config = vars(parse_args())

PAM = config["pam"].split("|")[0].upper()
assert PAM in ["NGG", "NGA", "TTTV"], "pam {} is not supported, please use NGG, NGA or TTTV".format(PAM)

#load scoring matrices
def get_mm_pam_scores():
    """
    """
    dataDir = Path(__file__).parent.parent / "utils" / "crisporScores" / "bin" / "CFD_Scoring"
    mm_scores = pickle.load(open(str(dataDir / "mismatch_score.pkl"), "rb"))
    pam_scores = pickle.load(open(str(dataDir / "pam_scores.pkl"), "rb"))
    return (mm_scores, pam_scores)

def get_mm_pam_scores_cas12a():
    """
    """
    dataDir = Path(__file__).parent.parent / "utils" / "crisporScores" / "bin" / "CFD_Scoring_Cas12a"
    df = pd.read_csv(dataDir / "off_targ_enCas12a.csv")
    df['mm_pos'] = df["MM"] + "," + df["Pos"].astype(str)
    mm_scores = dict(zip(df['mm_pos'], df['avg_percent_active']))
    return mm_scores

global mm_scores_cas12a
mm_scores_cas12a = get_mm_pam_scores_cas12a()
global mm_scores
global pam_scores
mm_scores,pam_scores = get_mm_pam_scores()
global hitScoreM
hitScoreM = [
    0,
    0,
    0.014,
    0,
    0,
    0.395,
    0.317,
    0,
    0.389,
    0.079,
    0.445,
    0.508,
    0.613,
    0.851,
    0.732,
    0.828,
    0.615,
    0.804,
    0.685,
    0.583,
]


#####################
##      main       ##
#####################
def main():
    try:

        #check input
        if (config["gzfile"] is None or config["gzfile"] == ""):
            log.error("please specify --gzfile")
            sys.exit("Please fix the error(s) above and rerun the script")

        starttime = datetime.datetime.now()
        gRNA_count = 0
        #keep a dictionary of multi-mapping guides for current file
        gRNA_match = dict()

        #search pattern
        gap = re.compile("-+")

        #make output dir
        outdir_path =config["gzdir"] + ".scored"
        if os.path.isdir(outdir_path):
            #shutil.rmtree(outdir_path)  # remove existing dir
            pass
        else:
            os.makedirs(outdir_path)

        #process gRNA ()

        file = os.path.join(config["gzdir"],config["gzfile"])      
        Chr = config["gzfile"].rstrip('.tab.mapped.gz').split(".")[0].split(r'.')[0]
        #print(Chr)
        file_gRNA_count = 0

        PAM_length = len(PAM)

        with gzip.open(file, "rt") as fh, gzip.open(os.path.join(outdir_path,f"{config['gzfile'].rstrip('.gz')}.scored.gz"), "wt") as wgfh:

            for line in fh:
                #print(line)
                #read gRNA from gzip file
                fields = line.split("\t")

                seq = fields[0]
                pam = fields[1]
                st = fields[2]
                en = fields[3]
                strand = fields[4].rstrip()
                leftFlank = fields[5]
                rightFlank = fields[6]
                gRNA_name = f"{Chr}_{st}_{en}_{strand}_{pam}"
                if strand == "+":
                    gRNA_name2 = f"{Chr}:{st}-{en}_{strand}"
                else:
                    gRNA_name2 = f"{Chr}:{en}-{st}_{strand}"  #adjust to match the format of bwa_mapping_coords (#bwa coords is always start <= end)

                #####################
                #calculate cfd score#
                #####################

                if len(fields[7])<=2:       # somes gRNAs does not have a bwa mapping, for example:  nnnnnnnnnnnnnngattca      tgg     257653  257675  + 
                    bwa_mapping = []
                else:
                    try:
                        bwa_mapping = ast.literal_eval(fields[7])
                    except (ValueError, SyntaxError):
                        # Handle the error or assign an empty list or a default value if the conversion fails
                        bwa_mapping = []

                #get bwa_mapping seqs
                bwa_mapping_seqs = [i.split('|')[2] for i in bwa_mapping]
                bwa_mapping_coords = [i.split('|')[0] for i in bwa_mapping]   #bwa coords is always start <= end
                
                #remove on-target seq from bwa_mapping
                if (gRNA_name2 in bwa_mapping_coords):
                    idx = bwa_mapping_coords.index(gRNA_name2)
                    del(bwa_mapping_coords[idx])
                    del(bwa_mapping_seqs[idx])
                else:
                    print(f"warning: bwa hits does not contain the ontarget of {gRNA_name2}, bwa hit coords: {bwa_mapping_coords}")

                if PAM in ["NGG", "NGA"]:
                    # MIT score must *not* include the PAM
                    guideNoPam = seq.upper()
                    pam_end_check = pam[-2:].upper() == "GG"
                    
                    mit_Ot_scores = [
                        (lambda mitScore: mitScore * 0.2 if pam_end_check and otSeq[-2:].upper() != "GG" else mitScore)(
                            calcHitScore(guideNoPam, otSeq[:-PAM_length].upper())
                        )
                        for otSeq in bwa_mapping_seqs
                        ]
                    
                    
                    # CFD score must include the PAM
                    # cfd_Ot_scores = [calcCfdScore(f"{seq}{pam}", otseq) for otseq in bwa_mapping_seqs]

                    #print out the cfd scores together with alignments
                    # print(gRNA_name2)
                    # for idx, otseq in enumerate(bwa_mapping_seqs):
                    #     score = calcCfdScore(f"{seq}{pam}", otseq)
                    #     coord = bwa_mapping_coords[idx]
                    #     print(f"{seq}{pam}\t{gRNA_name2}\n{otseq}\t{coord} \tCFD: {score}\n")
                    
                    #calc MIT specificity score for the current guide
                    mit_Ot_scores = [s for s in mit_Ot_scores if s]  #remove entries where cfd score failed to calcuate. For example, there is an offending "N" in NGTGGTGGTGGTAGAGGTGGTGG at position 6:167591372-167591394_-
                    guideMITScore = calcMitGuideScore(sum(mit_Ot_scores))

                    #calc CFD specificity score for the current guide
                    #cfd_Ot_scores = [s for s in cfd_Ot_scores if s]  #remove entries where cfd score failed to calcuate. For example, there is an offending "N" in NGTGGTGGTGGTAGAGGTGGTGG at position 6:167591372-167591394_-
                    #guideCfdScore = calcMitGuideScore(sum(cfd_Ot_scores))

                    #calc CFD specificity score v2
                    #guideCfdScorev2 = calcMitGuideScore_v2(sum(cfd_Ot_scores))

                    #calc CFD specificity score v3
                    #cfd_Ot_score_mod = [s if s<1 else 100 for s in cfd_Ot_scores]
                    #guideCfdScorev3 = calcMitGuideScore(sum(cfd_Ot_score_mod))

                    #print("seq:{}\tpam:{}\tscore:{}\tPAM:{}".format(seq,pam,guideMITScore,PAM))
                
                elif PAM in ["TTTV"]:

                    cfd_Ot_scores = [calcCfdScore_cas12a(seq, otseq[PAM_length:]) for otseq in bwa_mapping_seqs] # CFD score for TTTV PAM must not include PAM
                    cfd_Ot_scores = [s*100 for s in cfd_Ot_scores if s] 
                    guideCFDscore = calcMitGuideScore(sum(cfd_Ot_scores))

                    #print("seq:{}\tpam:{}\tscore:{}\tPAM:{}".format(seq,pam,guideCFDscore,PAM))
                
                
                #print(guideCfdScore)
                ######################
                #calculate eff scores#
                ######################
                scores_flattened=""

                if not config["skip_eff_score"]:
                    #input for calcAllScores() is  100bp sequences (50bp 5' of PAM, 50bp 3' of PAM) calculate all efficiency scores
                    len_L = 50 - len(seq)
                    len_R = 50 - len(pam)
                    gRNA_with_flank = leftFlank[-len_L:] + seq + pam + rightFlank[:len_R]
                    #print(len(gRNA_with_flank))
                    #test 
                    #print(sorted(calcAllScores(["CCACGTCTCCACACATCAGCACAACTACGCAGCGCCTCCCTCCACTCGGAAGGACTATCCTGCTGCCAAGAGGGTCAAGTTGGACAGTGTCAGAGTCCTG"]).items()))
                    m = re.search(gap, gRNA_with_flank)
                    
                    if m: #gap exists in gRNA_with_flank
                        scores_flattened = "|N.A. due to gap(s)|"
                    else: #no gaps, safe to calculate efficiency scores
                        scores = calcAllScores([gRNA_with_flank])
                        scores_flattened = flatten_score_dict(scores)
                    #print(scores_flattened)

                ################################################
                #write to scores (with gRNA info) to a new file#
                ################################################
                #print(f"{seq}\t{pam}\t{st}\t{en}\t{strand}\t{guideMITScore}\t{guideCfdScore}\t{guideCfdScorev2}\t{guideCfdScorev3}\t{scores_flattened}\n")
                if PAM in ["NGG", "NGA"]:
                    wgfh.write(f"{seq}\t{pam}\t{st}\t{en}\t{strand}\t{guideMITScore}\t{scores_flattened}\n")
                elif PAM in ["TTTV"]:
                    wgfh.write(f"{seq}\t{pam}\t{st}\t{en}\t{strand}\t{guideCFDscore}\t{scores_flattened}\n")
                
                gRNA_count += 1

                if gRNA_count%1000 == 0:
                    endtime = datetime.datetime.now()
                    elapsed_sec = endtime - starttime
                    elapsed_min = elapsed_sec.seconds / 60
                    print(f"elapsed {elapsed_min:.2f} min, scored {gRNA_count} gRNAs in file {file}", flush=True)

        endtime = datetime.datetime.now()
        elapsed_sec = endtime - starttime
        elapsed_min = elapsed_sec.seconds / 60
        log.info(f"finished in {elapsed_min:.2f} min, scored {gRNA_count} gRNAs in file {file}")
        print(f"finished in {elapsed_min:.2f} min, scored {gRNA_count} gRNAs in file {file}", flush=True)

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

def calcCfdScore_cas12a(guideSeq, otSeq):
    """ PAM sequence must not be included
    >>> calcCfdScore("GGGGGGGGGGGGGGGGGGGGGGG", "GGGGGGGGGGGGGGGGGAAAGGG")
    0.4635989007074176
    >>> calcCfdScore("GGGGGGGGGGGGGGGGGGGGGGG", "GGGGGGGGGGGGGGGGGGGGGGG")
    1.0
    >>> calcCfdScore("GGGGGGGGGGGGGGGGGGGGGGG", "aaaaGaGaGGGGGGGGGGGGGGG")
    0.5140384614450001
    # mismatches:      *               !!
    >>> calcCfdScore("ATGGTCGGACTCCCTGCCAGAGG", "ATGGTGGGACTCCCTGCCAGAGG")
    0.5

    # mismatches:    *  ** *
    >>> calcCfdScore("ATGGTCGGACTCCCTGCCAGAGG", "ATGATCCAAATCCCTGCCAGAGG")
    0.53625000020625

    >>> calcCfdScore("ATGTGGAGATTGCCACCTACCGG", "ATCTGGAGATTGCCACCTACAGG")

    """
    wt = guideSeq.upper()
    off = otSeq.upper()
    m_wt = re.search("[^ATCG]", wt)
    m_off = re.search("[^ATCG]", off)
    if (m_wt is None) and (m_off is None):
        sg = off[:20]
        wt = wt[:20]
        cfd_score = calc_cfd_cas12a(wt, sg)
        return cfd_score
    
def calc_cfd_cas12a(wt, sg):
    score = 1
    sg = sg.replace("T", "U")
    wt = wt.replace("T", "U")
    s_list = list(sg)
    wt_list = list(wt)
    for i, sl in enumerate(s_list):
        if wt_list[i] == sl:
            score *= 1
        else:
            key = "r" + wt_list[i] + ":d" + revcom(sl) + "," + str(i + 1)
            score *= mm_scores_cas12a[key]
    return score

def revcom(s):
    basecomp = {"A": "T", "C": "G", "G": "C", "T": "A", "U": "A"}
    letters = list(s[::-1])
    letters = [basecomp[base] for base in letters]
    return "".join(letters)

def calcHitScore(string1, string2):
    " see 'Scores of single hits' on http://crispr.mit.edu/about "
    # The Patrick Hsu weighting scheme
    # S. aureus requires 21bp long guides. We fudge by using only last 20bp
    matrixStart = 0
    maxDist = 19

    assert string1[0].isupper()
    assert len(string1) == len(string2)
    # for nmCas9 and a few others with longer guides, we limit ourselves to 20bp
    if len(string1) > 20:
        string1 = string1[-20:]
        string2 = string2[-20:]
    # for 19bp guides, we fudge a little, but first pos has no weight anyways
    elif len(string1) == 19:
        string1 = "A" + string1
        string2 = "A" + string2
    # for shorter guides, I'm not sure if this score makes sense anymore, we force things
    elif len(string1) < 19:
        matrixStart = 20 - len(string1)
        maxDist = len(string1) - 1

    assert len(string1) == len(string2)

    dists = []  # distances between mismatches, for part 2
    mmCount = 0  # number of mismatches, for part 3
    lastMmPos = None  # position of last mismatch, used to calculate distance

    score1 = 1.0
    for pos in range(matrixStart, len(string1)):
        if string1[pos] != string2[pos]:
            mmCount += 1
            if lastMmPos != None:
                dists.append(pos - lastMmPos)
            score1 *= 1 - hitScoreM[pos]
            lastMmPos = pos
    # 2nd part of the score
    if mmCount < 2:  # special case, not shown in the paper
        score2 = 1.0
    else:
        avgDist = sum(dists) / len(dists)
        score2 = 1.0 / (((maxDist - avgDist) / float(maxDist)) * 4 + 1)
    # 3rd part of the score
    if mmCount == 0:  # special case, not shown in the paper
        score3 = 1.0
    else:
        score3 = 1.0 / (mmCount ** 2)
    score = score1 * score2 * score3 * 100
    return score

def calcMitGuideScore(hitSum):
    """ Sguide defined on http://crispr.mit.edu/about
    Input is the sum of all off-target hit scores. Returns the specificity of the guide.
    """
    score = 100 / (100 + hitSum)
    score = int(round(score * 100))
    return score

def flatten_score_dict(score_dict):
    score_string = ""
    if 'fusi' in score_dict.keys():
        score_string = score_string + "|" + str(score_dict['fusi'][0])
    else:
        score_string = score_string + "|" + "NA"
    if 'fusiOld' in score_dict.keys():
        score_string = score_string + "|" + str(score_dict['fusiOld'][0])
    else:
        score_string = score_string + "|" + "NA"
    if 'chariRank' in score_dict.keys():
        score_string = score_string + "|" + str(score_dict['chariRank'][0])
    else:
        score_string = score_string + "|" + "NA"
    if 'ssc' in score_dict.keys():
        score_string = score_string + "|" + str(score_dict['ssc'][0])
    else:
        score_string = score_string + "|" + "NA"
    if 'wuCrispr' in score_dict.keys():
        score_string = score_string + "|" + str(score_dict['wuCrispr'][0])
    else:
        score_string = score_string + "|" + "NA"
    if 'doench' in score_dict.keys():
        score_string = score_string + "|" + str(score_dict['doench'][0])
    else:
        score_string = score_string + "|" + "NA"
    if 'wang' in score_dict.keys():
        score_string = score_string + "|" + str(score_dict['wang'][0])
    else:
        score_string = score_string + "|" + "NA"
    if 'crisprScan' in score_dict.keys(): #Moreno-Mateos
        score_string = score_string + "|" + str(score_dict['crisprScan'][0])
    else:
        score_string = score_string + "|" + "NA"
    if 'aziInVitro' in score_dict.keys():
        score_string = score_string + "|" + str(score_dict['aziInVitro'][0])
    else:
        score_string = score_string + "|" + "NA"
    if 'ccTop' in score_dict.keys():
        score_string = score_string + "|" + str(score_dict['ccTop'][0])
    else:
        score_string = score_string + "|" + "NA"
    if 'chariRaw' in score_dict.keys():
        score_string = score_string + "|" + str(score_dict['chariRaw'][0])
    else:
        score_string = score_string + "|" + "NA"
    if 'housden' in score_dict.keys():
        score_string = score_string + "|" + str(score_dict['housden'][0])
    else:
        score_string = score_string + "|" + "NA"                

    return score_string

if __name__ == "__main__": main()

