from Bio.Seq import Seq
import os.path
import pandas as pd
from Bio.Seq import reverse_complement
import argparse
import sys
import math
import pickle



import logging
#################
# custom logging #
#################
from protoSpaceJAM.util.hdr import HDR_flank


class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write("error: %s\n" % message)
        self.print_help()
        sys.exit(2)


BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE = range(8)
# The background is set with 40 plus the number of the color, and the foreground with 30

# These are the sequences need to get colored ouput
RESET_SEQ = "\033[0m"
COLOR_SEQ = "\033[1;%dm"
BOLD_SEQ = "\033[1m"


def formatter_message(message, use_color=True):
    if use_color:
        message = message.replace("$RESET", RESET_SEQ).replace("$BOLD", BOLD_SEQ)
    else:
        message = message.replace("$RESET", "").replace("$BOLD", "")
    return message


COLORS = {
    "WARNING": YELLOW,
    "INFO": WHITE,
    "DEBUG": BLUE,
    "CRITICAL": YELLOW,
    "ERROR": RED,
}


class ColoredFormatter(logging.Formatter):
    def __init__(self, msg, use_color=True):
        logging.Formatter.__init__(self, msg)
        self.use_color = use_color

    def format(self, record):
        levelname = record.levelname
        if self.use_color and levelname in COLORS:
            levelname_color = (
                COLOR_SEQ % (30 + COLORS[levelname]) + levelname + RESET_SEQ
            )
            record.levelname = levelname_color
        return logging.Formatter.format(self, record)


# Custom logger class with multiple destinations
class ColoredLogger(logging.Logger):
    FORMAT = "[$BOLD%(name)-1s$RESET][%(levelname)-1s]  %(message)s "  # ($BOLD%(filename)s$RESET:%(lineno)d)
    COLOR_FORMAT = formatter_message(FORMAT, True)

    def __init__(self, name):
        logging.Logger.__init__(self, name, logging.DEBUG)

        color_formatter = ColoredFormatter(self.COLOR_FORMAT)

        console = logging.StreamHandler(stream=sys.stdout)
        console.setFormatter(color_formatter)

        self.addHandler(console)
        return


logging.setLoggerClass(ColoredLogger)
# logging.basicConfig()
log = logging.getLogger("ProtospaceX")
log.propagate = False
log.setLevel(logging.INFO)


def get_phase_in_codon(Chr, Pos, ENST_ID, ENST_PhaseInCodon):
    """
    get the phase in codon for the current position and flanking 2 bp
    returns a dictionary:
    {0:[-1/-2/-3/1/2/3/0] 0 stands for not coding sequence (current position is not in a codon)
    -1:[-1/-2/-3/1/2/3/0]
    -2:[-1/-2/-3/1/2/3/0]
    +1:[-1/-2/-3/1/2/3/0]
    +2:[-1/-2/-3/1/2/3/0]}
    """
    mydict = {-2: 0, -1: 0, 0: 0, 1: 0, 2: 0}
    if Chr in ENST_PhaseInCodon.keys():
        Chr_dict = ENST_PhaseInCodon[Chr]
        if Pos in Chr_dict.keys():
            if ENST_ID in Chr_dict[Pos].keys():
                mydict[0] = Chr_dict[Pos][ENST_ID]
        if Pos + 1 in Chr_dict.keys():
            if ENST_ID in Chr_dict[Pos + 1].keys():
                mydict[+1] = Chr_dict[Pos + 1][ENST_ID]
        if Pos + 2 in Chr_dict.keys():
            if ENST_ID in Chr_dict[Pos + 2].keys():
                mydict[+2] = Chr_dict[Pos + 2][ENST_ID]
        if Pos - 1 in Chr_dict.keys():
            if ENST_ID in Chr_dict[Pos - 1].keys():
                mydict[-1] = Chr_dict[Pos - 1][ENST_ID]
        if Pos - 2 in Chr_dict.keys():
            if ENST_ID in Chr_dict[Pos - 2].keys():
                mydict[-2] = Chr_dict[Pos - 2][ENST_ID]
    return mydict


def get_phase_in_codon0(Chr, Pos, ENST_ID, ENST_PhaseInCodon):
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


def get_range(start, end):  # TODO phase should be ENST specific
    """
    return a list of number, start to end, step size = 1, -1 (if start > end)
    """
    if start <= end:
        return list(range(start, end + 1, 1))
    else:
        return list(range(start, end - 1, -1))


def get_HDR_template(
    df,
    ENST_info,
    type,
    ENST_PhaseInCodon,
    HDR_arm_len,
    genome_ver,
    tag,
    loc2posType,
    ssODN_max_size,
    recoding_args,
    Donor_type,
    Strand_choice,
    syn_check_args,
):
    for index, row in df.iterrows():
        ENST_ID = row["ID"]
        ENST_genename = ENST_info[ENST_ID].name
        ENST_strand = ENST_info[ENST_ID].features[0].strand
        Chr = row["chr"]
        InsPos = row[
            "Insert_pos"
        ]  # InsPos is the first letter of stop codon "T"AA or the last letter of the start codon AT"G"
        gStart = row["start"]
        gStrand = convert_strand(row["strand"])  # gRNA strand
        CutPos = get_cut_pos(gStart, gStrand)
        Cut2Ins_dist = row["Cut2Ins_dist"]

        ##########################
        # important debug info
        # print(f"ENST: {ENST_ID} ENST_strand: {ENST_strand} chr {Chr} type {type}-tagging InsPos {InsPos} gStrand {gStrand} CutPos {CutPos} gStart {gStart} ")
        ##########################
        # CutPos_phase = get_phase_in_codon(Chr=Chr, Pos=InsPos, ENST_PhaseInCodon=ENST_PhaseInCodon)
        # print(f"CutPos phase: {CutPos_phase}")170

        # get target_seq
        # target_seq = get_target_seq(Chr= Chr, InsPos = InsPos, gRNAstrand = gStrand, CutPos = CutPos , type = type, ENST_ID = ENST_ID, ENST_info = ENST_info)

        leftArm, rightArm, left_start, left_end, right_start, right_end = get_HDR_arms(
            loc=[Chr, InsPos, ENST_strand],
            half_len=HDR_arm_len,
            type=type,
            genome_ver=genome_ver,
        )  # start>end is possible
        left_Arm_Phases = [
            get_phase_in_codon0(
                Chr=Chr, Pos=i, ENST_ID=ENST_ID, ENST_PhaseInCodon=ENST_PhaseInCodon
            )
            for i in get_range(left_start, left_end)
        ]
        right_Arm_Phases = [
            get_phase_in_codon0(
                Chr=Chr, Pos=i, ENST_ID=ENST_ID, ENST_PhaseInCodon=ENST_PhaseInCodon
            )
            for i in get_range(right_start, right_end)
        ]

        myflank = HDR_flank(
            left_flk_seq=leftArm,
            right_flk_seq=rightArm,
            left_flk_coord_lst=[left_start, left_end],
            right_flk_coord_lst=[right_start, right_end],
            left_flk_phases=left_Arm_Phases,
            right_flk_phases=right_Arm_Phases,
            type=type,
            ENST_ID=ENST_ID,
            name=ENST_genename,
            ENST_strand=ENST_strand,
            ENST_chr=Chr,
            gStart=gStart,
            gStrand=gStrand,
            InsPos=InsPos,
            CutPos=CutPos,
            Cut2Ins_dist=Cut2Ins_dist,
            tag=tag,
            loc2posType=loc2posType,
            ssODN_max_size=ssODN_max_size,
            Donor_type=Donor_type,
            Strand_choice=Strand_choice,
            recoding_args=recoding_args,
            syn_check_args=syn_check_args,
        )
        return myflank
        # log IDs whose gRNA is not in the default-size HDR arms
        # if myflank.entire_gRNA_in_HDR_arms == False:
        #    gRNA_out_of_arms[type][ENST_ID] = False


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
    # convert strand from +/- to 1/-1
    if strand == "+":
        return 1
    elif strand == "-":
        return -1
    else:
        return f"input strand:{strand} needs to be +/-"


def get_HDR_arms(loc, half_len, type, genome_ver):
    """
    input:  loc         [chr,pos,strand]  #start < end , strand is the coding strand
            half_len      length of the HDR arm (one sided)
    return: HDR arms -> in coding strand <-
            [5'arm, 3'arm]
    """
    Chr, Pos, Strand = loc
    # get arms
    vanilla_left_arm, vanilla_right_arm = None, None
    if type == "start":
        if Strand == 1:
            vanilla_left_arm = get_seq(
                chr=Chr,
                start=Pos - half_len + 1,
                end=Pos + 1,
                strand=1,
                genome_ver=genome_ver,
            )
            vanilla_right_arm = get_seq(
                chr=Chr,
                start=Pos + 1,
                end=Pos + half_len + 1,
                strand=1,
                genome_ver=genome_ver,
            )
            return [
                vanilla_left_arm,
                vanilla_right_arm,
                Pos - half_len + 1,
                Pos,
                Pos + 1,
                Pos + half_len,
            ]
        elif Strand == -1:
            vanilla_left_arm = get_seq(
                chr=Chr, start=Pos - half_len, end=Pos, strand=1, genome_ver=genome_ver
            )
            vanilla_right_arm = get_seq(
                chr=Chr, start=Pos, end=Pos + half_len, strand=1, genome_ver=genome_ver
            )
            return [
                reverse_complement(vanilla_right_arm),
                reverse_complement(vanilla_left_arm),
                Pos + half_len - 1,
                Pos,
                Pos - 1,
                Pos - half_len,
            ]
        else:
            sys.exit(f"unknown strand: {Strand}, acceptable values are -1 and 1")
    elif type == "stop":
        if Strand == 1:
            vanilla_left_arm = get_seq(
                chr=Chr, start=Pos - half_len, end=Pos, strand=1, genome_ver=genome_ver
            )
            vanilla_right_arm = get_seq(
                chr=Chr, start=Pos, end=Pos + half_len, strand=1, genome_ver=genome_ver
            )
            return [
                vanilla_left_arm,
                vanilla_right_arm,
                Pos - half_len,
                Pos - 1,
                Pos,
                Pos + half_len - 1,
            ]
        elif Strand == -1:
            vanilla_left_arm = get_seq(
                chr=Chr,
                start=Pos - half_len + 1,
                end=Pos + 1,
                strand=1,
                genome_ver=genome_ver,
            )
            vanilla_right_arm = get_seq(
                chr=Chr,
                start=Pos + 1,
                end=Pos + half_len + 1,
                strand=1,
                genome_ver=genome_ver,
            )
            return [
                reverse_complement(vanilla_right_arm),
                reverse_complement(vanilla_left_arm),
                Pos + half_len,
                Pos + 1,
                Pos,
                Pos - half_len + 1,
            ]
    else:
        sys.exit("unknown type {type}, acceptable values: start, stop")

    # rev_com ajustment
    if Strand == 1:
        return [vanilla_left_arm, vanilla_right_arm]
    elif Strand == -1:
        return [
            reverse_complement(vanilla_right_arm),
            reverse_complement(vanilla_left_arm),
        ]
    else:
        sys.exit(f"unknown strand: {Strand}, acceptable values are -1 and 1")


def get_gRNAs(
    ENST_ID,
    ENST_info,
    freq_dict,
    loc2file_index,
    loc2posType,
    genome_ver,
    spec_score_flavor,
    dist=50,
):
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
    # get ENST strand
    ENST_strand = ENST_info[ENST_ID].features[0].strand
    ##################################
    # get gRNAs around the start codon#
    ##################################
    # get the start codon chromosomal location
    log.debug(f"ATG_loc: {ATG_loc}")
    end_of_ATG_loc = get_end_pos_of_ATG(ATG_loc)  # [chr, pos,strand]
    log.debug(f"end of the ATG: {end_of_ATG_loc}")
    # get gRNA around the chromosomeal location (near ATG)
    df_gRNAs_ATG = get_gRNAs_near_loc(
        loc=end_of_ATG_loc,
        dist=dist,
        loc2file_index=loc2file_index,
        genome_ver=genome_ver,
    )
    # rank gRNAs
    ranked_df_gRNAs_ATG = rank_gRNAs_for_tagging(
        loc=end_of_ATG_loc,
        gRNA_df=df_gRNAs_ATG,
        loc2posType=loc2posType,
        ENST_ID=ENST_ID,
        ENST_strand=ENST_strand,
        type="start",
        spec_score_flavor=spec_score_flavor,
    )
    ranked_df_gRNAs_ATG = ranked_df_gRNAs_ATG.sort_values(
        "final_weight", ascending=False
    )  # sort descending on final weight

    ##################################
    # get gRNAs around the stop  codon#
    ##################################
    # get gRNAs around the stop codon
    log.debug(f"stop_loc: {stop_loc}")
    start_of_stop_loc = get_start_pos_of_stop(stop_loc)
    log.debug(f"start of stop: {start_of_stop_loc}")  # [chr, pos,strand]
    # get gRNA around the chromosomeal location (near stop location)
    df_gRNAs_stop = get_gRNAs_near_loc(
        loc=start_of_stop_loc,
        dist=dist,
        loc2file_index=loc2file_index,
        genome_ver=genome_ver,
    )
    # rank gRNAs
    ranked_df_gRNAs_stop = rank_gRNAs_for_tagging(
        loc=start_of_stop_loc,
        gRNA_df=df_gRNAs_stop,
        loc2posType=loc2posType,
        ENST_ID=ENST_ID,
        ENST_strand=ENST_strand,
        type="stop",
        spec_score_flavor=spec_score_flavor,
    )
    ranked_df_gRNAs_stop = ranked_df_gRNAs_stop.sort_values(
        "final_weight", ascending=False
    )  # sort descending on final weight

    return [ranked_df_gRNAs_ATG, ranked_df_gRNAs_stop]


def rank_gRNAs_for_tagging(
    loc, gRNA_df, loc2posType, ENST_ID, ENST_strand, type, spec_score_flavor, alpha=1
):
    """
    input:  loc         [chr,pos,strand]  #start < end
            gRNA_df     pandas dataframe, *unranked*   columns: "seq","pam","start","end", "strand", "guideMITScore","guideCfdScore","guideCfdScorev2","guideCfdScorev3", "Eff_scores"  !! neg strand: start > end
            alpha       specificity weight is raised to the power of alpha
            type        "start" or "stop:
    output: gRNA_df     pandas dataframe *ranked*      columns: "seq","pam","start","end", "strand", "guideMITScore","guideCfdScore","guideCfdScorev2","guideCfdScorev3", "Eff_scores"  !! neg strand: start > end
    """
    insPos = loc[
        1
    ]  # InsPos is the first letter of stop codon "T"AA or the last letter of the start codon AT"G"
    Chr = loc[0]

    col_spec_weight = []
    col_dist_weight = []
    col_pos_weight = []
    col_final_weight = []
    Chrs = []
    ENSTs = []
    InsertPos = []
    cut2insDist_list = []
    # assign a score to each gRNA
    for index, row in gRNA_df.iterrows():
        start = row[2]
        end = row[3]
        strand = row[4]
        # Get cut to insert distance
        cutPos = get_cut_pos(start, strand)
        cut2insDist = cutPos - insPos

        # adjust cut2insDist
        if type == "start" and ENST_strand == -1:
            cut2insDist += 1
        if type == "stop" and ENST_strand == 1:
            cut2insDist += 1

        # calc. specificity_weight
        CSS = row[spec_score_flavor]
        specificity_weight = _specificity_weight(CSS)
        col_spec_weight.append(specificity_weight)

        # calc. distance_weight
        distance_weight = _dist_weight(hdr_dist=cut2insDist)
        col_dist_weight.append(distance_weight)

        # get position_weight
        position_type = _get_position_type(
            chr=Chr, ID=ENST_ID, pos=cutPos, loc2posType=loc2posType
        )
        position_weight = _position_weight(position_type)

        position_type_nextbp = _get_position_type(
            chr=Chr, ID=ENST_ID, pos=cutPos + 1, loc2posType=loc2posType
        )
        position_weight_nextbp = _position_weight(position_type_nextbp)

        position_weight = min([position_weight, position_weight_nextbp])

        col_pos_weight.append(position_weight)

        # add info to the df
        Chrs.append(Chr)
        ENSTs.append(ENST_ID)
        InsertPos.append(insPos)
        cut2insDist_list.append(cut2insDist)

        log.debug(
            f"strand {strand} {start}-{end} cutPos {cutPos} insert_loc {loc} cut2insDist {cut2insDist} distance_weight {distance_weight:.2f} CFD_score {CSS} specificity_weight {specificity_weight} pos_type {position_type} position_weight {position_weight}"
        )

        final_score = (
            float(pow(specificity_weight, alpha))
            * float(distance_weight)
            * float(position_weight)
        )
        col_final_weight.append(final_score)

    # add info and weight columns to the df
    gRNA_df["chr"] = Chrs
    gRNA_df["ID"] = ENSTs
    gRNA_df["Insert_pos"] = InsertPos
    gRNA_df["Cut2Ins_dist"] = cut2insDist_list
    gRNA_df["spec_weight"] = col_spec_weight
    gRNA_df["dist_weight"] = col_dist_weight
    gRNA_df["pos_weight"] = col_pos_weight
    gRNA_df["final_weight"] = col_final_weight

    # rank gRNAs based on the score
    gRNA_df["final_pct_rank"] = gRNA_df["final_weight"].rank(pct=True)

    return gRNA_df


def get_cut_pos(start, strand):
    """
    start:gRNA start
    strand:gRNA strand
    """
    if strand == "+" or strand == "1" or strand == 1:
        cutPos = start + 16
    else:
        cutPos = (
            start - 16 - 1
        )  # -1 because we want the cutsite to be behind the pos (viewed in the +1 strand)
    return cutPos


def _get_position_type(chr, ID, pos, loc2posType):
    """
    #input: mostly self-explanatory, loc2posType is a dictionary that translates location into types (e.g. exon/intron junctions etc)
    #return a list of types for the input position/ID combination
    >>> loc2posType = read_pickle_files(os.path.join("..","genome_files","parsed_gff3", "GRCh38","loc2posType.pickle"))

    #intron - exon junction
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50399367, loc2posType = loc2posType) #-4bp
    >>> print(postype)
    ['3N4bp_up_of_intron_exon_junction']
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50399368, loc2posType = loc2posType) #-3bp
    >>> print(postype)
    ['within_3bp_of_intron_exon_junction', '3N4bp_up_of_intron_exon_junction']
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50399369, loc2posType = loc2posType) #-2bp
    >>> print(postype)
    ['within_2bp_of_intron_exon_junction', 'within_3bp_of_intron_exon_junction']
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50399370, loc2posType = loc2posType) #-1bp
    >>> print(postype)
    ['within_2bp_of_intron_exon_junction', 'within_3bp_of_intron_exon_junction']
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50399371, loc2posType = loc2posType) #1bp
    >>> print(postype)
    ['cds', 'within_2bp_of_intron_exon_junction', 'within_3bp_of_intron_exon_junction']
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50399372, loc2posType = loc2posType) #2bp
    >>> print(postype)
    ['cds', 'within_2bp_of_intron_exon_junction', 'within_3bp_of_intron_exon_junction']
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50399373, loc2posType = loc2posType) #3bp
    >>> print(postype)
    ['cds', 'within_3bp_of_intron_exon_junction', '3N4bp_down_of_intron_exon_junction']
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50399374, loc2posType = loc2posType) #4bp
    >>> print(postype)
    ['cds', '3N4bp_down_of_intron_exon_junction']

    #intron - exon junction
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50399481, loc2posType = loc2posType) #-4bp
    >>> print(postype)
    ['cds', '3N4bp_up_of_exon_intron_junction']
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50399482, loc2posType = loc2posType) #-3bp
    >>> print(postype)
    ['cds', 'within_3bp_of_exon_intron_junction', '3N4bp_up_of_exon_intron_junction']
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50399483, loc2posType = loc2posType) #-2bp
    >>> print(postype)
    ['cds', 'within_2bp_of_exon_intron_junction', 'within_3bp_of_exon_intron_junction']
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50399484, loc2posType = loc2posType) #-1bp
    >>> print(postype)
    ['cds', 'within_2bp_of_exon_intron_junction', 'within_3bp_of_exon_intron_junction']
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50399485, loc2posType = loc2posType) #1bp
    >>> print(postype)
    ['within_2bp_of_exon_intron_junction', 'within_3bp_of_exon_intron_junction']
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50399486, loc2posType = loc2posType) #2bp
    >>> print(postype)
    ['within_2bp_of_exon_intron_junction', 'within_3bp_of_exon_intron_junction']
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50399487, loc2posType = loc2posType) #3bp
    >>> print(postype)
    ['within_3bp_of_exon_intron_junction', '3_to_6bp_down_of_exon_intron_junction']
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50399488, loc2posType = loc2posType) #4bp
    >>> print(postype)
    ['3_to_6bp_down_of_exon_intron_junction']
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50399489, loc2posType = loc2posType) #5bp
    >>> print(postype)
    ['3_to_6bp_down_of_exon_intron_junction']
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50399490, loc2posType = loc2posType) #6bp
    >>> print(postype)
    ['3_to_6bp_down_of_exon_intron_junction']

    #non cds exon
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50398847, loc2posType = loc2posType) #-4bp
    >>> print(postype)
    ['3N4bp_up_of_intron_exon_junction']
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50398848, loc2posType = loc2posType) #-3bp
    >>> print(postype)
    ['within_3bp_of_intron_exon_junction', '3N4bp_up_of_intron_exon_junction']
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50398849, loc2posType = loc2posType) #-2bp
    >>> print(postype)
    ['within_2bp_of_intron_exon_junction', 'within_3bp_of_intron_exon_junction']
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50398850, loc2posType = loc2posType) #-1bp
    >>> print(postype)
    ['within_2bp_of_intron_exon_junction', 'within_3bp_of_intron_exon_junction']
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50398851, loc2posType = loc2posType) #1bp
    >>> print(postype)
    ['5UTR', 'within_2bp_of_intron_exon_junction', 'within_3bp_of_intron_exon_junction']
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50398852, loc2posType = loc2posType) #2bp
    >>> print(postype)
    ['cds', 'within_2bp_of_intron_exon_junction', 'within_3bp_of_intron_exon_junction']
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50398853, loc2posType = loc2posType) #3bp
    >>> print(postype)
    ['cds', 'within_3bp_of_intron_exon_junction', '3N4bp_down_of_intron_exon_junction']
    >>> postype = _get_position_type(chr="19", ID="ENST00000440232", pos=50398854, loc2posType = loc2posType) #4bp
    >>> print(postype)
    ['cds', '3N4bp_down_of_intron_exon_junction']
    """
    # print(f"{type(chr)} {ID} {type(pos)}")
    chr_dict = loc2posType[chr]
    if not ID in chr_dict.keys():
        return []
    else:
        types = []
        mapping_dict = chr_dict[ID]
        for key in mapping_dict.keys():
            if in_interval_leftrightInclusive(pos, key):
                types.append(mapping_dict[key])
        return types


def _position_weight(types):
    """
    input: a list of types
    output: the lowest weight among all the types
    """
    mapping_dict = {
        "5UTR": 0.4,
        "3UTR": 1,
        "cds": 1,
        "within_2bp_of_exon_intron_junction": 0.01,
        "within_2bp_of_intron_exon_junction": 0.01,
        "3N4bp_up_of_exon_intron_junction": 0.1,
        "3_to_6bp_down_of_exon_intron_junction": 0.1,
        "3N4bp_up_of_intron_exon_junction": 0.1,
        "3N4bp_down_of_intron_exon_junction": 0.5,
    }
    lowest_weight = 1
    for type in types:
        if type in mapping_dict.keys():
            weight = mapping_dict[type]
            if weight < lowest_weight:
                lowest_weight = weight
        else:
            pass
            # sys.exit(f"unexpected position type: {type}")
    return lowest_weight


def _dist_weight(hdr_dist: int, _dist_weight_variance=55) -> float:
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
    weight = math.exp((-1 * hdr_dist ** 2) / (2 * variance))
    assert weight >= 0 and weight <= 1
    return weight


def _specificity_weight(
    specificity_score: float, _specificity_weight_low=45, _specificity_weight_high=65
):
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
        return [ATG_loc[0], ATG_loc[2], ATG_loc[3]]
    else:
        return [ATG_loc[0], ATG_loc[1], ATG_loc[3]]


def get_start_pos_of_stop(stop_loc):
    """
    input: stop_loc                                 [chr,start,end,strand]  #start < end
    output: the pos of first base in the stop codon [chr,pos,strand]        #start < end
    """
    strand = stop_loc[3]
    if str(strand) == "+" or str(strand) == "1":
        return [stop_loc[0], stop_loc[1], stop_loc[3]]
    else:
        return [stop_loc[0], stop_loc[2], stop_loc[3]]


def get_start_stop_loc(ENST_ID, ENST_info):
    """
    Get the chromosomal location of start and stop codons
    input: ENST_ID, ENST_info
    output: a list of three items
            ATG_loc: [chr,start,end,strand]  #start < end
            stop_loc: [chr,start,end,strand] #start < end
            Exon_end_ATG: Bool
    """
    my_transcript = ENST_info[ENST_ID]  # get the seq record
    # constructing the list of cds
    cdsList = [feat for feat in my_transcript.features if feat.type == "CDS"]
    CDS_first = cdsList[0]
    CDS_last = cdsList[len(cdsList) - 1]
    # check if ATG is at the end of the first exon
    if check_ATG_at_exonEnd(my_transcript):
        CDS_first = cdsList[
            1
        ]  # use the second cds if ATG is at the end of the first exon
    # get start codon location
    if CDS_first.strand == 1:
        ATG_loc = [
            CDS_first.location.ref,
            CDS_first.location.start + 0,
            CDS_first.location.start + 2,
            1,
        ]  # format [start, end, strand]
    else:
        stop_loc = [
            CDS_first.location.ref,
            CDS_first.location.start + 0,
            CDS_first.location.start + 2,
            -1,
        ]
    # get stop codon location
    if CDS_last.strand == 1:
        stop_loc = [
            CDS_last.location.ref,
            CDS_last.location.end - 2,
            CDS_last.location.end + 0,
            1,
        ]
    else:
        ATG_loc = [
            CDS_last.location.ref,
            CDS_last.location.end - 2,
            CDS_last.location.end + 0,
            -1,
        ]

    return [ATG_loc, stop_loc]


def get_gRNAs_near_loc(loc, dist, loc2file_index, genome_ver):
    """
    input
        loc: [chr,pos,strand]
        dist: max cut to loc distance
    return:
        a dataframe of guide RNAs which cuts <[dist] to the loc, the columns are "seq","pam","start","end", "strand", "guideMITScore","guideCfdScore","guideCfdScorev2","guideCfdScorev3", "Eff_scores"
    """
    chr = loc[0]
    pos = loc[1]
    chr_dict = loc2file_index[chr]
    target_files = []  # a list of file names containing gRNAs near loc
    # lookup the file
    for key in chr_dict.keys():
        file_start = key.split("-")[0]
        file_end = key.split("-")[1]
        interval = [file_start, file_end]
        if (
            in_interval(pos, interval)
            or in_interval(pos - 1000, interval)
            or in_interval(pos + 1000, interval)
        ):
            target_files.append(chr_dict[key])
    # print(target_files)
    # load the gRNAs in the file
    dfs = []
    for file in target_files:
        file_path = os.path.join(
            "precomputed_gRNAs",
            f"gRNA_{genome_ver}",
            "gRNA.tab.gz.split.BwaMapped.scored",
            file,
        )
        df_tmp = pd.read_csv(
            file_path,
            sep="\t",
            compression="infer",
            header=None,
            names=[
                "seq",
                "pam",
                "start",
                "end",
                "strand",
                "guideMITScore",
                "guideCfdScore",
                "guideCfdScorev2",
                "guideCfdScorev3",
                "Eff_scores",
            ],
        )
        dfs.append(df_tmp)
    df_gRNA = pd.concat(dfs)

    # subset gRNA based on strand  !ATTN: start > end when strand is '-'
    df_gRNA_on_sense = df_gRNA[(df_gRNA["strand"] == "+")]
    df_gRNA_on_antisense = df_gRNA[(df_gRNA["strand"] == "-")]

    # subset gRNAs and retain those cuts <[dist] to the loc
    df_gRNA_on_sense = df_gRNA_on_sense[
        (df_gRNA_on_sense["start"] > (pos - 17 - dist))
        & (df_gRNA_on_sense["start"] < (pos - 17 + dist))
    ]
    df_gRNA_on_antisense = df_gRNA_on_antisense[
        (df_gRNA_on_antisense["start"] > (pos + 17 - dist))
        & (df_gRNA_on_antisense["start"] < (pos + 17 + dist))
    ]

    return pd.concat([df_gRNA_on_sense, df_gRNA_on_antisense])


def in_interval(pos, interval):
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
        # log.debug(f"{pos} is in {interval}")
        return True
    else:
        return False


def in_interval_leftrightInclusive(pos, interval):
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
        # log.debug(f"{pos} is in {interval}")
        return True
    else:
        return False


def update_dict_count(
    key, dict
):  # update the dictionary that keeps the count of each string (key)
    if key in dict.keys():
        dict[key] += 1
    else:
        dict[key] = 1
    return dict


def get_seq(chr, start, end, strand, genome_ver):
    """
    chr
    start (1-indexed)
    end (the end position is not included
    strand 1 or -1 (str)
    """
    chr_file_path = os.path.join("genome_files", "fa_pickle", genome_ver, f"{chr}.pk")
    log.debug(f"opening file {chr_file_path}")
    if os.path.isfile(chr_file_path):
        # read file
        chr_seqrecord = read_pickle_files(chr_file_path)
        subseq = str(chr_seqrecord.seq)[
            (start - 1) : (end - 1)
        ]  # use -1 to convert 1-index to 0-index
        if strand == "-1" or strand == -1:
            return reverse_complement(subseq)
        else:
            return subseq
    else:
        sys.exit(f"ERROR: file not found: {chr_file_path}")


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


def count_ATG_at_exonEnd(ENST_info):
    """
    return a list of two items:
        count of number of ENST_IDs with ATG at the end of the exon
        list of such ENST_ID
    """
    count = 0
    list = []
    for ENST_ID in ENST_info.keys():
        my_transcript = ENST_info[ENST_ID]  # get the seq record
        if check_ATG_at_exonEnd(my_transcript):
            list.append(ENST_ID)
            count += 1
    return [count, list]


def check_ATG_at_exonEnd(my_transcript):
    """
    input: transcript object
    output: Bool
    """
    transcript_type = my_transcript.description.split("|")[1]
    if transcript_type == "protein_coding":  # only look at protein-coding transcripts
        # constructing the list of cds
        cdsList = [feat for feat in my_transcript.features if feat.type == "CDS"]
        if len(cdsList) >= 1:  # has more than 1 cds
            CDS_first = cdsList[0]
            cds_len = (
                abs(CDS_first.location.start - CDS_first.location.end) + 1
            )  # for ATG to be the end of the exon, the first exon length is 3bp
            if cds_len == 3:
                return True
            else:
                return False
        else:
            return False
    else:
        return False


def get_cds_seq_in_transcript(mytranscript):
    """
    input: Bio.SeqRecord  (such as that from function fetch_ensembl_transcript)
    return the cds sequence as a string
    """
    wholeSeq = str(mytranscript.seq)
    wholeSeq_rc = str(Seq(wholeSeq).reverse_complement())
    cds_seqs = []
    negative_strand_flag = False
    for feat in mytranscript.features:
        if feat.type == "cds":
            st = feat.location.start.position
            en = feat.location.end.position
            if feat.strand == -1:  # neg strand
                cds_seqs.append(
                    str(Seq(wholeSeq_rc[st:en]).reverse_complement())
                )  # the coord are respective to the revcom of the retrieved seq (weird)
                negative_strand_flag = True
            else:  # pos strand
                cds_seqs.append(wholeSeq[st:en])
    if negative_strand_flag == True:
        return "".join(cds_seqs[::-1])
    else:
        return "".join(cds_seqs)


def get_cds_seqNflank(transcriptObj, which_cds, cds_flank_len):
    """
    returns the n-th cds and its flanking sequences
    the returned sequence will be in the coding strand
    """
    total_cds_num = transcriptObj.num_cds
    if which_cds > total_cds_num:
        raise ValueError(
            f"Requested coding exon is {which_cds}, but the transcript only has {total_cds_num} coding exons"
        )

    wholeSeq = str(transcriptObj.seq)
    wholeSeq_rc = str(Seq(wholeSeq).reverse_complement())

    cds_list = [feat for feat in transcriptObj.features if feat.type == "cds"]
    strand = list(
        set([feat.strand for feat in transcriptObj.features if feat.type == "cds"])
    )[0]
    which_cds -= 1

    if strand == 1:
        which_cds = which_cds
    else:
        which_cds = len(cds_list) - which_cds - 1

    # get cds seq with flank
    target_cds = cds_list[which_cds]
    cds_st = target_cds.location.start.position
    cds_en = target_cds.location.end.position
    if strand == 1:
        cds_seq = wholeSeq[cds_st:cds_en]
        Lflank_st = max([i for i in [cds_st - cds_flank_len, 0] if i >= 0])
        RFlank_en = min(
            [i for i in [cds_en + cds_flank_len, len(wholeSeq)] if i <= len(wholeSeq)]
        )
        Lflank = wholeSeq[Lflank_st:cds_st]
        Rflank = wholeSeq[cds_en:RFlank_en]
    else:
        cds_seq = str(
            Seq(wholeSeq_rc[cds_st:cds_en]).reverse_complement()
        )  # the coord are respective to the revcom of the retrieved seq (weird)
        Lflank_st = max([i for i in [cds_st - cds_flank_len, 0] if i >= 0])
        RFlank_en = min(
            [
                i
                for i in [cds_en + cds_flank_len, len(wholeSeq_rc)]
                if i <= len(wholeSeq_rc)
            ]
        )
        Lflank = str(Seq(wholeSeq_rc[Lflank_st:cds_st]).reverse_complement())
        Rflank = str(Seq(wholeSeq_rc[cds_en:RFlank_en]).reverse_complement())
        Lflank, Rflank = Rflank, Lflank  # for -1 strand, the Lflank is the Rflank

    return {"cds_seq": cds_seq, "Lflank": Lflank, "Rflank": Rflank}


def get_exon_concat_with_flank(transcriptObj, which_cds, cds_flank_len):
    """
    get the nth cds segment from the transcript and concat with flanks
    the returned sequence will be in the coding strand
    """
    cds_seq_N_flank = get_cds_seqNflank(
        transcriptObj=transcriptObj, which_cds=which_cds, cds_flank_len=cds_flank_len
    )
    cds_with_flank = "".join(
        [
            cds_seq_N_flank["Lflank"],
            cds_seq_N_flank["cds_seq"],
            cds_seq_N_flank["Rflank"],
        ]
    )
    return cds_with_flank


def get_cutsite_in_gene(transcriptObj, listOfgRNAObj, which_cds, cds_flank_len):
    """
    find the coordinate (respective to the gene) of the cutsite in each gRNA

    Parameters
    ----------
    transcriptObj
    listOfgRNAObj
    which_cds: the i-th cds used to extract the gRNA
    cds_flank_len: the len of flank added to the cds when extracting gRNA sequences

    Returns
    -------
    list of gRNA object (updated with the cutsite)
    """
    cds_w_flank = get_exon_concat_with_flank(
        transcriptObj=transcriptObj, which_cds=which_cds, cds_flank_len=cds_flank_len
    )

    # find cds start and end
    wholeSeq = str(transcriptObj.seq)
    wholeSeq_rc = str(Seq(wholeSeq).reverse_complement())
    cds_list = [feat for feat in transcriptObj.features if feat.type == "cds"]
    cds_strand = list(
        set([feat.strand for feat in transcriptObj.features if feat.type == "cds"])
    )[0]
    cds_human_readable = which_cds
    which_cds -= 1
    if cds_strand == 1:
        which_cds = which_cds
    else:
        which_cds = len(cds_list) - which_cds - 1

    target_cds = cds_list[which_cds]
    cds_st = target_cds.location.start.position
    cds_en = target_cds.location.end.position
    cds_len = cds_en - cds_st

    # update the cutsite for each gRNA
    for idx, gRNAObj in enumerate(listOfgRNAObj):
        # print(f"{gRNAObj.protospacer} {gRNAObj.pam} {gRNAObj.g_strand} {gRNAObj.g_st} {gRNAObj.g_en}")

        # find gRNA cut site in cds
        gRNACut_in_cds = (
            gRNAObj.g_st - cds_flank_len + 17
        )  # the num of bp before cutsite

        if gRNAObj.g_strand == "-":
            gRNACut_in_cds = (len(cds_w_flank) - 2 * cds_flank_len) - gRNACut_in_cds

        # find gRNA cut site in gene, this is relative to the reference_left_index
        gRNACut_in_gene = gRNACut_in_cds + cds_st
        if (
            cds_strand == -1
        ):  # RNACut_in_cds needs to be from the right if strand == -1 (making it relative to the reference_left_index)
            tmp_gRNACut_in_cds = (len(cds_w_flank) - 2 * cds_flank_len) - gRNACut_in_cds
            gRNACut_in_gene = tmp_gRNACut_in_cds + cds_st

        # find gRNA cut site in chromosome
        gRNACut_in_chr = (
            gRNACut_in_gene + transcriptObj.annotations["reference_left_index"]
        )

        listOfgRNAObj[idx].gRNACut_in_cds = gRNACut_in_cds
        listOfgRNAObj[idx].gRNACut_in_gene = gRNACut_in_gene
        listOfgRNAObj[idx].gRNACut_in_chr = gRNACut_in_chr
        listOfgRNAObj[idx].cds_len = cds_len
        listOfgRNAObj[idx].cds = cds_human_readable
        listOfgRNAObj[idx].Ensemble_ID = transcriptObj.id
        listOfgRNAObj[idx].Ensemble_ref = transcriptObj.annotations["reference_species"]
        listOfgRNAObj[idx].Ensemble_chr = transcriptObj.annotations[
            "reference_chromosome_number"
        ]
        listOfgRNAObj[idx].Ensemble_chr_left_idx = transcriptObj.annotations[
            "reference_left_index"
        ]
        listOfgRNAObj[idx].Ensemble_chr_right_idx = transcriptObj.annotations[
            "reference_right_index"
        ]
        listOfgRNAObj[idx].Ensemble_transcript_strand = transcriptObj.annotations[
            "transcript_strand"
        ]

    return listOfgRNAObj


def get_HDR_flank(transcriptObj, listOfgRNAObj, HDR_flank_len):
    """
    Parameters
    ----------
    transcriptObj
    listOfgRNAObj
    HDR_flank_len: [int] the desired length of the HDR length
    which_cds: the i-th cds used to extract the gRNA (this is used for adjusting the frame)
    adjust_frame: [bool] True/False

    Returns
    -------
    list of gRNA object (updated with the HDR flanks, and whether if the HDR flank is shorter than expected, flags: HDR_Rflank_short and HDR_Lflank_short)
    """
    cds_strand = list(
        set([feat.strand for feat in transcriptObj.features if feat.type == "cds"])
    )[0]
    wholeSeq = str(transcriptObj.seq)
    wholeSeq_rc = str(Seq(wholeSeq).reverse_complement())
    HDR_flank_len = int(HDR_flank_len)
    gene_len = len(wholeSeq)

    for idx, gRNAObj in enumerate(listOfgRNAObj):
        listOfgRNAObj[idx].HDR_Lflank_short = 0
        listOfgRNAObj[idx].HDR_Rflank_short = 0
        gRNACut_in_gene = int(gRNAObj.gRNACut_in_gene)
        HDR_Lflank = ""
        HDR_Rflank = ""
        if cds_strand == 1:
            L_st = gRNACut_in_gene - HDR_flank_len
            L_en = gRNACut_in_gene
            R_st = gRNACut_in_gene
            R_en = gRNACut_in_gene + HDR_flank_len

            if L_st < 0:  # check out of bounds
                L_st = 0
                listOfgRNAObj[idx].HDR_Lflank_short = 1
            if R_st < 0:
                R_st = 0
                listOfgRNAObj[idx].HDR_Rflank_short = 1
            if L_en >= gene_len:
                L_en = gene_len - 1
                listOfgRNAObj[idx].HDR_Lflank_short = 1
            if R_en >= gene_len:
                R_en = gene_len - 1
                listOfgRNAObj[idx].HDR_Rflank_short = 1

            HDR_Lflank = wholeSeq[L_st:L_en]
            HDR_Rflank = wholeSeq[R_st:R_en]

        else:  # strand = -1
            L_st = gRNACut_in_gene
            L_en = gRNACut_in_gene + HDR_flank_len
            R_st = gRNACut_in_gene - HDR_flank_len
            R_en = gRNACut_in_gene

            if L_st < 0:  # check out of bounds
                L_st = 0
                listOfgRNAObj[idx].HDR_Lflank_short = 1
            if R_st < 0:
                R_st = 0
                listOfgRNAObj[idx].HDR_Rflank_short = 1
            if L_en >= gene_len:
                L_en = gene_len - 1
                listOfgRNAObj[idx].HDR_Lflank_short = 1
            if R_en >= gene_len:
                R_en = gene_len - 1
                listOfgRNAObj[idx].HDR_Rflank_short = 1

            HDR_Lflank = str(Seq(wholeSeq_rc[L_st:L_en]).reverse_complement())
            HDR_Rflank = str(
                Seq(wholeSeq_rc[R_st:R_en]).reverse_complement()
            )  # the coord are respective to the revcom of the retrieved seq (weird)

        listOfgRNAObj[idx].HDR_Lflank = HDR_Lflank
        listOfgRNAObj[idx].HDR_Rflank = HDR_Rflank

    return listOfgRNAObj


def calculate_offset_cutsite(
    gRNACut_in_cds, leading_nonTriplet
):  # calculate the offset so the edit will be in frame
    """
    Parameters
    ----------
    leading_nonTriplet

    Returns
    -------
    adjusted L_en
    """
    adj = 0
    if gRNACut_in_cds == 0:
        if leading_nonTriplet == 0:
            adj = 0
        elif leading_nonTriplet == 1:
            adj = 1
        elif leading_nonTriplet == 2:
            adj = -1
    elif gRNACut_in_cds == 1:
        if leading_nonTriplet == 0:
            adj = -1
        elif leading_nonTriplet == 1:
            adj = 0
        elif leading_nonTriplet == 2:
            adj = 1
    elif gRNACut_in_cds == 2:
        if leading_nonTriplet == 0:
            adj = 1
        elif leading_nonTriplet == 1:
            adj = -1
        elif leading_nonTriplet == 2:
            adj = 0
    elif gRNACut_in_cds >= 3:
        trailing_nonTriplet = (gRNACut_in_cds - leading_nonTriplet) % 3
        if trailing_nonTriplet == 0:
            adj = 0
        elif trailing_nonTriplet == 1:
            adj = -1
        elif trailing_nonTriplet == 2:
            adj = 1

    return adj


def nudge_cutsite_inframe(transcriptObj, listOfgRNAObj, HDR_flank_len, which_cds):
    """
    nudge the cutsite so that it is in frame (although it doesn't reflect the true cutsite now)
    Parameters
    ----------
    L_en: Lflank end
    leading_nonTriplet

    Returns
    -------
    adjusted L_en
    """
    cds_strand = list(
        set([feat.strand for feat in transcriptObj.features if feat.type == "cds"])
    )[0]
    wholeSeq = str(transcriptObj.seq)
    wholeSeq_rc = str(Seq(wholeSeq).reverse_complement())
    HDR_flank_len = int(HDR_flank_len)
    gene_len = len(wholeSeq)

    # calculate the number of leading nonTriplet at the 5' of the cds/exon
    frame = get_cds_frame(transcriptObj, which_cds)
    if frame == 1:
        leading_nonTriplet = 0
    if frame == 2:
        leading_nonTriplet = 1
    if frame == 3:
        leading_nonTriplet = 2

    for idx, gRNAObj in enumerate(listOfgRNAObj):
        gRNACut_in_cds = int(gRNAObj.gRNACut_in_cds)
        adj = calculate_offset_cutsite(gRNACut_in_cds, leading_nonTriplet)
        # nudge cutsite (respective to the gene)
        if cds_strand == -1:
            adj = 0 - adj
        listOfgRNAObj[idx].gRNACut_in_gene = listOfgRNAObj[idx].gRNACut_in_gene + adj
        listOfgRNAObj[idx].cds_leading_nonTriplet = leading_nonTriplet
    return listOfgRNAObj


def get_cds_frame(mytranscript, which_cds):
    """
    for the n-th cds, calculate the frame (using coding sequence from 1 to n-1 th cds)
    """
    wholeSeq = str(mytranscript.seq)
    wholeSeq_rc = str(Seq(wholeSeq).reverse_complement())

    cds_list = [feat for feat in mytranscript.features if feat.type == "cds"]
    strand = list(
        set([feat.strand for feat in mytranscript.features if feat.type == "cds"])
    )[0]
    which_cds -= 1
    cds_until = ""

    if strand == 1:
        which_cds = which_cds
        for i in range(0, which_cds):
            st = cds_list[i].location.start.position
            en = cds_list[i].location.end.position
            cds_until = cds_until + wholeSeq[st:en]
    else:
        which_cds = len(cds_list) - which_cds - 1
        for i in range(len(cds_list) - 1, which_cds, -1):
            st = cds_list[i].location.start.position
            en = cds_list[i].location.end.position
            cds_until = cds_until + str(Seq(wholeSeq_rc[st:en]).reverse_complement())

    # calculate frame
    frame = len(cds_until) % 3
    if frame == 0:
        frame = 1
    elif frame == 1:
        frame = 3
    elif frame == 2:
        frame = 2
    return frame


if __name__ == "__main__":
    import doctest

    doctest.testmod()
