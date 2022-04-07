from Bio.Seq import Seq



def get_cds_seq_in_transcript(mytranscript):
    '''
    input: Bio.SeqRecord  (such as that from function fetch_ensembl_transcript)
    return the cds sequence as a string
    '''
    wholeSeq = str(mytranscript.seq)
    wholeSeq_rc = str(Seq(wholeSeq).reverse_complement())
    cds_seqs = []
    negative_strand_flag = False
    for feat in mytranscript.features:
        if feat.type == "cds":
            st = feat.location.start.position
            en = feat.location.end.position
            if (feat.strand==-1): # neg strand
                cds_seqs.append(str(Seq(wholeSeq_rc[st:en]).reverse_complement())) # the coord are respective to the revcom of the retrieved seq (weird)
                negative_strand_flag = True
            else:  # pos strand
                cds_seqs.append(wholeSeq[st:en])
    if negative_strand_flag == True:
        return(''.join(cds_seqs[::-1]))
    else:
        return(''.join(cds_seqs))


# get cds seq and flank
def get_cds_seqNflank(transcriptObj, which_cds, cds_flank_len):
    '''
    returns the n-th cds and its flanking sequences
    the returned sequence will be in the coding strand
    '''
    total_cds_num = transcriptObj.num_cds
    if which_cds > total_cds_num:
        raise ValueError(f"Requested coding exon is {which_cds}, but the transcript only has {total_cds_num} coding exons")

    wholeSeq = str(transcriptObj.seq)
    wholeSeq_rc = str(Seq(wholeSeq).reverse_complement())

    cds_list = [feat for feat in transcriptObj.features if feat.type == 'cds']
    strand = list(set([feat.strand for feat in transcriptObj.features if feat.type == 'cds']))[0]
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
        RFlank_en = min([i for i in [cds_en + cds_flank_len, len(wholeSeq)] if i <= len(wholeSeq)])
        Lflank = wholeSeq[Lflank_st:cds_st]
        Rflank = wholeSeq[cds_en:RFlank_en]
    else:
        cds_seq = str(Seq(wholeSeq_rc[
                          cds_st:cds_en]).reverse_complement())  # the coord are respective to the revcom of the retrieved seq (weird)
        Lflank_st = max([i for i in [cds_st - cds_flank_len, 0] if i >= 0])
        RFlank_en = min([i for i in [cds_en + cds_flank_len, len(wholeSeq_rc)] if i <= len(wholeSeq_rc)])
        Lflank = str(Seq(wholeSeq_rc[Lflank_st:cds_st]).reverse_complement())
        Rflank = str(Seq(wholeSeq_rc[cds_en:RFlank_en]).reverse_complement())
        Lflank, Rflank = Rflank, Lflank  # for -1 strand, the Lflank is the Rflank

    return ({"cds_seq": cds_seq,
             "Lflank": Lflank,
             "Rflank": Rflank})

def get_exon_concat_with_flank(transcriptObj, which_cds, cds_flank_len):
    '''
    get the nth cds segment from the transcript and concat with flanks
    the returned sequence will be in the coding strand
    '''
    cds_seq_N_flank = get_cds_seqNflank(transcriptObj=transcriptObj, which_cds=which_cds, cds_flank_len=cds_flank_len)
    cds_with_flank = "".join([cds_seq_N_flank["Lflank"], cds_seq_N_flank["cds_seq"], cds_seq_N_flank["Rflank"]])
    return cds_with_flank

def get_cutsite_in_gene(transcriptObj, listOfgRNAObj, which_cds, cds_flank_len):
    '''
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
    '''
    cds_w_flank = get_exon_concat_with_flank(transcriptObj=transcriptObj, which_cds=which_cds, cds_flank_len=cds_flank_len)

    #find cds start and end
    wholeSeq = str(transcriptObj.seq)
    wholeSeq_rc = str(Seq(wholeSeq).reverse_complement())
    cds_list = [feat for feat in transcriptObj.features if feat.type == 'cds']
    cds_strand = list(set([feat.strand for feat in transcriptObj.features if feat.type == 'cds']))[0]
    cds_human_readable = which_cds
    which_cds -= 1
    if cds_strand == 1:
        which_cds = which_cds
    else:
        which_cds = len(cds_list) - which_cds - 1

    target_cds = cds_list[which_cds]
    cds_st = target_cds.location.start.position
    cds_en = target_cds.location.end.position
    cds_len = cds_en-cds_st

    #update the cutsite for each gRNA
    for idx, gRNAObj in enumerate(listOfgRNAObj):
        #print(f"{gRNAObj.protospacer} {gRNAObj.pam} {gRNAObj.g_strand} {gRNAObj.g_st} {gRNAObj.g_en}")

        #find gRNA cut site in cds
        gRNACut_in_cds = gRNAObj.g_st - cds_flank_len + 17 # the num of bp before cutsite

        if gRNAObj.g_strand == "-":
            gRNACut_in_cds = (len(cds_w_flank)-2*cds_flank_len) - gRNACut_in_cds

        #find gRNA cut site in gene, this is relative to the reference_left_index
        gRNACut_in_gene = gRNACut_in_cds + cds_st
        if cds_strand ==  -1: # RNACut_in_cds needs to be from the right if strand == -1 (making it relative to the reference_left_index)
            tmp_gRNACut_in_cds = (len(cds_w_flank)-2*cds_flank_len) - gRNACut_in_cds
            gRNACut_in_gene = tmp_gRNACut_in_cds + cds_st

        #find gRNA cut site in chromosome
        gRNACut_in_chr = gRNACut_in_gene + transcriptObj.annotations["reference_left_index"]

        listOfgRNAObj[idx].gRNACut_in_cds = gRNACut_in_cds
        listOfgRNAObj[idx].gRNACut_in_gene = gRNACut_in_gene
        listOfgRNAObj[idx].gRNACut_in_chr = gRNACut_in_chr
        listOfgRNAObj[idx].cds_len = cds_len
        listOfgRNAObj[idx].cds = cds_human_readable
        listOfgRNAObj[idx].Ensemble_ID = transcriptObj.id
        listOfgRNAObj[idx].Ensemble_ref = transcriptObj.annotations["reference_species"]
        listOfgRNAObj[idx].Ensemble_chr = transcriptObj.annotations["reference_chromosome_number"]
        listOfgRNAObj[idx].Ensemble_chr_left_idx = transcriptObj.annotations["reference_left_index"]
        listOfgRNAObj[idx].Ensemble_chr_right_idx = transcriptObj.annotations["reference_right_index"]
        listOfgRNAObj[idx].Ensemble_transcript_strand = transcriptObj.annotations["transcript_strand"]


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
    cds_strand = list(set([feat.strand for feat in transcriptObj.features if feat.type == 'cds']))[0]
    wholeSeq = str(transcriptObj.seq)
    wholeSeq_rc = str(Seq(wholeSeq).reverse_complement())
    HDR_flank_len = int(HDR_flank_len)
    gene_len = len(wholeSeq)

    for idx, gRNAObj in enumerate(listOfgRNAObj):
        listOfgRNAObj[idx].HDR_Lflank_short = 0
        listOfgRNAObj[idx].HDR_Rflank_short = 0
        gRNACut_in_gene = int(gRNAObj.gRNACut_in_gene)
        HDR_Lflank=""
        HDR_Rflank=""
        if cds_strand == 1:
            L_st = gRNACut_in_gene-HDR_flank_len
            L_en = gRNACut_in_gene
            R_st = gRNACut_in_gene
            R_en = gRNACut_in_gene + HDR_flank_len

            if L_st < 0: #check out of bounds
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

        else: #strand = -1
            L_st= gRNACut_in_gene
            L_en= gRNACut_in_gene + HDR_flank_len
            R_st= gRNACut_in_gene-HDR_flank_len
            R_en= gRNACut_in_gene

            if L_st < 0: #check out of bounds
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
            HDR_Rflank = str(Seq(wholeSeq_rc[R_st:R_en]).reverse_complement())  # the coord are respective to the revcom of the retrieved seq (weird)

        listOfgRNAObj[idx].HDR_Lflank = HDR_Lflank
        listOfgRNAObj[idx].HDR_Rflank = HDR_Rflank

    return listOfgRNAObj

def calculate_offset_cutsite(gRNACut_in_cds , leading_nonTriplet): # calculate the offset so the edit will be in frame
    """
    Parameters
    ----------
    leading_nonTriplet

    Returns
    -------
    adjusted L_en
    """
    adj=0
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
            adj = - 1
        elif leading_nonTriplet == 2:
            adj = 0
    elif gRNACut_in_cds >=3:
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
    cds_strand = list(set([feat.strand for feat in transcriptObj.features if feat.type == 'cds']))[0]
    wholeSeq = str(transcriptObj.seq)
    wholeSeq_rc = str(Seq(wholeSeq).reverse_complement())
    HDR_flank_len = int(HDR_flank_len)
    gene_len = len(wholeSeq)

    # calculate the number of leading nonTriplet at the 5' of the cds/exon
    frame = get_cds_frame(transcriptObj, which_cds)
    if frame == 1: leading_nonTriplet = 0
    if frame == 2: leading_nonTriplet = 1
    if frame == 3: leading_nonTriplet = 2

    for idx, gRNAObj in enumerate(listOfgRNAObj):
        gRNACut_in_cds = int(gRNAObj.gRNACut_in_cds)
        adj = calculate_offset_cutsite(gRNACut_in_cds, leading_nonTriplet)
        #nudge cutsite (respective to the gene)
        if cds_strand == -1:
            adj = 0 - adj
        listOfgRNAObj[idx].gRNACut_in_gene = listOfgRNAObj[idx].gRNACut_in_gene + adj
        listOfgRNAObj[idx].cds_leading_nonTriplet = leading_nonTriplet
    return listOfgRNAObj

# calcualte frame
def get_cds_frame(mytranscript, which_cds):
    """
    for the n-th cds, calculate the frame (using coding sequence from 1 to n-1 th cds)
    """
    wholeSeq = str(mytranscript.seq)
    wholeSeq_rc = str(Seq(wholeSeq).reverse_complement())

    cds_list = [feat for feat in mytranscript.features if feat.type == 'cds']
    strand = list(set([feat.strand for feat in mytranscript.features if feat.type == 'cds']))[0]
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
    return (frame)



#################
#custom logging #
#################
import logging

BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE = range(8)
#The background is set with 40 plus the number of the color, and the foreground with 30

#These are the sequences need to get colored ouput
RESET_SEQ = "\033[0m"
COLOR_SEQ = "\033[1;%dm"
BOLD_SEQ = "\033[1m"

def formatter_message(message, use_color = True):
    if use_color:
        message = message.replace("$RESET", RESET_SEQ).replace("$BOLD", BOLD_SEQ)
    else:
        message = message.replace("$RESET", "").replace("$BOLD", "")
    return message

COLORS = {
    'WARNING': YELLOW,
    'INFO': WHITE,
    'DEBUG': BLUE,
    'CRITICAL': YELLOW,
    'ERROR': RED
}

class ColoredFormatter(logging.Formatter):
    def __init__(self, msg, use_color = True):
        logging.Formatter.__init__(self, msg)
        self.use_color = use_color

    def format(self, record):
        levelname = record.levelname
        if self.use_color and levelname in COLORS:
            levelname_color = COLOR_SEQ % (30 + COLORS[levelname]) + levelname + RESET_SEQ
            record.levelname = levelname_color
        return logging.Formatter.format(self, record)

# Custom logger class with multiple destinations
class ColoredLogger(logging.Logger):
    FORMAT = "[$BOLD%(name)-1s$RESET][%(levelname)-1s]  %(message)s " #($BOLD%(filename)s$RESET:%(lineno)d)
    COLOR_FORMAT = formatter_message(FORMAT, True)
    def __init__(self, name):
        logging.Logger.__init__(self, name, logging.DEBUG)

        color_formatter = ColoredFormatter(self.COLOR_FORMAT)

        console = logging.StreamHandler()
        console.setFormatter(color_formatter)

        self.addHandler(console)
        return