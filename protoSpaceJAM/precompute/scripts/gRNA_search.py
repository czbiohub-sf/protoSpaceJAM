import re
from Bio.Seq import Seq
from utils import *
import gc

class gRNA:
    """
    Stores one gRNA
    """

    def __init__(
            self,
            name: str = '',
            protospacer: str = '',
            pam: str = 'NGG',
            pos_strand: bool = True,
            g_st: int = "",
            g_en: int = "",
            g_strand: str = "",
            prsp_len: int = '',
            protospacer5p_flank: str = '',
            pam3p_flank: str = '') -> None:
        
        self.prsp_len = prsp_len
        self.name = name
        self.protospacer = protospacer
        self.pam = pam
        self.pos_strand = pos_strand
        self.g_st = g_st
        self.g_en = g_en
        self.g_strand = g_strand
        self.protospacer5p_flank = protospacer5p_flank
        self.pam3p_flank = pam3p_flank
        self.offTargets = []

def search_gRNA(protosp_len, PAM, search_in, flanksize=100):
    """
    Features:
        search gRNA in both positive and negative strand of the input
        supports IUPAC nucleotide code

    Return:
        a list of gRNA objects
    """
    protosp_len = protosp_len
    gPAM = PAM
    gPAM = re.sub(
        "N", "[ACGT]", gPAM, flags=re.IGNORECASE
    )  # IUPAC nucleotide code http://www.bioinformatics.org/sms/iupac.html
    gPAM = re.sub("R", "[AG]", gPAM, flags=re.IGNORECASE)  # IUPAC nucleotide code
    gPAM = re.sub("Y", "[CT]", gPAM, flags=re.IGNORECASE)  # IUPAC nucleotide code
    gPAM = re.sub("S", "[GC]", gPAM, flags=re.IGNORECASE)  # IUPAC nucleotide code
    gPAM = re.sub("W", "[AT]", gPAM, flags=re.IGNORECASE)  # IUPAC nucleotide code
    gPAM = re.sub("K", "[GT]", gPAM, flags=re.IGNORECASE)  # IUPAC nucleotide code
    gPAM = re.sub("M", "[AC]", gPAM, flags=re.IGNORECASE)  # IUPAC nucleotide code
    gPAM = re.sub("B", "[CGT]", gPAM, flags=re.IGNORECASE)  # IUPAC nucleotide code
    gPAM = re.sub("D", "[AGT]", gPAM, flags=re.IGNORECASE)  # IUPAC nucleotide code
    gPAM = re.sub("H", "[ACT]", gPAM, flags=re.IGNORECASE)  # IUPAC nucleotide code
    gPAM = re.sub("V", "[ACG]", gPAM, flags=re.IGNORECASE)  # IUPAC nucleotide code
    gRNA_regex = re.compile(
        r"(?=(.{" + re.escape(str(protosp_len)) + r"})" + r"(" + gPAM + r"))",
        re.IGNORECASE,
    )  # use lookahead assertion (?=...) to find overlapping gRNAs
    # print(gRNA_regex)
    gRNA_list = []

    # search in positive strand
    for m in gRNA_regex.finditer(search_in):
        if m:
            if len(m.groups()) >= 2:
                protosp = m.group(1)
                gPAM = m.group(2)
                st, en = m.span()
                end = st + len(protosp) + len(gPAM)
                st = st + 1  # change from 0-index to 1-indexed

                # get flanks, st and end are 1-indexed
                lf = get_leftflank(
                    seq=search_in, pos=st, flanksize=flanksize
                )  # pos always is in respect to the 5'
                rf = get_rightflank(
                    seq=search_in, pos=end, flanksize=flanksize
                )  # pos always is in respect to the 5'

                gRNA_list.append(
                    gRNA(
                        name=f"gRNA_+_{st}_{end}",
                        protospacer=protosp,
                        pam=gPAM,
                        g_strand="+",
                        g_st=st,
                        g_en=end,
                        prsp_len=len(protosp),
                        protospacer5p_flank=lf,
                        pam3p_flank=rf,
                    )
                )

    # search in negative strand
    search_in_revcom = str(Seq(search_in).reverse_complement())  # revcom
    seq_len = len(search_in_revcom)
    for m in gRNA_regex.finditer(search_in_revcom):
        if m:
            if len(m.groups()) >= 2:
                protosp = m.group(1)
                gPAM = m.group(2)
                st, en = m.span()
                end = st + len(protosp) + len(gPAM)

                # get flanks, need to convert to 1-indexed
                lf = get_leftflank(
                    seq=search_in_revcom, pos=st + 1, flanksize=flanksize
                )  # pos always is in respect to the 5'. We need flanks on the - strand, thus before changing to indexing from the 5' of the + strand
                rf = get_rightflank(seq=search_in_revcom, pos=end, flanksize=flanksize)

                st = (
                    seq_len - st
                )  # change to indexing from the 5' of the + strand, and also change from 0-index to 1-indexed
                end = (
                    seq_len - end + 1
                )  # change to indexing from the 5' of the + strand,

                gRNA_list.append(
                    gRNA(
                        name=f"gRNA_-_{st}_{end}",
                        protospacer=protosp,
                        pam=gPAM,
                        g_strand="-",
                        g_st=st,
                        g_en=end,
                        prsp_len=len(protosp),
                        protospacer5p_flank=lf,
                        pam3p_flank=rf,
                    )
                )
    del search_in_revcom
    del search_in
    gc.collect()
    return gRNA_list


def get_rightflank(seq, pos, flanksize):
    """
    Parameters
    ----------
    seq
    pos: position (1-indexed), will start extracting seq immediately after this position
    flanksize

    Returns
    -------
    the flank sequence to the *right* of pos
    """
    if flanksize < 0 or pos < 0:
        raise ValueError(
            f"flanksize must >=0 and pos must >=0; flanksize is set to {flanksize}, pos is set to {pos}"
        )
    if pos > len(seq):
        raise ValueError(f"pos must <= length of seq, pos is set to {pos}")

    remaining_bases = len(seq) - pos
    if flanksize <= remaining_bases:  # enough length
        return seq[pos : (pos + flanksize)]
    else:  # not enough length
        return seq[pos:] + "-" * (flanksize - remaining_bases)


def get_leftflank(seq, pos, flanksize):
    """
    Parameters
    ----------
    seq
    pos: position (1-indexed), will start extracting seq immediately after this position
    flanksize

    Returns
    -------
    the flank sequence to the *left* of pos
    """
    if flanksize < 0 or pos < 0:
        raise ValueError(
            f"flanksize must >=0 and pos must >=0; flanksize is set to {flanksize}, pos is set to {pos}"
        )
    if pos > len(seq):
        raise ValueError(f"pos must <= length of seq, pos is set to {pos}")

    remaining_bases = pos - 1
    if flanksize <= remaining_bases:  # enough length
        return seq[(pos - flanksize - 1) : pos - 1]
    elif pos == 0:
        print(f"[WARNING]: pos was set to 0, while it should be 1-indexed ")
        return "-" * flanksize
    else:  # not enough length
        return "-" * (flanksize - remaining_bases) + seq[: pos - 1]


def search_gRNA_in_cdsList(
    cdsList, my_transcript, cds_flank_len, PAM, protosp_len, adjust_frame, HDR_flank_len
):
    """

    Parameters
    ----------
    cdsList
    my_transcript
    cds_flank_len
    PAM
    protosp_len
    adjust_frame
    HDR_flank_len

    Returns
    -------
    a list of gRNA objects
    """
    res_gRNA_list = []
    # process cds sequentially
    for current_cds in cdsList:

        which_cds = current_cds + 1  # which_cds is 1 indexed

        # print(f"working on cds #{which_cds}")

        # get exon sequence from transcript object
        cds_w_flank = get_exon_concat_with_flank(
            transcriptObj=my_transcript,
            which_cds=which_cds,
            cds_flank_len=cds_flank_len,
        )

        # search gRNA and update scores
        gRNA_list = search_gRNA(protosp_len=protosp_len, PAM=PAM, search_in=cds_w_flank)

        # update gRNA objects with cut site
        gRNA_list = get_cutsite_in_gene(
            transcriptObj=my_transcript,
            listOfgRNAObj=gRNA_list,
            which_cds=which_cds,
            cds_flank_len=cds_flank_len,
        )

        # adjust frame
        if adjust_frame == True:
            gRNA_list = nudge_cutsite_inframe(
                transcriptObj=my_transcript,
                listOfgRNAObj=gRNA_list,
                HDR_flank_len=HDR_flank_len,
                which_cds=which_cds,
            )

        # update gRNA objects with flanks
        gRNA_list = get_HDR_flank(
            transcriptObj=my_transcript,
            listOfgRNAObj=gRNA_list,
            HDR_flank_len=HDR_flank_len,
        )

        # for gRNA in gRNA_list:
        #    print(f"{gRNA.protospacer} {gRNA.pam} {gRNA.g_strand} {gRNA.g_st} {gRNA.g_en}")

        res_gRNA_list.extend(gRNA_list)

    return res_gRNA_list
