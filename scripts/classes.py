
from scoring_routine import gcContent
from scoring_routine import manu_score
############
# gRNA class#
############

class gRNA:
    """
    Stores one gRNA

    Arguments:
        name:               str    #*required*
        protospacer:        str    #*required*
        pam:                str   #PAM sequence,          default = "NGG"
        pos_strand:         bool   #is gRNA on pos strand, default = True
        prsp_len:           int    #length of protospacer, default = 20

    optional:
        parent gene name:               str
        parent gene seq:                str
        g_start:                        int
        g_end:                          int
        g_strand:                       str  [+ or -]
        parent gene is transcript:      bool
        parent gene coding frame:       int  [1,2 or 3]

    Derived properties:
        CFD (Cutting Frequency Determination score)     float
        GC  (GC content)                                float 
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
            pam3p_flank: str = '',
            protospacer5p_flank: str = '') -> None:
        
        self.prsp_len = prsp_len
        self.name = name
        self.protospacer = protospacer
        self.pam = pam
        self.pos_strand = pos_strand
        self.g_st = g_st
        self.g_en = g_en
        self.g_strand = g_strand
        self.valid = self.__validate()
        self.offTargets = []
        self.pam3p_flank = pam3p_flank
        self.protospacer5p_flank = protospacer5p_flank

        # add calls to score calculations here, and update the scores accordingly
        self.gc_score = gcContent(protospacer)
        #self.manu_score = manu_score()

    def __validate(self) -> bool:
        len_flag = len(self.protospacer) >= 1
        return all([len_flag])

    def __new_off_target(self,
                       off_target: object) -> None:
        self.offTargets.append(off_target)


###################
# off-target class#
###################

class off_target:
    """
    Stores one off-target

    Arguments:
        parent_gRNA:        object #gRNA object            *required*
        contig_name: str    #contig/chromosome name *required*
        start:       int,   #match start            *required*
        end:         int    #match end              *required*
        mm_pos:      list   #mis-match positions    *required*

    """

    def __init__(
            self,
            parent_gRNA: object,
            contig_name: str,
            start: str,
            end: bool,
            mm_pos: int) -> None:
        self.parent_gRNA = parent_gRNA
        self.contig_name = contig_name
        self.start = start
        self.end = end
        self.MM_pos = mm_pos


#############
# gene class#
#############

class gene:
    """
    Stores one gene

    Arguments:
        gid:         str    #*required*
        gid_type:    str    #*required*
        seq:        str
    """

    def __init__(
            self,
            gid: str = '',
            gid_type: str = '',
            seq: str = '') -> None:
        self.gid = gid
        self.gid_type = gid_type
        self.seq = seq

    def _validate(self) -> bool:
        flag = (self.gid != '') and (self.gid_type != '')
        return flag
