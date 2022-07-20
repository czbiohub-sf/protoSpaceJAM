"""
Taken from https://github.com/czbiohub/crispycrunch/blob/master/utils/hdr.py
Transformations of genome sequences for HDR.
"""
import functools
import itertools
import logging
from typing import List
import sys
from typing import Iterator
import copy
from Bio.Seq import Seq

try:
    from . import cfdscore, mitscore
except ImportError:
    import cfdscore  # type: ignore
    import mitscore  # type: ignore

logger = logging.getLogger(__name__)

class seq_w_phase:
    '''
    seqs with phases and start + end info
    '''
    def __init__(self,
                 seq:str,
                 phases:str,
                 start:int,
                 end:int)->None:
        self.seq = seq
        self.phases = phases
        self.start = start
        self.end = end

#by DP
class HDR_flank:

    def __init__(
            self,
            left_flk_seq: str,
            right_flk_seq: str,
            left_flk_coord_lst: List[int],
            right_flk_coord_lst: List[int],
            left_flk_phases: List[int],
            right_flk_phases: List[int],
            type:str,
            ENST_ID: str,
            ENST_strand : int,
            gStart : int,
            gStrand : int,
            InsPos:int,
            CutPos:int,
            Cut2Ins_dist:int,
            tag:str) -> None:

        self.left_flk_seq = left_flk_seq
        self.right_flk_seq = right_flk_seq
        self.ENST_ID = ENST_ID
        self.ENST_strand = ENST_strand
        self.gStrand = gStrand
        self.gStart = gStart
        self.InsPos = InsPos # InsPos is the first letter of stop codon "T"AA or the last letter of the start codon AT"G"
        self.CutPos = CutPos
        self.Cut2Ins_dist = Cut2Ins_dist
        self.tag = tag

        assert len(left_flk_coord_lst) > 1
        self.left_flk_coord_lst = left_flk_coord_lst

        assert len(right_flk_coord_lst) > 1
        self.right_flk_coord_lst = right_flk_coord_lst

        assert type in ('start', 'stop')
        self.type = type

        assert len(left_flk_phases) > 1
        self.left_flk_phases = self.join_int_list(left_flk_phases)

        assert len(right_flk_phases) > 1
        self.right_flk_phases = self.join_int_list(right_flk_phases)

        # adjust insPos, for stop-tagging ,the insertion site is now before the first base of the stop codon
        if type == "stop":
            if self.ENST_strand == 1 or self.ENST_strand == "1" or self.ENST_strand == "+":
                self.InsPos = int(self.InsPos) - 1
            else:
                self.InsPos = int(self.InsPos) + 1
        # adjust cutPos for -1 strands genes (the gRNA cut site was based on +1 strand), we need to view the gene in coding sequence
        if self.ENST_strand == -1 or self.ENST_strand == "-1" or self.ENST_strand == "-":
            self.CutPos = self.CutPos + 1


        #TODO: extend the flank to include the gRNA (de-prioritized b/c with 100bp arm and 50bp max-cut-to-insert-distance, gRNA will never be outside the HDR arm)
        ATG_at_end_of_exon = False
        self.cutPos, self.gPAM_end = self.get_gRNA_pos(self.gStart, self.gStrand)
        #check if the whole 23nt gRNA is in the HDR flank
        self.entire_gRNA_in_HDR_arms, self.gStart_in_HDR_arms, self.gPAM_end_in_HDR_arms = self.check_gRNA_in_HDR_arms(gStart=self.gStart, gPAM_end=self.gPAM_end, left_flk_coord_lst=self.left_flk_coord_lst, right_flk_coord_lst=self.right_flk_coord_lst)

        #TODO: check if ATG is at the end of exon

        #make gRNA lowercase
        self.gRNA_lc_Larm, self.gRNA_lc_Rarm = self.make_gRNA_lowercase()

        #get seq and phase between insertion and cut site
        self.ins2cut = self.get_ins2cut_seq()

        #check if the ins2cut_seq is the same length as the Cut2Ins_dist calculated elsewhere
        if len(self.ins2cut.seq)>0 and (len(self.ins2cut.seq)-4 != abs(self.Cut2Ins_dist)):
            sys.exit(f"ins2cut_seq:{self.ins2cut.seq} is not the same length as reported: Cut2Ins_dist={self.Cut2Ins_dist}")

        #trim insert-to-cut into frame
        self.ins2cut_Ltrimed = self.trim_left_into_frame(self.ins2cut)
        self.ins2cut_LRtrimed = self.trim_right_into_frame(self.ins2cut_Ltrimed)

        # mutate insert-to-cut sequence
        mutated_subseq = ""
        for mutated_subseq in mutate_silently(self.ins2cut_LRtrimed.seq.upper()): # get the last item in the generator
            pass
        #put_silent_mutation_subseq_back into HDR flank
        left_flk_seq_CodonMut, right_flk_seq_CodonMut = self.put_silent_mutation_subseq_back(mutated_subseq = mutated_subseq, start = self.ins2cut_LRtrimed.start, end = self.ins2cut_LRtrimed.end)

        #selected strand
        ssODN = self.select_ssODN_strand(left_flk_seq_CodonMut + right_flk_seq_CodonMut)

        #check if gRNA is affected by insertion
        self.gRNA_seq, Null = self.get_post_integration_gRNA(self.left_flk_seq,self.right_flk_seq) #get the original gRNA
        Null, self.post_mut_ins_gRNA = self.get_post_integration_gRNA(left_flk_seq_CodonMut,right_flk_seq_CodonMut) #get the post mutation and insertion gRNA


        #print for debug purposes
        print(f"{self.ENST_ID}\tstrand:{self.ENST_strand}\ttype:{type}-tagging\tInsPos:{self.InsPos}\tgRNA:{self.gStart}-{self.gPAM_end}\tstrand:{self.gStrand}\tCutPos:{self.CutPos}\tCut2Ins-dist:{self.Cut2Ins_dist}")
        print(f"1. left | right arms: {self.gRNA_lc_Larm}|{self.gRNA_lc_Rarm}\n"
              f"2. Phases           : {self.left_flk_phases}|{self.right_flk_phases}\n"
              f"3. Coordinates      :\t{self.left_flk_coord_lst[0]}-{self.left_flk_coord_lst[1]} | {self.right_flk_coord_lst[0]}-{self.right_flk_coord_lst[1]}\n"
              f"---------------------------------------------\n"
              f"4. cut2insert (with 2bp padding in lowercase)\n"
              f"5. seq         :{self.ins2cut.seq}\n"
              f"6. Phases      :{self.ins2cut.phases}\n"
              f"7. Coordinates :{self.ins2cut.start}-{self.ins2cut.end}\n"
              f"---------------------------------------------\n"
              f"8. cut2insert (trimmed into frame):\n"
              f"9. seq         :{self.ins2cut_LRtrimed.seq}\n"
              f"10.Phases      :{self.ins2cut_LRtrimed.phases}\n"
              f"11.Coordinates :{self.ins2cut_LRtrimed.start}-{self.ins2cut_LRtrimed.end}\n"
              f"12.mutated seq :{mutated_subseq}\n"
              f"---------------------------------------------\n"
              f"13.updated arms :{left_flk_seq_CodonMut}|{right_flk_seq_CodonMut}\n"
              f"14.original arms:{left_flk_seq}|{right_flk_seq}\n"
              f"15.ssODN strand :{ssODN}\n"
              f"16.                                  gRNA:{self.gRNA_seq}\n"
              f"17.gRNA post codon mutation and insertion:{self.post_mut_ins_gRNA}\n")

    #TODO: check CFD score, if too high get gRNA codon and mutate
    #TODO:


    #END OF INIT
    def get_post_integration_gRNA(self, leftArm, rightArm):
        """
        input: leftArm, rightArm (based on which the gRNA and insert-disrupted gRNA will be extracted
        return: gRNA, chimeric_gRNA
        """
        Lstart = self.left_flk_coord_lst[0]
        Rend = self.right_flk_coord_lst[1]

        if (Lstart<Rend and self.ENST_strand==-1) or (Lstart>Rend and self.ENST_strand==1):
            sys.exit("start<end while strand==-1 or start>end while strand==1")

        newLarm = leftArm
        newRarm = rightArm
        whole_arm = newLarm + newRarm

        # get the left, right coordinate of the gRNA
        if Lstart<Rend:
            newgRNAstart = min([self.gStart - Lstart , self.gPAM_end - Lstart])
            newgRNAend = max([self.gStart - Lstart , self.gPAM_end - Lstart])
            gRNAleft=newgRNAstart
            gRNAright=newgRNAend
        else:
            newgRNAstart = min(self.gPAM_end - Rend, self.gStart - Rend) #the ENSTID is on the -1 strand, using Rend as reference coordinate
            newgRNAend = max(self.gPAM_end - Rend, self.gStart - Rend)
            gRNAleft = len(whole_arm) - newgRNAend - 1
            gRNAright = len(whole_arm) - newgRNAstart - 1
        #get the gRNA (as in the coding strand)
        gRNA = whole_arm[gRNAleft:gRNAright+1]

        #get truncate gRNA:
        chimeric_gRNA = gRNA
        if (gRNAleft < int(len(whole_arm)/2) < gRNAright) and (int(self.ENST_strand) * int(self.gStrand) > 0): #gRNA is truncated, and gRNA is on the coding strand
            chimeric_gRNA = whole_arm[int(len(whole_arm)/2):gRNAright+1]
            trunc_len = int(len(whole_arm)/2) - gRNAleft
            chimeric_gRNA = self.tag[-trunc_len:] + chimeric_gRNA #make the chimeric gRNA
        elif (gRNAleft < int(len(whole_arm)/2) < gRNAright) and (int(self.ENST_strand) * int(self.gStrand) < 0): #gRNA is truncated, and gRNA is NOT on the coding strand
            chimeric_gRNA = whole_arm[gRNAleft:int(len(whole_arm)/2)+1]
            trunc_len = gRNAright - int(len(whole_arm)/2)
            chimeric_gRNA = chimeric_gRNA + self.tag[:trunc_len] #make the chimeric gRNA

        #reverse complement gRNA, if gRNA is on a different strand compared to the
        if (int(self.ENST_strand) * int(self.gStrand) < 0):
            #revcom gRNA and trunc_gRNA
            gRNA = str(Seq(gRNA).reverse_complement())
            chimeric_gRNA =  str(Seq(chimeric_gRNA).reverse_complement())

        return gRNA, chimeric_gRNA

    def select_ssODN_strand(self, seq):
        """
        select strand
        input: seq
        output: seq in the preferred strand
        """
        if self.InsPos <= self.CutPos:
            return seq
        else:
            return str(Seq(seq).reverse_complement())

    def put_silent_mutation_subseq_back(self, mutated_subseq, start, end):
        if start == end:
            return self.left_flk_seq, self.right_flk_seq
        else:
            whole_arm = self.left_flk_seq + self.right_flk_seq
            whole_arm = whole_arm[0:start] + mutated_subseq + whole_arm[end+1:]
            arm_len = int(len(whole_arm)/2)
            return whole_arm[:arm_len], whole_arm[arm_len:]

    def trim_left_into_frame(self, in_obj):
        obj = copy.copy(in_obj)
        if obj.seq == "":
            return obj
        for i in copy.copy(obj.phases):
            if i != "1":
                obj.seq = obj.seq[1:]
                obj.phases = obj.phases[1:]
                obj.start = obj.start + 1
            else:
                return(obj)
        if obj.seq == "": #trimmed all of the sequence
            obj.start = obj.end
        return(obj)

    def trim_right_into_frame(self, in_obj):
        obj = copy.copy(in_obj)
        if obj.seq == "":
            return obj
        for i in reversed(copy.copy(obj.phases)):
            if i != "3":
                obj.seq = obj.seq[:-1]
                obj.phases = obj.phases[:-1]
                obj.end = obj.end - 1
            else:
                return(obj)
        if obj.seq == "": #trimmed all of the sequence
            obj.end = obj.start
        return (obj)

    def to_0_index(self, start, end, strand, seq=""):
        if (start<end and strand==-1) or (start>end and strand==1):
            sys.exit("start<end while strand==-1 or start>end while strand==1")
        if start<end:
            newStart = start - start
            newEnd = end - start
            return(seq,newStart,newEnd)
        else:
            seq = seq[::-1]
            newStart=end - end
            newEnd=start - end
            return(seq,newStart,newEnd)

    def get_ins2cut_seq(self):
        '''
        get the sequence between insertion and cut sites
        return [seq, phases, start, end]
        Note, the sequence is padded with 2bp on each side (to avoid truncating codons)
        '''
        if self.InsPos == self.CutPos: #cut and insertion are at the same site
            ins2cut = seq_w_phase(seq = "", phases = "", start = self.InsPos, end = self.CutPos)
            return ins2cut

        #convert into 0-index
        Lstart = self.left_flk_coord_lst[0]
        Rend = self.right_flk_coord_lst[1]
        if (Lstart<Rend and self.ENST_strand==-1) or (Lstart>Rend and self.ENST_strand==1):
            sys.exit("start<end while strand==-1 or start>end while strand==1")

        if Lstart<Rend:
            newLstart = Lstart - Lstart
            newRend = Rend -Lstart
            ins2cutStart = min([self.InsPos - Lstart + 1, self.CutPos  - Lstart + 1]) # exclude the base after which the cut/insertion is made
            ins2cutEnd = max([self.InsPos - Lstart , self.CutPos - Lstart])
            newLarm = self.left_flk_seq
            newRarm = self.right_flk_seq
            newLarm_phases = self.left_flk_phases
            newRarm_phases = self.right_flk_phases
            whole_arm = newLarm + newRarm
            whole_arm_phases = newLarm_phases + newRarm_phases
        else:
            newLstart = Rend - Rend
            newRend = Lstart - Rend
            ins2cutStart = min(self.InsPos - Rend, self.CutPos - Rend)
            ins2cutEnd = max(self.InsPos - Rend -1, self.CutPos - Rend -1) # exclude the base after which the cut/insertion is made
            newLarm = self.left_flk_seq[::-1]
            newRarm = self.right_flk_seq[::-1]
            newLarm_phases = self.left_flk_phases[::-1]
            newRarm_phases = self.right_flk_phases[::-1]
            whole_arm = newRarm + newLarm
            whole_arm_phases = newRarm_phases + newLarm_phases

        #get the ins2cut_seq
        ins2cut_seq = whole_arm[ins2cutStart-2:ins2cutStart].lower() + whole_arm[ins2cutStart:ins2cutEnd+1] + whole_arm[ins2cutEnd+1:ins2cutEnd+3].lower()
        ins2cut_phases = whole_arm_phases[ins2cutStart-2:ins2cutStart] + whole_arm_phases[ins2cutStart:ins2cutEnd+1] + whole_arm_phases[ins2cutEnd+1:ins2cutEnd+3]

        if Lstart>Rend:
            ins2cut_seq = ins2cut_seq[::-1] #reverse
            ins2cut_phases = ins2cut_phases[::-1] #reverse

        ins2cut = seq_w_phase(seq = ins2cut_seq, phases = self.join_int_list(ins2cut_phases), start = ins2cutStart-2, end = ins2cutEnd+2)
        return ins2cut

    def make_gRNA_lowercase(self):
        #convert into 0-index
        Lstart = self.left_flk_coord_lst[0]
        Rend = self.right_flk_coord_lst[1]
        if (Lstart<Rend and self.ENST_strand==-1) or (Lstart>Rend and self.ENST_strand==1):
            sys.exit("start<end while strand==-1 or start>end while strand==1")

        if Lstart<Rend:
            newLstart = Lstart - Lstart
            newRend = Rend -Lstart
            newgRNAstart = min([self.gStart - Lstart , self.gPAM_end - Lstart])
            newgRNAend = max([self.gStart - Lstart , self.gPAM_end - Lstart])
            newLarm = self.left_flk_seq
            newRarm = self.right_flk_seq
            whole_arm = newLarm + newRarm
        else:
            newLstart = Rend - Rend
            newRend = Lstart - Rend
            newgRNAstart = min(self.gPAM_end - Rend, self.gStart - Rend)
            newgRNAend = max(self.gPAM_end - Rend, self.gStart - Rend)
            newLarm = self.left_flk_seq[::-1]
            newRarm = self.right_flk_seq[::-1]
            whole_arm = newRarm + newLarm
        #get the gRNA
        gRNA_seq = whole_arm[newgRNAstart:newgRNAend+1]
        #make gRNA lowercase
        newWhole_arm = whole_arm[0:newgRNAstart].upper() + gRNA_seq.lower() + whole_arm[newgRNAend+1:].upper()
        if Lstart>Rend:
            newWhole_arm = newWhole_arm[::-1] #reverse
        newLarm = newWhole_arm[0:int((newRend+1)/2)]
        newRarm = newWhole_arm[int((newRend+1)/2):]
        #print for debug
        #print(f"left | right arms: {newLarm}|{newRarm}\n")
        return([newLarm, newRarm])

    #TODO: check the gRNA CFD scores

    def join_int_list(self, mylist):
        return ''.join([str(abs(int(i))) for i in mylist])

    def check_gRNA_in_HDR_arms(self, gStart, gPAM_end, left_flk_coord_lst, right_flk_coord_lst):
        """
        return a list of three bools: [entire_gRNA_in_HDR_arms, gStart_in_HDR_arms, gPAM_end_in_HDR_arms]
        """
        entire_gRNA_in_HDR_arms = False
        max_arm_pos = max(left_flk_coord_lst + right_flk_coord_lst)
        min_arm_pos = min(left_flk_coord_lst + right_flk_coord_lst)
        gStart_in_HDR_arms = min_arm_pos <= gStart <= max_arm_pos
        gPAM_end_in_HDR_arms = min_arm_pos <= gPAM_end <= max_arm_pos

        if gStart_in_HDR_arms and gPAM_end_in_HDR_arms:
            entire_gRNA_in_HDR_arms = True
        #elif not gPAM_end_in_HDR_arms:
        return([entire_gRNA_in_HDR_arms, gStart_in_HDR_arms, gPAM_end_in_HDR_arms])

    def get_gRNA_pos(self, start, strand):
        """
        start:gRNA start
        strand:gRNA strand
        return [cutPos,PAM_endPos]
        """
        if strand == "+" or strand == "1" or strand == 1:
            cutPos = start + 16
            gPAM_end = start + 22
        else:
            cutPos = start - 16
            gPAM_end = start - 22
        return([cutPos,gPAM_end])


#taken from CRISPYcrunch
class MutatedSeq(str):

    # For mypy
    max_score: float
    max_seq: str

    def __new__(cls, val, **attrs):
        inst = super().__new__(cls, val)  # type: ignore
        for a, v in attrs.items():
            setattr(inst, a, v)
        return inst

    # For mypy
    def __init__(self, val, **attrs) -> None:
        pass


class HDR:
    """
    Encapsulates all the HDR transformations of sequences described in
    https://czi.quip.com/YbAhAbOV4aXi/ . Get a mutated HDR inserted that varies
    depending on start or stop codon, the cut-to-insert distance, the
    strandedness of the guide, and the amount of mutation desired.
    The target sequence should be in the direction of the gene. Reading from
    left to right, it should have either a ATG or one of TAG, TGA, or TAA.
    The target sequence must be codon aligned so the target codon can be found!
    The intron/exon junctions should be denoted by lowercase in the target seq.
    target_mutation_score is the minimum MIT score needed to stop silent mutation.
    guide_strand_same refers to strand of target_seq.
    """

    # TODO (gdingle): review whether we still need all these configs
    # TODO (gdingle): write tests that don't rely on defaults which change
    guide_seq_aligned_length = 27
    # Instead of MIT score
    use_cfd_score = True
    # Try mutate all codons, not just linearly.
    # TODO (gdingle): change to default True
    mutate_all_permutations = True

    # Default based on analysis of
    # https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1012-2
    target_mutation_score = 0.1

    def __init__(
            self,
            target_seq: str,
            hdr_seq: str = '',
            hdr_tag: str = 'start_codon',
            hdr_dist: int = 0,
            guide_strand_same: bool = None,
            codon_at: int = -1) -> None:

        _validate_seq(target_seq)
        self.target_seq = target_seq
        _validate_seq(hdr_seq)
        self.hdr_seq = hdr_seq

        assert hdr_tag in ('start_codon', 'stop_codon')
        self.hdr_tag = hdr_tag

        assert abs(hdr_dist) < len(target_seq)
        self.hdr_dist = hdr_dist

        # TODO (gdingle): refactor _target_codon_at
        assert codon_at < len(target_seq)
        self._codon_at = codon_at
        if hdr_tag == 'start_codon':
            self.boundary_codons = set(['ATG'])
            # just after start codon
            self.insert_at = self._target_codon_at() + 3
        else:
            self.boundary_codons = set(['TAG', 'TGA', 'TAA'])
            # just before stop codon
            self.insert_at = self._target_codon_at()

        if guide_strand_same is not None:
            assert guide_strand_same in (True, False)
            self.guide_strand_same = guide_strand_same
        else:
            self.guide_strand_same = self._guide_strand_same()

    def __repr__(self):
        return "HDR('{}', '{}', '{}', {}, {}, {})".format(
            self.target_seq,
            self.hdr_seq,
            self.hdr_tag,
            self.hdr_dist,
            self.guide_strand_same,
            self._codon_at,
        )

    def _guide_strand_same(self) -> bool:
        """
        Infer guide direction by expected PAM locations.
        We try both directions because we don't know guide direction yet.
        There is a small chance that there could be PAMs equidistant in both
        directions.
        See get_guide_cut_to_insert.
        >>> hdr = HDR('GCCATGGCTGAGCTGGATCCGTTCGGC', hdr_dist=14)
        >>> hdr._guide_strand_same()
        True
        >>> hdr = HDR('GCCATGGCTGAGCTGGATCCGTTCGGC', hdr_dist=1)
        >>> hdr._guide_strand_same()
        False
        """

        cut_at = self.cut_at
        pam1 = self.target_seq[cut_at + 3:cut_at + 6].upper()
        pam2 = self.target_seq[cut_at - 6:cut_at - 3].upper()
        is_for = pam1.endswith('GG')
        is_rev = pam2.startswith('CC')
        assert is_for or is_rev, (pam1, pam2)
        assert not (is_for and is_rev)
        return True if is_for else False

    @functools.lru_cache(1024 * 1024)
    def _target_codon_at(self) -> int:
        # If codon position explicitly passed in, use that.
        if self._codon_at != -1:
            return self._codon_at
        for i, codon in enumerate(_left_to_right_codons(self.target_seq)):
            if codon.upper() in self.boundary_codons:
                return i * 3

        assert False

    @property
    def cut_at(self):
        cut_at = self.insert_at + self.hdr_dist
        assert cut_at >= 0, (self.insert_at, self.hdr_dist)
        return cut_at

    @property
    def cut_in_junction(self) -> bool:
        """
        Determines whether the cut location is inside an intron/exon junction,
        as previously marked out by lowercasing.
        Cut in junction.
        >>> hdr = HDR('GCCATGGCTGAGCTGGAtccgttCGGC', hdr_dist=14)
        >>> (hdr.cut_at, hdr.cut_in_junction)
        (20, True)
        Cut just after junction.
        >>> hdr = HDR('CCNNNNtaannnnnn', hdr_dist=0, hdr_tag='stop_codon', guide_strand_same=True)
        >>> (hdr.cut_at, hdr.cut_in_junction)
        (6, False)
        No junction to cut.
        >>> hdr = HDR('ATGNGG', hdr_dist=-3)
        >>> hdr.cut_in_junction
        False
        """
        return self.target_seq[self.cut_at - 1].islower()

    @property
    def guide_seq(self):
        """
        Returns 23bp guide sequence that includes PAM.
        >>> hdr = HDR('GCCATGGCTGAGCTGGATCCGTTCGGC', hdr_dist=14)
        >>> hdr.guide_seq
        'ATGGCTGAGCTGGATCCGTTCGG'
        >>> hdr = HDR('GCCATGGCTGAGCTGGATCCGTTCGGC', hdr_dist=1)
        >>> hdr.guide_seq
        'CCATGGCTGAGCTGGATCCGTTC'
        """
        cut_at = self.cut_at
        if self.guide_strand_same:
            guide_seq = self.target_seq[cut_at - 17:cut_at + 6]
        else:
            guide_seq = self.target_seq[cut_at - 6:cut_at + 17]
        assert len(guide_seq) == 23, cut_at
        return guide_seq

    def anchor_seq(self, size: int = 9) -> str:
        """
        Returns a codon-aligned sequence around the codon for finding the codon
        again in a primer sequence, which should itself be somewhat centered
        around the codon.
        >>> hdr = HDR('CCTTGGCTGATGTGGATCCGTTCGGC', hdr_dist=-12)
        >>> hdr.anchor_seq()
        'CCTTGGCTGATGTGGATC'
        """
        codon_at = self._target_codon_at()
        return self.target_seq[codon_at - size:codon_at + size]

    @property
    def inserted(self) -> str:
        """
        >>> hdr = HDR('GCCATGGCTGAGCTGGATCCGTTCGGC', 'NNN', hdr_dist=14)
        >>> hdr.inserted
        'GCCATGnnnGCTGAGCTGGATCCGTTCGGC'
        """
        return (
            self.target_seq[:self.insert_at] +
            self.hdr_seq.lower() +
            self.target_seq[self.insert_at:])

    @property
    def inserted_mutated(self) -> MutatedSeq:
        """
        >>> hdr = HDR('GCCATGGCTGAGCTGGATCCGTTCGGC', 'TTT', hdr_dist=14)
        >>> hdr.target_mutation_score = 0.01
        >>> hdr.inserted_mutated
        'GCCATGTTcGCcGAatTaGAcCCcTTtGGC'
        """
        return self._mutate(do_insert=True)

    @property
    def mutated(self) -> str:
        """
        Mutates target sequence. If the guide PAM is outside the coding region,
        the PAM is mutated in place. Otherwise, some codons in the guide are
        mutated silently.
        >>> hdr = HDR('GCCATGGCTGAGCTGGATCCGTTCGGC', hdr_dist=14)
        >>> hdr.target_mutation_score = 0.5
        >>> hdr.mutated
        'GCCATGGCTGAGCTGGATCCcTTCGGC'
        PAM is outside.
        >>> hdr = HDR('CCTTGGCTGATGTGGATCCGTTCGGC', hdr_dist=-12)
        >>> hdr.mutated
        'CCTTccCTGATGTGGATCCGTTCGGC'
        >>> hdr = HDR('TGACCTAGAGATTGCAAGGGCGGG', hdr_dist=9, guide_strand_same=False, hdr_tag='stop_codon')
        >>> hdr.pam_outside_cds
        True
        >>> hdr.mutated
        'TGAggTAGAGATTGCAAGGGCGGG'
        >>> hdr = HDR('TGATCCCAAATTTGTCCATAGCTGAAG', hdr_dist=10, guide_strand_same=False, hdr_tag='stop_codon')
        >>> hdr.pam_outside_cds
        True
        >>> hdr.mutated
        'TGATggCAAATTTGTCCATAGCTGAAG'
        """
        return self._mutate(do_insert=False)

    @property
    def _pam_mutated(self) -> str:
        """
        Target seq with 3bp PAM mutated inside it.
        # TODO (gdingle): fixme
        >> hdr = HDR('CCTTGGCTGATGTGGATCCGTTCGGC', hdr_dist=-12)
        >> hdr._pam_mutated
        'CCTTccCTGATGTGGATCCGTTCGGC'
        >> hdr = HDR('ATGCCTTGGCTGATATGGATCCGT', hdr_dist=6, guide_strand_same=False)
        >> hdr._pam_mutated
        'ATGggTTGGCTGATATGGATCCGT'
        """
        before, pam, after = (
            self.target_seq[:self.pam_at],
            self.target_seq[self.pam_at:self.pam_at + 3].upper(),
            self.target_seq[self.pam_at + 3:]
        )
        assert len(pam) == 3
        if self.guide_strand_same:
            assert 'GG' in pam, pam
            pam_mutated = pam.replace('GG', 'cc')
        else:
            assert 'CC' in pam, pam
            pam_mutated = pam.replace('CC', 'gg')
        combined = before + pam_mutated + after
        assert len(combined) == len(self.target_seq)
        return combined

    def _mutate(self, do_insert: bool = True) -> MutatedSeq:
        """
        Returns the target_seq with inserted hdr_seq after optimal mutation.
        >>> hdr = HDR('GCCATGGCTGAGCTGGATCCGTTCGGCTAT', 'aaa', hdr_dist=14)
        >>> hdr.guide_seq
        'ATGGCTGAGCTGGATCCGTTCGG'
        >>> hdr.target_mutation_score = 0.5
        >>> hdr.inserted
        'GCCATGaaaGCTGAGCTGGATCCGTTCGGCTAT'
        >>> hdr._mutate()
        'GCCATGAAAGCTGAGCTGGATCCcTTCGGCTAT'
        >>> hdr.target_mutation_score = 0.1
        >>> hdr.inserted
        'GCCATGaaaGCTGAGCTGGATCCGTTCGGCTAT'
        >>> hdr._mutate()
        'GCCATGAAAGCTGAGCTGGAcCCcTTCGGCTAT'
        >> hdr.inserted_mutated # legacy function result
        'GCCATGaaaGCTGAGCTGGATCCcTTtGGCTAT'
        # No insert
        >>> hdr = HDR('GCCATGGCTGAGCTGGATCCGTTCGGCTAT', '', hdr_dist=14)
        >>> hdr.target_mutation_score = 0.1
        >>> hdr._mutate()
        'GCCATGGCTGAaCTGGAcCCcTTCGGCTAT'
        # Mutates insert
        >>> hdr = HDR('GCCGCTGAGCTGGATCCGATGTTCGG', 'TTCGG', hdr_dist=-1)
        >>> hdr.target_mutation_score = 0.1
        >>> hdr.guide_seq
        'GCTGAGCTGGATCCGATGTTCGG'
        >>> hdr.inserted
        'GCCGCTGAGCTGGATCCGATGttcggTTCGG'
        >>> hdr._mutate()
        'GCCGCcGAGtTaGATCCcATGTTCGGTTCGG'
        # Test fallback, 21bp
        >>> hdr = HDR('GCCATGGCTGAGCTGGATCCGTTCGG', hdr_dist=14)
        >>> hdr.target_mutation_score = 0.01
        >>> hdr._mutate()
        'GCCATGGCcGAatTaGAcCCcTTtGG'
        # Mutate PAM only
        >>> hdr = HDR('CTCAGAAGATGATGACTGAAAGGGACTCGGGACT', 'atg', 'stop_codon', 9, True, 15)
        >>> hdr.guide_seq
        'GATGATGACTGAAAGGGACTCGG'
        >>> hdr.target_mutation_score = 0.03
        >>> hdr.should_mutate
        True
        >>> hdr.pam_outside_cds
        True
        >>> hdr._mutate()
        'CTCAGAAGATGATGAATGCTGAAAGGGACTCccGACT'
        """
        if self.pam_outside_cds and self.should_mutate:
            # Skip other kinds of mutations because PAM mutation is enough
            seq = self._pam_mutated
            if do_insert:
                seq = (seq[:self.insert_at]
                       + self.hdr_seq.upper()
                       + seq[self.insert_at:])
            return MutatedSeq(seq, max_score=0.0, max_seq='')

        length = self.guide_seq_aligned_length

        ret_seq = list(self.target_seq.upper())  # To mutate string in place

        if do_insert:
            ret_seq = (ret_seq[:self.insert_at]
                       + list(self.hdr_seq.upper())
                       + ret_seq[self.insert_at:])

        if len(ret_seq) < length:
            logging.warning('Cannot find {}bp for mutation in target_seq {}. Falling back to 21bp.'.format(
                self.guide_seq_aligned_length,
                self.target_seq,))
            length = 21
        test_length = min(len(self.guide_seq), length)
        assert test_length >= 20, self.guide_seq

        # TODO (gdingle): remove MIT score?
        if self.use_cfd_score:
            hit_score_func = functools.partial(
                cfdscore.cfd_score,
                wt=self.guide_seq.upper(),
                guide_strand_same=self.guide_strand_same)
        else:
            hit_score_func = functools.partial(
                mitscore.mit_hit_score,
                self.guide_seq.upper(),
                guide_strand_same=self.guide_strand_same,
                include_pam=test_length == 23)

        max_score = 0.0
        max_seq = ''

        # Iterate over codons
        for left in range(0, len(ret_seq) - length + 1, 3):
            right = left + length
            # Iterate within 27bp
            for start in range(0, 5):
                mutate_seq = ''.join(ret_seq[left:right])
                end = start + test_length
                if test_length == 23:
                    # mask outside of 23bp of guide seq
                    mutate_seq = (
                        mutate_seq[:start].lower() +
                        mutate_seq[start:end] +
                        mutate_seq[end:].lower()
                    )
                assert len(mutate_seq) == length, len(mutate_seq)

                mutated, score = _best_mutation(
                    mutate_seq,
                    self.guide_strand_same,
                    self.mutate_all_permutations,
                    start,
                    end,
                    hit_score_func,
                    self.target_mutation_score)

                for k, c in enumerate(mutated):
                    if c.upper() != self.inserted[left + k].upper():
                        ret_seq[left + k] = c.lower()

                max_score = max(score, max_score)
                if score == max_score:
                    max_seq = mutated[start:end]

        # TODO (gdingle): this doesn't always make sense
        if max_score > self.target_mutation_score:
            logger.warning('Unable to mutate enough. Max score {}, target score {}, max seq {}, guide seq {}'.format(
                max_score, self.target_mutation_score, max_seq, self.guide_seq))

        return MutatedSeq(''.join(ret_seq), max_score=max_score, max_seq=max_seq)

    @property
    def _mutated_score(self) -> float:
        """
        >>> hdr = HDR('GCCATGGCTGAGCTGGATCCGTTCGGCTAT', '', hdr_dist=14)
        >>> hdr.target_mutation_score = 0.1
        >>> hdr.mutated
        'GCCATGGCTGAaCTGGAcCCcTTCGGCTAT'
        >>> hdr._mutated_score
        0.08348794069944342
        """
        if self.pam_outside_cds and self.should_mutate:
            return 0

        return self._mutate().max_score

    @property
    def mutation_in_junction(self) -> bool:
        """
        Determines whether there is a mutation inside an intron/exon junction.
        # TODO (gdingle): use mutation masking here instead of warning?
        # Then we would need to preserve lowercasing in guide_seq_aligned
        # See https://trello.com/c/HoEcAlVj/54-filter-out-mutations-in-intron-exon-junction
        Mutation in junction.
        >>> hdr = HDR('CATATGatccggagCCCGCCCCGCCCCCGAGCCGCAT', hdr_dist=8, guide_strand_same=False)
        >>> hdr.guide_seq
        'ccggagCCCGCCCCGCCCCCGAG'
        >>> hdr.mutated
        'CATATGATtCGGAGCCCGCCCCGCCCCCGAGCCGCAT'
        >>> hdr.mutation_in_junction
        True
        No junction.
        >>> hdr = HDR('CATATGATCCGGAGCCCGCCCCGCCCCCGAGCCGCAT', hdr_dist=8, guide_strand_same=False)
        >>> hdr.mutation_in_junction
        False
        """
        if all(u.isupper() for u in self.target_seq):
            return False
        # we want only the lowercase intron/exons
        inserted = self.inserted.replace(self.hdr_seq.lower(), self.hdr_seq.upper())
        # TODO (gdingle): another perf problem!!! :(
        for i, c in enumerate(self.inserted_mutated):
            u = inserted[i]
            if u.islower() and u.upper() != c.upper():
                return True
        return False

    @property
    def should_mutate(self) -> bool:
        # TODO (gdingle): remove me if no longer needed because of cdf score
        """
        Determines whether a guide should be mutated depending on the cut to
        insert distance and the guide orientation. The rule is: mutate if more
        than Xbp of the PAM-side of protospacer will be intact after insertion.
        1a. 14bp or more intact on PAM side, positive guide.
        GCC|ATG|GCTGAGCTGGATCC|GTT|CGG|C
            codon              cut pam
        >>> hdr = HDR('GCCATGGCTGAGCTGGATCCGTTCGGC', hdr_dist=14)
        >>> hdr.should_mutate
        True
        1b. Less than 14bp intact on PAM side, positive guide.
        GCC|ATG|GCTGA|GTT|CGG|C
            codon           cut pam
        >>> hdr = HDR('GCCATGGAGCTGTTCGGC', hdr_dist=5)
        >>> hdr.should_mutate
        False
        2a. 14bp or more intact on PAM side, negative guide.
        |CCA|CGA|GCGGCGGCGGCG|ATG|
         pam cut              codon
        >>> hdr = HDR('CCACGAGCGGCGGCGGCGATG', hdr_dist=-15, guide_strand_same=False)
        >>> hdr.should_mutate
        True
        2b. Less than 14bp intact on PAM side, negative guide.
        |CCA|CGA|GCG|ATG|GCTGAGCTGGATCCG
         pam cut     codon
        >>> hdr = HDR('CCACGAGCGATGGCTGAGCTGGATCCG', hdr_dist=-6, guide_strand_same=False)
        >>> hdr.should_mutate
        False
        3a. Insert is outside of guide, positive guide.
        |CCT|TGG|CTG|ATG|TGGATCCGTTCGGC
         cut pam     codon
        >>> hdr = HDR('CCTTGGCTGATGTGGATCCGTTCGGC', hdr_dist=-12)
        >>> hdr.should_mutate
        True
        3b. Insert is outside of guide, negative guide.
        |ATG|CCT|TGG|CTGATATGGATCCGT
         cod pam cut
        >>> hdr = HDR('ATGCCTTGGCTGATATGGATCCGT', hdr_dist=6, guide_strand_same=False)
        >>> hdr.should_mutate
        True
        """
        if self.guide_strand_same:
            guide_right = self.cut_at + 3
            intact = guide_right - self.insert_at
        else:
            guide_left = self.cut_at - 3
            intact = self.insert_at - guide_left

        # intact <= -3 means the insert is outside the guide + pam
        # return intact <= -3 or intact >= 14
        # modified to be more lenient
        return intact <= -3 or intact >= 11

    @property
    def pam_at(self) -> int:
        """
        >>> hdr = HDR('ATGCCTTGGCTGATATGGATCCGT', hdr_dist=6, guide_strand_same=False)
        >>> hdr.pam_at
        3
        >> hdr = HDR('CCTTGGCTGATGTGGATCCGTTCGGC', hdr_dist=-12)
        >> hdr.pam_at
        3
        >>> hdr = HDR('TGACCTAGAGATTGCAAGGGCGGG', hdr_dist=9, guide_strand_same=False, hdr_tag='stop_codon')
        >>> hdr.pam_at
        3
        """
        if self.guide_strand_same:
            return self.cut_at + 3
        else:
            return self.cut_at - 6

    @property
    def pam_outside_cds(self) -> bool:
        """
        >>> hdr = HDR('ATGCCTTGGCTGATATGGATCCGT', hdr_dist=6, guide_strand_same=False)
        >>> hdr.pam_outside_cds
        False
        >>> hdr = HDR('CCTTGGCTGATGTGGATCCGTTCGGC', hdr_dist=-12)
        >>> hdr.pam_outside_cds
        True
        >>> hdr = HDR('TGACCTAGAGATTGCAAGGGCGGG', hdr_dist=9, guide_strand_same=False, hdr_tag='stop_codon')
        >>> hdr.pam_outside_cds
        True
        """
        if self.hdr_tag == 'start_codon':
            return self.pam_at <= self._target_codon_at() - 3
        else:
            return self.pam_at >= self._target_codon_at() + 3


def mutate_silently(
        guide_seq: str,
        guide_strand_same: bool=False,
        skip_start_stop_codon: bool=True,
        all_permutations: bool=False) -> Iterator[str]:
    """
    Generator that silently mutates input sequence by substituing a different
    codon that encodes the same amino acid. Changes one codon per iteration.
    Direction is from PAM inwards, unless all_permutations is True. The new
    codon is the selected by frequency in the human genome.
    Data from http://biopython.org/DIST/docs/api/Bio.SeqUtils.CodonUsage-pysrc.html
    The input is assumed to a multiple of 3bp codons.
    By default, does not mutate stop codons, because such mutations are not
    always silent.
    If all_permutations is True, all possible orderings of mutations will be
    returned.
    Lowercase letters in guide_seq will never be mutated. This is useful for
    masking inputs. The letters will be uppercased on return, to avoid
    ambiguity with mutated letters.
    >>> it = mutate_silently('TGTTGCGATGAC')
    >>> next(it)
    'TGcTGCGATGAC'
    >>> next(it)
    'TGcTGtGATGAC'
    No possible synonyms.
    >>> next(mutate_silently('ATG'))
    'ATG'
    Right to left.
    >>> it = mutate_silently('TGTTGCGATGAC', True)
    >>> next(it)
    'TGTTGCGATGAt'
    >>> next(it)
    'TGTTGCGAcGAt'
    Skip stop codon.
    >>> it = mutate_silently('TAG')
    >>> next(it)
    'TAG'
    >>> it = mutate_silently('TAG', skip_stop_codon=False)
    >>> next(it)
    'Tga'
    all_permutations
    >>> it = mutate_silently('TGTTGCGATGAC', all_permutations=True)
    >>> next(it) # first one no effect when all_permutations
    'TGTTGCGATGAC'
    >>> next(it)
    'TGcTGCGATGAC'
    >>> next(it)
    'TGTTGtGATGAC'
    >>> next(it)
    'TGTTGCGAcGAC'
    >>> next(it)
    'TGTTGCGATGAt'
    all_permutations, other direction
    >>> it = mutate_silently('TGTTGCGATGAC', all_permutations=True, guide_strand_same=True)
    >>> next(it)
    'TGTTGCGATGAC'
    >>> next(it)
    'TGTTGCGATGAt'
    >>> next(it)
    'TGTTGCGAcGAC'
    >>> next(it)
    'TGTTGtGATGAC'
    Lowercase masking.
    >>> it = mutate_silently('TGtTGTTgT')
    >>> next(it)
    'TGTTGTTGT'
    >>> next(it)
    'TGTTGcTGT'
    >>> next(it)
    'TGTTGcTGc'
    Problem case.
    >>> it = mutate_silently('tCCAGGTAGTGCCGCGCTGCCTGCacc', all_permutations=True)
    >>> next(it)
    'TCCAGGTAGTGCCGCGCTGCCTGCACC'
    """
    synonymous = {
        'CYS': ['TGT', 'TGC'],
        'ASP': ['GAT', 'GAC'],
        'SER': ['TCT', 'TCG', 'TCA', 'TCC', 'AGC', 'AGT'],
        'GLN': ['CAA', 'CAG'],
        'MET': ['ATG'],
        'ASN': ['AAC', 'AAT'],
        'PRO': ['CCT', 'CCG', 'CCA', 'CCC'],
        'LYS': ['AAG', 'AAA'],
        'STOP': ['TAG', 'TGA', 'TAA'],
        'THR': ['ACC', 'ACA', 'ACG', 'ACT'],
        'PHE': ['TTT', 'TTC'],
        'ALA': ['GCA', 'GCC', 'GCG', 'GCT'],
        'GLY': ['GGT', 'GGG', 'GGA', 'GGC'],
        'ILE': ['ATC', 'ATA', 'ATT'],
        'LEU': ['TTA', 'TTG', 'CTC', 'CTT', 'CTG', 'CTA'],
        'HIS': ['CAT', 'CAC'],
        'ARG': ['CGA', 'CGC', 'CGG', 'CGT', 'AGG', 'AGA'],
        'TRP': ['TGG'],
        'VAL': ['GTA', 'GTC', 'GTG', 'GTT'],
        'GLU': ['GAG', 'GAA'],
        'TYR': ['TAT', 'TAC'],
    }
    # Fraction of occurences among synonyms in Human genome.
    # see https://www.genscript.com/tools/codon-frequency-table
    syn_fractions = {
        'ATG': 1, 'TGG': 1, 'CAG': 0.75, 'CAC': 0.59, 'AAG': 0.58, 'GAG': 0.58,
        'TAC': 0.57, 'TTC': 0.55, 'TGC': 0.55, 'AAC': 0.54, 'GAC': 0.54, 'TGA':
        0.52, 'ATC': 0.48, 'GTG': 0.47, 'AAT': 0.46, 'GAT': 0.46, 'TTT': 0.45,
        'TGT': 0.45, 'TAT': 0.43, 'AAA': 0.42, 'GAA': 0.42, 'CTG': 0.41, 'CAT':
        0.41, 'GCC': 0.4, 'ATT': 0.36, 'ACC': 0.36, 'GGC': 0.34, 'CCC': 0.33,
        'TAA': 0.28, 'CCT': 0.28, 'ACA': 0.28, 'CCA': 0.27, 'GCT': 0.26, 'CAA':
        0.25, 'GGA': 0.25, 'GGG': 0.25, 'GTC': 0.24, 'ACT': 0.24, 'AGC': 0.27,
        'GCA': 0.23, 'TCC': 0.22, 'CGG': 0.21, 'TAG': 0.2, 'CTC': 0.2, 'AGA':
        0.2, 'AGG': 0.2, 'CGC': 0.19, 'GTT': 0.18, 'TCT': 0.18, 'ATA': 0.16,
        'GGT': 0.16, 'TCA': 0.15, 'AGT': 0.15, 'TTG': 0.13, 'CTT': 0.13, 'ACG':
        0.12, 'GTA': 0.11, 'CCG': 0.11, 'CGA': 0.11, 'GCG': 0.11, 'CGT': 0.08,
        'TTA': 0.07, 'CTA': 0.07, 'TCG': 0.06,
    }
    # Rare codons in human sequences (defined by: frequency less than 7.0e-3 AND
    # less than half of median usage for that amino acid):
    blacklist = {
        'TCG': 'SER',
        'CCG': 'PRO',
        'ACG': 'THR',
        'GCG': 'ALA',
        'CGT': 'ARG',
        'ATA': 'ILE',
    }
    synonymous_index = dict(
        (codon, aa)
        for aa, codons in synonymous.items()
        for codon in codons
    )
    _validate_seq(guide_seq)

    def _mutate_codons(codons) -> Iterator:
        mutated_codons = []
        for codon in codons:
            mutated_codons.append(_mutate_codon(codon))
            yield mutated_codons

    @functools.lru_cache(maxsize=1024 * 1024)
    def _mutate_codon(codon: str) -> str:
        # Make copy and remove current codon
        syns = list(synonymous[synonymous_index[codon.upper()]])
        # Remove self codon
        syns.remove(codon.upper())
        # Filter out blacklist
        syns = [syn for syn in syns if not blacklist.get(syn)]

        # Skip mutations that affect lowercase masked base pairs
        for syn in syns.copy():
            for i, c in enumerate(codon):
                # if masked and mutated, skip
                if c != c.upper() and c.upper() != syn[i] and syn in syns:
                    syns.remove(syn)

        if skip_start_stop_codon and codon in ['TAG', 'TGA', 'TAA', 'ATG']:
            return codon
        elif len(syns):
            top = _select_syn(codon, syns)
            lowered = ''.join([
                c.upper() if c.upper() == top[i] else top[i].lower()
                for i, c in enumerate(codon)
            ])
            return lowered
        else:
            return codon.upper()  # erase lowercase masking

    # TODO (gdingle): make this work with hashable input
    # @functools.lru_cache(maxsize=1024 * 1024)
    def _select_syn(codon: str, syns: list) -> str:
        """
        Selects the most different by base pairs, or the most frequent in
        the genome.
        >>> syns = ['TCT', 'TCG', 'TCA', 'TCC', 'AGC', 'AGT']
        >>> _select_syn('TCT', syns)
        'AGC'
        """
        codon = codon.upper()
        diffs = [(sum([
            codon[0] != syn[0],
            codon[1] != syn[1],
            codon[2] != syn[2],
        ]), syn) for syn in syns]

        if all(diffs[0][0] == s for s, syn in diffs):
            # no differences, use old method
            fractions = tuple((syn_fractions[syn], syn) for syn in syns)
            top = max(fractions)[1]
        else:
            top = max(diffs)[1]

        return top

    def _all_permutations(guide_seq) -> Iterator:
        """This will return increasing numbers of mutations, from right to left,
        or left to right depending on strand."""
        codons = list(_left_to_right_codons(guide_seq))
        masks = itertools.product([False, True], repeat=len(codons))
        # sort to ensure strictly increasing number of mutations

        for mask in sorted(masks, key=lambda m: sum(m)):
            assert len(mask) == len(codons)
            if not guide_strand_same:
                mask = mask[::-1]
            new_codons = []
            for i, do_mutate in enumerate(mask):
                if do_mutate:
                    new_codons.append(_mutate_codon(codons[i]))
                else:
                    new_codons.append(codons[i].upper())  # erase lowercase masking
            assert len(new_codons) == len(codons)
            yield ''.join(new_codons)

    def _pam_inwards(guide_seq) -> Iterator:
        if guide_strand_same:
            codons = _right_to_left_codons(guide_seq)
        else:
            codons = _left_to_right_codons(guide_seq)

        for mutated_codons in _mutate_codons(codons):
            guide_seq = guide_seq.upper()  # erase lowercase masking
            if guide_strand_same:
                new_guide_str = ''.join(mutated_codons[::-1])
                combined = guide_seq[:-len(new_guide_str)] + new_guide_str
            else:
                new_guide_str = ''.join(mutated_codons)
                combined = new_guide_str + guide_seq[len(new_guide_str):]

            assert len(combined) == len(guide_seq), (combined, guide_seq)
            yield combined

    if all_permutations:
        yield from _all_permutations(guide_seq)
    else:
        yield from _pam_inwards(guide_seq)


@functools.lru_cache(maxsize=1024 * 1024)
def _validate_seq(seq: str):
    assert all(b.upper() in 'AGCTN' for b in seq), seq
    if seq != '':
        assert len(seq) >= 3, seq


def _right_to_left_codons(seq: str) -> Iterator[str]:
    """
    >>> it = _right_to_left_codons('TGTTGCGATGAC')
    >>> next(it)
    'GAC'
    >>> next(it)
    'GAT'
    >>> next(it)
    'TGC'
    >>> next(it)
    'TGT'
    >>> next(it)
    Traceback (most recent call last):
    ...
    StopIteration
    """
    for i in range(len(seq), 0, -3):
        codon = seq[i - 3:i]
        yield codon


def _left_to_right_codons(seq: str) -> Iterator[str]:
    """
    >>> it = _left_to_right_codons('TGTTGCGATGAC')
    >>> next(it)
    'TGT'
    >>> next(it)
    'TGC'
    """
    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3]
        yield codon


@functools.lru_cache(1024 * 1024)
def _best_mutation(
    mutate_seq: str,
    guide_strand_same: bool,
    all_permutations: bool,
    start: int,
    end: int,
    hit_score_func: functools.partial,
    target_mutation_score: float,
) -> tuple:
    for mutated in mutate_silently(
        mutate_seq,
        guide_strand_same,
        True,
        all_permutations
    ):
        # note: no change on first pass
        assert len(mutated) == len(mutate_seq), len(mutated)
        if end - start == 23:
            mutated_test_seq = mutated[start:end].upper()
        else:
            mutated_test_seq = mutated.upper()

        score = hit_score_func(sg=mutated_test_seq)

        # TODO (gdingle): also check rev comp?
        # rmutated_test_seq = cfdscore._revcom(mutated_test_seq)

        if score <= target_mutation_score:
            return mutated, score

    return mutated, score


if __name__ == '__main__':
    import doctest
    doctest.testmod()