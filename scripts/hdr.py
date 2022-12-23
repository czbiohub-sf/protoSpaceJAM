import functools
import itertools
import logging
from typing import List
import sys
from typing import Iterator
import copy
from Bio.Seq import Seq
from Bio import Restriction
from Bio.SeqUtils import GC
import re
from scripts.cfdscore import *
from itertools import islice
import math


try:
    from . import cfdscore, mitscore
except ImportError:
    import cfdscore  # type: ignore
    import mitscore  # type: ignore

logger = logging.getLogger(__name__)

#by DP
class seq_w_phase:
    '''
    seqs with phases and start + end info
    '''
    def __init__(self,
                 seq:str, phases:str, start:int, end:int)->None:
        self.seq = seq
        self.phases = phases
        self.start = start
        self.end = end

class seq_w_phase_cfd:
    '''
    seqs with phases and start + end info
    '''
    def __init__(self,
                 seq:str, phases:str, cfd:float)->None:
        self.seq = seq
        self.phases = phases
        self.cfd = cfd

class chimeric_gRNA:
    '''
    gRNA   .seq .cfd .phases
    recut  .seq .cfd .phases
    recut_rc .seq .cfd .phases
    '''
    def __init__(self,
                 gRNA_seq:str, gRNA_phases:str, recut_seq:str, recut_cfd:float,recut_phases:str)->None:
        self.gRNA = seq_w_phase_cfd(seq = gRNA_seq, phases = gRNA_phases, cfd = 1)
        self.recut = seq_w_phase_cfd(seq = recut_seq, phases = recut_phases, cfd = recut_cfd)

    def update(self,newSeq):
        """
        update newSeq only if cfd score is lower than the current one
        """
        newCfd = cfd_score(self.gRNA.seq, newSeq)
        if newCfd < self.recut.cfd:
            self.recut.seq = newSeq
            self.recut.cfd = newCfd

class HDR_flank:
    """
    useful objects:
    .gRNA: instance of the class seq_w_phase, gRNA is in the format of 20nt-NGG
    .chimeric_gRNA: str, chimeric gRNA
    .chimeric_gRNA_phases: str, phases of the chimeric gRNA
    """
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
            ENST_chr: str,
            gStart : int,
            gStrand : int,
            InsPos:int,
            CutPos:int,
            Cut2Ins_dist:int,
            tag:str,
            name:str,
            ssDNA_max_size,
            loc2posType,
            recoding_args,
            Donor_type,
            Strand_choice,
            syn_check_args) -> None:

        self.left_flk_seq = left_flk_seq
        self.right_flk_seq = right_flk_seq
        self.ENST_ID = ENST_ID
        self.ENST_strand = ENST_strand
        self.ENST_chr = ENST_chr
        self.gStrand = gStrand
        self.gStart = gStart
        self.InsPos = InsPos # InsPos is the first letter of stop codon "T"AA or the last letter of the start codon AT"G"
        self.CutPos = CutPos
        self.Cut2Ins_dist = Cut2Ins_dist
        self.tag = tag
        self.loc2posType = loc2posType
        self.name = name
        self.ssDNA_max_size = ssDNA_max_size
        self.recoding_args = recoding_args
        self.cfdThres = recoding_args["cfdThres"]
        self.synFlags=[]
        self.Donor_type=Donor_type
        self.Strand_choice=Strand_choice
        self.enzymes2check = syn_check_args["check_enzymes"]
        self.CustomSeq2Avoid = syn_check_args["CustomSeq2Avoid"]
        self.recode_order = recoding_args["recode_order"]

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

        if gStrand * ENST_strand >= 0:
            self.gRNA_in_coding_strand = True
        else:
            self.gRNA_in_coding_strand = False

        use_junction_masked_phase = True

        #print(f"{self.left_flk_phases}||{self.right_flk_phases}\n{self.left_flk_coord_lst[0]}-{self.left_flk_coord_lst[1]}||{self.right_flk_coord_lst[0]}-{self.right_flk_coord_lst[1]}")

        #for positions within 3bp to exon/intron junctions, change phase from "0" to "9", for cds, change from "1/2/3" to "8"
        #for position in 5UTR, change phase from "0" to "5"
        coord1, coord2 = self.left_flk_coord_lst[0], self.left_flk_coord_lst[1]
        self.left_flk_phases_masked = self.mask_phase_within_3bp_exon_intron_junction(coord1, coord2, self.left_flk_phases)
        self.left_flk_phases_masked = self.mask_phase_in_5UTR(coord1, coord2, self.left_flk_phases_masked)
        coord1, coord2 = self.right_flk_coord_lst[0], self.right_flk_coord_lst[1]
        self.right_flk_phases_masked = self.mask_phase_within_3bp_exon_intron_junction(coord1, coord2, self.right_flk_phases)
        self.right_flk_phases_masked = self.mask_phase_in_5UTR(coord1, coord2, self.right_flk_phases_masked)


        # adjust insPos, for stop-tagging ,the insertion site is now before the first base of the stop codon
        if type == "stop":
            if self.ENST_strand == 1 or self.ENST_strand == "1" or self.ENST_strand == "+":
                self.InsPos = int(self.InsPos) - 1
            else:
                self.InsPos = int(self.InsPos) + 1
        # adjust cutPos for -1 strands genes (the gRNA cut site was based on +1 strand), we need to view the gene in coding sequence
        if self.ENST_strand == -1 or self.ENST_strand == "-1" or self.ENST_strand == "-":
            self.CutPos = self.CutPos + 1

        #TODO: check this line
        self.cutPos, self.gPAM_end = self.get_gRNA_pos(self.gStart, self.gStrand)

        #check if the whole 23nt gRNA is in the HDR flank
        self.entire_gRNA_in_HDR_arms, self.gStart_in_HDR_arms, self.gPAM_end_in_HDR_arms = self.check_gRNA_in_HDR_arms(gStart=self.gStart, gPAM_end=self.gPAM_end, left_flk_coord_lst=self.left_flk_coord_lst, right_flk_coord_lst=self.right_flk_coord_lst)

        #make gRNA lowercase
        self.gRNA_lc_Larm, self.gRNA_lc_Rarm = self.make_gRNA_lowercase()

        #use masked phases (5UTR, 3bp exon-intron junction)
        #print(f"{self.left_flk_phases_masked}||{self.right_flk_phases_masked}\n{self.gRNA_lc_Larm}||{self.gRNA_lc_Rarm}\n{self.left_flk_coord_lst[0]}-{self.left_flk_coord_lst[1]}||{self.right_flk_coord_lst[0]}-{self.right_flk_coord_lst[1]}")
        if use_junction_masked_phase == True:
            self.left_flk_phases = self.left_flk_phases_masked
            self.right_flk_phases = self.right_flk_phases_masked

        #get gRNA left/right coord respective to the whole arm (gRNA orientation can be both -1 and +1)
        null, self.g_leftcoord, self.g_rightcoord = self.get_trunc_gRNA(self.left_flk_seq, self.right_flk_seq)

        #get seq and phase between insertion and cut site
        self.ins2cut = self.get_ins2cut_seq()

        #check if the ins2cut_seq is the same length as the Cut2Ins_dist calculated elsewhere
        if len(self.ins2cut.seq)>0 and (len(self.ins2cut.seq)-4 != abs(self.Cut2Ins_dist)):
            sys.exit(f"ins2cut_seq:{self.ins2cut.seq} is not the same length as reported: Cut2Ins_dist={self.Cut2Ins_dist}")

        # log information #
        self.info = f"\n#################\n#{self.ENST_ID}#\n#################\n {self.ENST_ID}\t{self.name}\tstrand:{self.ENST_strand}\tgRNA_strand:{self.gStrand}\t{type}-tagging\tCut2Ins-dist:{self.Cut2Ins_dist}\ngRNA:{self.gStart}-{self.gPAM_end} ({self.g_leftcoord}-{self.g_rightcoord})\tCutPos:{self.CutPos}\tInsPos:{self.InsPos}\n"
        self.info_arm = "".join(
            f"--------------------HDR arms (in coding strand)-----------------------------------------------------------------------------------\n"
            f"1. left | right arms:{self.gRNA_lc_Larm}||{self.gRNA_lc_Rarm}\n"
            f"2. Phases           :{self.left_flk_phases}||{self.right_flk_phases}\n"
            f"3. Coordinates      :\t{self.left_flk_coord_lst[0]}-{self.left_flk_coord_lst[1]} || {self.right_flk_coord_lst[0]}-{self.right_flk_coord_lst[1]}\n")

        ####################################################################################
        #start recoding
        #phase 1: Silently mutate sequence between insert and cut
        #         check CFD post payload integration
        #         Arm: left_flk_seq_CodonMut right_flk_seq_CodonMut
        #         gRNA: self.post_mut_ins_gRNA_seq (the returned gRNA is in strand same as the PAM)

        #phase 2: if cfd > self.cfdThres, mutate PAM and protospacer if in 3UTR/intron  (tunable order PAM <=> protospacer)
        #         Arm: left_flk_seq_CodonMut3 right_flk_seq_CodonMut2
        #         gRNA: self.post_mut2_gRNA_seq (the returned gRNA is in strand same as the PAM)

        #phase 3: if cfd > self.cfdThres, silently mutate gRNA seq not covered between insert and cut
        #         Arm: left_flk_seq_CodonMut2 right_flk_seq_CodonMut3
        #         gRNA: self.post_mut3_gRNA_seq

        #phase 4: if cfd > self.cfdThres, mutate protospacer and PAM if in 5UTR  (tunable order PAM <=> protospacer)
        #         Arm: left_flk_seq_CodonMut4 right_flk_seq_CodonMut4
        #         gRNA: self.post_mut4_gRNA_seq (the returned gRNA is in strand same as the PAM)

        #phase 5: sliding window check of recutting

        #####################################################################################
        self.info_p1=''
        self.info_p2=''
        self.info_p3=''
        self.info_p4=''
        self.info_p5=''

        if not self.recoding_args["recoding_off"]:
            self.left_flk_seq_CodonMut = self.left_flk_seq #initialize for later use ( the definition may be skipped in phase 1)
            self.right_flk_seq_CodonMut = self.right_flk_seq #initialize for later use ( the definition may be skipped in phase 1)
            if not self.recoding_args["recoding_stop_recut_only"]:
                #########
                #Phase 1#
                #########
                #trim insert-to-cut into frame
                self.ins2cut_Ltrimed = self.trim_left_into_frame(self.ins2cut)
                self.ins2cut_LRtrimed = self.trim_right_into_frame(self.ins2cut_Ltrimed)

                # mutate insert-to-cut sequence
                mutated_subseq = self.get_silent_mutations(self.ins2cut_LRtrimed.seq.upper()) # mutated_subseq is always in coding strand

                #put_silent_mutation_subseq_back into HDR flank
                self.left_flk_seq_CodonMut, self.right_flk_seq_CodonMut = self.put_silent_mutation_subseq_back(L_arm=self.left_flk_seq, R_arm=self.right_flk_seq,
                                                                                                               mutated_subseq = mutated_subseq,                                      # mutated_subseq is in the coding strand
                                                                                                               start = self.ins2cut_LRtrimed.start, end = self.ins2cut_LRtrimed.end) #|-> start, end are local to the whole arm = L_arm + R_arm <-|#
            #check if gRNA is affected by insertion
            self.gRNA_seq, Null, self.gRNA_seq_phases, Null = self.get_post_integration_gRNA(self.left_flk_seq,self.right_flk_seq) #get the original gRNA
            Null, self.post_mut_ins_gRNA_seq, Null, self.post_mut_ins_gRNA_seq_phases = self.get_post_integration_gRNA(self.left_flk_seq_CodonMut, self.right_flk_seq_CodonMut) #get the post mutation and insertion gRNA
            #print(f"caculating CFD scores using: {self.gRNA_seq} {self.post_mut_ins_gRNA_seq}")

            ###############################################################
            #only uncomment this part to debut the short-HDR-crash problem#
            ###############################################################

            # self.info_p1 = "".join(
            #         f"--------------------phase 1 mutate seq between cut to insert----------------------------------------------------------------------\n"
            #         f"phase1.cut-to-insert + 2bp padding on both sides (extend to full codons)\n"
            #         f"phase1.seq                             :{self.ins2cut.seq}\n"
            #         f"phase1.Phases                          :{self.ins2cut.phases}\n"
            #         f"phase1.Coordinates                     :{self.ins2cut.start}-{self.ins2cut.end}\n"
            #         f"phase1.seq         (trimmed into frame):{self.ins2cut_LRtrimed.seq}\n"
            #         f"phase1.Phases      (trimmed into frame):{self.ins2cut_LRtrimed.phases}\n"
            #         f"phase1.Coordinates (trimmed into frame):{self.ins2cut_LRtrimed.start}-{self.ins2cut_LRtrimed.end}\n"
            #         f"phase1.mutated seq (trimmed into frame):{mutated_subseq}\n"
            #         f"--> display gRNA and check disruption <--\n"
            #         f"phase1.                      gRNA:{self.gRNA_seq}\n"
            #         f"phase1.                    Phases:{self.gRNA_seq_phases}\n"
            #         f"phase1.after mutation and payload:{self.post_mut_ins_gRNA_seq}\n"
            #         f"phase1.                    Phases:{self.post_mut_ins_gRNA_seq_phases}\n"
            #         #f"phase1.                       CFD:{self.cdf_score_post_mut_ins:.4f}\n"
            #         f"--> display mutation in HDR arms <--\n"
            #         f"phase1.original arms   :{self.gRNA_lc_Larm}||{self.gRNA_lc_Rarm}\n"
            #         f"phase1.cut2insert mut  :{self.left_flk_seq_CodonMut}||{self.right_flk_seq_CodonMut}\n"
            #         f"phase1.arms with tag   :{self.left_flk_seq_CodonMut}|{self.tag}|{self.right_flk_seq_CodonMut}\n")
            #print(self.info)
            #print(self.info_arm)
            #print(self.info_p1)

            self.cdf_score_post_mut_ins  = cfd_score(self.gRNA_seq, self.post_mut_ins_gRNA_seq)

            #########
            #phase 2#
            #########
            #TODO test(implement adjustable order of mutating PAM <=> protospacer)
            if self.cdf_score_post_mut_ins > self.cfdThres: # mutate protospacer if in 3' UTR or intron #

                # mutate PAM if PAM first
                left, right, null, seq, phases = self.get_uptodate_mut() #get up-to-date gRNA seq and phases
                self.post_mut2_gRNA_seq = seq
                self.post_mut2_gRNA_seq_phases = phases
                if self.recode_order == "PAM_first":
                    #get PAM position types (3UTR etc)
                    PAM_coords = self.get_pam_loc()
                    PAM_pos_types = []
                    for pos in PAM_coords:
                        position_type = _get_position_type(chr = self.ENST_chr, ID = self.ENST_ID, pos = pos, loc2posType = self.loc2posType)
                        PAM_pos_types = PAM_pos_types + position_type
                    left, right, cdf, seq, phases = self.get_uptodate_mut()  # get up-to-date gRNA seq and phases
                    self.post_mut2_gRNA_seq = seq
                    self.post_mut2_gRNA_seq_phases = phases
                    #mutate PAM if in 3UTR/intron
                    if (phases[-1:] == "0" or phases[-2:-1] == "0"): #   unnecessary to double-check for 3UTR/intron set(PAM_pos_types) == set(["3UTR"])
                        # mutate PAM if it's in 3 UTR (phase == 0)
                        if phases[-1:] == "0": #last position phase = 0 (the second G in NGG)
                            self.post_mut2_gRNA_seq = seq[:-1] + "c"
                        if phases[-2:-1] == "0": #second position phase = 0 (the first G in NGG)
                            self.post_mut2_gRNA_seq = self.post_mut2_gRNA_seq[:-2] + "c" + self.post_mut2_gRNA_seq[-1:]
                        self.cdf_score_post_mut2 = cfd_score(self.gRNA_seq, self.post_mut2_gRNA_seq)

                        #put disrupted-and-mutated seq back to arms
                        ssODN = f"{left}{self.tag}{right}"
                        ssODN = ssODN.replace(seq,self.post_mut2_gRNA_seq)
                        ssODN = ssODN.replace(self.revcom(seq),self.revcom(self.post_mut2_gRNA_seq))
                        self.left_flk_seq_CodonMut2 = ssODN[0:len(self.left_flk_seq)]
                        self.right_flk_seq_CodonMut2 =ssODN[-len(self.left_flk_seq):]

                #mutate protospacer
                left, right, null, seq, phases = self.get_uptodate_mut() #get up-to-date gRNA seq and phases
                self.post_mut2_gRNA_seq = seq
                self.post_mut2_gRNA_seq_phases = phases
                counter = -1
                latest_cfd = self.cdf_score_post_mut_ins
                if hasattr(self,"cdf_score_post_mut2"):
                    latest_cfd = self.cdf_score_post_mut2
                if latest_cfd > self.cfdThres:
                    for idx,item in reversed(list(enumerate(phases))):
                        if (idx == (len(phases) - 1) or idx == (len(phases) - 2)):
                            continue #skip PAM
                        if item == "0": # "0" means 3'UTR (5'UTR is labeled "5")
                            counter += 1
                            if counter % 3 != 0:
                                continue #  mutate 1 in every 3 bp
                            base = seq[idx]
                            mutbase = self.single_base_muation(base)
                            self.post_mut2_gRNA_seq = self.post_mut2_gRNA_seq[:idx] + mutbase + self.post_mut2_gRNA_seq[idx+1:]
                            self.cdf_score_post_mut2 = cfd_score(self.gRNA_seq, self.post_mut2_gRNA_seq)
                            #early stop if CFD goes below self.cfdThres
                            if self.cdf_score_post_mut2<self.cfdThres:
                                break
                        #put mutated seq back to arms
                        ssODN = f"{left}{self.tag}{right}"
                        ssODN = ssODN.replace(seq,self.post_mut2_gRNA_seq)
                        ssODN = ssODN.replace(self.revcom(seq),self.revcom(self.post_mut2_gRNA_seq))
                        self.left_flk_seq_CodonMut2 = ssODN[0:len(self.left_flk_seq)]
                        self.right_flk_seq_CodonMut2 =ssODN[-len(self.left_flk_seq):]

                # mutate PAM afterwards (if protospacer first)
                latest_cfd = self.cdf_score_post_mut_ins
                if hasattr(self,"cdf_score_post_mut2"):
                    latest_cfd = self.cdf_score_post_mut2
                if latest_cfd > self.cfdThres and self.recode_order == "protospacer_first": # mutate PAM if it was skipped earlier
                    #get PAM position types (3UTR etc)
                    PAM_coords = self.get_pam_loc()
                    PAM_pos_types = []
                    for pos in PAM_coords:
                        position_type = _get_position_type(chr = self.ENST_chr, ID = self.ENST_ID, pos = pos, loc2posType = self.loc2posType)
                        PAM_pos_types = PAM_pos_types + position_type
                    left, right, cdf, seq, phases = self.get_uptodate_mut()  # get up-to-date gRNA seq and phases
                    self.post_mut2_gRNA_seq = seq
                    self.post_mut2_gRNA_seq_phases = phases
                    #mutate PAM if in 3UTR/intron
                    if (phases[-1:] == "0" or phases[-2:-1] == "0"): #   unnecessary to double-check for 3UTR/intron set(PAM_pos_types) == set(["3UTR"])
                        # mutate PAM if it's in 3 UTR (phase == 0)
                        if phases[-1:] == "0": #last position phase = 0 (the second G in NGG)
                            self.post_mut2_gRNA_seq = seq[:-1] + "c"
                        if phases[-2:-1] == "0": #second position phase = 0 (the first G in NGG)
                            self.post_mut2_gRNA_seq = self.post_mut2_gRNA_seq[:-2] + "c" + self.post_mut2_gRNA_seq[-1:]
                        self.cdf_score_post_mut2 = cfd_score(self.gRNA_seq, self.post_mut2_gRNA_seq)

                        #put disrupted-and-mutated seq back to arms
                        ssODN = f"{left}{self.tag}{right}"
                        ssODN = ssODN.replace(seq,self.post_mut2_gRNA_seq)
                        ssODN = ssODN.replace(self.revcom(seq),self.revcom(self.post_mut2_gRNA_seq))
                        self.left_flk_seq_CodonMut2 = ssODN[0:len(self.left_flk_seq)]
                        self.right_flk_seq_CodonMut2 =ssODN[-len(self.left_flk_seq):]

            left,right,cfd,seq,phases = self.get_uptodate_mut()
            self.postphase2ODN = left + tag + right
            #########
            #phase 3#
            #########
            left, right, cdf, seq, phases = self.get_uptodate_mut()  # get up-to-date gRNA seq and phases
            if cdf > self.cfdThres:   # CDS>=self.cfdThres, try mutating gRNA silently(codons) in parts not covered between insert and cut

                #print(f"phase 3: gRNA in coding {self.gRNA_in_coding_strand}\n{seq}\n{phases}")

                #get gRNA in coding strand
                if self.gRNA_in_coding_strand == True:
                    seq_obj = seq_w_phase(seq = seq, phases =phases, start = 1, end = 23)
                    seq_obj_ltrim = self.trim_left_into_frame(seq_obj)
                    seq_obj_lrtrim = self.trim_right_into_frame(seq_obj_ltrim)
                    #print(f"phase 3 trimmed:\n{seq_obj_lrtrim.seq}\n{seq_obj_lrtrim.phases}")
                    if len(seq_obj_lrtrim.seq) >= 3:
                        mutated = self.get_silent_mutations(seq_obj_lrtrim.seq)
                    else:
                        mutated = seq_obj_lrtrim.seq
                    untrimmed = seq_obj.seq.replace(seq_obj_lrtrim.seq, mutated) # untrim: replace trimmed part with the mutated part
                    #print(f"phase 3 mutated:\n{mutated}\nuntrimmed:\n{untrimmed}")
                else:
                    seq_obj = seq_w_phase(seq = self.revcom(seq), phases = phases[::-1], start = 1, end = 23)
                    seq_obj_ltrim = self.trim_left_into_frame(seq_obj)
                    seq_obj_lrtrim = self.trim_right_into_frame(seq_obj_ltrim)
                    #print(f"phase 3 trimmed:\n{seq_obj_lrtrim.seq}\n{seq_obj_lrtrim.phases}")
                    if len(seq_obj_lrtrim.seq) >= 3:
                        mutated = self.get_silent_mutations(seq_obj_lrtrim.seq)
                    else:
                        mutated = seq_obj_lrtrim.seq
                    untrimmed = seq_obj.seq.replace(seq_obj_lrtrim.seq, mutated) # untrim: replace trimmed part with the mutated part
                    untrimmed = self.revcom(untrimmed)
                    #print(f"phase 3 mutated:\n{mutated}\nuntrimmed:\n{untrimmed}")

                # # get the sequence between PAM and cut site
                # self.trunc_gRNA, null, null = self.get_trunc_gRNA(left, right)
                # # trim insert-to-cut into frame
                # self.trunc_gRNA_Ltrimed = self.trim_left_into_frame(self.trunc_gRNA)
                # self.trunc_gRNA_LRtrimed = self.trim_right_into_frame(self.trunc_gRNA_Ltrimed)
                #
                # if len(self.trunc_gRNA_LRtrimed.seq)>=3:
                #     self.mutated_trunc_gRNA = self.get_silent_mutations(self.trunc_gRNA_LRtrimed.seq)
                # else:
                #     self.mutated_trunc_gRNA = self.trunc_gRNA_LRtrimed.seq

                #put mutated seq back to arms
                ssODN = f"{left}{self.tag}{right}"
                ssODN = ssODN.replace(seq,untrimmed)
                ssODN = ssODN.replace(self.revcom(seq),self.revcom(untrimmed))
                self.left_flk_seq_CodonMut3 = ssODN[0:len(self.left_flk_seq)]
                self.right_flk_seq_CodonMut3 =ssODN[-len(self.left_flk_seq):]
                #put mutated seq back to arms
                # self.left_flk_seq_CodonMut3, self.right_flk_seq_CodonMut3 = self.put_silent_mutation_subseq_back(L_arm=self.left_flk_seq_CodonMut, R_arm=self.right_flk_seq_CodonMut,
                #                                                                                                  mutated_subseq = self.mutated_trunc_gRNA,                                   # mutated_subseq is already in coding strand,
                #                                                                                                  start = self.trunc_gRNA_LRtrimed.start, end = self.trunc_gRNA_LRtrimed.end) #|-> start, end are local to the whole arm = L_arm + R_arm <-| #TODO: fix
                #
                # extract gRNA
                Null, self.post_mut3_gRNA_seq, Null, self.post_mut3_gRNA_seq_phases = self.get_post_integration_gRNA(self.left_flk_seq_CodonMut3, self.right_flk_seq_CodonMut3) #get the gRNA after mutating sequence
                #check CFD
                self.cdf_score_post_mut3  = cfd_score(self.gRNA_seq, self.post_mut3_gRNA_seq)
            left,right,cfd,seq,phases = self.get_uptodate_mut()
            self.postphase3ODN = left + tag + right
            #########
            #phase 4#
            #########
            #TODO test(implement adjustable order of mutating PAM <=> protospacer)
            #If CFD>self.cfdThres, mutate PAM & protospacer in 5’ UTR
            left, right, cfd, seq, phases = self.get_uptodate_mut()  # get up-to-date gRNA seq and phases
            if cfd >= self.cfdThres: #

                #mutate PAM if PAM first
                latest_cfd = cfd
                if hasattr(self,"cdf_score_post_mut4"):
                    latest_cfd = self.cdf_score_post_mut4
                if latest_cfd > self.cfdThres and self.recode_order == "PAM_first":
                    left, right, cdf, seq, phases = self.get_uptodate_mut()  # get up-to-date gRNA seq and phases
                    self.post_mut4_gRNA_seq = seq
                    self.post_mut4_gRNA_seq_phases = phases
                    #mutate PAM if in 5 UTR/intron
                    if (phases[-1:] == "5" or phases[-2:-1] == "5"):
                        # mutate PAM if it's in 5 UTR (phase == 5)
                        if phases[-1:] == "5": #last position phase = 5 (the second G in NGG)
                            self.post_mut4_gRNA_seq = seq[:-1] + "c"
                        if phases[-2:-1] == "5": #second position phase = 5 (the first G in NGG)
                            self.post_mut4_gRNA_seq = self.post_mut4_gRNA_seq[:-2] + "c" + self.post_mut4_gRNA_seq[-1:]
                        self.cdf_score_post_mut4 = cfd_score(self.gRNA_seq, self.post_mut4_gRNA_seq)

                        #put disrupted-and-mutated seq back to arms
                        ssODN = f"{left}{self.tag}{right}"
                        ssODN = ssODN.replace(seq,self.post_mut4_gRNA_seq)
                        ssODN = ssODN.replace(self.revcom(seq),self.revcom(self.post_mut4_gRNA_seq))
                        self.left_flk_seq_CodonMut4 = ssODN[0:len(self.left_flk_seq)]
                        self.right_flk_seq_CodonMut4 =ssODN[-len(self.left_flk_seq):]
                        print(self.post_mut4_gRNA_seq)

                #mutate protospacer
                left, right, null, seq, phases = self.get_uptodate_mut() #get up-to-date gRNA seq and phases
                self.post_mut4_gRNA_seq = seq
                self.post_mut4_gRNA_seq_phases = phases
                latest_cfd = cfd
                if hasattr(self,"cdf_score_post_mut4"):
                    latest_cfd = self.cdf_score_post_mut4
                counter = -1
                if latest_cfd>self.cfdThres:
                    for idx,item in list(enumerate(phases)): # go through protospacer prior to PAM
                        if (idx == (len(phases) - 1) or idx == (len(phases) - 2)): #skip PAM (in this protospacer block of code)
                            continue
                        if item == "5": # 0=3'UTR, 5=5'UTR (minus sites that are 3-bp dist to exon-intro junctions) #TODO check remove of item=="0"
                            counter += 1
                            if counter % 3 != 0:
                                continue #  mutate 1 in every 3 bp
                            base = seq[idx]
                            mutbase = self.single_base_muation(base)
                            self.post_mut4_gRNA_seq = self.post_mut4_gRNA_seq[:idx] + mutbase + self.post_mut4_gRNA_seq[idx+1:]
                            self.cdf_score_post_mut4 = cfd_score(self.gRNA_seq, self.post_mut4_gRNA_seq)
                            #early stop if CFD goes below self.cfdThres
                            if self.cdf_score_post_mut4<self.cfdThres:
                                break
                    #put mutated seq back to arms
                    ssODN = f"{left}{self.tag}{right}"
                    ssODN = ssODN.replace(seq,self.post_mut4_gRNA_seq)
                    ssODN = ssODN.replace(self.revcom(seq),self.revcom(self.post_mut4_gRNA_seq))
                    self.left_flk_seq_CodonMut4 = ssODN[0:len(self.left_flk_seq)]
                    self.right_flk_seq_CodonMut4 =ssODN[-len(self.left_flk_seq):]
                    print(self.post_mut4_gRNA_seq)

                # mutate PAM afterwards (if skipped earlier)
                latest_cfd = cfd
                if hasattr(self,"cdf_score_post_mut4"):
                    latest_cfd = self.cdf_score_post_mut4
                if latest_cfd > self.cfdThres and self.recode_order == "protospacer_first":
                    left, right, cdf, seq, phases = self.get_uptodate_mut()  # get up-to-date gRNA seq and phases
                    self.post_mut4_gRNA_seq = seq
                    self.post_mut4_gRNA_seq_phases = phases
                    #mutate PAM if in 5 UTR/intron
                    if (phases[-1:] == "5" or phases[-2:-1] == "5"):
                        # mutate PAM if it's in 5 UTR (phase == 5)
                        if phases[-1:] == "5": #last position phase = 5 (the second G in NGG)
                            self.post_mut4_gRNA_seq = seq[:-1] + "c"
                        if phases[-2:-1] == "5": #second position phase = 5 (the first G in NGG)
                            self.post_mut4_gRNA_seq = self.post_mut4_gRNA_seq[:-2] + "c" + self.post_mut4_gRNA_seq[-1:]
                        self.cdf_score_post_mut4 = cfd_score(self.gRNA_seq, self.post_mut4_gRNA_seq)

                        #put disrupted-and-mutated seq back to arms
                        ssODN = f"{left}{self.tag}{right}"
                        ssODN = ssODN.replace(seq,self.post_mut4_gRNA_seq)
                        ssODN = ssODN.replace(self.revcom(seq),self.revcom(self.post_mut4_gRNA_seq))
                        self.left_flk_seq_CodonMut4 = ssODN[0:len(self.left_flk_seq)]
                        self.right_flk_seq_CodonMut4 =ssODN[-len(self.left_flk_seq):]
                        print(self.post_mut4_gRNA_seq)

                if hasattr(self,"cdf_score_post_mut4"):
                    self.info_phase4_5UTR=[cfd, self.cdf_score_post_mut4]

            left,right,cfd,seq,phases = self.get_uptodate_mut()
            self.postphase4ODN = left + tag + right

            #########
            #phase 5#
            #########
            left,right,cfd,seq,phases = self.get_uptodate_mut()
            self.Donor_phases = self.left_flk_phases + "X"*len(self.tag) + self.right_flk_phases
            self.Donor_postMut = left + self.tag + right
            self.Donor_prephase5 = self.Donor_postMut

            #slide windows scan for new cut sites
            self.info_p5 = ""
            #print(f"{self.Donor_postMut}\n{self.Donor_phases}")

            #scan for maximum chimeric cfd (no recoding yet)
            fwd_scan_highest_cfd_no_recode = self.slide_win_cfd_coding(self.Donor_postMut, self.Donor_phases) #this is the chimeric cfd if no recoding
            rev_scan_highest_cfd_no_recode = self.slide_win_cfd_noncoding(str(Seq(self.Donor_postMut).reverse_complement()), self.Donor_phases[::-1])
            self.scan_highest_cfd_no_recode = max([fwd_scan_highest_cfd_no_recode, rev_scan_highest_cfd_no_recode])

            #scan and recode
            self.phase5recoded = False
            fwd_scan_highest_cfd = self.slide_win_mutation_coding(self.Donor_postMut, self.Donor_phases) # this is the chimeric cfd after recoding
            #scan and recode revcom
            #print(f"{str(Seq(self.Donor_postMut).reverse_complement())}\n{self.Donor_phases[::-1]}")
            rev_scan_highest_cfd = self.slide_win_mutation_noncoding(str(Seq(self.Donor_postMut).reverse_complement()), self.Donor_phases[::-1]) # this function modifies self.Donor_postMut


            #get the highest cfd among all possible cutsites
            scan_highest_cfd = max([fwd_scan_highest_cfd,rev_scan_highest_cfd])

            self.cdf_score_highest_in_win_scan = scan_highest_cfd

            ############################
            # log recoding information #
            ############################
            self.info_p1 = ""
            if not self.recoding_args["recoding_stop_recut_only"]:
                self.info_p1 = "".join(
                    f"--------------------phase 1 mutate seq between cut to insert----------------------------------------------------------------------\n"
                    f"phase1.cut-to-insert + 2bp padding on both sides (extend to full codons)\n"
                    f"phase1.seq                             :{self.ins2cut.seq}\n"
                    f"phase1.Phases                          :{self.ins2cut.phases}\n"
                    f"phase1.Coordinates                     :{self.ins2cut.start}-{self.ins2cut.end}\n"
                    f"phase1.seq         (trimmed into frame):{self.ins2cut_LRtrimed.seq}\n"
                    f"phase1.Phases      (trimmed into frame):{self.ins2cut_LRtrimed.phases}\n"
                    f"phase1.Coordinates (trimmed into frame):{self.ins2cut_LRtrimed.start}-{self.ins2cut_LRtrimed.end}\n"
                    f"phase1.mutated seq (trimmed into frame):{mutated_subseq}\n"
                    f"--> display gRNA and check disruption <--\n"
                    f"phase1.                      gRNA:{self.gRNA_seq}\n"
                    f"phase1.                    Phases:{self.gRNA_seq_phases}\n"
                    f"phase1.after mutation and payload:{self.post_mut_ins_gRNA_seq}\n"
                    f"phase1.                    Phases:{self.post_mut_ins_gRNA_seq_phases}\n"
                    f"phase1.                       CFD:{self.cdf_score_post_mut_ins:.4f}\n"
                    f"--> display mutation in HDR arms <--\n"
                    f"phase1.original arms   :{self.gRNA_lc_Larm}||{self.gRNA_lc_Rarm}\n"
                    f"phase1.cut2insert mut  :{self.left_flk_seq_CodonMut}||{self.right_flk_seq_CodonMut}\n"
                    f"phase1.arms with tag   :{self.left_flk_seq_CodonMut}|{self.tag}|{self.right_flk_seq_CodonMut}\n")

            if hasattr(self,"cdf_score_post_mut2"):
                self.info_p2 = "".join(
                     f"--------------------phase 2: if cfd > self.cfdThres, mutating PAM and protospacer if in 3UTR-----------------------------------------------\n"
                     f"phase 2.                          gRNA:{self.gRNA_seq}\n"
                     f"phase 2.                        Phases:{self.gRNA_seq_phases}\n"
                     f"phase 2.             gRNA post phase 2:{self.post_mut2_gRNA_seq}\n"
                     f"phase 2.                        Phases:{self.post_mut2_gRNA_seq_phases}\n"
                     f"phase 2.                           CFD:{self.cdf_score_post_mut2:.4f}\n"
                     f"phase 2.       post phase 2 ssODN:{self.postphase2ODN}\n"
                )
            else:
                self.info_p2 = \
                     f"--------------------phase 2 skipped-----------------------------------------------------------------------------------------------\n"
            if hasattr(self,"cdf_score_post_mut3"):
                self.info_p3 = "".join(
                     f"--------------------phase 3: if cfd > self.cfdThres, silently mutate gRNA seq not covered between insert and cut---------------------------\n"
                     f"phase 3.                     gRNA:{self.gRNA_seq}\n"
                     f"phase 3.                   Phases:{self.gRNA_seq_phases}\n"
                     f"phase 3.        gRNA post phase 3:{self.post_mut3_gRNA_seq}\n"
                     f"phase 3.                   Phases:{self.post_mut3_gRNA_seq_phases}\n"
                     f"phase 3.                      CFD:{self.cdf_score_post_mut3:.4f}\n"
                     f"phase 3.       post phase 3 ssODN:{self.postphase3ODN}\n"
                    )
            else:
                self.info_p3 = \
                     f"--------------------phase 3 skipped-----------------------------------------------------------------------------------------------\n"
            if hasattr(self,"cdf_score_post_mut4"):
                self.info_p4 = "".join(
                     f"--------------------phase 4: if cfd > self.cfdThres, mutating gRNA if in 5UTR--------------------------------------------------------------\n"
                     f"phase 4.                     gRNA:{self.gRNA_seq}\n"
                     f"phase 4.                   Phases:{self.gRNA_seq_phases}\n"
                     f"phase 4.        gRNA post phase 4:{self.post_mut4_gRNA_seq}\n"
                     f"phase 4.                   Phases:{self.post_mut4_gRNA_seq_phases}\n"
                     f"phase 4.                      CFD:{self.cdf_score_post_mut4:.4f}\n"
                     f"phase 4.       post phase 4 ssODN:{self.postphase4ODN}\n"
                )
            else:
                self.info_p4 = \
                    f"--------------------phase 4 skipped-----------------------------------------------------------------------------------------------\n"
            if self.phase5recoded == True:
                self.info_p5 = "".join(
                    f"--------------------phase 5: sliding window check of recutting--------------------------------------------------------------------\n"
                    f"phase 5.  ssODN pre-phase 5:{self.Donor_prephase5}\n"
                    f"phase 5.ssODN after-phase 5:{self.Donor_postMut}\n"
                    f"phase 5.             phases:{self.Donor_phases}\n"
                    f"phase 5. maximum cfd in window scan analysis(after recoding): {self.cdf_score_highest_in_win_scan:.6f}\n") + self.info_p5
            else:
                self.info_p5 = "".join(
                    f"--------------------phase 5: sliding window check of recutting--------------------------------------------------------------------\n"
                    f"phase 5.  no recoding was performed b/c maximum cfd in window scan analysis is {self.cdf_score_highest_in_win_scan:.6f} \n") + self.info_p5
        #########################
        # entirely skip recoding#
        #########################
        else:
            #scan to check the highest cfd
            self.gRNA_seq, Null, self.gRNA_seq_phases, Null = self.get_post_integration_gRNA(self.left_flk_seq,self.right_flk_seq) #get the original gRNA
            ODN = f"{self.gRNA_lc_Larm}{self.tag}{self.gRNA_lc_Rarm}"
            Donor_phases = self.gRNA_lc_Larm + "X"*len(self.tag) + self.gRNA_lc_Rarm

            fwd_scan_highest_cfd = self.slide_win_cfd_coding(ODN, Donor_phases) # this function modifies self.Donor_postMut
            #scan the revcom
            rev_scan_highest_cfd = self.slide_win_cfd_noncoding(str(Seq(ODN).reverse_complement()), Donor_phases[::-1]) # this function modifies self.Donor_postMut

            #get the highest cfd among all possible cutsites
            scan_highest_cfd = max([fwd_scan_highest_cfd,rev_scan_highest_cfd])

            self.cdf_score_post_ins = scan_highest_cfd
            self.Donor_postMut = "recoding turned off"

        ##############################
        #finished or skipped recoding#
        ##############################
        left,right,cfd,seq,phases = self.get_uptodate_mut() # not including slide window scan and mutation
        self.Donor_vanillia = f"{self.gRNA_lc_Larm}{self.tag}{self.gRNA_lc_Rarm}"
        if self.Donor_postMut == "recoding turned off":
            self.final_cfd = self.cdf_score_post_ins
        else:
            self.final_cfd = max(cfd,scan_highest_cfd) #this should be the highest cfd from all phases,  cfd= phase 1-4, scan_highest_cfd = phase5

        ################
        #dsDNA donor   #
        ################
        if self.Donor_type == "dsDNA":
            self.Donor_final = self.Donor_vanillia
            #get recoded donor
            if self.Donor_postMut != "recoding turned off":
                self.Donor_final = self.Donor_postMut
            #set effective HA length
            self.effective_HA_len="N/A for dsDNA"

            ################################
            #check synthesis considerations#
            ################################
            #restriction cuts
            #print(Restriction.BsaI.site)
            enzyme_list = self.enzymes2check.split("|")
            if len(enzyme_list)>=1:
                for enzyme in enzyme_list:
                    if hasattr(Restriction, enzyme):
                        RE_object = getattr(Restriction, enzyme)
                        enzyme_cutsites = [str(pos) for pos in RE_object.search(Seq(self.Donor_final))] # will find cutsites on both strands
                        if len(enzyme_cutsites) > 0:
                            enzyme_cutPos = ";".join(enzyme_cutsites)
                            self.synFlags.append(f"Cut by {enzyme} @{enzyme_cutPos}")

            #custom sequences to avoid
            if self.CustomSeq2Avoid!="":
                seqs2avoid = self.CustomSeq2Avoid.split("|")
                if len(seqs2avoid)>=1:
                    for seq in seqs2avoid:
                        locations = []
                        for m in re.finditer(seq, self.Donor_final,flags=re.IGNORECASE):
                            locations.append(str(m.start()+1))
                        for m in re.finditer(seq, self.revcom(self.Donor_final),flags=re.IGNORECASE):
                            start = m.start() + 1
                            locations.append(str(len(self.Donor_final) - start + 1))
                        if len(locations) > 0: # current seq found
                            locs = "@".join(locations)
                            self.synFlags.append(f"{seq}@{locs}")

            #GC content
            seq_noAmbiguous = re.sub(r'[^ATCGatcg]', '', self.Donor_final)
            global_GC = GC(seq_noAmbiguous)
            #print(global_GC)
            if global_GC < 25:
                self.synFlags.append(f"global GC content {global_GC:.2f}% < 25%")
            if global_GC > 65:
                self.synFlags.append(f"global GC content {global_GC:.2f}% > 65%")

            #GC content skew (slide windown analysis
            win_GC = self.slide_win_GC_content(seq=seq_noAmbiguous, win_size=50)
            #print(win_GC)
            max_diff = max(win_GC) - min(win_GC)
            #print(max_diff)
            if max_diff > 52:
                self.synFlags.append(f"Max difference of slide window GC content {max_diff:.2f}% > 52%")

            #homopolyer
            hp_res = [(m.group(), m.start()+1) for m in re.finditer(r'([ACGT])\1{9,}', seq_noAmbiguous.upper())]
            if len(hp_res)>0:
                hp_res_display = [f"({t[0]}@{t[1]})" for t in hp_res]
                hp_res_display = "".join(hp_res_display)
                self.synFlags.append(f"Homopolymer > 10bp (sequence@start): {hp_res_display}")

            if len(self.synFlags) == 0:
                self.synFlags = "None"
            else:
                self.synFlags = "; ".join(self.synFlags)

        ################
        #ssDNA donor   #
        ################
        if self.Donor_type == "ssDNA":
            self.Donor_final = self.Donor_vanillia
            #get recoded donor
            if not self.recoding_args["recoding_off"]: #recoding is on
                self.Donor_final = self.Donor_postMut
            self.synFlags = "N/A for ssDNA"
            self.effective_HA_len = "N/A if not enforcing max donor length"
            ###################################################
            #Enforce max payload size and centering           #
            #This handles both recoded and non-recoded donor  #
            ###################################################
            if self.ssDNA_max_size is not None:
                self.effective_HA_len = len(self.gRNA_lc_Larm) #initizalize effective HA length with the maximum value
                #print(f"Centering")
                #print(f"Donor_vanillia {self.gRNA_lc_Larm}{self.tag}{self.gRNA_lc_Rarm}\n"
                #      f"Donor_postmut  {self.Donor_postMut}")
                if not self.recoding_args["recoding_off"]: #recoding is on
                    diff_loc = self.get_diff_loc(str1 = f"{self.gRNA_lc_Larm}{self.tag}{self.gRNA_lc_Rarm}", str2 =f"{self.Donor_postMut}")
                else:
                    diff_loc = [] # skip get_diff_loc if recoding is off
                #print(f"{diff_loc}")
                #print(f"lengths:{len(self.gRNA_lc_Larm)}|{len(self.tag)}|{len(self.gRNA_lc_Rarm)}")

                if len(diff_loc) == 0: #no recoding (including recoding turned off and no-recoding with recoding turned on)
                    _HA_len = (self.ssDNA_max_size - len(self.tag))/2
                    start = len(self.gRNA_lc_Larm) - _HA_len - 1
                    end = start + _HA_len + len (self.tag) + _HA_len
                    #print(f"no recoding, start={start}\tend={end}")
                else: #with recoding
                    recoding_left = min(diff_loc) # 0-indexed
                    recoding_right = max(diff_loc) # 0-indexed
                    tag_start = len(self.gRNA_lc_Larm) # 0-indexed
                    tag_end = len(self.gRNA_lc_Larm) + len(tag) - 1 # 0-indexed
                    centerpiece_start = min([recoding_left,recoding_right,tag_start,tag_end]) #center piece is the payload + recoded region
                    centerpiece_end = max([recoding_left,recoding_right,tag_start,tag_end])
                    centerpiece_len = centerpiece_end - centerpiece_start + 1
                    _HA_len = math.floor((self.ssDNA_max_size - centerpiece_len)/2)
                    start = centerpiece_start - _HA_len
                    end = centerpiece_end + _HA_len + 1 # need to get the base at position:end
                    _len= end - start
                    #print(f"recoding, start={start}\tend={end}\t centerpiece:{centerpiece_start}-{centerpiece_end} len={_len} HA_len={_HA_len}")
                if hasattr(self,"Donor_postMut") and not self.recoding_args["recoding_off"]:  # donor is recoded
                    self.Donor_final = self.Donor_postMut[int(start):int(end)]
                else: # donor is not recoded
                    #print(f"{start}-{end}")
                    self.Donor_final = self.Donor_final[int(start):int(end)]
                self.effective_HA_len = _HA_len
                #print(f"effective_HA_len {_HA_len}")
            ###################
            #strand selection #
            ###################
            #check PAM-less cutting
            fwd_scan_highest_PAMless_cfd = self.slide_win_PAMLESScfd_coding(self.Donor_final) # this function modifies self.Donor_postMut
            rev_scan_highest_PAMless_cfd = self.slide_win_PAMLESScfd_noncoding(str(Seq(self.Donor_final).reverse_complement())) # this function modifies self.Donor_postMut
            PAMless_cfd = max([fwd_scan_highest_PAMless_cfd,rev_scan_highest_PAMless_cfd])
            #print(f"PAMless_cfd:{PAMless_cfd}")

            if PAMless_cfd > self.cfdThres: # PAM-independent cutting can happen, choose gRNA strand
                if not self.ENST_strand * self.gStrand:
                    self.Donor_final = str(Seq(self.Donor_final).reverse_complement()) # take revcom if gRNA is not on the same strand as the ENST (coding)
            else:   # PAM-independent cutting canNOT happen, choose Manu strand
                self.Donor_final = self.select_Manu_strand(self.Donor_final)


    #############
    #END OF INIT#
    #############
    #TODO scan recoded sequence for PAM-less recut
    def slide_win_PAMLESScfd_noncoding(self,seq):
        _arm_len = int((len(seq)-len(self.tag))/2)
        highest_cfd = 0
        #get sliding window as an iterator
        it = self.sliding_window(seq,23)
        #go through sliding windows #NOTE: Donor_postMut is always in the coding straind
        n_window = 0
        for n_mer in it:
            #check for "N"s in the sequence
            n_mer = "".join(n_mer)
            if "N" in n_mer or "n" in n_mer:
                n_window+=1
                continue
            #replace the last 3-mer with the PAM from gRNA
            n_mer_addPAM = n_mer[:-3] + self.gRNA_seq[-3:]
            #print(f"{self.gRNA_seq}\n{n_mer_addPAM}\n")
            cfd = cfd_score(self.gRNA_seq, n_mer_addPAM)
            #print(f"DIAG cfd_noncoding: {cfd:.6f}")
            if cfd>=highest_cfd:
                highest_cfd = cfd
            n_window+=1
        return highest_cfd

    def slide_win_PAMLESScfd_coding(self,seq):
        highest_cfd = 0
         #get sliding window as an iterator
        it = self.sliding_window(seq,23)
        #go through sliding windows #NOTE: Donor_postMut is always in the coding straind
        n_window = 0
        for n_mer in it:
            #check for "N"s in the sequence
            n_mer = "".join(n_mer)
            if "N" in n_mer or "n" in n_mer:
                n_window+=1
                continue
            #replace the last 3-mer with the PAM from gRNA
            n_mer_addPAM = n_mer[:-3] + self.gRNA_seq[-3:]
            #print(f"{self.gRNA_seq}\n{n_mer_addPAM}\n")
            cfd = cfd_score(self.gRNA_seq, n_mer_addPAM)
            #print(f"DIAG cfd_coding: {cfd:.6f}")
            if cfd>=highest_cfd:
                highest_cfd = cfd
            n_window+=1
        return highest_cfd

    def slide_win_GC_content(self,seq,win_size):
        '''
        returns a list of GC contents of sliding windows
        '''
        GC_list=[]
        it = self.sliding_window(seq,win_size)
        for n_mer in it:
            n_mer = "".join(n_mer)
            #print(n_mer)
            GC_list.append(GC(n_mer))
        return GC_list

    def get_diff_loc(self, str1,str2):
        """
        get the locations where two strings differ
        """
        assert len(str1) == len(str2)
        diff_loc = []
        for idx,val in enumerate(str1):
            char1 = val.upper()
            char2 = str2[idx].upper()
            if char1 != char2:
                diff_loc.append(idx)
        return(diff_loc)

    def revcom(self,seq):
        return str(Seq(seq).reverse_complement())

    def slide_win_mutation_noncoding(self,seq,phases):
        '''
        for noncoding strand
        use with untrimmed donor
        '''
        highest_cfd = 0
        #get sliding window as an iterator
        it = self.sliding_window(seq,23)
        #go through sliding windows #NOTE: Donor_postMut is always in the coding straind
        n_window = 0
        for n_mer in it:
            #skip the homology arms (while not skipping the payload)
            if 0 <= n_window <= (len(self.left_flk_seq) - 23 - 1):
                #print(f"skipping window: {n_window}")
                n_window+=1
                continue
            if n_window >= len(self.left_flk_seq) + len(self.tag) - 1:
                n_window+=1
                #print(f"skipping window: {n_window}")
                continue
            #check for "N"s in the sequence
            n_mer = "".join(n_mer)
            if "N" in n_mer or "n" in n_mer:
                n_window+=1
                continue
            cfd = cfd_score(self.gRNA_seq, n_mer)
            n_mer_phases = phases[n_window: n_window + 23]
            o = chimeric_gRNA(gRNA_seq=self.gRNA_seq, gRNA_phases="",
                              recut_seq=n_mer, recut_cfd = cfd, recut_phases=n_mer_phases)
            if cfd>self.cfdThres:
                self.info_p5 = self.info_p5 + "".join(f"win-     original gRNA:{self.gRNA_seq}\n"
                                                      f"         chimeric site:{o.recut.seq}\tcfd: {o.recut.cfd}\n"
                                                      f"       chimeric phases:{o.recut.phases}\n")
                # print(f"win-     original gRNA:{self.gRNA_seq}\n"
                #       f"         chimeric site:{o.recut.seq}\tcfd: {o.recut.cfd}\n"
                #       f"       chimeric phases:{o.recut.phases}")
                #start to mutate (it is sufficient to only work with the coding strand
                #try silently mutate **the coding strand**
                seq_obj = seq_w_phase(seq = self.revcom(o.recut.seq), phases = o.recut.phases[::-1], start = 1, end = 23)
                seq_obj_ltrim = self.trim_left_into_frame(seq_obj)
                seq_obj_lrtrim = self.trim_right_into_frame(seq_obj_ltrim)

                if len(seq_obj_lrtrim.seq) >= 3:
                    #print(f"mutating:\n{seq_obj_lrtrim.seq}\n{seq_obj_lrtrim.phases}")
                    mutated = self.get_silent_mutations(seq_obj_lrtrim.seq.upper())
                    self.phase5recoded = True
                else:
                    mutated = seq_obj_lrtrim.seq
                untrimmed = seq_obj.seq.replace(seq_obj_lrtrim.seq, mutated) # untrim: replace trimmed part with the mutated part
                o.update(newSeq=self.revcom(untrimmed)) #update cfd, rc seq , rc cfd according to mutated seq
                self.info_p5 = self.info_p5 + "".join(
                      f"syn          before mut:{n_mer}\n"
                      f"              after mut:{o.recut.seq}\t cfd:{o.recut.cfd}\n"
                      f"                 phases:{n_mer_phases}\n")
                # print(f"syn          before mut:{n_mer}\n"
                #       f"              after mut:{o.recut.seq}\t cfd:{o.recut.cfd}\n"
                #       f"                 phases:{n_mer_phases}")
            #try mutate in 3UTR
            if o.recut.cfd>self.cfdThres:
                # mutate PAM if it's in 3' UTR or intron (phase == 0)
                if n_mer_phases[-1:] == "0": #last position phase = 0 (the second G in NGG)
                    o.update(newSeq=o.recut.seq[:-1] + "c")
                    self.phase5recoded = True
                if n_mer_phases[-2:-1] == "0": #second position phase = 0 (the first G in NGG)
                    o.update(newSeq=o.recut.seq[:-2] + "c" + o.recut.seq[-1:])
                    self.phase5recoded = True

                if o.recut.cfd>self.cfdThres:
                    for idx,item in reversed(list(enumerate(n_mer_phases))):
                        if idx == (len(n_mer_phases) - 1) or idx == (len(n_mer_phases) - 2):
                            continue #skip PAM
                        if item == "0": # "0" means 3'UTR (5'UTR is labeled "5")
                            self.phase5recoded = True
                            base = o.recut.seq[idx]
                            mutbase = self.single_base_muation(base)
                            o.update(newSeq=o.recut.seq[:idx] + mutbase + o.recut.seq[idx+1:])
                            #early stop if CFD goes below self.cfdThres
                            if o.recut.cfd<self.cfdThres:
                                break
                self.info_p5 = self.info_p5 + "".join(f"3UTR/intron  before mut:{n_mer}\n"
                                                      f"              after mut:{o.recut.seq}\t cfd:{o.recut.cfd}\n"
                                                      f"                 phases:{n_mer_phases}\n")
                # print(f"3UTR/intron  before mut:{n_mer}\n"
                #       f"              after mut:{o.recut.seq}\t cfd:{o.recut.cfd}\n"
                #       f"                 phases:{n_mer_phases}")
            #try mutate in 5UTR
            if o.recut.cfd>self.cfdThres:
                self.info_phase5_5UTR=[o.recut.cfd, ""]
                # mutate PAM if it's in 5' UTR (phase == 5)
                if n_mer_phases[-1:] == "5": #last position phase = 5 (the second G in NGG)
                    o.update(newSeq=o.recut.seq[:-1] + "c")
                    self.phase5recoded = True
                if n_mer_phases[-2:-1] == "5": #second position phase = 5 (the first G in NGG)
                    o.update(newSeq=o.recut.seq[:-2] + "c" + o.recut.seq[-1:])
                    self.phase5recoded = True

                self.info_p5 = self.info_p5 + "".join(f"5UTR         before mut:{n_mer}\n"
                                      f"              after mut:{o.recut.seq}\t cfd:{o.recut.cfd}\n"
                                      f"                 phases:{n_mer_phases}\n")
                if o.recut.cfd>self.cfdThres:
                    #mutate protospacer if in 5UTR
                    for idx,item in reversed(list(enumerate(n_mer_phases))):
                        if idx == (len(n_mer_phases) - 1) or idx == (len(n_mer_phases) - 2):
                            continue #skip PAM
                        if item == "5": # 5'UTR is labeled "5"
                            self.phase5recoded = True
                            base = o.recut.seq[idx]
                            mutbase = self.single_base_muation(base)
                            o.update(newSeq=o.recut.seq[:idx] + mutbase + o.recut.seq[idx+1:])
                            #early stop if CFD goes below self.cfdThres
                            if o.recut.cfd<self.cfdThres:
                                break
                    self.info_p5 = self.info_p5 + "".join(f"5UTR         before mut:{n_mer}\n"
                                                          f"              after mut:{o.recut.seq}\t cfd:{o.recut.cfd}\n"
                                                          f"                 phases:{n_mer_phases}\n")
                    # print(f"5UTR         before mut:{n_mer}\n"
                    #       f"              after mut:{o.recut.seq}\t cfd:{o.recut.cfd}\n"
                    #       f"                 phases:{n_mer_phases}")

                self.info_phase5_5UTR[1] = o.recut.cfd
            #print(f"DIAG mut_noncoding: {o.recut.cfd:.6f}")
            if o.recut.seq != n_mer and o.recut.cfd <= cfd: # mutations been made and it decreases cfd
                self.info_p5 = self.info_p5 + "".join(f"phase 5 final\n"
                                                      f"  orig gRNA: {o.gRNA.seq}\n"
                                                      f"pre-mut seq: {n_mer} cfd:{cfd}\n"
                                                      f"   post-mut: {o.recut.seq} cfd:{o.recut.cfd}\n")
                # print(f"phase 5 final for current window\n"
                #       f"  orig gRNA: {o.gRNA.seq}\n"
                #       f"pre-mut seq: {n_mer} cfd:{cfd}\n"
                #       f"   post-mut: {o.recut.seq} cfd:{o.recut.cfd}")
                #put mut seq back into ssODN, need revcom here b/c the input/seq is revcomed from self.Donor_postMut
                self.Donor_postMut = self.revcom(seq[0:n_window] + o.recut.seq + seq[n_window + 23:])
                highest_cfd = max([highest_cfd, o.recut.cfd]) #update highest cfd
            else:
                highest_cfd = max([highest_cfd, cfd]) #update highest cfd
            n_window+=1
        return highest_cfd

    def slide_win_cfd_noncoding(self,seq,phases):
        '''
        scan the donor and return the highest cfd
        '''
        highest_cfd = 0
        #get sliding window as an iterator
        it = self.sliding_window(seq,23)
        #go through sliding windows #NOTE: Donor_postMut is always in the coding straind
        n_window = 0
        for n_mer in it:
            # #skip non-chimeric part of the homology arm #!!!=> we should not skip non-chimeric parts b/c the payload may contain cutsites
            # if 0 <= n_window <= (len(self.left_flk_seq) - 23 - 1):
            #     #print(f"skipping window: {n_window}")
            #     n_window+=1
            #     continue
            # if n_window >= len(self.left_flk_seq) + len(self.tag) - 1:
            #     n_window+=1
            #     #print(f"skipping window: {n_window}")
            #    continue
            #check for "N"s in the sequence
            n_mer = "".join(n_mer)
            if "N" in n_mer or "n" in n_mer:
                n_window+=1
                continue
            cfd = cfd_score(self.gRNA_seq, n_mer)
            #print(f"DIAG cfd_noncoding: {cfd:.6f}")
            if cfd>=highest_cfd:
                highest_cfd = cfd
            n_window+=1
        return highest_cfd

    def slide_win_mutation_coding(self,seq,phases):
        '''
        for coding strand
        '''
        highest_cfd = 0
         #get sliding window as an iterator
        it = self.sliding_window(seq,23)
        #go through sliding windows #NOTE: Donor_postMut is always in the coding straind
        n_window = 0
        for n_mer in it:
            #skip the homology arms (while not skipping the payload)
            if 0 <= n_window <= (len(self.left_flk_seq) - 23 - 1):
                #print(f"skipping window: {n_window}")
                n_window+=1
                continue
            if n_window >= len(self.left_flk_seq) + len(self.tag) - 1:
                n_window+=1
                #print(f"skipping window: {n_window}")
                continue
            #check for "N"s in the sequence
            n_mer = "".join(n_mer)
            if "N" in n_mer or "n" in n_mer:
                n_window+=1
                continue
            cfd = cfd_score(self.gRNA_seq, n_mer)
            n_mer_phases = phases[n_window: n_window + 23]
            o = chimeric_gRNA(gRNA_seq=self.gRNA_seq, gRNA_phases="",
                              recut_seq=n_mer, recut_cfd = cfd, recut_phases=n_mer_phases)
            if cfd>self.cfdThres:
                self.info_p5 = self.info_p5 + "".join(f"win+     original gRNA:{self.gRNA_seq}\n"
                                      f"         chimeric site:{o.recut.seq}\tcfd: {o.recut.cfd}\n"
                                      f"       chimeric phases:{o.recut.phases}\n")
                # print(f"win+     original gRNA:{self.gRNA_seq}\n"
                #       f"         chimeric site:{o.recut.seq}\tcfd: {o.recut.cfd}\n"
                #       f"       chimeric phases:{o.recut.phases}")
                #start to mutate (it is sufficient to only work with the coding strand
                #try silently mutate the coding strand
                seq_obj = seq_w_phase(seq = o.recut.seq, phases = o.recut.phases, start = 1, end = 23)
                seq_obj_ltrim = self.trim_left_into_frame(seq_obj)
                seq_obj_lrtrim = self.trim_right_into_frame(seq_obj_ltrim)

                if len(seq_obj_lrtrim.seq)>=3:
                    #print(f"mutating:\n{seq_obj_lrtrim.seq}\n{seq_obj_lrtrim.phases}")
                    mutated = self.get_silent_mutations(seq_obj_lrtrim.seq.upper())
                    self.phase5recoded = True
                else:
                    mutated = seq_obj_lrtrim.seq
                o.update(newSeq=o.recut.seq.replace(seq_obj_lrtrim.seq, mutated)) #update cfd, rc seq , rc cfd according to mutated seq
                self.info_p5 = self.info_p5 + "".join(
                      f"syn          before mut:{n_mer}\n"
                      f"              after mut:{o.recut.seq}\t cfd:{o.recut.cfd}\n"
                      f"                 phases:{n_mer_phases}\n")
                # print(f"syn          before mut:{n_mer}\n"
                #       f"              after mut:{o.recut.seq}\t cfd:{o.recut.cfd}\n"
                #       f"                 phases:{n_mer_phases}")
            #try mutate in 3UTR/intron
            if o.recut.cfd>self.cfdThres:
                # mutate PAM if it's in 3' UTR (phase == 0)
                if n_mer_phases[-1:] == "0": #last position phase = 0 (the second G in NGG)
                    o.update(newSeq=o.recut.seq[:-1] + "c")
                    self.phase5recoded = True
                if n_mer_phases[-2:-1] == "0": #second position phase = 0 (the first G in NGG)
                    o.update(newSeq=o.recut.seq[:-2] + "c" + o.recut.seq[-1:])
                    self.phase5recoded = True
                if o.recut.cfd>self.cfdThres:
                    #mutate protospacer if in 3UTR/intron
                    for idx,item in reversed(list(enumerate(n_mer_phases))):
                        if idx == (len(n_mer_phases) - 1) or idx == (len(n_mer_phases) - 2):
                            continue #skip PAM
                        if item == "0": # "0" means 3'UTR (5'UTR is labeled "5")
                            self.phase5recoded = True
                            base = o.recut.seq[idx]
                            mutbase = self.single_base_muation(base)
                            o.update(newSeq=o.recut.seq[:idx] + mutbase + o.recut.seq[idx+1:])
                            #early stop if CFD goes below self.cfdThres
                            if o.recut.cfd<self.cfdThres:
                                break
                self.info_p5 = self.info_p5 + "".join(f"3UTR/intron  before mut:{n_mer}\n"
                                                      f"              after mut:{o.recut.seq}\t cfd:{o.recut.cfd}\n"
                                                      f"                 phases:{n_mer_phases}\n")
                # print(f"3UTR/intron  before mut:{n_mer}\n"
                #       f"              after mut:{o.recut.seq}\t cfd:{o.recut.cfd}\n"
                #       f"                 phases:{n_mer_phases}")
            #try mutate in 5UTR
            if o.recut.cfd>self.cfdThres:
                self.info_phase5_5UTR=[o.recut.cfd, ""]
                # mutate PAM if it's in 5' UTR (phase == 5)
                if n_mer_phases[-1:] == "5": #last position phase = 5 (the second G in NGG)
                    o.update(newSeq=o.recut.seq[:-1] + "c")
                    self.phase5recoded = True
                if n_mer_phases[-2:-1] == "5": #second position phase = 5 (the first G in NGG)
                    o.update(newSeq=o.recut.seq[:-2] + "c" + o.recut.seq[-1:])
                    self.phase5recoded = True

                self.info_p5 = self.info_p5 + "".join(f"5UTR         before mut:{n_mer}\n"
                                                      f"              after mut:{o.recut.seq}\t cfd:{o.recut.cfd}\n"
                                                      f"                 phases:{n_mer_phases}\n")
                if o.recut.cfd>self.cfdThres:
                    #mutate protospacer if in 5UTR
                    for idx,item in reversed(list(enumerate(n_mer_phases))):
                        if idx == (len(n_mer_phases) - 1) or idx == (len(n_mer_phases) - 2):
                            continue #skip PAM
                        if item == "5": # 5'UTR is labeled "5"
                            self.phase5recoded = True
                            base = o.recut.seq[idx]
                            mutbase = self.single_base_muation(base)
                            o.update(newSeq=o.recut.seq[:idx] + mutbase + o.recut.seq[idx+1:])
                            #early stop if CFD goes below self.cfdThres
                            if o.recut.cfd<self.cfdThres:
                                break
                    self.info_p5 = self.info_p5 + "".join(f"5UTR         before mut:{n_mer}\n"
                                                          f"              after mut:{o.recut.seq}\t cfd:{o.recut.cfd}\n"
                                                          f"                 phases:{n_mer_phases}\n")
                    # print(f"5UTR         before mut:{n_mer}\n"
                    #       f"              after mut:{o.recut.seq}\t cfd:{o.recut.cfd}\n"
                    #       f"                 phases:{n_mer_phases}")

                self.info_phase5_5UTR[1] = o.recut.cfd
            #print(f"DIAG_mut_coding: {o.recut.cfd:.6f}")
            if o.recut.seq != n_mer and o.recut.cfd <= cfd: # mutations been made and it decreases cfd
                n_mer_confirm = self.Donor_postMut[n_window: n_window + 23]
                self.info_p5 = self.info_p5 + "".join(f"phase 5 final\n"
                                                      f" orig gRNA: {o.gRNA.seq}\n"
                                                      f"premut seq: {n_mer} cfd:{cfd}\n"
                                                      f"   postmut: {o.recut.seq} cfd:{o.recut.cfd}\n")
                # print(f"phase 5 final for current window\n"
                #       f" orig gRNA: {o.gRNA.seq}\n"
                #       f"premut seq: {n_mer} cfd:{cfd}\n"
                #       f"   postmut: {o.recut.seq} cfd:{o.recut.cfd}")
                #put mut seq back into ssODN
                self.Donor_postMut = self.Donor_postMut[0:n_window] + o.recut.seq + self.Donor_postMut[n_window + 23:]
                highest_cfd = max([highest_cfd, o.recut.cfd]) #update highest cfd
            else:
                highest_cfd = max([highest_cfd, cfd]) #update highest cfd
            n_window+=1
        return highest_cfd

    def slide_win_cfd_coding(self,seq,phases):
        '''
        scan the donor and return the highest cfd
        for coding strand only
        '''
        highest_cfd = 0
         #get sliding window as an iterator
        it = self.sliding_window(seq,23)
        #go through sliding windows #NOTE: Donor_postMut is always in the coding straind
        n_window = 0
        for n_mer in it:
            # #skip non-chimeric part of the homology arm #!!!=> we should not skip non-chimeric parts b/c the payload may contain cutsites
            # if 0 <= n_window <= (len(self.left_flk_seq) - 23 - 1):
            #     #print(f"skipping window: {n_window}")
            #     n_window+=1
            #     continue
            # if n_window >= len(self.left_flk_seq) + len(self.tag) - 1:
            #     n_window+=1
            #     #print(f"skipping window: {n_window}")
            #     continue
            #check for "N"s in the sequence
            n_mer = "".join(n_mer)
            if "N" in n_mer or "n" in n_mer:
                n_window+=1
                continue
            cfd = cfd_score(self.gRNA_seq, n_mer)
            #print(f"DIAG cfd_coding: {cfd:.6f}")
            if cfd>=highest_cfd:
                highest_cfd = cfd
            n_window+=1
        return highest_cfd

    def sliding_window(self, seq, n):
        "Returns a sliding window (of width n) over data from the iterable"
        "   s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...                   "
        it = iter(seq)
        result = tuple(islice(it, n))
        if len(result) == n:
            yield result
        for elem in it:
            result = result[1:] + (elem,)
            yield result

    def mask_phase_in_5UTR(self, coord1, coord2, phases):
        """
        similar to mask_phase_within_3bp_exon_intron_junction
        change 5UTR phase from "0" to "5"
        """
        idx = 0
        tmp_lst = list(phases)
        if coord1<coord2:
            for c in range(coord1,coord2+1):
                flag = check_5UTR(chr = self.ENST_chr, ID = self.ENST_ID, pos = c, loc2posType = self.loc2posType)
                if flag:
                    tmp_lst[idx] = "5"
                idx+=1
        else:
            for c in reversed(range(coord2,coord1+1)):
                flag = check_5UTR(chr = self.ENST_chr, ID = self.ENST_ID, pos = c, loc2posType = self.loc2posType) ##
                if flag:
                    tmp_lst[idx] = "5"
                idx+=1
        return "".join(tmp_lst)

    def mask_phase_within_3bp_exon_intron_junction(self, coord1, coord2, phases):
        """
        example input:
        phases = 0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000123
        coord1 = 24813282 (left of phases)
        coord1 = 24813183 (right of phases) Note: left can > right
        output: phases ("0" changed to "9" whereas it's within 3bp to exon/intron junctions)
                phases ("1/2/3" changed to "8" whereas it's within 3bp to exon/intron junctions)
        """
        idx = 0
        tmp_lst = list(phases)
        if coord1<coord2:
            for c in range(coord1,coord2+1):
                flag1 = check_within_3bp_exon_intron_junction(chr = self.ENST_chr, ID = self.ENST_ID, pos = c, loc2posType = self.loc2posType)
                flag2 = check_cds(chr = self.ENST_chr, ID = self.ENST_ID, pos = c, loc2posType = self.loc2posType)
                if flag2 and flag1:
                    tmp_lst[idx] = "8"
                elif flag1:
                    tmp_lst[idx] = "9"
                idx+=1
        else:
            for c in reversed(range(coord2,coord1+1)):
                flag1 = check_within_3bp_exon_intron_junction(chr = self.ENST_chr, ID = self.ENST_ID, pos = c, loc2posType = self.loc2posType)
                flag2 = check_cds(chr = self.ENST_chr, ID = self.ENST_ID, pos = c, loc2posType = self.loc2posType)
                if flag2 and flag1:
                    tmp_lst[idx] = "8" # idx+1 is caused by reversing (range(coord2,coord1))
                elif flag1:
                    tmp_lst[idx] = "9"
                idx+=1
        return "".join(tmp_lst)

    def get_pam_loc(self):
        """
        return the genomic coordinates of the GG in PAM
        """
        if self.ENST_strand == 1:
            return ([self.gStart+21, self.gStart+22])
        else:
            return ([self.gStart-21, self.gStart-22])


    def to_coding_strand(self, gRNA_seq):
        """
        input:  gRNA sequence in the PAM strand
        output: gRNA sequence in coding strand
        """
        if int(self.ENST_strand) * int(self.gStrand) > 0:
            return gRNA_seq
        else:
            return (str(Seq(gRNA_seq).reverse_complement()))

    def get_uptodate_mut(self):
        """
        return the latest left, right arms, cfd, gRNA seq and phases
        NOTE: gRNA is not consistently on the PAM strand
        """
        if hasattr(self,"cdf_score_post_mut4"):
            left = self.left_flk_seq_CodonMut4
            right = self.right_flk_seq_CodonMut4
            cfd = self.cdf_score_post_mut4
            seq = self.post_mut4_gRNA_seq
            phases = self.post_mut4_gRNA_seq_phases
        elif hasattr(self,"cdf_score_post_mut3"):
            left = self.left_flk_seq_CodonMut3
            right = self.right_flk_seq_CodonMut3
            cfd = self.cdf_score_post_mut3
            seq = self.post_mut3_gRNA_seq
            phases = self.post_mut3_gRNA_seq_phases
        elif hasattr(self,"cdf_score_post_mut2"):
            left = self.left_flk_seq_CodonMut2
            right = self.right_flk_seq_CodonMut2
            cfd = self.cdf_score_post_mut2
            seq = self.post_mut2_gRNA_seq
            phases = self.post_mut2_gRNA_seq_phases
        elif hasattr(self, "cdf_score_post_mut_ins"):
            left = self.left_flk_seq_CodonMut
            right = self.right_flk_seq_CodonMut
            cfd = self.cdf_score_post_mut_ins
            seq = self.post_mut_ins_gRNA_seq
            phases = self.post_mut_ins_gRNA_seq_phases
        else:
            left = self.left_flk_seq
            right = self.right_flk_seq
            cfd = self.cdf_score_post_ins
            seq = self.gRNA_seq
            phases = self.gRNA_seq_phases
        return([left,right,cfd,seq,phases])
    def single_base_muation(self,base):
        mapping = {"A":"t","a":"t",
                   "C":"g","c":"g",
                   "G":"c","g":"c",
                   "T":"a","t":"a"}
        return mapping[base]

    def check_overlap(self,obj1,obj2):
        """
        return true if the overlap part is longer than 50% of either obj1 or obj2
        """
        if min([obj1.start,obj1.end]) <= max([obj2.start,obj2.end]):
            overlap_size = abs(obj1.end - obj2.start)
        if min([obj2.start,obj2.end]) <= max([obj1.start,obj1.end]):
            overlap_size = abs(obj2.end - obj1.start)
        if overlap_size >= 0.5 * abs(obj1.start - obj1.end) or overlap_size >= 0.5 * abs(obj2.start - obj1.end):
            return True
        else:
            return False

    def get_silent_mutations(self, seq):
        """
        input: a stretch of codons
        """
        mutated_subseq = ""
        #mutate
        for mutated_subseq in mutate_silently(seq): # get the last item in the generator
            pass
        #preserve lowercase in input seq
        tmp = list(mutated_subseq)
        for idx, item in enumerate(seq):
            if item.islower():
                tmp[idx] = tmp[idx].lower()
        mutated_subseq = "".join(tmp)

        return mutated_subseq

    def get_trunc_gRNA(self, leftArm, rightArm):
        """
        returns the gRNA sequence between PAM and insert, in coding strand (for later mutations)
        returns the left, right coordinate of the gRNA respective to the whole arm (whole arm is in coding strand)
        all returns are in coding strand
        """
        Lstart = self.left_flk_coord_lst[0]
        Rend = self.right_flk_coord_lst[1]

        if (Lstart<Rend and self.ENST_strand==-1) or (Lstart>Rend and self.ENST_strand==1):
            sys.exit("start<end while strand==-1 or start>end while strand==1")

        newLarm = leftArm
        newRarm = rightArm
        whole_arm = newLarm + newRarm
        whole_arm_ph = self.left_flk_phases + self.right_flk_phases

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
        gRNA = whole_arm[gRNAleft-2:gRNAleft].lower() + whole_arm[gRNAleft:gRNAright+1] + whole_arm[gRNAright+1:gRNAright+1+2].lower()
        gRNA_ph = whole_arm_ph[gRNAleft-2:gRNAleft] + whole_arm_ph[gRNAleft:gRNAright+1] + whole_arm_ph[gRNAright+1:gRNAright+1+2] #get phases
        start = gRNAleft -2 #pad 2bp
        end = gRNAright +2 #pad 2bp

        #get seq between gRNA and cut site:
        trunc_gRNA = gRNA
        trunc_gRNA_ph = gRNA_ph
        if (gRNAleft < int(len(whole_arm)/2) < gRNAright) and (int(self.ENST_strand) * int(self.gStrand) > 0): #gRNA is truncated, and gRNA is on the coding strand
            trunc_gRNA = whole_arm[int(len(whole_arm)/2)-2:int(len(whole_arm)/2)].lower() + whole_arm[int(len(whole_arm)/2):gRNAright+1] + whole_arm[gRNAright+1:gRNAright+1+2].lower() #pad 2bp
            trunc_gRNA_ph = whole_arm_ph[int(len(whole_arm)/2)-2:int(len(whole_arm)/2)] + whole_arm_ph[int(len(whole_arm)/2):gRNAright+1] + whole_arm_ph[gRNAright+1:gRNAright+1+2] #pad 2bp
            start = int(len(whole_arm)/2) -2 #pad 2bp
            end = gRNAright +2 #pad 2bp
        elif (gRNAleft < int(len(whole_arm)/2) < gRNAright) and (int(self.ENST_strand) * int(self.gStrand) < 0): #gRNA is truncated, and gRNA is NOT on the coding strand
            trunc_gRNA = whole_arm[gRNAleft-2:gRNAleft].lower() + whole_arm[gRNAleft:int(len(whole_arm)/2)] + whole_arm[int(len(whole_arm)/2):int(len(whole_arm)/2)+2].lower() #pad 2bp
            trunc_gRNA_ph = whole_arm_ph[gRNAleft-2:gRNAleft] + whole_arm_ph[gRNAleft:int(len(whole_arm)/2)] + whole_arm_ph[int(len(whole_arm)/2):int(len(whole_arm)/2)+2] #pad 2bp
            start = gRNAleft -2 #pad 2bp
            end = int(len(whole_arm)/2) -1 +2 #pad 2bp

        #skipped revcom, b/c the truncated gRNA needs to be in the coding frame for mutation

        return [seq_w_phase(seq = trunc_gRNA, phases = trunc_gRNA_ph, start = start, end = end),gRNAleft,gRNAright]

    def get_post_integration_gRNA(self, leftArm, rightArm):
        """
        input: leftArm, rightArm (based on which the gRNA and insert-disrupted gRNA will be extracted
        return: gRNA, chimeric_gRNA, gRNA_phases, chimeric_gRNA_phases (phases of the insert are denoted with X)
        NOTE: the returned gRNA is in strand same as the PAM for CFD calcualte (the strand is not necessarily the coding strand!)
        #TODO bugfix: crashes when the payload is less than the gRNA length
        """
        Lstart = self.left_flk_coord_lst[0]
        Rend = self.right_flk_coord_lst[1]

        if (Lstart<Rend and self.ENST_strand==-1) or (Lstart>Rend and self.ENST_strand==1):
            sys.exit("start<end while strand==-1 or start>end while strand==1")

        newLarm = leftArm
        newRarm = rightArm
        whole_arm = newLarm + newRarm
        whole_arm_ph = self.left_flk_phases + self.right_flk_phases

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
        gRNA_ph = whole_arm_ph[gRNAleft:gRNAright+1] #get phases

        #get chimeric gRNA:
        chimeric_gRNA = gRNA
        chimeric_gRNA_ph = gRNA_ph

        #print(f"gRNAleft {gRNAleft} {int(len(whole_arm)/2)} gRNAright{gRNAright}")
        if (gRNAleft < int(len(whole_arm)/2) <= gRNAright) and (int(self.ENST_strand) * int(self.gStrand) > 0): #gRNA is truncated, and gRNA is on the coding strand
            chimeric_gRNA = whole_arm[int(len(whole_arm)/2):gRNAright+1]
            trunc_len = int(len(whole_arm)/2) - gRNAleft
            chimeric_gRNA = self.tag[-trunc_len:] + chimeric_gRNA #make the chimeric gRNA
            #chimeric phases
            chimeric_gRNA_ph = whole_arm_ph[int(len(whole_arm)/2):gRNAright+1]
            chimeric_gRNA_ph = "X"*trunc_len + chimeric_gRNA_ph  # make the chimeric gRNA  #TODO missing the left side of Xs

        elif (gRNAleft < int(len(whole_arm)/2) <= gRNAright) and (int(self.ENST_strand) * int(self.gStrand) < 0): #gRNA is truncated, and gRNA is NOT on the coding strand
            chimeric_gRNA = whole_arm[gRNAleft:int(len(whole_arm)/2)]
            trunc_len = gRNAright - int(len(whole_arm)/2) + 1
            chimeric_gRNA = chimeric_gRNA + self.tag[:trunc_len] #make the chimeric gRNA
            #chimeric phases
            chimeric_gRNA_ph = whole_arm_ph[gRNAleft:int(len(whole_arm)/2)]
            chimeric_gRNA_ph = chimeric_gRNA_ph + "X"*trunc_len  #TODO missing the right side of Xs

        #reverse complement gRNA, if needed
        if (int(self.ENST_strand) * int(self.gStrand) < 0):
            #revcom gRNA and trunc_gRNA
            gRNA = str(Seq(gRNA).reverse_complement())
            chimeric_gRNA =  str(Seq(chimeric_gRNA).reverse_complement())
            #reverse phase
            gRNA_ph = gRNA_ph[::-1]
            chimeric_gRNA_ph = chimeric_gRNA_ph[::-1]

        return gRNA, chimeric_gRNA, gRNA_ph, chimeric_gRNA_ph #cannot use the seq_with_phase object b/c not always in the coding strand

    def select_Manu_strand(self, seq):
        """
        select strand:
            if cut-to-insert is 0: use same strand as gRNA
            if recoding is off:
         input: seq
        output: seq in the preferred strand
        """
        if self.Cut2Ins_dist == 0: # use same strand as gRNA,  note the ssODN is initialized to be on the coding strand (same as ENST_strand)
            if self.ENST_strand * self.gStrand > 0:
                return seq
            elif self.ENST_strand * self.gStrand < 0:
                return str(Seq(seq).reverse_complement())

        if self.InsPos <= self.CutPos:
            return seq
        else:
            return str(Seq(seq).reverse_complement())

    def put_silent_mutation_subseq_back(self,L_arm, R_arm, mutated_subseq, start, end):
        """
        start, end are local to the whole arm = L_arm + R_arm
        mutated_subseq has to be in coding strand
        """
        if start == end:
            return L_arm, R_arm
        else:
            if start > end:
                start, end = end, start
            whole_arm = L_arm + R_arm
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
        return [seq, phases, start, end] (the returned seq is always in the coding strand)
        Note, the sequence is padded with 2bp on each side (to avoid truncating codons)
        Note 2: returns two sets of coordinates, both with respective to the whole arm, one is adjusted that the start is always the gRNA start
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
        ins2cut_seq = whole_arm[ins2cutStart-2:ins2cutStart].lower() + whole_arm[ins2cutStart:ins2cutEnd+1] + whole_arm[ins2cutEnd+1:ins2cutEnd+3].lower() #pad 2bp on each side
        ins2cut_phases = whole_arm_phases[ins2cutStart-2:ins2cutStart] + whole_arm_phases[ins2cutStart:ins2cutEnd+1] + whole_arm_phases[ins2cutEnd+1:ins2cutEnd+3] #pad 2bp on each side

        if Lstart>Rend:
            ins2cut_seq = ins2cut_seq[::-1] #reverse
            ins2cut_phases = ins2cut_phases[::-1] #reverse
            #ins2cutStart, ins2cutEnd = ins2cutEnd, ins2cutStart  #swap start and end
            ins2cutStart, ins2cutEnd = (len(whole_arm) - ins2cutEnd - 1), (len(whole_arm) - ins2cutStart - 1 )  #swap start and end

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
        #TODO finish gRNA_obj
        #gRNA_obj = seq_w_phase(seq = gRNA, phases = gRNA_ph, start = gRNAleft, end = gRNAright)

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


def check_within_3bp_exon_intron_junction(chr, ID, pos, loc2posType):
    """
    return True if the position is within_3bp_of_intron_exon_junction or within_3bp_of_exon_intron_junction
    """
    pos_types = _get_position_type(chr, ID, pos, loc2posType)
    if 'within_3bp_of_intron_exon_junction' in set(pos_types) or 'within_3bp_of_exon_intron_junction' in set(pos_types):
        return True
    else:
        return False

def check_5UTR(chr, ID, pos, loc2posType):
    """
    return True if the position is within_3bp_of_intron_exon_junction or within_3bp_of_exon_intron_junction
    """
    pos_types = _get_position_type(chr, ID, pos, loc2posType)
    if '5UTR' in set(pos_types):
        return True
    else:
        return False

def check_cds(chr, ID, pos, loc2posType):
    """
    return True if the position is within_3bp_of_intron_exon_junction or within_3bp_of_exon_intron_junction
    """
    pos_types = _get_position_type(chr, ID, pos, loc2posType)
    if 'cds' in set(pos_types):
        return True
    else:
        return False

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

def check_RE_site(RE,seq):
    """
    RE: restriction site object
    seq: plain sequence
    return a list of restriction sites if found
    ***does not check revcom***
    >>>check_RE_site(Restriction.BsaI,"AAAAAAAAAAAAAA")
    []
    >>>check_RE_site(Restriction.BsaI,"AAAAAGGTCTCAAA")
    [13]
    >>>check_RE_site(Restriction.BsaI,"AAAAAGAGACCAAA")
    []
    """
    return RE.search(Seq(seq))
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
        >>> hdr.target_mutation_score = self.cfdThres
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