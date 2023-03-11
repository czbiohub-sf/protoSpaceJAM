# Finds OffTargets using bwa
# input: fasta file
# output: bwa_wd_BRCA2_[fasta file name]/bwa.out.bed

import os
from subprocess import Popen
import subprocess
import shutil
import datetime
from itertools import product
from Bio import SeqIO
import gzip
import logging
import argparse
import sys

from protoSpaceJAM.util.utils import MyParser

logging.basicConfig()  # setup the default configuration set
log = logging.getLogger("FinOfftargetBwa")  # Logger name
log.propagate = False
log.setLevel(logging.INFO)  # set the level of warning displayed


def parse_args():
    parser = MyParser(description="This script does XXX")
    parser.add_argument(
        "--fa", default="", type=str, help="name of the input fa", metavar=""
    )
    parser.add_argument(
        "--fa_dir",
        default="",
        type=str,
        help="path to the input fa directory",
        metavar="",
    )
    parser.add_argument(
        "--bin_dir", default="", type=str, help="path to the bin directory", metavar=""
    )
    parser.add_argument(
        "--script_dir",
        default="",
        type=str,
        help="path to the script directory",
        metavar="",
    )
    parser.add_argument(
        "--bwa_idx", default="", type=str, help="path to bwa index", metavar=""
    )
    parser.add_argument(
        "--genome_fa", default="", type=str, help="name of genome fasta", metavar=""
    )
    parser.add_argument(
        "--guideLen", default="", type=str, help="guide length", metavar=""
    )
    parser.add_argument(
        "--thread", default="", type=str, help="num of thread to use", metavar=""
    )

    config = parser.parse_args()
    if len(sys.argv) == 1:  # print help message if arguments are not valid
        parser.print_help()
        sys.exit(1)
    return config


config = vars(parse_args())

GRNA_FA = config["fa"]
GRNA_FA_PATH = os.path.join(config["fa_dir"], config["fa"])
BIN = config["bin_dir"]
SCRIPT = config["script_dir"]
# GENOME = config["genome_fa"]
GENOME_IDX_BWA = config["bwa_idx"]
GUIDELEN = config["guideLen"]
thread2use = config["thread"]

###########################################
# ATTENTION:################################
# Need to copy the GENOME file into /dev/shm

guide_count = 0
with open(GRNA_FA_PATH, "r", encoding="utf-8") as handle:
    fasta_sequences = SeqIO.parse(handle, "fasta")
    for entry in fasta_sequences:
        name, desc, seq = entry.id, entry.description, str(entry.seq)
        guide_count += 1
# print(f"total guide count: {guide_count}\n")


# for some PAMs, we allow other alternative motifs when searching for offtargets
# MIT and eCrisp do that, they use the motif NGG + NAG, we add one more, based on the
# on the guideSeq results in Tsai et al, Nat Biot 2014
# The NGA -> NGG rule was described by Kleinstiver...Young 2015 "Improved Cas9 Specificity..."
# NNGTRRT rule for S. aureus is in the new protocol "SaCas9 User manual"
# ! the length of the alternate PAM has to be the same as the original PAM!
offtargetPams = {
    "NGG": ["NAG", "NGA"],
    # "NGN" : ["GAW"],
    "NGK": ["GAW"],
    "NGA": ["NGG"],
    "NNGRRT": ["NNGRRN"],
    "TTTV": ["TTTN"],
    "ATTN": ["TTTN", "GTTN"],
    "TTYN": ["VTTV", "TRTV"],
}
DEFAULTPAM = "NGG"

# minimum off-target score of standard off-targets (those that end with NGG)
# This should probably be based on the CFD score these days
# But for now, I'll let the user do the filtering
MINSCORE = 0.0
# minimum off-target score for alternative PAM off-targets
# There is not a lot of data to support this cutoff, but it seems
# reasonable to have at least some cutoff, as otherwise we would show
# NAG and NGA like NGG and the data shows clearly that the alternative
# PAMs are not recognized as well as the main NGG PAM.
# so for now, I just filter out very degenerative ones. the best solution
# would be to have a special penalty on the CFD score, but CFS does not
# support non-NGG PAMs (is this actually true?)
ALTPAMMINSCORE = 1.0

# labels and descriptions of eff. scores
scoreDescs = {
    "doench": (
        "Doench '14",
        "Range: 0-100. Linear regression model trained on 880 guides transfected into human MOLM13/NB4/TF1 cells (three genes) and mouse cells (six genes). Delivery: lentivirus. The Fusi score can be considered an updated version this score, as their training data overlaps a lot. See <a target='_blank' href='http://www.nature.com/nbt/journal/v32/n12/full/nbt.3026.html'>Doench et al.</a>",
    ),
    "wuCrispr": (
        "Wu-Crispr",
        "Range 0-100. Aka 'Wong score'. SVM model trained on previously published data. The aim is to identify only a subset of efficient guides, many guides will have a score of 0. Takes into account RNA structure. See <a target='_blank' href='https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0784-0'>Wong et al., Gen Biol 2015</a>",
    ),
    "ssc": (
        "Xu",
        "Range ~ -2 - +2. Aka 'SSC score'. Linear regression model trained on data from &gt;1000 genes in human KBM7/HL60 cells (Wang et al) and mouse (Koike-Yusa et al.). Delivery: lentivirus. Ranges mostly -2 to +2. See <a target='_blank' href='http://genome.cshlp.org/content/early/2015/06/10/gr.191452.115'>Xu et al.</a>",
    ),
    "crisprScan": [
        "Moreno-Mateos",
        "Also called 'CrisprScan'. Range: mostly 0-100. Linear regression model, trained on data from 1000 guides on &gt;100 genes, from zebrafish 1-cell stage embryos injected with mRNA. See <a target=_blank href='http://www.nature.com/nmeth/journal/v12/n10/full/nmeth.3543.html'>Moreno-Mateos et al.</a>. Recommended for guides transcribed <i>in-vitro</i> (T7 promoter). Click to sort by this score. Note that under 'Show all scores', you can find a Doench2016 model trained on Zebrafish scores, Azimuth in-vitro, which should be slightly better than this model for zebrafish.",
    ],
    "wang": (
        "Wang",
        "Range: 0-100. SVM model trained on human cell culture data on guides from &gt;1000 genes. The Xu score can be considered an updated version of this score, as the training data overlaps a lot. Delivery: lentivirus. See <a target='_blank' href='http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3972032/'>Wang et al.</a>",
    ),
    "chariRank": (
        "Chari",
        "Range: 0-100. Support Vector Machine, converted to rank-percent, trained on data from 1235 guides targeting sequences that were also transfected with a lentivirus into human 293T cells. See <a target='_blank' href='http://www.nature.com/nmeth/journal/v12/n9/abs/nmeth.3473.html'>Chari et al.</a>",
    ),
    "fusi": (
        "Doench '16",
        "Aka the 'Fusi-Score', since V4.4 using the version 'Azimuth', scores are slightly different than before April 2018 but very similar (click 'show all' to see the old scores). Range: 0-100. Boosted Regression Tree model, trained on data produced by Doench et al (881 guides, MOLM13/NB4/TF1 cells + unpublished additional data). Delivery: lentivirus. See <a target='_blank' href='http://biorxiv.org/content/early/2015/06/26/021568'>Fusi et al. 2015</a> and <a target='_blank' href='http://www.nature.com/nbt/journal/v34/n2/full/nbt.3437.html'>Doench et al. 2016</a> and <a target=_blank href='https://crispr.ml/'>crispr.ml</a>. Recommended for guides expressed in cells (U6 promoter). Click to sort the table by this score.",
    ),
    "fusiOld": (
        "OldDoench '16",
        "The original implementation of the Doench 2016 score, as received from John Doench. The scores are similar, but not exactly identical to the 'Azimuth' version of the Doench 2016 model that is currently the default on this site, since Apr 2018.",
    ),
    "najm": (
        "Najm 2018",
        "A modified version of the Doench 2016 score ('Azimuth'), by Mudra Hegde for S. aureus Cas9. Range 0-100. See <a target=_blank href='https://www.nature.com/articles/nbt.4048'>Najm et al 2018</a>.",
    ),
    "ccTop": ("CCTop", "The efficiency score used by CCTop, called 'crisprRank'."),
    "aziInVitro": (
        "Azimuth in-vitro",
        "The Doench 2016 model trained on the Moreno-Mateos zebrafish data. Unpublished model, gratefully provided by J. Listgarden. This should be better than Moreno-Mateos, but we have not found the time to evaluate it yet.",
    ),
    "housden": (
        "Housden",
        "Range: ~ 1-10. Weight matrix model trained on data from Drosophila mRNA injections. See <a target='_blank' href='http://stke.sciencemag.org/content/8/393/rs9.long'>Housden et al.</a>",
    ),
    "proxGc": ("ProxGCCount", "Number of GCs in the last 4pb before the PAM"),
    "seqDeepCpf1": (
        "DeepCpf1",
        "Range: ~ 0-100. Convolutional Neural Network trained on ~20k Cpf1 lentiviral guide results. This is the score without DNAse information, 'Seq-DeepCpf1' in the paper. See <a target='_blank' href='https://www.nature.com/articles/nbt.4061'>Kim et al. 2018</a>",
    ),
    "oof": (
        "Out-of-Frame",
        "Range: 0-100. Out-of-Frame score, only for deletions. Predicts the percentage of clones that will carry out-of-frame deletions, based on the micro-homology in the sequence flanking the target site. See <a target='_blank' href='http://www.nature.com/nmeth/journal/v11/n7/full/nmeth.3015.html'>Bae et al. 2014</a>. Click the score to show the predicted deletions.",
    ),
    "lindel": (
        "Lindel",
        "Wei Chen Frameshift ratio (0-100). Predicts probability of a frameshift caused by any type of insertion or deletion. See <a href='https://academic.oup.com/nar/article/47/15/7989/5511473'>Wei Chen et al, Bioinf 2018</a>. Click the score to see the most likely deletions and insertions.",
    ),
}


def parseFastaAsList(fileObj):
    " parse a fasta file, return list (id, seq) "
    seqs = []
    parts = []
    seqId = None
    for line in fileObj:
        line = line.rstrip("\n").rstrip("\r")
        if line.startswith(">"):
            if seqId != None:
                seqs.append((seqId, "".join(parts)))
            seqId = line.lstrip(">")
            parts = []
        else:
            parts.append(line)
    if len(parts) != 0:
        seqs.append((seqId, "".join(parts)))
    return seqs


def expandIupac(seq):
    """ expand all IUPAC characters to nucleotides, returns list.
    >>> expandIupac("NY")
    ['GC', 'GT', 'AC', 'AT', 'TC', 'TT', 'CC', 'CT']
    """
    # http://stackoverflow.com/questions/27551921/how-to-extend-ambiguous-dna-sequence
    d = {
        "A": "A",
        "C": "C",
        "B": "CGT",
        "D": "AGT",
        "G": "G",
        "H": "ACT",
        "K": "GT",
        "M": "AC",
        "N": "GATC",
        "S": "CG",
        "R": "AG",
        "T": "T",
        "W": "AT",
        "V": "ACG",
        "Y": "CT",
        "X": "GATC",
    }
    seqs = []
    for i in product(*[d[j] for j in seq]):
        seqs.append("".join(i))
    return seqs


def annotateBedWithPos(inBed, outBed):
    """
    given an input bed4 and an output bed filename, add an additional column 5 to the bed file
    that is a descriptive text of the chromosome pos (e.g. chr1:1.23 Mbp).
    """
    with open(outBed, "w") as ofh, open(inBed, "r") as fh:
        for line in open(inBed):
            chrom, start = line.split("\t")[:2]
            start = int(start)
            if start > 1000000:
                startStr = "%.2f Mbp" % (float(start) / 1000000)
            else:
                startStr = "%.2f Kbp" % (float(start) / 1000)
            desc = "%s %s" % (chrom, startStr)

            ofh.write(line.rstrip("\n"))
            ofh.write("\t")
            ofh.write(desc)
            ofh.write("\n")


#########
###bwa###
#########

pam = "NGG"

# BWA: allow up to X mismatches
maxMMs = 4

# maximum number of occurences in the genome to get flagged as repeats.
# This is used in bwa samse, when converting the same file
# and for warnings in the table output.
MAXOCC = 60000

# the BWA queue size is 2M by default. We derive the queue size from MAXOCC
MFAC = 2000000 / MAXOCC


# make output dir
bwa_wd_path = os.path.join(config["fa_dir"], f"{GRNA_FA}_bwa")
if os.path.isdir(bwa_wd_path):
    shutil.rmtree(bwa_wd_path)
os.makedirs(bwa_wd_path)

# print(f"wd {os.getcwd()}")
# print(f"workding dir {workding_dir}")
# print(f"bwa_wd_path {bwa_wd_path}")


seqLen = GUIDELEN
maxDiff = maxMMs
bwaM = MFAC * MAXOCC  # -m is queue size in bwa
maxOcc = MAXOCC

bwa_genome_idx = os.path.join(GENOME_IDX_BWA)
saFname = "sai.out"
bedFname = "bwa.out.bed"
matchesBedFname = "matches.bed"
filtMatchesBedFname = "filtMatches.bed"


saFname = os.path.join(bwa_wd_path, f"{saFname}")
bedFname = os.path.join(bwa_wd_path, f"{bedFname}")
matchesBedFname = os.path.join(bwa_wd_path, f"{matchesBedFname}")
filtMatchesBedFname = os.path.join(bwa_wd_path, f"{filtMatchesBedFname}")

starttime = datetime.datetime.now()

command = [
    f"{BIN}/bwa",
    "aln",
    "-o",
    "0",
    "-m",
    f"{bwaM}",
    "-n",
    f"{maxDiff}",  # 	Maximum edit distance if the value is INT, or the fraction of missing alignments given 2% uniform base error rate if FLOAT. In the latter case, the maximum edit distance is automatically chosen for different read lengths.
    "-k",
    f"{maxDiff}",  # Maximum edit distance in the seed
    "-N",  # -N Disable iterative search. All hits with no more than maxDiff differences will be found. This mode is much slower than the default.
    "-l",
    f"{seqLen}",  # Take the first INT subsequence as seed. If INT is larger than the query sequence, seeding will be disabled. For long reads, this option is typically ranged from 25 to 35 for ‘-k 2’. [inf]
    "-t",
    f"{thread2use}",
    f"{bwa_genome_idx}",
    f"{GRNA_FA_PATH}",
]
path_to_stderr_file = os.path.join(bwa_wd_path, f"bwa.aln.stderr.txt")
path_to_stdout_file = saFname
mystdput = open(path_to_stdout_file, "w+")
mystderr = open(path_to_stderr_file, "w+")
p = Popen(command, stdout=mystdput, stderr=mystderr)
p.communicate()  # wait for the commands to process

# samse
# execuated in a one-liner
# command = [f"{BIN}/bwa",
#             "samse",
#             "-n", f"{maxOcc}", #Maximum number of alignments to output in the XA tag for reads paired properly. If a read has more than INT hits, the XA tag will not be written. [3]
#             f"{bwa_genome_idx}",
#             f"{saFname}",
#             f"{GRNA_FA_PATH}"
#             ]
# path_to_stderr_file = os.path.join(bwa_wd_path, f"bwa.samse.stderr.txt")
# path_to_stdout_file = os.path.join(bwa_wd_path, f"samse.out")
# mystdput = open(path_to_stdout_file, 'w+')
# mystderr = open(path_to_stderr_file, 'w+')
# p = Popen(command, stdout=mystdput, stderr=mystderr, universal_newlines=True)
# p.communicate()  # wait for the commands to process

# next 4 blocks
# extract matches and convert to bed format
# execuated in a one-liner
# command = f"cat samse.out | {SCRIPT}/xa2multi.pl > xa2multi.out"
# path_to_stderr_file = os.path.join(bwa_wd_path, f"xa2multi.stderr.txt")
# mystderr = open(path_to_stderr_file, 'w+')
# p = Popen(command,  stdout=subprocess.PIPE, shell=True, stderr=mystderr, universal_newlines=True)
# p.communicate()  # wait for the commands to process

# command = f"cat xa2multi.out | python {SCRIPT}/samToBed {pam} {seqLen} > samtobed.out"
# path_to_stderr_file = os.path.join(bwa_wd_path, f"samToBed.stderr.txt")
# mystderr = open(path_to_stderr_file, 'w+')
# p = Popen(command,  stdout=subprocess.PIPE, shell=True, stderr=mystderr, universal_newlines=True)
# p.communicate()  # wait for the commands to process

# command = f"cat samtobed.out |sort -k1,1 -k2,2n > sort.out"
# path_to_stderr_file = os.path.join(bwa_wd_path, f"sort.stderr.txt")
# mystderr = open(path_to_stderr_file, 'w+')
# p = Popen(command,  stdout=subprocess.PIPE, shell=True, stderr=mystderr, universal_newlines=True)
# p.communicate()  # wait for the commands to process

# command = f"cat sort.out | {BIN}/bedClip stdin {bwa_genome_idx}.sizes stdout >> {matchesBedFname} "
# path_to_stderr_file = os.path.join(bwa_wd_path, f"sorbedClip.stderr.txt")
# mystderr = open(path_to_stderr_file, 'w+')
# p = Popen(command,  stdout=subprocess.PIPE, shell=True, stderr=mystderr, universal_newlines=True)
# p.communicate()

# one-liner command that parses the matches and convert to bed format
path_to_stderr_file = os.path.join(bwa_wd_path, f"bwa.extract_matches_pipe.stderr.txt")
mystderr = open(path_to_stderr_file, "w+")
command = f"{BIN}/bwa samse -n {maxOcc} {bwa_genome_idx} {saFname} {GRNA_FA_PATH} | {SCRIPT}/xa2multi.pl | python {SCRIPT}/samToBed {pam} {seqLen} | sort -k1,1 -k2,2n | {BIN}/bedClip stdin {bwa_genome_idx}.sizes stdout >> {matchesBedFname} "
p = Popen(command, stdout=subprocess.PIPE, shell=True, stderr=mystderr)
p.communicate()  # wait for the commands to process

# queue.startStep(batchId, "filter", "Removing matches without a PAM motif")
altPats = ",".join(offtargetPams.get(pam, ["na"]))
altPamMinScore = str(ALTPAMMINSCORE)
# shmFaFname = os.path.join("/dev/shm", GENOME) #bwa_genome_idx is the name ends with .fa

# EXTRACTION OF SEQUENCES + ANNOTATION - big headache!!
# twoBitToFa was 15x slower than python's twobitreader, after markd's fix it is better
# but bedtools uses an fa.idx file and also mmap, so is a LOT faster
# arguments: guideSeq, mainPat, altPats, altScore, passTotalAlnCount
path_to_stderr_file = os.path.join(bwa_wd_path, f"bwa.seqExtraction.stderr.txt")
mystderr = open(path_to_stderr_file, "w+")

# if os.path.isfile(shmFaFname):
# logging.info("Using bedtools and genome fasta on ramdisk, %s" % shmFaFname)
# use genome fa in /dev/shm
# command = f"time {BIN}/bedtools getfasta -s -name -fi {shmFaFname} -bed {matchesBedFname} -fo /dev/stdout | python {SCRIPT}/filterFaToBed {GRNA_FA_PATH} {pam} {altPats} {altPamMinScore} > {filtMatchesBedFname}"
# don't use genome fa in /dev/shm
command = f"time {BIN}/bedtools getfasta -s -name -fi {GENOME_IDX_BWA} -bed {matchesBedFname} -fo /dev/stdout | python {SCRIPT}/filterFaToBed {GRNA_FA_PATH} {pam} {altPats} {altPamMinScore} > {filtMatchesBedFname}"

# else:
#    logging.info("Using twoBitfofa, %s" % shmFaFname)
#    command = f"time {BIN}/twoBitToFa {bwa_genome_idx}.2bit stdout -bed={matchesBedFname} | python {SCRIPT}/filterFaToBed {GRNA_FA_PATH} {pam} {altPats} {altPamMinScore} > {filtMatchesBedFname}"


# cmd = "$SCRIPT/twoBitToFaPython %(genomeDir)s/%(genome)s/%(genome)s.2bit %(matchesBedFname)s | $SCRIPT/filterFaToBed %(faFname)s %(pam)s %(altPats)s %(altPamMinScore)s %(maxOcc)d > %(filtMatchesBedFname)s" % locals()
p = Popen(command, stdout=subprocess.PIPE, shell=True, stderr=mystderr)
p.communicate()  # wait for the commands to process

# annotated Bed with chromosome position
bedFnameTmp = bedFname + ".tmp"
annotateBedWithPos(filtMatchesBedFname, bedFnameTmp)
shutil.move(bedFnameTmp, bedFname)

# remove the temporary files
tempFnames = [saFname, matchesBedFname, filtMatchesBedFname]
for tfn in tempFnames:
    if os.path.isfile(tfn):
        os.remove(tfn)

endtime = datetime.datetime.now()
elapsed = endtime - starttime
elapsedmin = elapsed.seconds / 60
print(
    f"FindOfftargetsBwa.py finished in {elapsedmin:.2f} min ({elapsed} sec)", flush=True
)
