import os
from subprocess import Popen
import shutil
import datetime
from itertools import product
from Bio import SeqIO

GRNA_FA = "BRCA2_guides.fa"
GRNA_FA_PATH = os.path.join("/home/duo.peng/github_repos/protospaceXS/test_off_target",GRNA_FA)
BIN="/home/duo.peng/github_repos/protospaceXS/test_off_target/Linux"
GENOME = "Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
GENOME_IDX_BWA = "/home/duo.peng/github_repos/protospaceXS/genome_files/indexes_bwa/" + GENOME
GENOME_IDX_BOWTIE = "/home/duo.peng/github_repos/protospaceXS/genome_files/indexes_bowtie/" + GENOME
GUIDELEN = 20


guide_count = 0
with open(GRNA_FA_PATH, "r", encoding="utf-8") as handle: 
    fasta_sequences = SeqIO.parse(handle,'fasta')
    for entry in fasta_sequences:
        name, desc,seq = entry.id, entry.description, str(entry.seq)
        guide_count +=1
print(f"total guide count: {guide_count}\n")


# for some PAMs, we allow other alternative motifs when searching for offtargets
# MIT and eCrisp do that, they use the motif NGG + NAG, we add one more, based on the
# on the guideSeq results in Tsai et al, Nat Biot 2014
# The NGA -> NGG rule was described by Kleinstiver...Young 2015 "Improved Cas9 Specificity..."
# NNGTRRT rule for S. aureus is in the new protocol "SaCas9 User manual"
# ! the length of the alternate PAM has to be the same as the original PAM!
offtargetPams = {
    "NGG" : ["NAG","NGA"],
    #"NGN" : ["GAW"],
    "NGK" : ["GAW"],
    "NGA" : ["NGG"],
    "NNGRRT" : ["NNGRRN"],
    "TTTV" : ["TTTN"],
    'ATTN' : ["TTTN", "GTTN"],
    "TTYN" : ["VTTV", "TRTV"]
}
DEFAULTPAM = 'NGG'

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
    "doench" : ("Doench '14", "Range: 0-100. Linear regression model trained on 880 guides transfected into human MOLM13/NB4/TF1 cells (three genes) and mouse cells (six genes). Delivery: lentivirus. The Fusi score can be considered an updated version this score, as their training data overlaps a lot. See <a target='_blank' href='http://www.nature.com/nbt/journal/v32/n12/full/nbt.3026.html'>Doench et al.</a>"),
    "wuCrispr" : ("Wu-Crispr", "Range 0-100. Aka 'Wong score'. SVM model trained on previously published data. The aim is to identify only a subset of efficient guides, many guides will have a score of 0. Takes into account RNA structure. See <a target='_blank' href='https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0784-0'>Wong et al., Gen Biol 2015</a>"),
    "ssc" : ("Xu", "Range ~ -2 - +2. Aka 'SSC score'. Linear regression model trained on data from &gt;1000 genes in human KBM7/HL60 cells (Wang et al) and mouse (Koike-Yusa et al.). Delivery: lentivirus. Ranges mostly -2 to +2. See <a target='_blank' href='http://genome.cshlp.org/content/early/2015/06/10/gr.191452.115'>Xu et al.</a>"),
    "crisprScan" : ["Moreno-Mateos", "Also called 'CrisprScan'. Range: mostly 0-100. Linear regression model, trained on data from 1000 guides on &gt;100 genes, from zebrafish 1-cell stage embryos injected with mRNA. See <a target=_blank href='http://www.nature.com/nmeth/journal/v12/n10/full/nmeth.3543.html'>Moreno-Mateos et al.</a>. Recommended for guides transcribed <i>in-vitro</i> (T7 promoter). Click to sort by this score. Note that under 'Show all scores', you can find a Doench2016 model trained on Zebrafish scores, Azimuth in-vitro, which should be slightly better than this model for zebrafish."],
    "wang" : ("Wang", "Range: 0-100. SVM model trained on human cell culture data on guides from &gt;1000 genes. The Xu score can be considered an updated version of this score, as the training data overlaps a lot. Delivery: lentivirus. See <a target='_blank' href='http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3972032/'>Wang et al.</a>"),
    "chariRank" : ("Chari", "Range: 0-100. Support Vector Machine, converted to rank-percent, trained on data from 1235 guides targeting sequences that were also transfected with a lentivirus into human 293T cells. See <a target='_blank' href='http://www.nature.com/nmeth/journal/v12/n9/abs/nmeth.3473.html'>Chari et al.</a>"),
    "fusi" : ("Doench '16", "Aka the 'Fusi-Score', since V4.4 using the version 'Azimuth', scores are slightly different than before April 2018 but very similar (click 'show all' to see the old scores). Range: 0-100. Boosted Regression Tree model, trained on data produced by Doench et al (881 guides, MOLM13/NB4/TF1 cells + unpublished additional data). Delivery: lentivirus. See <a target='_blank' href='http://biorxiv.org/content/early/2015/06/26/021568'>Fusi et al. 2015</a> and <a target='_blank' href='http://www.nature.com/nbt/journal/v34/n2/full/nbt.3437.html'>Doench et al. 2016</a> and <a target=_blank href='https://crispr.ml/'>crispr.ml</a>. Recommended for guides expressed in cells (U6 promoter). Click to sort the table by this score."),
    "fusiOld" : ("OldDoench '16", "The original implementation of the Doench 2016 score, as received from John Doench. The scores are similar, but not exactly identical to the 'Azimuth' version of the Doench 2016 model that is currently the default on this site, since Apr 2018."),
    "najm" : ("Najm 2018", "A modified version of the Doench 2016 score ('Azimuth'), by Mudra Hegde for S. aureus Cas9. Range 0-100. See <a target=_blank href='https://www.nature.com/articles/nbt.4048'>Najm et al 2018</a>."),
    "ccTop" : ("CCTop", "The efficiency score used by CCTop, called 'crisprRank'."),
    "aziInVitro" : ("Azimuth in-vitro", "The Doench 2016 model trained on the Moreno-Mateos zebrafish data. Unpublished model, gratefully provided by J. Listgarden. This should be better than Moreno-Mateos, but we have not found the time to evaluate it yet."),
    "housden" : ("Housden", "Range: ~ 1-10. Weight matrix model trained on data from Drosophila mRNA injections. See <a target='_blank' href='http://stke.sciencemag.org/content/8/393/rs9.long'>Housden et al.</a>"),
    "proxGc" : ("ProxGCCount", "Number of GCs in the last 4pb before the PAM"),
    "seqDeepCpf1" : ("DeepCpf1", "Range: ~ 0-100. Convolutional Neural Network trained on ~20k Cpf1 lentiviral guide results. This is the score without DNAse information, 'Seq-DeepCpf1' in the paper. See <a target='_blank' href='https://www.nature.com/articles/nbt.4061'>Kim et al. 2018</a>"),
    "oof" : ("Out-of-Frame", "Range: 0-100. Out-of-Frame score, only for deletions. Predicts the percentage of clones that will carry out-of-frame deletions, based on the micro-homology in the sequence flanking the target site. See <a target='_blank' href='http://www.nature.com/nmeth/journal/v11/n7/full/nmeth.3015.html'>Bae et al. 2014</a>. Click the score to show the predicted deletions."),
    "lindel": ("Lindel", "Wei Chen Frameshift ratio (0-100). Predicts probability of a frameshift caused by any type of insertion or deletion. See <a href='https://academic.oup.com/nar/article/47/15/7989/5511473'>Wei Chen et al, Bioinf 2018</a>. Click the score to see the most likely deletions and insertions.")
}

#########
###bwa###
#########

# BWA: allow up to X mismatches
maxMMs=4

# maximum number of occurences in the genome to get flagged as repeats.
# This is used in bwa samse, when converting the same file
# and for warnings in the table output.
MAXOCC = 60000

# the BWA queue size is 2M by default. We derive the queue size from MAXOCC
MFAC = 2000000/MAXOCC



def parseFastaAsList(fileObj):
    " parse a fasta file, return list (id, seq) "
    seqs = []
    parts = []
    seqId = None
    for line in fileObj:
        line = line.rstrip("\n").rstrip("\r")
        if line.startswith(">"):
            if seqId!=None:
                seqs.append( (seqId, "".join(parts)) )
            seqId = line.lstrip(">")
            parts = []
        else:
            parts.append(line)
    if len(parts)!=0:
        seqs.append( (seqId, "".join(parts)) )
    return seqs

def expandIupac(seq):
    """ expand all IUPAC characters to nucleotides, returns list.
    >>> expandIupac("NY")
    ['GC', 'GT', 'AC', 'AT', 'TC', 'TT', 'CC', 'CT']
    """
    # http://stackoverflow.com/questions/27551921/how-to-extend-ambiguous-dna-sequence
    d = {'A': 'A', 'C': 'C', 'B': 'CGT', 'D': 'AGT', 'G': 'G', \
        'H': 'ACT', 'K': 'GT', 'M': 'AC', 'N': 'GATC', 'S': 'CG', \
        'R': 'AG', 'T': 'T', 'W': 'AT', 'V': 'ACG', 'Y': 'CT', 'X': 'GATC'}
    seqs = []
    for i in product(*[d[j] for j in seq]):
       seqs.append("".join(i))
    return seqs

current_dir = os.getcwd()


bwa_path = os.path.join(current_dir,"bwa")
if os.path.isdir(bwa_path):
    shutil.rmtree(bwa_path)
os.makedirs(bwa_path)

os.chdir(bwa_path)

seqLen = GUIDELEN
maxDiff = maxMMs
bwaM = MFAC*MAXOCC # -m is queue size in bwa
maxOcc = MAXOCC

bwa_genome_idx = os.path.join(GENOME_IDX_BWA)
saFname = "sai.out"

starttime = datetime.datetime.now()

command = [f"{BIN}/bwa",
            "aln",
            "-o", "0",
            "-m", f"{bwaM}",
            "-n", f"{maxDiff}",
            "-k", f"{maxDiff}",
            "-N",
            "-l", f"{seqLen}",
            f"{bwa_genome_idx}",
            f"{GRNA_FA_PATH}"
            ]
path_to_stderr_file = os.path.join(f"bwa.aln.stderr.txt")
path_to_stdout_file = os.path.join(f"{saFname}")
mystdput = open(path_to_stdout_file, 'w+')
mystderr = open(path_to_stderr_file, 'w+')
p = Popen(command, stdout=mystdput, stderr=mystderr)
p.communicate()  # wait for the commands to process

command = [f"{BIN}/bwa",
            "samse",
            "-n", f"{maxOcc}",
            f"{bwa_genome_idx}",
            f"{saFname}",
            f"{GRNA_FA_PATH}"
            ]
path_to_stderr_file = os.path.join(f"bwa.samse.stderr.txt")
path_to_stdout_file = os.path.join(f"bwa.out")
mystdput = open(path_to_stdout_file, 'w+')
mystderr = open(path_to_stderr_file, 'w+')
p = Popen(command, stdout=mystdput, stderr=mystderr, universal_newlines=True)
p.communicate()  # wait for the commands to process


filtMatchesBedFname = batchBase+".filtMatches.bed"
#queue.startStep(batchId, "filter", "Removing matches without a PAM motif")
altPats = ",".join(offtargetPams.get(pam, ["na"]))
bedFnameTmp = bedFname+".tmp"
altPamMinScore = str(ALTPAMMINSCORE)
shmFaFname = join("/dev/shm", genome+".fa")





endtime = datetime.datetime.now()
elapsed = endtime - starttime
elapsedmin = elapsed.seconds/60
print(f"bwa finished in {elapsedmin:.2f} min ({elapsed} sec)")


########
#bowtie#
########
def makeVariants(seq):
    " generate all possible variants of sequence at 1bp-distance"
    seqs = []
    for i in range(0, len(seq)):
        for l in "ACTG":
            if l==seq[i]:
                continue
            newSeq = seq[:i]+l+seq[i+1:]
            seqs.append((i, seq[i], l, newSeq))
    return seqs

def writeBowtieSequences(inFaFname, outFname, pamPat):
    """ write the sequence and one-bp-distant-sequences + all possible PAM sequences to outFname
    Return dict querySeqId -> querySeq and a list of all
    possible PAMs, as nucleotide sequences (not IUPAC-patterns)
    """
    ofh = open(outFname, "w")
    outCount = 0
    inCount = 0
    guideSeqs = {} # 20mer guide sequences
    qSeqs = {} # 23mer query sequences for bowtie, produced by expanding guide sequences
    allPamSeqs = expandIupac(pamPat)
    for seqId, seq in parseFastaAsList(open(inFaFname)):
        inCount += 1
        guideSeqs[seqId] = seq
        for pamSeq in allPamSeqs:
            # the input sequence + the PAM
            newSeqId = "%s.%s" % (seqId, pamSeq)
            newFullSeq = seq+pamSeq
            ofh.write(">%s\n%s\n" % (newSeqId, newFullSeq))
            qSeqs[newSeqId] = newFullSeq

            # all one-bp mutations of the input sequence + the PAM
            for nPos, fromNucl, toNucl, newSeq in makeVariants(seq):
                newSeqId = "%s.%s.%d:%s>%s" % (seqId, pamSeq, nPos, fromNucl, toNucl)
                newFullSeq = newSeq+pamSeq
                ofh.write(">%s\n%s\n" % (newSeqId, newFullSeq))
                qSeqs[newSeqId] = newFullSeq
                outCount += 1
    ofh.close()
    #logging.debug("Wrote %d variants+expandedPam of %d sequences to %s" % (outCount, inCount, outFname))
    return guideSeqs, qSeqs, allPamSeqs
    
bowtie_path = os.path.join(current_dir,"bowtie")
if os.path.isdir(bowtie_path):
    shutil.rmtree(bowtie_path)
os.makedirs(bowtie_path)

bowtie_genome_idx = os.path.join(GENOME_IDX_BOWTIE)

# write out the input sequences for bowtie
pamPat = 'NGG'
bwFaFname = os.path.join(bowtie_path, "bowtieIn.fa")
guideSeqs, qSeqs, allPamSeqs = writeBowtieSequences(GRNA_FA_PATH, bwFaFname, pamPat)

os.chdir(bowtie_path)

command = [f"{BIN}/bowtie",
            f"-e", f"1000",
            f"{bowtie_genome_idx}",
            f"-f", f"{bwFaFname}",
            f"-v", f"3",
            f"-y", 
            f"-t",
            f"-k", f"{maxOcc}",
            f"-m", f"{maxOcc}",
            f"dummy",
            f"--max", "tooManyHits.txt",
            f"--mm",
            #f"--refout", #this option is deprecated. the default is not refout
            f"--maxbts=2000",
            f"-p", f"4",
]

starttime = datetime.datetime.now()

path_to_stderr_file = os.path.join(f"bowtie.stderr.txt")
path_to_stdout_file = os.path.join(f"bowtie.out")
mystdput = open(path_to_stdout_file, 'w+')
mystderr = open(path_to_stderr_file, 'w+')
p = Popen(command, stdout=mystdput, stderr=mystderr, universal_newlines=True)
p.communicate()  # wait for the commands to process


endtime = datetime.datetime.now()
elapsed = endtime - starttime
elapsedmin = elapsed.seconds/60
print(f"bowtie finished in {elapsedmin:.2f} min ({elapsed} sec)")





