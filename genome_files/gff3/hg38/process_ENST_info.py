import os.path
import pandas as pd
from subprocess import Popen
from Bio import SeqIO
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
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
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


class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

def parse_args():
    parser= MyParser(description='This scripts creates a mapping of ENST to chr from gff3')
    parser.add_argument('--gff3_gzfile', default="", type=str, help='path to the gzfile', metavar='')
    config = parser.parse_args()
    if len(sys.argv)==1: # print help message if arguments are not valid
        parser.print_help()
        sys.exit(1)
    return config

logging.setLoggerClass(ColoredLogger)
#logging.basicConfig()
log = logging.getLogger("get_scores")
log.propagate = False
log.setLevel(logging.DEBUG) #set the level of warning displayed
#log.setLevel(logging.INFO) #set the level of warning displayed

config = vars(parse_args())

#####################
##      main       ##
#####################
def main():
    try:
        starttime = datetime.datetime.now()

        file = config["gff3_gzfile"]
        outfile = file.rstrip(".gz") + ".ENST2Chr.gz"

        #dicts
        ENST_info = dict() #dict of seq records
        ENST_exon_dict = dict()  # a dict of dicts
        ENST_exon_dict2 = dict() # a dict of req records
        ENST_CDS_dict = dict()
        loc2exonID_dict = dict() # the same location may have different ENSE IDs, and need ENST to distinguish
        loc2posType = dict() # a dict that maps location to position types (e.g. 5UTR, exon intron, junction, 3UTR)

        #search patterns
        ENST_pattern = re.compile('(ENST.+?);')
        transcript_ID_pattern = re.compile('ID=transcript:(.+?);')
        Parent_pattern = re.compile('Parent=transcript:(.+?);')
        exon_id_pattern = re.compile('exon_id=(.+?);')
        CDS_id_pattern = re.compile('ID=CDS:(.+?);')
        rank_pattern = re.compile('rank=(.+?);')
        phase_pattern = re.compile('ensembl_phase=(.+?);')
        end_phase_pattern = re.compile('ensembl_end_phase=(.+?);')
        name_pattern = re.compile('Name=(.+?);')
        biotype_pattern = re.compile('biotype=(.+?);')
        version_pattern = re.compile('version=(.+?)')

        line_count = 0

        #go through gff3 file, store all ENST IDs
        with gzip.open(file, "rt") as fh:
            for line in fh:
                fields = line.rstrip().split("\t")
                if len(fields)>=2:
                    chr = fields[0]
                    type = fields[2]
                    m = re.search(transcript_ID_pattern, line) #require ID=ENSTXXXX
                    if m:
                        ENST_id=m.group(1)
                        if ENST_id in ENST_info:
                            sys.exit(f"{ENST_id} is already in ENST_info_dict") # ID=transcript:ENSTXXXX should be unique
                        #get transcript parent gene and it's name
                        biotype = ""
                        version = "0"
                        name = ""
                        m2 = re.search(biotype_pattern, line)
                        m3 = re.search(version_pattern, line)
                        m4 = re.search(name_pattern,line)
                        if m2: biotype = m2.group(1)
                        if m3: version = m3.group(1)
                        if m4: name = m4.group(1)
                        description = [f"{ENST_id}.{version}",biotype]
                        ENST_info[ENST_id] = SeqRecord("", id=ENST_id, description="|".join(description), name = name)
                        line_count +=1

        #go through gff3 file, sand store all exons in ENST_exon_dict
        with gzip.open(file, "rt") as fh:
            parent_ENST_id = ""
            for line in fh:
                # print(line.rstrip())
                fields = line.rstrip().split("\t")
                if len(fields) >= 2:
                    chr = fields[0]
                    type = fields[2]
                    m2 = re.search(Parent_pattern, line)
                    if m2 and type == "exon":
                        parent_ENST_id = m2.group(1)
                        if not parent_ENST_id in ENST_info.keys():
                            sys.exit(f"{parent_ENST_id} is not found in ENST_info.keys(), offending line: {line.rstrip()} ")
                        exon_id = exon_loc = f"{parent_ENST_id}__{chr}_{fields[3]}-{fields[4]}_{fields[6]}"  # chr_start-end_strand
                        m3 = re.search(exon_id_pattern, line)
                        if m3:  # there is an exon ID
                            exon_id = m3.group(1)
                            loc2exonID_dict[exon_loc] = exon_id
                        exon_start = int(fields[3])
                        exon_end = int(fields[4])
                        exon_strand = plus_minus_strand_to_numeric(fields[6])
                        #append exon to seq record
                        ENST_info[parent_ENST_id].features.append(SeqFeature(location=FeatureLocation(exon_start,exon_end,strand=exon_strand, ref=chr),type=type, id=exon_id))

        # parse GFF3 again, extracting CDS info, and referencing exon info
        with gzip.open(file, "rt") as fh:
            parent_ENST_id = ""
            for line in fh:
                # print(line.rstrip())
                fields = line.rstrip().split("\t")
                if len(fields)>=2:
                    chr = fields[0]
                    type = fields[2]
                    m2 = re.search(Parent_pattern, line)
                    if m2 and type == "CDS":
                        parent_ENST_id = m2.group(1)
                        if not parent_ENST_id in ENST_info.keys():
                            sys.exit(f"{parent_ENST_id} is not found in ENST_info.keys(), offending line: {line.rstrip()} ")

                        CDS_id = f"{parent_ENST_id}__{chr}_{fields[3]}-{fields[4]}_{fields[6]}"  # chr_start-end_strand
                        # determine if CDS superimposes with an exon
                        if CDS_id in loc2exonID_dict.keys():
                            exon_id = loc2exonID_dict[CDS_id]
                            CDS_id = exon_id

                        # populate CDS info
                        CDS_start = int(fields[3])
                        CDS_end = int(fields[4])
                        CDS_strand = plus_minus_strand_to_numeric(fields[6])
                        CDS_phase = fields[7]
                        #append exon to seq record
                        #print(f"{CDS_start} {CDS_end} {CDS_strand} {type} {CDS_id}")
                        ENST_info[parent_ENST_id].features.append(SeqFeature(location=FeatureLocation(CDS_start,CDS_end,strand=CDS_strand, ref=chr),type=type, id=CDS_id))
                        # copy over additional info from exon dict
                        # if parent_ENST_id in ENST_exon_dict.keys():
                        #     if exon_id in ENST_exon_dict[parent_ENST_id].keys():
                        #         ENST_CDS_dict[parent_ENST_id][CDS_id]["exon_rank"] = ENST_exon_dict[parent_ENST_id][exon_id]["rank"]
                        #         ENST_CDS_dict[parent_ENST_id][CDS_id]["exon_phase"] = ENST_exon_dict[parent_ENST_id][exon_id]["phase"]
                        #         ENST_CDS_dict[parent_ENST_id][CDS_id]["exon_end_phase"] = ENST_exon_dict[parent_ENST_id][exon_id]["end_phase"]

        # write dict to file
        with open('ENST_info.pickle', 'wb') as handle:
            pickle.dump(ENST_info, handle, protocol=pickle.HIGHEST_PROTOCOL)

        #populate loc2posType dict
        for ENST_ID in ENST_info.keys():
            UTR5p, UTR3p = get_UTR_loc(ENST_ID,ENST_info) #get UTR loc
            cds_loc = get_cds_loc(ENST_ID,ENST_info) #get cds loc
            log.debug(f"{ENST_info[ENST_ID].description}")
            log.debug(f"=======================")
            #mark UTRs
            #log.debug(f"ENST {ENST_ID} UTR5p {UTR5p} UTR3p {UTR3p}")

            #mark cds and exon/intron junctions
            for idx, loc in enumerate(cds_loc):
                chr, start, end, strand = loc
                #mark cds
                loc2posType = update_dictOfDict(mydict= loc2posType, key=chr, key2 = f"{start}-{end}", value = "cds")
                #mark junctions
                if idx != 0 and idx != (len(cds_loc)-1): #not first or last cds
                    if strand == 1 or strand == "1" or strand == "+": #pos strand
                        loc2posType = update_dictOfDict(mydict=loc2posType, key=chr, key2=f"{end-1}-{end+2}", value="within_2bp_of_exon_intron_junction")
                        loc2posType = update_dictOfDict(mydict=loc2posType, key=chr, key2=f"{end-3}-{end-2}", value="3N4bp_up_of_exon_intron_junction")
                        loc2posType = update_dictOfDict(mydict=loc2posType, key=chr, key2=f"{end+3}-{end+6}", value="3_to_6bp_down_of_exon_intron_junction")
                        loc2posType = update_dictOfDict(mydict=loc2posType, key=chr, key2=f"{start-2}-{start+1}", value="within_2bp_of_intron_exon_junction")
                        loc2posType = update_dictOfDict(mydict=loc2posType, key=chr, key2=f"{start-4}-{start-3}", value="3N4bp_up_of_intron_exon_junction")
                        loc2posType = update_dictOfDict(mydict=loc2posType, key=chr, key2=f"{start+2}-{start+3}", value="3N4bp_down_of_intron_exon_junction")
                    else: #neg strand
                        loc2posType = update_dictOfDict(mydict=loc2posType, key=chr, key2=f"{start-2}-{end+1}", value="within_2bp_of_exon_intron_junction")
                        loc2posType = update_dictOfDict(mydict=loc2posType, key=chr, key2=f"{start+2}-{start+3}", value="3N4bp_up_of_exon_intron_junction")
                        loc2posType = update_dictOfDict(mydict=loc2posType, key=chr, key2=f"{start-6}-{start-3}", value="3_to_6bp_down_of_exon_intron_junction")
                        loc2posType = update_dictOfDict(mydict=loc2posType, key=chr, key2=f"{end-1}-{end+2}", value="within_2bp_of_intron_exon_junction")
                        loc2posType = update_dictOfDict(mydict=loc2posType, key=chr, key2=f"{end+3}-{end+4}", value="3N4bp_up_of_intron_exon_junction")
                        loc2posType = update_dictOfDict(mydict=loc2posType, key=chr, key2=f"{end-3}-{end-2}" ,value="3N4bp_down_of_intron_exon_junction")
                elif idx == 0:
                    if strand == 1 or strand == "1" or strand == "+": #pos strand, first cds
                        loc2posType = update_dictOfDict(mydict=loc2posType, key=chr, key2=f"{end-1}-{end+2}", value="within_2bp_of_exon_intron_junction")
                        loc2posType = update_dictOfDict(mydict=loc2posType, key=chr, key2=f"{end-3}-{end-2}", value="3N4bp_up_of_exon_intron_junction")
                        loc2posType = update_dictOfDict(mydict=loc2posType, key=chr, key2=f"{end+3}-{end+6}", value="3_to_6bp_down_of_exon_intron_junction")
                    else:  # neg strand last cds
                        loc2posType = update_dictOfDict(mydict=loc2posType, key=chr, key2=f"{end-1}-{end+2}", value="within_2bp_of_intron_exon_junction")
                        loc2posType = update_dictOfDict(mydict=loc2posType, key=chr, key2=f"{end+3}-{end+4}", value="3N4bp_up_of_intron_exon_junction")
                        loc2posType = update_dictOfDict(mydict=loc2posType, key=chr, key2=f"{end-3}-{end-2}" ,value="3N4bp_down_of_intron_exon_junction")
                elif idx == (len(cds_loc)-1):
                    if strand == 1 or strand == "1" or strand == "+":  # pos strand, last cds
                        loc2posType = update_dictOfDict(mydict=loc2posType, key=chr, key2=f"{start-2}-{start+1}", value="within_2bp_of_intron_exon_junction")
                        loc2posType = update_dictOfDict(mydict=loc2posType, key=chr, key2=f"{start-4}-{start-3}", value="3N4bp_up_of_intron_exon_junction")
                        loc2posType = update_dictOfDict(mydict=loc2posType, key=chr, key2=f"{start+2}-{start+3}", value="3N4bp_down_of_intron_exon_junction")
                    else:  # neg strand, first cds
                        loc2posType = update_dictOfDict(mydict=loc2posType, key=chr, key2=f"{start-2}-{end+1}", value="within_2bp_of_exon_intron_junction")
                        loc2posType = update_dictOfDict(mydict=loc2posType, key=chr, key2=f"{start+2}-{start+3}", value="3N4bp_up_of_exon_intron_junction")
                        loc2posType = update_dictOfDict(mydict=loc2posType, key=chr, key2=f"{start-6}-{start-3}", value="3_to_6bp_down_of_exon_intron_junction")



        endtime = datetime.datetime.now()
        elapsed_sec = endtime - starttime
        elapsed_min = elapsed_sec.seconds / 60
        log.info(f"finished in {elapsed_min:.2f} min, processed {line_count} lines {file}")
        print(f"finished in {elapsed_min:.2f} min, processed {line_count} lines {file}", flush=True)

    except Exception as e:
        print("Unexpected error:", str(sys.exc_info()))
        print("additional information:", e)
        PrintException()

##########################
## function definitions ##
##########################

def update_dictOfDict(mydict, key, key2, value):
    """
    make the following update:  mydict[key][key2] = value
    """
    if key in mydict.keys():
        mydict[key][key2] = value
    else:
        mydict[key] = dict()
        mydict[key][key2] = value
    return mydict

def get_UTR_loc(ENST_ID, ENST_info):
    """
    return [UTR5p,UTR3p]  #UTR5p:[start,end]
    """
    locList = []
    my_transcript = ENST_info[ENST_ID]  # get the seq record
    # construct the list of cds
    cdsList = [feat for feat in my_transcript.features if feat.type == 'CDS']
    # construct the list of exon
    exonList = [feat for feat in my_transcript.features if feat.type == 'exon']
    UTR5p = None
    UTR3p = None

    for idx, exon in enumerate(exonList):
        superimpose_cds = find_superimpose_cds(exon,cdsList)
        overlap_cds = find_overlap_cds(exon,cdsList)
        strand = exon.location.strand
        if len(superimpose_cds)==0 and len(overlap_cds)==0: #untranscribed exon
            log.debug(f"ENST {ENST_ID} exon_idx {idx} strand {strand} untranscribed exon")
        elif len(superimpose_cds)==0 and len(overlap_cds)!=0:
            log.debug(f"ENST {ENST_ID} exon_idx {idx} strand {strand} partial-transcribed exon")
        elif len(superimpose_cds)!=0:
            log.debug(f"ENST {ENST_ID} exon_idx {idx} strand {strand} fully-transcribed exon")

    return([UTR5p,UTR3p])

def find_overlap_cds(exon, cdsList):
    """
    return a list of index of cds that overlap with input exon
    """
    overlap_cds = []
    for idx, cds in enumerate(cdsList):
        if cds.location.start <= exon.location.start <= cds.location.end or cds.location.start <= exon.location.end <= cds.location.end or exon.location.start <= cds.location.start <= exon.location.end or exon.location.start <= cds.location.end <= exon.location.end:
            overlap_cds.append(idx)
    return(overlap_cds)

def find_superimpose_cds(exon, cdsList):
    """
    return a list of index of cds that superimpose with input exon
    """
    superimpose_cds = []
    for idx, cds in enumerate(cdsList):
        if cds.location.start == exon.location.start and cds.location.end == exon.location.end:
            superimpose_cds.append(idx)
    return(superimpose_cds)

def get_nonoverlap(loc1, loc2):
    """
    compare two locations [start,end] loc2 [start,end] , one loc encapsulates the other
    return the non-overlapping interval
    """
    if loc1[0]==loc2[0]:
        if loc1[1] != loc2[1]:
            return([min([loc1[1],loc2[1]])+1,max([loc1[1],loc2[1]])])
    elif loc1[1]==loc2[1]:
        if loc1[0] != loc2[0]:
            return([min([loc1[0],loc2[0]]),max([loc1[0],loc2[0]])-1])
    else:
        return(None)

def get_cds_loc(ENST_ID,ENST_info):
    """
    Get the chromosomal location of CDS
    input: ENST_ID, ENST_info
    output: a list loc [chr,start,end,strand] #start < end
    """
    locList = []
    my_transcript = ENST_info[ENST_ID]  # get the seq record
    # constructing the list of cds
    cdsList = [feat for feat in my_transcript.features if feat.type == 'CDS']
    for cds in cdsList:
        locList.append([cds.location.ref, cds.location.start, cds.location.end, cds.strand])
    return(locList)

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
    if len(cdsList)==0:
        return [None, None]
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

def PrintException():
    exc_type, exc_obj, tb = sys.exc_info()
    f = tb.tb_frame
    lineno = tb.tb_lineno
    filename = f.f_code.co_filename
    linecache.checkcache(filename)
    line = linecache.getline(filename, lineno, f.f_globals)
    print('EXCEPTION IN ({}, LINE {} "{}"): {}'.format(filename, lineno, line.strip(), exc_obj))
def plus_minus_strand_to_numeric(input):
    if input == "+":
        return(1)
    elif input == "-":
        return(-1)
    else:
        return(0)
if __name__ == "__main__": main()




