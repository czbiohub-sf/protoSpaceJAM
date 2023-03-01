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
import itertools
import traceback
import os
import urllib
import zlib

#################
# custom logging #
#################
import logging

from protoSpaceJAM.util.utils import MyParser

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

        console = logging.StreamHandler()
        console.setFormatter(color_formatter)

        self.addHandler(console)
        return


def parse_args():
    parser = MyParser(
        description="This scripts processes the gff3 file, generates a dictionary files (pickle format) of (1) ENST transcript information, (2) loc2posType, (3) codon phase information. The dictionaries are input for main.py"
    )
    # parser.add_argument('--release', default="106", type=str, help='Ensembl release', metavar='')
    # parser.add_argument('--genome_ver', default="GRCh38", type=str, help='genome+version, choices are GRCh38, GRCz11, GRCm39', metavar='')
    parser.add_argument(
        "--gff3_gz",
        default="../genome_files/Homo_sapiens.GRCh38.107.gff3.gz",
        type=str,
        help="path to gff3.gz file",
        metavar="",
    )
    parser.add_argument(
        "--out_dir",
        default="../genome_files/parsed_gff3/GRCh38",
        type=str,
        help="path to outputdir",
        metavar="",
    )
    config = parser.parse_args()

    return config


logging.setLoggerClass(ColoredLogger)
# logging.basicConfig()
log = logging.getLogger("process_ENST_info")
log.propagate = False
# log.setLevel(logging.DEBUG) #set the level of warning displayed
log.setLevel(logging.INFO)  # set the level of warning displayed

config = vars(parse_args())

#####################
##      main       ##
#####################
# TODO: run zebrafish
# TODO: output to genome dir
def main():
    try:
        starttime = datetime.datetime.now()

        # genome_ver = config["genome_ver"]
        # release = config["release"]
        # prefix_mapping = {"GRCz": "Danio_rerio",
        #           "GRCh": "Homo_sapiens",
        #           "GRCm": "Mus_musculus"}
        # spp = prefix_mapping[genome_ver[0:4]]

        # file_directory = f"genome_files/gff3/{genome_ver}" # this is the output directory
        # file = spp + '.' + genome_ver + '.' + release + ".gff3.gz"
        file_path = config["gff3_gz"]
        out_dir = config["out_dir"]

        print(f"processing {file_path}")

        # prepare gff3.gz
        # if os.path.isfile(file_path):
        #     print(f"found {file} in directory {file_directory}\nchecking md5 checksums")
        #     #check md5
        #     url=f"http://ftp.ensembl.org/pub/release-{release}/gff3/" + spp.lower() + "/CHECKSUMS"
        #     outmd5 = f"{file_directory}/CHECKSUMS"
        #     # Download the file from `url` and save it locally under `file_name`:
        #     with urllib.request.urlopen(url) as response, open(outmd5, 'wb') as out_file:
        #         shutil.copyfileobj(response, out_file)
        #     # md5_pass = False
        #     # with open(outmd5, "r") as f:
        #     #     for line in f:
        #     #         checksum, fsize, fname = line.rstrip().split()
        #     #         if fname == file and checksum == crc32(file_path):
        #     #             print("md5 checksums matched")
        #     #             md5_pass==True
        # else:
        #     #download
        #     print(f"downloading {file} to directory {file_directory}")
        #     if make_output_dir(file_directory):
        #         url = f"http://ftp.ensembl.org/pub/release-{release}/gff3/" + spp.lower() + "/" + file
        #         # Download the file from `url` and save it locally under `file_name`:
        #         with urllib.request.urlopen(url) as response, open(file_path, 'wb') as out_file:
        #             shutil.copyfileobj(response, out_file)
        #     else:
        #         sys.exit(f"failed to create directory genome_files/gff3/{genome_ver}")

        file = file_path
        outfile = file.rstrip(".gz") + ".ENST2Chr.gz"

        # dicts
        ENST_info = dict()  # dict of seq records
        ENST_exon_dict = dict()  # a dict of dicts
        ENST_exon_dict2 = dict()  # a dict of req records
        ENST_CDS_dict = dict()
        ENST_codons_dict = dict()
        loc2exonID_dict = (
            dict()
        )  # the same location may have different ENSE IDs, and need ENST to distinguish
        loc2posType = (
            dict()
        )  # a dict that maps location to position types (e.g. 5UTR, exon intron, junction, 3UTR)

        # search patterns
        ENST_pattern = re.compile("(ENST.+?);")
        transcript_ID_pattern = re.compile("ID=transcript:(.+?);")
        Parent_pattern = re.compile("Parent=transcript:(.+?);")
        exon_id_pattern = re.compile("exon_id=(.+?);")
        CDS_id_pattern = re.compile("ID=CDS:(.+?);")
        rank_pattern = re.compile("rank=(.+?);")
        phase_pattern = re.compile("ensembl_phase=(.+?);")
        end_phase_pattern = re.compile("ensembl_end_phase=(.+?);")
        name_pattern = re.compile("Name=(.+?);")
        biotype_pattern = re.compile("biotype=(.+?);")
        version_pattern = re.compile("version=(.+?)")

        line_count = 0

        # go through gff3 file, store all ENST IDs
        print("-processing ENST IDs")
        with gzip.open(file, "rt") as fh:
            for line in fh:
                fields = line.rstrip().split("\t")
                if len(fields) >= 2:
                    chr = fields[0]
                    type = fields[2]
                    m = re.search(transcript_ID_pattern, line)  # require ID=ENSTXXXX
                    if m:
                        ENST_id = m.group(1)
                        if ENST_id in ENST_info:
                            sys.exit(
                                f"{ENST_id} is already in ENST_info_dict"
                            )  # ID=transcript:ENSTXXXX should be unique
                        # get transcript parent gene and it's name
                        biotype = ""
                        version = "0"
                        name = ""
                        m2 = re.search(biotype_pattern, line)
                        m3 = re.search(version_pattern, line)
                        m4 = re.search(name_pattern, line)
                        if m2:
                            biotype = m2.group(1)
                        if m3:
                            version = m3.group(1)
                        if m4:
                            name = m4.group(1)
                        description = [f"{ENST_id}.{version}", biotype]
                        ENST_info[ENST_id] = SeqRecord(
                            "", id=ENST_id, description="|".join(description), name=name
                        )
                        line_count += 1

        # go through gff3 file, and store all exons in ENST_exon_dict
        print("-processing exons")
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
                            sys.exit(
                                f"{parent_ENST_id} is not found in ENST_info.keys(), offending line: {line.rstrip()} "
                            )
                        exon_id = (
                            exon_loc
                        ) = f"{parent_ENST_id}__{chr}_{fields[3]}-{fields[4]}_{fields[6]}"  # chr_start-end_strand
                        m3 = re.search(exon_id_pattern, line)
                        if m3:  # there is an exon ID
                            exon_id = m3.group(1)
                            loc2exonID_dict[exon_loc] = exon_id
                        exon_start = int(fields[3])
                        exon_end = int(fields[4])
                        exon_strand = plus_minus_strand_to_numeric(fields[6])
                        # append exon to seq record
                        ENST_info[parent_ENST_id].features.append(
                            SeqFeature(
                                location=FeatureLocation(
                                    exon_start, exon_end, strand=exon_strand, ref=chr
                                ),
                                type=type,
                                id=exon_id,
                            )
                        )

        # parse GFF3 again, extracting CDS info, and referencing exon info
        print("-processing cds")
        with gzip.open(file, "rt") as fh:
            parent_ENST_id = ""
            for line in fh:
                # print(line.rstrip())
                fields = line.rstrip().split("\t")
                if len(fields) >= 2:
                    chr = fields[0]
                    type = fields[2]
                    m2 = re.search(Parent_pattern, line)
                    if m2 and type == "CDS":
                        parent_ENST_id = m2.group(1)
                        if not parent_ENST_id in ENST_info.keys():
                            sys.exit(
                                f"{parent_ENST_id} is not found in ENST_info.keys(), offending line: {line.rstrip()} "
                            )

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
                        # append exon to seq record
                        # print(f"{CDS_start} {CDS_end} {CDS_strand} {type} {CDS_id}")
                        ENST_info[parent_ENST_id].features.append(
                            SeqFeature(
                                location=FeatureLocation(
                                    CDS_start, CDS_end, strand=CDS_strand, ref=chr
                                ),
                                type=type,
                                id=CDS_id,
                            )
                        )
                        # copy over additional info from exon dict
                        # if parent_ENST_id in ENST_exon_dict.keys():
                        #     if exon_id in ENST_exon_dict[parent_ENST_id].keys():
                        #         ENST_CDS_dict[parent_ENST_id][CDS_id]["exon_rank"] = ENST_exon_dict[parent_ENST_id][exon_id]["rank"]
                        #         ENST_CDS_dict[parent_ENST_id][CDS_id]["exon_phase"] = ENST_exon_dict[parent_ENST_id][exon_id]["phase"]
                        #         ENST_CDS_dict[parent_ENST_id][CDS_id]["exon_end_phase"] = ENST_exon_dict[parent_ENST_id][exon_id]["end_phase"]

        # go through ENST_info, calculate span_start span_end for each ID
        for ID in ENST_info:
            if len(ENST_info[ID].features) > 0:
                coords = [i.location.nofuzzy_end for i in ENST_info[ID].features] + [
                    i.location.nofuzzy_start for i in ENST_info[ID].features
                ]
                ENST_info[ID].span_start = min(coords)
                ENST_info[ID].span_end = max(coords)
                ENST_info[ID].chr = ENST_info[ID].features[0].ref

        # write dict to file
        with open(f"{out_dir}/ENST_info.pickle", "wb") as handle:
            pickle.dump(ENST_info, handle, protocol=pickle.HIGHEST_PROTOCOL)

        # split ENST_into into fileparts
        # split the dict and write to pickle files
        ENST_info_splitfile = (
            {}
        )  # indexing dictionary, the values of this dict is an int (file part number)
        ENST_info_outdir = os.path.join(out_dir, "ENST_info")
        if os.path.exists(ENST_info_outdir):
            shutil.rmtree(ENST_info_outdir)
        os.makedirs(ENST_info_outdir)
        # initialization
        part = 0
        temp_dict = {}
        for ENST in ENST_info.keys():  # loop through all ENST
            # populate tmp dict
            temp_dict[ENST] = ENST_info[ENST]
            # update indexing dictionary
            ENST_info_splitfile[ENST] = part
            if len(temp_dict) >= 200:
                # write tmp dict to file
                with open(
                    f"{ENST_info_outdir}/ENST_info_part{part}.pickle", "wb"
                ) as handle:
                    pickle.dump(temp_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
                temp_dict = {}  # reset temp dict
                part += 1
        # write last chunk of the current chromosome
        if len(temp_dict) > 0:
            with open(
                f"{ENST_info_outdir}/ENST_info_part{part}.pickle", "wb"
            ) as handle:
                pickle.dump(temp_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
                part += 1
        # write indexing dictionary
        with open(f"{ENST_info_outdir}/ENST_info_index.pickle", "wb") as handle:
            pickle.dump(ENST_info_splitfile, handle, protocol=pickle.HIGHEST_PROTOCOL)

        # parse codons and assign phase to each (coding) chr position
        print("-parsing codons")
        with open(f"{out_dir}/debug.txt", "w") as wfh:
            for ENST_ID in ENST_info.keys():
                my_transcript = ENST_info[ENST_ID]
                transcript_type = my_transcript.description.split("|")[1]
                if transcript_type == "protein_coding":

                    # constructing the list of cds
                    cdsList = [
                        feat for feat in my_transcript.features if feat.type == "CDS"
                    ]
                    CDS_first = cdsList[0]
                    CDS_last = cdsList[len(cdsList) - 1]
                    #
                    strand = CDS_first.strand
                    Chr = CDS_first.location.ref
                    if not Chr in ENST_codons_dict.keys():
                        ENST_codons_dict[Chr] = dict()
                    # parse codons
                    if strand == 1:
                        start_phase = 1
                        count = 1
                        for cds in cdsList:
                            start = cds.location.start + 0
                            end = cds.location.end + 0
                            codon_assignment, end_phase = assign_codon_position(
                                start=start, end=end, start_phase=start_phase
                            )
                            start_phase = end_phase + 1
                            if start_phase == 4:
                                start_phase = 1
                            for key, val in codon_assignment.items():
                                if not key in ENST_codons_dict[Chr].keys():
                                    ENST_codons_dict[Chr][key] = dict()
                                ENST_codons_dict[Chr][key][ENST_ID] = val
                            if count == 1:  # debug
                                wfh.write(
                                    f"{ENST_ID} +1\nstrand first cds:{codon_assignment}\n"
                                )
                            if count == len(cdsList):
                                wfh.write(f"last cds:{codon_assignment}\n")
                            count += 1

                    else:  # -1 strand
                        start_phase = -1
                        count = 1
                        for cds in reversed(cdsList):
                            start = cds.location.start + 0
                            end = cds.location.end + 0
                            codon_assignment, end_phase = assign_codon_position(
                                start=end, end=start, start_phase=start_phase
                            )  # start needs to be the coordinate of the first nt in the first codon
                            start_phase = end_phase - 1
                            if start_phase == -4:
                                start_phase = -1
                            for key, val in codon_assignment.items():
                                if not key in ENST_codons_dict[Chr].keys():
                                    ENST_codons_dict[Chr][key] = dict()
                                ENST_codons_dict[Chr][key][ENST_ID] = val
                            if count == 1:  # debug
                                wfh.write(
                                    f"{ENST_ID} -1 strand\nfirst cds:{codon_assignment}\n"
                                )
                            if count == len(cdsList):
                                wfh.write(f"last cds:{codon_assignment}\n")
                            count += 1

        # write dict to file
        with open(f"{out_dir}/ENST_codonPhase.pickle", "wb") as handle:
            pickle.dump(ENST_codons_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

        # split the dict and write to pickle files
        ENST_to_codoPhase_splitfile = (
            {}
        )  # indexing dictionary, the values of this dict is a list of all the files that may contain the ENST codon phase info
        codonPhase_outdir = os.path.join(out_dir, "ENST_codonPhases")
        if os.path.exists(codonPhase_outdir):
            shutil.rmtree(codonPhase_outdir)
        os.makedirs(codonPhase_outdir)
        # loop through all chrs
        part = 0
        for Chr in ENST_codons_dict.keys():
            # initialization
            tmp_file_ENST_list = []
            temp_dict = {}
            for pos in ENST_codons_dict[Chr].keys():  # loop through all pos
                # populate tmp dict
                update_dict(
                    temp_dict, key=Chr, key2=pos, value=ENST_codons_dict[Chr][pos]
                )
                # update indexing dictionary
                for ID in ENST_codons_dict[Chr][pos].keys():
                    if ID in ENST_to_codoPhase_splitfile:
                        if not part in ENST_to_codoPhase_splitfile[ID]:
                            ENST_to_codoPhase_splitfile[ID].append(part)
                    else:
                        ENST_to_codoPhase_splitfile[ID] = [part]
                    if not ID in tmp_file_ENST_list:
                        tmp_file_ENST_list.append(ID)
                if len(tmp_file_ENST_list) >= 200:
                    # write tmp dict to file
                    with open(
                        f"{codonPhase_outdir}/ENST_codonPhase_part{part}.pickle", "wb"
                    ) as handle:
                        pickle.dump(temp_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
                    temp_dict = {}  # reset temp dict
                    part += 1
                    tmp_file_ENST_list = []  # reset temp ENST list
            # write last chunk of the current chromosome
            if len(tmp_file_ENST_list) > 0:
                with open(
                    f"{codonPhase_outdir}/ENST_codonPhase_part{part}.pickle", "wb"
                ) as handle:
                    pickle.dump(temp_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
                    part += 1
        # write indexing dictionary
        with open(f"{codonPhase_outdir}/ENST_codonPhase_index.pickle", "wb") as handle:
            pickle.dump(
                ENST_to_codoPhase_splitfile, handle, protocol=pickle.HIGHEST_PROTOCOL
            )

        # populate loc2posType dict
        print("-parsing location types(e.g. exon/intron junctions)")
        for ENST_ID in ENST_info.keys():
            UTR5p, UTR3p = get_UTR_loc(ENST_ID, ENST_info)  # get UTR loc
            cds_loc = get_cds_loc(ENST_ID, ENST_info)  # get cds loc
            exon_loc = get_exon_loc(ENST_ID, ENST_info)  # get cds loc
            log.debug(f"{ENST_info[ENST_ID].description}")
            log.debug(f"5UTR {UTR5p}")
            log.debug(f"3UTR {UTR3p}")
            log.debug(f"=======================")

            # mark UTRs
            for loc in UTR5p:
                loc2posType = update_dictOfDict(
                    mydict=loc2posType,
                    key=loc[2],
                    key2=ENST_ID,
                    key3=tuple([loc[0], loc[1]]),
                    value="5UTR",
                )
            for loc in UTR3p:
                loc2posType = update_dictOfDict(
                    mydict=loc2posType,
                    key=loc[2],
                    key2=ENST_ID,
                    key3=tuple([loc[0], loc[1]]),
                    value="3UTR",
                )

            # mark cds
            for idx, loc in enumerate(cds_loc):
                chr, start, end, strand = loc
                # mark cds
                loc2posType = update_dictOfDict(
                    mydict=loc2posType,
                    key=chr,
                    key2=ENST_ID,
                    key3=tuple([start + 0, end + 0]),
                    value="cds",
                )
            # mark exon/intron junctions
            for idx, loc in enumerate(exon_loc):
                chr, start, end, strand = loc
                if idx != 0 and idx != (len(exon_loc) - 1):  # not first or last cds
                    if strand == 1 or strand == "1" or strand == "+":  # pos strand
                        # exon-intron
                        loc2posType = update_dictOfDict(
                            mydict=loc2posType,
                            key=chr,
                            key2=ENST_ID,
                            key3=tuple([end - 1, end + 2]),
                            value="within_2bp_of_exon_intron_junction",
                        )
                        loc2posType = update_dictOfDict(
                            mydict=loc2posType,
                            key=chr,
                            key2=ENST_ID,
                            key3=tuple([end - 2, end + 3]),
                            value="within_3bp_of_exon_intron_junction",
                        )  # for recoding off limit (old)
                        loc2posType = update_dictOfDict(
                            mydict=loc2posType,
                            key=chr,
                            key2=ENST_ID,
                            key3=tuple([end - 3, end - 2]),
                            value="3N4bp_up_of_exon_intron_junction",
                        )
                        loc2posType = update_dictOfDict(
                            mydict=loc2posType,
                            key=chr,
                            key2=ENST_ID,
                            key3=tuple([end + 3, end + 6]),
                            value="3_to_6bp_down_of_exon_intron_junction",
                        )
                        loc2posType = update_dictOfDict(
                            mydict=loc2posType,
                            key=chr,
                            key2=ENST_ID,
                            key3=tuple([end - 2, end + 6]),
                            value="-3_to_+6bp_of_exon_intron_junction",
                        )  # for recoding off limit (new)
                        # intron-exon
                        loc2posType = update_dictOfDict(
                            mydict=loc2posType,
                            key=chr,
                            key2=ENST_ID,
                            key3=tuple([start - 2, start + 1]),
                            value="within_2bp_of_intron_exon_junction",
                        )
                        loc2posType = update_dictOfDict(
                            mydict=loc2posType,
                            key=chr,
                            key2=ENST_ID,
                            key3=tuple([start - 3, start + 2]),
                            value="within_3bp_of_intron_exon_junction",
                        )  # for recoding off limit (old)
                        loc2posType = update_dictOfDict(
                            mydict=loc2posType,
                            key=chr,
                            key2=ENST_ID,
                            key3=tuple([start - 4, start - 3]),
                            value="3N4bp_up_of_intron_exon_junction",
                        )
                        loc2posType = update_dictOfDict(
                            mydict=loc2posType,
                            key=chr,
                            key2=ENST_ID,
                            key3=tuple([start + 2, start + 3]),
                            value="3N4bp_down_of_intron_exon_junction",
                        )
                        loc2posType = update_dictOfDict(
                            mydict=loc2posType,
                            key=chr,
                            key2=ENST_ID,
                            key3=tuple([start - 3, start + 1]),
                            value="-3_to_+2bp_of_intron_exon_junction",
                        )  # for recoding off limit (new)
                    else:  # neg strand
                        # exon-intron
                        loc2posType = update_dictOfDict(
                            mydict=loc2posType,
                            key=chr,
                            key2=ENST_ID,
                            key3=tuple([start - 2, start + 1]),
                            value="within_2bp_of_exon_intron_junction",
                        )
                        loc2posType = update_dictOfDict(
                            mydict=loc2posType,
                            key=chr,
                            key2=ENST_ID,
                            key3=tuple([start - 3, start + 2]),
                            value="within_3bp_of_exon_intron_junction",
                        )  # for recoding off limit (old)
                        loc2posType = update_dictOfDict(
                            mydict=loc2posType,
                            key=chr,
                            key2=ENST_ID,
                            key3=tuple([start + 2, start + 3]),
                            value="3N4bp_up_of_exon_intron_junction",
                        )
                        loc2posType = update_dictOfDict(
                            mydict=loc2posType,
                            key=chr,
                            key2=ENST_ID,
                            key3=tuple([start - 6, start - 3]),
                            value="3_to_6bp_down_of_exon_intron_junction",
                        )
                        loc2posType = update_dictOfDict(
                            mydict=loc2posType,
                            key=chr,
                            key2=ENST_ID,
                            key3=tuple([start - 6, start + 2]),
                            value="-3_to_+6bp_of_exon_intron_junction",
                        )  # for recoding off limit (new)
                        # intron-exon
                        loc2posType = update_dictOfDict(
                            mydict=loc2posType,
                            key=chr,
                            key2=ENST_ID,
                            key3=tuple([end - 1, end + 2]),
                            value="within_2bp_of_intron_exon_junction",
                        )
                        loc2posType = update_dictOfDict(
                            mydict=loc2posType,
                            key=chr,
                            key2=ENST_ID,
                            key3=tuple([end - 2, end + 3]),
                            value="within_3bp_of_intron_exon_junction",
                        )  # for recoding off limit (old)
                        loc2posType = update_dictOfDict(
                            mydict=loc2posType,
                            key=chr,
                            key2=ENST_ID,
                            key3=tuple([end + 3, end + 4]),
                            value="3N4bp_up_of_intron_exon_junction",
                        )
                        loc2posType = update_dictOfDict(
                            mydict=loc2posType,
                            key=chr,
                            key2=ENST_ID,
                            key3=tuple([end - 3, end - 2]),
                            value="3N4bp_down_of_intron_exon_junction",
                        )
                        loc2posType = update_dictOfDict(
                            mydict=loc2posType,
                            key=chr,
                            key2=ENST_ID,
                            key3=tuple([end - 1, end + 3]),
                            value="-3_to_+2bp_of_intron_exon_junction",
                        )  # for recoding off limit (new)
                elif idx == 0 and idx != (
                    len(exon_loc) - 1
                ):  # first cds (not necessarily the one with the start codon)
                    if (
                        strand == 1 or strand == "1" or strand == "+"
                    ):  # pos strand, first cds
                        # exon-intron
                        loc2posType = update_dictOfDict(
                            mydict=loc2posType,
                            key=chr,
                            key2=ENST_ID,
                            key3=tuple([end - 1, end + 2]),
                            value="within_2bp_of_exon_intron_junction",
                        )
                        loc2posType = update_dictOfDict(
                            mydict=loc2posType,
                            key=chr,
                            key2=ENST_ID,
                            key3=tuple([end - 2, end + 3]),
                            value="within_3bp_of_exon_intron_junction",
                        )  # for recoding off limit (old)
                        loc2posType = update_dictOfDict(
                            mydict=loc2posType,
                            key=chr,
                            key2=ENST_ID,
                            key3=tuple([end - 3, end - 2]),
                            value="3N4bp_up_of_exon_intron_junction",
                        )
                        loc2posType = update_dictOfDict(
                            mydict=loc2posType,
                            key=chr,
                            key2=ENST_ID,
                            key3=tuple([end + 3, end + 6]),
                            value="3_to_6bp_down_of_exon_intron_junction",
                        )
                        loc2posType = update_dictOfDict(
                            mydict=loc2posType,
                            key=chr,
                            key2=ENST_ID,
                            key3=tuple([end - 2, end + 6]),
                            value="-3_to_+6bp_of_exon_intron_junction",
                        )  # for recoding off limit (new)
                    else:  # neg strand last cds
                        # intron-exon
                        loc2posType = update_dictOfDict(
                            mydict=loc2posType,
                            key=chr,
                            key2=ENST_ID,
                            key3=tuple([end - 1, end + 2]),
                            value="within_2bp_of_intron_exon_junction",
                        )
                        loc2posType = update_dictOfDict(
                            mydict=loc2posType,
                            key=chr,
                            key2=ENST_ID,
                            key3=tuple([end - 2, end + 3]),
                            value="within_3bp_of_intron_exon_junction",
                        )  # for recoding off limit (old)
                        loc2posType = update_dictOfDict(
                            mydict=loc2posType,
                            key=chr,
                            key2=ENST_ID,
                            key3=tuple([end + 3, end + 4]),
                            value="3N4bp_up_of_intron_exon_junction",
                        )
                        loc2posType = update_dictOfDict(
                            mydict=loc2posType,
                            key=chr,
                            key2=ENST_ID,
                            key3=tuple([end - 3, end - 2]),
                            value="3N4bp_down_of_intron_exon_junction",
                        )
                        loc2posType = update_dictOfDict(
                            mydict=loc2posType,
                            key=chr,
                            key2=ENST_ID,
                            key3=tuple([end - 1, end + 3]),
                            value="-3_to_+2bp_of_intron_exon_junction",
                        )  # for recoding off limit (new)
                elif (
                    idx == (len(exon_loc) - 1) and idx != 0
                ):  # last cds (not necessarily the one with the stop codon)
                    if (
                        strand == 1 or strand == "1" or strand == "+"
                    ):  # pos strand, last cds
                        # intron-exon
                        loc2posType = update_dictOfDict(
                            mydict=loc2posType,
                            key=chr,
                            key2=ENST_ID,
                            key3=tuple([start - 2, start + 1]),
                            value="within_2bp_of_intron_exon_junction",
                        )
                        loc2posType = update_dictOfDict(
                            mydict=loc2posType,
                            key=chr,
                            key2=ENST_ID,
                            key3=tuple([start - 3, start + 2]),
                            value="within_3bp_of_intron_exon_junction",
                        )  # for recoding off limit (old)
                        loc2posType = update_dictOfDict(
                            mydict=loc2posType,
                            key=chr,
                            key2=ENST_ID,
                            key3=tuple([start - 4, start - 3]),
                            value="3N4bp_up_of_intron_exon_junction",
                        )
                        loc2posType = update_dictOfDict(
                            mydict=loc2posType,
                            key=chr,
                            key2=ENST_ID,
                            key3=tuple([start + 2, start + 3]),
                            value="3N4bp_down_of_intron_exon_junction",
                        )
                        loc2posType = update_dictOfDict(
                            mydict=loc2posType,
                            key=chr,
                            key2=ENST_ID,
                            key3=tuple([start - 3, start + 1]),
                            value="-3_to_+2bp_of_intron_exon_junction",
                        )  # for recoding off limit (new)
                    else:  # neg strand, first cds
                        # exon-intron
                        loc2posType = update_dictOfDict(
                            mydict=loc2posType,
                            key=chr,
                            key2=ENST_ID,
                            key3=tuple([start - 2, start + 1]),
                            value="within_2bp_of_exon_intron_junction",
                        )
                        loc2posType = update_dictOfDict(
                            mydict=loc2posType,
                            key=chr,
                            key2=ENST_ID,
                            key3=tuple([start - 3, start + 2]),
                            value="within_3bp_of_exon_intron_junction",
                        )  # for recoding off limit (old)
                        loc2posType = update_dictOfDict(
                            mydict=loc2posType,
                            key=chr,
                            key2=ENST_ID,
                            key3=tuple([start + 2, start + 3]),
                            value="3N4bp_up_of_exon_intron_junction",
                        )
                        loc2posType = update_dictOfDict(
                            mydict=loc2posType,
                            key=chr,
                            key2=ENST_ID,
                            key3=tuple([start - 6, start - 3]),
                            value="3_to_6bp_down_of_exon_intron_junction",
                        )
                        loc2posType = update_dictOfDict(
                            mydict=loc2posType,
                            key=chr,
                            key2=ENST_ID,
                            key3=tuple([start - 6, start + 2]),
                            value="-3_to_+6bp_of_exon_intron_junction",
                        )  # for recoding off limit (new)

        # write dict to file
        with open(f"{out_dir}/loc2posType.pickle", "wb") as handle:
            pickle.dump(loc2posType, handle, protocol=pickle.HIGHEST_PROTOCOL)

        endtime = datetime.datetime.now()
        elapsed_sec = endtime - starttime
        elapsed_min = elapsed_sec.seconds / 60
        log.info(
            f"finished in {elapsed_min:.2f} min, processed {line_count} lines {file}"
        )
        print(
            f"finished in {elapsed_min:.2f} min, processed {line_count} lines {file}",
            flush=True,
        )

    except Exception as e:
        print("Unexpected error:", str(sys.exc_info()))
        traceback.print_exc()
        print("additional information:", e)
        PrintException()


##########################
## function definitions ##
##########################
def make_output_dir(path):
    if not os.path.exists(path):  # output doesn't exist
        os.makedirs(path)
        return True
    else:  # output dir exists
        inp = input(f"{path} already exists...\npress y to confirm overwritting")
        if inp == "y" or inp == "Y":
            shutil.rmtree(path, ignore_errors=True, onerror=onerror)
            os.makedirs(path, exist_ok=True)
            return True
        else:
            return False


def onerror(func, path, exc_info):
    """
    Error handler for ``shutil.rmtree``.

    If the error is due to an access error (read only file)
    it attempts to add write permission and then retries.

    If the error is for another reason it re-raises the error.

    Usage : ``shutil.rmtree(path, onerror=onerror)``
    """
    import stat

    # Is the error an access error?
    if not os.access(path, os.W_OK):
        os.chmod(path, stat.S_IWUSR)
        func(path)
    else:
        raise


def crc32(fileName):
    prev = 0
    for eachLine in open(fileName, "rb"):
        prev = zlib.crc32(eachLine, prev)
    return "%X" % (prev & 0xFFFFFFFF)


def assign_codon_position(start, end, start_phase):
    """
    phase = 1 , 2 or 3
    >>>assign_codon_position(1,10,1)
    ({1: 1, 2: 2, 3: 3, 4: 1, 5: 2, 6: 3, 7: 1, 8: 2, 9: 3, 10: 1}, 1)
    >>>assign_codon_position(1,10,2)
    ({1: 2, 2: 3, 3: 1, 4: 2, 5: 3, 6: 1, 7: 2, 8: 3, 9: 1, 10: 2}, 2)
    >>>assign_codon_position(1,10,3)
    ({1: 3, 2: 1, 3: 2, 4: 3, 5: 1, 6: 2, 7: 3, 8: 1, 9: 2, 10: 3}, 3)
     phase = -1 , -2 or -3
    >>>assign_codon_position(10,1,-1)
    ({10: -1, 9: -2, 8: -3, 7: -1, 6: -2, 5: -3, 4: -1, 3: -2, 2: -3, 1: -1}, -1)
    >>>assign_codon_position(10,1,-2)
    ({10: -2, 9: -3, 8: -1, 7: -2, 6: -3, 5: -1, 4: -2, 3: -3, 2: -1, 1: -2}, -2)
    >>>assign_codon_position(10,1,-3)
    ({10: -3, 9: -1, 8: -2, 7: -3, 6: -1, 5: -2, 4: -3, 3: -1, 2: -2, 1: -3}, -3)
    """
    if start <= end and start_phase > 0:  # positive strand
        length = end - start + 1
        codon_pos = list(
            itertools.islice(
                itertools.cycle([1, 2, 3]), start_phase - 1, length + start_phase - 1
            )
        )
        mydict = dict(zip(list(range(start, end + 1)), codon_pos))
        end_phase = codon_pos[-1]
        return mydict, end_phase
    elif start >= end and start_phase < 0:  # negative strand
        length = start - end + 1
        codon_pos = list(
            itertools.islice(
                itertools.cycle([-1, -2, -3]),
                (0 - start_phase) - 1,
                length + (0 - start_phase) - 1,
            )
        )
        mydict = dict(zip(list(range(start, end - 1, -1)), codon_pos))
        end_phase = codon_pos[-1]
        return mydict, end_phase
    else:
        sys.exit(
            "ERROR in assign_codon_position(start, end, start_phase), if start_phase>0, start must < end, vice versa"
        )


def update_dict(mydict, key, key2, value):
    """
    make the following update:  mydict[key][key2] = value
    """
    if not key in mydict.keys():
        mydict[key] = dict()
    mydict[key][key2] = value

    return mydict


def update_dictOfDict(mydict, key, key2, key3, value):
    """
    make the following update:  mydict[key][key2][key3] = value
    """
    if not key in mydict.keys():
        mydict[key] = dict()
    if not key2 in mydict[key].keys():
        mydict[key][key2] = dict()
    mydict[key][key2][key3] = value
    return mydict


def get_UTR_loc(ENST_ID, ENST_info):
    """
    get the UTR locations of a transcript
    return [UTR5p,UTR3p]  #UTR5p:[[start,end],...]
    """
    locList = []
    my_transcript = ENST_info[ENST_ID]  # get the seq record
    # construct the list of cds
    cdsList = [feat for feat in my_transcript.features if feat.type == "CDS"]
    # construct the list of exon
    exonList = [feat for feat in my_transcript.features if feat.type == "exon"]
    UTR5p = []
    UTR3p = []

    res = dict()  # {[start,end]:UTR_type} # temp storage of UTR results
    status = 0  # 0 1 #status indicator, this parameter signals the end of UTR (depending on the strand, it may be the end of 5p or 3p)
    status_dict = {
        1: {0: "5UTR", 1: "3UTR"},  # status_dict[strand][status]
        -1: {0: "3UTR", 1: "5UTR"},
    }

    for idx, exon in enumerate(exonList):  # go through exons
        superimpose_cds = find_superimpose_cds(
            exon, cdsList
        )  # find superimposed cds with the current exon
        overlap_cds = find_overlap_cds(
            exon, cdsList
        )  # find overlapping cds with the current exon
        strand = exon.location.strand
        chr = exon.location.ref
        if len(superimpose_cds) == 0 and len(overlap_cds) == 0:  # untranslated exon
            UTR_type = status_dict[strand][status]  # get the type of UTR, 5p or 3p
            interval = [exon.location.start + 0, exon.location.end + 0, chr]
            res[tuple(interval)] = UTR_type  # save UTR
            log.debug(
                f"{ENST_ID} exon_idx={idx} strand={strand} untranslated exon UTR_type={UTR_type}"
            )
        elif (
            len(superimpose_cds) == 0 and len(overlap_cds) != 0
        ):  # partial-transcribed exon
            UTR_type = status_dict[strand][status]  # get the type of UTR, 5p or 3p
            cds = cdsList[overlap_cds[0]]
            interval = get_nonoverlap(
                [exon.location.start + 0, exon.location.end + 0],
                [cds.location.start + 0, cds.location.end + 0],
            )  # will return None if both 5UTR and 3UTR are in the same exon
            # check if the partially untranslated part is 5UTR or 3UTR
            # only one UTR is in the exon
            if not interval is None:
                if strand == 1 or strand == "1" or strand == "+":
                    if exon.location.start == cds.location.start:
                        UTR_type = "3UTR"
                    else:
                        UTR_type = "5UTR"
                else:  # neg strand
                    if exon.location.start == cds.location.start:
                        UTR_type = "5UTR"
                    else:
                        UTR_type = "3UTR"
                interval.append(chr)
                res[tuple(interval)] = UTR_type  # save UTR
                log.debug(
                    f"{ENST_ID} exon_idx={idx} strand={strand} partial-translated exon UTR_type={UTR_type}"
                )
            else:  # both 5UTR and 3UTR are in the same exon
                interval_1 = [exon.location.start + 0, cds.location.start - 1, chr]
                res[tuple(interval_1)] = UTR_type  # save UTR
                status = 1  # maually update UTR type
                UTR_type = status_dict[strand][status]  # maually update UTR type
                interval_2 = [cds.location.end + 1, exon.location.end + 0, chr]
                res[tuple(interval_2)] = UTR_type  # save UTR
                log.debug(
                    f"{ENST_ID} exon_idx={idx} strand={strand} partial-translated exon UTR_type=5UTR+3UTR"
                )
            status = 1  # this parameter signals the end of UTR (depending on the strand, it may be the end of 5p or 3p)
        elif len(superimpose_cds) != 0:  # fully-transcribed exon
            log.debug(f"{ENST_ID} exon_idx={idx} strand={strand} fully-translated exon")
            status = 1  # this parameter signals the end of UTR (depending on the strand, it may be the end of 5p or 3p)

    for key, val in res.items():
        if val == "5UTR":
            UTR5p.append(key)
        if val == "3UTR":
            UTR3p.append(key)
    return [UTR5p, UTR3p]


def find_overlap_cds(exon, cdsList):
    """
    return a list of index of cds that overlap with input exon
    """
    overlap_cds = []
    for idx, cds in enumerate(cdsList):
        if (
            cds.location.start <= exon.location.start <= cds.location.end
            or cds.location.start <= exon.location.end <= cds.location.end
            or exon.location.start <= cds.location.start <= exon.location.end
            or exon.location.start <= cds.location.end <= exon.location.end
        ):
            overlap_cds.append(idx)
    return overlap_cds


def find_superimpose_cds(exon, cdsList):
    """
    return a list of index of cds that superimpose with input exon
    """
    superimpose_cds = []
    for idx, cds in enumerate(cdsList):
        if (
            cds.location.start == exon.location.start
            and cds.location.end == exon.location.end
        ):
            superimpose_cds.append(idx)
    return superimpose_cds


def get_nonoverlap(loc1, loc2):
    """
    compare two locations [start,end] loc2 [start,end] , one loc encapsulates the other
    return the non-overlapping interval
    """
    if loc1[0] == loc2[0]:
        if loc1[1] != loc2[1]:
            return [min([loc1[1], loc2[1]]) + 1, max([loc1[1], loc2[1]])]
    elif loc1[1] == loc2[1]:
        if loc1[0] != loc2[0]:
            return [min([loc1[0], loc2[0]]), max([loc1[0], loc2[0]]) - 1]
    else:
        return None


def get_cds_loc(ENST_ID, ENST_info):
    """
    Get the chromosomal location of CDS
    input: ENST_ID, ENST_info
    output: a list loc [chr,start,end,strand] #start < end
    """
    locList = []
    my_transcript = ENST_info[ENST_ID]  # get the seq record
    # constructing the list of cds
    cdsList = [feat for feat in my_transcript.features if feat.type == "CDS"]
    for cds in cdsList:
        locList.append(
            [cds.location.ref, cds.location.start, cds.location.end, cds.strand]
        )
    return locList


def get_exon_loc(ENST_ID, ENST_info):
    """
    Get the chromosomal location of exon
    input: ENST_ID, ENST_info
    output: a list loc [chr,start,end,strand] #start < end
    """
    locList = []
    my_transcript = ENST_info[ENST_ID]  # get the seq record
    # constructing the list of cds
    exonList = [feat for feat in my_transcript.features if feat.type == "exon"]
    for exon in exonList:
        locList.append(
            [exon.location.ref, exon.location.start, exon.location.end, exon.strand]
        )
    return locList


def get_start_stop_loc(ENST_ID, ENST_info):
    """
    Get the chromosomal location of start and stop codons
    input: ENST_ID, ENST_info
    output: a list of two items
            ATG_loc: [chr,start,end,strand]  #start < end
            stop_loc: [chr,start,end,strand] #start < end
    """
    my_transcript = ENST_info[ENST_ID]  # get the seq record
    # constructing the list of cds
    cdsList = [feat for feat in my_transcript.features if feat.type == "CDS"]
    if len(cdsList) == 0:
        return [None, None]
    CDS_first = cdsList[0]
    CDS_last = cdsList[len(cdsList) - 1]
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


def PrintException():
    exc_type, exc_obj, tb = sys.exc_info()
    f = tb.tb_frame
    lineno = tb.tb_lineno
    filename = f.f_code.co_filename
    linecache.checkcache(filename)
    line = linecache.getline(filename, lineno, f.f_globals)
    print(
        'EXCEPTION IN ({}, LINE {} "{}"): {}'.format(
            filename, lineno, line.strip(), exc_obj
        )
    )


def plus_minus_strand_to_numeric(input):
    if input == "+":
        return 1
    elif input == "-":
        return -1
    else:
        return 0


if __name__ == "__main__":
    main()
