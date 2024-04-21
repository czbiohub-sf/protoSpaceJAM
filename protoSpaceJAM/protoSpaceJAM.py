import datetime
import linecache
import logging
import os.path
import pickle
import sys
import traceback
import time
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.GenBank import Record


from protoSpaceJAM.util.utils import MyParser, ColoredLogger, read_pickle_files, cal_elapsed_time, get_gRNAs,get_gRNAs_target_coordinate, \
    get_HDR_template #uncomment this for pip installation

# from util.utils import MyParser, ColoredLogger, read_pickle_files, cal_elapsed_time, get_gRNAs, get_gRNAs_target_coordinate, \
#     get_HDR_template

def parse_args(test_mode=False):
    parser = MyParser(description="protoSpaceJAM: perfectionist CRISPR knock-in design at scale\n")
    IO = parser.add_argument_group('input/output')
    IO.add_argument(
        "--path2csv",
        default=os.path.join("input","test_input.csv"),
        type=str,
        help="Path to a csv file containing the input knock-in sites, see input/test_input.csv for an example\n *required columns*: 'Ensembl_ID' (specifying the transcript ID), and either 'Target_terminus' or 'Chromosome','Coordinate' (specifying the terminus of the transcript or a genomic coordinate in the transcript)",
        metavar="<PATH_TO_CSV>",
    )
    IO.add_argument(
        "--outdir",
        default=os.path.join("output","test"),
        type=str,
        metavar = "<PATH_TO_OUTPUT_DIRECTORY>",
        help="Path to the output directory"
    )
    genome = parser.add_argument_group('genome')
    genome.add_argument(
        "--genome_ver",
        default="GRCh38",
        type=str,
        help="Genome and version to use, possible values are GRCh38, GRCm39, and GRCz11",
        metavar="<string>",
    )
    gRNA = parser.add_argument_group('gRNA')
    gRNA.add_argument(
        "--pam",
        default="NGG",
        type=str,
        help="PAM sequence (default: NGG)",
        metavar="<string>",
    )
    gRNA.add_argument(
        "--num_gRNA_per_design",
        default=1,
        type=int,
        help="Number of gRNAs to return per site (default: 1)",
        metavar="<integer>",
    )
    gRNA.add_argument(
        "--no_regulatory_penalty",
        default=False,
        action="store_true",
        help="Turn off penalty for gRNAs cutting in UTRs or near splice junctions, default: penalty on",
    )
    gRNA.add_argument(
        "--alpha1",
        default=1.0,
        type=float,
        help="raise the specificity weight to the power of this number, default: 1.0, range: [0,1]",
        metavar="<float>",
    )
    gRNA.add_argument(
        "--alpha2",
        default=1.0,
        type=float,
        help="raise the insert dist. weight to the power of this number, default: 1.0, range: [0,1]",
        metavar="<float>",
    )
    gRNA.add_argument(
        "--alpha3",
        default=1.0,
        type=float,
        help="raise the position weight to the power of this number, default: 1.0, range: [0,1]",
        metavar="<float>",
    )
    payload = parser.add_argument_group('payload')
    payload.add_argument(
        "--payload",
        default="",
        type=str,
        help="Define the payload sequence for every site, regardless of terminus or coordinates, overrides --Npayload, --Cpayload, POSpayload, --Tag, --Linker",
        metavar="<string>",
    )
    payload.add_argument(
        "--Npayload",
        default="ACCGAGCTCAACTTCAAGGAGTGGCAAAAGGCCTTTACCGATATGATGGGTGGCGGATTGGAAGTTTTGTTTCAAGGTCCAGGAAGTGGT",
        type=str,
        help="Payload sequence to use at the N terminus (default: mNG11 + XTEN80): ACCGAGCTCAACTTCAAGGAGTGGCAAAAGGCCTTTACCGATATGATGGGTGGCGGATTGGAAGTTTTGTTTCAAGGTCCAGGAAGTGGT, overrides --Tag and --Linker",
        metavar="<string>",
    )
    payload.add_argument(
        "--Cpayload",
        default="GGTGGCGGATTGGAAGTTTTGTTTCAAGGTCCAGGAAGTGGTACCGAGCTCAACTTCAAGGAGTGGCAAAAGGCCTTTACCGATATGATG",
        type=str,
        help="Payload sequence to use at the C terminus (default: XTEN80 + mNG11): GGTGGCGGATTGGAAGTTTTGTTTCAAGGTCCAGGAAGTGGTACCGAGCTCAACTTCAAGGAGTGGCAAAAGGCCTTTACCGATATGATG, overrides --Tag and --Linker",
        metavar="<string>",
    )
    payload.add_argument(
        "--POSpayload",
        default="GGTGGCGGATTGGAAGTTTTGTTTCAAGGTCCAGGAAGTGGTACCGAGCTCAACTTCAAGGAGTGGCAAAAGGCCTTTACCGATATGATG",
        type=str,
        help="Payload sequence to use at the specific genomic coordinates (default: XTEN80 + mNG11): GGTGGCGGATTGGAAGTTTTGTTTCAAGGTCCAGGAAGTGGTACCGAGCTCAACTTCAAGGAGTGGCAAAAGGCCTTTACCGATATGATG, overrides --Tag and --Linker",
        metavar="<string>",
    )
    payload.add_argument(
        "--Tag",
        default="ACCGAGCTCAACTTCAAGGAGTGGCAAAAGGCCTTTACCGATATGATG",
        type=str,
        help="default is the mNG11 tag",
        metavar="<string>",
    )
    payload.add_argument(
        "--Linker",
        default="GGTGGCGGATTGGAAGTTTTGTTTCAAGGTCCAGGAAGTGGT",
        type=str,
        help="default is the XTEN80 linker",
        metavar="<string>",
    )

    donor = parser.add_argument_group('donor')
    donor.add_argument(
        "--Donor_type",
        default="ssODN",
        help="Set the type of donor, possible values are ssODN and dsDNA (default: ssODN)",
        type=str,
        metavar="<string>",
    )
    donor.add_argument(
        "--HA_len",
        default=500,
        help="[dsDNA] Length of the desired homology arm on each side (default: 500)",
        type=int,
        metavar="<integer>",
    )
    donor.add_argument(
        "--Strand_choice",
        default="auto",
        help="[ssODN] Strand choice of ssoODN, Possible values are 'auto', 'TargetStrand', 'NonTargetStrand', 'CodingStrand' and 'NonCodingStrand'",
        type=str,
        metavar="<string>",
    )
    donor.add_argument(
        "--ssODN_max_size",
        type=int,
        default=200,
        help="Enforce a length restraint on the the ssODN donor (default: 200), The ssODN donor will be centered on the payload and the recoded region",
        metavar="<int>",
    )
    donor.add_argument(
        "--CheckEnzymes",
        default="",
        help="[dsDNA] Name of Restriction digestion enzymes, separated by '|', to flag and trim, for example BsaI|EcoRI (default: None)",
        type=str,
        metavar="<string>",
    )
    donor.add_argument(
        "--CustomSeq2Avoid",
        default="",
        help="[dsDNA] Custom sequences, separated by '|', to flag and trim (default: None)",
        type=str,
        metavar="<string>",
    )
    donor.add_argument(
        "--MinArmLenPostTrim",
        default=0,
        help="[dsDNA] Minimum length of the homology arm after trimming. Set to 0 to turn off trimming (default: 0)",
        type=int,
        metavar="<integer>",
    )
    recoding = parser.add_argument_group('recoding')
    recoding.add_argument(
        "--recoding_off",
        default=False,
        action="store_true",
        help="Turn off *all* recoding",
    )
    recoding.add_argument(
        "--recoding_stop_recut_only",
        default=False,
        action="store_true",
        help="Recode the gRNA recognition site to prevent recut",
    )
    recoding.add_argument(
        "--recoding_full",
        default=False,
        action="store_true",
        help="Use full recoding: recode both the gRNA recognition site and the cut-to-insert region (default: on)",
    )
    recoding.add_argument(
        "--recoding_coding_region_only",
        default=False,
        action="store_true",
        help="Only recode the coding region",
    )
    recoding.add_argument(
        "--cfdThres",
        default=0.03,
        help="Threshold that protoSpaceJAM will attempt to lower the recut potential (measured by the CFD score) to (default: 0.03)",
        metavar="<float>",
    )
    recoding.add_argument(
        "--recode_order",
        default="PAM_first",
        help="Prioritize recoding in the PAM or in protospacer, possible values: protospacer_first, PAM_first (default: PAM_first)",
        metavar="<string>",
    )
    misc = parser.add_argument_group('misc.')
    misc.add_argument(
        "--test_mode",
        default=False,
        help="used by the unit tests, not user-oriented",
        metavar="<boolean>",
    )
    config = parser.parse_args()
    return config, parser


def main(custom_args=None):
    """
    main function
    custom_args: a dict of arguments to override the default arguments
    """
    try:
        #set up working directory
        if not os.path.exists(os.path.join("precomputed_gRNAs")):
            if not os.path.exists(os.path.join("protoSpaceJAM", "precomputed_gRNAs")):
                sys.exit("precomputed_gRNAs folder not found, please run the script from the repo's root directory")
            else:
                os.chdir("protoSpaceJAM")

        logging.setLoggerClass(ColoredLogger)
        # logging.basicConfig()
        log = logging.getLogger("protoSpaceJAM")
        log.propagate = False
        log.setLevel(logging.INFO)  # set the level of warning displayed
        # log.setLevel(logging.DEBUG) #set the level of warning displayed

        # configs
        config = vars(parse_args()[0])
        parser = parse_args()[1]

        # apply custom args
        if not custom_args is None and len(custom_args) > 0:
            for c_arg in custom_args:
                config[c_arg] = custom_args[c_arg]

        # Exit if no arguments provided and not in test mode
        if len(sys.argv)==1 and config["test_mode"] == False:
            print("[Message] Pleases provide the following arguments: --path2csv --outdir")
            print("[Message] To run a quick example: protoSpaceJAM --path2csv input/test_input.csv --outdir output/test\n")

            parser.print_help(sys.stderr)
            sys.exit(1)

        gRNA_num_out = config["num_gRNA_per_design"]
        max_cut2ins_dist = 50  # deprecated?
        HDR_arm_len = config["HA_len"]
        ssODN_max_size = config["ssODN_max_size"]
        spec_score_flavor = "guideMITScore"
        outdir = config["outdir"]
        reg_penalty = not config["no_regulatory_penalty"]
        syn_check_args = {
            "check_enzymes": config["CheckEnzymes"],
            "CustomSeq2Avoid": config["CustomSeq2Avoid"],
            "MinArmLenPostTrim": config["MinArmLenPostTrim"],
        }  # dictionary for multiple synthesis check arguments

        # check recoding args
        assert (
            config["recode_order"] == "protospacer_first"
            or config["recode_order"] == "PAM_first"
        )
        if config["recoding_full"] and any(
            [config["recoding_off"], config["recoding_stop_recut_only"]]
        ):
            sys.exit(
                f"Found conflicts in recoding arguments: --recoding_full cannot be used with --recoding_off or --recoding_stop_recut_only\nplease correct the issue and try again"
            )
        if config["recoding_off"] and config["recoding_stop_recut_only"]:
            sys.exit(
                f"Found conflicts in recoding arguments: --recoding_off cannot be used with --recoding_stop_recut_only\nplease correct the issue and try again"
            )

        # process recoding args
        if (not config["recoding_off"]) and (not config["recoding_stop_recut_only"]):
            config["recoding_full"] = True

        if config["recoding_off"] or config["recoding_stop_recut_only"]:
            config["recoding_full"] = False

        recoding_args = {
            "recoding_off": config["recoding_off"],
            "recoding_stop_recut_only": config["recoding_stop_recut_only"],
            "recoding_full": config["recoding_full"],
            "cfdThres": float(config["cfdThres"]),
            "recode_order": config["recode_order"],
            "recoding_coding_region_only": config["recoding_coding_region_only"],
        }

        # check donor args
        if not config["Donor_type"] in ["ssODN", "dsDNA"]:
            sys.exit(
                "Donor_type must be ssODN or dsDNA, offending value:"
                + config["Donor_type"]
                + ", please correct the issue and try again"
            )
        if not config["Strand_choice"] in [
            "auto",
            "TargetStrand",
            "NonTargetStrand",
            "CodingStrand",
            "NonCodingStrand",
        ]:
            sys.exit(
                "Strand_choice must be auto,TargetStrand,NonTargetStrand,CodingStrand or NonCodingStrand, offending value:"
                + config["Strand_choice"]
                + ", please correct the issue and try again"
            )

        # check pam
        if not config["pam"].upper() in ["NGG", "NGA", "TTTV"]:
            sys.exit("PAM must be NGG, NGA or TTTV, please correct the issue and try again")

        # parse payload
        Linker = config["Linker"]
        Tag = config["Tag"]

        if config["payload"] == "":  # no payload override
            if config["Npayload"] == "":  # no Npayload override
                config["Npayload"] = Tag + Linker
            if config["Cpayload"] == "":  # no Cpayload override
                config["Cpayload"] = Linker + Tag
        else:  # payload override
            config["Npayload"] = config["payload"]
            config["Cpayload"] = config["payload"]
            config["POSpayload"] = config["payload"]


        # check if HA_len is too short to satisfy ssODN_max_size
        if ssODN_max_size is not None and config["Donor_type"] == "ssODN":
            max_payload_size = max([len(config["Npayload"]), len(config["Cpayload"])])
            derived_HDR_arm_len = round((ssODN_max_size - max_payload_size) / 2)
            if derived_HDR_arm_len >= HDR_arm_len:
                print(
                    f"HA_len={HDR_arm_len} is to short to meet the requirement of ssODN_max_size={ssODN_max_size}, payload size={max_payload_size}\n ssODN_max_size={ssODN_max_size} requires HA_len = ssODN_max_size- max_payload_size / 2 = {derived_HDR_arm_len}"
                )
                HDR_arm_len = derived_HDR_arm_len + 100
                print(f"HA_len is adjusted to {HDR_arm_len}")

        # check alpha values
        if config["alpha1"] < 0 or config["alpha2"] < 0 or config["alpha3"] < 0:
            sys.exit("alpha values must >= 0, please correct the issue and try again")
        if config["alpha1"] >1 or config["alpha2"] >1  or config["alpha3"] >1 :
            sys.exit("alpha values must be <= 1, please correct the issue and try again")
        if config["alpha1"] == 0 and config["alpha2"] == 0 and config["alpha3"] == 0:
            sys.exit("At least one alpha value must be > 0, please correct the issue and try again")
        alphas = [config["alpha1"], config["alpha2"], config["alpha3"]] 

        # TODO: fix  the potential infinite loop
        # check memory requirement
        enough_mem = test_memory(4200)
        while not enough_mem:  # test if at least 4.2 GB memory is available
            time.sleep(5)  # retry in 5 seconds
            enough_mem = test_memory(4200)

        starttime = datetime.datetime.now()
        freq_dict = dict()

        # load gRNA info index (mapping of chromosomal location to file parts)
        log.info("loading precomputed gRNA **index to file parts**")
        loc2file_index = read_pickle_files(
            os.path.join(
                "precomputed_gRNAs",
                "gRNAs_" + config["pam"].upper(),
                "gRNA_" + config["genome_ver"],
                "gRNA.tab.gz.split.BwaMapped.scored",
                "loc2file_index.pickle",
            )
        )

        elapsed = cal_elapsed_time(starttime, datetime.datetime.now())
        log.info(f"finished loading in {elapsed[0]:.2f} min ({elapsed[1]} sec)")

        # load chr location to type (e.g. UTR, cds, exon/intron junction) mappings
        log.info(
            "loading the mapping of chromosomal location to type (e.g. UTR, cds, exon/intron junction)"
        )
        loc2posType = read_pickle_files(
            os.path.join(
                "genome_files",
                "parsed_gff3",
                config["genome_ver"],
                "loc2posType.pickle",
            )
        )

        elapsed = cal_elapsed_time(starttime, datetime.datetime.now())
        log.info(f"finished loading in {elapsed[0]:.2f} min ({elapsed[1]} sec)")

        # load gene model info
        # log.info("loading gene model info")
        # ENST_info = read_pickle_files(os.path.join("genome_files","parsed_gff3", config['genome_ver'],"ENST_info.pickle"))

        # load the mapping of ENST to the info file part
        log.info("loading gene model info **index to file parts**")
        ENST_info_index = read_pickle_files(
            os.path.join(
                "genome_files",
                "parsed_gff3",
                config["genome_ver"],
                "ENST_info",
                "ENST_info_index.pickle",
            )
        )

        elapsed = cal_elapsed_time(starttime, datetime.datetime.now())
        log.info(f"finished loading in {elapsed[0]:.2f} min ({elapsed[1]} sec)")

        # load codon phase index
        # log.info("loading codon phase info")
        # ENST_PhaseInCodon = read_pickle_files(os.path.join("genome_files","parsed_gff3", config['genome_ver'],"ENST_codonPhase.pickle"))

        log.info("loading codon phase info **index to file parts**")
        ENST_PhaseInCodon_index = read_pickle_files(
            os.path.join(
                "genome_files",
                "parsed_gff3",
                config["genome_ver"],
                "ENST_codonPhases",
                "ENST_codonPhase_index.pickle",
            )
        )

        # report time used
        elapsed = cal_elapsed_time(starttime, datetime.datetime.now())
        log.info(f"finished loading in {elapsed[0]:.2f} min ({elapsed[1]} sec)")

        # report the number of ENSTs which has ATG at the end of the exon
        # ExonEnd_ATG_count,ExonEnd_ATG_list = count_ATG_at_exonEnd(ENST_info)

        # open log files
        mkdir(outdir)
        mkdir(os.path.join(outdir, "genbank_files"))
        recut_CFD_all = open(os.path.join(outdir, "recut_CFD_all.txt"), "w")
        recut_CFD_fail = open(os.path.join(outdir, "recut_CFD_fail.txt"), "w")
        csvout_N = open(os.path.join(outdir, "out_Nterm_recut_cfd.csv"), "w")
        csvout_C = open(os.path.join(outdir, "out_Cterm_recut_cfd.csv"), "w")
        csvout_header = "ID,cfd1,cfd2,cfd3,cfd4,cfdScan,cfdScanNoRecode,cfd_max\n"
        csvout_N.write(csvout_header)
        csvout_C.write(csvout_header)
        fiveUTR_log = open(os.path.join(outdir, "fiveUTR.txt"), "w")

        # open result file and write header
        csvout_res = open(f"{outdir}/result.csv", "w")
        csvout_res.write(
            f"Entry,ID,chr,transcript_type,name,terminus,gRNA_name,gRNA_seq,PAM,gRNA_start,gRNA_end,gRNA_cut_pos,edit_pos,distance_between_cut_and_edit(cut_pos-insert_pos),specificity_score,specificity_weight,distance_weight,position_weight,final_weight,cfd_before_recoding,cfd_after_recoding,cfd_after_windowScan_and_recoding,max_recut_cfd,name_of_DNA_donor,DNA donor,name_of_trimmed_DNA_Donor,trimmed_DNA_donor,effective_HA_len,synthesis_problems,cutPos2nearestOffLimitJunc,strand(gene/gRNA/donor)\n"
        )   #"Entry,ID,chr,transcript_type,name,terminus,gRNA_seq,PAM,gRNA_start,gRNA_end,gRNA_cut_pos,edit_pos,distance_between_cut_and_edit(cut pos - insert pos),specificity_score,specificity_weight,distance_weight,position_weight,final_weight,cfd_before_recoding,cfd_after_recoding,cfd_after_windowScan_and_recoding,max_recut_cfd,DNA donor,effective_HA_len,synthesis_problems,cutPos2nearestOffLimitJunc,strand(gene/gRNA/donor)\n"

        # open result file2 for GenoPrimer input
        csvout_res2 = open(f"{outdir}/input_for_GenoPrimer.csv", "w")
        csvout_res2.write(f"Entry,ref,chr,coordinate,ID,geneSymbol\n")

        # dataframes to store best gRNAs
        best_start_gRNAs = pd.DataFrame()
        best_stop_gRNAs = pd.DataFrame()

        # logging cfd score, failed gRNAs etc
        start_info = info()
        stop_info = info()

        # load ENST list (the user input list or the whole transcriptome)
        if os.path.isfile(config["path2csv"]):
            log.info(
                f"begin processing user-supplied list of gene IDs in file {config['path2csv']}"
            )
            df = pd.read_csv(os.path.join(config["path2csv"]), dtype = str)
            # check csv columns
            keys2check = set(["Ensembl_ID"])
            if not keys2check.issubset(df.columns):
                log.error(
                    f'Missing columns in the input csv file\n Required columns:"Ensembl_ID"'
                )
                log.info(f"Please fix the input csv file and try again")
                sys.exit()
        else:
            sys.exit(f"ERROR: The input file {config['path2csv']} is not found")
            # log.warning(f"The input file {config['path2csv']} is not found, using the whole human transcriptome")
            # input("Press Enter to continue...")
            # df = pd.DataFrame(ENST_info.keys(), columns = ["Ensembl_ID"]) # create data frame from ENST_info

        # loop through each entry in the input csv file
        transcript_count = 0
        protein_coding_transcripts_count = 0
        target_terminus = "None"
        target_coordinate = "None"
        ENST_in_db = False
        ENST_design_counts = {} #used in the names of gRNA and donors
        Entry = 0 # Entry is defined by the portal, and if not using the portal, it is just the index of the row
        for index, row in df.iterrows():
            ENST_ID = row["Ensembl_ID"].rstrip().lstrip()
            if "Target_terminus" in df.columns and row.isnull()["Target_terminus"] == False:
                target_terminus = row["Target_terminus"].rstrip().lstrip().upper()

                if (
                    target_terminus != "N"
                    and target_terminus != "C"
                    and target_terminus != "ALL"
                    and target_terminus != ""
                ):
                    sys.exit(f"invalid target terminus: {target_terminus}")

            # determine if the input is ENST-based or coordinate-based
            ENST_based, coordinate_based = False, False
            if "Chromosome" in df.columns and "Coordinate" in df.columns and row.isnull()["Chromosome"] == False and row.isnull()["Coordinate"] == False:
                chrom = str(row["Chromosome"]).rstrip().lstrip()
                coordinate = str(row["Coordinate"]).rstrip().lstrip()
                if (
                    chrom != ""
                    and coordinate != ""
                    and coordinate.isdigit()
                ):
                    coordinate_based = True
            # if coordinate_based didn't check out, revert to ENST-based
            if coordinate_based == False:
                ENST_based = True
                if not target_terminus in ["N","C","ALL"]:
                    target_terminus = "ALL"

            if "Entry" in df.columns:
                try:
                    Entry = str(row["Entry"]).rstrip().lstrip()
                    ENST_in_db = True
                except:
                    ENST_in_db = False

            # check if ENST_ID is in the database
            if not ENST_ID in ENST_info_index.keys():
                if not ENST_in_db:
                    Entry += 1
                log.warning(
                    f"skipping {ENST_ID} b/c transcript is not in the annotated ENST collection (excluding those on chr_patch_hapl_scaff)"
                )
                genome_ver = config["genome_ver"]
                csvout_res.write(
                    f"{Entry},{ENST_ID},ERROR: this ID was not found in the genome {genome_ver}, most likely this ID was deprecated\n"
                )
                continue

            # check if codon phase info exists
            if not ENST_ID in ENST_PhaseInCodon_index.keys():
                if not ENST_in_db:
                    Entry += 1
                log.warning(
                    f"skipping {ENST_ID} b/c transcript has no codon phase information"
                )
                csvout_res.write(
                    f"{Entry},{ENST_ID},ERROR: this ID is either not protein-coding and/or has no codon phase information\n"
                )
                continue

            # load the ENST_info for current ID
            part = ENST_info_index[ENST_ID]
            ENST_info = read_pickle_files(
                os.path.join(
                    "genome_files",
                    "parsed_gff3",
                    config["genome_ver"],
                    "ENST_info",
                    f"ENST_info_part{part}.pickle",
                )
            )

            transcript_type = ENST_info[ENST_ID].description.split("|")[1]
            # if transcript_type == "protein_coding": # and ENST_ID == "ENST00000398165":
            # if not ENST_ID in ExonEnd_ATG_list: # only process edge cases in which genes with ATG are at the end of exons
            #     continue
            log.info(f"processing {ENST_ID}\ttranscript type: {transcript_type}")

            if hasattr(ENST_info[ENST_ID], "name"):
                name = ENST_info[ENST_ID].name
            else:
                name = ""
            row_prefix = f"{ENST_ID},{ENST_info[ENST_ID].chr},{transcript_type},{name}"

            # get codon_phase information for current ENST
            file_parts_list = ENST_PhaseInCodon_index[ENST_ID]
            ENST_PhaseInCodon = {}
            for part in file_parts_list:
                Codon_phase_dict = read_pickle_files(
                    os.path.join(
                        "genome_files",
                        "parsed_gff3",
                        config["genome_ver"],
                        "ENST_codonPhases",
                        f"ENST_codonPhase_part{str(part)}.pickle",
                    )
                )
                ENST_PhaseInCodon = deepmerge(
                    ENST_PhaseInCodon, Codon_phase_dict
                )  # merge file parts if ENST codon info is split among fileparts

            ######################################
            # best gRNA for a specific coordinate#
            ######################################
            if coordinate_based == True:
                if not ENST_in_db:
                    Entry += 1
                csvout_N.write(ENST_ID)
                #check if the coordinate is in the ENST
                if ENST_info[ENST_ID].chr == chrom and min(ENST_info[ENST_ID].span_start, ENST_info[ENST_ID].span_end) <= int(coordinate) <=  max(ENST_info[ENST_ID].span_start, ENST_info[ENST_ID].span_end):
                    
                    # get gRNAs
                    ranked_df_gRNAs_target_pos = get_gRNAs_target_coordinate(
                        ENST_ID=ENST_ID,
                        chrom = chrom,
                        pos = int(coordinate),
                        ENST_info=ENST_info,
                        freq_dict=freq_dict,
                        loc2file_index=loc2file_index,
                        loc2posType=loc2posType,
                        dist=max_cut2ins_dist,
                        genome_ver=config["genome_ver"],
                        pam=config["pam"],
                        spec_score_flavor=spec_score_flavor,
                        reg_penalty=reg_penalty,
                        alphas = alphas,
                    )

                    if ranked_df_gRNAs_target_pos.empty == True:
                        csvout_res.write(f"{Entry},{ENST_ID},ERROR: no suitable gRNAs found\n")

                    for i in range(0, min([gRNA_num_out, ranked_df_gRNAs_target_pos.shape[0]])):
                        current_gRNA = ranked_df_gRNAs_target_pos.iloc[[i]]

                        # get HDR template
                        try:
                            HDR_template = get_HDR_template(
                                df=current_gRNA,
                                ENST_info=ENST_info,
                                type="start", # this borrowed option specifies that the edit is immediately after the coordinate
                                ENST_PhaseInCodon=ENST_PhaseInCodon,
                                loc2posType=loc2posType,
                                genome_ver=config["genome_ver"],
                                HDR_arm_len=HDR_arm_len,
                                tag=config["POSpayload"],
                                ssODN_max_size=ssODN_max_size,
                                Donor_type=config["Donor_type"],
                                Strand_choice=config["Strand_choice"],
                                recoding_args=recoding_args,
                                syn_check_args=syn_check_args,
                            )
                        except Exception as e:
                            print("Unexpected error:", str(sys.exc_info()))
                            traceback.print_exc()
                            print("additional information:", e)
                            PrintException()

                        # append the best gRNA to the final df
                        if i == 0:
                            best_start_gRNAs = pd.concat([best_start_gRNAs, current_gRNA])

                        # append cfd score to list for plotting
                        pre_recoding_cfd_score = HDR_template.pre_recoding_cfd_score
                        cfd1 = ""
                        if hasattr(HDR_template, "cfd_score_post_mut_ins"):
                            cfd1 = HDR_template.cfd_score_post_mut_ins
                        if not hasattr(HDR_template, "cfd_score_post_mut2"):
                            cfd2 = cfd1
                        else:
                            cfd2 = HDR_template.cfd_score_post_mut2
                        if not hasattr(HDR_template, "cfd_score_post_mut3"):
                            cfd3 = cfd2
                        else:
                            cfd3 = HDR_template.cfd_score_post_mut3
                        if not hasattr(HDR_template, "cfd_score_post_mut4"):
                            cfd4 = cfd3
                        else:
                            cfd4 = HDR_template.cfd_score_post_mut4
                        cfd_scan = 0
                        cfd_scan_no_recode = 0
                        if hasattr(HDR_template, "cfd_score_highest_in_win_scan"):
                            cfd_scan = HDR_template.cfd_score_highest_in_win_scan
                            cfd_scan_no_recode = HDR_template.scan_highest_cfd_no_recode

                        cfdfinal = HDR_template.final_cfd

                        strands = f"{HDR_template.ENST_strand}/{HDR_template.gStrand}/{HDR_template.Donor_strand}"

                        # write csv
                        (
                            spec_score,
                            seq,
                            pam,
                            s,
                            e,
                            cut2ins_dist,
                            spec_weight,
                            dist_weight,
                            pos_weight,
                            final_weight,
                        ) = get_res(current_gRNA, spec_score_flavor)
                        donor = HDR_template.Donor_final
                        donor_trimmed = "N/A for ssODN"
                        if config["Donor_type"] == "dsDNA":
                            donor_trimmed = HDR_template.Donor_final
                            donor = HDR_template.Donor_pretrim
                            
                        # gRNA and donor names
                        ENST_design_counts[ENST_ID] = ENST_design_counts.get(ENST_ID, 0) + 1
                        gRNA_name = f"{ENST_ID}_gRNA_{ENST_design_counts[ENST_ID]}"
                        donor_name = f"{ENST_ID}_donor_{ENST_design_counts[ENST_ID]}"
                        if config["Donor_type"] == "dsDNA":
                            donor_trimmed_name = f"{ENST_ID}_donor_trimmed_{ENST_design_counts[ENST_ID]}"
                        else:
                            donor_trimmed_name = "N/A for ssODN"

                        gRNA_cut_pos = (
                            HDR_template.CutPos
                        )  # InsPos is the first letter of stop codon "T"AA or the last letter of the start codon AT"G"
                        insert_pos = HDR_template.InsPos
                        if config["recoding_off"]:
                            csvout_res.write(
                                f"{Entry},{row_prefix},-,{gRNA_name},{seq},{pam},{s},{e},{gRNA_cut_pos},{insert_pos},{cut2ins_dist},{spec_score},{ret_six_dec(spec_weight)},{ret_six_dec(dist_weight)},{ret_six_dec(pos_weight)},{ret_six_dec(final_weight)},{ret_six_dec(pre_recoding_cfd_score)},recoding turned off,,{ret_six_dec(cfdfinal)},{donor_name},{donor},{donor_trimmed_name},{donor_trimmed},{HDR_template.effective_HA_len},{HDR_template.synFlags},{HDR_template.cutPos2nearestOffLimitJunc},{strands}\n"
                            )
                            csvout_res2.write( f"{Entry},"+
                                config["genome_ver"]
                                + f",{HDR_template.ENST_chr},{insert_pos},{ENST_ID},{name}\n"
                            )
                        else:
                            if not isinstance(cfd4, float):
                                cfd4 = ""
                            csvout_res.write(
                                f"{Entry},{row_prefix},-,{gRNA_name},{seq},{pam},{s},{e},{gRNA_cut_pos},{insert_pos},{cut2ins_dist},{spec_score},{ret_six_dec(spec_weight)},{ret_six_dec(dist_weight)},{ret_six_dec(pos_weight)},{ret_six_dec(final_weight)},{ret_six_dec(pre_recoding_cfd_score)},{ret_six_dec(cfd4)},{ret_six_dec(cfd_scan)},{ret_six_dec(cfdfinal)},{donor_name},{donor},{donor_trimmed_name},{donor_trimmed},{HDR_template.effective_HA_len},{HDR_template.synFlags},{HDR_template.cutPos2nearestOffLimitJunc},{strands}\n"
                            )
                            csvout_res2.write( f"{Entry},"+
                                config["genome_ver"]
                                + f",{HDR_template.ENST_chr},{insert_pos},{ENST_ID},{name}\n"
                            )

                        # donor features
                        donor_features = HDR_template.Donor_features
                        # write genbank file
                        with open(os.path.join(outdir, "genbank_files", f"{donor_name}.gb"), "w") as gb_handle:
                            write_genbank(handle = gb_handle, data_obj = HDR_template, donor_name = donor_name, donor_type = config["Donor_type"])


                        # write log
                        this_log = f"{HDR_template.info}{HDR_template.info_arm}{HDR_template.info_p1}{HDR_template.info_p2}{HDR_template.info_p3}{HDR_template.info_p4}{HDR_template.info_p5}{HDR_template.info_p6}\n--------------------final CFD:{ret_six_dec(HDR_template.final_cfd)}\n    donor before any recoding:{HDR_template.Donor_vanillia}\n     donor after all recoding:{HDR_template.Donor_postMut}\ndonor centered(if applicable):{HDR_template.Donor_final}\n          donor (best strand):{HDR_template.Donor_final}\n\n"
                        this_log = f"{this_log}Donor features:\n{donor_features}\n\n"
                        recut_CFD_all.write(this_log)
                        if HDR_template.final_cfd > 0.03:
                            recut_CFD_fail.write(this_log)

                        if hasattr(HDR_template, "info_phase4_5UTR"):
                            fiveUTR_log.write(
                                f"phase4_UTR\t{HDR_template.info_phase4_5UTR[0]}\t{HDR_template.info_phase4_5UTR[1]}\n"
                            )
                        if hasattr(HDR_template, "info_phase5_5UTR"):
                            fiveUTR_log.write(
                                f"phase5_UTR\t{HDR_template.info_phase5_5UTR[0]}\t{HDR_template.info_phase5_5UTR[1]}\n"
                            )
                else:
                    csvout_res.write(f"{Entry},{ENST_ID},ERROR: provided genomic coordinates are not in the {ENST_ID}\n")

            ###################################
            # best start gRNA and HDR template#
            ###################################
            if ENST_based == True and (target_terminus == "ALL" or target_terminus == "N"):
                if not ENST_in_db:
                    Entry += 1
                csvout_N.write(ENST_ID)
                # get gRNAs
                ranked_df_gRNAs_ATG, ranked_df_gRNAs_stop = get_gRNAs(
                    ENST_ID=ENST_ID,
                    ENST_info=ENST_info,
                    freq_dict=freq_dict,
                    loc2file_index=loc2file_index,
                    loc2posType=loc2posType,
                    dist=max_cut2ins_dist,
                    genome_ver=config["genome_ver"],
                    pam=config["pam"],
                    spec_score_flavor=spec_score_flavor,
                    reg_penalty=reg_penalty,
                    alphas = alphas,
                )

                if ranked_df_gRNAs_ATG.empty == True:
                    start_info.failed.append(ENST_ID)
                    csvout_N.write(",,,,,\n")
                    csvout_res.write(f"{Entry},{ENST_ID},ERROR: no suitable gRNAs found\n")

                for i in range(0, min([gRNA_num_out, ranked_df_gRNAs_ATG.shape[0]])):
                    # if best_start_gRNA.shape[0] > 1: # multiple best scoring gRNA
                    #     best_start_gRNA = best_start_gRNA[best_start_gRNA["CSS"] == best_start_gRNA["CSS"].max()] # break the tie by CSS score
                    #     best_start_gRNA = best_start_gRNA.head(1) #get the first row in case of ties
                    current_gRNA = ranked_df_gRNAs_ATG.iloc[[i]]

                    # get HDR template
                    try:
                        HDR_template = get_HDR_template(
                            df=current_gRNA,
                            ENST_info=ENST_info,
                            type="start",
                            ENST_PhaseInCodon=ENST_PhaseInCodon,
                            loc2posType=loc2posType,
                            genome_ver=config["genome_ver"],
                            HDR_arm_len=HDR_arm_len,
                            tag=config["Npayload"],
                            ssODN_max_size=ssODN_max_size,
                            Donor_type=config["Donor_type"],
                            Strand_choice=config["Strand_choice"],
                            recoding_args=recoding_args,
                            syn_check_args=syn_check_args,
                        )
                    except Exception as e:
                        print("Unexpected error:", str(sys.exc_info()))
                        traceback.print_exc()
                        print("additional information:", e)
                        PrintException()

                    # append the best gRNA to the final df
                    if i == 0:
                        best_start_gRNAs = pd.concat([best_start_gRNAs, current_gRNA])

                    # append cfd score to list for plotting
                    pre_recoding_cfd_score = HDR_template.pre_recoding_cfd_score
                    cfd1 = ""
                    if hasattr(HDR_template, "cfd_score_post_mut_ins"):
                        cfd1 = HDR_template.cfd_score_post_mut_ins
                    if not hasattr(HDR_template, "cfd_score_post_mut2"):
                        cfd2 = cfd1
                    else:
                        cfd2 = HDR_template.cfd_score_post_mut2
                    if not hasattr(HDR_template, "cfd_score_post_mut3"):
                        cfd3 = cfd2
                    else:
                        cfd3 = HDR_template.cfd_score_post_mut3
                    if not hasattr(HDR_template, "cfd_score_post_mut4"):
                        cfd4 = cfd3
                    else:
                        cfd4 = HDR_template.cfd_score_post_mut4
                    cfd_scan = 0
                    cfd_scan_no_recode = 0
                    if hasattr(HDR_template, "cfd_score_highest_in_win_scan"):
                        cfd_scan = HDR_template.cfd_score_highest_in_win_scan
                        cfd_scan_no_recode = HDR_template.scan_highest_cfd_no_recode

                    cfdfinal = HDR_template.final_cfd

                    strands = f"{HDR_template.ENST_strand}/{HDR_template.gStrand}/{HDR_template.Donor_strand}"

                    # write csv
                    (
                        spec_score,
                        seq,
                        pam,
                        s,
                        e,
                        cut2ins_dist,
                        spec_weight,
                        dist_weight,
                        pos_weight,
                        final_weight,
                    ) = get_res(current_gRNA, spec_score_flavor)

                    donor = HDR_template.Donor_final
                    donor_trimmed = "N/A for ssODN"
                    if config["Donor_type"] == "dsDNA":
                        donor_trimmed = HDR_template.Donor_final
                        donor = HDR_template.Donor_pretrim
                        
                    # gRNA and donor names
                    ENST_design_counts[ENST_ID] = ENST_design_counts.get(ENST_ID, 0) + 1
                    gRNA_name = f"{ENST_ID}_gRNA_{ENST_design_counts[ENST_ID]}"
                    donor_name = f"{ENST_ID}_donor_{ENST_design_counts[ENST_ID]}"
                    if config["Donor_type"] == "dsDNA":
                        donor_trimmed_name = f"{ENST_ID}_donor_trimmed_{ENST_design_counts[ENST_ID]}"
                    else:
                        donor_trimmed_name = "N/A for ssODN"

                    gRNA_cut_pos = (
                        HDR_template.CutPos
                    )  # InsPos is the first letter of stop codon "T"AA or the last letter of the start codon AT"G"
                    insert_pos = HDR_template.InsPos
                    if config["recoding_off"]:
                        csvout_N.write(
                            f",{cfd1},{cfd2},{cfd3},{cfd4},{cfd_scan},{cfd_scan_no_recode},{cfdfinal}\n"
                        )
                        csvout_res.write(
                            f"{Entry},{row_prefix},N,{gRNA_name},{seq},{pam},{s},{e},{gRNA_cut_pos},{insert_pos},{cut2ins_dist},{spec_score},{ret_six_dec(spec_weight)},{ret_six_dec(dist_weight)},{ret_six_dec(pos_weight)},{ret_six_dec(final_weight)},{ret_six_dec(pre_recoding_cfd_score)},recoding turned off,,{ret_six_dec(cfdfinal)},{donor_name},{donor},{donor_trimmed_name},{donor_trimmed},{HDR_template.effective_HA_len},{HDR_template.synFlags},{HDR_template.cutPos2nearestOffLimitJunc},{strands}\n"
                        )
                        csvout_res2.write(f"{Entry},"+
                            config["genome_ver"]
                            + f",{HDR_template.ENST_chr},{insert_pos},{ENST_ID},{name}\n"
                        )
                    else:
                        csvout_N.write(
                            f",{cfd1},{cfd2},{cfd3},{cfd4},{cfd_scan},{cfd_scan_no_recode},{cfdfinal}\n"
                        )
                        if not isinstance(cfd4, float):
                            cfd4 = ""
                        csvout_res.write(
                            f"{Entry},{row_prefix},N,{gRNA_name},{seq},{pam},{s},{e},{gRNA_cut_pos},{insert_pos},{cut2ins_dist},{spec_score},{ret_six_dec(spec_weight)},{ret_six_dec(dist_weight)},{ret_six_dec(pos_weight)},{ret_six_dec(final_weight)},{ret_six_dec(pre_recoding_cfd_score)},{ret_six_dec(cfd4)},{ret_six_dec(cfd_scan)},{ret_six_dec(cfdfinal)},{donor_name},{donor},{donor_trimmed_name},{donor_trimmed},{HDR_template.effective_HA_len},{HDR_template.synFlags},{HDR_template.cutPos2nearestOffLimitJunc},{strands}\n"
                        )
                        csvout_res2.write(f"{Entry},"+
                            config["genome_ver"]
                            + f",{HDR_template.ENST_chr},{insert_pos},{ENST_ID},{name}\n"
                        )

                    # donor features
                    donor_features = HDR_template.Donor_features
                    # write genbank file
                    with open(os.path.join(outdir, "genbank_files", f"{donor_name}.gb"), "w") as gb_handle:
                        write_genbank(handle = gb_handle, data_obj = HDR_template, donor_name = donor_name, donor_type = config["Donor_type"])

                    # write log
                    this_log = f"{HDR_template.info}{HDR_template.info_arm}{HDR_template.info_p1}{HDR_template.info_p2}{HDR_template.info_p3}{HDR_template.info_p4}{HDR_template.info_p5}{HDR_template.info_p6}\n--------------------final CFD:{ret_six_dec(HDR_template.final_cfd)}\n    donor before any recoding:{HDR_template.Donor_vanillia}\n     donor after all recoding:{HDR_template.Donor_postMut}\ndonor centered(if applicable):{HDR_template.Donor_final}\n          donor (best strand):{HDR_template.Donor_final}\n\n"
                    this_log = f"{this_log}Donor features:\n{donor_features}\n\n"
                    recut_CFD_all.write(this_log)
                    if HDR_template.final_cfd > 0.03:
                        recut_CFD_fail.write(this_log)

                    if hasattr(HDR_template, "info_phase4_5UTR"):
                        fiveUTR_log.write(
                            f"phase4_UTR\t{HDR_template.info_phase4_5UTR[0]}\t{HDR_template.info_phase4_5UTR[1]}\n"
                        )
                    if hasattr(HDR_template, "info_phase5_5UTR"):
                        fiveUTR_log.write(
                            f"phase5_UTR\t{HDR_template.info_phase5_5UTR[0]}\t{HDR_template.info_phase5_5UTR[1]}\n"
                        )

            ##################################
            # best stop gRNA and HDR template#
            ##################################
            if ENST_based == True and (target_terminus == "ALL" or target_terminus == "C"):
                if not ENST_in_db:
                    Entry += 1
                csvout_C.write(ENST_ID)

                # get gRNAs
                ranked_df_gRNAs_ATG, ranked_df_gRNAs_stop = get_gRNAs(
                    ENST_ID=ENST_ID,
                    ENST_info=ENST_info,
                    freq_dict=freq_dict,
                    loc2file_index=loc2file_index,
                    loc2posType=loc2posType,
                    dist=max_cut2ins_dist,
                    genome_ver=config["genome_ver"],
                    pam=config["pam"],
                    spec_score_flavor=spec_score_flavor,
                    reg_penalty=reg_penalty,
                    alphas = alphas,
                )

                if ranked_df_gRNAs_stop.empty == True:
                    stop_info.failed.append(ENST_ID)
                    csvout_C.write(",,,,,\n")
                    csvout_res.write(f"{Entry},{ENST_ID},ERROR: no suitable gRNAs found\n")

                for i in range(0, min([gRNA_num_out, ranked_df_gRNAs_stop.shape[0]])):
                    # if best_stop_gRNA.shape[0] > 1: # multiple best scoring gRNA
                    #     best_stop_gRNA = best_stop_gRNA[best_stop_gRNA["CSS"] == best_stop_gRNA["CSS"].max()] # break the tie by CSS score
                    #     best_stop_gRNA = best_stop_gRNA.head(1) #get the first row in case of ties
                    current_gRNA = ranked_df_gRNAs_stop.iloc[[i]]

                    # get HDR template
                    try:
                        HDR_template = get_HDR_template(
                            df=current_gRNA,
                            ENST_info=ENST_info,
                            type="stop",
                            ENST_PhaseInCodon=ENST_PhaseInCodon,
                            loc2posType=loc2posType,
                            HDR_arm_len=HDR_arm_len,
                            genome_ver=config["genome_ver"],
                            tag=config["Cpayload"],
                            Donor_type=config["Donor_type"],
                            Strand_choice=config["Strand_choice"],
                            ssODN_max_size=ssODN_max_size,
                            recoding_args=recoding_args,
                            syn_check_args=syn_check_args,
                        )
                    except Exception as e:
                        print("Unexpected error:", str(sys.exc_info()))
                        traceback.print_exc()
                        print("additional information:", e)
                        PrintException()

                    # append the best gRNA to the final df
                    best_stop_gRNAs = pd.concat([best_stop_gRNAs, current_gRNA])

                    # append cfd score to list for plotting
                    pre_recoding_cfd_score = HDR_template.pre_recoding_cfd_score
                    cfd1 = ""
                    if hasattr(HDR_template, "cfd_score_post_mut_ins"):
                        cfd1 = HDR_template.cfd_score_post_mut_ins
                    if not hasattr(HDR_template, "cfd_score_post_mut2"):
                        cfd2 = cfd1
                    else:
                        cfd2 = HDR_template.cfd_score_post_mut2
                    if not hasattr(HDR_template, "cfd_score_post_mut3"):
                        cfd3 = cfd2
                    else:
                        cfd3 = HDR_template.cfd_score_post_mut3
                    if not hasattr(HDR_template, "cfd_score_post_mut4"):
                        cfd4 = cfd3
                    else:
                        cfd4 = HDR_template.cfd_score_post_mut4

                    cfd_scan = 0
                    cfd_scan_no_recode = 0
                    if hasattr(HDR_template, "cfd_score_highest_in_win_scan"):
                        cfd_scan = HDR_template.cfd_score_highest_in_win_scan
                        cfd_scan_no_recode = HDR_template.scan_highest_cfd_no_recode

                    cfdfinal = HDR_template.final_cfd

                    strands = f"{HDR_template.ENST_strand}/{HDR_template.gStrand}/{HDR_template.Donor_strand}"

                    # write csv
                    (
                        spec_score,
                        seq,
                        pam,
                        s,
                        e,
                        cut2ins_dist,
                        spec_weight,
                        dist_weight,
                        pos_weight,
                        final_weight,
                    ) = get_res(current_gRNA, spec_score_flavor)

                    donor = HDR_template.Donor_final
                    donor_trimmed = "N/A for ssODN"
                    if config["Donor_type"] == "dsDNA":
                        donor_trimmed = HDR_template.Donor_final
                        donor = HDR_template.Donor_pretrim
                        
                    # gRNA and donor names
                    ENST_design_counts[ENST_ID] = ENST_design_counts.get(ENST_ID, 0) + 1
                    gRNA_name = f"{ENST_ID}_gRNA_{ENST_design_counts[ENST_ID]}"
                    donor_name = f"{ENST_ID}_donor_{ENST_design_counts[ENST_ID]}"
                    if config["Donor_type"] == "dsDNA":
                        donor_trimmed_name = f"{ENST_ID}_donor_trimmed_{ENST_design_counts[ENST_ID]}"
                    else:
                        donor_trimmed_name = "N/A for ssODN"

                    gRNA_cut_pos = (
                        HDR_template.CutPos
                    )  # InsPos is the first letter of stop codon "T"AA or the last letter of the start codon AT"G"
                    insert_pos = HDR_template.InsPos
                    if config["recoding_off"]:
                        csvout_C.write(
                            f",{cfd1},{cfd2},{cfd3},{cfd4},{cfd_scan},{cfd_scan_no_recode},{cfdfinal}\n"
                        )
                        csvout_res.write(
                            f"{Entry},{row_prefix},C,{gRNA_name},{seq},{pam},{s},{e},{gRNA_cut_pos},{insert_pos},{cut2ins_dist},{spec_score},{ret_six_dec(spec_weight)},{ret_six_dec(dist_weight)},{ret_six_dec(pos_weight)},{ret_six_dec(final_weight)},{ret_six_dec(pre_recoding_cfd_score)},recoding turned off,,{ret_six_dec(cfdfinal)},{donor_name},{donor},{donor_trimmed_name},{donor_trimmed},{HDR_template.effective_HA_len},{HDR_template.synFlags},{HDR_template.cutPos2nearestOffLimitJunc},{strands}\n"
                        )
                        csvout_res2.write( f"{Entry},"+
                            config["genome_ver"]
                            + f",{HDR_template.ENST_chr},{insert_pos},{ENST_ID},{name}\n"
                        )
                    else:
                        csvout_C.write(
                            f",{cfd1},{cfd2},{cfd3},{cfd4},{cfd_scan},{cfd_scan_no_recode},{cfdfinal}\n"
                        )
                        if not isinstance(cfd4, float):
                            cfd4 = ""
                        csvout_res.write(
                            f"{Entry},{row_prefix},C,{gRNA_name},{seq},{pam},{s},{e},{gRNA_cut_pos},{insert_pos},{cut2ins_dist},{spec_score},{ret_six_dec(spec_weight)},{ret_six_dec(dist_weight)},{ret_six_dec(pos_weight)},{ret_six_dec(final_weight)},{ret_six_dec(pre_recoding_cfd_score)},{ret_six_dec(cfd4)},{ret_six_dec(cfd_scan)},{ret_six_dec(cfdfinal)},{donor_name},{donor},{donor_trimmed_name},{donor_trimmed},{HDR_template.effective_HA_len},{HDR_template.synFlags},{HDR_template.cutPos2nearestOffLimitJunc},{strands}\n"
                        )
                        csvout_res2.write(f"{Entry},"+
                            config["genome_ver"]
                            + f",{HDR_template.ENST_chr},{insert_pos},{ENST_ID},{name}\n"
                        )

                    # donor features
                    donor_features = HDR_template.Donor_features
                    # write genbank file
                    with open(os.path.join(outdir, "genbank_files", f"{donor_name}.gb"), "w") as gb_handle:
                        write_genbank(handle = gb_handle, data_obj = HDR_template, donor_name = donor_name, donor_type = config["Donor_type"])

                    # write log
                    this_log = f"{HDR_template.info}{HDR_template.info_arm}{HDR_template.info_p1}{HDR_template.info_p2}{HDR_template.info_p3}{HDR_template.info_p4}{HDR_template.info_p5}{HDR_template.info_p6}\n--------------------final CFD:{ret_six_dec(HDR_template.final_cfd)}\n   donor before any recoding:{HDR_template.Donor_vanillia}\n    donor after all recoding:{HDR_template.Donor_postMut}\n             donor centered:{HDR_template.Donor_final}\ndonor centered (best strand):{HDR_template.Donor_final}\n\n"
                    this_log = f"{this_log}Donor features:\n{donor_features}\n\n"
                    recut_CFD_all.write(this_log)
                    if HDR_template.final_cfd > 0.03:
                        recut_CFD_fail.write(this_log)

                    if hasattr(HDR_template, "info_phase4_5UTR"):
                        fiveUTR_log.write(
                            f"phase4_UTR\t{HDR_template.info_phase4_5UTR[0]}\t{HDR_template.info_phase4_5UTR[1]}\n"
                        )
                    if hasattr(HDR_template, "info_phase5_5UTR"):
                        fiveUTR_log.write(
                            f"phase5_UTR\t{HDR_template.info_phase5_5UTR[0]}\t{HDR_template.info_phase5_5UTR[1]}\n"
                        )

            protein_coding_transcripts_count += 1
            # else:
            #     log.info(f"skipping {ENST_ID} transcript type: {transcript_type} b/c transcript is not protein_coding")
            transcript_count += 1
            # report progress
            if (
                protein_coding_transcripts_count % 100 == 0
                and protein_coding_transcripts_count != 0
            ):
                endtime = datetime.datetime.now()
                elapsed_sec = endtime - starttime
                elapsed_min = elapsed_sec.seconds / 60
                log.info(
                    f"processed {protein_coding_transcripts_count}/{transcript_count} transcripts, elapsed time {elapsed_min:.2f} min ({elapsed_sec} sec)"
                )

        # write csv out
        endtime = datetime.datetime.now()
        elapsed_sec = endtime - starttime
        elapsed_min = elapsed_sec.seconds / 60

        if "num_to_process" in locals():
            pass
        else:
            num_to_process = "all"

        log.info(
            f"finished in {elapsed_min:.2f} min ({elapsed_sec} sec) , processed {protein_coding_transcripts_count}/{transcript_count} transcripts\nnonprotein-coding transcripts were skipped\n"
            f"results written to {config['outdir']}"
        )

        recut_CFD_all.close()
        recut_CFD_fail.close()
        csvout_N.close()
        csvout_C.close()
        fiveUTR_log.close()
        csvout_res.close()
        csvout_res2.close()

    except Exception as e:
        print("Unexpected error:", str(sys.exc_info()))
        traceback.print_exc()
        print("additional information:", e)
        PrintException()

## end of main()

def write_genbank(handle, data_obj, donor_name, donor_type):
    gb_record = Record.Record()
    gb_record.locus = donor_name
    gb_record.size = len(data_obj.Donor_final)
    gb_record.residue_type = 'DNA'
    gb_record.data_file_division = 'PLN'
    gb_record.definition = ''
    gb_record.accession = ['']
    gb_record.version = ''
    gb_record.keywords = ['DNA donor']
    gb_record.source = "synthetic DNA donor"
    gb_record.organism = 'synthetic DNA donor'
    gb_record.taxonomy = ['synthetic DNA donor']
    gb_record.sequence=data_obj.Donor_final

    # Convert the Bio.GenBank record to a SeqRecord (needed for writing with SeqIO)
    sequence = Seq(gb_record.sequence)
    seq_record = SeqRecord(sequence, id=gb_record.accession[0], name=donor_name, description= f"DNA donor type: {donor_type}")
    seq_record.annotations["date"] = get_current_date_formatted()
    seq_record.annotations["data_file_division"] = gb_record.data_file_division
    seq_record.annotations["organism"] = 'synthetic DNA donor'
    seq_record.annotations["molecule_type"] = "DNA"

    # Features 
    if "gRNA_coord" in data_obj.Donor_features:
        for feat in data_obj.Donor_features["gRNA_coord"]:
            feature = SeqFeature(FeatureLocation(start=feat[0], end=feat[1], strand=data_obj.Donor_features["gRNA_strand"]), type='gRNA+PAM', qualifiers={"label": "gRNA+PAM"})
            seq_record.features.append(feature)
    if "recoding_coord" in data_obj.Donor_features:
        for feat in data_obj.Donor_features["recoding_coord"]:
            feature = SeqFeature(FeatureLocation(start=feat[0], end=feat[1]), type='recode', qualifiers={"label": "recode"})
            seq_record.features.append(feature)

    feature = SeqFeature(FeatureLocation(start=data_obj.Donor_features["left_arm_coord"][0], end=data_obj.Donor_features["left_arm_coord"][1], strand=data_obj.Donor_features["HA_payload_strand"]), type='left HA', qualifiers={"label": "left HA"})
    seq_record.features.append(feature)

    feature = SeqFeature(FeatureLocation(start=data_obj.Donor_features["right_arm_coord"][0], end=data_obj.Donor_features["right_arm_coord"][1], strand=data_obj.Donor_features["HA_payload_strand"]), type='right HA', qualifiers={"label": "right HA"})
    seq_record.features.append(feature)

    feature = SeqFeature(FeatureLocation(start=data_obj.Donor_features["tag_coord"][0], end=data_obj.Donor_features["tag_coord"][1], strand=data_obj.Donor_features["HA_payload_strand"]), type='payload', qualifiers={"label": "payload"})
    seq_record.features.append(feature)


    # Write to a GenBank file using SeqIO for the actual file writing
    SeqIO.write([seq_record], handle, 'genbank')

def ret_six_dec(myvar):
    """
    retain six decimal points for printout
    """
    if type(myvar) == float:
        return f"{myvar:.6f}"
    elif type(myvar) == int:
        myvar = float(myvar)
        return f"{myvar:.6f}"
    else:
        return myvar


def mkdir(mypath):
    if not os.path.exists(mypath):
        os.makedirs(mypath)


def deepmerge(dict1, dict2):
    """
    merge two dictionary at the secondary key level
    """
    if len(dict1) == 0:
        return dict2
    if len(dict2) == 0:
        return dict1
    # start merging
    dictm = dict1
    for k in dict2:
        if k in dict1:  # shared key
            for k2 in dict2[k].keys():
                dictm[k][k2] = dict2[k][k2]
        else:
            dictm[k] = dict2[k]
    return dictm


class info:
    """
    info log class
    """
    def __init__(self) -> None:
        self.cfd1 = []
        self.cfd2 = []
        self.cfd3 = []
        self.cfd4 = []
        self.cfdfinal = []
        self.failed = []


def get_res(best_start_gRNA, spec_score_flavor):
    spec_score = best_start_gRNA[spec_score_flavor].values[0]
    seq = best_start_gRNA["seq"].values[0]
    pam = best_start_gRNA["pam"].values[0]
    s = best_start_gRNA["start"].values[0]
    e = best_start_gRNA["end"].values[0]
    cut2ins_dist = best_start_gRNA["Cut2Ins_dist"].values[0]
    spec_weight = best_start_gRNA["spec_weight"].values[0]
    dist_weight = best_start_gRNA["dist_weight"].values[0]
    pos_weight = best_start_gRNA["pos_weight"].values[0]
    final_weight = best_start_gRNA["final_weight"].values[0]
    return [
        spec_score,
        seq,
        pam,
        s,
        e,
        cut2ins_dist,
        spec_weight,
        dist_weight,
        pos_weight,
        final_weight,
    ]


def test_memory(n):
    """
    try allocate n MB of memory
    return true if can, and false otherwise
    """
    try:
        x = bytearray(1024 * 1000 * n)
        del x
        return True
    except:
        return False

def get_current_date_formatted():
    # Get the current date
    current_date = datetime.datetime.now()
    # Format the date as yyyy-mm-dd
    formatted_date = current_date.strftime("%d-%b-%Y").upper()
    return formatted_date

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


if __name__ == "__main__":
    main()
