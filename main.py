import os.path
from scripts.utils import *
import traceback
import time

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

def parse_args():
    parser= MyParser(description='ProtospaceX')
    parser.add_argument('--genome_ver', default="GRCh38", type=str, help='currently supports three genomes: GRCh38, GRCm39, GRCz11, ', metavar='')
    parser.add_argument('--path2csv',   default="input/mart_export_canonical_proteincoding.csv", type=str,help='path to a csv file containing ENST information\n *required columns*: Ensemble_ID',metavar='')

    #gRNA
    parser.add_argument('--num_gRNA_per_term',  default=1, type=int, help='payload at the N terminus', metavar='')

    #donor
    parser.add_argument('--HA_len',  default=500, help='length of the homology arm (on each side), will be the final arm length for dsDNA donors', type=int, metavar='')
    parser.add_argument('--Donor_type',  default="ssDNA", help='ssDNA(default) or dsDNA', type=str, metavar='')
    parser.add_argument('--Strand_choice',  default="auto", help='only applies when --Donor_type is set to ssDNA, possible values are auto,TargetStrand,NonTargetStrand,CodingStrand,NonCodingStrand', type=str, metavar='')
    parser.add_argument('--ssDNA_max_size', type=int, help='only applies when --Donor_type is set to ssDNA. Enforce a length restraint of the donor (both arms + payload), setting this option will center the ssODN with respect to the payload and the recoded region', metavar='')
    parser.add_argument('--CheckEnzymes',  default="", help='Restriction enzyme sites to check, separated by |, for example: BsaI|EcoRI', type=str, metavar='')
    parser.add_argument('--CustomSeq2Avoid',  default="", help='custom sequences to avoid, separated by |', type=str, metavar='')
    parser.add_argument('--MinArmLenPostTrim',  default=0, help='Minimum length the homology arm after trimming,  default is 0 (turning off trimming)', type=str, metavar='')

    #payload
    parser.add_argument('--payload',   default="", type=str, help='payload, overrides --Npayloadf and --Cpayload, --Tag, --Linker', metavar='')
    parser.add_argument('--Npayload',  default="", type=str, help='payload at the N terminus, overrides --Tag --Linker', metavar='')
    parser.add_argument('--Cpayload',  default="", type=str, help='payload at the N terminus, overrides --Tag --Linker', metavar='')
    parser.add_argument('--Tag',  default="ACCGAGCTCAACTTCAAGGAGTGGCAAAAGGCCTTTACCGATATGATG", type=str, help='default is the mNG11 tag', metavar='')
    parser.add_argument('--Linker',  default="GGTGGCGGATTGGAAGTTTTGTTTCAAGGTCCAGGAAGTGGT", type=str, help='default is the GS linker', metavar='')

    #recoding
    parser.add_argument('--recoding_off',             default = False, action='store_true', help='turn off *all* recoding')
    parser.add_argument('--recoding_stop_recut_only', default = False, action='store_true', help='use recoding to prevent recut')
    parser.add_argument('--recoding_full',            default = False, action='store_true', help='use recoding to prevent recut + recode region between insert and cut site')
    parser.add_argument('--cfdThres',                 default = 0.03, help='ProtospaceX will attempt to lower the recut cfd to this threshold (by recoding), cfd values lower than the threshold will be considered not suceptible to being recut anymore.')
    parser.add_argument('--recode_order',             default = "protospacer_first", help='possible values: protospacer_first, PAM_first')

    #output
    parser.add_argument('--outdir',   default="logs", type=str, help='output directory')

    config = parser.parse_args()
    return config

logging.setLoggerClass(ColoredLogger)
#logging.basicConfig()
log = logging.getLogger("ProtospaceX")
log.propagate = False
log.setLevel(logging.INFO) #set the level of warning displayed
#log.setLevel(logging.DEBUG) #set the level of warning displayed

#configs
config = vars(parse_args())
gRNA_num_out = config['num_gRNA_per_term']
max_cut2ins_dist = 50 #deprecated?
HDR_arm_len = config['HA_len']
ssDNA_max_size = config["ssDNA_max_size"]
spec_score_flavor = "guideMITScore"
outdir = config['outdir']
syn_check_args = {
                    "check_enzymes": config["CheckEnzymes"],
                    "CustomSeq2Avoid": config["CustomSeq2Avoid"],
                    "MinArmLenPostTrim": config["MinArmLenPostTrim"]
} # dictionary for multiple synthesis check arguments

#check recoding args
assert config["recode_order"] == "protospacer_first" or config["recode_order"] == "PAM_first"
if (config["recoding_full"] and any([config["recoding_off"],config["recoding_stop_recut_only"]])):
    sys.exit(f"Found conflicts in recoding arguments: --recoding_full cannot be used with --recoding_off or --recoding_stop_recut_only\nplease correct the issue and try again")
if config["recoding_off"] and config["recoding_stop_recut_only"]:
    sys.exit(f"Found conflicts in recoding arguments: --recoding_off cannot be used with --recoding_stop_recut_only\nplease correct the issue and try again")

#process recoding args
if (not config["recoding_off"]) and (not config["recoding_stop_recut_only"]):
    config["recoding_full"] = True

if config["recoding_off"] or config["recoding_stop_recut_only"]:
    config["recoding_full"] = False

recoding_args = {"recoding_off":config["recoding_off"],
                 "recoding_stop_recut_only":config["recoding_stop_recut_only"],
                 "recoding_full":config["recoding_full"],
                 "cfdThres":float(config['cfdThres']),
                 "recode_order":config['recode_order']}

#check donor args
if not config["Donor_type"] in ["ssDNA", "dsDNA"]:
    sys.exit("Donor_type must be ssDNA or dsDNA, offending value:" + config["Donor_type"] + ", please correct the issue and try again")
if not config["Strand_choice"] in ["auto","TargetStrand","NonTargetStrand","CodingStrand","NonCodingStrand"]:
    sys.exit("Strand_choice must be auto,TargetStrand,NonTargetStrand,CodingStrand or NonCodingStrand, offending value:" + config["Strand_choice"] + ", please correct the issue and try again")

#parse payload
Linker = config["Linker"]
Tag = config["Tag"]

if config["payload"] == "": #no payload override
    if config["Npayload"] == "": #no Npayload override
        config["Npayload"] = Tag + Linker
    if config["Cpayload"] == "": #no Cpayload override
        config["Cpayload"] = Linker + Tag
else: #payload override
    config["Npayload"] = config["payload"]
    config["Cpayload"] = config["payload"]


#check if HA_len is too short to satisfy ssDNA_max_size
if ssDNA_max_size is not None:
    max_payload_size = max([len(config["Npayload"]),len(config["Cpayload"])])
    derived_HDR_arm_len = ssDNA_max_size- max_payload_size / 2
    if derived_HDR_arm_len >= HDR_arm_len:
        print(f"HA_len={HDR_arm_len} is to short to meet the requirement of ssDNA_max_size={ssDNA_max_size}, payload size={max_payload_size}\n ssDNA_max_size={ssDNA_max_size} requires HA_len = ssDNA_max_size- max_payload_size / 2 = {derived_HDR_arm_len}")
        HDR_arm_len = derived_HDR_arm_len + 100
        print(f"HA_len is adjusted to {HDR_arm_len}")

#####################
##      main       ##
#####################
def main(outdir):
    try:
        #check memory requirement
        enough_mem = test_memory(4200)
        while not enough_mem: # test if at least 4.2 GB memory is available
            time.sleep(5) #retry in 5 seconds
            enough_mem = test_memory(4200)

        starttime = datetime.datetime.now()
        freq_dict = dict()

        #load gRNA info index (mapping of chromosomal location to file parts)
        log.info("loading precomputed gRNA **index to file parts**")
        loc2file_index = read_pickle_files(os.path.join("precomuted_gRNAs", "gRNA_" + config['genome_ver'],"gRNA.tab.gz.split.BwaMapped.scored","loc2file_index.pickle"))

        elapsed = cal_elapsed_time(starttime,datetime.datetime.now())
        log.info(f"finished loading in {elapsed[0]:.2f} min ({elapsed[1]} sec)")

        #load chr location to type (e.g. UTR, cds, exon/intron junction) mappings
        log.info("loading the mapping of chromosomal location to type (e.g. UTR, cds, exon/intron junction)")
        loc2posType = read_pickle_files(os.path.join("genome_files","parsed_gff3", config['genome_ver'],"loc2posType.pickle"))

        elapsed = cal_elapsed_time(starttime,datetime.datetime.now())
        log.info(f"finished loading in {elapsed[0]:.2f} min ({elapsed[1]} sec)")

        #load gene model info
        #log.info("loading gene model info")
        #ENST_info = read_pickle_files(os.path.join("genome_files","parsed_gff3", config['genome_ver'],"ENST_info.pickle"))

        #load the mapping of ENST to the info file part
        log.info("loading gene model info **index to file parts**")
        ENST_info_index = read_pickle_files(os.path.join("genome_files", "parsed_gff3", config['genome_ver'],"ENST_info","ENST_info_index.pickle"))

        elapsed = cal_elapsed_time(starttime,datetime.datetime.now())
        log.info(f"finished loading in {elapsed[0]:.2f} min ({elapsed[1]} sec)")

        #load codon phase index
        #log.info("loading codon phase info")
        #ENST_PhaseInCodon = read_pickle_files(os.path.join("genome_files","parsed_gff3", config['genome_ver'],"ENST_codonPhase.pickle"))

        log.info("loading codon phase info **index to file parts**")
        ENST_PhaseInCodon_index = read_pickle_files(os.path.join("genome_files","parsed_gff3", config['genome_ver'],"ENST_codonPhases","ENST_codonPhase_index.pickle"))

        #report time used
        elapsed = cal_elapsed_time(starttime,datetime.datetime.now())
        log.info(f"finished loading in {elapsed[0]:.2f} min ({elapsed[1]} sec)")

        #report the number of ENSTs which has ATG at the end of the exon
        #ExonEnd_ATG_count,ExonEnd_ATG_list = count_ATG_at_exonEnd(ENST_info)

        #open log files
        mkdir(outdir)
        recut_CFD_all = open(os.path.join(outdir,"recut_CFD_all.txt"), "w")
        recut_CFD_fail = open(os.path.join(outdir,"recut_CFD_fail.txt"), "w")
        csvout_N = open(os.path.join(outdir,"out_Nterm_recut_cfd.csv"), "w")
        csvout_C = open(os.path.join(outdir,"out_Cterm_recut_cfd.csv"), "w")
        csvout_header = "ID,cfd1,cfd2,cfd3,cfd4,cfdScan,cfdScanNoRecode,cfd_max\n"
        csvout_N.write(csvout_header)
        csvout_C.write(csvout_header)
        fiveUTR_log = open(os.path.join(outdir,"fiveUTR.txt"), "w")

        #open result file and write header
        csvout_res = open(f"{outdir}/result.csv", "w")
        csvout_res.write(f"ID,chr,transcript_type,name,terminus,gRNA_seq,PAM,gRNA_start,gRNA_end,gRNA_cut_pos,edit_pos,distance_between_cut_and_edit(cut pos - insert pos),specificity_score,specificity_weight,distance_weight,position_weight,final_weight,cfd_before_recoding,cfd_after_recoding,cfd_after_windowScan_and_recoding,max_recut_cfd,DNA donor,effective_HA_len,synthesis_problems\n")
        #open result file2 for GenoPrimer input
        csvout_res2 = open(f"{outdir}/input_for_GenoPrimer.csv", "w")
        csvout_res2.write(f"ref,chr,coordinate\n")


        #dataframes to store best gRNAs
        best_start_gRNAs = pd.DataFrame()
        best_stop_gRNAs = pd.DataFrame()

        #logging cfd score, failed gRNAs etc
        start_info = info()
        stop_info = info()

        #load ENST list (the user input list or the whole transcriptome)
        if os.path.isfile(config['path2csv']):
            log.info(f"begin processing user-supplied list of gene IDs in file {config['path2csv']}")
            df = pd.read_csv(os.path.join(config['path2csv']))
            # check csv columns
            keys2check = set(["Ensemble_ID"])
            if not keys2check.issubset(df.columns):
                log.error(f"Missing columns in the input csv file\n Required columns:\"Ensemble_ID\"")
                log.info(f"Please fix the input csv file and try again")
                sys.exit()
        else:
            sys.exit(f"ERROR: The input file {config['path2csv']} is not found")
            # log.warning(f"The input file {config['path2csv']} is not found, using the whole human transcriptome")
            # input("Press Enter to continue...")
            # df = pd.DataFrame(ENST_info.keys(), columns = ["Ensemble_ID"]) # create data frame from ENST_info

        #loop through each ENST
        transcript_count = 0
        protein_coding_transcripts_count = 0
        target_terminus = "all"
        for index, row in df.iterrows():
            ENST_ID = row["Ensemble_ID"].rstrip().lstrip()
            if "Target_terminus" in df.columns:
                target_terminus = row["Target_terminus"].rstrip().lstrip().upper()

                if target_terminus!="N" and target_terminus!="C" and target_terminus!="all":
                    sys.exit(f"invalid target terminus: {target_terminus}")

            #check if ENST_ID is in the database
            if not ENST_ID in ENST_info_index.keys():
                log.warning(f"skipping {ENST_ID} b/c transcript is not in the annotated ENST collection (excluding those on chr_patch_hapl_scaff)")
                genome_ver = config['genome_ver']
                csvout_res.write(f"{ENST_ID},ERROR: this ID was not found in the genome {genome_ver}, most likely this ID was deprecated\n")
                continue

            #check if codon phase info exists
            if not ENST_ID in ENST_PhaseInCodon_index.keys():
                log.warning(f"skipping {ENST_ID} b/c transcript has no codon phase information")
                csvout_res.write(f"{ENST_ID},ERROR: this ID is either not protein-coding and/or has no codon phase information\n")
                continue

            #load the ENST_info for current ID
            part = ENST_info_index[ENST_ID]
            ENST_info = read_pickle_files(os.path.join("genome_files", "parsed_gff3", config['genome_ver'],"ENST_info",f"ENST_info_part{part}.pickle"))

            transcript_type = ENST_info[ENST_ID].description.split("|")[1]
            # if transcript_type == "protein_coding": # and ENST_ID == "ENST00000398165":
            # if not ENST_ID in ExonEnd_ATG_list: # only process edge cases in which genes with ATG are at the end of exons
            #     continue
            log.info(f"processing {ENST_ID}\ttranscript type: {transcript_type}")

            if hasattr(ENST_info[ENST_ID],"name"):
                name = ENST_info[ENST_ID].name
            else:
                name = ""
            row_prefix = f"{ENST_ID},{ENST_info[ENST_ID].chr},{transcript_type},{name}"

            #get gRNAs
            ranked_df_gRNAs_ATG, ranked_df_gRNAs_stop = get_gRNAs(ENST_ID = ENST_ID, ENST_info= ENST_info, freq_dict = freq_dict, loc2file_index= loc2file_index, loc2posType = loc2posType, dist = max_cut2ins_dist, genome_ver=config["genome_ver"], spec_score_flavor = spec_score_flavor)

            #get codon_phase information for current ENST
            file_parts_list = ENST_PhaseInCodon_index[ENST_ID]
            ENST_PhaseInCodon = {}
            for part in file_parts_list:
                Codon_phase_dict = read_pickle_files(os.path.join("genome_files","parsed_gff3", config['genome_ver'],"ENST_codonPhases",f"ENST_codonPhase_part{str(part)}.pickle"))
                ENST_PhaseInCodon = deepmerge(ENST_PhaseInCodon,Codon_phase_dict) #merge file parts if ENST codon info is split among fileparts


            ##################################
            #best start gRNA and HDR template#
            ##################################
            if target_terminus=="all" or target_terminus=="N":
                csvout_N.write(ENST_ID)
                if ranked_df_gRNAs_ATG.empty == True:
                    start_info.failed.append(ENST_ID)
                    csvout_N.write(",,,,,\n")
                    csvout_res.write(f"{ENST_ID},ERROR: no suitable gRNAs found\n")

                for i in range(0,min([gRNA_num_out, ranked_df_gRNAs_ATG.shape[0]])):
                    # if best_start_gRNA.shape[0] > 1: # multiple best scoring gRNA
                    #     best_start_gRNA = best_start_gRNA[best_start_gRNA["CSS"] == best_start_gRNA["CSS"].max()] # break the tie by CSS score
                    #     best_start_gRNA = best_start_gRNA.head(1) #get the first row in case of ties
                    current_gRNA = ranked_df_gRNAs_ATG.iloc[[i]]

                    #get HDR template
                    try:
                        HDR_template = get_HDR_template(df = current_gRNA, ENST_info = ENST_info, type = "start", ENST_PhaseInCodon = ENST_PhaseInCodon, loc2posType = loc2posType, genome_ver=config["genome_ver"],
                                                    HDR_arm_len=HDR_arm_len, tag = config["Npayload"],  ssDNA_max_size = ssDNA_max_size, Donor_type = config["Donor_type"] ,Strand_choice= config['Strand_choice'],
                                                    recoding_args = recoding_args, syn_check_args = syn_check_args)
                    except:
                        pass

                    # append the best gRNA to the final df
                    if i==0:
                        best_start_gRNAs = pd.concat([best_start_gRNAs, current_gRNA])

                    #append cfd score to list for plotting
                    pre_recoding_cfd_score = HDR_template.pre_recoding_cfd_score
                    cfd1 = ""
                    if hasattr(HDR_template,"cdf_score_post_mut_ins"):
                        cfd1 = HDR_template.cdf_score_post_mut_ins
                    if not hasattr(HDR_template,"cdf_score_post_mut2"):
                        cfd2 = cfd1
                    else:
                        cfd2 = HDR_template.cdf_score_post_mut2
                    if not hasattr(HDR_template,"cdf_score_post_mut3"):
                        cfd3 = cfd2
                    else:
                        cfd3 = HDR_template.cdf_score_post_mut3
                    if not hasattr(HDR_template,"cdf_score_post_mut4"):
                        cfd4 = cfd3
                    else:
                        cfd4 = HDR_template.cdf_score_post_mut4
                    start_info.cfd4.append(cfd4)
                    cfd_scan = 0
                    cfd_scan_no_recode = 0
                    if hasattr(HDR_template,"cdf_score_highest_in_win_scan"):
                        cfd_scan = HDR_template.cdf_score_highest_in_win_scan
                        cfd_scan_no_recode = HDR_template.scan_highest_cfd_no_recode

                    cfdfinal = HDR_template.final_cfd

                    #write csv
                    spec_score, seq, pam, s, e, cut2ins_dist, spec_weight, dist_weight, pos_weight, final_weight = get_res(current_gRNA)
                    donor = HDR_template.Donor_final
                    gRNA_cut_pos = HDR_template.CutPos # InsPos is the first letter of stop codon "T"AA or the last letter of the start codon AT"G"
                    insert_pos = HDR_template.InsPos
                    if config["recoding_off"]:
                        csvout_N.write(f",{cfd1},{cfd2},{cfd3},{cfd4},{cfd_scan},{cfd_scan_no_recode},{cfdfinal}\n")
                        csvout_res.write(f"{row_prefix},N,{seq},{pam},{s},{e},{gRNA_cut_pos},{insert_pos},{cut2ins_dist},{spec_score},{ret_six_dec(spec_weight)},{ret_six_dec(dist_weight)},{ret_six_dec(pos_weight)},{ret_six_dec(final_weight)},{ret_six_dec(pre_recoding_cfd_score)},recoding turned off,,{ret_six_dec(cfdfinal)},{donor},{HDR_template.effective_HA_len},{HDR_template.synFlags}\n")
                        csvout_res2.write(config["genome_ver"] + f",{HDR_template.ENST_chr},{insert_pos}\n")
                    else:
                        csvout_N.write(f",{cfd1},{cfd2},{cfd3},{cfd4},{cfd_scan},{cfd_scan_no_recode},{cfdfinal}\n")
                        if not isinstance(cfd4, float):
                            cfd4=""
                        csvout_res.write(f"{row_prefix},N,{seq},{pam},{s},{e},{gRNA_cut_pos},{insert_pos},{cut2ins_dist},{spec_score},{ret_six_dec(spec_weight)},{ret_six_dec(dist_weight)},{ret_six_dec(pos_weight)},{ret_six_dec(final_weight)},{ret_six_dec(pre_recoding_cfd_score)},{ret_six_dec(cfd4)},{ret_six_dec(cfd_scan)},{ret_six_dec(cfdfinal)},{donor},{HDR_template.effective_HA_len},{HDR_template.synFlags}\n")
                        csvout_res2.write(config["genome_ver"] + f",{HDR_template.ENST_chr},{insert_pos}\n")

                    #write log
                    this_log = f"{HDR_template.info}{HDR_template.info_arm}{HDR_template.info_p1}{HDR_template.info_p2}{HDR_template.info_p3}{HDR_template.info_p4}{HDR_template.info_p5}--------------------final CFD:{ret_six_dec(HDR_template.final_cfd)}\n    donor before any recoding:{HDR_template.Donor_vanillia}\n     donor after all recoding:{HDR_template.Donor_postMut}\ndonor centered(if applicable):{HDR_template.Donor_final}\n          donor (best strand):{HDR_template.Donor_final}\n\n"
                    recut_CFD_all.write(this_log)
                    if HDR_template.final_cfd > 0.03:
                        recut_CFD_fail.write(this_log)

                    if hasattr(HDR_template,"info_phase4_5UTR"):
                        fiveUTR_log.write(f"phase4_UTR\t{HDR_template.info_phase4_5UTR[0]}\t{HDR_template.info_phase4_5UTR[1]}\n")
                    if hasattr(HDR_template,"info_phase5_5UTR"):
                        fiveUTR_log.write(f"phase5_UTR\t{HDR_template.info_phase5_5UTR[0]}\t{HDR_template.info_phase5_5UTR[1]}\n")

            #################################
            #best stop gRNA and HDR template#
            #################################
            if target_terminus=="all" or target_terminus=="C":
                csvout_C.write(ENST_ID)
                if ranked_df_gRNAs_stop.empty == True:
                    stop_info.failed.append(ENST_ID)
                    csvout_C.write(",,,,,\n")
                    csvout_res.write(f"{ENST_ID},ERROR: no suitable gRNAs found\n")

                for i in range(0,min([gRNA_num_out, ranked_df_gRNAs_stop.shape[0]])):
                    # if best_stop_gRNA.shape[0] > 1: # multiple best scoring gRNA
                    #     best_stop_gRNA = best_stop_gRNA[best_stop_gRNA["CSS"] == best_stop_gRNA["CSS"].max()] # break the tie by CSS score
                    #     best_stop_gRNA = best_stop_gRNA.head(1) #get the first row in case of ties
                    current_gRNA = ranked_df_gRNAs_stop.iloc[[i]]

                    #get HDR template
                    try:
                        HDR_template = get_HDR_template(df=current_gRNA, ENST_info=ENST_info, type="stop", ENST_PhaseInCodon = ENST_PhaseInCodon, loc2posType = loc2posType,
                                                    HDR_arm_len = HDR_arm_len, genome_ver=config["genome_ver"], tag = config["Cpayload"], Donor_type = config["Donor_type"] ,Strand_choice= config['Strand_choice'], ssDNA_max_size = ssDNA_max_size,
                                                    recoding_args = recoding_args, syn_check_args = syn_check_args)
                    except:
                        pass

                    # append the best gRNA to the final df
                    best_stop_gRNAs = pd.concat([best_stop_gRNAs, current_gRNA])

                    #append cfd score to list for plotting
                    pre_recoding_cfd_score = HDR_template.pre_recoding_cfd_score
                    cfd1 = ""
                    if hasattr(HDR_template,"cdf_score_post_mut_ins"):
                        cfd1 = HDR_template.cdf_score_post_mut_ins
                    if not hasattr(HDR_template,"cdf_score_post_mut2"):
                        cfd2 = cfd1
                    else:
                        cfd2 = HDR_template.cdf_score_post_mut2
                    if not hasattr(HDR_template,"cdf_score_post_mut3"):
                        cfd3 = cfd2
                    else:
                        cfd3 = HDR_template.cdf_score_post_mut3
                    if not hasattr(HDR_template,"cdf_score_post_mut4"):
                        cfd4 = cfd3
                    else:
                        cfd4 = HDR_template.cdf_score_post_mut4

                    cfd_scan = 0
                    cfd_scan_no_recode = 0
                    if hasattr(HDR_template,"cdf_score_highest_in_win_scan"):
                        cfd_scan = HDR_template.cdf_score_highest_in_win_scan
                        cfd_scan_no_recode = HDR_template.scan_highest_cfd_no_recode

                    cfdfinal = HDR_template.final_cfd

                    #write csv
                    spec_score, seq, pam, s, e, cut2ins_dist, spec_weight, dist_weight, pos_weight, final_weight = get_res(current_gRNA)
                    donor = HDR_template.Donor_final
                    gRNA_cut_pos = HDR_template.CutPos # InsPos is the first letter of stop codon "T"AA or the last letter of the start codon AT"G"
                    insert_pos = HDR_template.InsPos
                    if config["recoding_off"]:
                        csvout_C.write(f",{cfd1},{cfd2},{cfd3},{cfd4},{cfd_scan},{cfd_scan_no_recode},{cfdfinal}\n")
                        csvout_res.write(f"{row_prefix},C,{seq},{pam},{s},{e},{gRNA_cut_pos},{insert_pos},{cut2ins_dist},{spec_score},{ret_six_dec(spec_weight)},{ret_six_dec(dist_weight)},{ret_six_dec(pos_weight)},{ret_six_dec(final_weight)},{ret_six_dec(pre_recoding_cfd_score)},recoding turned off,,{ret_six_dec(cfdfinal)},{donor},{HDR_template.effective_HA_len},{HDR_template.synFlags}\n")
                        csvout_res2.write(config["genome_ver"] + f",{HDR_template.ENST_chr},{insert_pos}\n")
                    else:
                        csvout_C.write(f",{cfd1},{cfd2},{cfd3},{cfd4},{cfd_scan},{cfd_scan_no_recode},{cfdfinal}\n")
                        if not isinstance(cfd4, float):
                            cfd4=""
                        csvout_res.write(f"{row_prefix},C,{seq},{pam},{s},{e},{gRNA_cut_pos},{insert_pos},{cut2ins_dist},{spec_score},{ret_six_dec(spec_weight)},{ret_six_dec(dist_weight)},{ret_six_dec(pos_weight)},{ret_six_dec(final_weight)},{ret_six_dec(pre_recoding_cfd_score)},{ret_six_dec(cfd4)},{ret_six_dec(cfd_scan)},{ret_six_dec(cfdfinal)},{donor},{HDR_template.effective_HA_len},{HDR_template.synFlags}\n")
                        csvout_res2.write(config["genome_ver"] + f",{HDR_template.ENST_chr},{insert_pos}\n")
                        #print(f"{row_prefix},C,{seq},{pam},{s},{e},{cut2ins_dist},{spec_score},{spec_weight:.6f},{dist_weight:.6f},{pos_weight:.6f},{final_weight:.6f},{cfd4},{cfd_scan},{cfdfinal},{donor},{HDR_template.effective_HA_len}\n")

                    #write log
                    this_log = f"{HDR_template.info}{HDR_template.info_arm}{HDR_template.info_p1}{HDR_template.info_p2}{HDR_template.info_p3}{HDR_template.info_p4}{HDR_template.info_p5}--------------------final CFD:{ret_six_dec(HDR_template.final_cfd)}\n   donor before any recoding:{HDR_template.Donor_vanillia}\n    donor after all recoding:{HDR_template.Donor_postMut}\n             donor centered:{HDR_template.Donor_final}\ndonor centered (best strand):{HDR_template.Donor_final}\n\n"
                    recut_CFD_all.write(this_log)
                    if HDR_template.final_cfd > 0.03:
                        recut_CFD_fail.write(this_log)

                    if hasattr(HDR_template,"info_phase4_5UTR"):
                        fiveUTR_log.write(f"phase4_UTR\t{HDR_template.info_phase4_5UTR[0]}\t{HDR_template.info_phase4_5UTR[1]}\n")
                    if hasattr(HDR_template,"info_phase5_5UTR"):
                        fiveUTR_log.write(f"phase5_UTR\t{HDR_template.info_phase5_5UTR[0]}\t{HDR_template.info_phase5_5UTR[1]}\n")


            protein_coding_transcripts_count +=1
            # else:
            #     log.info(f"skipping {ENST_ID} transcript type: {transcript_type} b/c transcript is not protein_coding")
            transcript_count +=1
            #report progress
            if protein_coding_transcripts_count%10 == 0 and protein_coding_transcripts_count != 0:
                endtime = datetime.datetime.now()
                elapsed_sec = endtime - starttime
                elapsed_min = elapsed_sec.seconds / 60
                log.info(f"processed {protein_coding_transcripts_count}/{transcript_count} transcripts, elapsed time {elapsed_min:.2f} min ({elapsed_sec} sec)")
                #gnum = len(gRNA_out_of_arms["start"]) + len(gRNA_out_of_arms["stop"])
                #log.info(f"number of gRNAs outside the HDR arm: {gnum}")

            ################
            #early stopping#
            ################
            # if ENST_ID == "ENST00000562221":
            #     sys.exit()
            # num_to_process = 100
            # if protein_coding_transcripts_count >=num_to_process:
            #     break

        #write csv out
        endtime = datetime.datetime.now()
        elapsed_sec = endtime - starttime
        elapsed_min = elapsed_sec.seconds / 60

        if 'num_to_process' in locals():
            pass
        else:
            num_to_process = "all"

        # write best gRNA dfs to file
        outdir = "pickles"
        mkdir(outdir)
        with open(f"{outdir}/best_start_gRNAs_of_{num_to_process}_genes.pickle", 'wb') as handle:
            pickle.dump(best_start_gRNAs, handle, protocol=pickle.HIGHEST_PROTOCOL)
        with open(f"{outdir}/best_stop_gRNAs_of_{num_to_process}_genes.pickle", 'wb') as handle:
            pickle.dump(best_stop_gRNAs, handle, protocol=pickle.HIGHEST_PROTOCOL)

        # write failed ENSTs to file
        with open(f"{outdir}/start_failed_IDs_of_{num_to_process}_genes.pickle", 'wb') as handle:
            pickle.dump(start_info.failed, handle, protocol=pickle.HIGHEST_PROTOCOL)
        with open(f"{outdir}/stop_failed_IDs_of_{num_to_process}_genes.pickle", 'wb') as handle:
            pickle.dump(stop_info.failed, handle, protocol=pickle.HIGHEST_PROTOCOL)

        # write ENSTs (whose gRNA is outside of the default HDR arm) to file
        #with open(f"pickles/gRNA_out_of_arms_{num_to_process}_genes.pickle", 'wb') as handle:
        #    pickle.dump(gRNA_out_of_arms, handle, protocol=pickle.HIGHEST_PROTOCOL)

        log.info(f"finished in {elapsed_min:.2f} min ({elapsed_sec} sec) , processed {protein_coding_transcripts_count}/{transcript_count} transcripts\nnonprotein-coding transcripts were skipped")

        recut_CFD_all.close()
        recut_CFD_fail.close()
        csvout_N.close()
        csvout_C.close()

    except Exception as e:
        print("Unexpected error:", str(sys.exc_info()))
        traceback.print_exc()
        print("additional information:", e)
        PrintException()

##########################
## function definitions ##
##########################
def ret_six_dec(myvar):
    '''
    retain six decimal points for printout
    '''

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
    '''
    merge two dictionary at the secondary key level
    '''
    if len(dict1) == 0:
        return dict2
    if len(dict2) == 0:
        return dict1
    #start merging
    dictm = dict1
    for k in dict2:
        if k in dict1: #shared key
            for k2 in dict2[k].keys():
                dictm[k][k2] = dict2[k][k2]
        else:
            dictm[k]=dict2[k]
    return dictm



class info:
    '''
    info log class
    '''
    def __init__(self)->None:
        self.cfd1 = []
        self.cfd2 = []
        self.cfd3 = []
        self.cfd4 = []
        self.cfdfinal = []
        self.failed = []

def get_res(best_start_gRNA):
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
    return([spec_score,seq,pam,s,e,cut2ins_dist,spec_weight,dist_weight,pos_weight,final_weight])

def test_memory(n):
    '''
    try allocate n MB of memory
    return true if can, and false otherwise
    '''
    try:
        x = bytearray(1024*1000*n)
        del(x)
        return True
    except:
        return False

def PrintException():
    exc_type, exc_obj, tb = sys.exc_info()
    f = tb.tb_frame
    lineno = tb.tb_lineno
    filename = f.f_code.co_filename
    linecache.checkcache(filename)
    line = linecache.getline(filename, lineno, f.f_globals)
    print('EXCEPTION IN ({}, LINE {} "{}"): {}'.format(filename, lineno, line.strip(), exc_obj))

if __name__ == "__main__": main(outdir)




