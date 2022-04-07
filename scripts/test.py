
import fetch_ensembl
from utils import *
from gRNA_search import *
from Bio.Seq import Seq
import csv

logging.setLoggerClass(ColoredLogger)
#logging.basicConfig()
log = logging.getLogger()
log.setLevel(logging.INFO) #set the level of warning displayed
###################################
# search gRNA in the transcript   #
###################################
#my_transcript = fetch_ensembl.fetch_ensembl_transcript(ensembl_transcript_id ="ENST00000620157", exon_annot = True) # reverse strand
my_transcript = fetch_ensembl.fetch_ensembl_transcript(ensembl_transcript_id ="ENST00000252108", exon_annot = True) # forward strand
#print(my_transcript.features)

cdsList =  "all"  # all cds
#cdsList = [2,3] # 2nd and 3rd cds
cds_flank_len = 3 # length of flanking intro sequences added so that the gRNA can cut at the edge of exons
protosp_len = 20
PAM = "NGG"
HDR_flank_len = 50
adjust_frame = True
res_gRNA_list = []

#constructing the list of cds to process
if cdsList == "all": #update cdsList with all
    tmp_list = [feat for feat in my_transcript.features if feat.type == 'cds']
    cdsList = [*range(0, len(tmp_list))]
elif all(isinstance(x, int) for x in cdsList) :
    pass
else:
    log.error(f"invalid value for cdsList, acceptable values: (1) a list of integers or (2)\"all\"\n Offending value: {cdsList}")
    quit()

#processing cds sequentially
for current_cds in cdsList:

    which_cds = current_cds + 1 # which_cds is 1 indexed

    print(f"working on cds #{which_cds}")

    #get exon sequence from transcript object
    cds_w_flank = get_exon_concat_with_flank(transcriptObj = my_transcript,
                                             which_cds = which_cds,
                                             cds_flank_len = cds_flank_len)

    # search gRNA and update scores
    gRNA_list = search_gRNA(protosp_len=protosp_len,
                            PAM=PAM,search_in=cds_w_flank)

    #update gRNA objects with cut site
    gRNA_list = get_cutsite_in_gene(transcriptObj=my_transcript,
                                    listOfgRNAObj=gRNA_list,
                                    which_cds=which_cds,
                                    cds_flank_len = cds_flank_len)

    #adjust frame
    if (adjust_frame == True):
        gRNA_list = nudge_cutsite_inframe(transcriptObj=my_transcript,
                                          listOfgRNAObj=gRNA_list,
                                          HDR_flank_len=HDR_flank_len,
                                          which_cds=which_cds)

    #update gRNA objects with flanks
    gRNA_list = get_HDR_flank(transcriptObj=my_transcript,
                              listOfgRNAObj=gRNA_list,
                              HDR_flank_len=HDR_flank_len)

    for gRNA in gRNA_list:
        print(f"{gRNA.protospacer} {gRNA.pam} {gRNA.g_strand} {gRNA.g_st} {gRNA.g_en}")

    res_gRNA_list.extend(gRNA_list)



#write gRNAs to csv
attrNames = [attr for attr in dir(res_gRNA_list[0]) if not callable(getattr(res_gRNA_list[0], attr)) and not attr.startswith("__")]
with open("gRNA.csv",'w', newline="") as resultFile:
    wr = csv.writer(resultFile, dialect='excel')
    wr.writerow(attrNames)
    for gRNA in res_gRNA_list:
        attr2write = [getattr(gRNA,attr) for attr in attrNames]
        wr.writerow(attr2write)



gRNA_list = search_gRNA(protosp_len=protosp_len,
                        PAM=PAM, search_in="TTTTTTTTTTAAAAAAAAAAAAAAAAAAAAAGGCCAAAAAAAAAAAAAAAAAAAAATTTTTTTTTT", flanksize=11)

#test getting flank
a = get_rightflank(seq = "123456789", pos=8, flanksize=4)
print(a)

a = get_leftflank(seq = "123456789", pos=8, flanksize=4)
print(a)