import os
import gzip
import re
import sys
import argparse

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write("error: %s\n" % message)
        self.print_help()
        sys.exit(2)

def parse_args():
    parser = MyParser(description='')
    parser.add_argument('--pam', default="", type=str, help='path to the gzfile', metavar='')
    parser.add_argument('--infile', default="", type=str, help='path to the dir containing gzfile', metavar='')
    parser.add_argument('--outfile', default="", type=str, help='name of genome fasta file', metavar='')
    
    config = parser.parse_args()
    if len(sys.argv)==1: # print help message if arguments are not valid
        parser.print_help()
        sys.exit(1)
    return config

config = vars(parse_args())

pam = config['pam'].split("|")[0]
infile = config['infile']
outfile = config['outfile']

def check_pam(seq, gPAM):
    seq = seq.upper()
    gPAM = gPAM.upper()
    gPAM = re.sub(
        "N", "[ACGT]", gPAM, flags=re.IGNORECASE
    )  # IUPAC nucleotide code http://www.bioinformatics.org/sms/iupac.html
    gPAM = re.sub("R", "[AG]", gPAM, flags=re.IGNORECASE)  # IUPAC nucleotide code
    gPAM = re.sub("Y", "[CT]", gPAM, flags=re.IGNORECASE)  # IUPAC nucleotide code
    gPAM = re.sub("S", "[GC]", gPAM, flags=re.IGNORECASE)  # IUPAC nucleotide code
    gPAM = re.sub("W", "[AT]", gPAM, flags=re.IGNORECASE)  # IUPAC nucleotide code
    gPAM = re.sub("K", "[GT]", gPAM, flags=re.IGNORECASE)  # IUPAC nucleotide code
    gPAM = re.sub("M", "[AC]", gPAM, flags=re.IGNORECASE)  # IUPAC nucleotide code
    gPAM = re.sub("B", "[CGT]", gPAM, flags=re.IGNORECASE)  # IUPAC nucleotide code
    gPAM = re.sub("D", "[AGT]", gPAM, flags=re.IGNORECASE)  # IUPAC nucleotide code
    gPAM = re.sub("H", "[ACT]", gPAM, flags=re.IGNORECASE)  # IUPAC nucleotide code
    gPAM = re.sub("V", "[ACG]", gPAM, flags=re.IGNORECASE)  # IUPAC nucleotide code
    if re.search(gPAM, seq):
        return True
    else:
        return False
    
def main(pam, infile, outfile):
    illegal_pams = 0
   
    print("Scanning folder {}\nUsing PAM: {}".format(infile, pam))

    if infile.endswith(".gz"):
        with gzip.open(infile, "rt") as f, gzip.open(outfile, "wt") as o:
            for line in f:
                fields = line.strip().split("\t")
                if check_pam(fields[1], pam):
                    o.write(line)
                else:
                    illegal_pams += 1
    print("Scanning complete. Found {} illegal PAMs.".format(illegal_pams))


if __name__ == "__main__":
    main(pam, infile, outfile)