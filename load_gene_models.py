import os
genome_ver = "hg38"
pk_file = "ENST_info.pickle"
pk_file_path = os.path.join("genome_files", "gff3", version, pk_file)

def read_pickle_files(file):
    with open(file, 'rb') as handle:
        mydict = pickle.load(handle)
    return mydict

#read dict from file
ENST_info = read_pickle_files(pk_file)
