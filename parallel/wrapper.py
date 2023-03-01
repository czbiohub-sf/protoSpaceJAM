import os.path
import subprocess
import datetime
import traceback
import time
import os
import shutil
import argparse

import pandas as pd
import multiprocessing

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write("error: %s\n" % message)
        self.print_help()
        sys.exit(2)

def parse_args():
    parser = MyParser(description="parallel wrapper for protoSpaceJAM, all arguments are the same as the main script (and passed to the main script), except for --path2csv, and --num_workers, which is required")
    parser.add_argument(
        "--num_workers",
        default=1,
        type=int,
        help="number of parallel worker processes",
        metavar="",
    )
    parser.add_argument(
        "--genome_ver",
        default="GRCh38",
        type=str,
        help="currently supports three genomes: GRCh38, GRCm39, GRCz11, ",
        metavar="",
    )
    parser.add_argument(
        "--path2csv",
        default="input/mart_export_canonical_proteincoding.csv",
        type=str,
        help="path to a csv file containing ENST information\n *required columns*: Ensemble_ID",
        metavar="",
    )
    # gRNA
    parser.add_argument(
        "--num_gRNA_per_term",
        default=1,
        type=int,
        help="payload at the N terminus",
        metavar="",
    )
    # donor
    parser.add_argument(
        "--HA_len",
        default=500,
        help="length of the homology arm (on each side), will be the final arm length for dsDNA donors",
        type=int,
        metavar="",
    )
    parser.add_argument(
        "--Donor_type",
        default="ssODN",
        help="ssODN(default) or dsDNA",
        type=str,
        metavar="",
    )
    parser.add_argument(
        "--Strand_choice",
        default="NonTargetStrand",
        help="only applies when --Donor_type is set to ssODN, possible values are auto,TargetStrand,NonTargetStrand,CodingStrand,NonCodingStrand",
        type=str,
        metavar="",
    )
    parser.add_argument(
        "--ssODN_max_size",
        type=int,
        help="only applies when --Donor_type is set to ssODN. Enforce a length restraint of the donor (both arms + payload), setting this option will center the ssODN with respect to the payload and the recoded region",
        metavar="",
    )
    parser.add_argument(
        "--CheckEnzymes",
        default="",
        help="Restriction enzyme sites to check, separated by |, for example: BsaI|EcoRI",
        type=str,
        metavar="",
    )
    parser.add_argument(
        "--CustomSeq2Avoid",
        default="",
        help="custom sequences to avoid, separated by |",
        type=str,
        metavar="",
    )
    parser.add_argument(
        "--MinArmLenPostTrim",
        default=0,
        help="Minimum length the homology arm after trimming,  default is 0 (turning off trimming)",
        type=str,
        metavar="",
    )

    # payload
    parser.add_argument(
        "--payload",
        default="",
        type=str,
        help="payload, overrides --Npayloadf and --Cpayload, --Tag, --Linker",
        metavar="",
    )
    parser.add_argument(
        "--Npayload",
        default="",
        type=str,
        help="payload at the N terminus, overrides --Tag --Linker",
        metavar="",
    )
    parser.add_argument(
        "--Cpayload",
        default="",
        type=str,
        help="payload at the N terminus, overrides --Tag --Linker",
        metavar="",
    )
    parser.add_argument(
        "--Tag",
        default="ACCGAGCTCAACTTCAAGGAGTGGCAAAAGGCCTTTACCGATATGATG",
        type=str,
        help="default is the mNG11 tag",
        metavar="",
    )
    parser.add_argument(
        "--Linker",
        default="GGTGGCGGATTGGAAGTTTTGTTTCAAGGTCCAGGAAGTGGT",
        type=str,
        help="default is the GS linker",
        metavar="",
    )

    # recoding
    parser.add_argument(
        "--recoding_off",
        default=False,
        action="store_true",
        help="turn off *all* recoding",
    )
    parser.add_argument(
        "--recoding_stop_recut_only",
        default=False,
        action="store_true",
        help="use recoding to prevent recut",
    )
    parser.add_argument(
        "--recoding_full",
        default=False,
        action="store_true",
        help="use recoding to prevent recut + recode region between insert and cut site",
    )
    parser.add_argument(
        "--cfdThres",
        default=0.03,
        help="protoSpaceJAM will attempt to lower the recut cfd to this threshold (by recoding), cfd values lower than the threshold will be considered not suceptible to being recut anymore.",
    )
    parser.add_argument(
        "--recode_order",
        default="PAM_first",
        help="possible values: protospacer_first, PAM_first",
    )

    # output
    parser.add_argument("--outdir", default="logs", type=str, help="output directory")

    config = parser.parse_args()
    return config


def pwrapper():
    try:
        starttime = datetime.datetime.now()

        config = vars(parse_args())

        ##############################################################
        #split the csv file into chunks equal to the number of workers
        input_file_basename = os.path.basename(config["path2csv"])
        output_dir = config["outdir"]
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        else:
            shutil.rmtree(output_dir)
            os.makedirs(output_dir)
        chunk_size = len(pd.read_csv(config["path2csv"])) // (config["num_workers"] - 1)
        split_csv_by_rows(input_file_path = config["path2csv"], chunk_size = chunk_size, output_file_prefix = os.path.join(output_dir, input_file_basename + "_chunk_"))

        ##################################
        #run the main script on each chunk
        processes = []
        parent_conns = []
        events = []
        results = []

        for i in range(config["num_workers"]):
            parent_conn, child_conn = multiprocessing.Pipe()
            event = multiprocessing.Event()
            process = multiprocessing.Process(target=worker, args=(i, event, child_conn, config))
            processes.append(process)
            parent_conns.append(parent_conn)
            events.append(event)
            process.start()
            #print(f"process {i+1} started")

        for i in range(config["num_workers"]):
            events[i].wait()
            results.append(parent_conns[i].recv())
            processes[i].join()

        ##################
        #merge the results
        with open(os.path.join(output_dir, "result.csv"), "w") as outfile:
            for i in range(config["num_workers"]):
                chunk = i + 1
                chunk_res = os.path.join(output_dir, "outdir_chunk_" + str(chunk), "result.csv")
                with open(chunk_res, "r") as infile:
                    if i != 0:
                        infile.readline() #skip header
                    for line in infile:
                        outfile.write(line)

        with open(os.path.join(output_dir, "out_Cterm_recut_cfd.csv"), "w") as outfile:
            for i in range(config["num_workers"]):
                chunk = i + 1
                chunk_res = os.path.join(output_dir, "outdir_chunk_" + str(chunk), "out_Cterm_recut_cfd.csv")
                with open(chunk_res, "r") as infile:
                    if i != 0:
                        infile.readline() #skip header
                    for line in infile:
                        outfile.write(line)

        with open(os.path.join(output_dir, "out_Nterm_recut_cfd.csv"), "w") as outfile:
            for i in range(config["num_workers"]):
                chunk = i + 1
                chunk_res = os.path.join(output_dir, "outdir_chunk_" + str(chunk), "out_Nterm_recut_cfd.csv")
                with open(chunk_res, "r") as infile:
                    if i != 0:
                        infile.readline() #skip header
                    for line in infile:
                        outfile.write(line)


        elapsed = cal_elapsed_time(starttime, datetime.datetime.now())
        print(f"finished in {elapsed[0]:.2f} min ({elapsed[1]} sec)")


    except:
        traceback.print_exc()
        raise

def split_csv_by_rows(input_file_path, chunk_size, output_file_prefix):
    """
    Splits a CSV file into smaller files with a specified number of rows per file.
    :param input_file_path: The path to the input CSV file.
    :param chunk_size: The number of rows to include in each output file.
    :param output_file_prefix: The prefix to use for the output file names.
    """
    # Load the input CSV file into a Pandas DataFrame.
    input_df = pd.read_csv(input_file_path)

    # Determine the total number of chunks to create.
    total_chunks = (len(input_df) + chunk_size - 1) // chunk_size

    # Split the input DataFrame into chunks and write each chunk to a separate CSV file.
    for i in range(total_chunks):
        output_file_path = f"{output_file_prefix}{i+1}.csv"
        start_index = i * chunk_size
        end_index = start_index + chunk_size
        chunk_df = input_df[start_index:end_index]
        chunk_df.to_csv(output_file_path, index=False)

def worker(num, event, conn, workerConfig):
    """Worker function"""
    print(f"Worker {num+1} is running on chunk {num+1}")

    workerConfig = workerConfig.copy()
    output_dir = workerConfig["outdir"]
    input_file_basename = os.path.basename(workerConfig["path2csv"])
    del workerConfig["num_workers"]

    #no input file or input file contains less than 2 rows
    if not os.path.isfile(workerConfig["path2csv"]) or len(pd.read_csv(workerConfig["path2csv"])) < 2:
        result = 0
        conn.send(result)
        conn.close()
        event.set()


    workerConfig["path2csv"] = os.path.join(output_dir, input_file_basename + "_chunk_" + str(num+1) + ".csv")
    workerConfig["outdir"] = os.path.join(output_dir, "outdir_chunk_" + str(num+1))
    if not os.path.exists(workerConfig["outdir"]):
        os.makedirs(workerConfig["outdir"])

    #constrcut the command line arguments
    cmd = ["python", "main.py"]
    for key, val in workerConfig.items():

        if key in ["recoding_full", "recoding_stop_recut_only", "recoding_off"]:
            if val == "True" or val == True:
                cmd.append("--recoding_full")
            elif val == "True" or val == True:
                cmd.append("--recoding_stop_recut_only")
            elif val == "True" or val == True:
                cmd.append("--recoding_off")
            continue

        if val != "" and val is not None and val != 'None':
            cmd.append(f"--{key}")
            cmd.append(str(val))


    #print(cmd)
    with open(os.path.join(workerConfig["outdir"], "stdout.txt"), "w") as stdout, open(os.path.join(workerConfig["outdir"], "stderr.txt"), "w") as stderr:
        subprocess.call(cmd, stdout=stdout, stderr=stderr)

    print(f"Worker {num+1} finished on chunk {num+1}")
    result = 0
    conn.send(result)
    conn.close()
    event.set()

def cal_elapsed_time(starttime, endtime):
    """
    output: [elapsed_min,elapsed_sec]
    """
    elapsed_sec = endtime - starttime
    elapsed_min = elapsed_sec.seconds / 60
    return [elapsed_min, elapsed_sec]

if __name__ == "__main__":
    pwrapper()
