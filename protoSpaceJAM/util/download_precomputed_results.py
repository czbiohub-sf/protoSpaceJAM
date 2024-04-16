import os
import tarfile
import gdown
import sys
import argparse

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write("error: %s\n" % message)
        self.print_help()
        sys.exit(2)

def parse_args():
    parser = MyParser(description="download precomputed results\n")
    parser.add_argument(
        "--SpCas9_only",
        default=False,
        action="store_true",
        help
        ="if set, download precomputed results for SpCas9 only. Otherwise, download precomputed results for SpCas9, VQR-Cas9 and enAsCas12a.",
    )
    parser.add_argument(
        "--cleanup",
        default=False,
        action="store_true",
        help="[default: False] if set, remove the tarball file after unzipping.",
    )
    parser.add_argument(
        "--use_local_file",
        default=False,
        action="store_true",
        help="[default: False] if set, use the local file instead of downloading it from the Google Drive.\n The path to the  protoSpaceJAM/precomputed_results.tar",
    )
    parser.add_argument(
        "--local_file_path",
        default="protoSpaceJAM/precomputed_results.tar",
        help="the path to the precomputed results tarball file. [default: protoSpaceJAM/precomputed_results.tar] ",
    )

    parser.add_argument
    config = parser.parse_args()
    return config

def download_from_gdrive(
    id, name, dest_folder=".", overwrite=False, unzip=False, cleanup=False, use_local_file=False, local_file_path="protoSpaceJAM/precomputed_results.tar"
):

    try:
        os.makedirs(dest_folder)
    except Exception:
        pass

    url = f'https://drive.google.com/uc?id={id}'
    output_path = os.path.join(dest_folder, name)

    if use_local_file:
        print(f"Using local file {local_file_path}.")
        output_path = local_file_path
        if not os.path.isfile(output_path):
            sys.exit(f"Local file {output_path} does not exist.")

    elif overwrite or not os.path.isfile(output_path) :
        print(f"Downloading file {output_path} as it does not exist yet.")
        try:
            gdown.download(url, output_path, quiet=False, resume=True)
        except Exception as e:
            sys.exit(f"Failed to download file from {url} with error: {e}")

    if not os.path.isfile(output_path):
        sys.exit(f"Error reading file {output_path}, either download failed or the local file is missing.")

    if unzip:
        print(f"Unzipping file {output_path}...")
        tFile = tarfile.open(output_path)
        tFile.extractall(dest_folder)
        tFile.close()
        # zip_ref = zipfile.ZipFile(output_path, 'r')
        # zip_ref.extractall(dest_folder)
        # zip_ref.close()
        # os.remove(output_path)

    
    # remove the tarfile
    if cleanup:
        os.remove(output_path)
  
        return output_path

    return None


if __name__ == '__main__':
    if not os.path.exists("protoSpaceJAM"):
        sys.exit("Please run this script from the repo's root directory.")
    config = vars(parse_args())

    id = "1E2BlwrPjHXQgGAiZFU9YYdPIR3Ovslna"
    filename = "precomputed_results.tar"
    if config["SpCas9_only"]:
            id="13uPuSPXZTP0F8oK1DzziE9CcuonvtkFd"
            filename = "precomputed_results_SpCas9_only.tar"

    download_from_gdrive(
        id=id,
        dest_folder=os.path.join("protoSpaceJAM"),
        name=filename,
        unzip=True,
        overwrite=True,
        cleanup=config["cleanup"],
        use_local_file=config["use_local_file"],
        local_file_path=config["local_file_path"],
    )
