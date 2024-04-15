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
    config = parser.parse_args()
    if len(sys.argv) == 1:  # print help message if arguments are not valid
        parser.print_help()
        sys.exit(1)
    return config

def download_from_gdrive(
    id, name, dest_folder=".", overwrite=False, unzip=False
):

    try:
        os.makedirs(dest_folder)
    except Exception:
        pass

    url = f'https://drive.google.com/uc?id={id}'
    output_path = os.path.join(dest_folder, name)
    if overwrite or not os.path.isfile(output_path):
        print(f"Downloading file {output_path} as it does not exist yet.")
        gdown.download(url, output_path, quiet=False)

        if unzip:
            print(f"Unzipping file {output_path}...")
            tFile = tarfile.open(output_path)
            tFile.extractall(dest_folder)
            tFile.close()
            # zip_ref = zipfile.ZipFile(output_path, 'r')
            # zip_ref.extractall(dest_folder)
            # zip_ref.close()
            # os.remove(output_path)
  
        return output_path
    else:
        print(f"Not downloading file {output_path} as it already exists.")
        return None


if __name__ == '__main__':
    if not os.path.exists("protoSpaceJAM"):
        sys.exit("Please run this script from the repo's root directory.")
    config = vars(parse_args())

    id = "1BNq9j2yBW-a1oCcKZkNVSDFil8_KpQyt"
    if config["SpCas9_only"]:
            id="1AunvrrNNJydBYQVI9eyniMXaKULDCF8K"

    download_from_gdrive(
        id=id,
        dest_folder=os.path.join("protoSpaceJAM"),
        name="precomputed_results.tar.gz",
        unzip=True,
        overwrite=True
    )
