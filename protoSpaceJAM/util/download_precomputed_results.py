import os
import tarfile
import gdown
import sys

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
    download_from_gdrive(
        id="1kFl22tFgCyr9Ipnm7h6Y-fu0Usn-4wsa",
        dest_folder=os.path.join("."),
        name="precomputed_results.tar.gz",
        unzip=True,
        overwrite=True
    )
