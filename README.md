![Python version](https://img.shields.io/badge/python-3.9%20%7C%203.10-blue)
![Webserver dependencies-0-brightgreen](https://user-images.githubusercontent.com/4129442/198696112-92ecc372-f3b5-4498-8cd9-4a01de0f851b.svg)
![status-active development-blueviolet](https://user-images.githubusercontent.com/4129442/198695999-a70bcd5f-c52e-4895-a1e7-d6b0da132812.svg)

# protoSpaceJAM 
A standalone program to design guide RNA and repair donors for CRISPR knock-in experiments  

## Key features:  
- Fully standalone, no calling to other bioinformatics servers
- Sophisticated guide RNA ranking system
  - Specificity weight
  - Penalize cuts near exon-intron junctions etc.
  - Penalize cuts far away from the payload insertion site
- Sophisticated DNA donor design
  - Recode to prevent recut
  - Recode to facilitate payload insertion and prevent recut
  - Center the DNA donor around the region containing the payload and recoded bases. 
  - Enforce maximum length of DNA donor
  - Automated selection of of ssDNA donor for maximum chance of payload insertion and prevent recut
  - Scan and trim hard-to-synthesis motifs (coming soon)
- Use/generate pre-computed genome-wide guide RNAs and their properties for fast runtime


## Usage

### Clone the repository
```
git clone https://github.com/czbiohub/protoSpaceJAM
```
### Go the repository directory, switch the branch if running branch other than master
```
cd protoSpaceJAM
git checkout <branch you'd like to run>
```
### Create conda environment and activate it
```
conda create -y -n protospacejam python=3.9 && conda activate protospacejam
```
### Install required packages
```
pip install -r requirements.txt
```
### Download and unzip pre-computed data
```
python ./scripts/download_precomputed_results.py
```
PS: Default dependencies will only allow users to run protoSpaceJam 
with existing pre-computed data. If you like to execute pre-computation
steps please follow [here](). TODO: fill  this link

### Run protoSpaceJAM
```
conda activate protospacejam
python main.py --path2csv input/test_protospaceX.csv --ssODN_max_size 200 --recoding_all
```

