![Python version](https://img.shields.io/badge/python-3.9%20%7C%203.10-blue)
![Webserver dependencies-0-brightgreen](https://user-images.githubusercontent.com/4129442/198696112-92ecc372-f3b5-4498-8cd9-4a01de0f851b.svg)
![status-active development-blueviolet](https://user-images.githubusercontent.com/4129442/198695999-a70bcd5f-c52e-4895-a1e7-d6b0da132812.svg)

# protoSpaceJAM 
A standalone program to design guide RNA and repair donors for CRISPR knock-in experiments  
A web-portal is available at http://protospacejam.czbiohub.org/

## Key features:  
- Fully standalone, no calling to other bioinformatics servers
- Sophisticated guide RNA ranking system ([details](https://czbiohub.github.io/protoSpaceJAM/algorithm.html#grna-scoring))
  - Specificity weight
  - Penalize cuts near exon-intron junctions etc.
  - Penalize cuts far away from the payload insertion site
- Sophisticated DNA donor design ([details](https://czbiohub.github.io/protoSpaceJAM/algorithm.html#recoding-strategy))
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
git checkout main
```
### Create conda environment and activate it
```
conda create -y -n protospacejam python=3.9 && conda activate protospacejam
```
### Install protoSpaceJAM
```
pip install .
```
### Download and unzip pre-computed data
```
python ./protoSpaceJAM/util/download_precomputed_results.py
```
PS: Default dependencies will only allow users to run protoSpaceJam 
with existing pre-computed data. If you like to execute pre-computation
steps please follow [here](https://github.com/czbiohub/protoSpaceJAM/tree/main/precomputed_gRNAs)

### Run a quick test to verify installation
```
python protoSpaceJAM/tests/run_quick_test_pJAM.py
```
A successful test will have a printout similar to `Ran 2 tests in 14.644s   OK` at the end.

### Run protoSpaceJAM
```
conda activate protospacejam
python main.py --path2csv input/test_protoSpaceJAM.csv --outdir logs/test
```

## License
Distributed under the terms of the BSD-3 license, "protoSpaceJAM" is free and open source software
