import argparse
import requests
import logging
import time
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry

log = logging.getLogger(__name__)
logging.basicConfig(level=logging.WARNING)

def fetch_ensembl_transcript(
    ensembl_transcript_id, 
    exon_annot = False,
    is_retry: bool = False):

    '''
    Fetch the requested Ensembl transcript.
    Get the requested Ensembl transcript, together with exon and
    coding region (CDS) boundaries.
    Parameters
    ----------
    ensembl_transcript_id : str
      the ensembl transcript id, of the form ENST...
    
    Returns
    -------
    `Bio.SeqRecord`
      The requested transcript sequence, in 5' -> 3' order, together
      with exon and CDS features. The coordinates of exons and CDS
      features are relative to the sequence fragment.
    '''

    #retry strategy
    retry_strategy = Retry(
        total=7,
        backoff_factor=2, #how long the processes will sleep = {backoff factor} * (2 ** ({number of total retries} - 1));  2 seconds = 1, 2, 4, 8, 16, 32, 64, 128, 256, 512
        status_forcelist=[429, 500, 502, 503, 504],
        method_whitelist=["HEAD", "GET", "OPTIONS"]
    )
    adapter = HTTPAdapter(max_retries=retry_strategy)
    http = requests.Session()
    http.mount("https://", adapter)
    http.mount("http://", adapter)

    base_url = "http://rest.ensembl.org"

    # First, fetch the transcript sequence
    #url = base_url + f"/sequence/id/{ensembl_transcript_id}"
    #response = requests.get(url, {"type": "genomic", "content-type": "application/json"})

    log.info(f"Querying Ensembl for sequence of {ensembl_transcript_id}")

    url = base_url + f"/sequence/id/{ensembl_transcript_id}?type=genomic"
    response = http.get(url, headers= { "content-type": "application/json" })

    try:
        response.raise_for_status()
    except requests.HTTPError:
        log.error("Ensembl sequence REST query returned error "
                  "{}".format(response.text))
        raise ValueError(response.text)

    response_data = response.json()

    try:
        description = response_data['desc'].split(':')
        species = description[1]
        try:
            chromosome_number = int(description[2])
        except:
            chromosome_number = str(description[2])
        sequence_left = int(description[3])
        sequence_right = int(description[4])
        transcript_strand = int(description[5])

        if sequence_left > sequence_right:
            raise ValueError(f"Expected left sequence boundary {sequence_left} "
                             f"<= right sequence boundary {sequence_right}: did "
                             "the format of the Ensembl REST response change?")
        
        sequence_id = response_data['id']
        
        seq_str = response_data['seq']

        log.info(f"Retrieved sequence {response_data['desc']} of length "
                 f"{sequence_right - sequence_left} for species {species} on "
                 f"strand {transcript_strand}")
    except (KeyError, ValueError) as e:
        log.error(e)
        log.error('Error parsing sequence metadata from Ensembl REST response - '
                  'did the format of the response change?')
        raise ValueError(e)
    
    seq = Seq(seq_str)
    
    record = SeqRecord(seq, id=sequence_id,
                       description=":".join(description))
    if exon_annot:
        adapter = HTTPAdapter(max_retries=retry_strategy)
        http = requests.Session()
        http.mount("https://", adapter)
        http.mount("http://", adapter)

        #url = base_url + f"/overlap/id/{ensembl_transcript_id}"
        #response = requests.get(url, {"feature": ["cds", "exon"], "content-type": "application/json"})

        log.info(f"Querying Ensembl for overlaps of {ensembl_transcript_id}")

        url = base_url + f"/overlap/id/{ensembl_transcript_id}?feature=cds&exon"
        response = http.get(url, headers = {  "content-type": "application/json" })
        try:
            response.raise_for_status()
        except requests.HTTPError:
            log.error("Ensembl sequence REST query returned error "
                      "{}".format(response.text))
            raise ValueError(response.text)

        response_data = response.json()

        try:
            # Handle the unlikely event of a single piece of information overlapping a lonely transcript

            if not hasattr(response_data, '__iter__'):
                response_data = [response_data]

            for response_datum in response_data:
                if response_datum['Parent'] != ensembl_transcript_id:
                    continue

                if response_datum['assembly_name'] != species:
                    continue

                # Store feature locations 0-indexed from the left-most sequence boundary

                record.features.append(SeqFeature(
                    location=FeatureLocation(
                        int(response_datum['start']) - sequence_left,
                        int(response_datum['end']) - sequence_left + 1,
                        strand=int(response_datum['strand'])),
                    type=response_datum['feature_type']))
            num_exon_boundaries = len([f for f in record.features
                                       if f.type == 'exon'])

            num_cds_boundaries = len([f for f in record.features
                                      if f.type == 'cds'])

            log.info(f"Retrieved {num_exon_boundaries} exons and "
                     f"{num_cds_boundaries} coding regions for transcript "
                     f"{ensembl_transcript_id}")
        except (KeyError, ValueError) as e:
            log.error(e)
            log.error('Error parsing overlap metadata from Ensembl REST response - '
                      'did the format of the response change?')
            raise ValueError(e)

    record.annotations['reference_species'] = species
    record.annotations['reference_chromosome_number'] = chromosome_number
    record.annotations['reference_left_index'] = sequence_left
    record.annotations['reference_right_index'] = sequence_right
    record.annotations['transcript_strand'] = transcript_strand

    # Finally, sort features by their start locations
    record.features.sort(key=lambda f: f.location.start)

    # Count exon numbers
    cds_list = [feat for feat in record.features if feat.type == 'cds']
    num_cds = len(cds_list)
    record.num_cds = num_cds

    return record

def fetch_ensembl_sequence(
    chromosome, 
    region_left, 
    region_right, 
    expand = 200):

    '''
    Returns genome sequence based on chromosome range. The sequence is expanded by a flat amount on both the 5 and
    3 termini. 
    '''
    
    base_url = "http://rest.ensembl.org"
    ext = f"/sequence/region/human/{chromosome}:{region_left}..{region_right}:1?expand_5prime={expand};expand_3prime={expand}"
    r = requests.get(base_url + ext, headers = {"Content-Type":"text/plain"})
    
    if not r.ok:
        r.raise_for_status()
    
    #sequence = Seq(r.text, IUPACUnambiguousDNA())
    sequence = Seq(r.text)
    return sequence