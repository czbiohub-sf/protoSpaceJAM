# 1. Calculates the Cutting Frequency Determination score
# 
#   Copied and adapted from
#   https://github.com/maximilianh/crisporWebsite/blob/master/crispor.py

#   Input: 1. 23mer WT sgRNA sequence: str
#          2. 23mer Off-target sgRNA sequence: str
#
#   Output: CFD scores (Doench score, MIT score, and Manu's score)
# 
# 2. Calculates the GC content
#     
#   Input: sequence: str
#
#   Output: GC content in %

from functools import lru_cache
import math

_specificity_weight_low = 45
_specificity_weight_high = 65
_dist_weight_variance = 55
_before_start_codon_penalty = 0.2

# Loaded from pam_scores.pkl
pam_scores = {
    'AA': 0.0,
    'AC': 0.0,
    'GT': 0.016129031999999998,
    'AG': 0.25925925899999996,
    'CC': 0.0,
    'CA': 0.0,
    'CG': 0.107142857,
    'TT': 0.0,
    'GG': 1.0,
    'GC': 0.022222222000000003,
    'AT': 0.0,
    'GA': 0.06944444400000001,
    'TG': 0.038961038999999996,
    'TA': 0.0,
    'TC': 0.0,
    'CT': 0.0
}

# Loaded from mismatch_score.pkl
mm_scores = {
    'rU:dT,12': 0.8,
    'rU:dT,13': 0.692307692,
    'rU:dC,5': 0.64,
    'rG:dA,14': 0.26666666699999997,
    'rG:dG,19': 0.448275862,
    'rG:dG,18': 0.47619047600000003,
    'rG:dG,15': 0.272727273,
    'rG:dG,14': 0.428571429,
    'rG:dG,17': 0.235294118,
    'rG:dG,16': 0.0,
    'rC:dC,20': 0.058823529000000006,
    'rG:dT,20': 0.9375,
    'rG:dG,13': 0.42105263200000004,
    'rG:dG,12': 0.529411765,
    'rU:dC,6': 0.571428571,
    'rU:dG,14': 0.28571428600000004,
    'rU:dT,18': 0.666666667,
    'rA:dG,13': 0.21052631600000002,
    'rA:dG,12': 0.263157895,
    'rA:dG,11': 0.4,
    'rA:dG,10': 0.333333333,
    'rA:dA,19': 0.538461538,
    'rA:dA,18': 0.5,
    'rA:dG,15': 0.272727273,
    'rA:dG,14': 0.214285714,
    'rA:dA,15': 0.2,
    'rA:dA,14': 0.533333333,
    'rA:dA,17': 0.133333333,
    'rA:dA,16': 0.0,
    'rA:dA,11': 0.307692308,
    'rA:dA,10': 0.882352941,
    'rA:dA,13': 0.3,
    'rA:dA,12': 0.333333333,
    'rG:dA,13': 0.3,
    'rG:dA,12': 0.384615385,
    'rG:dA,11': 0.384615385,
    'rG:dA,10': 0.8125,
    'rG:dA,17': 0.25,
    'rG:dA,16': 0.0,
    'rG:dA,15': 0.14285714300000002,
    'rG:dA,6': 0.666666667,
    'rG:dG,20': 0.428571429,
    'rG:dA,19': 0.666666667,
    'rG:dA,18': 0.666666667,
    'rU:dC,4': 0.625,
    'rG:dT,12': 0.933333333,
    'rG:dT,13': 0.923076923,
    'rU:dG,11': 0.666666667,
    'rC:dA,3': 0.6875,
    'rC:dA,2': 0.9090909090000001,
    'rC:dA,1': 1.0,
    'rC:dA,7': 0.8125,
    'rC:dA,6': 0.9285714290000001,
    'rC:dA,5': 0.636363636,
    'rC:dA,4': 0.8,
    'rC:dA,9': 0.875,
    'rC:dA,8': 0.875,
    'rU:dT,6': 0.8666666670000001,
    'rA:dG,20': 0.22727272699999998,
    'rG:dT,18': 0.692307692,
    'rU:dG,10': 0.533333333,
    'rG:dT,19': 0.7142857140000001,
    'rG:dA,20': 0.7,
    'rC:dT,20': 0.5,
    'rU:dC,2': 0.84,
    'rG:dG,10': 0.4,
    'rC:dA,17': 0.46666666700000003,
    'rC:dA,16': 0.307692308,
    'rC:dA,15': 0.066666667,
    'rC:dA,14': 0.7333333329999999,
    'rC:dA,13': 0.7,
    'rC:dA,12': 0.538461538,
    'rC:dA,11': 0.307692308,
    'rC:dA,10': 0.9411764709999999,
    'rG:dG,11': 0.428571429,
    'rU:dC,20': 0.176470588,
    'rG:dG,3': 0.384615385,
    'rC:dA,19': 0.46153846200000004,
    'rC:dA,18': 0.642857143,
    'rU:dG,17': 0.705882353,
    'rU:dG,16': 0.666666667,
    'rU:dG,15': 0.272727273,
    'rG:dG,2': 0.692307692,
    'rU:dG,13': 0.7894736840000001,
    'rU:dG,12': 0.947368421,
    'rG:dA,9': 0.533333333,
    'rG:dA,8': 0.625,
    'rG:dA,7': 0.571428571,
    'rG:dG,5': 0.7857142859999999,
    'rG:dA,5': 0.3,
    'rG:dA,4': 0.363636364,
    'rG:dA,3': 0.5,
    'rG:dA,2': 0.636363636,
    'rG:dA,1': 1.0,
    'rG:dG,4': 0.529411765,
    'rG:dG,1': 0.7142857140000001,
    'rA:dC,9': 0.666666667,
    'rG:dG,7': 0.6875,
    'rG:dT,5': 0.8666666670000001,
    'rU:dT,20': 0.5625,
    'rC:dC,15': 0.05,
    'rC:dC,14': 0.0,
    'rC:dC,17': 0.058823529000000006,
    'rC:dC,16': 0.153846154,
    'rC:dC,11': 0.25,
    'rC:dC,10': 0.38888888899999996,
    'rC:dC,13': 0.13636363599999998,
    'rC:dC,12': 0.444444444,
    'rC:dA,20': 0.3,
    'rC:dC,19': 0.125,
    'rC:dC,18': 0.133333333,
    'rA:dA,1': 1.0,
    'rA:dA,3': 0.705882353,
    'rA:dA,2': 0.727272727,
    'rA:dA,5': 0.363636364,
    'rA:dA,4': 0.636363636,
    'rA:dA,7': 0.4375,
    'rA:dA,6': 0.7142857140000001,
    'rA:dA,9': 0.6,
    'rA:dA,8': 0.428571429,
    'rU:dG,20': 0.090909091,
    'rC:dC,9': 0.6190476189999999,
    'rC:dC,8': 0.642857143,
    'rU:dT,10': 0.857142857,
    'rU:dT,11': 0.75,
    'rU:dT,16': 0.9090909090000001,
    'rU:dT,17': 0.533333333,
    'rU:dT,14': 0.6190476189999999,
    'rU:dT,15': 0.578947368,
    'rC:dC,1': 0.913043478,
    'rU:dT,3': 0.7142857140000001,
    'rC:dC,3': 0.5,
    'rC:dC,2': 0.695652174,
    'rC:dC,5': 0.6,
    'rC:dC,4': 0.5,
    'rC:dC,7': 0.470588235,
    'rC:dC,6': 0.5,
    'rU:dT,4': 0.47619047600000003,
    'rU:dT,8': 0.8,
    'rU:dT,9': 0.9285714290000001,
    'rA:dC,19': 0.375,
    'rA:dC,18': 0.4,
    'rA:dC,17': 0.176470588,
    'rA:dC,16': 0.192307692,
    'rA:dC,15': 0.65,
    'rA:dC,14': 0.46666666700000003,
    'rA:dC,13': 0.6521739129999999,
    'rA:dC,12': 0.7222222220000001,
    'rA:dC,11': 0.65,
    'rA:dC,10': 0.5555555560000001,
    'rU:dC,7': 0.588235294,
    'rC:dT,8': 0.65,
    'rC:dT,9': 0.857142857,
    'rC:dT,6': 0.9285714290000001,
    'rC:dT,7': 0.75,
    'rC:dT,4': 0.842105263,
    'rC:dT,5': 0.571428571,
    'rC:dT,2': 0.727272727,
    'rC:dT,3': 0.8666666670000001,
    'rC:dT,1': 1.0,
    'rA:dC,8': 0.7333333329999999,
    'rU:dT,1': 1.0,
    'rU:dC,3': 0.5,
    'rU:dC,1': 0.956521739,
    'rU:dT,2': 0.846153846,
    'rU:dG,19': 0.275862069,
    'rG:dT,14': 0.75,
    'rG:dT,15': 0.9411764709999999,
    'rG:dT,16': 1.0,
    'rG:dT,17': 0.933333333,
    'rG:dT,10': 0.933333333,
    'rG:dT,11': 1.0,
    'rA:dG,9': 0.571428571,
    'rA:dG,8': 0.428571429,
    'rA:dG,7': 0.4375,
    'rA:dG,6': 0.454545455,
    'rA:dG,5': 0.5,
    'rA:dG,4': 0.352941176,
    'rA:dG,3': 0.428571429,
    'rA:dG,2': 0.7857142859999999,
    'rA:dG,1': 0.857142857,
    'rU:dT,5': 0.5,
    'rG:dT,2': 0.846153846,
    'rA:dC,3': 0.611111111,
    'rA:dC,20': 0.764705882,
    'rG:dT,1': 0.9,
    'rG:dT,6': 1.0,
    'rG:dT,7': 1.0,
    'rG:dT,4': 0.9,
    'rC:dT,19': 0.428571429,
    'rG:dG,9': 0.538461538,
    'rG:dG,8': 0.615384615,
    'rG:dT,8': 1.0,
    'rG:dT,9': 0.642857143,
    'rU:dG,18': 0.428571429,
    'rU:dT,7': 0.875,
    'rG:dG,6': 0.681818182,
    'rA:dA,20': 0.6,
    'rU:dC,9': 0.6190476189999999,
    'rA:dG,17': 0.176470588,
    'rU:dC,8': 0.7333333329999999,
    'rA:dG,16': 0.0,
    'rA:dG,19': 0.20689655199999998,
    'rG:dT,3': 0.75,
    'rU:dG,3': 0.428571429,
    'rU:dG,2': 0.857142857,
    'rU:dG,1': 0.857142857,
    'rA:dG,18': 0.19047619,
    'rU:dG,7': 0.6875,
    'rU:dG,6': 0.9090909090000001,
    'rU:dG,5': 1.0,
    'rU:dG,4': 0.647058824,
    'rU:dG,9': 0.923076923,
    'rU:dG,8': 1.0,
    'rU:dC,19': 0.25,
    'rU:dC,18': 0.333333333,
    'rU:dC,13': 0.260869565,
    'rU:dC,12': 0.5,
    'rU:dC,11': 0.4,
    'rU:dC,10': 0.5,
    'rU:dC,17': 0.117647059,
    'rU:dC,16': 0.346153846,
    'rU:dC,15': 0.05,
    'rU:dC,14': 0.0,
    'rC:dT,10': 0.8666666670000001,
    'rC:dT,11': 0.75,
    'rC:dT,12': 0.7142857140000001,
    'rC:dT,13': 0.384615385,
    'rC:dT,14': 0.35,
    'rC:dT,15': 0.222222222,
    'rC:dT,16': 1.0,
    'rC:dT,17': 0.46666666700000003,
    'rC:dT,18': 0.538461538,
    'rA:dC,2': 0.8,
    'rA:dC,1': 1.0,
    'rA:dC,7': 0.705882353,
    'rA:dC,6': 0.7142857140000001,
    'rA:dC,5': 0.72,
    'rA:dC,4': 0.625,
    'rU:dT,19': 0.28571428600000004
}


@lru_cache(maxsize=1024 * 1024 * 1024)
def _revcom(s: str) -> str:
    """
    Reverse complements a given string
    """
    basecomp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U': 'A', 'N': 'N'}
    letters = list(s[::-1])
    letters = [basecomp[base] for base in letters]
    return ''.join(letters)


@lru_cache(maxsize=1024 * 1024 * 1024)
def calc_cfd(wt: str, sg: str, pam: str) -> float:
    """
    Calculates CFD score for NGG Cas9 guides.
    wt: WT 23mer sgRNA sequence
    sg: Off-target 20mer sgRNA sequence
    pam: 2bp PAM
    """
    # Disabled for perf. # TODO (gdingle): re-enable
    # assert len(pam) == 2, 'last two chars of PAM only'
    assert len(sg) == 20, len(sg)
    assert len(wt) == 20, len(wt)
    score = 1.0
    score *= pam_scores[pam]

    if score == 0:
        # early exit
        return score

    # Perf optimization
    mm = mm_scores
    for i, sl in enumerate(sg):
        if wt[i] != sl:
            score *= mm[_key(wt[i], sl, i)]

    return score


@lru_cache(maxsize=1024 * 1024 * 1024)
def _key(r, sl, i) -> str:
    return 'r' + r.replace('T', 'U') + ':d' \
        + _revcom(sl.replace('T', 'U')) + ',' + str(i + 1)


@lru_cache(maxsize=1024 * 1024 * 1024)
def cfd_score(
    wt: str, 
    sg: str, 
    pam='NGG', 
    guide_strand_same=True) -> float:
    """
    Alternate calling of calc_cfd with common inputs.
    Extremes.
    >>> cfd_score('AAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAA')
    1.0
    >>> cfd_score('AAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAG')
    0.764705882
    >>> cfd_score('GGGGGGGGGGGGGGGGGGGG', 'AAAAAAAAAAAAAAAAAAAA')
    0.0703372084129132
    >>> cfd_score('TTTTTTTTTTTTTTTTTTTT', 'AAAAAAAAAAAAAAAAAAAA')
    0.000537370080345105
    >>> cfd_score('CCCCCCCCCCCCCCCCCCCC', 'AAAAAAAAAAAAAAAAAAAA')
    8.801223840655041e-05
    Not important switch.
    'rG:dA,1': 1.0,
    >>> cfd_score('AAAAAAAAAAAAAAAAAAAA', 'GAAAAAAAAAAAAAAAAAAA')
    1.0
    Super important switch.
    'rA:dG,16': 1.0,
    >>> cfd_score('AAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAACAAAA')
    0.0
    Realistic.
    >>> cfd_score('AAGGCCAACCGGCGCCGCGC', 'GCGCGGCGCCGGTTGGCCTT')
    2.705350026063259e-06
    >>> cfd_score('GAAGGCCAACCGGCGCCGCG', 'CGCGGCGCCGGTTGGCCTTC')
    0.0
    Same examples as in https://github.com/maximilianh/crisporWebsite/blob/master/crispor.py#L1902.
    Based on source code provided by John Doench
    >>> cfd_score("GGGGGGGGGGGGGGGGGGGG", "GGGGGGGGGGGGGGGGGAAA")
    0.4635989007074176
    >>> cfd_score("GGGGGGGGGGGGGGGGGGGG", "GGGGGGGGGGGGGGGGGGGG")
    1.0
    >>> cfd_score("GGGGGGGGGGGGGGGGGGGG", "aaaaGaGaGGGGGGGGGGGG")
    0.5140384614450001
    >>> cfd_score("ATGGTCGGACTCCCTGCCAG", "ATGGTGGGACTCCCTGCCAG")
    0.5
    >>> cfd_score("ATGGTCGGACTCCCTGCCAG", "ATGATCCAAATCCCTGCCAG")
    0.53625000020625
    Change PAM.
    >>> cfd_score("AAAAAAAAAAAAAAAAAAAA", "AAAAAAAAAAAAAAAAAAAA", 'NGG')
    1.0
    >>> cfd_score("AAAAAAAAAAAAAAAAAAAA", "AAAAAAAAAAAAAAAAAAAA", 'NCC')
    0.0
    >>> cfd_score("AAAAAAAAAAAAAAAAAAAA", "AAAAAAAAAAAAAAAAAAAA", 'NGC')
    0.022222222000000003
    Reverse complements.
    >>> cfd_score("AAAAAAAAAAAAAAAAAAAA", "AAAAAAAAAAAAAAAAAAAA", 'CCN', guide_strand_same=False)
    1.0
    23bp calls.
    >>> cfd_score("AAAAAAAAAAAAAAAAAAAA", "CCNAAAAAAAAAAAAAAAAAAAA", guide_strand_same=False)
    1.0
    >>> cfd_score("AAAAAAAAAAAAAAAAAAAA", "TTNAAAAAAAAAAAAAAAAAAAA", guide_strand_same=False)
    0.0
    Non-symmetric.
    >>> cfd_score('CCAAAAAGGAGAGATGGTGCTGG', 'CCAAAAAGGAGAAATGGTCCTGG')
    0.41379310335013264
    >>> cfd_score('CCAAAAAGGAGAAATGGTCCTGG', 'CCAAAAAGGAGAGATGGTGCTGG')
    0.08152173912499999
    WT PAM bug
    >>> cfd_score('CCAGGTAGTGCCGCGCTGCCTGC', 'CCAGGggtggcggattggaagtt'.upper(), guide_strand_same=False)
    0.0001810095144455036
    >>> cfd_score('CCAGGTAGTGCCGCGCTGCCTGC', 'tgatgTAGTGCCGCGCTGCCTGC'.upper(), guide_strand_same=False)
    0.0
    >>> cfd_score('GCAGGCAGCGCGGCACTACC', 'TTTGAAGGTTAGGCGGTGGG')
    5.716801869405846e-05
    >>> cfd_score('CCAGCCAGTCCCACTCCAGCTCC', 'CCAGGTAGTGCCGCGCTGCCTGC', guide_strand_same=False)
    0.01913357577630766
    >>> cfd_score('CCAGGTAGTGCCGCGCTGCCTGC', 'CCAGCCAGTCCCACTCCAGCTCC', guide_strand_same=False)
    0.0
    problem case ENST00000368809
    >>> cfd_score('GTCGTCCACATGAAGCAGAAGGG', 'GTAGTACACATGAAGCAGAA')
    0.8047619054428573
    >>> cfd_score('CATCCTCCTGGACTCAATCA', 'CATCCTCCTGGACTCAATCANGT')
    0.016129031999999998
    >>> cfd_score('CATCCTCCTGGACTCAATCANGT', 'CATCCTCCTGGACTCAATCA')
    Traceback (most recent call last):
    ...
    ValueError: wild type should end in NGG
    >>> cfd_score('CCAATATGAAAAGGCCTAGTAAG', 'CCGATATGATGGGTGGCGGATTG', guide_strand_same=False)
    0.032817109477196474
    """
    if guide_strand_same is False:
        wt = _revcom(wt)
        sg = _revcom(sg)
        pam = _revcom(pam)

    # assumed to end in NGG
    if len(wt) >= 23 and not wt[20:23].endswith('GG'):
        raise ValueError('wild type should end in NGG')

    wt = wt[:20]

    if len(sg) == 23:
        sg, pam = sg[:20], sg[20:]
    else:
        sg = sg[:20]

    return calc_cfd(wt.upper(), sg.upper(), pam[-2:].upper())

@lru_cache(maxsize=1024 * 1024 * 1024)
def mit_hit_score(
    seq1: str,
    seq2: str,
    guide_strand_same=True,
    include_pam=False) -> float:

    """Compute MIT mismatch score between two 20-mers or 23-mers.
    See 'Scores of single hits' on http://crispr.mit.edu/about
    See calcHitScore in
    https://github.com/maximilianh/crisporWebsite/blob/master/crispor.py
    Parameters
    ----------
    seq1, seq2 : sequence
        two 20-mers to compare
    guide_strand_same : optional direction for starting with PAM
    include_pam : optional include extra 3bp for PAM.
    Returns
    -------
    float
        MIT mismatch score between the two sequences
    Extremes.
    >>> mit_hit_score('AAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAA')
    100.0
    >>> mit_hit_score('AAAAAAAAAAAAAAAAAAAA', 'GAAAAAAAAAAAAAAAAAAA')
    100.0
    >>> mit_hit_score('AAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAG')
    41.7
    >>> mit_hit_score('ZZZZZZZZZZZZZZZZZZZZ', 'AAAAAAAAAAAAAAAAAAAA')
    8.609700038185587e-08
    Realistic.
    >>> mit_hit_score('AAGGCCAACCGGCGCCGCGC', 'GCGCGGCGCCGGTTGGCCTT')
    6.039504885480631e-06
    >>> mit_hit_score('GAAGGCCAACCGGCGCCGCG', 'CGCGGCGCCGGTTGGCCTTC')
    1.6703747039472636e-05
    Other direction.
    >>> mit_hit_score('AAAAAAAAAAAAAAAAAAAA', 'GAAAAAAAAAAAAAAAAAAA', False)
    41.7
    >>> mit_hit_score('AAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAG', False)
    100.0
    Real case.
    >>> seq1 = 'CTAAGAGCATTTACACAATACA'[::-1]
    >>> seq2 = 'ctgAGAGCATTTACACAATACA'[::-1]
    >>> mit_hit_score(seq1, seq2)
    0.05972723076923077
    Include PAM.
    >>> mit_hit_score('AAAAAAAAAAAAAAAAAAAAGGG', 'AAAAAAAAAAAAAAAAAAAATGG', include_pam=True)
    100.0
    >>> mit_hit_score('AAAAAAAAAAAAAAAAAAAAAGG', 'AAAAAAAAAAAAAAAAAAAAATT', include_pam=True)
    0.20754716981132063
    """
    # aka Matrix "M"
    hit_score_m = [0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079, 0.445, 0.508,
                   0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583]
    if include_pam:
        # Add some high values, determined intuitively.
        hit_score_m += [0, 0.8, 0.8]

    # Go towards PAM
    if guide_strand_same is False:
        seq1 = seq1[::-1]
        seq2 = seq2[::-1]

    if include_pam:
        assert len(seq1) == 23
        max_dist = 22
    else:
        # Use most important 20bp only
        seq1 = seq1[-20:]
        seq2 = seq2[-20:]
        assert(len(seq1) == 20)
        max_dist = 19

    dists = []  # distances between mismatches, for part 2
    mm_count = 0  # number of mismatches, for part 3
    last_mm_pos = None  # position of last mismatch, used to calculate distance

    score1 = 1.0
    for pos in range(0, len(seq1)):
        if seq1[pos] != seq2[pos]:
            mm_count += 1
            if last_mm_pos != None:
                dists.append(pos - last_mm_pos)  # type: ignore
            score1 *= 1 - hit_score_m[pos]
            last_mm_pos = pos
    # 2nd part of the score
    if mm_count < 2:  # special case, not shown in the paper
        score2 = 1.0
    else:
        avg_dist = sum(dists) / len(dists)
        score2 = 1.0 / (((max_dist - avg_dist) / float(max_dist)) * 4 + 1)
    # 3rd part of the score
    if mm_count == 0:  # special case, not shown in the paper
        score3 = 1.0
    else:
        score3 = 1.0 / (mm_count**2)

    return score1 * score2 * score3 * 100

@lru_cache(maxsize=1024 * 1024 * 1024)
def manu_score(
    specificity_score: float, 
    hdr_dist: int, 
    codon='stop_codon') -> float:
    """
    Composite score to optimize guide selection by Manuel Leonetti.
    See https://goo.gl/KsUYJa.
    specificity_score as reported by Crispor.
    hdr_dist is cut-to-insert distance.
    codon was added to avoid regulatory region before or inside of start codon.
    >>> manu_score(100, 0)
    1.0
    >>> manu_score(0, 0)
    0.0
    >>> manu_score(0, 100)
    0.0
    >>> manu_score(60, 2)
    0.7232171842235433
    >>> manu_score(100, 0, 'start_codon')
    1.0
    >>> manu_score(100, -1, 'start_codon')
    0.19819005765761588
    >>> manu_score(60, -2, 'start_codon')
    0.14464343684470868
    """
    assert codon in ('start_codon', 'stop_codon'), codon
    score = _specificity_weight(specificity_score) * _dist_weight(hdr_dist)
    if codon == 'start_codon' and hdr_dist < -1:
        score *= _before_start_codon_penalty
    assert score >= 0 and score <= 1
    return score

@lru_cache(maxsize=1024 * 1024 * 1024)
def _specificity_weight(specificity_score: float):
    """
    >>> _specificity_weight(20)
    0
    >>> _specificity_weight(60)
    0.75
    >>> _specificity_weight(80)
    1
    """
    low = _specificity_weight_low
    high = _specificity_weight_high

    assert specificity_score >= 0 and specificity_score <= 100
    if specificity_score <= low:
        return 0
    elif specificity_score >= high:
        return 1
    else:
        return 1 / (high - low) * (specificity_score - low)

@lru_cache(maxsize=1024 * 1024 * 1024)
def _dist_weight(hdr_dist: int) -> float:
    """
    >>> _dist_weight(0)
    1.0
    >>> _dist_weight(5)
    0.7967034698934616
    >>> _dist_weight(10)
    0.402890321529133
    >>> _dist_weight(-20)
    0.026347980814448734
    """
    variance = _dist_weight_variance

    hdr_dist = abs(hdr_dist)  # make symmetric
    assert hdr_dist >= 0 and hdr_dist <= 100  # 100 is resonable upper bound

    # Returns a gaussian
    weight = math.exp((-1 * hdr_dist**2) / (2 * variance))
    assert weight >= 0 and weight <= 1
    return weight

@lru_cache(maxsize=1024 * 1024 * 1024)
def gcContent(seq):
    " return GC content as a float "
    c = 0
    for x in seq:
        if x in ["G","C"]:
            c+= 1
    return (float(c)/len(seq))

'''
if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.FAIL_FAST | doctest.ELLIPSIS)
'''

if __name__ == '__main__':
    import doctest
    doctest.testmod()

    import sys
    if len(sys.argv) == 3:
        print(cfd_score(sys.argv[1], sys.argv[2]))