# taken from https://github.com/czbiohub/crispycrunch/blob/master/utils/mitscore.py
from functools import lru_cache


@lru_cache(maxsize=1024 * 1024)
def mit_hit_score(
    seq1: str, seq2: str, guide_strand_same=True, include_pam=False
) -> float:
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
    hit_score_m = [
        0,
        0,
        0.014,
        0,
        0,
        0.395,
        0.317,
        0,
        0.389,
        0.079,
        0.445,
        0.508,
        0.613,
        0.851,
        0.732,
        0.828,
        0.615,
        0.804,
        0.685,
        0.583,
    ]
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
        assert len(seq1) == 20
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
        score3 = 1.0 / (mm_count ** 2)

    return score1 * score2 * score3 * 100


if __name__ == "__main__":
    import doctest

    doctest.testmod()
