# === SOURCE CODE cfd-score-calculator.py provided by John Doench =====
# The CFD score is an improved specificity score
import os
import pickle
import re
import sys
from os.path import dirname
import pandas as pd

def calcMitGuideScore(hitSum):
    """ Sguide defined on http://crispr.mit.edu/about
    Input is the sum of all off-target hit scores. Returns the specificity of the guide.
    """
    score = 100 / (100 + hitSum)
    score = int(round(score * 100))
    return score


def calcMitGuideScore_v2(hitSum):
    """ Sguide defined on http://crispr.mit.edu/about
    Input is the sum of all off-target hit scores. Returns the specificity of the guide.
    """
    score = 1 / (1 + hitSum)
    score = int(round(score * 100))
    return score


def get_mm_pam_scores():
    """
    """
    dataDir = os.path.join(dirname(__file__), "bin", "CFD_Scoring")
    mm_scores = pickle.load(open(os.path.join(dataDir, "mismatch_score.pkl"), "rb"))
    pam_scores = pickle.load(open(os.path.join(dataDir, "pam_scores.pkl"), "rb"))
    return (mm_scores, pam_scores)

def get_mm_pam_scores_cas12a():
    """
    """
    dataDir = os.path.join(dirname(__file__), "bin", "CFD_Scoring_Cas12a")
    df = pd.read_csv(os.path.join(dataDir, "off_targ_enCas12a.csv"))
    df['mm_pos'] = df["MM"] + "," + df["Pos"].astype(str)
    mm_scores = dict(zip(df['mm_pos'], df['avg_percent_active']))
    return mm_scores


# Reverse complements a given string
def revcom(s):
    basecomp = {"A": "T", "C": "G", "G": "C", "T": "A", "U": "A"}
    letters = list(s[::-1])
    letters = [basecomp[base] for base in letters]
    return "".join(letters)


# Calculates CFD score
def calc_cfd(wt, sg, pam):
    # mm_scores,pam_scores = get_mm_pam_scores()
    score = 1
    sg = sg.replace("T", "U")
    wt = wt.replace("T", "U")
    s_list = list(sg)
    wt_list = list(wt)
    for i, sl in enumerate(s_list):
        if wt_list[i] == sl:
            score *= 1
        else:
            key = "r" + wt_list[i] + ":d" + revcom(sl) + "," + str(i + 1)
            score *= mm_scores[key]
    score *= pam_scores[pam]
    return score

def calc_cfd_cas12a(wt, sg):
    mm_scores_cas12a = get_mm_pam_scores_cas12a()
    score = 1
    sg = sg.replace("T", "U")
    wt = wt.replace("T", "U")
    s_list = list(sg)
    wt_list = list(wt)
    for i, sl in enumerate(s_list):
        if wt_list[i] == sl:
            score *= 1
        else:
            key = "r" + wt_list[i] + ":d" + revcom(sl) + "," + str(i + 1)
            score *= mm_scores_cas12a[key]
    return score


mm_scores, pam_scores = None, None
mm_scores_cas12a = None


def calcCfdScore(guideSeq, otSeq):
    """ based on source code provided by John Doench
    >>> calcCfdScore("GGGGGGGGGGGGGGGGGGGGGGG", "GGGGGGGGGGGGGGGGGAAAGGG")
    0.4635989007074176
    >>> calcCfdScore("GGGGGGGGGGGGGGGGGGGGGGG", "GGGGGGGGGGGGGGGGGGGGGGG")
    1.0
    >>> calcCfdScore("GGGGGGGGGGGGGGGGGGGGGGG", "aaaaGaGaGGGGGGGGGGGGGGG")
    0.5140384614450001
    # mismatches:      *               !!
    >>> calcCfdScore("ATGGTCGGACTCCCTGCCAGAGG", "ATGGTGGGACTCCCTGCCAGAGG")
    0.5

    # mismatches:    *  ** *
    >>> calcCfdScore("ATGGTCGGACTCCCTGCCAGAGG", "ATGATCCAAATCCCTGCCAGAGG")
    0.53625000020625

    >>> calcCfdScore("ATGTGGAGATTGCCACCTACCGG", "ATCTGGAGATTGCCACCTACAGG")

    """
    global mm_scores, pam_scores
    if mm_scores is None:
        mm_scores, pam_scores = get_mm_pam_scores()
    wt = guideSeq.upper()
    off = otSeq.upper()
    m_wt = re.search("[^ATCG]", wt)
    m_off = re.search("[^ATCG]", off)
    if (m_wt is None) and (m_off is None):
        pam = off[-2:]
        sg = off[:20]
        cfd_score = calc_cfd(wt, sg, pam)
        return cfd_score


def calcCfdScore_cas12a(guideSeq, otSeq):
    """ PAM sequence must not be included
    >>> calcCfdScore("GGGGGGGGGGGGGGGGGGGGGGG", "GGGGGGGGGGGGGGGGGAAAGGG")
    0.4635989007074176
    >>> calcCfdScore("GGGGGGGGGGGGGGGGGGGGGGG", "GGGGGGGGGGGGGGGGGGGGGGG")
    1.0
    >>> calcCfdScore("GGGGGGGGGGGGGGGGGGGGGGG", "aaaaGaGaGGGGGGGGGGGGGGG")
    0.5140384614450001
    # mismatches:      *               !!
    >>> calcCfdScore("ATGGTCGGACTCCCTGCCAGAGG", "ATGGTGGGACTCCCTGCCAGAGG")
    0.5

    # mismatches:    *  ** *
    >>> calcCfdScore("ATGGTCGGACTCCCTGCCAGAGG", "ATGATCCAAATCCCTGCCAGAGG")
    0.53625000020625

    >>> calcCfdScore("ATGTGGAGATTGCCACCTACCGG", "ATCTGGAGATTGCCACCTACAGG")

    """
    global mm_scores_cas12a
    if mm_scores_cas12a is None:
        mm_scores_cas12a = get_mm_pam_scores_cas12a()
    wt = guideSeq.upper()
    off = otSeq.upper()
    m_wt = re.search("[^ATCG]", wt)
    m_off = re.search("[^ATCG]", off)
    if (m_wt is None) and (m_off is None):
        sg = off[:20]
        wt = wt[:20]
        cfd_score = calc_cfd_cas12a(wt, sg)
        return cfd_score

# ==== END CFD score source provided by John Doench

# from https://github.com/maximilianh/crisporWebsite
hitScoreM = [
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


def calcHitScore(string1, string2):
    " see 'Scores of single hits' on http://crispr.mit.edu/about "
    # The Patrick Hsu weighting scheme
    # S. aureus requires 21bp long guides. We fudge by using only last 20bp
    matrixStart = 0
    maxDist = 19

    assert string1[0].isupper()
    assert len(string1) == len(string2)
    # for nmCas9 and a few others with longer guides, we limit ourselves to 20bp
    if len(string1) > 20:
        string1 = string1[-20:]
        string2 = string2[-20:]
    # for 19bp guides, we fudge a little, but first pos has no weight anyways
    elif len(string1) == 19:
        string1 = "A" + string1
        string2 = "A" + string2
    # for shorter guides, I'm not sure if this score makes sense anymore, we force things
    elif len(string1) < 19:
        matrixStart = 20 - len(string1)
        maxDist = len(string1) - 1

    assert len(string1) == len(string2)

    dists = []  # distances between mismatches, for part 2
    mmCount = 0  # number of mismatches, for part 3
    lastMmPos = None  # position of last mismatch, used to calculate distance

    score1 = 1.0
    for pos in range(matrixStart, len(string1)):
        if string1[pos] != string2[pos]:
            mmCount += 1
            if lastMmPos != None:
                dists.append(pos - lastMmPos)
            score1 *= 1 - hitScoreM[pos]
            lastMmPos = pos
    # 2nd part of the score
    if mmCount < 2:  # special case, not shown in the paper
        score2 = 1.0
    else:
        avgDist = sum(dists) / len(dists)
        score2 = 1.0 / (((maxDist - avgDist) / float(maxDist)) * 4 + 1)
    # 3rd part of the score
    if mmCount == 0:  # special case, not shown in the paper
        score3 = 1.0
    else:
        score3 = 1.0 / (mmCount ** 2)
    score = score1 * score2 * score3 * 100
    return score


if __name__ == "__main__":
    print(calcCfdScore("GGGGGGGGGGGGGGGGGGGGGGG", "GGGGGGGGGGGGGGGGGGGCGGG"))
    print(calcHitScore("ATGTGGAGATTGCCACCTAC", "ATCTGGAGATTGCCACCTAC"))
    print(calcCfdScore_cas12a("GGGGGGGGGGGGGGGGGGGGGGG", "GGGGGGGGGGGGGGGGGGGCGGG"))
    print(calcHitScore("ATGTGGAGATTGCCACCTAC", "ATCTGGAGATTGCCACCTAC"))
