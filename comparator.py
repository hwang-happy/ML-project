from __future__ import division
import collections
import re

## alphabet_structure: List containing SSP secondary structure symbols
alphabet_structure = ['C', 'H', 'E', 'T']
## alphabet_proteins: List containing aminoacid IUPAC symbols
alphabet_proteins = ['A','R','N','D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

def compare_structures(list_of_ss):
    """
    This function return basic statistics of each seqence from
    from input list, icluding percentage of structure. After
    little modification it could return the most popular
    secondary structure in each protein.
    :param list_of_ss: list of strings containing seqences
    :return: list of structures` percentage.
    """

    list_of_dict = []
    common = []
    for ss in list_of_ss:
        list_of_dict.append(collections.Counter(ss))
        common.append(collections.Counter(ss).most_common(1))

    statistic = []
    for elem in list_of_dict:
        elem_stats = []
        sum_of_elem = sum(elem.values())
        for key in ('H', 'E', 'T', 'C'):
            if key not in elem.keys():
                elem_stats.append(float(0))
            else:
            #print dict(elem)[key]
                elem_stats.append(float(dict(elem)[key])/float(sum_of_elem))
        statistic.append(elem_stats)

    # print common
    out = []
    for seq in statistic:
        seq_score = []
        for score in seq:
            seq_score.append(round(score, 2))
        out.append(seq_score)
    return out

def HMM_consensus(list_of_seq):
    """
    This function compares 3 or more secondary structure
    seqences and returns consensus secondary structure
    :param list_of_seq:
    :return: consensus secondary structure (string)
    """
    consensus = ''
    alphabet = ['H', 'E', 'T', 'C']
    for i in range(len(list_of_seq[0])):
        count_ss = [0, 0, 0, 0]
        for ss in list_of_seq:
            if ss[i] == 'H':
                count_ss[0] += 1
            elif ss[i] == 'E':
                count_ss[1] += 1
            elif ss[i] == 'T':
                count_ss[2] += 1
            else:
                count_ss[3] +=1

        consensus += alphabet[count_ss.index(max(count_ss))]
    sec1 = re.sub('([^H])HH([^H])', r'\1\1\2\2', consensus)
    sec2 = re.sub('([^C])C([^C])', r'\1\1\2', sec1)
    sec3 = re.sub('([^H])HHH([^H])', r'\1HHHH', sec2)
    sec4 = re.sub('([^T])T([^T])', r'\1\1\2', sec3)
    sec5 = list(re.sub('([^H])H([^H])', r'\1\1\2', sec4))
    return ''.join(sec5)
#najpierw HMM z DSSP, potem CF z DSSP, uruchamiane gdy jest PDB potem HMM z CF

def identities(two_seq):
    """
    This function compares two seqences and checks how many identic
    elements are presented in each position of this seqences
    :param two_seq: list of two sequences (string)
    :return: percentage of identities (float)
    """
    identity = 0

    countBadSS = 0
    countBadAA = 0
    for literka in two_seq[0]:
        if literka not in alphabet_structure:
            countBadSS += 1
    for literka in two_seq[1]:
        if literka not in alphabet_proteins:
            countBadAA += 1
    if countBadSS != 0 or countBadAA != 0 or len(two_seq[0]) != len(two_seq[1]):
            return
    for i in range(len(two_seq[0])):
        if two_seq[0][i] == two_seq[1][i]:
            identity += 1

    out = round(identity/len(two_seq[0]), 2)
    return out

"""
Test sequences
"""
# list_of_ss = ['HHHHHHHHHHHHHHHHHHHHHHHHHEEEEEEEEEEEEEEEEHHHHCCCCCCCCCCCCTTTTTTTTTTTHHHHHHHHHHHHH',
#               'HHHHHHHHHHHHHEEEEEEEEEEEEEEEEEEEEEEEEEEHHHHHHHHHHHHHHHHHHHCCCCCCCCCCCCTTTTTTTTTTT']
#
# print compare_structures(list_of_ss=list_of_ss)
# print HMM_consensus(list_of_ss)
# print identities(list_of_ss)
