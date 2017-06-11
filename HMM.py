# -*- coding: utf-8 -*-
"""
@author: Paulina Kania
"""
from collections import Counter
from collections import defaultdict
import pickle
from ghmm import *

##proteins is used to store all possible amino acids sequences: list of strings; it's taken from dssp files
proteins = []
##structures is used to store all possible structures sequences: list of strings; its taken from dssp files
structures = []

#max let you can set the number of files you want to use for training the HMM
max = float('inf')
i = 0

##path is the file for the training of the HMM
path = "output5597.out"
with open(path) as f:
    for line in f:
        splited = line.strip().split(' ')
        proteins.append(splited[1].upper())
        structures.append(splited[2].upper())
        i += 1
        if i > max:
            break

alphabet_structure = ['C', 'H', 'E', 'T']
alphabet_proteins = ['A','R','N','D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

def create_mapping(alphabet):
    """
    Change list of chars to list of ints, taking sequencial natural numbers
    :param alphabet: list of char
    :return: dictionary with keys that are letters from alphabet and ints as values
    """
    mapping = {}
    for (letter, i) in zip(alphabet, range(len(alphabet))):
        mapping[letter] = i
    return mapping

protein_mapping = create_mapping(alphabet_proteins)
structure_mapping = create_mapping(alphabet_structure)

# print protein_mapping
# print structure_mapping

"""
Matrix B and vector pi creation
Matrix B contains the probabilities of the presence of each amino acid with
each type of structure
"""

#Counter is a dictionary with letters as keys
structure_counter = Counter()
pi = Counter()
proteins_counter = Counter()
B = defaultdict(lambda: Counter())

alphabet_structure_set = set(alphabet_structure)
alphabet_proteins_set = set(alphabet_proteins)

bad_letters_proteins_all = []
bad_letters_structures_all = []

for (structure, protein) in zip(structures, proteins):
    ##it counts which letter of structure appears the most
    bad_letters_struct = set(structure) - alphabet_structure_set
    if bad_letters_struct:
        bad_letters_structures_all.append(bad_letters_struct)
        continue

    bad_letters_proteins = set(protein) - alphabet_proteins_set
    if bad_letters_proteins:
        bad_letters_proteins_all.append(bad_letters_proteins)
        continue

    structure_counter += Counter(structure)
    pi[structure[0]] = pi[structure[0]] + 1
    for leter_structure in alphabet_structure:
        #search for positions in which each type of structure is present
        indices = [i for i, x in enumerate(structure) if x == leter_structure]
        #we search for amino acids which are in the indices found above
        secondary_by_indices_counter = Counter([protein[i] for i in indices])
        #sum of appearance for each letter
        B[leter_structure] = B[leter_structure] + secondary_by_indices_counter

print "Bad proteins letters:", bad_letters_proteins_all
print "Bad structures letters:", bad_letters_structures_all

##normalization of the pi vector - we want to have the probability instead of the exact number
# pi vector is the path = sequence of states
for k in pi:
    pi[k] /= float(len(structures))

# print pi

"""
Matrix A creation
Matrix A contains the probabilities of passing consequently
from one structure's type to another
"""
A = defaultdict(lambda: defaultdict(lambda: 0))
for structure in structures:
    bad_letters_struct = set(structure) - alphabet_structure_set
    if bad_letters_struct:
        continue
    for (structure_begin, structure_ends) in zip(structure[:-1], structure[1:]):
        A[structure_begin][structure_ends] += 1

""" matrix A, B normalization"""
for matrix in [A, B]:
    for structure_letter, counter_structure_letter in matrix.items():
        structure_counter_sum = sum(counter_structure_letter.values(), 0.0)
        for protein_letter in counter_structure_letter:
            counter_structure_letter[protein_letter] /= structure_counter_sum

def dict_to_list(my_dict, mapping):
    """
    having a dictionary, it returns a list
    :param my_dict: dictionary to be mapped
    :param mapping: alphabet for mapping
    :return: list of ints
    """
    my_list = [0] * len(mapping.keys())
    for key in my_dict:
        my_list[mapping[key]] = my_dict[key]
    return my_list

def dict_dict_to_list_list(my_dict, column_mapping, row_mapping):
    """
    this function takes a dictionary of dictionaries and returns a list of lists (matrix)
    :param my_dict: dictionary to being mapped
    :param column_mapping: dictionary containing the alphabet for the mapping
    :param row_mapping: dictionary containing the alphabet for the mapping
    :return: matrix with probabilities
    """
    matrix = [[0] * len(column_mapping.keys()) ] * len(row_mapping) #rows creation
    for key in my_dict:
        matrix[row_mapping[key]] = dict_to_list(my_dict[key], column_mapping)
    return matrix

pi_hmm = dict_to_list(pi, structure_mapping)
B_hmm = dict_dict_to_list_list(B, protein_mapping, structure_mapping)
A_hmm = dict_dict_to_list_list(A, structure_mapping, structure_mapping)

##sigma is our alphabet in hmm format
sigma = Alphabet(alphabet_proteins)
##m=HMM object
m = HMMFromMatrices(sigma,DiscreteDistribution(sigma),A_hmm,B_hmm,pi_hmm)

# index = 5500
# seqFasta = [ char for char in proteins[index]]
# seqFasta = [ char for char in "GTAGKVIKCKAAVLWEQKQPFSIEEIEVAPPKTKEVRIKILATGICRTDDHVIKGTMVSKFPVIVGHEATGIVESIGEGVTTVKPGDKVIPLFLPQCRECNACRNPDGNLCIRSDITGRGVLADGTTRFTCKGKPVHHFLNTSTFTEYTVVDESSVAKIDDAAPPEKVCLIGCGFSTGYGAAVKTGKVKPGSTCVVFGLGGVGLSVIMGCKSAGASRIIGIDLNKDKFEKAMAVGATECISPKDSTKPISEVLSEMTGNNVGYTFEVIGHLETMIDALASCHMNYGTSVVVGVPPSAKMLTYDPMLLFTGRTWKGCVFGGLKSRDDVPKLVTEFLAKKFDLDQLITHVLPFKKISEGFELLNSGQSIRTVLTF"]
# str_dssp = "CCTTCCEEEEEEEECCTTCCCEEEEEEECCCCTTEEEEEEEEEECCHHHHHHHHTCCCCCCCECCCCEEEEEEEEECTTCCCCCTTCEEEECCCCCCCCCHHHHCTTCCCCTTCCCCCCCCCTTCCCCEECCCCEEECCTTTCCCECEEEEEHHHEEECCTTCCHHHHHHHHTHHHHHHHHHHHHCCCCTTCEEEEECCCHHHHHHHHHHHHTTCCEEEEECCCHHHHHHHHHHTCCEEECHHHCCCCHHHHHHHHHTCCCCEEEECCCCHHHHHHHHTTCCTTTCEEEECCCCCTTCCEEECTHHHHTTCEEEECCHHHCCHHHHHHHHHHHHTTTCCCCHHHEEEEEEHHHHHHHHHHHHTTCCCEEEEEC"

def hmm_structure(method, seqFasta):
    """
    This function returns the mapped predicted structure according to the chosen method
    :param method: decodification method; it can be 'viterbi', 'forward' or 'backward'
    :return: predicted mapped structure sequence
    """
    seqFasta = [char for char in seqFasta]
    if method=='viterbi':
        ss_hmm = ''.join(map(str, m.viterbi(EmissionSequence(sigma, seqFasta))[0]))
        return ss_hmm
    elif method=='forward':
        forward = m.forward(EmissionSequence(sigma,seqFasta))
        #ss_f contains the secondary structure predicted by the forward algorithm (mapped - values 0,1,2,3)
        ss_f = ""
        for lista in forward[0]:
        #we add to the ss_f string the index with the highest probability
            ss_f += str(lista.index((sorted(lista))[-1]))
        return ss_f
    elif method=='backward':
        backward = m.backward(EmissionSequence(sigma, seqFasta),m.forward(EmissionSequence(sigma,seqFasta))[-1])
        #m.forward(EmissionSequence(sigma,seqFasta))[-1] is the scaling vector from the forward algorithm
        #ss_b contains the secondary structure predicted by the backward algorithm (mapped - values 0,1,2,3)
        ss_b = ""
        for lista in backward:
            ss_b += str(lista.index((sorted(lista))[-1]))
        return ss_b

structure_mapping_numbers = {0: 'C', 1: 'H', 2: 'E', 3: 'T'}
def mapping_num_to_letters(sss):
    """
    This function takes as input a string containing as key the number of a structure type
    and returns the corresponding letter in SS
    :param sss: predicted mapped secondary structure - sequence of numbers
    :return: string containing the unmapped secondary structure
    H = helix
    E = sheet
    T = turn
    C = other
    """
    ss_out = ""
    for num in sss:
        for key in structure_mapping_numbers.keys():
            if int(num) == key:
                ss_out += structure_mapping_numbers[key]
    return ss_out

# print "A", A
# print "B", B
# print "pi", pi
# print "HMM     ", mapping_num_to_letters(hmm_structure('viterbi', seqFasta))
# print "Forward ", mapping_num_to_letters(hmm_structure('forward', seqFasta))
# print "Backward", mapping_num_to_letters(hmm_structure('backward', seqFasta))
# print "DSSP    ", str_dssp


data_dict = {
    "A": A_hmm,
    "B": B_hmm,
    "pi": pi_hmm
}

#writing our matrixes to a file
model_path = "hmm.dssp.model"
pickle.dump( data_dict, open( model_path, "wb" ) )

data_dict_loaded = pickle.load(open( model_path, "rb" ) )
# print data_dict_loaded["A"]
