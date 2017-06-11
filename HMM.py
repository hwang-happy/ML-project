__author__ = 'Paulina Kania'
import pickle
from ghmm import *

##proteins is used to store all possible amino acids sequences: list of strings; it's taken from dssp files
proteins = []
##structures is used to store all possible structures sequences: list of strings; its taken from dssp files
structures = []

i = 0

alphabet_structure = ['C', 'H', 'E', 'T']
alphabet_proteins = ['A','R','N','D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

##model_path is the file containing the matrixes A, B and pi, which were prepaired in the file HMM_training.py
model_path = "hmm.dssp.model"

data_dict_loaded = pickle.load(open( model_path, "rb" ) )

A_hmm = data_dict_loaded["A"]
B_hmm = data_dict_loaded["B"]
pi_hmm = data_dict_loaded["pi"]

##sigma is our alphabet in hmm format
sigma = Alphabet(alphabet_proteins)
##m = HMM object
m = HMMFromMatrices(sigma,DiscreteDistribution(sigma),A_hmm,B_hmm,pi_hmm)

#index = 5600
#seqFasta = [ char for char in proteins[index]]
#seqFasta = [ char for char in "GTAGKVIKCKAAVLWEQKQPFSIEEIEVAPPKTKEVRIKILATGICRTDDHVIKGTMVSKFPVIVGHEATGIVESIGEGVTTVKPGDKVIPLFLPQCRECNACRNPDGNLCIRSDITGRGVLADGTTRFTCKGKPVHHFLNTSTFTEYTVVDESSVAKIDDAAPPEKVCLIGCGFSTGYGAAVKTGKVKPGSTCVVFGLGGVGLSVIMGCKSAGASRIIGIDLNKDKFEKAMAVGATECISPKDSTKPISEVLSEMTGNNVGYTFEVIGHLETMIDALASCHMNYGTSVVVGVPPSAKMLTYDPMLLFTGRTWKGCVFGGLKSRDDVPKLVTEFLAKKFDLDQLITHVLPFKKISEGFELLNSGQSIRTVLTF"]
#str_dssp = "CCTTCCEEEEEEEECCTTCCCEEEEEEECCCCTTEEEEEEEEEECCHHHHHHHHTCCCCCCCECCCCEEEEEEEEECTTCCCCCTTCEEEECCCCCCCCCHHHHCTTCCCCTTCCCCCCCCCTTCCCCEECCCCEEECCTTTCCCECEEEEEHHHEEECCTTCCHHHHHHHHTHHHHHHHHHHHHCCCCTTCEEEEECCCHHHHHHHHHHHHTTCCEEEEECCCHHHHHHHHHHTCCEEECHHHCCCCHHHHHHHHHTCCCCEEEECCCCHHHHHHHHTTCCTTTCEEEECCCCCTTCCEEECTHHHHTTCEEEECCHHHCCHHHHHHHHHHHHTTTCCCCHHHEEEEEEHHHHHHHHHHHHTTCCCEEEEEC"

def hmm_structure(method, seqFasta):
    """
    This function returns the mapped predicted structure according to the chosen method
    :param method: decodification method; it can be 'viterbi', 'forward' or 'backward'
    :return: predicted mapped structure sequence
    """
    seqFasta = [char for char in seqFasta]
    if method=='forward':
        forward = m.forward(EmissionSequence(sigma,seqFasta))
        ## ss_f contains the secondary structure predicted by the forward algorithm (mapped - values 0,1,2,3)
        ss_f = ""
        for lista in forward[0]:
        #we add to the ss_f string the index with the highest probability
            ss_f += str(lista.index((sorted(lista))[-1]))
        return ss_f
    elif method=='backward':
        backward = m.backward(EmissionSequence(sigma, seqFasta),m.forward(EmissionSequence(sigma,seqFasta))[-1])
        #m.forward(EmissionSequence(sigma,seqFasta))[-1] is the scaling vector from the forward algorithm
        ## ss_b contains the secondary structure predicted by the backward algorithm (mapped - values 0,1,2,3)
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
# print "HMM     ", mapping_num_to_letters(hmm_structure('viterbi'))
# print "Forward ", mapping_num_to_letters(hmm_structure('forward'))
# print "Backward", mapping_num_to_letters(hmm_structure('backward'))
# print "DSSP    ", str_dssp
