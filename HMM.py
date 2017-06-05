# -*- coding: utf-8 -*-
"""
@author: Paulina Kania
"""

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

print protein_mapping
print structure_mapping
