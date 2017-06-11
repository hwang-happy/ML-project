# import Bio.PDB
# from Bio.PDB import PDBList
# from Bio.PDB.PDBParser import PDBParser
# from Bio.PDB import DSSP
# from Bio import SeqIO

from Bio.PDB import *

def dssp_ss(fileName):
    """
    This fuction retrieves secondary structure from DSSP database
    Parametry:
    :fileName: - PDB ID (String)
    :return: secondary structure (String)
    """
    path = "output102361.out"
    with open(path) as f:
        found = False
        for line in f:
            if fileName in line:
                found = True
                return line.split(' ')[2]
        if found == False:
            raise ValueError('Could not find id %s' % fileName)

