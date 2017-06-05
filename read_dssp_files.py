# -*- coding: utf-8 -*-
"""
@author: Paulina Kania
"""

import Bio
from Bio.PDB import PDBParser
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import os

files = 5597
dssp_path = "/home/paulina/GIT/github/ML/ML-project/dssp/"


def readingFiles():
    last_file_name = ""
    try:
        listOfFiles = os.listdir(dssp_path)
        output = open("/home/paulina/GIT/github/ML/ML-project/output" + str(files) + ".out", "w")
        for name in listOfFiles[:files]:
            last_file_name = name
            f = open(dssp_path + str(name))
            dssp_start = defining_start_dssp(f)
            struct = ""
            iter = 0
            seq = ""
            for letter in sequence(f, dssp_start):
                if letter != '!':
                    iter += 1
                    seq += letter
                else:
                    break
            for symbol in structure(f, dssp_start)[:iter]:
                if symbol == 'G' or symbol == 'I' :
                    struct += "H"
                elif symbol == 'B':
                    struct += 'E'
                elif symbol in ["-", "S"]:
                    struct += "C"
                else:
                    struct += symbol
            output.write(getID(f) + " " + seq + " " + struct + "\n")
            f.close()
        output.close()
    except TypeError:
        print last_file_name

def defining_start_dssp(file):
    '''
    :param file: input file
    :return: next line after "#"
    '''
    file.seek(0,0)
    i = 1
    for line in file:
        if "#" in line:
            return i
        else:
            i += 1
            continue

def structure(file, dssp_start):
    '''
    :param file: input file
    :param dssp_start: result of the function defining_start_dssp(file)
    :return: secondary structure predicted by DSSP
    '''
    file.seek(0,0)
    out = []
    for line in file:
        out.append(line[16])
    out = ''.join(out[dssp_start:])
    out = out.replace(" ", "-")
    return out

def sequence(file, dssp_start):
    '''
    :param file: input file
    :param dssp_start: result of the function defining_start_dssp(file)
    :return: aminoacid sequence read from PDB database, included in DSSP file
    '''
    file.seek(0,0)
    out = []
    i = 0
    for line in file:
        out.append(line[13])
    out = ''.join(out[dssp_start:])
    return out

def getID(file):
    """
    This function takes the PDB ID from DSSP file
    :param file: input file in DSSP format
    :return: PDB ID as string
    """
    file.seek(0,0)
    splited = []
    for line in file:
        splited = line.split()
        if "HEADER" in line:
            return splited[-2]


readingFiles()
