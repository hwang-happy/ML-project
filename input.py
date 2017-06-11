__author__ = 'Paulina Kania'

import tkMessageBox

def check_input(fasta):
    """
    This function checks is the input is in the correct format (fasta protein sequence)
    :param fasta: input in fasta format
    :return: protein sequence as string
    """
    onlySEQ = ""
    pdbID = ""
    for sign in fasta:
        if sign[0] == '>':
            break
        else:
            tkMessageBox.showinfo("Info", "This is not a fasta format")
            return
    ## countNA is the counter for nucleotides; used for checking if the code provided is not RNA or DNA
    countNA = 0
    ## countAll is the counter for all the aminoacids in the sequence
    countAll = 0
    countEndLine = 0
    if fasta[-1] == '*':
        fasta[-1] = ""
    for sign in fasta:
        if sign != unichr(10) and countEndLine == 0:
            continue
        else:
            if sign in ('A','C','T','G','U'):
                countNA += 1
                countAll += 1
                onlySEQ += sign
            elif sign == unichr(10):
                countEndLine += 1
            else:
                countAll += 1
                onlySEQ += sign
    if countNA == countAll or countNA==countAll-1:
        tkMessageBox.showinfo("Info", "Please, insert a protein sequence")
    else:
        pdbID = fasta[1:5]
        return pdbID, onlySEQ

