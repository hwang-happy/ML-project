from __future__ import division
import re

##Algorithm parameters.
## TURN_PRODUCT is cutoff value of Turns
TURN_PRODUCT = 0.75E-4
## ALPHA_CUTis cutoff value of Helix
ALPHA_CUT = 100
## BETA_CUT is cutoff value of Sheet
BETA_CUT = 100

##dataDict contains a Chou Fasman probabilities of amino acid occurences in secondary structures.
# For more infomration in proteinData description or in publication

dataDictFromPub = {'A' : [142,139, 83, 79, 66, 0.060,0.076,0.035,0.058,    'H','H','i','i',0,0,0],
            'R' : [98,100, 93, 94, 95,  0.070,0.106,0.099,0.085,    'i','h','i','i',0,0,0],
            'N' : [67, 78, 89, 66,156,  0.161,0.083,0.191,0.091,    'b','i','i','b',0,0,0],
            'D' : [101,106, 54, 66,146, 0.147,0.110,0.179,0.081,    'I','h','B','b',0,0,0],
            'C' : [70, 95,119,107,119,  0.149,0.053,0.117,0.128,    'i','i','h','h',0,0,0],
            'Q' : [111,112,110,100, 98, 0.074,0.098,0.037,0.098,    'h','h','h','I',0,0,0],
            'E' : [151,144, 37, 51, 74, 0.056,0.060,0.077,0.064,    'H','H','B','b',0,0,0],
            'G' : [57, 64, 75, 87,156,  0.102,0.085,0.190,0.152,    'B','B','b','i',0,0,0],
            'H' : [100,112, 87, 83, 95, 0.140,0.047,0.093,0.054,    'I','h','i','i',0,0,0],
            'I' : [108, 99,160,157, 47, 0.043,0.034,0.013,0.056,    'h','i','H','H',0,0,0],
            'L' : [121,130,130,117, 59, 0.061,0.025,0.036,0.070,    'H','H','H','h',0,0,0],
            'K' : [116,121, 74, 73,101, 0.055,0.115,0.072,0.095,    'h','h','b','b',0,0,0],
            'M' : [145,132,105,101, 60, 0.088,0.082,0.014,0.055,    'H','H','h','I',0,0,0],
            'F' : [113,111,138,123, 60, 0.059,0.041,0.065,0.065,    'h','h','h','h',0,0,0],
            'P' : [57, 55, 55, 62,152,  0.102,0.301,0.034,0.068,    'B','B','B','B',0,0,0],
            'S' : [77, 72, 75, 94,143,  0.120,0.139,0.125,0.106,    'i','b','b','i',0,0,0],
            'T' : [83, 78,119,133, 96,  0.086,0.108,0.065,0.079,    'i','i','h','h',0,0,0],
            'W' : [108,103,137,124, 96, 0.017,0.013,0.064,0.167,    'h','I','h','h',0,0,0],
            'Y' : [69, 73,147,131,114,  0.082,0.065,0.114,0.125,    'b','b','H','h',0,0,0],
            'V' : [106, 97,170,164, 50, 0.062,0.048,0.028,0.053,    'h','i','H','H',0,0,0]
            }

class proteinData():
    """
    A Chou Fasman data structure connected to dataDict.
    """
    ## code: (char) IUPAC amino acid code
    code = 0
    ## p_alpha: (float) Probability of helix for amino acid (dataDict[code][0])
    p_alpha = 0
    ## p_beta:(float) Probability of sheet for amino acid (dataDict[code][2])
    p_beta = 0
    ## p_turn: (float) Probability of turn for amino acid (dataDict[code][4])
    p_turn = 0
    ## p_turn: (float) Statistical list of params (dataDict[code][5:8])
    bend = []
    ## alpha_class: (char) Other param from Chou-Fasman publication
    alpha_class = 0
    ## beta_class: (char) Other param from Chou-Fasman publication
    beta_class = 0
    ## p4_alpha: (float) Probability of helix for four aminoacids
    p4_alpha = 0
    ## p4_beta: (float) Probability of sheet for four aminoacids
    p4_beta = 0
    ## p4_turn: (float) Probability of turn for four aminoacids
    p4_turn = 0
    ## turn_prod: (float) turn cutoff value
    turn_prod = 0
    ## helix_nucl: (boolean) probably helix (* in publication)
    helix_nucl = False
    ## sheet_nucl: (boolean) probably sheet (* in publication)
    sheet_nucl = False
    ## turn_prop: (boolean) probably turn (* in publication)
    turn_prop = False

    def __init__(self, char):
        self.code = char

    def __getitem__(self, key):
        return key

class ChouFasman():
    def __init__(self):
        self.dataDict = dataDictFromPub
        ##self.aaData contains Chou Fasman structures of each aminoacid in protein
        self.aaData = []

    def main(self, seqq):
        """
        The main function creates a list of Chou Fasman data structures.
        :param seqq: input sequence
        """
        for aa in seqq:
            aU = aa.upper()
            self.aaData.append(proteinData(aU))

    def print_it(self):
        """
        The print_it function prints the Chou Fasman data structure.
        """
        for i in range(0, len(self.aaData)):
            print i, self.aaData[i].code, self.aaData[i].p_alpha, self.aaData[i].p_beta, self.aaData[i].p_turn,\
                  self.aaData[i].alpha_class, self.aaData[i].beta_class,\
                  self.aaData[i].p4_alpha, self.aaData[i].p4_beta, self.aaData[i].p4_turn, self.aaData[i].turn_prod,\
                self.aaData[i].helix_nucl, self.aaData[i].sheet_nucl, self.aaData[i].turn_prop

    def get_probability(self):
        """
        The get_probability function counts probabilities of secondary structure at amino acid position.
        """
        for aa in self.aaData:
            aa.p_alpha = self.dataDict[aa.code][0]
            aa.p_beta = self.dataDict[aa.code][2]
            aa.p_turn = self.dataDict[aa.code][4]
            aa.bend = [self.dataDict[aa.code][5], self.dataDict[aa.code][6], self.dataDict[aa.code][7],
                       self.dataDict[aa.code][8]]
            aa.alpha_class = self.dataDict[aa.code][9]
            aa.beta_class = self.dataDict[aa.code][10]

    def tetra_ave(self, seqq):
        """
        The tetra_ave function counts probabilities of secondary structure at four aminoacids positions.
        :param seqq: input sequence
        """
        for i in range (1, len(seqq)- 3):
            asum = bsum = tsum = 0
            tprod = 1.0
            for j in range(4):
                asum += self.aaData[i+j].p_alpha
                bsum += self.aaData[i+j].p_beta
                tsum += self.aaData[i+j].p_turn
                tprod *= float(self.aaData[i+j].bend[j])
            self.aaData[i].p4_alpha = asum/4
            self.aaData[i].p4_beta = bsum/4
            self.aaData[i].p4_turn = tsum/4
            self.aaData[i].turn_prod = tprod


    def propagation(self):
        """
        The propagation function compares probabilities of secondary structure variants and complements
        Chou Fasman data structure by secondary structure candidates.
        """
        #Helix
        helix_iterator = 0
        for i in range(len(self.aaData)-3):
            if self.aaData[i].p4_alpha >= 100 and self.aaData[i+1].p4_alpha >= 100 and self.aaData[i+2].p4_alpha >= 100 \
                    and self.aaData[i+3].p4_alpha >= 100:
                self.aaData[i].helix_nucl = True
                self.aaData[i+1].helix_nucl = True
                self.aaData[i+2].helix_nucl = True
                self.aaData[i+3].helix_nucl = True
                i += 4
                helix_iterator += 4
                if self.aaData[i].p4_alpha >= 100 and self.aaData[i].code != 'P' \
                        and (self.aaData[i].alpha_class == 'H' or  self.aaData[i].alpha_class == 'h'):
                    self.aaData[i].helix_nucl = True
                    helix_iterator += 1
                elif self.aaData[i].p4_alpha >= 103 and self.aaData[i].code != 'P' and self.aaData[i-1].helix_nucl == True \
                        and self.aaData[i].alpha_class == 'i' and (helix_iterator + 1) % 4 == 0:
                    self.aaData[i].helix_nucl = True
                    helix_iterator += 1

                elif self.aaData[i].code == 'P':
                    helix_iterator == 0
                    i += 1
                else:
                    helix_iterator == 0
        beta_iterator = 0

        #Sheet
        for i in range(len(self.aaData)-3):
            if self.aaData[i].p4_beta >= 100 and self.aaData[i+1].p4_beta >= 100 and self.aaData[i+2].p4_alpha >= 100:
                self.aaData[i].sheet_nucl = True
                self.aaData[i+1].sheet_nucl = True
                self.aaData[i+2].sheet_nucl = True
                i += 3
                beta_iterator += 3
                if self.aaData[i].p4_beta >= 100 and (self.aaData[i].beta_class == 'H'
                                                      or  self.aaData[i].beta_class == 'h'):
                    self.aaData[i].sheet_nucl = True
                    beta_iterator += 1
                    i += 1
                elif self.aaData[i].p4_beta >= 100 and self.aaData[i].beta_class == 'i':
                    p4_beta_SUM = 0
                    for j in range(i - beta_iterator, i):
                        p4_beta_SUM += self.aaData[j].p4_beta

                    if p4_beta_SUM / len(range(i - beta_iterator, i)) >= 103:
                        self.aaData[i].sheet_nucl = True
                        beta_iterator += 1
                        i += 1
                else:
                    beta_iterator == 0
        #Turns
        for i in range(len(self.aaData)-3):
            if self.aaData[i].p4_turn > 100 and self.aaData[i].p4_turn > self.aaData[i].p4_alpha \
                    and self.aaData[i].p4_turn > self.aaData[i].p4_beta and self.aaData[i].turn_prod > TURN_PRODUCT:
                self.aaData[i].turn_prop = True

    def secondary_structure(self, seqq):
        """
        The secondary_structure function assigns secondary structure
        symbol (E, H, T, C)->(Sheet, Helix, Turn, Random) for each aminoacid.
        :param seqq: input aminoacidic sequence
        :return: secondary structure sequence
        """
        secondary = ['C']*len(seqq)
        for i in range(len(seqq)-3):
            if self.aaData[i].helix_nucl == True and self.aaData[i].sheet_nucl == False \
                    and self.aaData[i].turn_prop == False:
                secondary[i] = 'H'
            elif self.aaData[i].helix_nucl == False and self.aaData[i].sheet_nucl == True \
                    and self.aaData[i].turn_prop == False:
                secondary[i] = 'E'
            elif self.aaData[i].helix_nucl == False and self.aaData[i].sheet_nucl == False \
                    and self.aaData[i].turn_prop == True:
                secondary[i] = 'T'
            elif self.aaData[i].helix_nucl == True and self.aaData[i].sheet_nucl == True \
                    and self.aaData[i].turn_prop == False:
                if self.aaData[i].p4_alpha >=  self.aaData[i].p4_beta:   ## >= tu jest slaby element...
                    secondary[i] = 'H'
                else:
                    secondary[i] = 'E'
    
        # secondary_structure cleaning
        for i in range(len(secondary) - 2):
            if secondary[i-1] == secondary[i+1] and secondary[i] != secondary[i-1]:
                secondary[i] = secondary[i-1]
            elif secondary[i] != secondary [i - 1] and secondary[i] != secondary[i + 1]:
                if secondary[i - 1] == 'C' and secondary[i + 1] in ('T', 'E'):
                    secondary[i] = secondary[i + 1]
                elif secondary[i + 1] == 'C' and secondary[i - 1] in ('T', 'E'):
                    secondary[i] = secondary[i - 1]
        # short helices removing
        sec = ''.join(secondary)
        sec1 = re.sub('([^H])HH([^H])', r'\1\1\2\2', sec)
        sec2 = re.sub('([^C])C([^C])', r'\1\1\2', sec1)
        sec3 = re.sub('([^H])HHH([^H])', r'\1HHHH', sec2)
        sec4 = list(re.sub('([^T])T([^T])', r'\1\1\2',sec3))
        return ''.join(sec4)
    
    def engine(self, seq_fasta):
        """
        This function execute Chou Fasman's algorithm
        :param seq_fasta: input sequence in fasta format
        :return: predicted secondary structure (string)
        """
        self.main(seq_fasta)
        self.get_probability()
        self.tetra_ave(seq_fasta)
        self.propagation()
        return self.secondary_structure(seq_fasta)

# seq = "TKYHDLKNKVIVIAGGIKNLGALTAKTFALESVNLVLHYHQAKDSDTANKLKDELEDQGAKVALYQSDLSNEEEVAKLFDFAEKEFGKVDIAINTVGKVLKKPIVETSEAEFDAMDTINNKVAYFFIKQAAKHMNPNGHIITIATSLLAAYTGFYST"
# engine(seq)
