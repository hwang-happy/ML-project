import HMM
from ChouFasman import ChouFasman
import comparator

import matplotlib.pyplot as plt
import numpy as np

# path = "output102361.out"
path = "output49561.out"
ids, sequences, structures, skipped = [], [], [], []
scores_hmm_f_dssp, scores_hmm_back_dssp, scores_cf_dssp, scores_consensus_dssp, percentages = [], [], [], [], []

## alphabet_structure: List containing SSP secondary structure symbols
alphabet_structure = ['C', 'H', 'E', 'T']
## alphabet_proteins: List containing aminoacid IUPAC symbols
alphabet_proteins = ['A','R','N','D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

bad_1 = "B"
bad_2 = "X"
bad_3 = "J"
max = 20000

def run():
    i = 0
    with open(path) as f:
        for line in f:
            if i < max:
                if bad_1 in line.split(' ')[1].upper() or bad_2 in line.split(' ')[1].upper() or bad_3 in line.split(' ')[1].upper():
                    skipped.append(line.split(' ')[0])
                else:
                    ids.append(line.split(' ')[0])
                    sequences.append(line.split(' ')[1])
                    structures.append(line.split(' ')[2])
            i += 1

    length = len(ids)
    # print length
    for i in range(length):
        pdbID = ids[i]
        #print pdbID
        seqInput = sequences[i]
        #print i, pdbID
        CF = ChouFasman()
        ss_forward = ""
        ss_backward = ""
        ss_chou = ""
        consensus = ""

        ss_dssp = structures[i]
        ss_dssp_out = ""
        seqInput_out = ""
        diff = len(ss_dssp) - len(seqInput)
        # CHECKING LENGTH OF INPUT
        if diff > 0:
            for i in range(len(seqInput)):
                ss_dssp_out += ss_dssp[i]
                seqInput_out += seqInput[i]
        elif diff < 0:
            for i in range(len(ss_dssp)-1):
                ss_dssp_out += ss_dssp[i]
                seqInput_out += seqInput[i]

        else:
            ss_dssp_out = ss_dssp
            seqInput_out = seqInput

        ss_dssp_out = ss_dssp_out.upper()
        seqInput_out = seqInput_out.upper()

        ss_forward = HMM.mapping_num_to_letters(HMM.hmm_structure('forward', seqInput_out))
        ss_backward = HMM.mapping_num_to_letters(HMM.hmm_structure('backward', seqInput_out))
        ss_chou = ChouFasman.engine(CF, seq_fasta=seqInput_out)
        consensus = comparator.HMM_consensus([ss_forward,ss_backward,ss_chou])
        scores_hmm_f_dssp.append(comparator.identities([ss_forward, ss_dssp_out]))
        scores_hmm_back_dssp.append(comparator.identities([ss_backward, ss_dssp_out]))
        scores_cf_dssp.append(comparator.identities([ss_chou, ss_dssp_out]))
        scores_consensus_dssp.append(comparator.identities([consensus, ss_dssp_out]))
        percentages.append(comparator.compare_structures([ss_chou, ss_backward, ss_forward, ss_dssp_out, consensus]))

run()

# print scores_hmm_f_dssp
#
# import matplotlib.pyplot as plt
# import numpy as np
#
# x = np.arange(760)
#
# plt.plot(scores_hmm_f_dssp)
# plt.plot(scores_hmm_back_dssp)
# plt.plot(scores_cf_dssp)
# plt.plot(scores_consensus_dssp)
#
# plt.legend(['HMM forward', 'HMM backward', 'Chou Fasman', 'Consensus'], loc='upper left')
#
# plt.show()


# import matplotlib.pyplot as plt
# import numpy as np
#
# # Learn about API authentication here: https://plot.ly/python/getting-started
# # Find your api_key here: https://plot.ly/settings/api
#
# # gaussian_numbers = np.random.randn(1000)
# plt.hist(scores_hmm_f_dssp, edgecolor='None', alpha = 0.5, color= 'b')
# plt.hist(scores_hmm_back_dssp, edgecolor='None', alpha = 0.5, color= 'r')
# plt.hist(scores_cf_dssp, edgecolor='None', alpha = 0.5, color= 'k')
# plt.hist(scores_consensus_dssp, edgecolor='None', alpha = 0.5, color= 'g')
# plt.title("Precision by method")
# plt.xlabel("Precision")
# plt.ylabel("")
# plt.legend(['HMM forward', 'HMM backward', 'Chou Fasman', 'Consensus'], loc='best')
#
# fig = plt.gcf()
# plt.show()


def plot_methods_precisions():
    y = [np.mean(scores_hmm_f_dssp), np.mean(scores_hmm_back_dssp),
         np.mean(scores_cf_dssp), np.mean(scores_consensus_dssp)]
    N = len(y)

    ind = np.arange(N)  # the x locations for the groups
    width = 0.35       # the width of the bars

    fig, ax = plt.subplots()
    rects1 = ax.bar(ind, y, width, color='r')

    # add some text for labels, title and axes ticks
    ax.set_ylabel('Precision')
    ax.set_title('Mean precision by method')
    ax.set_xticks(ind + width / 2)
    ax.set_xticklabels(('HMM forward', 'HMM backward', 'Chou Fasman', 'Consensus'))
    for rect in rects1:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2., 1.01*height, np.around(height, 2),
                ha='center', va='bottom')
    plt.show()

plot_methods_precisions()
