### run with singlecell conda env
### make sure to change gRNA_pilot_screen_file and pickle_out

import itertools
import numpy as np
import pandas as pd

gRNA_pilot_screen_file = "/gpfs/commons/home/imascio/lab/APA_pilot/PerturbSci/APA_1SRSF7_sgRNA_info_table.txt"

gRNA_0 = pd.read_csv(gRNA_pilot_screen_file, sep="\t", header=0)
gRNA = list(gRNA_0["gRNA_seq"])
print("Generating the dictionary for gRNA barcodes...")

gRNA_len = len(gRNA[0])
n_gRNA = len(gRNA)

gRNA_dict = dict()
gRNA_correction_black_list = set()

for each_gRNA in gRNA:
    ###add gRNA themselves to the dictionaries first
    gRNA_dict[each_gRNA] = each_gRNA
    
    ###then consider 1bp correction
    for each_letter in range(gRNA_len):
        letter_list = list(each_gRNA)
        for each_alternative in ["A", "G", "C", "T", "N"]:
            letter_list[each_letter] = each_alternative
            dist_seq = "".join(letter_list)
            
            ###only consider seqs which need corrections
            if dist_seq != each_gRNA:
                ###if the corrected gRNA sequence collides with any other sequence, just don't use any of them
                if dist_seq not in gRNA_dict.keys() and dist_seq not in gRNA_correction_black_list:
                    gRNA_dict[dist_seq] = each_gRNA
                else:
                    gRNA_correction_black_list.add(dist_seq)
                    if dist_seq in gRNA_dict.keys():
                        del gRNA_dict[dist_seq]

print(len(gRNA_dict))
print(len(gRNA_correction_black_list))
            
import pickle
pickle_out = open("/gpfs/commons/home/imascio/lab/APA_pilot/PerturbSci/APA_1SRSF7_sgRNAseq.pickle2", "wb")
pickle.dump(gRNA_dict, pickle_out, 2)
pickle_out.close()
