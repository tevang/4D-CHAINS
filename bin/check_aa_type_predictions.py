#!/usr/bin/env python2.7
# 4D-CHAINS software is a property of Thomas Evangelidis and Konstantinos Tripsianes. The code is licensed under the Attribution-NonCommercial-NoDerivatives 4.0 International (CC BY-# NC-ND 4.0). You are free to:
# * Share - copy and redistribute the material in any medium or format.
# * The licensor cannot revoke these freedoms as long as you follow the license terms.
# Under the following terms:
# * Attribution - You must give appropriate credit, provide a link to the license, and indicate if changes were made. You may do so in any reasonable manner, but not in any way 
#   that suggests the licensor endorses you or your use.
# * NonCommercial - You may not use the material for commercial purposes.
# * NoDerivatives - If you remix, transform, or build upon the material, you may not distribute the modified material.
# * No additional restrictions - You may not apply legal terms or technological measures that legally restrict others from doing anything the license permits.
# To view a full copy of this license, visit https://creativecommons.org/licenses/by-nc-nd/4.0/legalcode.



import sys, re

code3 = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
code1 = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
aa1to3_dict = dict((c1,c3) for c1,c3 in zip(code1,code3))
aa3to1_dict = dict((c3,c1) for c3,c1 in zip(code3,code1))


predictions_file = sys.argv[1]
fasta = sys.argv[2]
start = int(sys.argv[3])

template_sequence_string = ""
with open(fasta, 'r') as f:
    for line in f:
        if not re.match("^>", line):
            template_sequence_string += re.sub(r"[^A-Z]", "", line)

resid2AAtype_dict = {}
for i, AA in enumerate(template_sequence_string):
    resid2AAtype_dict[start+i] = aa1to3_dict[AA]

missing_predictions_list = []
with open(predictions_file, 'r') as f:
    for line in f:
        mo = re.search("^([A-Z])([0-9]+)\s+", line)
        if mo:
            aa_type =aa1to3_dict[mo.group(1)]
            resid = int(mo.group(2))
            correct_prediction = resid2AAtype_dict[resid-1]
            no = re.search("\[(.*)\]", line)
            prediction_triplets = no.group(1).split(', ')
            CORRECT = False
            for triplet in prediction_triplets:
                if correct_prediction in triplet:
                    CORRECT = True
                    break
            if CORRECT == False:
                missing_predictions_list.append(aa_type+str(resid)+"-->"+resid2AAtype_dict[resid-1])


print "There are ", len(missing_predictions_list), " missing predictions: ", "\t".join(missing_predictions_list)