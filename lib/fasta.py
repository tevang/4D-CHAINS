#!/usr/bin/env python
# 4D-CHAINS software is a property of is a property of Masaryk university and the authors are
# Thomas Evangelidis and Konstantinos Tripsianes. The code is licensed under the Attribution-NonCommercial-NoDerivatives 4.0
# International (CC BY-# NC-ND 4.0). You are free to:
# * Share - copy and redistribute the material in any medium or format.
# * The licensor cannot revoke these freedoms as long as you follow the license terms.
# Under the following terms:
# * Attribution - You must give appropriate credit, provide a link to the license, and indicate if changes were made. You may do so in any reasonable manner, but not in any way 
#   that suggests the licensor endorses you or your use.
# * NonCommercial - You may not use the material for commercial purposes.
# * NoDerivatives - If you remix, transform, or build upon the material, you may not distribute the modified material.
# * No additional restrictions - You may not apply legal terms or technological measures that legally restrict others from doing anything the license permits.
# To view a full copy of this license, visit https://creativecommons.org/licenses/by-nc-nd/4.0/legalcode.


import re
from .global_func import *
from lib.global_vars import *


class FASTA():

    def __init__(self, fasta_file, full_seq=False):
        """

        :param fasta_file:
        :param full_seq: do not omit the 1st residue (not visible in 4D-spectra)
        """
        # Load the fasta sequence to a list
        protseq_str = ""
        self.protname = ""
        with open(fasta_file, 'r') as f:
            for line in f:
                if re.match("^>", line):
                    self.protname = line.replace(">", "")
                else:
                    protseq_str += re.sub(r"[^A-Z]", "", line)

        self.protname = self.protname.rstrip()
        protseq_str = protseq_str.rstrip()

        print(protseq_str)
        if full_seq:
            self.protseq_list = list(protseq_str)
            self.protseq_str = protseq_str
        else:
            self.protseq_list = list(protseq_str[1:])   # WARNING: may be dangerous to omit the first residue.
            self.protseq_str = protseq_str[1:]

    def count_AAIGs(self):
        """
        HN-N Peaks (AAIGs) expected from sequence:
        XXX backbone (all AA minus PRO minus FIRST)
        XXX SIDECHAIN pairs (ASN+GLN)
        XXX SIDECHAIN (TRP)

        :return:
        """

        PRO = self.protseq_list.count('P')
        AA = len(self.protseq_list) - PRO
        ASN = self.protseq_list.count('N')
        GLN = self.protseq_list.count('Q')
        TRP = self.protseq_list.count('W')
        self.BB_AAIG_num = AA
        self.SC_AAIG_num = 2*(ASN + GLN) + TRP
        self.SC_NQ_num = 2*(ASN + GLN)  # each sidechain has NH2 (two protons)
        self.SC_W_num = TRP             # each sidechain has NH (one proton)

    def count_AA_frequencies(self):
        """

        :return: aatype_P_dict: the percentage of each AA type in the protein sequence
        """
        aatype_P_dict = {}  # dict with the P[aatype(i-1)] for every aatype(i-1): the percentage of each aatype in the
                            # template sequence
        for aa_type in list(aatype_maxH_C_pairs_dict.keys()):
            aatype_P_dict[aa_type] = self.protseq_list.count(aa3to1_dict[aa_type]) / \
                                     float(len(self.protseq_list))
        return aatype_P_dict