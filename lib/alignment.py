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


#!/usr/bin/env python
import os
import pickle
import re
import shutil
import sys
from collections import OrderedDict
from operator import itemgetter

from ete3 import Tree
from lib.alignment_v2 import print_results_summary
from lib.fasta import FASTA
from lib.global_func import Debuginfo, tree, remove_NH_suffix, save_pickle
from scoop import shared
from tabulate import tabulate
import numpy as np
import pandas as pd

class Alignment():

    def __init__(self,
                 ABSOLUTE_MATCHES_FILE,
                 FIRST_RESIDUE_NUMBER,
                 spectrum_combo="TOCSY-HCNH"):
        self.ABSOLUTE_MATCHES_FILE, self.FIRST_RESIDUE_NUMBER, self.spectrum_combo = \
            ABSOLUTE_MATCHES_FILE, FIRST_RESIDUE_NUMBER, spectrum_combo
        self.absolute_AAIGmatches_alignment_list = []   # alignment of AAIG signatures to the protein sequence
        self.protein_alignment_list = []    # the protein sequence in a list form
        self.absolute_matches_alignment_list = []   # alignment of real residue names, not AAIG signatures

    @staticmethod
    def read_NHmap_file(absolute_matches_file, get_protein_alignment=False):
        ## READ FILE WITH ABSOLUTE MATCHES
        protein_alignment_list = []
        absolute_matches_alignment_list = []
        with open(absolute_matches_file, 'r') as f:
            for line in f:
                if line[0:53] == "CONSENSUS CONFIDENCE TOP-SCORED AND ABSOLUTE MATCHES:" or \
                        line[0:36] == "ABSOLUTE MATCHES WITH CHAIN LINKERS:":
                    try:
                        while True:
                            line = next(f)
                            while line[0] == '\n':
                                line = next(f)
                            if line.startswith("STATISTICS:"):  # this is where the NH mapping stops, don't read in more
                                break
                            line = next(f)
                            # print "DEBUG: Protein:", line
                            protein_alignment_list.extend(line.replace('│', '').split())
                            line = next(f)
                            line = next(f)
                            # print "DEBUG: Matches:", line
                            absolute_matches_alignment_list.extend(line.replace('│', '').replace('*', '-').split())
                            line = next(f)
                            line = next(f)
                            line = next(f)
                            line = next(f)
                    except StopIteration:
                        break

        # protein_alignment_list = protein_alignment_list[:-1]    # remove the last element which is 'N/A'

        if get_protein_alignment:
            return absolute_matches_alignment_list, protein_alignment_list  # include the last element: 'N/A'
        else:
            return absolute_matches_alignment_list

    def create_absolute_AAIGmatches_alignment(self):
        """
        Method to READ FILE WITH ABSOLUTE MATCHES (NHmap)
        :param ABSOLUTE_MATCHES_FILE:
        :populate: self.absolute_AAIGmatches_alignment_list:     contains the absolutely matched root index groups (AAIG),
                                                                namely the names given to each N, HN peak in the HSQC.
                                                                It does not necessarily contain real residue names!
        """

        self.absolute_AAIGmatches_alignment_list, self.protein_alignment_list = \
            Alignment.read_NHmap_file(self.ABSOLUTE_MATCHES_FILE, get_protein_alignment=True)
        # Remove terminal 'N/A' otherwise they will cause you trouble later
        NA_indices = []
        if self.protein_alignment_list[0] == 'N/A' or self.absolute_AAIGmatches_alignment_list[0] == 'N/A':
            NA_indices.append(0)
            self.spectrum_combo = "TOCSY-HCNH"
        else:
            self.spectrum_combo = "HCNH-HCNH"
        if self.protein_alignment_list[-1] == 'N/A' or self.absolute_AAIGmatches_alignment_list[-1] == 'N/A':
            NA_indices.append(len(self.protein_alignment_list) - 1)
        # remove the 'N/A' elements from self.absolute_AAIGmatches_alignment_list and self.protein_alignment_list
        self.protein_alignment_list = [self.protein_alignment_list[i] for i in range(len(self.protein_alignment_list))
                                  if i not in NA_indices]
        self.absolute_AAIGmatches_alignment_list = [self.absolute_AAIGmatches_alignment_list[i]
                                                   for i in range(len(self.absolute_AAIGmatches_alignment_list))
                                                   if i not in NA_indices]

    def create_absolute_matches_alignment(self):
        """
        Method to be called after absolute_AAIGmatches_alignment() in order to rename the abstract AAIG signatures in
        the self.absolute_AAIGmatches_alignment_list to correct residue names according to the protein sequence and
        starting resid.

        :param spectrum_combo:
        :return:    populates self.absolute_matches_alignment_list
        """
        if self.spectrum_combo == "TOCSY-HCNH":
            resid_offset = 1
        elif self.spectrum_combo == "HCNH-HCNH":
            resid_offset = 0
        print("DEBUG: resid_offset=", resid_offset)
        print("DEBUG: self.FIRST_RESIDUE_NUMBER=", self.FIRST_RESIDUE_NUMBER)
        resid = self.FIRST_RESIDUE_NUMBER  # start counting from the 1st residue number specified by the user
        self.absolute_matches_alignment_list = ['-'] * len(self.absolute_AAIGmatches_alignment_list)  # contains the absolutely matched residue names,
                                                      # namely real residue names not AAIGs!
        for position in range(len(self.protein_alignment_list)):
            resname = self.protein_alignment_list[position]
            AAIG = self.absolute_AAIGmatches_alignment_list[position]
            if AAIG == '-':
                self.absolute_matches_alignment_list[position] = '-'
            elif AAIG == 'N/A':
                self.absolute_matches_alignment_list[position] = 'N/A'
            else:
                residue = resname + str(resid + resid_offset)
                self.absolute_matches_alignment_list[position] = residue
            resid += 1

    def write_alignment_to_file(self,
                                outfname,
                                protein_alignment_list,
                                new_absolute_AAIGmatches_alignment_list,
                                i_iminus1_dict={},
                                patched_TIGs_list=[]):

        # CREATE THE CONNECTIVITIES ALIGNMENT
        print("DEBUG: i_iminus1_dict=", i_iminus1_dict)
        connectivities_alignment = []
        if i_iminus1_dict:  # if connectivities have been provided
            for i in range(len(new_absolute_AAIGmatches_alignment_list) - 1):
                TIG_i = new_absolute_AAIGmatches_alignment_list[i]
                TIG_iplus1 = new_absolute_AAIGmatches_alignment_list[i + 1]
                if TIG_i != '-' and TIG_iplus1 != '-':
                    try:
                        connectivity = [c for c in i_iminus1_dict[TIG_iplus1] if c[0] == TIG_i][
                            0]  # connectivity = (TIG_i, occupancy, TOCSY_resonnum)
                        print("DEBUG: connectivity=", connectivity)
                        if TIG_iplus1 in patched_TIGs_list:  # if it has been patched, write "0/0" to distinguish it
                            connectivities_alignment.append("0/0")
                        else:
                            connectivities_alignment.append(str(connectivity[1]) + "/" + str(connectivity[2]))
                    except (IndexError, KeyError):
                        connectivities_alignment.append('-')
                else:
                    connectivities_alignment.append('-')
            connectivities_alignment.append('-')
        connectivities_alignment.insert(0, '-')
        connectivities_alignment[0] = 'N/A'

        # CALCULATE NH-MAPPING STATISTICS
        CORRECT, WRONG, COVERAGE = 0, 0, 0
        wrong_AAIG_list = []
        for position in range(1, len(protein_alignment_list)-1):
            aa_type = protein_alignment_list[position]
            pred_aa_type = new_absolute_AAIGmatches_alignment_list[position][0]
            if pred_aa_type != '-':
                if pred_aa_type == aa_type:
                    CORRECT += 1
                else:
                    WRONG += 1
                    wrong_AAIG_list.append(new_absolute_AAIGmatches_alignment_list[position])
                COVERAGE += 1
        COVERAGE = 100.0*COVERAGE/len(protein_alignment_list)


        with open(outfname, 'w') as f:
            f.write("\n\n\nABSOLUTE MATCHES WITH CHAIN LINKERS:\n\n\n" + "\n")
            print("\n\n\nABSOLUTE MATCHES WITH CHAIN LINKERS:\n\n\n")
            # print "DEBUG: new_absolute_matches_alignment_list=", new_absolute_matches_alignment_list
            # print "DEBUG: protein_alignment_list=", protein_alignment_list
            for start, end in zip(list(range(0, len(new_absolute_AAIGmatches_alignment_list), 17)),
                                  list(range(17, len(new_absolute_AAIGmatches_alignment_list), 17))):
                headers = protein_alignment_list[start:end]
                table = [new_absolute_AAIGmatches_alignment_list[start:end]]
                table.append(connectivities_alignment[start:end])
                print(tabulate(table, headers, tablefmt="fancy_grid"))
                f.write(tabulate(table, headers, tablefmt="fancy_grid") +
                        "\n")  # convert the \n to the same format (binary), too.
                print("\n\n")
                f.write("\n\n" + "\n")
            headers = protein_alignment_list[end:]
            table = [new_absolute_AAIGmatches_alignment_list[end:]]
            table.append(connectivities_alignment[end:])
            print(tabulate(table, headers, tablefmt="fancy_grid"))
            f.write(tabulate(table, headers, tablefmt="fancy_grid") +
                    "\n")  # convert the \n to the same format (binary), too.
            # WRITE THE STATISTICS OF NH-MAPPING AT THE END OF THE FILE
            f.write("\n\n\nSTATISTICS:\n\n")
            print("COVERAGE = " + str(COVERAGE) + " %")
            f.write("COVERAGE = " + str(COVERAGE) + " %\n")
            print("CORRECT = " + str(CORRECT))
            f.write("CORRECT = " + str(CORRECT) + "\n")
            print("WRONG = " + str(WRONG) + " (" + ", ".join(wrong_AAIG_list) + ")")
            f.write("WRONG = " + str(WRONG) + " (" + ", ".join(wrong_AAIG_list) + ")\n")

def get_alignments_with_unique_AAIG(consensus_overlappingChains_alignment_list, max_peptide_length):
    """
        FUNCTION to break long consensus sequence alignments with multiple occurrences of one or more TOCSY indices
        into smaller ones with no duplications.
    """
    print("DEBUG: consensus_overlappingChains_alignment_list=", consensus_overlappingChains_alignment_list)
    ungapped_consensus_overlappingChains_alignment_list = [a for a in consensus_overlappingChains_alignment_list if a != '-']
    
    all_subsequences_list = []
    subsequence = []
    for item in ungapped_consensus_overlappingChains_alignment_list:
        if item in subsequence:
            if len(subsequence) >= max_peptide_length:
                #print subsequence
                all_subsequences_list.append(subsequence)
            index = subsequence.index(item)+1
            subsequence = subsequence[index:]
        subsequence.append(item)
    if len(subsequence) >= max_peptide_length:
        #print subsequence
        all_subsequences_list.append(subsequence)
    
    #print all_subsequences_list
    subsequence_alignments_list = []
    for subsequence in all_subsequences_list:
        subsequence_alignment = ["-"] * len(consensus_overlappingChains_alignment_list)
        subseq_length = len(subsequence)
        for index in range(len(consensus_overlappingChains_alignment_list)):
            num_of_matches = 0
            for offset in range(subseq_length):
                if consensus_overlappingChains_alignment_list[index + offset] == subsequence[offset]:
                    num_of_matches += 1
            if num_of_matches == subseq_length:
                start = index
                end = index + offset
                for position, TAAIG in zip(list(range(start, end +1)), subsequence):
                    subsequence_alignment[position] = TAAIG
                break
        subsequence_alignments_list.append(subsequence_alignment)
    
    return subsequence_alignments_list

    
def get_alignment_with_absolute_matches(consensus_sequence_alignments_set, args):
    
    print("DEBUG: consensus_sequence_alignments_set=", consensus_sequence_alignments_set)
    alignment_length = len(next(iter(consensus_sequence_alignments_set)))
    
    # IF WE HAVE RESTRAINTS, FILL IN THE POSITIONS THAT HAVE GAPS IN THE absolute_matches_alignment WITH THE RESPECTIVE AAIG OF THE RESTRAINTS ALIGNMENT FILE
    if args.ABSOLUTE_MATCHES_FILE:
        raw_restraints_alignment = Alignment.read_NHmap_file(args.ABSOLUTE_MATCHES_FILE)
        restraints_alignment = ['*' if x == '-' else x for x in raw_restraints_alignment[1:]]   # remove 'N/A' at the beginning and convert '-' to '*'
        del raw_restraints_alignment
        #print "DEBUG: restraints_alignment =", restraints_alignment
    
    # PLACE THE CONSENSUS TAAIG/gap OF ALL CONTIGS AT EVERY POSITION OF THE absolute_matches_alignment
    absolute_matches_alignment = ["*"] * alignment_length
    for position in range(alignment_length):
        try:
            AAIG_list = []
            for chain_alignment_list in consensus_sequence_alignments_set:
                AAIG_list.append(chain_alignment_list[position])
            AAIG_list = [a for a in AAIG_list if a != '-']   # remove gaps from the list
            if len(set(AAIG_list)) == 1:    # if we found a unique TAAIG at this position, save it into absolute_matches_alignment
                absolute_matches_alignment[position] = AAIG_list[0]
            elif len(set(AAIG_list)) == 0 and args.ABSOLUTE_MATCHES_FILE:  # if according the our contigs there is a gap there, place whatever the restraint alignment has at that position (TAAIG or '-')
                absolute_matches_alignment[position] = restraints_alignment[position]
        except IndexError:
            print("DEBUG: position=", position)
            print("DEBUG: absolute_matches_alignment=", absolute_matches_alignment)
            if args.ABSOLUTE_MATCHES_FILE: print("DEBUG: restraints_alignment =", restraints_alignment)
            sys.exit(1)
        
    #print "DEBUG: raw absolute_matches_alignment = ", absolute_matches_alignment
    
    ## CLEAN THE ABSOLUTE MATCHES ALIGNMENT FROM DUPLICATE TOCSY INDICES
    Tdindex2position_dict = {}
    duplicateTAAIG_indices_dict = {}
    for TAAIG, position in zip(absolute_matches_alignment, list(range(len(absolute_matches_alignment)))):
        if TAAIG == '*':
            continue
        if TAAIG in list(Tdindex2position_dict.keys()) and position != Tdindex2position_dict[TAAIG]:
            indices = [i for i, x in enumerate(absolute_matches_alignment) if x == TAAIG]
            duplicateTAAIG_indices_dict[TAAIG] = indices
        else:
            Tdindex2position_dict[TAAIG] = position
    
    # NOW CLEAN THE GROUPS OF DUPLICATE TOCSY INDICES FROM OVERHANGING ENDS AND KEEP ONLY THOSE THAT CONTAIN ONE NON-OVERHANGING AAIG
    #print "DEBUG: absolute_matches_alignment before=", absolute_matches_alignment
    for TAAIG, index_group in list(duplicateTAAIG_indices_dict.items()):
        #print "DEBUG: index_group=", index_group
        index_locationSet_dict = {}
        for index in index_group:
            index_locationSet_dict[index] = set()
            for chain_alignment_list in consensus_sequence_alignments_set:
                if chain_alignment_list[index] == TAAIG:   # just found a consensus sequence alignment with the current TAAIG
                    if index == 0:  # in case it is N-term of C-term it's still an overhanging end ????
                        index_locationSet_dict[index].add("N-term")
                    elif index == len(chain_alignment_list)-1:
                        index_locationSet_dict[index].add("C-term")
                    elif chain_alignment_list[index-1] != "-" and chain_alignment_list[index+1] != "-":
                        print(chain_alignment_list)
                        index_locationSet_dict[index].add("middle")
                    elif chain_alignment_list[index-1] == "-" and chain_alignment_list[index+1] != "-":
                        index_locationSet_dict[index].add("N-term")
                    elif chain_alignment_list[index-1] != "-" and chain_alignment_list[index+1] == "-":
                        index_locationSet_dict[index].add("C-term")
        # SAVE THE INDICES THAT ARE OVERHANGING ENDS
        overhanging_ends_list = []
        for index, location_set in list(index_locationSet_dict.items()):
            #print "DEBUG: TAAIG=", TAAIG, "index=", index, "location_set=", location_set
            if len(location_set) == 1 and not "middle" in location_set:
                overhanging_ends_list.append(index)
        #print "DEBUG: overhanging_ends_list=", overhanging_ends_list
        if len(overhanging_ends_list) == len(index_group)-1:   # if there was only one position that was not overhanging end, keep the AAIG only at that position
            for index in overhanging_ends_list:
                absolute_matches_alignment[index] = "*"
        else:   # if there were multiple positions that were not overhanging ends, remove all instances of that AAIG
            for index in index_group:
                absolute_matches_alignment[index] = "*"
    
    #print "DEBUG: absolute_matches_alignment  after=", absolute_matches_alignment
    ### RECOVER THE CHAIN ALIGNMENT STRING
    #chain_line = '     '
    #for word in absolute_matches_alignment:
    #    if word == '*':
    #        chain_line += '  *  '
    #    else:
    #        chain_line += "%12s" % word
    #print ''.join(chain_line)
    #return ''.join(chain_line) + "\n"
    
    ## RECOVER THE CHAIN ALIGNMENT LIST
    chain_list = [' N/A ']
    for word in absolute_matches_alignment:
        if word == '*':
            chain_list.append('*')
        else:
            chain_list.append("%12s" % word)
    return chain_list


def is_valid_alignment(peptide_aln, chain_aln, max_peptide_length):
    """
        FUNCTION to check for some exceptions in overlapping alignment filtering for TAAIG prediction scores. 
    """
    peptide_aln_list = peptide_aln.split() # this was string, convert it to list like chain_aln
    peptide_list = [a for a in peptide_aln_list if a != '-']
    chain_list = [a for a in chain_aln if a != '-']
    #if "I145" in chain_aln and "C146" in chain_aln and "E147" in chain_aln:
    #    print "DEBUG: peptide_aln_list=", peptide_aln_list
    #    print "DEBUG: chain_aln=", chain_aln
    #    sys.exit(1)
    if peptide_list[0] == 'P':  # if Proline is the N-term of the peptide
        return True
    #elif len(peptide_list) <= max_peptide_length and (peptide_aln_list[0] != '-' or chain_aln[-2] != '-'):
    #    return True
        
    return False
    

def create_consensus_alignemnt(peptide_identity_dict,
                               peptide_score_dict,
                               peptideName2chain_dict,
                               peptide_alignment_dict,
                               peptideName_chainScoreList_dict,
                               args):
    """
        FUNCTION to process the raw data that resulted from loading the alignment files and to create various alignment forms
        in useful format.
        
        ARGUMENTS:
        peptide_identity_dict:          peptide name --> seq. identity to the protein sequence
        peptide_score_dict:             peptide name --> alignment score
        peptideName2chain_dict:         peptide name --> source TAAIG chain
        peptide_alignment_dict:         peptide_name --> (template_alignment string, peptide_alignment string)
        peptideName_chainScoreList_dict: peptide_name --> list of confidence scores for each position of the peptide, e.g. [5.6381008762e-05, 5.61539400189e-05, 3.28593502441e-05]
        
        RETURNS:
        ordered_consensus_peptide_alignment_list:       list of strings of the form:  '-        -     ...   G    L    S         -           -     '
        ordered_consensus_chain_alignment_list:         list of lists of the form: ['-', '-', 'T252NH', 'S267NH', 'L268NH', '-', '-']
        ordered_consensus_TAAIGScore_alignment_list:   list of lists of the form: [('', 0.0), ('', 0.0), ... , ('T252NH', 0.01318430523215), ('S267NH', 7.53408792136e-10), ('L268', 0.008264462809928), ('', 0.0), ('', 0.0)]
        template_sequence_alignment_string:             string of the protein in the form: "   A    N    I    V    G    G    I ...  L    V    T    G  "
        
    """

    ## CREATE THE CONSENSUS ALIGNMENT
    consensus_peptideName_list = []     # a list to keep record of the peptide names in the order they occur in the consensus alignment
    consensus_peptide_alignment_list = []
    consensus_chain_alignment_list = []
    consensus_TAAIGScore_alignment_list = []
    #print '  ' + '  '.join(peptide_alignment_dict.values()[0][0]) + '  '
    template_sequence_alignment_string = ('  ' + '    '.join(list(peptide_alignment_dict.values())[0][0].replace('-', '')) + '  ')
    predictionScoreTuples_list = []    # alignment position -> list of tuples of the form (aa type, score)
    #print "DEBUG: identity",len(peptide_identity_dict)
    #print "DEBUG: peptide_score_dict=",len(peptide_score_dict)
    #print "DEBUG: peptide_alignment_dict=", len(peptide_alignment_dict)
    for peptide_name, alignments in list(peptide_alignment_dict.items()):
        peptide_length = int(peptide_name.split('_')[1])
        chainscore_list = peptideName_chainScoreList_dict[peptide_name]
        peptide_line = ""
        chain_line = "     "
        TAAIGScoreTuple_list = [] # list of tuples of the form (AAIG, confidence score) that has the same alignment as chain alignment
        DISCARD_PEPTIDE = False
        if peptide_identity_dict[peptide_name] >= args.IDENTITY_CUTOFF:                  # <== CHANGE ME
            chain = peptideName2chain_dict[peptide_name].split('-')
            chain_index = -peptide_length
            for aln_index in range(len(alignments[1])):
                try:
                    template_char = alignments[0][aln_index+1]
                except IndexError:
                    template_char = ""
                peptide_char = alignments[1][aln_index]
                if peptide_char == '-':
                    peptide_line += '      -     '
                    chain_line += '      -     '
                    TAAIGScoreTuple_list.append(("", 0.0))
                elif template_char == 'P':     # because a real proline cannot give the aa type of the previous residue (NOT NECESSARY FILTER)
                    DISCARD_PEPTIDE = True
                    break
                else:
                    peptide_line += '  '+peptide_char+'  '
                    #print "DEBUG: chain=", chain
                    #print "DEBUG: chain_index=", chain_index
                    #print "DEBUG: peptide_name=", peptide_name
                    chain_line += "%12s" % chain[chain_index]
                    ## WEIGHT THE CHAIN SCORE BY THE ALIGNMENT SCORE OF THE PEPTIDE (LONGER PEPTIDES HAVE HIGHER WEIGHTS)
                    TAAIGScoreTuple_list.append((chain[chain_index], chainscore_list[chain_index] * peptide_score_dict[peptide_name]))
                    predictionScoreTuples_list
                    chain_index += 1
            
            print(peptide_line)
            print(chain_line)
            IS_PEPTIDE_DUPLICATE = False
            for pepindex in range(len(consensus_peptide_alignment_list)):
                if peptide_line == consensus_peptide_alignment_list[pepindex] and chain_line == consensus_chain_alignment_list[pepindex]:
                    IS_PEPTIDE_DUPLICATE = True
                    break
            if not IS_PEPTIDE_DUPLICATE and not DISCARD_PEPTIDE:
                consensus_peptideName_list.append(peptide_name)
                consensus_peptide_alignment_list.append(peptide_line)
                consensus_chain_alignment_list.append(chain_line)
                consensus_TAAIGScore_alignment_list.append(TAAIGScoreTuple_list)
    
    ## CLEAN THE peptide_score_dict FROM MISMATCHED PEPTIDES
    for peptide_name in list(peptide_score_dict.keys()):
        if not peptide_name in consensus_peptideName_list:
            del peptide_score_dict[peptide_name]
    
    ## REORDER THE CHAIN ALIGNMENT TO VISUALIZE OVERLAPS EASIER
    global ordered_consensus_peptide_alignment_list, ordered_consensus_chain_alignment_list, ordered_consensus_TAAIGScore_alignment_list    
    ordered_consensus_peptide_alignment_list = []
    ordered_consensus_chain_alignment_list = []
    ordered_consensus_TAAIGScore_alignment_list = []
    #print "DEBUG: consensus_chain_alignment_list[0]=",consensus_chain_alignment_list[0]
    #print "DEBUG: template_sequence_alignment_string =", len(template_sequence_alignment_string.split())
    used_indices_set = set()    # set with saved alignments to avoid saving them in multiple start_column
    for start_column in range(len(consensus_chain_alignment_list[0].split())-2):  # the minimum peptide length is 2 so end at position -2
        #print "DEBUG: start_column = ", start_column)
        index = 0
        for chain_alignment_string, peptide_name, peptide_line, TAAIGScoreTuple_list in zip(consensus_chain_alignment_list, consensus_peptideName_list, consensus_peptide_alignment_list, consensus_TAAIGScore_alignment_list):
            chain_alignment = chain_alignment_string.split()
            peptide_length = int(peptide_name.split('_')[1])
            if not '-' in chain_alignment[start_column:(start_column+peptide_length)] and not index in used_indices_set:    # append every alignment only once
                #print "DEBUG: saving chain", chain_alignment[start_column:(start_column+peptide_length)]
                #print "DEBUG: saving chain_alignment", chain_alignment
                ordered_consensus_peptide_alignment_list.append(peptide_line)
                ordered_consensus_chain_alignment_list.append(chain_alignment)
                ordered_consensus_TAAIGScore_alignment_list.append(TAAIGScoreTuple_list)
                used_indices_set.add(index)
            index += 1
    
    #print "DEBUG: len(ordered_consensus_chain_alignment_list)",len(ordered_consensus_chain_alignment_list)," == len(consensus_chain_alignment_list)",len(consensus_chain_alignment_list)
    return ordered_consensus_peptide_alignment_list, ordered_consensus_chain_alignment_list, ordered_consensus_TAAIGScore_alignment_list, template_sequence_alignment_string


def write_consensus_alignment(ordered_consensus_peptide_alignment_list, ordered_consensus_chain_alignment_list, ordered_consensus_TAAIGScore_alignment_list,
                              template_sequence_alignment_string, args, iteration=1):
    """
        FUNCTION to write the ordered chain alignments into file consensus_alignment.iteration[1-7].
    """
    
    ## NOW PRINT THE CONSENSUS ALIGNMENT
    with open(args.consensus_alignment_file+".iteration"+str(iteration), 'w') as f:
        total_chars_per_line = len(template_sequence_alignment_string)
        start = 0
        while (start+args.chars_per_line-1) <= total_chars_per_line:
            print("\n\n")
            print(template_sequence_alignment_string[start:(start+args.chars_per_line-1)])
            f.write(template_sequence_alignment_string[start:(start+args.chars_per_line-1)] + "\n")
            for peptide_line, chain_alignment_list, TAAIGScoreTuple_list in zip(ordered_consensus_peptide_alignment_list, ordered_consensus_chain_alignment_list, ordered_consensus_TAAIGScore_alignment_list):
                #print peptide_line[start:(start+args.chars_per_line-1)]
                #f.write(peptide_line[start:(start+args.chars_per_line-1)] + "\n")
                ## RECOVER THE CHAIN ALIGNMENT STRING
                chain_line = '     '
                for word in chain_alignment_list:
                    if word == '-':
                        chain_line += '      -     '
                    else:
                        chain_line += "%12s" % word
                print(''.join(chain_line[start:(start+args.chars_per_line-1)]))
                f.write(''.join(chain_line[start:(start+args.chars_per_line-1)]) + "\n")
                #print TAAIGScoreTuple_list[start:(start+args.chars_per_line-1)])
                #f.write(TAAIGScoreTuple_list[start:(start+args.chars_per_line-1)] + "\n")
            start += args.chars_per_line-1
        print("\n\n")
        f.write("\n\n\n")
        print(template_sequence_alignment_string[start:])
        f.write(template_sequence_alignment_string[start:] + "\n")
        for peptide_line, chain_alignment_list, TAAIGScoreTuple_list in zip(ordered_consensus_peptide_alignment_list, ordered_consensus_chain_alignment_list, ordered_consensus_TAAIGScore_alignment_list):
            #print peptide_line[start:])
            #f.write(peptide_line[start:] + "\n")
            ## RECOVER THE CHAIN ALIGNMENT STRING
            chain_line = '     '
            for word in chain_alignment_list:
                if word == '-':
                    chain_line += '      -     '
                else:
                    chain_line += "%12s" % word
            print(''.join(chain_line[start:]))
            f.write(''.join(chain_line[start:]) + "\n")
            #print TAAIGScoreTuple_list[start:])
            #f.write(TAAIGScoreTuple_list[start:] + "\n")
    
    
## BUILD TREES OF OVERLAPPING CHAIN ALIGNMENTS
def populate_leaves(Alignment_Tree, refchain_overlappingChainsList_dict):
    """
        INNER FUNCTION that adds new branches to the leaves of the Tree.
        ARGUMENTS:
        Alignment_Tree:    The Tree structure with connectivities
        RETURNS:
        (Alignment_Tree, BOOLEAN):    A tuple with elements the input Tree structure with new branches (if applicable), and a BOOLEAN value which is True if the function added
                                       new leaves to the Tree, or False otherwise
    """
    
    number_of_new_leaves = 0
    # ATTENTION: never use Assignment_Tree.iter_leaf_names(), it doesn't return the names in the order
    # ATTENTION: corresponding to Assignment_Tree.get_leaves()!!!
    for leaf in Alignment_Tree.get_leaves():
        try:
            #for child_tuple in refchain_overlappingChainsList_dict[name]:
            for query_aligned_chain in refchain_overlappingChainsList_dict[tuple(leaf.name)]:    # convert the chain alignment lists to tuples
                query_chain = [a for a in query_aligned_chain if a != '-']
                #ancestors_list = [ancestor.name for ancestor in leaf.get_ancestors()]
                #if query_aligned_chain in ancestors_list:     # if the current NOESY AAIG is already a node or leaf in the Tree, continue to the next
                #    continue
                new_child = leaf.add_child(name=query_aligned_chain, dist=1) # add a new brach
                new_child.add_features(chain=query_chain)
                number_of_new_leaves += 1
                #print "DEBUG: adding connection: ",name,"-->",NOESYaaindex
        except KeyError:
            continue
    
    #print Alignment_Tree.get_ascii(show_internal=True, compact=False, attributes=["name", "dist", "occupancy", "numOfResonances"])
    #print Alignment_Tree.get_ascii(show_internal=True, compact=False)
    if number_of_new_leaves > 0:
        return (Alignment_Tree, True)
    else:
        return (Alignment_Tree, False)


def build_Alignment_Tree(reference_aligned_chain, refchain_overlappingChainsList_dict, ordered_consensus_protein_alignment_list,
                            ordered_consensus_chain_alignment_list):
    """
    FUNCTION to create contigs of overlapping chain alignments.
    
    ARGUMENTS:
    reference_aligned_chain:                one chain alignment (to be used as the root of the tree) from which it will start forming contigs of overlapping chain alignments.
    refchain_overlappingChainsList_dict:    ordereddict with keys the reference chain alignment that have at least 1 overlapping chain alignment, and values tuples of the
                                            overlapping alignments
    """
    
    
    expand_tree = True
    Alignment_Tree = Tree()
    Root = Alignment_Tree.get_tree_root()
    reference_chain = [a for a in reference_aligned_chain if a != '-']
    Root.add_features(name=reference_aligned_chain, chain=reference_chain)
    level = 1
    sys.stdout.write("Expanding tree from level ")
    while expand_tree:
        sys.stdout.write(str(level)+" ")
        sys.stdout.flush()
        Alignment_Tree, expand_tree = populate_leaves(Alignment_Tree, refchain_overlappingChainsList_dict)
        level += 1
        #if level == args.MAX_PEPTIDE_LENGTH:
        #    break
    # Print the Tree
    #print Alignment_Tree.get_ascii(show_internal=True, compact=False)
    #print Alignment_Tree.get_ascii(show_internal=True, compact=False, attributes=["name", "dist", "occupancy", "numOfResonances"])
    
    print("\nSaving overlapping chain alignments from Tree...")
    
    all_overlapped_chain_alignments_list = []   # this is a local list to this function and is irrelevant to the all_overlapped_chain_alignments_list in create_overlapped_chain_alignments()
    for leaf in Alignment_Tree.get_leaves():
        overlapped_chain_alignments_list = []
        overlapped_chain_alignments_list.append(leaf.name)
# DEACTIVATED            #consensus_overlapping_sequence = []     # store the consensus overlapping sequence so far to prevent the occurrence of duplicate AAIG
# DEACTIVATED            #consensus_overlapping_sequence.extend(leaf.chain)
        for ancestor in leaf.get_ancestors():
            # TODO: check if including the Root to the overlapped_chain_alignments_list is correct.
# DEACTIVATED                #if ancestor.chain[0] in consensus_overlapping_sequence:   # if the overhanging N-term is in the consensus sequence already (recall that we go the opposite direction here from leaf->root)
# DEACTIVATED                #    print "DEBUG: first element of ancestor.chain", ancestor.chain, "is in consensus_overlapping_sequence",consensus_overlapping_sequence)
# DEACTIVATED                #    break
# DEACTIVATED                #consensus_overlapping_sequence.extend(ancestor.chain[-1])
            overlapped_chain_alignments_list.append(ancestor.name)
        if len(overlapped_chain_alignments_list) > 1:   # save only if we have at least 2 overlapping chains
            overlapped_chain_alignments_tuple = tuple(tuple(x) for x in reversed(overlapped_chain_alignments_list)) # reverse the order to the correct
            all_overlapped_chain_alignments_list.append(overlapped_chain_alignments_tuple)
        elif len(overlapped_chain_alignments_list) == 1: # check if a Proline is in the N-term of the respective peptide, otherwise don't save the chain
            chain_alignment_list = overlapped_chain_alignments_list[0]
            position = ordered_consensus_chain_alignment_list.index(list(chain_alignment_list))
            peptide_alignment_list = ordered_consensus_protein_alignment_list[position]
            peptide_sequence_list = [a for a in peptide_alignment_list.split() if a != '-']
            print(peptide_sequence_list)
            if peptide_sequence_list[0] == "P": # if the peptide starts from Pro there is not way to find an overlapping peptide, so save the respective chain
                overlapped_chain_alignments_tuple = tuple(tuple(x) for x in reversed(overlapped_chain_alignments_list)) # reverse the order to the correct
                all_overlapped_chain_alignments_list.append(overlapped_chain_alignments_tuple)
    
    return all_overlapped_chain_alignments_list


def check_if_subgroup(overlappingChainsGroup1_tuple, all_overlapped_chain_alignments_list):
    """
        FUNCTION to that finds if a contig of overlapping chain alignments is a subgroup of another bigger contig.
    """
    num_of_overlapping_chains = len(overlappingChainsGroup1_tuple)   # the number of overlapping chain alignments in the current group
    # NOW ITERATE OVER ALL OVERLAPPING CHAIN ALIGNMENT GROUPS AND CHECK IF ALL CHAIN ALIGNMENTS OF THE CURRENT GROUP ARE CONTAINED ALSO ELSEWHERE
    for overlappingChainsGroup2_tuple in all_overlapped_chain_alignments_list:
        IS_SUBGROUP = False
        num_of_matches = 0
        for chain_alignment_tuple1 in overlappingChainsGroup1_tuple:
            if chain_alignment_tuple1 in overlappingChainsGroup2_tuple and overlappingChainsGroup1_tuple != overlappingChainsGroup2_tuple:
                num_of_matches += 1
        if num_of_matches == num_of_overlapping_chains:
            IS_SUBGROUP = True
            break

    return IS_SUBGROUP


def save_overlappingChainsGroups(overlapped_chain_alignments_list, tuple_index_list, iterationNo_list):
    """
        FUNCTION to save each contig of overlapped chain alignments into a pickle file to save memory for subsequence calculations.
    """
    #print "DEBUG: Entered function ..."
    all_overlapped_chain_alignments_list = shared.getConst('ALL_OVERLAPPED_CHAIN_ALIGNMENTS_LIST_'+iterationNo_list[0])
    progbar = shared.getConst('PROGBAR')
    overlappingChainsGroup_num = len(overlapped_chain_alignments_list)
    pointer = 0
    for overlappingChainsGroup1_tuple, tuple_index in zip(overlapped_chain_alignments_list, tuple_index_list):
        IS_SUBGROUP = check_if_subgroup(overlappingChainsGroup1_tuple, all_overlapped_chain_alignments_list)
        #with open('tmp_overlapped_chain_alignments_folder/progress.txt', 'a') as f:
        #    f.write(str(tuple_index)+"\n")
        pointer += 1
        progbar.set_progress((float(pointer)/overlappingChainsGroup_num))
        if not IS_SUBGROUP:
            #unique_overlapped_chain_alignments_set.add(overlappingChainsGroup1_tuple)
            # save the contig of overlapping chains into a pickle file to save memory
            #print "DEBUG: saving overlappingChainsGroup number ", tuple_index
            pickle.dump(overlappingChainsGroup1_tuple, open('tmp_overlapped_chain_alignments_folder/overlappingChainsGroup_'+str(tuple_index)+'.pickle', 'wb'))


def create_overlapped_chain_alignments(ordered_consensus_protein_alignment_list,
                                       ordered_consensus_chain_alignment_list,
                                       args):
    """
        FUNCTION to find overlapping chain alignments and form contigs of overlapping chain alignments.
    """

    ## FIND OVERLAPPING CHAIN ALIGNMENTS AND FORM CONTIGS OF OVERLAPPING CHAIN ALIGNMENTS
    refchain_overlappingChainsList_dict = OrderedDict() # ordereddict with keys the reference chain alignment that have at least 1 overlapping chain alignment, and values
                                                        # tuples of the overlapping alignments
    start = 0
    for row_index in range(len(ordered_consensus_chain_alignment_list)):    # iterate over all chain alignments
        #try:
        reference_aligned_chain = ordered_consensus_chain_alignment_list[row_index]
        while reference_aligned_chain[start] == '-':
            start += 1
        #print "DEBUG: start=", start
        #print "DEBUG: reference_aligned_chain=", reference_aligned_chain
        num_of_potential_overlaps = 0
        overlapping_chain_segment_list = []  # list of the AAIG of the reference chain that must match with the overlapping chains
        consensus_sequence_list = [reference_aligned_chain[start]]
        for rTAAIG in reference_aligned_chain[(start+1):]:  # the overlap must start from the second TAAIG of the chain
            if rTAAIG == '-':
                break
            overlapping_chain_segment_list.append(rTAAIG)
            consensus_sequence_list.append(rTAAIG)

        if len(overlapping_chain_segment_list) < args.OVERLAP_CUTOFF:
            continue

        overlappingChains_list = []
        for query_aligned_chain in ordered_consensus_chain_alignment_list[row_index+1:]: # iterate over all remaining aligned chains to find overlaps
            # CRITERIA FOR PASS:
            # 1. if the maximum number of overlapping AAIG has been reached
            # 2. if the 2nd chain does not start from the 1st TAAIG of the 1st chain
            # 3. if the 2nd chain is not a subchain of the 1st chain
            # 4. if the C-term overhanging end of the 2nd chain is not already in the 1st chain
            try:
                if query_aligned_chain[start] == '-' and overlapping_chain_segment_list == query_aligned_chain[start+1:(start+1+len(overlapping_chain_segment_list))] and query_aligned_chain[start+1+len(overlapping_chain_segment_list)] != '-' and not query_aligned_chain[start+1+len(overlapping_chain_segment_list)] in consensus_sequence_list:
                    overlappingChains_list.append(query_aligned_chain)                   # save the overlaping chain
                    #ordered_consensus_chain_alignment_list.remove(query_aligned_chain)  # and remove it from the ordered_consensus_chain_alignment_list
            except IndexError:
                # if the C-term of reference_aligned_chain is not '-' and hence there is not C-term overhanging end at the query_aligned_chain, so the last condition does not apply here
                if len(query_aligned_chain) == (start+1+len(overlapping_chain_segment_list)) and query_aligned_chain[start] == '-' and overlapping_chain_segment_list == query_aligned_chain[start+1:(start+1+len(overlapping_chain_segment_list))]:
                    overlappingChains_list.append(query_aligned_chain)                   # save the overlaping chain
                    #ordered_consensus_chain_alignment_list.remove(query_aligned_chain)  # and remove it from the ordered_consensus_chain_alignment_list
                    continue
                else:
                    print("DEBUG: query_aligned_chain=", query_aligned_chain)
                    print("DEBUG: len(overlapping_chain_segment_list)=", len(overlapping_chain_segment_list))
                    print("DEBUG: overlapping_chain_segment_list=", overlapping_chain_segment_list)
                    #print "DEBUG: ordered_consensus_chain_alignment_list=", ordered_consensus_chain_alignment_list
                    print("DEBUG: consensus_sequence_list=", consensus_sequence_list)
                    sys.exit(1)

        #ordered_consensus_chain_alignment_list.remove(reference_aligned_chain)  # finally remove the reference chain alignment from the ordered_consensus_chain_alignment_list
        #if len(overlappingChains_list) > 1: # and save the list of the contiguous overlapping chains if they are more than 1
        refchain_overlappingChainsList_dict[tuple(reference_aligned_chain)] = tuple(overlappingChains_list)
        #except IndexError:  # it means than all aligned chains were analyzed and removed from ordered_consensus_chain_alignment_list
        #    break


    all_overlapped_chain_alignments_list = []
    for reference_aligned_chain in list(refchain_overlappingChainsList_dict.keys()):
        all_overlapped_chain_alignments_list.extend(build_Alignment_Tree(reference_aligned_chain, refchain_overlappingChainsList_dict, ordered_consensus_protein_alignment_list, ordered_consensus_chain_alignment_list))

    #unique_overlapped_chain_alignments_set = set()  # set with all unique alignments of chains that have at least one other overlapping chain
    if os.path.exists('tmp_overlapped_chain_alignments_folder/'):
        # if the folder exist from a previous run, remove it to create a clean one
        shutil.rmtree('tmp_overlapped_chain_alignments_folder/', ignore_errors=True)
    os.makedirs('tmp_overlapped_chain_alignments_folder/')

    overlappingChainsGroup_num = len(all_overlapped_chain_alignments_list)
    print("DEBUG: len(all_overlapped_chain_alignments_list)=", overlappingChainsGroup_num)
    if overlappingChainsGroup_num == 0:
        print("ATTENTION: NO OVERLAPPING PEPTIDES WERE CREATED! PROBABLY THERE IS SOMETHING WRONG WITH YOUR DATA. LOOK AT THE connectivities_all AND amino_acid_type_prediction_probabilities FILES.")
        sys.exit(1)

    return refchain_overlappingChainsList_dict, all_overlapped_chain_alignments_list


def get_unique_overlapped_chain_alignments_set():
    unique_overlapped_chain_alignments_set = set()  # set with all unique alignments of chains that have at least one other overlapping chain
    # LOAD ALL PREVIOUSLY WRITEN PEPTIDE FILES
    fnames=os.listdir('tmp_overlapped_chain_alignments_folder/')
    fpattern = re.compile('overlappingChainsGroup_[0-9]+.pickle')
    overlappingChainsGroups_list = list(filter(fpattern.search, fnames))
    for pickle_file in overlappingChainsGroups_list:
        overlappingChainsGroup_tuple = pickle.load( open( 'tmp_overlapped_chain_alignments_folder/' + pickle_file, "rb" ) )
        unique_overlapped_chain_alignments_set.add(overlappingChainsGroup_tuple)

    return unique_overlapped_chain_alignments_set



def write_overlapped_chains_alignment(template_sequence_alignment_string, args, iteration=1):
    """
        FUNCTION to write the various contigs of overlapping chain alignments into the files consensus_alignment.overlapped_chains.iteration[1-7].
    """

    unique_overlapped_chain_alignments_set = get_unique_overlapped_chain_alignments_set()
    ## PRINT THE OVERLAPPING CHAIN ALIGNMENTS
    with open(args.consensus_alignment_file+".overlapped_chains.iteration"+str(iteration), 'w') as f:
        total_chars_per_line = len(template_sequence_alignment_string)

        for overlappingChains_tuple in unique_overlapped_chain_alignments_set:
            #print "\n\n"    # change 3 lines to separate successive parts of the protein sequence
            f.write("\n\n\n")
            start = 0
            while (start+args.chars_per_line-1) <= total_chars_per_line:
                #print "\n\n"
                #print template_sequence_alignment_string[start:(start+args.chars_per_line-1)]
                f.write(template_sequence_alignment_string[start:(start+args.chars_per_line-1)] + "\n")
                for chain_alignment_tuple in overlappingChains_tuple:
                    ## RECOVER THE CHAIN ALIGNMENT STRING
                    chain_line = '     '
                    for word in chain_alignment_tuple:
                        if word == '-':
                            chain_line += '      -     '
                        else:
                            chain_line += "%12s" % word
                    #print chain_line[start:(start+args.chars_per_line-1)]
                    f.write(chain_line[start:(start+args.chars_per_line-1)] + "\n")
                #print ""    # change line to separate the groups of overlapping chains)
                f.write("\n")
                start += args.chars_per_line-1
        #for overlappingChains_tuple in unique_overlapped_chain_alignments_set:

            #print template_sequence_alignment_string[start:]
            f.write(template_sequence_alignment_string[start:] + "\n")
            for chain_alignment_tuple in overlappingChains_tuple:
                ## RECOVER THE CHAIN ALIGNMENT STRING
                chain_line = '     '
                for word in chain_alignment_tuple:
                    if word == '-':
                        chain_line += '      -     '
                    else:
                        chain_line += "%12s" % word
                #print chain_line[start:]
                f.write(chain_line[start:] + "\n")
            #print ""    # change line to separate the groups of overlapping chains)
            f.write("\n")


def compare_2_alignments(start, end, aln1, aln2):
    """
        FUNCTION to compare if a seuqence starting from position "start" and ending at position "end" is common between two contigs.
    """
    for i in range(start, end+1):
        if aln1[i] != aln2[i]:
            return False
    return True


def trim_overhanging_ends(consensus_sequence_alignments_set):
    """
        FUNCTION to trim the N- or C-term overhanging ends given the condition that there is at least one other contig which extends at least 2
        positions in the respective direction and that it conflicts with the N- or C-term of the currect contig.
    """

    # sort the contigs acoording to length
    alignmentLength_list = []
    for alignment in consensus_sequence_alignments_set:
        alignmentLength_list.append((list(alignment), len([c for c in alignment if not c in ["-", "N/A"]])))
    seq_length = len(alignmentLength_list[0][0])
    alignmentLength_list.sort(key=itemgetter(1))    # sort from shortest to longest
    contigs_list = [[a for a in aln if a != '-'] for aln,l in alignmentLength_list]
    for contig,duplet in zip(contigs_list, alignmentLength_list):
        aln = duplet[0]
        N_TIG = contig[0]
        N_index = aln.index(N_TIG)
        C_TIG = contig[-1]
        C_index = aln.index(C_TIG)
        if N_index < 2 or C_index >= seq_length-2:
            continue
        for other_aln in consensus_sequence_alignments_set:
            if compare_2_alignments(N_index+1, C_index-1, aln, other_aln) == False: # if they share the same intermediate region
                continue
            if other_aln[N_index-1] != "-" and other_aln[N_index-2] != "-" and other_aln[N_index] != N_TIG: # trim the N-term
                aln[N_index] = "-"
                break
        for other_aln in consensus_sequence_alignments_set:
            if compare_2_alignments(N_index+1, C_index-1, aln, other_aln) == False: # if they share the same intermediate region
                continue
            if other_aln[C_index+1] != "-" and other_aln[C_index+2] != "-" and other_aln[C_index] != C_TIG: # trim the N-term
                aln[C_index] = "-"
                break

    new_consensus_sequence_alignments_list = [a[0] for a in alignmentLength_list]
    return new_consensus_sequence_alignments_list


def write_consensus_alignment_and_absolute_macthes(template_sequence_alignment_string,
                                                   refchain_overlappingChainsList_dict,
                                                   max_peptide_length,
                                                   ordered_consensus_peptide_alignment_list,
                                                   ordered_consensus_chain_alignment_list,
                                                   ordered_consensus_TAAIGScore_alignment_list,
                                                   args,
                                                   position_TAAIG_confidenceScore_mdict=None,
                                                   iteration=1):
    """
        ARGUMENTS:
        template_sequence_alignment_string:         the alignment (without gaps) of the protein sequence in string form
        refchain_overlappingChainsList_dict:        ordereddict with keys the reference chain alignment that have at least 1 overlapping chain alignment, and values
                                                    tuples of the overlapping alignments
        max_peptide_length:                         the maximum length of the input peptides

    """

    unique_overlapped_chain_alignments_set = get_unique_overlapped_chain_alignments_set()   # get a set of tuples, each one containing one cluster of overlapping chain alignments as tuples
    ## GET THE ALIGNMENT OF THE CONSENSUS SEQUENCE FROM EACH CONTIG OF OVERLAPPING CHAIN ALIGNMENTS
    #print "DEBUG: unique_overlapped_chain_alignments_set=", unique_overlapped_chain_alignments_set
    consensus_sequence_alignments_set = set() # set with the alignment of the consensus sequence of every cluster of overlapping chains
    for overlappingChainsGroup_tuple in unique_overlapped_chain_alignments_set:
        alignment_length = len(overlappingChainsGroup_tuple[0])
        consensus_overlappingChains_alignment_list = ['-'] * alignment_length
        for chain_alignment_tuple in overlappingChainsGroup_tuple:
            for index in range(alignment_length):
                if consensus_overlappingChains_alignment_list[index] == '-' and chain_alignment_tuple[index] != '-':
                    consensus_overlappingChains_alignment_list[index] = chain_alignment_tuple[index]
        ## TRIM THE N-TERM OVERHANGING END
        #for forward_index in range(alignment_length):
        #    if consensus_overlappingChains_alignment_list[forward_index] != '-':
        #        consensus_overlappingChains_alignment_list[forward_index] = '-'
        #        break
        ## TRIM THE C-TERM OVERHANGING END
        #for backward_index in reversed(range(alignment_length)):
        #    if consensus_overlappingChains_alignment_list[backward_index] != '-':
        #        consensus_overlappingChains_alignment_list[backward_index] = '-'
        #        break
        ## TAKE CARE OF DUPLICATIONS IN THE ALIGNMENT
        #print "DEBUG: consensus_overlappingChains_alignment_list=", consensus_overlappingChains_alignment_list
        for aln in get_alignments_with_unique_AAIG(consensus_overlappingChains_alignment_list, max_peptide_length):
            consensus_sequence_alignments_set.add(tuple(aln))
    # print "DEBUG: consensus_sequence_alignments_set=", consensus_sequence_alignments_set
    if iteration >= 5 and args.TRIM_ENDS == True:
        consensus_sequence_alignments_list = trim_overhanging_ends(consensus_sequence_alignments_set)
    else:
        consensus_sequence_alignments_list = [list(a) for a in consensus_sequence_alignments_set]

    assert len(consensus_sequence_alignments_list) > 0, Debuginfo("FAIL: consensus sequence alignment is empty! Relax your mratio,"
            " mcutoff, zacutoff, zmcutoff criteria to form more chains/peptides.", fail=True)

    print("DEBUG: args=", args)
    ## NOW PRINT THE CONSENSUS SEQUENCE ALIGNMENT OF EVERY GROUP OF OVERLAPPING CHAINS
    with open(args.consensus_alignment_file+".overlapped_chains_common_sequence.iteration"+str(iteration), 'w') as f:
        total_chars_per_line = len(template_sequence_alignment_string)
        start = 0
        while (start+args.chars_per_line-1) <= total_chars_per_line:
            #print "\n\n"
            #print template_sequence_alignment_string[start:(start+args.chars_per_line-1)]
            f.write(template_sequence_alignment_string[start:(start+args.chars_per_line-1)] + "\n")
            for chain_alignment_list in consensus_sequence_alignments_list:
                ## RECOVER THE CHAIN ALIGNMENT STRING
                chain_line = '     '
                for word in chain_alignment_list:
                    if word == '-':
                        chain_line += '      -     '
                    else:
                        chain_line += "%12s" % word
                #print ''.join(chain_line[start:(start+args.chars_per_line-1)])
                f.write(''.join(chain_line[start:(start+args.chars_per_line-1)]) + "\n")
            start += args.chars_per_line-1
        #print "\n\n"
        f.write("\n\n\n")
        #print template_sequence_alignment_string[start:]
        f.write(template_sequence_alignment_string[start:] + "\n")
        for chain_alignment_list in consensus_sequence_alignments_list:
            ## RECOVER THE CHAIN ALIGNMENT STRING
            chain_line = '     '
            for word in chain_alignment_list:
                if word == '-':
                    chain_line += '      -     '
                else:
                    chain_line += "%12s" % word
            #print ''.join(chain_line[start:])
            f.write(''.join(chain_line[start:]) + "\n")

        ## NOW PRINT BELOW EACH POSTION OF THE PROTEIN SEQUENCE THE AAIG THAT WAS PRESENT IN ALL PREDICTED CONSENSUS
        # SEQUENCES (ABSOLUTE MATCHES)
        #print "DEBUG: consensus_sequence_alignments_list=", consensus_sequence_alignments_list
        absolute_matches_alignment = get_alignment_with_absolute_matches(consensus_sequence_alignments_list, args)
        template_sequence_alignment_string = template_sequence_alignment_string.split()
        if len(absolute_matches_alignment) == len(template_sequence_alignment_string) -1:
            # Debuginfo('ERROR: the protein sequence and the absolute AIIG matches in the alignment do not'
            #           ' have the same length.', fail=True)
            absolute_matches_alignment.insert(0, 'N/A')
        #f.write("\n\n\n")
        #f.write("ABSOLUTE MATCHES:\n")
        f.write("\n\n\nABSOLUTE MATCHES:\n\n\n" + "\n")
        print("\n\n\nABSOLUTE MATCHES:\n\n\n")
        #print "DEBUG: absolute_matches_alignment=", absolute_matches_alignment
        #print "DEBUG: template_sequence_alignment_string=", template_sequence_alignment_string
        for start, end in zip(list(range(0,len(absolute_matches_alignment),17)), list(range(17,len(absolute_matches_alignment),17))):
            headers = template_sequence_alignment_string[start:end]
            table = [ absolute_matches_alignment[start:end] ]
            print("DEBUG: headers=", headers)
            print("DEBUG: table=", table)
            print(tabulate(table, headers, tablefmt="fancy_grid"))
            f.write(tabulate(table, headers, tablefmt="fancy_grid") +
                    "\n")  # convert the \n to the same format (binary), too.
            print("\n\n")
            f.write("\n\n" + "\n")
        headers = template_sequence_alignment_string[end:]
        headers.append('N/A')
        table = [ absolute_matches_alignment[end:] ]  # the last element is a gap
        print(tabulate(table, headers, tablefmt="fancy_grid"))
        f.write(tabulate(table, headers, tablefmt="fancy_grid") +
                "\n")  # convert the \n to the same format (binary), too.
        #f.write(template_sequence_alignment_string + "\n\n")
        #f.write(get_alignment_with_absolute_matches(consensus_sequence_alignments_list))


    ## CALCULATE CONSENSUS SCORE FOR EACH POSITION OF THE ALIGNMENT
    if args.USE_ALL_PEPTIDES_FOR_PREDICTION:    ## use all peptides to score the AAIG prediction
        filtered_ordered_consensus_chain_alignment_list = ordered_consensus_chain_alignment_list
        filtered_ordered_consensus_TAAIGScore_alignment_list = ordered_consensus_TAAIGScore_alignment_list
    else:   ## use only overlapping peptides to score the AAIG prediction
        all_overlapping_chain_alignments_set = set()  # set with all the alignments of all the chains that have at least one ovelapping chain, in tuple form
        for ref_aln in list(refchain_overlappingChainsList_dict.keys()):
            all_overlapping_chain_alignments_set.add(ref_aln)
            for query_aln in refchain_overlappingChainsList_dict[ref_aln]:
                all_overlapping_chain_alignments_set.add(tuple(query_aln))
        ## check for exceptions in the "overlapping chain rule"
        for peptide_aln, chain_aln in zip(ordered_consensus_peptide_alignment_list, ordered_consensus_chain_alignment_list):
            if is_valid_alignment(peptide_aln, chain_aln, max_peptide_length): # add terminal chains or chains interrupted by Prolines
                all_overlapping_chain_alignments_set.add(tuple(chain_aln))

        filtered_ordered_consensus_chain_alignment_list = []    # list of lists of the form: ['-', '-', ..., 'T252', 'S267', 'L268', '-', '-']
        filtered_ordered_consensus_TAAIGScore_alignment_list = []  # list of lists of the form: [('', 0.0), ('', 0.0), ..., ('T252', 0.01318430523215), ('S267', 7.53408792136e-10), ('L268', 0.008264462809928), ('', 0.0), ('', 0.0)]
                                                                    # which correspond to the ordered chain alignments
        for chain_alignment_list, TAAIGScoreTuple_list in zip(ordered_consensus_chain_alignment_list, ordered_consensus_TAAIGScore_alignment_list):
            if tuple(chain_alignment_list) in all_overlapping_chain_alignments_set:
                filtered_ordered_consensus_chain_alignment_list.append(chain_alignment_list)
                filtered_ordered_consensus_TAAIGScore_alignment_list.append(TAAIGScoreTuple_list)

    print("DEBUG: ordered_consensus_chain_alignment_list=", len(ordered_consensus_chain_alignment_list))
    print("DEBUG: filtered_ordered_consensus_chain_alignment_list=", len(filtered_ordered_consensus_chain_alignment_list))
    ## Now sum the confidence scores at each position of the protein sequence
    if position_TAAIG_confidenceScore_mdict == None:
        position_TAAIG_confidenceScore_mdict = tree()
        for position in range(len(filtered_ordered_consensus_chain_alignment_list[0])):
            #print "\nDEBUG: position=",position
            for chain_alignment_list, TAAIGScoreTuple_list in zip(filtered_ordered_consensus_chain_alignment_list, filtered_ordered_consensus_TAAIGScore_alignment_list):
                    TAAIG = TAAIGScoreTuple_list[position][0]
                    if len(TAAIG) == 0:
                        continue
                    Score = TAAIGScoreTuple_list[position][1]
                    #sys.stdout.write(str(TAAIG)+"->"+str(Score)+" ")
                    try:
                        position_TAAIG_confidenceScore_mdict[position][TAAIG] += Score
                    except TypeError:
                        position_TAAIG_confidenceScore_mdict[position][TAAIG] = Score

    ## NOW WRITE THE AAIG PREDICTIONS AND CONSENSUS SCORES FOR EACH POSITION OF THE ALIGNMENT
    max_num_of_predictions = 0
    for position in range(len(filtered_ordered_consensus_chain_alignment_list[0])):
        if max_num_of_predictions < len(list(position_TAAIG_confidenceScore_mdict[position].keys())):
            max_num_of_predictions = len(list(position_TAAIG_confidenceScore_mdict[position].keys()))

    consensus_score_alignments_list = [[""] * len(filtered_ordered_consensus_chain_alignment_list[0])] * max_num_of_predictions
    consensus_score_alignments_TAAIG_array = np.array(consensus_score_alignments_list, dtype=object)   # use dtype=object to be able to change the elements to every string length
    consensus_score_alignments_Score_array = np.zeros( (max_num_of_predictions, len(filtered_ordered_consensus_chain_alignment_list[0])) )
    for position in range(len(filtered_ordered_consensus_chain_alignment_list[0])):
        row_index = 0
        TAAIGConsScoreTuples_list = list(position_TAAIG_confidenceScore_mdict[position].items())
        if len(TAAIGConsScoreTuples_list) > 0:
            TAAIGConsScoreTuples_list.sort(key=itemgetter(1), reverse=True)  # sort TAAIG predictions by consensus score
        for TAAIGConsScore_tuple in TAAIGConsScoreTuples_list:
            TAAIG = TAAIGConsScore_tuple[0]
            ConsScore = TAAIGConsScore_tuple[1]
            consensus_score_alignments_TAAIG_array[row_index][position] = TAAIG
            consensus_score_alignments_Score_array[row_index][position] = ConsScore
            row_index += 1

    with open(args.consensus_alignment_file+".prediction_scores.iteration"+str(iteration), 'w') as f:
        f.write("TOCSY INDICES AND PR"
                "EDICTION SCORES:\n" + "\n")
        print("AAIGs AND PREDICTION SCORES:\n")
        for start, end in zip(list(range(0,len(consensus_score_alignments_list[0]),5)), list(range(5,len(consensus_score_alignments_list[0]),5))):
            headers = template_sequence_alignment_string[start:end]
            table = []
            for row_index in range(len(consensus_score_alignments_list)):
                new_row = [re.sub(r"^\s0\.0", "", str(TAAIG)+" "+str(ConsScore)) for TAAIG, ConsScore in zip(consensus_score_alignments_TAAIG_array[row_index], consensus_score_alignments_Score_array[row_index]) ]
                new_row.insert(0, '')
                table.append(new_row[start:end])
            #print "DEBUG: table=", table
            #print "DEBUG: headers=", headers
            print(tabulate(table, headers, tablefmt="fancy_grid"))
            f.write(tabulate(table, headers, tablefmt="fancy_grid") +
                    "\n")  # convert the \n to the same format (binary), too.
            print("\n\n")
            f.write("\n\n\n")
        headers = template_sequence_alignment_string[end:]
        headers.append('N/A')
        table = []
        for row_index in range(len(consensus_score_alignments_list)):
            new_row = [re.sub(r"^\s0\.0", "", str(TAAIG)+" "+str(ConsScore)) for TAAIG, ConsScore in zip(consensus_score_alignments_TAAIG_array[row_index], consensus_score_alignments_Score_array[row_index]) ]
            new_row.insert(0, '')
            table.append(new_row[end:])   # the last element is a gap
        print(tabulate(table, headers, tablefmt="fancy_grid"))
        f.write(tabulate(table, headers, tablefmt="fancy_grid") +
                "\n")  # convert the \n to the same format (binary), too.

    return absolute_matches_alignment, position_TAAIG_confidenceScore_mdict


def filter_alignments_using_absolute_matches(absolute_matches_alignment, ordered_consensus_peptide_alignment_list, ordered_consensus_chain_alignment_list,
                                            ordered_consensus_TAAIGScore_alignment_list, position_TAAIG_confidenceScore_mdict):
    """
        FUNCTION to filter
        absolute_matches_alignment= [u' N/A ', u'*', u' Q136', u' N137', u'*', u'*', u'*', u'*', u'*', u'*', u'*', u'  S99', u' A100', u' I101', u' W102', u' S103', u'*', u'*',
        u'*', u'*', u'*', u'*', u'*', u'*', u'*', u' I113', u'*', u' S115', u' I116', u'*', u'*', u'*', u'*', u'*', u'*', u'*', u'*', u'*', u' V126', u'*', u'*', u'*', u' Y130',
        u' G131', u' N132', u' R133', u' E134', u'*', u'*', u'*', u'*', u'*', u'*', u' L141', u'*', u'*', u'*', u'*', u'*', u'*', u'*']
    """

    #print "DEBUG: position_TAAIG_confidenceScore_mdict=", position_TAAIG_confidenceScore_mdict
    # make the format of absolute_matches_alignment as that of ordered_consensus_peptide_alignment_list
    #absolute_matches_alignment[0] = '-'
    absolute_matches_alignment = absolute_matches_alignment[1:] # remove the first element: ' N/A '
    for position in range(len(absolute_matches_alignment)):
        if absolute_matches_alignment[position] == '*':
            absolute_matches_alignment[position] = '-'
        else:
            absolute_matches_alignment[position] = absolute_matches_alignment[position].replace(" ", "")
    #print "DEBUG: processed absolute_matches_alignment=", absolute_matches_alignment

    ## KEEP ONLY THE ABSOLUTE MATCHES THAT HAVE THE HIGHEST CONFIDENCE SCORE AT THEIR POSITION IN THE ALIGNMENT
    for position in range(len(absolute_matches_alignment)):
        if absolute_matches_alignment[position] == '-':
            continue
        # otherwise check if there is an agreement at that position between the absolute match and the confidence scores methods
        TAAIGConfscoreTuple_list = list(position_TAAIG_confidenceScore_mdict[position].items())
        if len(TAAIGConfscoreTuple_list) == 0: # if no confidence score is available for this position don't change the absolute_matches_alignment
            if absolute_matches_alignment[position] != '-': # if there is a TAAIG at this position (obviously from the applied restraints)
                rst_TAAIG = absolute_matches_alignment[position]
                position_TAAIG_confidenceScore_mdict[position][rst_TAAIG] = 10000 # instead of conf score place 10000 the distinguish that it is a restraint
                continue
            else:
                continue
        TAAIGConfscoreTuple_list.sort(key=itemgetter(1), reverse=True)
        try:
            if TAAIGConfscoreTuple_list[0][0] != absolute_matches_alignment[position]:
                absolute_matches_alignment[position] = '-'  # if there is a discrepancy, remove the absolute match !
            else:   # if there is an agreement, clean the alternative choices in position_TAAIG_confidenceScore_mdict
                for TAAIGConfscore_tuple in TAAIGConfscoreTuple_list[1:]:
                    TAAIG = TAAIGConfscore_tuple[0]
                    del position_TAAIG_confidenceScore_mdict[position][TAAIG]
        except IndexError:
            print("DEBUG: position=", position)
            print("DEBUG: TAAIGConfscoreTuple_list=", TAAIGConfscoreTuple_list)
            sys.exit(1)

    #print "DEBUG: absolute_matches_alignment=", absolute_matches_alignment)
    for peptide_line, chain_alignment_list, TAAIGScoreTuple_list in zip(ordered_consensus_peptide_alignment_list, ordered_consensus_chain_alignment_list, ordered_consensus_TAAIGScore_alignment_list):
        for position in range(len(chain_alignment_list)):
            try:
                if chain_alignment_list[position] == '-':
                    continue
                elif absolute_matches_alignment[position] != '-' and chain_alignment_list[position] == absolute_matches_alignment[position]:    # if there must be a different TAAIG at that position, discard it
                    continue
                elif (not chain_alignment_list[position] in absolute_matches_alignment) and absolute_matches_alignment[position] == '-':
                    continue
                elif chain_alignment_list[position] in absolute_matches_alignment and absolute_matches_alignment[position] == chain_alignment_list[position]:
                    continue
                else:
                    #try:
                    #print "DEBUG: deleting chain_alignment_list=", chain_alignment_list
                    ordered_consensus_peptide_alignment_list.remove(peptide_line)
                    ordered_consensus_chain_alignment_list.remove(chain_alignment_list)
                    ordered_consensus_TAAIGScore_alignment_list.remove(TAAIGScoreTuple_list)
                    break
            except IndexError:
                print("DEBUG: position=", position)
                print("DEBUG: absolute_matches_alignment=", absolute_matches_alignment)
                print("DEBUG: peptide_line=", peptide_line)
                print("DEBUG: chain_alignment_list=", chain_alignment_list)
                sys.exit(1)
        #print "DEBUG: keeping chain_alignment_list=", chain_alignment_list
    #print "DEBUG: absolute_matches_alignment=", absolute_matches_alignment

    ## REMOVE AAIG that are absolute matches but also top scored in their position in the alignment, from other positions in the alignment that are irrelevant.
    # E.g. we found S115 to be an consensus match at position 115. Remove it from all other positions, like position 101!

    # print the contents of position_TAAIG_confidenceScore_mdict after filtering
    #for pos in range(len(position_TAAIG_confidenceScore_mdict)):
    #    print "DEBUG: before pos",pos,"->",position_TAAIG_confidenceScore_mdict[pos].keys()

    for absmatch_position in range(len(absolute_matches_alignment)):
        TAAIG = absolute_matches_alignment[absmatch_position]  # these are the absolute matches that agree with cofnidence scores!
        if TAAIG != '-':
            for confscore_position in range(len(position_TAAIG_confidenceScore_mdict)):
                if confscore_position != absmatch_position:
                    #if len(position_TAAIG_confidenceScore_mdict[confscore_position]) == 1:
                    #    print "ERROR: "
                    # delete the consensus confscore/absolute match from every other position in position_TAAIG_confidenceScore_mdict
                    try:
                        del position_TAAIG_confidenceScore_mdict[confscore_position][TAAIG]
                    except KeyError:
                        continue

    # print the contents of position_TAAIG_confidenceScore_mdict after filtering
    #for pos in range(len(position_TAAIG_confidenceScore_mdict)):
    #    print "DEBUG: after pos",pos,"->",position_TAAIG_confidenceScore_mdict[pos].keys()

    return ordered_consensus_peptide_alignment_list, \
           ordered_consensus_chain_alignment_list, \
           ordered_consensus_TAAIGScore_alignment_list,\
           position_TAAIG_confidenceScore_mdict, \
           absolute_matches_alignment


def print_results_summary(template_sequence_alignment_string,
                          consensus_absolute_matches_alignment,
                          position_TAAIG_confidenceScore_mdict,
                          ordered_consensus_chain_alignment_list,
                          i_iminus1_dict=None,
                          SPECTRUM_COMBO="HCNH-HCNH",
                          out_fname="results_summary",
                          silent=False):
    """

    :param template_sequence_alignment_string: string of the protein sequence in the form: "   A    N    I    V    G    G    I ...  L    V    T    G  "
    :param consensus_absolute_matches_alignment: list of the form ['X12NXHX', 'X15NXHX', 'X16NXHX', 'X17NXHX', 'X18NXHX', '-', '-', ...]
    :param position_TAAIG_confidenceScore_mdict: otherwise position (int) -> AAIG signature -> confidence score (float)
    :param ordered_consensus_chain_alignment_list:  list of the chains aligned to the sequence. E.g.
    [['-', '-', 'F26NH', 'V27NH', 'A28NH', 'G29NH', ...], [...]]
    :param i_iminus1_dict:
    :param SPECTRUM_COMBO: controls the shift of the aligned peptides. E.g. in TOCSY-HCNH we had a
                            shift of -1, but in NOESY we have shift 0.
    :return:
    """

    ## NOW WRITE THE AAIG PREDICTIONS AND CONSENSUS SCORES FOR EACH POSITION OF THE ALIGNMENT
    max_num_of_predictions = 0
    for position in range(len(ordered_consensus_chain_alignment_list[0])):
        if max_num_of_predictions < len(list(position_TAAIG_confidenceScore_mdict[position].keys())):
            max_num_of_predictions = len(list(position_TAAIG_confidenceScore_mdict[position].keys()))

    consensus_score_alignments_list = [[""] * len(ordered_consensus_chain_alignment_list[0])] * max_num_of_predictions
    consensus_score_alignments_TAAIG_array = np.array(consensus_score_alignments_list,
                                                      dtype=object)  # use dtype=object to be able to change the elements to every string length
    consensus_score_alignments_Score_array = np.zeros(
        (max_num_of_predictions, len(ordered_consensus_chain_alignment_list[0])))
    for position in range(len(ordered_consensus_chain_alignment_list[0])):
        row_index = 0
        TAAIGConsScoreTuples_list = list(position_TAAIG_confidenceScore_mdict[position].items())
        if len(TAAIGConsScoreTuples_list) > 0:
            TAAIGConsScoreTuples_list.sort(key=itemgetter(1), reverse=True)  # sort TAAIG predictions by consensus score
        for TAAIGConsScore_tuple in TAAIGConsScoreTuples_list:
            TAAIG = TAAIGConsScore_tuple[0]
            ConsScore = TAAIGConsScore_tuple[1]
            consensus_score_alignments_TAAIG_array[row_index][position] = TAAIG
            consensus_score_alignments_Score_array[row_index][position] = ConsScore
            row_index += 1

    with open(out_fname, 'w') as f:
        template_sequence_alignment_string = template_sequence_alignment_string.split()
        f.write("AAIGs AND PREDICTION SCORES:\n" + "\n")
        if not silent:  print("AAIGs AND PREDICTION SCORES:\n")
        for start, end in zip(list(range(0, len(consensus_score_alignments_list[0]), 5)),
                              list(range(5, len(consensus_score_alignments_list[0]), 5))):
            headers = template_sequence_alignment_string[start:end]
            table = []
            for row_index in range(len(consensus_score_alignments_list)):
                new_row = [re.sub(r"^\s0\.0", "", str(TAAIG) + " " + str(ConsScore)) for TAAIG, ConsScore in
                           zip(consensus_score_alignments_TAAIG_array[row_index],
                               consensus_score_alignments_Score_array[row_index])]
                if SPECTRUM_COMBO == "TOCSY-HCNH":
                    new_row.insert(0, '')
                table.append(new_row[start:end])
            # print "DEBUG: table=", table
            # print "DEBUG: headers=", headers
            if not silent:  print(tabulate(table, headers, tablefmt="fancy_grid"))
            f.write(tabulate(table, headers, tablefmt="fancy_grid") +
                    "\n")  # convert the \n to the same format (binary), too.
            if not silent:  print("\n\n")
            f.write("\n\n\n")
        headers = template_sequence_alignment_string[end:]
        headers.append('N/A')
        table = []
        for row_index in range(len(consensus_score_alignments_list)):
            new_row = [re.sub(r"^\s0\.0", "", str(TAAIG) + " " + str(ConsScore)) for TAAIG, ConsScore in
                       zip(consensus_score_alignments_TAAIG_array[row_index],
                           consensus_score_alignments_Score_array[row_index])]
            if SPECTRUM_COMBO == "TOCSY-HCNH":
                new_row.insert(0, '')
            table.append(new_row[end:])  # the last element matches with N/A
        if not silent:  print(tabulate(table, headers, tablefmt="fancy_grid"))
        f.write(tabulate(table, headers, tablefmt="fancy_grid") +
                "\n")  # convert the \n to the same format (binary), too.

        ## NOW PRINT BELOW EACH POSITION OF THE PROTEIN SEQUENCE THE AAIG THAT WAS PRESENT IN ALL
        ## PREDICTED CONSENSUS SEQUENCES (ABSOLUTE MATCHES) & WAS RANKED 1st AT THAT POSITION BY
        ## THE CONFIDENCE SCORE
        connectivities_alignment = []
        if i_iminus1_dict != None:
            for i in range(len(consensus_absolute_matches_alignment) - 1):
                TIG_i = consensus_absolute_matches_alignment[i]
                TIG_iplus1 = consensus_absolute_matches_alignment[i + 1]
                if TIG_i != '-' and TIG_iplus1 != '-':
                    try:
                        connectivity = [c for c in i_iminus1_dict[TIG_iplus1] if c[0] == TIG_i][
                            0]  # connectivity = (TIG_i, occupancy, TOCSY_resonnum, [intersection], [multi])
                        print("DEBUG: connectivity=", connectivity)
                        connectivities_alignment.append(str(connectivity[1]) + "/" + str(connectivity[2]))
                    except (IndexError, KeyError):
                        connectivities_alignment.append('-')
                else:
                    connectivities_alignment.append('-')
            connectivities_alignment.append('-')

        connectivities_alignment.insert(0, '-')  # TODO: should I include it in the if statement???
        if SPECTRUM_COMBO == "TOCSY-HCNH":
            consensus_absolute_matches_alignment.insert(0, 'N/A')
            connectivities_alignment.insert(0, 'N/A')
        # consensus_absolute_matches_alignment = consensus_absolute_matches_alignment[:-1]    # fix alignment with protein sequence
        # f.write("\n\n\n")
        # f.write("ABSOLUTE MATCHES:\n")
        f.write("\n\n\nCONSENSUS CONFIDENCE TOP-SCORED AND ABSOLUTE MATCHES:\n\n\n" + "\n")
        if not silent:  print("\n\n\nCONSENSUS CONFIDENCE TOP-SCORED AND ABSOLUTE MATCHES:\n\n\n")
        # print("DEBUG: consensus_absolute_matches_alignment=", consensus_absolute_matches_alignment)
        # print "DEBUG: template_sequence_alignment_string=", template_sequence_alignment_string
        full_headers, full_NHmappings = [], []
        for start, end in zip(list(range(0, len(consensus_absolute_matches_alignment), 17)),
                              list(range(17, len(consensus_absolute_matches_alignment), 17))):
            # print "DEBUG: table=", table
            # print "DEBUG: headers=", headers
            headers = template_sequence_alignment_string[start:end]
            full_headers.extend(headers)    # for proof-reading
            table = [consensus_absolute_matches_alignment[start:end]]
            full_NHmappings.extend(consensus_absolute_matches_alignment[start:end]) # for proof-reading
            if i_iminus1_dict != None:
                table.append(connectivities_alignment[start:end])
            if not silent:  print(tabulate(table, headers, tablefmt="fancy_grid"))
            f.write(tabulate(table, headers, tablefmt="fancy_grid") +
                    "\n")  # convert the \n to the same format (binary), too.
            if not silent:  print("\n\n")
            f.write("\n\n" + "\n")
        headers = template_sequence_alignment_string[end:]
        full_headers.extend(headers)  # for proof-reading
        headers.append('N/A')
        table = [consensus_absolute_matches_alignment[end:]]  # the last element is a gap
        full_NHmappings.extend(consensus_absolute_matches_alignment[end:])  # for proof-reading
        if i_iminus1_dict != None:
            table.append(connectivities_alignment[end:])
        # print "DEBUG: table=", table
        # print "DEBUG: headers=", headersdf = pd
        if not silent:  print(tabulate(table, headers, tablefmt="fancy_grid"))
        f.write(tabulate(table, headers, tablefmt="fancy_grid") +
                "\n")  # convert the \n to the same format (binary), too.
        # f.write(template_sequence_alignment_string + "\n\n")
        # f.write(get_alignment_with_absolute_matches(consensus_sequence_alignments_set))

        # WRITE THE STATISTICS OF NH-MAPPING AT THE END OF THE FILE
        df = pd.DataFrame([full_headers, list(map(remove_NH_suffix, full_NHmappings))]) \
            .T.rename(columns={0: 'header', 1: 'AAIGsign'})
        df = df.join(df['AAIGsign'].str.extract('^([A-Z])([0-9]+)$') \
                     .rename({0: 'resname', 1: 'resid'}, axis=1) \
                     .fillna(9999999) \
                     .astype({'resname': str, 'resid': int}, errors='ignore'))
        df['resid'] = df['resid'] - df['resid'].min()
        df.index = df.index.values - df[df['resid'] == 0].index[0]
        CORRECT = df[(df['header']==df['resname']) & (df.index==df['resid'])].shape[0]
        WRONG = df[(~df['AAIGsign'].isin(['N/A', '-'])) & ((df['header']!=df['resname']) | (df.index!=df['resid']))].shape[0]
        wrong_AIIGsign_list = df[(~df['AAIGsign'].isin(['N/A', '-'])) & ((df['header']!=df['resname']) | (df.index!=df['resid']))].values.tolist()
        UNASSIGNED = df[df['AAIGsign'].isin(['N/A', '-'])].shape[0]
        f.write("\n\n\nSTATISTICS:\n\n")
        f.write("CORRECT = " + str(CORRECT) + "\n")
        f.write("WRONG = " + str(WRONG) + " (" + ", ".join(wrong_AIIGsign_list) + ")\n")
        f.write("UNASSIGNED = %i\n" % UNASSIGNED)
        if not silent:
            print("CORRECT = " + str(CORRECT))
            print("WRONG = " + str(WRONG) + " (" + ", ".join(wrong_AIIGsign_list) + ")")

def print_graph_results_summary(graph_output, fasta, out_fname, first_resid=1, silent=False):
    """
    Method to write the NH-mapping generated by the Graph algorithm in 4D-CHAINS format.

template_sequence_alignment_string,\
consensus_absolute_matches_alignment,\
final_position_TAAIG_confidenceScore_mdict,\
final_ordered_consensus_chain_alignment_list,\
i_iminus1_dict,\
ALIGNMENT_SHIFT = load_pickle("print_results_summary_input.pkl")


template_sequence_alignment_string: remains the same
consensus_absolute_matches_alignment: form it from graph_output and ALIGNMENT_SHIFT?
final_position_TAAIG_confidenceScore_mdict: creat dummy
final_ordered_consensus_chain_alignment_list: think!
i_iminus1_dict: ignore
    :return:
    """
    # # Example input
    # graph_output = "K82NH M2NH K3NH L4NH V5NH K6NH X21NH R8NH K9NH G10NH DUMMY58 N51ND2HD22 N51ND2HD21 " \
    #                "X86NH L15NH R16NH L17NH A18NH G19NH DUMMY40 N21ND2HD21 D22NH V23NH G85NH DUMMY5 F26NH " \
    #                "V27NH A28NH G29NH V30NH L31NH E32NH D33NH S34NH PRO35 A36NH A37NH K38NH E39NH G40NH " \
    #                "L41NH E42NH E43NH G44NH D45NH Q46NH I47NH L48NH R49NH DUMMY59 N89ND2HD21 DUMMY53 N51NH " \
    #                "V50NH F55NH T56NH N57NH I58NH I59NH R60NH E61NH E62NH A63NH V64NH L65NH F66NH L67NH " \
    #                "L68NH D69NH L70NH PRO107 K72NH G73NH E74NH E75NH V76NH T77NH I78NH L79NH A80NH Q81NH " \
    #                "DUMMY55 Q46NE2HE22 Q46NE2HE21 I25NH G24NH G20NH G84NH N89ND2HD22 L90NH Y91NH I92NH " \
    #                "Q93NH W94NH L95NH K96NH D97NH G98NH G99NH PRO71 V53NH D54NH G103NH R104NH PRO100 PRO106 " \
    #                "PRO105 S108NH"
    # Create template_sequence_alignment_string
    fasta = FASTA(fasta, full_seq=True)
    protresseq_list = ["%s%i" % (aa, first_resid+i) for i, aa in enumerate(fasta.protseq_list)]
    template_sequence_alignment_string = ('  ' + '    '.join(protresseq_list) + '  ')

    # Create consensus_absolute_matches_alignment
    consensus_absolute_matches_alignment = []
    for line in open(graph_output, 'r'):
        for word in line.split():
            if word.startswith("DUMMY"):
                word = "-"
            consensus_absolute_matches_alignment.append(word)

    # Create dummy final_position_TAAIG_confidenceScore_mdict
    position_TAAIG_confidenceScore_mdict = tree()
    for pos, AAIGsign in enumerate(consensus_absolute_matches_alignment):
        if AAIGsign == "-":
            continue
        position_TAAIG_confidenceScore_mdict[pos][AAIGsign] = 1.0     # default confidence score

    # Create dummy ordered_consensus_chain_alignment_list
    ordered_consensus_chain_alignment_list = [consensus_absolute_matches_alignment] # Graph does not produce aligned chains

    # Print and write the NH-mapping in 4D-CHAINS format
    print_results_summary(template_sequence_alignment_string=template_sequence_alignment_string,
                          consensus_absolute_matches_alignment=consensus_absolute_matches_alignment,
                          position_TAAIG_confidenceScore_mdict=position_TAAIG_confidenceScore_mdict,
                          ordered_consensus_chain_alignment_list=ordered_consensus_chain_alignment_list,
                          i_iminus1_dict=None,
                          SPECTRUM_COMBO='HCNH-HCNH',
                          out_fname=out_fname,
                          silent=silent)

def is_graph_NHmap(NHmap_file):
    for line in open(NHmap_file):
        if '│' in line:
            return False
    return True
