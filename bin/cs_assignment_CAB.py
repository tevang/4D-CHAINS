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


#!/usr/bin/env python2.7

import sys, re, os
from operator import itemgetter
from itertools import combinations, permutations
from ordereddict import OrderedDict
from argparse import ArgumentParser, SUPPRESS
from tabulate import tabulate
from ete3 import Tree
import gc, math
import collections, copy
import numpy as np
from cluster import HierarchicalClustering
def tree(): # function to create multidimensional dictionaries
    return collections.defaultdict(tree)
CHAINS_BIN_DIR = os.path.dirname(os.path.realpath(__file__))
CHAINS_LIB_DIR = CHAINS_BIN_DIR[:-3] + "lib"
sys.path.append(CHAINS_BIN_DIR)
sys.path.append(CHAINS_LIB_DIR)
from global_vars import *
from spectrum_processing import *
from chain_formation import *
from aa_prediction import *
from peptides import *
from global_func import *
from cs import *



def cmdlineparse():
    parser = ArgumentParser(description="command line arguments",
        epilog="EXAMPLE: cs_assignment.py -absfile results_summary.3mers_round1_rst.chainlinkers.patch -rstart 200 -tocsy tocsyHCNH.21.4.2016num.list -probprod -2dhist -noesy noesyHCNH.21.4.2016num.list")
    parser.add_argument("-absfile", dest="ABSOLUTE_MATCHES_FILE", required=True, help="file with absolute matches from previous run", metavar="<absolute matches file>")
    parser.add_argument("-rstart", dest="FIRST_RESIDUE_NUMBER", required=True, default=1, type=int, 
                        help="the number of the first residue in the protein (default: 1)", metavar="<first residue number>")
    parser.add_argument("-tocsy", dest="TOCSY_fname", required=True, help="4D TOCSY (HCTOCSYNH) file", metavar="<4D TOCSY input file>")
    parser.add_argument("-noesy", dest="NOESY_fname", required=False, help="optionally the 4D NOESY file in case the -patch argument was used when running chain_linker.py", metavar="<4D NOESY input file>")
    parser.add_argument("-wH", dest="H_weight", required=False, type=float, default=0.1, help="weight in (0.0-1.0] to apply on aa type prediction from aliphatic H resonances (default: 0.1)", metavar="<H weight>")
    parser.add_argument("-wC", dest="C_weight", required=False, type=float, default=1.0, help="weight in (0.0, 1.0] to apply on aa type prediction from aliphatic C resonances (default: 1.0)", metavar="<C weight>")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s version 0.9')
    parser.add_argument("-probprod", dest="PROB_PRODUCT", required=False, action='store_true', default=False,
                        help="select the best C-H type assignment combination based on the product of probabilities of the individual C-H assignments")
    parser.add_argument("-probmode", dest="CONSENSUS_CGROUP_PROB_MODE", required=False, default=4, type=int,
                        help="""The way to calculate the total score of a set of chemical shift assignment (default: 4). Can be:
                    0: just multiply all the probabilities of the individual peaks
                    The following values control how to calculate the consensus probability of each C-group. The total score will be the
                    product of this consensus C-group probabilities.
                    1: average;
                    2: sqrt(prob1*prob2)    ; geometric mean
                    3: (prob2-prob1)/(log(prob2)-log(prob1))    ; logarithmic mean
                    4: (prob1*prob2)/((prob1+prob2)/2.) ;composite average 
                    5: (log(prob1)+log(prob2))/2
                        """, metavar="<way to calculate the cs assignment score>")
    parser.add_argument("-2dhist", dest="USE_2D_HISTOGRAMS", required=False, action='store_true', default=False,
                        help="use 2D BMRB histograms for aa type prediction")
    args=parser.parse_args()
    return args


args = cmdlineparse()

args.ABSOLUTE_MATCHES_FILE = os.path.abspath(args.ABSOLUTE_MATCHES_FILE)
args.TOCSY_fname = os.path.abspath(args.TOCSY_fname)
args.NOESY_fname = os.path.abspath(args.NOESY_fname)



# It does not necessarily contain real residue names!
protein_alignment_list, absolute_RIGmatches_alignment_list  = read_NHmap(args.ABSOLUTE_MATCHES_FILE)


residue2Tindex_dict = OrderedDict()     # OrderedDict with keys: residue (resname{3 letters}+resid) -> assigned TOCSY index group (TIG)
residue2position_dict = OrderedDict()   # OrderedDict with keys: residue (resname+resid) -> position in protein sequence (starting from 0)
resid = args.FIRST_RESIDUE_NUMBER
for position in range(1, len(protein_alignment_list)):
    resname = protein_alignment_list[position-1]    # recall that the TOCSY resonances correspond to residue i-1
    Tindex = absolute_RIGmatches_alignment_list[position]
    if resname == 'N/A':
        continue
    residue = aa1to3_dict[resname] + str(resid)
    residue2position_dict[residue] = position
    if Tindex == '-' or Tindex == 'N/A':
        residue2Tindex_dict[residue] = None
    else:
        residue2Tindex_dict[residue] = Tindex 
    resid += 1

print "DEBUG: residue2Tindex_dict=", residue2Tindex_dict
def read_spectrum_file(query_fname):
    
    with open(query_fname, 'r') as f:
        tmp_query_contents=f.readlines()    # contents of original query_fname (4D TOCSY or 4D NOESY) in 5 column format (name H C N HN)
    query_contents=[]
    for line in tmp_query_contents: # CLEANING THE SPECTRUM FROM IRRELEVANT LINES
        word_list = line.split()
        try:
            float(word_list[1]); float(word_list[2]); float(word_list[3]); float(word_list[4]); # checking if the line contains numbers (resonances)
            if word_list[0][0:7] != "?-?-?-?":
                print "DEBUG: appending word_list:", word_list
                query_contents.append(" ".join(word_list[:5])+"\n")
        except (IndexError, ValueError):
            print "WARNING: Discarding TOCSY line:", line
    
    lines2remove_set = set()
    for qline1 in query_contents:
        try:
            counter = 0
            query1_words_list = qline1.split()
            q1_w1=float(query1_words_list[1])
            q1_w2=float(query1_words_list[2])
            q1_w3=float(query1_words_list[3])
            q1_w4=float(query1_words_list[4])
            for qline2 in query_contents:
                query2_words_list = qline2.split()
                q2_w1=float(query2_words_list[1])
                q2_w2=float(query2_words_list[2])
                q2_w3=float(query2_words_list[3])
                q2_w4=float(query2_words_list[4])
                if approx_equal_proportional(q1_w1, q2_w1) and approx_equal_proportional(q1_w2, q2_w2) and approx_equal_proportional(q1_w3, q2_w3) and approx_equal_proportional(q1_w4, q2_w4):
                    counter += 1
                    if counter > 1:
                        lines2remove_set.add(qline2)
        except (ValueError, IndexError):
            continue
    for qline in lines2remove_set:
        query_contents.remove(qline)
    
    query_lineLists_list = []   # list of the lines of query frame in list form not in string
    for qline in query_contents:
        qline_list = qline.split()
        query_lineLists_list.append(qline_list)
    sorted_query_lineLists_list = sorted(query_lineLists_list, key=itemgetter(0))   # sort the spectrum lines by the assigned TOCSY index
    Tindex_CSlist_dict = {}
    previous_resid = None
    for qline_list in sorted_query_lineLists_list:
        current_resid = qline_list[0].replace('?-?-', '').replace('N-H', '')
        if current_resid != previous_resid:
            Tindex_CSlist_dict[current_resid] = []
        Tindex_CSlist_dict[current_resid].append((qline_list[1:]))
        previous_resid = current_resid
    
    return Tindex_CSlist_dict

Tindex_CSlist_dict = read_spectrum_file(args.TOCSY_fname)
residue_CSlist_dict = OrderedDict() # ordereddict with keys the resname+resid -> the list of the associated H-C-N-HN resonances. E.g.
patched_residues_list = []  # list with the residues that were added by -patch option in the chain_linker.py (they have no TOCSY peaks but have NOESY)
for residue in residue2Tindex_dict.keys():
    print "DEBUG: residue=", residue 
    if residue2Tindex_dict[residue] == None:
        print "DEBUG: residue2Tindex_dict[residue]=", residue2Tindex_dict[residue]
        continue
    try:
        residue_CSlist_dict[residue] = Tindex_CSlist_dict[residue2Tindex_dict[residue]]
    except KeyError:    # in case this residue from the alignment has no TOCSY peaks (was added using -patch in the chain_linker.py), save it
        patched_residues_list.append(residue)
        continue

print "DEBUG: residue_CSlist_dict =", residue_CSlist_dict
print "DEBUG: patched_residues_list=", patched_residues_list
NOESY_residue_CSlist_dict = OrderedDict() # ordereddict with keys the resname+reside -> the list of the associated H-C-N-HN resonances. E.g.
if args.NOESY_fname:
    NOESY_Tindex_CSlist_dict = read_spectrum_file(args.NOESY_fname)  # same as Tindex_CSlist_dict, but for the NOESY; it includes the right peaks as well as peaks of nearby residues
    for residue in residue2Tindex_dict.keys():
        if residue2Tindex_dict[residue] == None:
            continue
        elif residue in patched_residues_list:
            try:
                NOESY_residue_CSlist_dict[residue] = NOESY_Tindex_CSlist_dict[residue2Tindex_dict[residue]]
            except KeyError:    # in case this residue from the alignment has no TOCSY and no NOESY peaks (was added using -patch in the chain_linker.py), print error!
                print "ERROR: residue ", residue, " was patched but has to NOESY peaks!!!"
                sys.exit(1)

print "DEBUG: NOESY_residue_CSlist_dict=", NOESY_residue_CSlist_dict



print "\nLoading BMRB chemical shift histograms..."
aa_carbon_binDensityList_multidict = tree()   # multidict with amino acid name -> Carbon atom type -> [array of bin limits, array of probability density]
aa_hydrogen_binDensityList_multidict = tree()   # multidict with amino acid name -> Hydrogen atom type -> [array of bin limits, array of probability density]
fnames = os.listdir(CHAINS_BIN_DIR+"/../databases/histograms/")
fpattern = re.compile("[A-Z]{3}_[A-Z0-9]+_hist.txt$")
hist_files_list = filter(fpattern.search, fnames)
for hist_file in hist_files_list:
    aa = hist_file.split("_")[0]
    atom = hist_file.split("_")[1]
    if atom in allowed_aa_atoms_dict[aa]:   # load histograms of the allowed atom types only
        bin_list = []
        density_list = []
        with open(CHAINS_BIN_DIR+"/../databases/histograms/"+hist_file, 'r') as f:
            for line in f:
                word_list = line.split()
                bin_list.append(float(word_list[0]))
                density_list.append(float(word_list[1]))
        bin_array = np.array(bin_list)
        density_array = np.array(density_list)/sum(density_list)
        if atom[0] == "C":
            aa_carbon_binDensityList_multidict[aa][atom] = [bin_array, density_array]
        elif atom[0] == "H":
            aa_hydrogen_binDensityList_multidict[aa][atom] = [bin_array, density_array]

if args.USE_2D_HISTOGRAMS == True:
    print "\nLoading 2D BMRB chemical shift histograms..."
    aa_CHpair_binProbabilityList_multidict = tree()   # multidict with amino acid name -> Carbon atom type -> [array of bin limits, array of probability density]
    fnames = os.listdir(CHAINS_BIN_DIR+"/../databases/histograms/")
    fpattern = re.compile("[A-Z]{3}_[A-Z0-9-]+_correlated_2Dhist.smoothed.txt$")
    hist_files_list = filter(fpattern.search, fnames)
    for hist_file in hist_files_list:
        aa = hist_file.split("_")[0]
        CH_pair = hist_file.split("_")[1]
        if CH_pair in aa_CHpair_2Dhist_multidict[aa].keys():   # load histograms of the allowed atom types only
            x_bin_list = []     # for H
            y_bin_list = []     # for C
            probability_list = []
            with open(CHAINS_BIN_DIR+"/../databases/histograms/"+hist_file, 'r') as f:
                for line in f:
                    word_list = line.split()
                    x_bin_list.append(float(word_list[0]))
                    y_bin_list.append(float(word_list[1]))
                    probability_list.append(float(word_list[2]))
            x_bin_array = np.array(x_bin_list)
            y_bin_array = np.array(y_bin_list)
            probability_array = np.array(probability_list)
            aa_CHpair_binProbabilityList_multidict[aa][CH_pair] = [x_bin_array, y_bin_array, probability_array]


    



residue_residueTypesProbTupleList_dict = OrderedDict() # ordereddict with keys the aa index of residue i and values lists of tuples of the form (residue i-1 matching aa type, average probability)
Num_of_TOCSY_resonances = 0 # save here the number of TOCSY resonances for a particular aa index
residue = None
possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list = [] # list of lists of the form (aa type, Carbon name, Hydrogen name, weighted average probability, TOCSY_reson_index, H resonance, C resonance, carbon group)
previous_aaindex = ""
atom_index = 1
xeasy_fout = open(args.ABSOLUTE_MATCHES_FILE+".xeasy", 'w')
sparky_fout = open(args.ABSOLUTE_MATCHES_FILE+".sparky", 'w')
residue_list = residue_CSlist_dict.keys()
carbon_groups_list = []  # list of lists of the form [H resonance, C resonance, N resonance, HN resonance, carbon group]
for res_index, residue in enumerate(residue_list):
    print "Assigning chemical shifts to residue ", residue
    aa_type = residue[0:3]  # aa type in 3-letter code
    TOCSYindex_ResonancesTuple_dict = {}
    for TOCSY_H_resonance,TOCSY_C_resonance,TOCSY_N_resonance,TOCSY_HN_resonance in residue_CSlist_dict[residue]:
        Num_of_TOCSY_resonances += 1
        if args.USE_2D_HISTOGRAMS == True:
            valid_matches_list = get_all_assignments_from_H_C_resonpair_2Dhist(aa_type, float(TOCSY_H_resonance), float(TOCSY_C_resonance),
                                                Num_of_TOCSY_resonances, aa_CHpair_binProbabilityList_multidict, aa_carbon_binDensityList_multidict,
                                                useonlyCarbon=False)
        else:
            valid_matches_list = get_aatypes_from_H_C_resonpair(aa_type, float(TOCSY_H_resonance), float(TOCSY_C_resonance), Num_of_TOCSY_resonances)
        TOCSYindex_ResonancesTuple_dict[Num_of_TOCSY_resonances] = (TOCSY_H_resonance, TOCSY_C_resonance, TOCSY_N_resonance, TOCSY_HN_resonance)
        possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list.extend(valid_matches_list)
        print "point 0"
        carbon_groups_list.append([TOCSY_H_resonance, TOCSY_C_resonance, TOCSY_N_resonance, TOCSY_HN_resonance, None])
        print "point 1"
        
    print "point 2"
    new_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list = group_carbons(possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list, iteration=None)
    print "point 3"
    aatype_CHnucleiType_presenceProbSum_multidict, aatype_CHnucleiType_TresonIndex_multidict = get_presenceprob_from_all_H_C_resonpairs(new_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list,
                                                                                                                    Num_of_TOCSY_resonances, args)
    if aatype_CHnucleiType_presenceProbSum_multidict == {} and aatype_CHnucleiType_TresonIndex_multidict == {}: # if no atom assignment could be made, skip this aa! (TEMPORARY FIX)
        continue
    print "point 4"
    possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list = [] # list of tuples of the form (aa type, probability, H_resonance, C_resonance, TOCSY_reson_index) containing all matching sets of H,C resonances
    
    print "point 5"
    Num_of_TOCSY_resonances = 0

    
    xeasy_lines_list, sparky_lines_list = [], []    # lines to save the xeasy and sparky lines for modification
    HN_resonances_list, N_resonances_list = [], []
    Creson_Cname_lines_list = [] # list with the names of the printed carbons to avoid double printing (e.g. CB, CG, CD, CE)
    print "DEBUG: aatype_CHnucleiType_TresonIndex_multidict[aa_type].keys() =", aatype_CHnucleiType_TresonIndex_multidict[aa_type].keys()
    for CHnucleiType, TresonIndex  in aatype_CHnucleiType_TresonIndex_multidict[aa_type].items():
        ResonancesTuple = TOCSYindex_ResonancesTuple_dict[TresonIndex]
        Carbon_name = CHnucleiType.split('_')[0]
        Hydrogen_name = CHnucleiType.split('_')[1]
        if len(aa_carbonBondedHydrogensDict_dict[aa_type][Carbon_name]) > 1:    # if this carbon has geminal protons
            print "DEBUG: aa_type=", aa_type, "Carbon_name=", Carbon_name
            print "DEBUG: aatype_CHnucleiType_TresonIndex_multidict[aa_type].keys()=", aatype_CHnucleiType_TresonIndex_multidict[aa_type].keys()
            print "DEBUG: aa_carbonBondedHydrogensDict_dict[aa_type][Carbon_name]=", aa_carbonBondedHydrogensDict_dict[aa_type][Carbon_name]
            if [Carbon_name+"_" in x for x in aatype_CHnucleiType_TresonIndex_multidict[aa_type].keys()].count(True) < len(aa_carbonBondedHydrogensDict_dict[aa_type][Carbon_name]):
                print "DEBUG entered condition 1!"
                if aa_type in aatype_carbon_methylHydrogens_multidict.keys() and Carbon_name in aatype_carbon_methylHydrogens_multidict[aa_type].keys():
                    print "DEBUG entered condition 2!"
                    new_Hydrogen_name = aatype_carbon_methylHydrogens_multidict[aa_type][Carbon_name]
                    Hydrogen_name = new_Hydrogen_name
                #    print "DEBUG entered condition 3!"
        if aa_type in aatype_carbon_methylHydrogens_multidict.keys() and Carbon_name in aatype_carbon_methylHydrogens_multidict[aa_type].keys() and (not aa_type in aa_CarbonListwithDegenerateH_dict.keys() or Carbon_name in aa_CarbonListwithDegenerateH_dict[aa_type]):
            print "DEBUG entered condition 4! aa_type = ", aa_type, "Carbon_name=", Carbon_name
            new_Hydrogen_name = aatype_carbon_methylHydrogens_multidict[aa_type][Carbon_name]
            Hydrogen_name = new_Hydrogen_name
        Carbon_resonance = ResonancesTuple[1]
        Hydrogen_resonance = ResonancesTuple[0]
        N_resonance = ResonancesTuple[2]
        HN_resonance = ResonancesTuple[3]
        HN_resonances_list.append(float(ResonancesTuple[3]))
        N_resonances_list.append(float(ResonancesTuple[2]))
        
        Creson_Cname_lines_list.append([float(Carbon_resonance), Carbon_name])  
            
        xeasy_lines_list.append( [atom_index, float(Hydrogen_resonance), 0.02, Hydrogen_name, int(residue[3:]), aa_type] )  # append the Hydrogen resonance line
        atom_index += 1
        try:
            sparky_lines_list.append( [aa3to1_dict[aa_type] + residue[3:] + Hydrogen_name, Carbon_name, protein_alignment_list[residue2position_dict[residue]] + str(int(residue[3:])+1) + "N",float(Hydrogen_resonance), float(Carbon_resonance), float(N_resonance), float(HN_resonance)] )
        except IndexError:
            sparky_lines_list.append( [aa3to1_dict[aa_type] + residue[3:] + Hydrogen_name, Carbon_name, str(int(residue[3:])+1) + "N", float(Hydrogen_resonance), float(Carbon_resonance), float(N_resonance), float(HN_resonance)] )
        
    
    Cnames_list = [cl[1] for cl in Creson_Cname_lines_list]
    Cnames_set = set(Cnames_list)
    for Carbon_name in Cnames_set:
        if Cnames_list.count(Carbon_name) > 1:  # if this carbon occurs in multiple C-H pairs, calculate its average resonance for the XEASY format
            Carbon_resonance = np.average([cl[0] for cl in Creson_Cname_lines_list if cl[1]==Carbon_name])
        else:
            Carbon_resonance = np.average([cl[0] for cl in Creson_Cname_lines_list if cl[1]==Carbon_name])
        xeasy_lines_list.append( [atom_index, float(Carbon_resonance), 0.2, Carbon_name, int(residue[3:]), aa_type] )
        atom_index += 1
    
    average_HN_resonance = np.mean(HN_resonances_list)
    stdev_HN_resonance = np.std(HN_resonances_list)
    average_N_resonance = np.mean(N_resonances_list)
    stdev_N_resonance = np.std(N_resonances_list)
    
    new_xeasy_lines_list, new_sparky_lines_list = rename_protons(xeasy_lines_list, sparky_lines_list)
    for xeasy_line in new_xeasy_lines_list:
        xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s\n" % (xeasy_line[0], xeasy_line[1], xeasy_line[2], xeasy_line[3], xeasy_line[4], xeasy_line[5]) )
    for sparky_line in new_sparky_lines_list:
        sparky_fout.write("\t%s-%s-%s-H\t\t%.3f\t%.3f\t%.3f\t%.3f\t\n" % (sparky_line[0], sparky_line[1], sparky_line[2], sparky_line[3], sparky_line[4], sparky_line[5], sparky_line[6]))
    
    
    try:
        xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s\n" % (atom_index, float(average_HN_resonance), 0.02, "H", int(residue[3:])+1, aa1to3_dict[protein_alignment_list[residue2position_dict[residue]]])) # The N-HN of i sees the C-H of i-1
        atom_index += 1
        xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s\n" % (atom_index, average_N_resonance, stdev_N_resonance, "N", int(residue[3:])+1, aa1to3_dict[protein_alignment_list[residue2position_dict[residue]]]))   # The N-HN of i sees the C-H of i-1
        atom_index += 1
    except IndexError:
        xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s\n" % (atom_index, float(average_HN_resonance), 0.02, "H", int(residue[3:])+1, aa1to3_dict[protein_alignment_list[residue2position_dict[residue]]])) # The N-HN of i sees the C-H of i-1
        atom_index += 1
        xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s\n" % (atom_index, average_N_resonance, stdev_N_resonance, "N", int(residue[3:])+1, aa1to3_dict[protein_alignment_list[residue2position_dict[residue]]]))   # The N-HN of i sees the C-H of i-1
        atom_index += 1
    
    


residue_residueTypesProbTupleList_dict = OrderedDict() # ordereddict with keys the aa index of residue i and values lists of tuples of the form (residue i-1 matching aa type, average probability)
Num_of_NOESY_resonances = 0 # save here the number of TOCSY resonances for a particular aa index
residue = None
possible_aatype_prob_C_H_resonpair_NOESYindex_list_list = [] # list of lists of the form (aa type, Carbon name, Hydrogen name, weighted average probability, NOESY_reson_index, H resonance, C resonance, carbon group)
previous_aaindex = ""
residue_list = NOESY_residue_CSlist_dict.keys()
carbon_groups_list = []  # list of lists of the form [H resonance, C resonance, N resonance, HN resonance, carbon group]
for res_index, residue in enumerate(residue_list):
        print "Assigning chemical shifts to residue ", residue
        aa_type = residue[0:3]  # aa type in 3-letter code
        NOESYindex_ResonancesTuple_dict = {}
        for NOESY_H_resonance,NOESY_C_resonance,NOESY_N_resonance,NOESY_HN_resonance in NOESY_residue_CSlist_dict[residue]:
            Num_of_NOESY_resonances += 1
            NOESYindex_ResonancesTuple_dict[Num_of_NOESY_resonances] = (NOESY_H_resonance, NOESY_C_resonance, NOESY_N_resonance, NOESY_HN_resonance)
        Num_of_NOESY_resonances = 0
    
        
        HN_resonances_list, N_resonances_list = [], []
        for TresonIndex  in NOESYindex_ResonancesTuple_dict.keys():
            ResonancesTuple = NOESYindex_ResonancesTuple_dict[TresonIndex]
            N_resonance = ResonancesTuple[2]
            HN_resonance = ResonancesTuple[3]
            HN_resonances_list.append(float(ResonancesTuple[3]))
            N_resonances_list.append(float(ResonancesTuple[2]))
        
        average_HN_resonance = np.mean(HN_resonances_list)
        stdev_HN_resonance = np.std(HN_resonances_list)
        average_N_resonance = np.mean(N_resonances_list)
        stdev_N_resonance = np.std(N_resonances_list)
        
        try:
            xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s\n" % (atom_index, float(average_HN_resonance), 0.02, "H", int(residue[3:])+1, aa1to3_dict[protein_alignment_list[residue2position_dict[residue]]])) # The N-HN of i sees the C-H of i-1
            atom_index += 1
            print "\t"+str(atom_index)+"\t"+str(float(average_HN_resonance))+"\t0.02\tH\t"+str(int(residue[3:])+1)
            print "\t"+str(atom_index)+"\t"+str(float(average_N_resonance))+"\t0.2\tN\t"+str(int(residue[3:])+1)
            xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s\n" % (atom_index, average_N_resonance, stdev_N_resonance, "N", int(residue[3:])+1, aa1to3_dict[protein_alignment_list[residue2position_dict[residue]]]))   # The N-HN of i sees the C-H of i-1
            atom_index += 1
        except IndexError:
            xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s\n" % (atom_index, float(average_HN_resonance), 0.02, "H", int(residue[3:])+1, aa1to3_dict[protein_alignment_list[residue2position_dict[residue]]])) # The N-HN of i sees the C-H of i-1
            atom_index += 1
            print "\t"+str(atom_index)+"\t"+str(float(average_HN_resonance))+"\t0.02\tH\t"+str(int(residue[3:])+1)
            print "\t"+str(atom_index)+"\t"+str(float(average_N_resonance))+"\t0.2\tN\t"+str(int(residue[3:])+1)
            xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s\n" % (atom_index, average_N_resonance, stdev_N_resonance, "N", int(residue[3:])+1, aa1to3_dict[protein_alignment_list[residue2position_dict[residue]]]))   # The N-HN of i sees the C-H of i-1
            atom_index += 1

xeasy_fout.close()
sparky_fout.close()


def read_whole_spectrum_file(query_fname):
    """ Function to load all TOCSY file contents, including lines ?-?-?-?. """
    with open(query_fname, 'r') as f:
        tmp_query_contents=f.readlines()    # contents of original query_fname (4D TOCSY or 4D NOESY) in 5 column format (name H C N HN)
    query_contents=[]   # list with the line_lists of TOCSY file
    for line in tmp_query_contents: # CLEANING THE SPECTRUM FROM IRRELEVANT LINES
        word_list = line.split()
        try:
            float(word_list[1]); float(word_list[2]); float(word_list[3]); float(word_list[4]); # checking if the line contains numbers (resonances)
            query_contents.append(line)
        except (IndexError, ValueError):
            print "WARNING: Discarding TOCSY line:", line
    
    lines2remove_set = set()
    for qline1 in query_contents:
        try:
            counter = 0
            query1_words_list = qline1.split()
            q1_w1=float(query1_words_list[1])
            q1_w2=float(query1_words_list[2])
            q1_w3=float(query1_words_list[3])
            q1_w4=float(query1_words_list[4])
            for qline2 in query_contents:
                query2_words_list = qline2.split()
                q2_w1=float(query2_words_list[1])
                q2_w2=float(query2_words_list[2])
                q2_w3=float(query2_words_list[3])
                q2_w4=float(query2_words_list[4])
                if approx_equal_proportional(q1_w1, q2_w1) and approx_equal_proportional(q1_w2, q2_w2) and approx_equal_proportional(q1_w3, q2_w3) and approx_equal_proportional(q1_w4, q2_w4):
                    counter += 1
                    if counter > 1:
                        lines2remove_set.add(qline2)
        except (ValueError, IndexError):
            continue
    for qline in lines2remove_set:
        query_contents.remove(qline)
    
    return query_contents


TOCSY_line_list = read_whole_spectrum_file(args.TOCSY_fname)    # read the original TOCSY, including the "?-?-?-?" lines
sparky_fout = open(args.ABSOLUTE_MATCHES_FILE+".sparky", 'a+')   # append to the annotated Sparky the unlabeled lines
sparky_lines = sparky_fout.readlines()
for TOCSY_line in TOCSY_line_list:
    TOCSY_lineList = TOCSY_line.split()
    print "DEBUG: TOCSY_lineList=", TOCSY_lineList
    FOUND = False
    TOCSY_label = TOCSY_lineList[0]
    TOCSY_Hreson = float(TOCSY_lineList[1])
    TOCSY_Creson = float(TOCSY_lineList[2])
    TOCSY_Nreson = float(TOCSY_lineList[3])
    TOCSY_HNreson = float(TOCSY_lineList[4])
    for line in sparky_lines:
        print "DEBUG: line=", line
        word_list = line.split()
        Hreson = float(word_list[1])
        Creson = float(word_list[2])
        Nreson = float(word_list[3])
        HNreson = float(word_list[4])
        if Hreson == TOCSY_Hreson and Creson == TOCSY_Creson and Nreson == TOCSY_Nreson and HNreson == TOCSY_HNreson:
            print "DEBUG: matching lines:"
            print line,
            print TOCSY_line
            FOUND = True
            break
    
    if FOUND == False:  # if this line hasn't been write into the annotated TOCSY Sparky file, append it
        sparky_fout.write("\t%s\t\t%.3f\t%.3f\t%.3f\t%.3f\t\n" % (TOCSY_label, TOCSY_Hreson, TOCSY_Creson, TOCSY_Nreson, TOCSY_HNreson))


actual_patched_residues_list = []   # patched_residues_list has the i-1 residues, actual_patched_residues_list has the i. I.e. N212 and A213,
for residue in patched_residues_list:
    resid = int(residue[3:])    # recall that we made residue variable by concatenating the aa type (3-letter) and the resid
    actual_patched_residues_list.append(protein_alignment_list[resid+1 -args.FIRST_RESIDUE_NUMBER] + str(resid+1))
sparky_fout.write("# PATCHED RESIDUES: " + " ".join(actual_patched_residues_list) + "\n")
sparky_fout.close()

