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
        epilog="EXAMPLE: cs_assignment_NOESY.py -absfile results_summary.3mers_round1_rst.chainlinkers.patch -rstart 200 -atocsy results_summary.3mers_round1_rst.chainlinkers.patch.sparky -noesy noesyHCNH.21.4.2016num.list -probprod -2dhist -wcthres 100 -bcthres 25 -percentile 0.9")
    parser.add_argument("-absfile", dest="ABSOLUTE_MATCHES_FILE", required=False, help="file with absolute matches from previous run", metavar="<absolute matches file>")
    parser.add_argument("-rstart", dest="FIRST_RESIDUE_NUMBER", required=False, default=1, type=int, 
                        help="the number of the first residue in the protein (default: 1)", metavar="<first residue number>")
    parser.add_argument("-noesy", dest="NOESY_fname", required=True, help="the 4D NOESY file (*num.list) produced by 4D_assignment_parallel.py", metavar="<4D NOESY input file>")
    parser.add_argument("-wH", dest="H_weight", required=False, type=float, default=0.1, help="weight in (0.0-1.0] to apply on aa type prediction from aliphatic H resonances (default: 0.1)", metavar="<H weight>")
    parser.add_argument("-wC", dest="C_weight", required=False, type=float, default=1.0, help="weight in (0.0, 1.0] to apply on aa type prediction from aliphatic C resonances (default: 1.0)", metavar="<C weight>")
    parser.add_argument("-o", dest="OUT_fname", required=True, default=None,
                        help="output file name (assigned 4D TOCSY in sparky format)", metavar="<output file>")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s version 0.9')
    parser.add_argument("-probprod", dest="PROB_PRODUCT", required=False, action='store_true', default=False,
                        help="select the best C-H type assignment combination based on the product of probabilities of the individual C-H assignments")
    parser.add_argument("-2dhist", dest="USE_2D_HISTOGRAMS", required=False, action='store_true', default=False,
                        help="use 2D BMRB histograms for aa type prediction")
    parser.add_argument("-wcthres1", dest="WITHIN_CGROUP_THRESHOLD_ITER1", required=False, default=0.0, type=float,
                        help="The first Carbon type must have probability greater this ratio from the second one of the SAME CARBON GROUP in order to be kept in the final NOESY assignments")
    parser.add_argument("-bcthres1", dest="BETWEEN_CGROUP_THRESHOLD_ITER1", required=False, default=10.0, type=float,
                        help="The same Carbon type of a Carbon group must have probability greater this ratio from the probability of the SAME CARBON TYPE in any other Carbon group in order to be kept in the final NOESY assignments")
    parser.add_argument("-wcthres2", dest="WITHIN_CGROUP_THRESHOLD_ITER2", required=False, default=0.0, type=float,
                        help="same as -wcthres1, but for the 2nd iteration")
    parser.add_argument("-bcthres2", dest="BETWEEN_CGROUP_THRESHOLD_ITER2", required=False, default=10.0, type=float,
                        help="same as -bcthres1, but for the 2nd iteration")
    parser.add_argument("-wcthres3", dest="WITHIN_CGROUP_THRESHOLD_ITER3", required=False, default=0.0, type=float,
                        help="same as -wcthres1, but for the 3rd iteration")
    parser.add_argument("-bcthres3", dest="BETWEEN_CGROUP_THRESHOLD_ITER3", required=False, default=0.0, type=float,
                        help="same as -bcthres1, but for the 3rd iteration")
    parser.add_argument("-wcthres4", dest="WITHIN_CGROUP_THRESHOLD_ITER4", required=False, default=0.0, type=float,
                        help="same as -wcthres1, but for the 4th iteration")
    parser.add_argument("-bcthres4", dest="BETWEEN_CGROUP_THRESHOLD_ITER4", required=False, default=0.0, type=float,
                        help="same as -bcthres1, but for the 4th iteration")
    parser.add_argument("-wcthres5", dest="WITHIN_CGROUP_THRESHOLD_ITER5", required=False, default=0.0, type=float,
                        help="same as -wcthres1, but for the 5th iteration")
    parser.add_argument("-bcthres5", dest="BETWEEN_CGROUP_THRESHOLD_ITER5", required=False, default=0.0, type=float,
                        help="same as -bcthres1, but for the 5th iteration")
    parser.add_argument("-percentile1", dest="PERCENTILE_ITER1", required=False, default=0.9, type=float,
                        help="This arguments is used as a theshold to discard low probability C-H type predictions in iteration 1. Can be 0.9, 0.85, 0.8, 0.75, or 0.0 (no filtering). A 0.9 percentile of value X means that the top 10 percent of the probability density values were above X.")
    parser.add_argument("-percentile2", dest="PERCENTILE_ITER2", required=False, default=0.8, type=float,
                        help="This arguments is used as a theshold to discard low probability C-H type predictions in iteration 2. Can be 0.9, 0.85, 0.8, 0.75, or 0.0 (no filtering). A 0.9 percentile of value X means that the top 10 percent of the probability density values were above X.")
    parser.add_argument("-percentile3", dest="PERCENTILE_ITER3", required=False, default=0.8, type=float,
                        help="This arguments is used as a theshold to discard low probability C-H type predictions in iteration 3. Can be 0.9, 0.85, 0.8, 0.75, or 0.0 (no filtering). A 0.9 percentile of value X means that the top 10 percent of the probability density values were above X.")
    parser.add_argument("-percentile4", dest="PERCENTILE_ITER4", required=False, default=0.8, type=float,
                        help="This arguments is used as a theshold to discard low probability C-H type predictions in iteration 4. Can be 0.9, 0.85, 0.8, 0.75, or 0.0 (no filtering). A 0.9 percentile of value X means that the top 10 percent of the probability density values were above X.")
    parser.add_argument("-percentile5", dest="PERCENTILE_ITER5", required=False, default=0.8, type=float,
                        help="This arguments is used as a theshold to discard low probability C-H type predictions in iteration 5. Can be 0.9, 0.85, 0.8, 0.75, or 0.0 (no filtering). A 0.9 percentile of value X means that the top 10 percent of the probability density values were above X.")
    parser.add_argument("-usertocsy", dest="user_TOCSY_fname", required=False, help="4D TOCSY (HCTOCSYNH) file in Sparky format with atom assignments made by the user", metavar="<Sparky 4D TOCSY input file with user-made atom assignments>")
    parser.add_argument("-usernoesy", dest="user_NOESY_fname", required=False, help="4D NOESY (HCTOCSYNH) file in Sparky format with atom assignments made by the user", metavar="<Sparky 4D NOESY input file with user-made atom assignments>")
    parser.add_argument("-int1", dest="USE_INTENSITIES_ITERATION1", required=False, action='store_true', default=False,
                        help="use peak intensities in iteration 1")
    parser.add_argument("-int2", dest="USE_INTENSITIES_ITERATION2", required=False, action='store_true', default=False,
                        help="use peak intensities in iteration 2")
    parser.add_argument("-int3", dest="USE_INTENSITIES_ITERATION3", required=False, action='store_true', default=False,
                        help="use peak intensities in iteration 3")
    parser.add_argument("-int4", dest="USE_INTENSITIES_ITERATION4", required=False, action='store_true', default=False,
                        help="use peak intensities in iteration 4")
    parser.add_argument("-int5", dest="USE_INTENSITIES_ITERATION5", required=False, action='store_true', default=False,
                        help="use peak intensities in iteration 5")
    parser.add_argument("-ithres1", dest="INTENSITY_THRESHOLD_ITERATION1", required=False, type=float, default=0.1,
                        help="use relative peak intensity threshold")
    parser.add_argument("-ithres2", dest="INTENSITY_THRESHOLD_ITERATION2", required=False, type=float, default=0.1,
                        help="use relative peak intensity threshold")
    parser.add_argument("-ithres3", dest="INTENSITY_THRESHOLD_ITERATION3", required=False, type=float, default=0.0,
                        help="use relative peak intensity threshold")
    parser.add_argument("-ithres4", dest="INTENSITY_THRESHOLD_ITERATION4", required=False, type=float, default=0.0,
                        help="use relative peak intensity threshold")
    parser.add_argument("-ithres5", dest="INTENSITY_THRESHOLD_ITERATION5", required=False, type=float, default=0.0,
                        help="use relative peak intensity threshold")
    parser.add_argument("-itrans1", dest="INTENSITY_TRANSFORM_TYPE_ITER1", required=False, type=int, default=4,
                        help="transform the intensities in iteration1. Can be 1: do nothing; 2: (default: 1)")
    parser.add_argument("-itrans2", dest="INTENSITY_TRANSFORM_TYPE_ITER2", required=False, type=int, default=1,
                        help="transform the intensities in iteration2. Can be 1: do nothing; 2: (default: 1)")
    parser.add_argument("-itrans3", dest="INTENSITY_TRANSFORM_TYPE_ITER3", required=False, type=int, default=4,
                        help="transform the intensities in iteration3. Can be 1: do nothing; 2: (default: 1)")
    parser.add_argument("-itrans4", dest="INTENSITY_TRANSFORM_TYPE_ITER4", required=False, type=int, default=1,
                        help="transform the intensities in iteration4. Can be 1: do nothing; 2: (default: 1)")
    parser.add_argument("-itrans5", dest="INTENSITY_TRANSFORM_TYPE_ITER5", required=False, type=int, default=1,
                        help="transform the intensities in iteration5. Can be 1: do nothing; 2: (default: 1)")
    parser.add_argument("-probmodel", dest="PROBABILITY_MODEL", required=False, type=int, default=2,
                        help="If '1' the probability of each peak will be given by [wH*1Dhist(H)+wC*1Dhist(C)]/(wH+wC). If '2' then by 1Dhist(H)*1Dhist(C)")
    parser.add_argument("-probmode", dest="CONSENSUS_CGROUP_PROB_MODE", required=False, default=1, type=int,
                        help="""The way to calculate the total score of a set of chemical shift assignment (default: 3).
                    The following values control how to calculate the consensus probability of each C-group. The total score will be the
                    product of this consensus C-group probabilities.
                    1: average;
                    2: sqrt(prob1*prob2)    ; geometric mean
                    3: (prob2-prob1)/(log(prob2)-log(prob1))    ; logarithmic mean
                    4: (prob1*prob2)/((prob1+prob2)/2.) ;composite average 
                    5: prob1 + prob2    ; simple sum
                        """, metavar="<way to calculate the cs assignment score>")
    parser.add_argument("-ithresA", dest="ithresA", type=float, default=0.1, help=SUPPRESS) # intensity threshold for ALA CB-QB methyl for both iteration 1 & 3
    parser.add_argument("-ithresT", dest="ithresT", type=float, default=0.1, help=SUPPRESS) # intensity threshold for THR CG2-QG2 methyl for both iteration 1 & 3
    parser.add_argument("-ithresV", dest="ithresV", type=float, default=0.1, help=SUPPRESS) # intensity threshold for VAL CG1-QG1 & CG2-QG2 methyls for both iteration 1 & 3
    parser.add_argument("-ithresI", dest="ithresI", type=float, default=0.1, help=SUPPRESS) # intensity threshold for ILE CD1-QD1 & CG2-QG2 methyl for both iteration 1 & 3
    parser.add_argument("-ithresL", dest="ithresL", type=float, default=None, help=SUPPRESS) # intensity threshold for LEU CD1-QD1 & CD2-QD2 methyls for both iteration 1 & 3
    parser.add_argument("-pthresA", dest="pthresA", type=float, default=0.8, help=SUPPRESS) # intensity threshold for ALA CB-QB methyl for both iteration 1 & 3
    parser.add_argument("-pthresT", dest="pthresT", type=float, default=0.8, help=SUPPRESS) # intensity threshold for THR CG2-QG2 methyl for both iteration 1 & 3
    parser.add_argument("-pthresV", dest="pthresV", type=float, default=0.8, help=SUPPRESS) # intensity threshold for VAL CG1-QG1 & CG2-QG2 methyls for both iteration 1 & 3
    parser.add_argument("-pthresI", dest="pthresI", type=float, default=0.8, help=SUPPRESS) # intensity threshold for ILE CD1-QD1 & CG2-QG2 methyl for both iteration 1 & 3
    parser.add_argument("-pthresL", dest="pthresL", type=float, default=None, help=SUPPRESS) # intensity threshold for LEU CD1-QD1 & CD2-QD2 methyls for both iteration 1 & 3
    
    parser.add_argument("-ithresMethyl", dest="ithresMethyl", type=float, default=None, help=SUPPRESS) # intensity threshold for all Methyls
    parser.add_argument("-pthresMethyl", dest="pthresMethyl", type=float, default=None, help=SUPPRESS) # intensity threshold for all Methyls
    
    args=parser.parse_args()
    return args

args = cmdlineparse()


if args.ithresMethyl:
    args.ithresA, args.ithresT, args.ithresV, args.ithresI, args.ithresL = args.ithresMethyl, args.ithresMethyl, args.ithresMethyl, args.ithresMethyl, args.ithresMethyl
if args.pthresMethyl:
    args.pthresA, args.pthresT, args.pthresV, args.pthresI, args.pthresL = args.pthresMethyl, args.pthresMethyl, args.pthresMethyl, args.pthresMethyl, args.pthresMethyl
    
args.ABSOLUTE_MATCHES_FILE = os.path.abspath(args.ABSOLUTE_MATCHES_FILE)
args.NOESY_fname = os.path.abspath(args.NOESY_fname)



protein_alignment_list, absolute_RIGmatches_alignment_list  = read_NHmap(args.ABSOLUTE_MATCHES_FILE)
NA_indices = []
if protein_alignment_list[0] == 'N/A' or absolute_RIGmatches_alignment_list[0] == 'N/A':
    NA_indices.append(0)
if protein_alignment_list[-1] == 'N/A' or absolute_RIGmatches_alignment_list[-1] == 'N/A':
    NA_indices.append(len(protein_alignment_list)-1)
protein_alignment_list = [protein_alignment_list[i] for i in range(len(protein_alignment_list)) if not i in NA_indices]
absolute_RIGmatches_alignment_list = [absolute_RIGmatches_alignment_list[i] for i in range(len(absolute_RIGmatches_alignment_list)) if not i in NA_indices]
print "DEBUG: protein_alignment_list=", protein_alignment_list
print "DEBUG: absolute_RIGmatches_alignment_list=", absolute_RIGmatches_alignment_list
sys.exit(1)


resid = args.FIRST_RESIDUE_NUMBER        # start counting from the 1st residue number specified by the user
absolute_matches_alignment_list = ['-']*len(absolute_RIGmatches_alignment_list)     # contains the absolutely matched residue names, namely real residue names not RIGs!
for position in range(len(protein_alignment_list)):
    resname = protein_alignment_list[position]
    RIG = absolute_RIGmatches_alignment_list[position]
    if RIG == '-':
        absolute_matches_alignment_list[position] = '-'
    elif RIG == 'N/A':
        absolute_matches_alignment_list[position] = 'N/A'
    else:
        residue = resname + str(resid)
        absolute_matches_alignment_list[position] = residue
    resid += 1

print "DEBUG: absolute_RIGmatches_alignment_list=", absolute_RIGmatches_alignment_list
print "DEBUG: protein_alignment_list=", protein_alignment_list

print "DEBUG: absolute_matches_alignment_list=", absolute_matches_alignment_list


print "\nLoading BMRB chemical shift histograms..."
aa_carbon_binDensityList_multidict = tree()   # multidict with amino acid name -> Carbon atom type -> [array of bin limits, array of probability density]
aa_hydrogen_binDensityList_multidict = tree()   # multidict with amino acid name -> Hydrogen atom type -> [array of bin limits, array of probability density]
fnames = os.listdir(CHAINS_BIN_DIR+"/../databases/histograms/")
fpattern = re.compile("[A-Z]{3}_[A-Z0-9]+_hist.txt$")
hist_files_list = filter(fpattern.search, fnames)
for hist_file in hist_files_list:
    aa = hist_file.split("_")[0]
    atom = hist_file.split("_")[1]
    atoms2load = allowed_aa_atoms_dict[aa]
    if aa == "PHE":
        atoms2load.extend(["HD1", "CD1"])
    elif aa == "TYR":
        atoms2load.extend(["HD1", "HE1", "CD1", "CE1"])
    atoms2load.extend([])
    if atom in atoms2load:   # load histograms of the allowed atom types only
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

aa_CHpair_binProbabilityList_multidict = tree()   # multidict with amino acid name -> Carbon atom type -> [array of bin limits, array of probability density]
print "\nLoading 2D BMRB chemical shift histograms..."
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




iter1_allowed_atoms_dict = {}   # resname->list of allowed carbons in iteration1
iter1_allowed_atoms_dict["LEU"] = ["CD1", "CD2", "CG"]
iter1_allowed_atoms_dict["VAL"] = ["CG1", "CG2"]
iter1_allowed_atoms_dict["ILE"] = ["CD1", "CG2"]
iter1_allowed_atoms_dict["THR"] = ["CG2"]
iter1_allowed_atoms_dict["ALA"] = ["CB"]

aa2ithres_dict = {}
if type(args.ithresA) == float:    aa2ithres_dict["ALA"] = args.ithresA
if type(args.ithresT) == float:    aa2ithres_dict["THR"] = args.ithresT
if type(args.ithresV) == float:    aa2ithres_dict["VAL"] = args.ithresV
if type(args.ithresI) == float:    aa2ithres_dict["ILE"] = args.ithresI
if type(args.ithresL) == float:    aa2ithres_dict["LEU"] = args.ithresL
aa2pthres_dict = {}
print "DEBUG: args.pthresA=", args.pthresA
if type(args.pthresA) == float:    aa2pthres_dict["ALA"] = args.pthresA
if type(args.pthresT) == float:    aa2pthres_dict["THR"] = args.pthresT
if type(args.pthresV) == float:    aa2pthres_dict["VAL"] = args.pthresV
if type(args.pthresI) == float:    aa2pthres_dict["ILE"] = args.pthresI
if type(args.pthresL) == float:    aa2pthres_dict["LEU"] = args.pthresL


TOCSY_residue_assignments_dict, TOCSY_residue_NHresonances_dict, i_to_iminus1_dict, patched_residues_list, TOCSY_residue_peak_intensity_multidict = {}, {}, {}, {}, {}



NOESY_residue_assignments_dict, NOESY_residue_NHresonances_dict, original_NOESY_peaks_dict, NOESY_residue_peak_intensity_multidict = read_assigned_TOCSY_file(args.NOESY_fname, 'NOESY', absolute_RIGmatches_alignment_list, absolute_matches_alignment_list)
print "DEBUG: NOESY_residue_assignments_dict=", NOESY_residue_assignments_dict
print "DEBUG: NOESY_residue_NHresonances_dict=", NOESY_residue_NHresonances_dict
print "DEBUG: original_NOESY_peaks_dict=", original_NOESY_peaks_dict
print "DEBUG: NOESY_residue_peak_intensity_multidict=",
for residue in NOESY_residue_peak_intensity_multidict.keys():
    for peak in NOESY_residue_peak_intensity_multidict[residue].keys():
        norm_intensity = NOESY_residue_peak_intensity_multidict[residue][peak]
        print residue, peak[1], peak[0], norm_intensity

matched_NOESY_residue_assignments_dict = copy.deepcopy(NOESY_residue_assignments_dict)
residue_unmatched_TOCSY_peaks_dict = {}

print "DEBUG iteration1 point1: matched_NOESY_residue_assignments_dict=", matched_NOESY_residue_assignments_dict
for k,v in matched_NOESY_residue_assignments_dict.items():
    print k, "-->", v
print "DEBUG: residue_unmatched_TOCSY_peaks_dict=", residue_unmatched_TOCSY_peaks_dict

print "DEBUG: point 2 TOCSY_residue_assignments_dict=", TOCSY_residue_assignments_dict


checked_residues_set = set()
Cterm_residue_set = set()


xeasy_fout = open(args.ABSOLUTE_MATCHES_FILE+".NOESY.iter1.xeasy", 'w')
sparky_fout = open(args.ABSOLUTE_MATCHES_FILE+".NOESY.iter1.sparky", 'w')
atom_index = 1
TOCSY_assigned_residue_list = [] # list of the residues that have their peaks already assigned using TOCSY & NOESY or only TOCSY
residues_with_written_NH_list = [] # list of the residues for which the N & HN resonances have been written in the xeasy file
revordered_residue_keys = [(k[0], int(k[1:])) for k in matched_NOESY_residue_assignments_dict.keys() if k in absolute_matches_alignment_list]
revordered_residue_keys.sort(key=itemgetter(1), reverse=True)  # reverse ordered keys, namely they start from the last residue in the protein sequence
matched_i_to_iplus1_peaks_multidict = tree()    # residue -> i peak -> matched i+1 peak
residues_with_written_TOCSY_peaks = []


selected_revordered_residue_keys = []
for aa_type  in ["ALA", "THR", "VAL", "ILE", "LEU"]:
    for i_aa_resid in revordered_residue_keys:
        i_residue = i_aa_resid[0]+str(i_aa_resid[1])
        aa = aa1to3_dict[get_aa_type_from_residue(i_residue)]
        if aa == aa_type:
            selected_revordered_residue_keys.append(i_aa_resid)
print "DEBUG: selected_revordered_residue_keys=", selected_revordered_residue_keys



    

for i_aa_resid in selected_revordered_residue_keys:    # both the keys and the C,H assignments correspond to residue i-1
    i_residue = i_aa_resid[0]+str(i_aa_resid[1])
    aa_type = aa1to3_dict[get_aa_type_from_residue(i_residue)]
    if not aa_type in iter1_allowed_atoms_dict.keys():
        continue
    if aa_type in aa2ithres_dict.keys():
        istart = int(10 * aa2ithres_dict[aa_type])
    else:
        istart = 5
    print "DEBUG: aa_type=", aa_type, "istart=", istart
    for I_THRESHOLD in reversed(range(1, istart+1, 1)):
        I_THRESHOLD *= 0.1
        all_possible_assignments_list = []  # list with all the peaks for this residue (TOCSY-NOESY, NOESY) and the respective possible C-H assignments
        print "ITERATION 1: Matching residue i ", i_residue," TOCSY PEAKS to residue's i NOESY peaks..., I_THRESHOLD=", I_THRESHOLD
        checked_residues_set.add(i_residue)
        all_allowed_carbons_list = [x for x in allowed_aa_atoms_dict[aa_type] if x[0] == "C"]   # all allowed carbons for this aa type
        print "DEBUG: all_allowed_carbons_list=", all_allowed_carbons_list
        tmp_assigned_carbons_list = [x[1] for x in matched_NOESY_residue_assignments_dict[i_residue] if x[0]==i_residue]   # these are the assignments of residue i
        print "DEBUG: tmp_assigned_carbons_list=", tmp_assigned_carbons_list
        assigned_carbons_list = []  # list of carbons which have ALL their protons assigned in TOCSY
        partly_assigned_carbons_list = []   # list of methylene carbons for which only one proton frequency was found in TOCSY
        for Cname in tmp_assigned_carbons_list:
            if Cname in aatype_carbon_nondegenerateHlist_multidict[aa_type].keys() and tmp_assigned_carbons_list.count(Cname) == 2:
                assigned_carbons_list.append(Cname)
            if Cname in aatype_carbon_nondegenerateHlist_multidict[aa_type].keys() and tmp_assigned_carbons_list.count(Cname) == 1:
                partly_assigned_carbons_list.append(Cname)
            elif not Cname in aatype_carbon_nondegenerateHlist_multidict[aa_type].keys():
                assigned_carbons_list.append(Cname)
        print "DEBUG: assigned_carbons_list=", assigned_carbons_list
        print "DEBUG: partly_assigned_carbons_list=", partly_assigned_carbons_list
        tmp_missing_carbons_list = []
        for carbon in all_allowed_carbons_list:     # keep only the carbon types of this aa type that were not assigned to peaks
            if not carbon in assigned_carbons_list:
                tmp_missing_carbons_list.append(carbon)
        print "DEBUG: dict key i_residue", i_residue, "tmp_missing_carbons_list=", tmp_missing_carbons_list
        if match_NOESY_i_iplus1_peaks(i_residue, absolute_matches_alignment_list, matched_NOESY_residue_assignments_dict, matched_i_to_iplus1_peaks_multidict) == False:  # add labels to the i->i+1 matched peaks. If i_residue is a C-term residue or presends a gap, skip it from this iteration
            continue
        else:   # otherwise proceed to assigned but not from the extreme carbons in exclude_from_i_iplus_matching_dict
            missing_carbons_list = [c for c in tmp_missing_carbons_list if c in iter1_allowed_atoms_dict[aa_type]]
            print "DEBUG: dict key i_residue", i_residue, "missing_carbons_list=", missing_carbons_list
            if len(missing_carbons_list) == 0:
                update_NOESY_peaks(i_residue, matched_NOESY_residue_assignments_dict, protein_alignment_list, absolute_matches_alignment_list, args, sparky_lines_list=[])   # remove "i+1" and "i-1" labels
                continue
        unassigned_NOESY_peaks_list = []
        if len(missing_carbons_list) != 0:   # IF THERE ARE MISSING CARBONS IN TOCSY, TRY TO FIND THEM IN NOESY (EXCLUDE THE ASSIGNED OR MATCHED NOESY PEAKS)
            print "Residue ", i_residue, " was not completely assigned in TOCSY. Using NOESY to assigned resonances to the rest Carbon types."
            try:    # DANGEROUS but works for the N-term residue (A200 in aLP)
                unassigned_NOESY_peaks_list = [peak for peak in matched_NOESY_residue_assignments_dict[i_residue] if peak[0]=="i+1" and peak[1]=="?" and peak[3]=="?"]
            except KeyError:
                continue
            print "DEBUG: unassigned_NOESY_peaks_list=", unassigned_NOESY_peaks_list
            for NOESY_peak in unassigned_NOESY_peaks_list:  # make C-H type predictions for each unassigned NOESY peak
                NOESY_Creson = NOESY_peak[2]
                NOESY_Hreson = NOESY_peak[4]
                all_possible_assignments_list.extend(get_probabilities_from_H_C_resonpair_2Dhist(i_residue, NOESY_Hreson, NOESY_Creson, missing_carbons_list,
                                                                partly_assigned_carbons_list, args, aa_CHpair_binProbabilityList_multidict, aa_carbon_binDensityList_multidict,
                                                                NOESY_residue_peak_intensity_multidict, matched_i_to_iplus1_peaks_multidict, aa2pthres_dict,
                                                                aa2ithres_dict, protein_alignment_list, absolute_matches_alignment_list, 
                                                                aa_hydrogen_binDensityList_multidict, iteration=1, I_THRESHOLD=I_THRESHOLD))
            [assignment.append('NOESY') for assignment in all_possible_assignments_list]    # we don't need to created a new list, the original all_possible_assignments_list will be updated
            print "DEBUG: i_residue=", i_residue, "all_possible_assignments_list=", all_possible_assignments_list
            written_nucleiNames_list = []   # list of nuclei that have been written already, to avoid duplication
            if len(all_possible_assignments_list) > 0:  # IF WE HAVE FOUND POSSIBLE NOESY PEAK ASSIGNMENTS OF THE MISSING CARBONS
                written_nucleiNames_list = select_best_NOESY_peak_combination(i_residue, all_possible_assignments_list, partly_assigned_carbons_list, atom_index, xeasy_fout, sparky_fout,
                                       protein_alignment_list, absolute_matches_alignment_list, args, NOESY_residue_NHresonances_dict,
                                       matched_NOESY_residue_assignments_dict, TOCSY_assigned_residue_list, iter1_allowed_atoms_dict,  iteration=1)
        else:   # IF ALL THE CARBONS WERE ASSIGNED IN TOCSY, FIND & SAVE PARTLY ASSIGNED CARBONS IN NOESY, SAVE THE MATCHED NOESY RESONANCES, AND IF NOT MATCHED, SAVE THE TOCSY RESONANCES
            print "Residue ", i_residue, " was completely assigned in TOCSY."
            try:
                aa_type = aa1to3_dict[get_aa_type_from_residue(i_residue)]
            except (IndexError, KeyError):    # skip residue i
                print "EXITING: COULD NOT FIND AA_TYPE OF RESIDUE", i_residue
                sys.exit(1)
            for TOCSY_assignment in TOCSY_residue_assignments_dict[i_residue]:
                spectrum_type = 'TOCSY (unmatched)'
                Cname = TOCSY_assignment[1]
                Creson = TOCSY_assignment[2]    # use the TOCSY resonance by default
                Hname = TOCSY_assignment[3]
                Hreson = TOCSY_assignment[4]
                if i_residue in matched_NOESY_residue_assignments_dict.keys():  # first check if this residue has NOESY peaks
                    matched_NOESY_peak = [peak for peak in matched_NOESY_residue_assignments_dict[i_residue] if peak[0]==i_residue and peak[1]==Cname and peak[3]==Hname]
                    if len(matched_NOESY_peak) == 1:    # if this TOCSY peak has been matched with a NOESY peak, use the NOESY resonances (more accurate)
                        Creson = matched_NOESY_peak[0][2]
                        Hreson = matched_NOESY_peak[0][4]
                        spectrum_type = 'TOCSY-NOESY'
                    elif len(matched_NOESY_peak) > 1:
                        print "ERROR: residue's ", i_residue, "TOCSY peak ", TOCSY_assignment, ", has been matched with multiple NOESY peaks", matched_NOESY_peak
                        sys.exit(1)
                all_possible_assignments_list.append([aa_type, Cname, Hname, 10000000.0, Hreson, Creson, spectrum_type])  # set the probability to 1 because we don't know it
                print "DEBUG: i_residue=", i_residue, "aa_type=", aa_type, "all_possible_assignments_list=", all_possible_assignments_list
            written_nucleiNames_list = []   # list of nuclei that have been written already, to avoid duplication
            if len(all_possible_assignments_list) > 0:  # IF WE HAVE FOUND POSSIBLE NOESY PEAK ASSIGNMENTS OF THE MISSING CARBONS
                written_nucleiNames_list = select_best_NOESY_peak_combination(i_residue, all_possible_assignments_list, partly_assigned_carbons_list, atom_index, xeasy_fout, sparky_fout,
                                       protein_alignment_list, absolute_matches_alignment_list, args, NOESY_residue_NHresonances_dict,
                                       matched_NOESY_residue_assignments_dict, TOCSY_assigned_residue_list, iter1_allowed_atoms_dict, iteration=1)


print "DEBUG iteration1 point2: matched_NOESY_residue_assignments_dict=", matched_NOESY_residue_assignments_dict

for residue in absolute_matches_alignment_list:
    if not residue in checked_residues_set and residue != '-':
        print "WARNING: Residue ", residue, " was not checked!!!"


write_NH_of_residues_without_CH(absolute_matches_alignment_list, NOESY_residue_NHresonances_dict, TOCSY_residue_NHresonances_dict,
                                    residues_with_written_NH_list, atom_index, xeasy_fout, sparky_fout, protein_alignment_list, args, 
                                    Cterm_residue_set, iteration=1)

xeasy_fout.close()
sparky_fout.close()


clean_resid_assignments_dict = clean_xeasy_file(args.ABSOLUTE_MATCHES_FILE+".NOESY.iter1.xeasy")  
for resid, peak_list in clean_resid_assignments_dict.items():
    peak_list.sort(key=itemgetter(3))
    clean_resid_assignments_dict[resid] = peak_list

print "DEBUG: ITERATION 1 point 1 clean_resid_assignments_dict=", clean_resid_assignments_dict
xeasy_out = open(args.ABSOLUTE_MATCHES_FILE+".NOESY.cleaned.iter1.xeasy", 'w')
atom_index = 0  # reinitialize atom index
for resid in clean_resid_assignments_dict.keys():
    for peak in clean_resid_assignments_dict[resid]:
        atom_index += 1
        peak[0] = atom_index
        line = "\t".join([str(p) for p in peak])
        xeasy_out.write(line + "\n")
xeasy_out.close()



print "DEBUG: args.user_TOCSY_fname=", args.user_TOCSY_fname, "args.user_NOESY_fname=", args.user_NOESY_fname
if args.user_TOCSY_fname and args.user_NOESY_fname:
    print "STARTING PROOF-READING NOESY CHEMICAL SHIFT ASSIGNMENTS"
    aa_equivalentCarbons_dict = {}
    aa_equivalentCarbons_dict["LEU"] = [["CD1", "CD2"]]
    aa_equivalentCarbons_dict["VAL"] = [["CG1", "CG2"]]
    
    user_TOCSY_residue_assignments_dict, user_TOCSY_residue_NHresonances_dict, user_i_to_iminus1_dict, user_patched_residues_list, user_TOCSY_residue_peak_intensity_multidict = read_assigned_TOCSY_file(args.user_TOCSY_fname, "TOCSY", absolute_RIGmatches_alignment_list, absolute_matches_alignment_list)
    print "DEBUG: user_TOCSY_residue_assignments_dict=", user_TOCSY_residue_assignments_dict
    print "DEBUG: user_TOCSY_residue_NHresonances_dict=", user_TOCSY_residue_NHresonances_dict
    for k,v in user_TOCSY_residue_assignments_dict.items():
        print k, "-->", v
        
    user_NOESY_residue_assignments_dict, user_NOESY_residue_NHresonances_dict, user_original_NOESY_peaks_dict, user_NOESY_residue_peak_instensity_dict = read_assigned_TOCSY_file(args.user_NOESY_fname, 'NOESY', absolute_RIGmatches_alignment_list, absolute_matches_alignment_list)
    print "DEBUG: user_NOESY_residue_assignments_dict=", user_NOESY_residue_assignments_dict
    print "DEBUG: user_NOESY_residue_NHresonances_dict=", user_NOESY_residue_NHresonances_dict
    print "DEBUG: user_original_NOESY_peaks_dict=", user_original_NOESY_peaks_dict
    for residue in user_NOESY_residue_assignments_dict:
        if residue[0] in ['Y', 'F']:
            for assignment in user_NOESY_residue_assignments_dict[residue]:
                if assignment[1] == "CD":
                    assignment[1] = "CD1"
                elif assignment[1] == "QD":
                    assignment[1] = "HD1"
                elif assignment[1] == "CE":
                    assignment[1] = "CE1"
                elif assignment[1] == "QE":
                    assignment[1] = "HE1"
    for residue in user_TOCSY_residue_assignments_dict:
        if residue[0] in ['Y', 'F']:
            for assignment in user_TOCSY_residue_assignments_dict[residue]:
                if assignment[1] == "CD":
                    assignment[1] = "CD1"
                elif assignment[1] == "QD":
                    assignment[1] = "HD1"
                elif assignment[1] == "CE":
                    assignment[1] = "CE1"
                elif assignment[1] == "QE":
                    assignment[1] = "HE1"
    for k,v in user_NOESY_residue_assignments_dict.items():
        print k, "-->", v
    
    user_NOESY_resid_assignments_dict = {}  # same as user_NOESY_resid_assignments_dict, but with keys resids instead of residues
    for k,v in user_NOESY_residue_assignments_dict.items():
        user_NOESY_resid_assignments_dict[int(k[1:])] = v
    user_TOCSY_resid_assignments_dict = {}  # same as user_TOCSY_resid_assignments_dict, but with keys resids instead of residues
    for k,v in user_TOCSY_residue_assignments_dict.items():
        user_TOCSY_resid_assignments_dict[int(k[1:])] = v
    
    print "DEBUG: user_NOESY_resid_assignments_dict=", user_NOESY_resid_assignments_dict
    print "DEBUG: ITERATION 1 point 2 clean_resid_assignments_dict=", clean_resid_assignments_dict
    
    user_VAL_LEU_NOESY_resid_assignments_dict, user_VAL_LEU_TOCSY_resid_assignments_dict = create_equivalentCarbons_assignments_dict(user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict, protein_alignment_list,
                                              absolute_matches_alignment_list, aa_equivalentCarbons_dict, args)
    print "DEBUG: user_VAL_LEU_NOESY_resid_assignments_dict=", user_VAL_LEU_NOESY_resid_assignments_dict
    print "DEBUG: user_VAL_LEU_TOCSY_resid_assignments_dict=", user_VAL_LEU_TOCSY_resid_assignments_dict
    
    
    proofread_clean_resid_assignments_dict = copy.deepcopy(clean_resid_assignments_dict)   # same as clean_resid_assignments_dict but with comments CORRECT or WRONG at the carbon lines
    for resid in proofread_clean_resid_assignments_dict.keys():
        print "DEBUG ITERATION1: proof reading resid=", resid
        resname = proofread_clean_resid_assignments_dict[resid][0][6]
        if resname in aa_equivalentCarbons_dict.keys(): # If this is a LEU or VAL
            Carbon_pair = aa_equivalentCarbons_dict[resname][0]
            if len([a for a in proofread_clean_resid_assignments_dict[resid] if a[3] in Carbon_pair]) > 0: # AND if the methyl carbons have been assigned
                Proton_pair = ["H"+C[1:] for C in Carbon_pair]   # protons of the equivalent Carbons
                program_assignments_list = [assignment for assignment in proofread_clean_resid_assignments_dict[resid] if assignment[3] in Carbon_pair+Proton_pair]   # keep only the resonances assigned by the program to the equivalent methyl carbons and protons
                proof_read_equivalentMethyls(program_assignments_list, resid, user_VAL_LEU_NOESY_resid_assignments_dict, user_VAL_LEU_TOCSY_resid_assignments_dict, Carbon_pair, Proton_pair)   # it doesn't matter if 1 or 2 of the peaks have been proof-read in iteration 1 & 2
        for assignment in proofread_clean_resid_assignments_dict[resid]:
            print "DEBUG ITERATION1: proof reading assignment=", assignment
            spectrum_type = re.sub(r' iter[1-9]$', '', assignment[7])    # remove the "iter[1-9]" suffix to be able to distinguish the type of assignment (NOESY, NOESY+TOCSY, etc.)
            comment = None
            aa_type = assignment[6]
            nucleus = assignment[3] # the nucleus name
            if aa_type in aa_equivalentCarbons_dict.keys() and nucleus in aa_equivalentCarbons_dict[aa_type][0]:   # if this peak has been assigned by the program to LEU or VAL CD1, CD2 or CG1,CG2 respectively, skip it
                continue
            
            if aa_type in aatype_carbon_nongeminalHname_multidict.keys() and nucleus in aatype_carbon_nongeminalHname_multidict[aa_type].keys():
                C_assignment = assignment # the program-made carbon assignment
                H_assignment = None # default value for the program-made proton assignment
                print "DEBUG: resid=", resid, "proofread_clean_resid_assignments_dict[resid]=", proofread_clean_resid_assignments_dict[resid]
                for a in proofread_clean_resid_assignments_dict[resid]:
                    if ((type(aatype_carbon_nongeminalHname_multidict[aa_type][nucleus]) == list and a[3] in aatype_carbon_nongeminalHname_multidict[aa_type][nucleus]) or
                        (type(aatype_carbon_nongeminalHname_multidict[aa_type][nucleus]) == str and a[3] == aatype_carbon_nongeminalHname_multidict[aa_type][nucleus])):
                        H_assignment = a
                        break
                print "DEBUG: C_assignment=", C_assignment, "H_assignment=", H_assignment
                comment = ""
                if spectrum_type in ['NOESY', 'NOESY (C-term hanging residue)', 'NOESY + TOCSY-NOESY', 'TOCSY-NOESY', 'TOCSY-NOESY + NOESY', 'TOCSY (unmatched) + TOCSY-NOESY', 'TOCSY-NOESY + TOCSY (unmatched)']:    # if this resonance comes only from NOESY
                    user_assignment = get_Carbon_resonance(resid, nucleus, 'NOESY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict,
                         aa_equivalentCarbons_dict, protein_alignment_list, absolute_matches_alignment_list, args, get_assignment=True)
                    if user_assignment == False:    # if no user assignment was made for this Carbon, add no extra comment
                        print "WARNING: no user assignment made in NOESY for resid", resid, "and Carbon", nucleus,". Checking TOCSY..."
                        user_assignment = get_Carbon_resonance(resid, nucleus, 'TOCSY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict,
                         aa_equivalentCarbons_dict, protein_alignment_list, absolute_matches_alignment_list, args, get_assignment=True)
                        if user_assignment == False:
                            print "WARNING: no user assignment made also in TOCSY for resid", resid, "and Carbon", nucleus,". Skipping nucleus."
                            assignment[7] += " <NOT ASSIGNED>"
                            continue
                    user_Creson = user_assignment[2]
                    user_Hreson = user_assignment[4]
                    print "DEBUG: user_Creson=", user_Creson, "user_Hreson=", user_Hreson
                    if approx_equal(C_assignment[1], user_Creson, 0.4) and approx_equal(H_assignment[1], user_Hreson, 0.04): # if the user-made assignment agrees with ours, add to the comment "CORRECT"
                        print "DEBUG: CORRECT C_assignment found:", C_assignment
                        print "DEBUG: CORRECT H_assignment found:", H_assignment
                        comment = " <CORRECT>"
                    else:
                        print "DEBUG: WRONG C_assignment found:", C_assignment
                        print "DEBUG: WRONG H_assignment found:", H_assignment
                        comment = " <WRONG>"
                elif spectrum_type in ['TOCSY (average)', 'TOCSY (unmatched)']:
                    user_assignment = get_Carbon_resonance(resid, nucleus, 'TOCSY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict,
                         aa_equivalentCarbons_dict, protein_alignment_list, absolute_matches_alignment_list, args, get_assignment=True)
                    if user_assignment == False:    # if no user assignment was made for this Carbon, add no extra comment
                        print "WARNING: no user assignment made in TOCSY for resid", resid, "and Carbon", nucleus,". Checking NOESY..."
                        user_assignment = get_Carbon_resonance(resid, nucleus, 'NOESY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict,
                         aa_equivalentCarbons_dict, protein_alignment_list, absolute_matches_alignment_list, args, get_assignment=True)
                        if user_assignment == False:
                            print "WARNING: no user assignment made also in NOESY for resid", resid, "and Carbon", nucleus,". Skipping nucleus."
                            assignment[7] += " <NOT ASSIGNED>"
                            continue
                    user_Creson = user_assignment[2]
                    user_Hreson = user_assignment[4]
                    if approx_equal(C_assignment[1], user_Creson, 0.4) and approx_equal(H_assignment[1], user_Hreson, 0.04): # if the user-made assignment agrees with ours, add to the comment "CORRECT"
                        print "DEBUG: CORRECT C_assignment found:", C_assignment
                        print "DEBUG: CORRECT H_assignment found:", H_assignment
                        comment = " <CORRECT>"
                    else:
                        print "DEBUG: WRONG C_assignment found:", C_assignment
                        print "DEBUG: WRONG H_assignment found:", H_assignment
                        comment = " <WRONG>"
                        
                if comment == "":
                    comment = " <NOT ASSIGNED>"
                C_assignment[7] += comment
                H_assignment[7] += comment
                comment = None  # a trick to add no more comments to "assignment" at the end of this block of code
            
            elif nucleus[0] == 'C':   # only if the assignment line belongs to a Carbon then check if the user has made the same assignment
                comment = ""
                if spectrum_type in ['NOESY', 'NOESY (C-term hanging residue)', 'NOESY + TOCSY-NOESY', 'TOCSY-NOESY', 'TOCSY-NOESY + NOESY', 'TOCSY (unmatched) + TOCSY-NOESY', 'TOCSY-NOESY + TOCSY (unmatched)']:    # if this resonance comes only from NOESY
                    user_assignment = get_Carbon_resonance(resid, nucleus, 'NOESY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict,
                         aa_equivalentCarbons_dict, protein_alignment_list, absolute_matches_alignment_list, args)
                    if user_assignment == False:    # if no user assignment was made for this Carbon, add no extra comment
                        print "WARNING: no user assignment made in NOESY for resid", resid, "and Carbon", nucleus,". Checking TOCSY..."
                        user_assignment = get_Carbon_resonance(resid, nucleus, 'TOCSY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict,
                         aa_equivalentCarbons_dict, protein_alignment_list, absolute_matches_alignment_list, args)
                        if user_assignment == False:
                            print "WARNING: no user assignment made also in TOCSY for resid", resid, "and Carbon", nucleus,". Skipping nucleus."
                            assignment[7] += " <NOT ASSIGNED>"
                            continue
                    if approx_equal(assignment[1], user_assignment, 0.4): # if the user-made assignment agrees with ours, add to the comment "CORRECT"
                        print "DEBUG: CORRECT assignment found:", assignment
                        comment = " <CORRECT>"
                    else:
                        print "DEBUG: WRONG assignment found:", assignment
                        comment = " <WRONG>"
                elif spectrum_type in ['TOCSY (average)', 'TOCSY (unmatched)']:
                    user_assignment = get_Carbon_resonance(resid, nucleus, 'TOCSY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict,
                         aa_equivalentCarbons_dict, protein_alignment_list, absolute_matches_alignment_list, args)
                    if user_assignment == False:    # if no user assignment was made for this Carbon, add no extra comment
                        print "WARNING: no user assignment made in TOCSY for resid", resid, "and Carbon", nucleus,". Checking NOESY..."
                        user_assignment = get_Carbon_resonance(resid, nucleus, 'NOESY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict,
                         aa_equivalentCarbons_dict, protein_alignment_list, absolute_matches_alignment_list, args)
                        if user_assignment == False:
                            print "WARNING: no user assignment made also in NOESY for resid", resid, "and Carbon", nucleus,". Skipping nucleus."
                            assignment[7] += " <NOT ASSIGNED>"
                            continue
                    if approx_equal(assignment[1], user_assignment, 0.4): # if the user-made assignment agrees with ours, add to the comment "CORRECT"
                        print "DEBUG: CORRECT assignment found:", assignment
                        comment = " <CORRECT>"
                    else:
                        print "DEBUG: WRONG assignment found:", assignment
                        comment = " <WRONG>"
                
            if comment == None:     # this is a proton or a methyl Carbon or a Carbon with only one proton
                continue
            elif comment == "":
                comment = " <NOT ASSIGNED>"
            assignment[7] += comment
    
    
    print "DEBUG: proofread_clean_resid_assignments_dict=", proofread_clean_resid_assignments_dict
    
    xeasy_out = open(args.ABSOLUTE_MATCHES_FILE+".NOESY.proofread.iter1.xeasy", 'w')
    atom_index = 0  # reinitialize atom index
    for resid in proofread_clean_resid_assignments_dict.keys():
        for peak in proofread_clean_resid_assignments_dict[resid]:
            atom_index += 1
            peak[0] = atom_index
            line = "\t".join([str(p) for p in peak])
            xeasy_out.write(line + "\n")
    xeasy_out.close()



print "DEBUG: ENTERING ITERATION 2 matched_NOESY_residue_assignments_dict=", matched_NOESY_residue_assignments_dict
print "DEBUG: ENTERING ITERATION 2 clean_resid_assignments_dict=", clean_resid_assignments_dict

xeasy_fout = open(args.ABSOLUTE_MATCHES_FILE+".NOESY.iter2.xeasy", 'w')
sparky_fout = open(args.ABSOLUTE_MATCHES_FILE+".NOESY.iter2.sparky", 'w')
atom_index = 1
revordered_residue_keys = [(k[0], int(k[1:])) for k in matched_NOESY_residue_assignments_dict.keys() if k in absolute_matches_alignment_list]
revordered_residue_keys.sort(key=itemgetter(1), reverse=True)  # reverse ordered keys, namely they start from the last residue in the protein sequence
matched_i_to_iplus1_peaks_multidict = tree()    # residue -> i peak -> matched i+1 peak
residues_with_written_prevAssigned_peaks = []

for i_aa_resid in revordered_residue_keys:    # both the keys and the C,H assignments correspond to residue i-1
    i_residue = i_aa_resid[0]+str(i_aa_resid[1])
    all_possible_assignments_list = []  # list with all the peaks for this residue (TOCSY-NOESY, NOESY) and the respective possible C-H assignments
    print "ITERATION 2: Matching residue i ", i_residue," TOCSY PEAKS to residue's i NOESY peaks..."
    i_resid = get_resid_from_residue(i_residue)
    checked_residues_set.add(i_residue)
    aa_type = aa1to3_dict[get_aa_type_from_residue(i_residue)]
    all_allowed_carbons_list = [x for x in allowed_aa_atoms_dict[aa_type] if x[0] == "C" and not (aa_type in iter1_allowed_atoms_dict.keys() and x in iter1_allowed_atoms_dict[aa_type])]   # all allowed carbons for this aa type (exclude methyls)
    print "DEBUG: all_allowed_carbons_list=", all_allowed_carbons_list
    if match_NOESY_i_iplus1_peaks(i_residue, absolute_matches_alignment_list, matched_NOESY_residue_assignments_dict, matched_i_to_iplus1_peaks_multidict) == False:  # add labels to the i->i+1 matched peaks. If i_residue is a C-term residue or presends a gap, skip it from this iteration
        continue
    else:   # otherwise proceed to assigned but not from the extreme carbons in exclude_from_i_iplus_matching_dict
        assigned_carbons_list = []  # list of carbons which have ALL their protons assigned in TOCSY
        partly_assigned_carbons_list = []   # list of methylene carbons for which only one proton frequency was found in TOCSY
        missing_carbons_list = []
        tmp_missing_carbons_list, tmp_partly_assigned_carbons_list = get_missing_carbons_from_xeasy_dict(aa_type, i_resid, clean_resid_assignments_dict, matched_NOESY_residue_assignments_dict)
        print "DEBUG: tmp_missing_carbons_list=", tmp_missing_carbons_list, "tmp_partly_assigned_carbons_list=", tmp_partly_assigned_carbons_list
        for carbon in tmp_missing_carbons_list:
            if aa_type in iter1_allowed_atoms_dict.keys() and carbon in iter1_allowed_atoms_dict[aa_type]:     # to exclude methyls from this iteration
                continue
            if not aa_type in exclude_from_i_iplus_matching_dict.keys() or (not carbon in exclude_from_i_iplus_matching_dict[aa_type]):
                missing_carbons_list.append(carbon)
        for carbon in tmp_partly_assigned_carbons_list:
            if not aa_type in exclude_from_i_iplus_matching_dict.keys() or (not carbon in exclude_from_i_iplus_matching_dict[aa_type]):
                partly_assigned_carbons_list.append(carbon)
    print "DEBUG: dict key i_residue", i_residue, "missing_carbons_list=", missing_carbons_list
    
    unassigned_NOESY_peaks_list = []    # ini
    if len(missing_carbons_list) != 0:   # IF THERE ARE MISSING CARBONS, TRY TO FIND THEM (EXCLUDE THE ASSIGNED OR MATCHED NOESY PEAKS)
        print "Residue ", i_residue, " was not completely assigned in ITERATION 1. Using NOESY to assign resonances to the rest Carbon types."
        try:    # DANGEROUS but works for the N-term residue (A200 in aLP)
            unassigned_NOESY_peaks_list = [peak for peak in matched_NOESY_residue_assignments_dict[i_residue] if peak[0]=="i+1" and peak[1]=="?" and peak[3]=="?"]
        except KeyError:
            continue
        print "DEBUG: unassigned_NOESY_peaks_list=", unassigned_NOESY_peaks_list
        for NOESY_peak in unassigned_NOESY_peaks_list:  # make C-H type predictions for each unassigned NOESY peak
            NOESY_Creson = NOESY_peak[2]
            NOESY_Hreson = NOESY_peak[4]
            all_possible_assignments_list.extend(get_probabilities_from_H_C_resonpair_2Dhist(i_residue, NOESY_Hreson, NOESY_Creson, missing_carbons_list,
                                                            partly_assigned_carbons_list, args, aa_CHpair_binProbabilityList_multidict, aa_carbon_binDensityList_multidict,
                                                            NOESY_residue_peak_intensity_multidict, matched_i_to_iplus1_peaks_multidict, aa2pthres_dict,
                                                            aa2ithres_dict, protein_alignment_list, absolute_matches_alignment_list, 
                                                            aa_hydrogen_binDensityList_multidict, iteration=2))
        [assignment.append('NOESY') for assignment in all_possible_assignments_list]    # we don't need to created a new list, the original all_possible_assignments_list will be updated
        print "DEBUG: i_residue=", i_residue, "all_possible_assignments_list=", all_possible_assignments_list
        written_nucleiNames_list = []
        if i_resid in clean_resid_assignments_dict.keys():
            written_nucleiNames_list = [line[3]  for line in clean_resid_assignments_dict[i_resid]]   # list of nuclei that have been written already, to avoid duplication
        if len(all_possible_assignments_list) > 0:  # IF WE HAVE FOUND POSSIBLE NOESY PEAK ASSIGNMENTS OF THE MISSING CARBONS
            written_nucleiNames_list.extend(select_best_NOESY_peak_combination(i_residue, all_possible_assignments_list, partly_assigned_carbons_list, atom_index, xeasy_fout, sparky_fout,
                                       protein_alignment_list, absolute_matches_alignment_list, args, NOESY_residue_NHresonances_dict,
                                       matched_NOESY_residue_assignments_dict, TOCSY_assigned_residue_list, iter1_allowed_atoms_dict, iteration=2))
        residues_with_written_prevAssigned_peaks.append(i_residue)
        write_peaks_assigned_in_previous_iteration(i_residue, clean_resid_assignments_dict, atom_index, xeasy_fout, 
                                               NOESY_residue_NHresonances_dict, residues_with_written_NH_list)
    
    elif len(missing_carbons_list) == 0 and len(partly_assigned_carbons_list) == 0: # NOTHING TO DO FOR THIS RESIDUE
        residues_with_written_prevAssigned_peaks.append(i_residue)
        write_peaks_assigned_in_previous_iteration(i_residue, clean_resid_assignments_dict, atom_index, xeasy_fout, 
                                               NOESY_residue_NHresonances_dict, residues_with_written_NH_list)
    
    else:   # IF ALL THE CARBONS WERE ASSIGNED, FIND & SAVE PARTLY ASSIGNED CARBONS, SAVE THE MATCHED NOESY RESONANCES, AND IF NOT MATCHED, SAVE THE TOCSY RESONANCES
        print "Residue ", i_residue, " was completely assigned in ITERATION 1."
        try:
            aa_type = aa1to3_dict[get_aa_type_from_residue(i_residue)]
        except (IndexError, KeyError):    # skip residue i
            print "EXITING: COULD NOT FIND AA_TYPE OF RESIDUE", i_residue
            sys.exit(1)
        
        assigned_peaks_list = get_peaks_from_xeasy(i_resid, clean_resid_assignments_dict, atom_index, xeasy_fout, 
                                               NOESY_residue_NHresonances_dict, residues_with_written_NH_list)     # get all assigned peaks of this residue from xeasy format
        for assigned_peak in assigned_peaks_list:
            spectrum_type = assigned_peak[5]
            Cname = assigned_peak[1]
            Creson = assigned_peak[2]    # use the TOCSY resonance by default
            Hname = assigned_peak[3]
            Hreson = assigned_peak[4]
            all_possible_assignments_list.append([aa_type, Cname, Hname, 10000000.0, Hreson, Creson, spectrum_type])  # set the probability to 1 cause it has been already assigned in iteration 1
            print "DEBUG: i_residue=", i_residue, "aa_type=", aa_type, "all_possible_assignments_list=", all_possible_assignments_list
        written_nucleiNames_list = []
        if i_resid in clean_resid_assignments_dict.keys():
            written_nucleiNames_list = [line[3]  for line in clean_resid_assignments_dict[i_resid]]   # list of nuclei that have been written already, to avoid duplication
        if len(all_possible_assignments_list) > 0:  # IF WE HAVE FOUND POSSIBLE NOESY PEAK ASSIGNMENTS OF THE MISSING CARBONS
            written_nucleiNames_list = select_best_NOESY_peak_combination(i_residue, all_possible_assignments_list, partly_assigned_carbons_list, atom_index, xeasy_fout, sparky_fout,
                                       protein_alignment_list, absolute_matches_alignment_list, args, NOESY_residue_NHresonances_dict,
                                       matched_NOESY_residue_assignments_dict, TOCSY_assigned_residue_list, iter1_allowed_atoms_dict, iteration=2)
    remove_iplus1_labels(i_residue, matched_NOESY_residue_assignments_dict, protein_alignment_list, absolute_matches_alignment_list, args)     # remove all the 'i+1', 'i-1' labels from this residue and residue i+1 (if applicable)
        
print "DEBUG iteration2 point2: matched_NOESY_residue_assignments_dict=", matched_NOESY_residue_assignments_dict

for i_aa_resid in revordered_residue_keys:
    i_residue = i_aa_resid[0]+str(i_aa_resid[1])
    aa_type = aa1to3_dict[i_aa_resid[0]]
    if not i_residue in residues_with_written_prevAssigned_peaks:
        write_peaks_assigned_in_previous_iteration(i_residue, clean_resid_assignments_dict, atom_index, xeasy_fout, 
                                               NOESY_residue_NHresonances_dict, residues_with_written_NH_list)


for residue in absolute_matches_alignment_list:
    if not residue in checked_residues_set and residue != '-':
        print "WARNING: Residue ", residue, " was not checked!!!"


write_NH_of_residues_without_CH(absolute_matches_alignment_list, NOESY_residue_NHresonances_dict, TOCSY_residue_NHresonances_dict,
                                    residues_with_written_NH_list, atom_index, xeasy_fout, sparky_fout, protein_alignment_list, args, 
                                    Cterm_residue_set, iteration=2)

xeasy_fout.close()
sparky_fout.close()


clean_resid_assignments_dict = clean_xeasy_file(args.ABSOLUTE_MATCHES_FILE+".NOESY.iter2.xeasy")  
for resid, peak_list in clean_resid_assignments_dict.items():
    peak_list.sort(key=itemgetter(3))
    clean_resid_assignments_dict[resid] = peak_list

print "DEBUG: ITERATION 2 point 1 clean_resid_assignments_dict=", clean_resid_assignments_dict
xeasy_out = open(args.ABSOLUTE_MATCHES_FILE+".NOESY.cleaned.iter2.xeasy", 'w')
atom_index = 0  # reinitialize atom index
for resid in clean_resid_assignments_dict.keys():
    for peak in clean_resid_assignments_dict[resid]:
        atom_index += 1
        peak[0] = atom_index
        line = "\t".join([str(p) for p in peak])
        xeasy_out.write(line + "\n")
xeasy_out.close()


print "DEBUG: args.user_TOCSY_fname=", args.user_TOCSY_fname, "args.user_NOESY_fname=", args.user_NOESY_fname
if args.user_TOCSY_fname and args.user_NOESY_fname:
    print "STARTING PROOF-READING NOESY CHEMICAL SHIFT ASSIGNMENTS"
    aa_equivalentCarbons_dict = {}
    aa_equivalentCarbons_dict["LEU"] = [["CD1", "CD2"]]
    aa_equivalentCarbons_dict["VAL"] = [["CG1", "CG2"]]
    
    user_TOCSY_residue_assignments_dict, user_TOCSY_residue_NHresonances_dict, user_i_to_iminus1_dict, user_patched_residues_list, user_TOCSY_residue_peak_intensity_multidict = read_assigned_TOCSY_file(args.user_TOCSY_fname, "TOCSY", absolute_RIGmatches_alignment_list, absolute_matches_alignment_list)
    print "DEBUG: user_TOCSY_residue_assignments_dict=", user_TOCSY_residue_assignments_dict
    print "DEBUG: user_TOCSY_residue_NHresonances_dict=", user_TOCSY_residue_NHresonances_dict
    for k,v in user_TOCSY_residue_assignments_dict.items():
        print k, "-->", v
        
    user_NOESY_residue_assignments_dict, user_NOESY_residue_NHresonances_dict, user_original_NOESY_peaks_dict, user_NOESY_residue_peak_instensity_dict = read_assigned_TOCSY_file(args.user_NOESY_fname, 'NOESY', absolute_RIGmatches_alignment_list, absolute_matches_alignment_list)
    print "DEBUG: user_NOESY_residue_assignments_dict=", user_NOESY_residue_assignments_dict
    print "DEBUG: user_NOESY_residue_NHresonances_dict=", user_NOESY_residue_NHresonances_dict
    print "DEBUG: user_original_NOESY_peaks_dict=", user_original_NOESY_peaks_dict
    for residue in user_NOESY_residue_assignments_dict:
        if residue[0] in ['Y', 'F']:
            for assignment in user_NOESY_residue_assignments_dict[residue]:
                if assignment[1] == "CD":
                    assignment[1] = "CD1"
                elif assignment[1] == "QD":
                    assignment[1] = "HD1"
                elif assignment[1] == "CE":
                    assignment[1] = "CE1"
                elif assignment[1] == "QE":
                    assignment[1] = "HE1"
    for residue in user_TOCSY_residue_assignments_dict:
        if residue[0] in ['Y', 'F']:
            for assignment in user_TOCSY_residue_assignments_dict[residue]:
                if assignment[1] == "CD":
                    assignment[1] = "CD1"
                elif assignment[1] == "QD":
                    assignment[1] = "HD1"
                elif assignment[1] == "CE":
                    assignment[1] = "CE1"
                elif assignment[1] == "QE":
                    assignment[1] = "HE1"
    for k,v in user_NOESY_residue_assignments_dict.items():
        print k, "-->", v
    
    user_NOESY_resid_assignments_dict = {}  # same as user_NOESY_resid_assignments_dict, but with keys resids instead of residues
    for k,v in user_NOESY_residue_assignments_dict.items():
        user_NOESY_resid_assignments_dict[int(k[1:])] = v
    user_TOCSY_resid_assignments_dict = {}  # same as user_TOCSY_resid_assignments_dict, but with keys resids instead of residues
    for k,v in user_TOCSY_residue_assignments_dict.items():
        user_TOCSY_resid_assignments_dict[int(k[1:])] = v
    
    print "DEBUG: user_NOESY_resid_assignments_dict=", user_NOESY_resid_assignments_dict
    print "DEBUG: ITERATION 2 point 2 clean_resid_assignments_dict=", clean_resid_assignments_dict
    
    user_VAL_LEU_NOESY_resid_assignments_dict, user_VAL_LEU_TOCSY_resid_assignments_dict = create_equivalentCarbons_assignments_dict(user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict, protein_alignment_list,
                                              absolute_matches_alignment_list, aa_equivalentCarbons_dict, args)
    print "DEBUG: user_VAL_LEU_NOESY_resid_assignments_dict=", user_VAL_LEU_NOESY_resid_assignments_dict
    print "DEBUG: user_VAL_LEU_TOCSY_resid_assignments_dict=", user_VAL_LEU_TOCSY_resid_assignments_dict
    
    
    proofread_clean_resid_assignments_dict = copy.deepcopy(clean_resid_assignments_dict)   # same as clean_resid_assignments_dict but with comments CORRECT or WRONG at the carbon lines
    for resid in proofread_clean_resid_assignments_dict.keys():
        print "DEBUG ITERATION2: proof reading resid=", resid
        resname = proofread_clean_resid_assignments_dict[resid][0][6]
        if resname in aa_equivalentCarbons_dict.keys(): # If this is a LEU or VAL
            Carbon_pair = aa_equivalentCarbons_dict[resname][0]
            if len([a for a in proofread_clean_resid_assignments_dict[resid] if a[3] in Carbon_pair]) > 0: # AND if the methyl carbons have been assigned
                Proton_pair = ["H"+C[1:] for C in Carbon_pair]   # protons of the equivalent Carbons
                program_assignments_list = [assignment for assignment in proofread_clean_resid_assignments_dict[resid] if assignment[3] in Carbon_pair+Proton_pair]   # keep only the resonances assigned by the program to the equivalent methyl carbons and protons
                proof_read_equivalentMethyls(program_assignments_list, resid, user_VAL_LEU_NOESY_resid_assignments_dict, user_VAL_LEU_TOCSY_resid_assignments_dict,
                                             Carbon_pair, Proton_pair)   # it doesn't matter if 1 or 2 of the peaks have been proof-read in iteration 1 & 2
        for assignment in proofread_clean_resid_assignments_dict[resid]:
            print "DEBUG ITERATION2: proof reading assignment=", assignment
            spectrum_type = re.sub(r' iter[1-9]$', '', assignment[7])    # remove the "iter[1-9]" suffix to be able to distinguish the type of assignment (NOESY, NOESY+TOCSY, etc.)
            comment = None
            aa_type = assignment[6]
            nucleus = assignment[3] # the nucleus name
            if aa_type in aa_equivalentCarbons_dict.keys() and nucleus in aa_equivalentCarbons_dict[aa_type][0]:   # if this peak has been assigned by the program to LEU or VAL CD1, CD2 or CG1,CG2 respectively, skip it
                continue
            
            if aa_type in aatype_carbon_nongeminalHname_multidict.keys() and nucleus in aatype_carbon_nongeminalHname_multidict[aa_type].keys():
                C_assignment = assignment # the program-made carbon assignment
                H_assignment = None # default value for the program-made proton assignment
                print "DEBUG: resid=", resid, "proofread_clean_resid_assignments_dict[resid]=", proofread_clean_resid_assignments_dict[resid]
                for a in proofread_clean_resid_assignments_dict[resid]:
                    if ((type(aatype_carbon_nongeminalHname_multidict[aa_type][nucleus]) == list and a[3] in aatype_carbon_nongeminalHname_multidict[aa_type][nucleus]) or
                        (type(aatype_carbon_nongeminalHname_multidict[aa_type][nucleus]) == str and a[3] == aatype_carbon_nongeminalHname_multidict[aa_type][nucleus])):
                        H_assignment = a
                        break
                print "DEBUG: C_assignment=", C_assignment, "H_assignment=", H_assignment
                comment = ""
                if spectrum_type in ['NOESY', 'NOESY (C-term hanging residue)', 'NOESY + TOCSY-NOESY', 'TOCSY-NOESY', 'TOCSY-NOESY + NOESY', 'TOCSY (unmatched) + TOCSY-NOESY', 'TOCSY-NOESY + TOCSY (unmatched)']:    # if this resonance comes only from NOESY
                    user_assignment = get_Carbon_resonance(resid, nucleus, 'NOESY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict,
                         aa_equivalentCarbons_dict, protein_alignment_list, absolute_matches_alignment_list, args, get_assignment=True)
                    if user_assignment == False:    # if no user assignment was made for this Carbon, add no extra comment
                        print "WARNING: no user assignment made in NOESY for resid", resid, "and Carbon", nucleus,". Checking TOCSY..."
                        user_assignment = get_Carbon_resonance(resid, nucleus, 'TOCSY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict,
                         aa_equivalentCarbons_dict, protein_alignment_list, absolute_matches_alignment_list, args, get_assignment=True)
                        if user_assignment == False:
                            print "WARNING: no user assignment made also in TOCSY for resid", resid, "and Carbon", nucleus,". Skipping nucleus."
                            assignment[7] += " <NOT ASSIGNED>"
                            continue
                    user_Creson = user_assignment[2]
                    user_Hreson = user_assignment[4]
                    print "DEBUG: user_Creson=", user_Creson, "user_Hreson=", user_Hreson
                    if approx_equal(C_assignment[1], user_Creson, 0.4) and approx_equal(H_assignment[1], user_Hreson, 0.04): # if the user-made assignment agrees with ours, add to the comment "CORRECT"
                        print "DEBUG: CORRECT C_assignment found:", C_assignment
                        print "DEBUG: CORRECT H_assignment found:", H_assignment
                        comment = " <CORRECT>"
                    else:
                        print "DEBUG: WRONG C_assignment found:", C_assignment
                        print "DEBUG: WRONG H_assignment found:", H_assignment
                        comment = " <WRONG>"
                elif spectrum_type in ['TOCSY (average)', 'TOCSY (unmatched)']:
                    user_assignment = get_Carbon_resonance(resid, nucleus, 'TOCSY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict,
                         aa_equivalentCarbons_dict, protein_alignment_list, absolute_matches_alignment_list, args, get_assignment=True)
                    if user_assignment == False:    # if no user assignment was made for this Carbon, add no extra comment
                        print "WARNING: no user assignment made in TOCSY for resid", resid, "and Carbon", nucleus,". Checking NOESY..."
                        user_assignment = get_Carbon_resonance(resid, nucleus, 'NOESY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict,
                         aa_equivalentCarbons_dict, protein_alignment_list, absolute_matches_alignment_list, args, get_assignment=True)
                        if user_assignment == False:
                            print "WARNING: no user assignment made also in NOESY for resid", resid, "and Carbon", nucleus,". Skipping nucleus."
                            assignment[7] += " <NOT ASSIGNED>"
                            continue
                    user_Creson = user_assignment[2]
                    user_Hreson = user_assignment[4]
                    if approx_equal(C_assignment[1], user_Creson, 0.4) and approx_equal(H_assignment[1], user_Hreson, 0.04): # if the user-made assignment agrees with ours, add to the comment "CORRECT"
                        print "DEBUG: CORRECT C_assignment found:", C_assignment
                        print "DEBUG: CORRECT H_assignment found:", H_assignment
                        comment = " <CORRECT>"
                    else:
                        print "DEBUG: WRONG C_assignment found:", C_assignment
                        print "DEBUG: WRONG H_assignment found:", H_assignment
                        comment = " <WRONG>"
                        
                if comment == "":
                    comment = " <NOT ASSIGNED>"
                C_assignment[7] += comment
                H_assignment[7] += comment
                comment = None  # a trick to add no more comments to "assignment" at the end of this block of code
            
            elif nucleus[0] == 'C':   # only if the assignment line belongs to a Carbon then check if the user has made the same assignment
                comment = ""
                if spectrum_type in ['NOESY', 'NOESY (C-term hanging residue)', 'NOESY + TOCSY-NOESY', 'TOCSY-NOESY', 'TOCSY-NOESY + NOESY', 'TOCSY (unmatched) + TOCSY-NOESY', 'TOCSY-NOESY + TOCSY (unmatched)']:    # if this resonance comes only from NOESY
                    user_assignment = get_Carbon_resonance(resid, nucleus, 'NOESY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict,
                         aa_equivalentCarbons_dict, protein_alignment_list, absolute_matches_alignment_list, args)
                    if user_assignment == False:    # if no user assignment was made for this Carbon, add no extra comment
                        print "WARNING: no user assignment made in NOESY for resid", resid, "and Carbon", nucleus,". Checking TOCSY..."
                        user_assignment = get_Carbon_resonance(resid, nucleus, 'TOCSY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict,
                         aa_equivalentCarbons_dict, protein_alignment_list, absolute_matches_alignment_list, args)
                        if user_assignment == False:
                            print "WARNING: no user assignment made also in TOCSY for resid", resid, "and Carbon", nucleus,". Skipping nucleus."
                            assignment[7] += " <NOT ASSIGNED>"
                            continue
                    if approx_equal(assignment[1], user_assignment, 0.4): # if the user-made assignment agrees with ours, add to the comment "CORRECT"
                        print "DEBUG: CORRECT assignment found:", assignment
                        comment = " <CORRECT>"
                    else:
                        print "DEBUG: WRONG assignment found:", assignment
                        comment = " <WRONG>"
                elif spectrum_type in ['TOCSY (average)', 'TOCSY (unmatched)']:
                    user_assignment = get_Carbon_resonance(resid, nucleus, 'TOCSY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict,
                         aa_equivalentCarbons_dict, protein_alignment_list, absolute_matches_alignment_list, args)
                    if user_assignment == False:    # if no user assignment was made for this Carbon, add no extra comment
                        print "WARNING: no user assignment made in TOCSY for resid", resid, "and Carbon", nucleus,". Checking NOESY..."
                        user_assignment = get_Carbon_resonance(resid, nucleus, 'NOESY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict,
                         aa_equivalentCarbons_dict, protein_alignment_list, absolute_matches_alignment_list, args)
                        if user_assignment == False:
                            print "WARNING: no user assignment made also in NOESY for resid", resid, "and Carbon", nucleus,". Skipping nucleus."
                            assignment[7] += " <NOT ASSIGNED>"
                            continue
                    if approx_equal(assignment[1], user_assignment, 0.4): # if the user-made assignment agrees with ours, add to the comment "CORRECT"
                        print "DEBUG: CORRECT assignment found:", assignment
                        comment = " <CORRECT>"
                    else:
                        print "DEBUG: WRONG assignment found:", assignment
                        comment = " <WRONG>"
                
            if comment == None:     # this is a proton or a methyl Carbon or a Carbon with only one proton
                continue
            elif comment == "":
                comment = " <NOT ASSIGNED>"
            assignment[7] += comment
    
    
    print "DEBUG: proofread_clean_resid_assignments_dict=", proofread_clean_resid_assignments_dict
    
    xeasy_out = open(args.ABSOLUTE_MATCHES_FILE+".NOESY.proofread.iter2.xeasy", 'w')
    atom_index = 0  # reinitialize atom index
    for resid in proofread_clean_resid_assignments_dict.keys():
        for peak in proofread_clean_resid_assignments_dict[resid]:
            atom_index += 1
            peak[0] = atom_index
            line = "\t".join([str(p) for p in peak])
            xeasy_out.write(line + "\n")
    xeasy_out.close()


print "DEBUG: ENTERING ITERATION 3 matched_NOESY_residue_assignments_dict=", matched_NOESY_residue_assignments_dict
print "DEBUG: ENTERING ITERATION 3 clean_resid_assignments_dict=", clean_resid_assignments_dict

xeasy_fout = open(args.ABSOLUTE_MATCHES_FILE+".NOESY.iter3.xeasy", 'w')
sparky_fout = open(args.ABSOLUTE_MATCHES_FILE+".NOESY.iter3.sparky", 'w')
atom_index = 1
revordered_residue_keys = [(k[0], int(k[1:])) for k in matched_NOESY_residue_assignments_dict.keys() if k in absolute_matches_alignment_list]
revordered_residue_keys.sort(key=itemgetter(1), reverse=True)  # reverse ordered keys, namely they start from the last residue in the protein sequence
residues_with_written_prevAssigned_peaks = []

selected_revordered_residue_keys = []
for aa_type  in ["ALA", "THR", "VAL", "ILE", "LEU"]:
    for i_aa_resid in revordered_residue_keys:
        i_residue = i_aa_resid[0]+str(i_aa_resid[1])
        aa = aa1to3_dict[get_aa_type_from_residue(i_residue)]
        if aa == aa_type:
            selected_revordered_residue_keys.append(i_aa_resid)
print "DEBUG: selected_revordered_residue_keys=", selected_revordered_residue_keys

for i_aa_resid in selected_revordered_residue_keys:    # both the keys and the C,H assignments correspond to residue i-1
    if aa_type in aa2ithres_dict.keys():
        istart = int(10 * aa2ithres_dict[aa_type])
    else:
        istart = 5
    for I_THRESHOLD in reversed(range(1, istart+1, 1)):
        I_THRESHOLD *= 0.1
        i_residue = i_aa_resid[0]+str(i_aa_resid[1])
        aa_type = aa1to3_dict[i_aa_resid[0]]
        i_resid = get_resid_from_residue(i_residue)
        if not aa_type in iter1_allowed_atoms_dict.keys():
            continue
        all_possible_assignments_list = []  # list with all the peaks for this residue (TOCSY-NOESY, NOESY) and the respective possible C-H assignments
        print "ITERATION 3: Matching residue i ", i_residue," TOCSY PEAKS to residue's i NOESY peaks... I_THRESHOLD=", I_THRESHOLD
        checked_residues_set.add(i_residue)
        tmp_missing_carbons_list, tmp_partly_assigned_carbons_list = get_missing_carbons_from_xeasy_dict(aa_type, i_resid, clean_resid_assignments_dict, matched_NOESY_residue_assignments_dict)
        missing_carbons_list = [c for c in tmp_missing_carbons_list if c in iter1_allowed_atoms_dict[aa_type]]
        partly_assigned_carbons_list = [c for c in tmp_partly_assigned_carbons_list if c in iter1_allowed_atoms_dict[aa_type]]
        print "DEBUG: missing_carbons_list=", missing_carbons_list, "partly_assigned_carbons_list=", partly_assigned_carbons_list
        unassigned_NOESY_peaks_list = []
        if len(missing_carbons_list) != 0:   # IF THERE ARE MISSING CARBONS IN ITERATION 1, TRY TO FIND THEM IN NOESY (EXCLUDE THE ASSIGNED OR MATCHED NOESY PEAKS)
            print "Residue ", i_residue, " was not completely assigned in iteration 1. Using NOESY to assigned resonances to the rest Carbon types."
            try:    # DANGEROUS but works for the N-term residue (A200 in aLP)
                unassigned_NOESY_peaks_list = [peak for peak in matched_NOESY_residue_assignments_dict[i_residue] if peak[0]=="?" and peak[1]=="?" and peak[3]=="?"]
            except KeyError:    # if this residue is not in the alignment (not in matched_NOESY_residue_assignments_dict), then we cannot assign the partly assigned Carbons and hence we move on to the next
                print "Residue ", i_residue, " is either not in the alignment or doesn't have NOESY. Therefore it is impossible to assign partly assigned Carbon peaks."
                continue
            for NOESY_peak in unassigned_NOESY_peaks_list:  # make C-H type predictions for each unassigned NOESY peak
                NOESY_Creson = NOESY_peak[2]
                NOESY_Hreson = NOESY_peak[4]
                all_possible_assignments_list.extend(get_probabilities_from_H_C_resonpair_2Dhist(i_residue, NOESY_Hreson, NOESY_Creson, missing_carbons_list,
                                                                partly_assigned_carbons_list, args, aa_CHpair_binProbabilityList_multidict, aa_carbon_binDensityList_multidict,
                                                                NOESY_residue_peak_intensity_multidict, matched_i_to_iplus1_peaks_multidict, aa2pthres_dict,
                                                                aa2ithres_dict, protein_alignment_list, absolute_matches_alignment_list, 
                                                                aa_hydrogen_binDensityList_multidict,  iteration=3, I_THRESHOLD=I_THRESHOLD))
            [assignment.append('NOESY') for assignment in all_possible_assignments_list]    # we don't need to created a new list, the original all_possible_assignments_list will be updated
            print "DEBUG: i_residue=", i_residue, "all_possible_assignments_list=", all_possible_assignments_list
            if i_resid in clean_resid_assignments_dict.keys():  # if this is a C-term residue that is analyzed for the first time in iter3
                written_nucleiNames_list = [line[3]  for line in clean_resid_assignments_dict[i_resid]]   # list of nuclei that have been written already, to avoid duplication
            else:
                written_nucleiNames_list = []
            if len(all_possible_assignments_list) > 0:  # IF WE HAVE FOUND POSSIBLE NOESY PEAK ASSIGNMENTS OF THE MISSING CARBONS
                written_nucleiNames_list.extend(select_best_NOESY_peak_combination(i_residue, all_possible_assignments_list, partly_assigned_carbons_list, atom_index, xeasy_fout, sparky_fout,
                                       protein_alignment_list, absolute_matches_alignment_list, args, NOESY_residue_NHresonances_dict,
                                       matched_NOESY_residue_assignments_dict, TOCSY_assigned_residue_list, iter1_allowed_atoms_dict, iteration=3))
            residues_with_written_prevAssigned_peaks.append(i_residue)
            write_peaks_assigned_in_previous_iteration(i_residue, clean_resid_assignments_dict, atom_index, xeasy_fout, 
                                               NOESY_residue_NHresonances_dict, residues_with_written_NH_list)
        elif len(missing_carbons_list) == 0 and len(partly_assigned_carbons_list) == 0: # NOTHING TO DO FOR THIS RESIDUE
            residues_with_written_prevAssigned_peaks.append(i_residue)
            write_peaks_assigned_in_previous_iteration(i_residue, clean_resid_assignments_dict, atom_index, xeasy_fout, 
                                               NOESY_residue_NHresonances_dict, residues_with_written_NH_list)
        else:   # IF ALL THE CARBONS WERE ASSIGNED IN ITERATION 1, SAVE THEM, BUT ALSO FIND IN NOESY AND SAVE PARTLY ASSIGNED CARBONS
            print "Residue ", i_residue, " was completely assigned in iteration 1."
            if not i_residue in matched_NOESY_residue_assignments_dict.keys():
                print "Residue ", i_residue, " is either not in the alignment or doesn't have NOESY. Therefore it is impossible to assign partly assigned Carbon peaks."
                continue
            try:
                aa_type = aa1to3_dict[get_aa_type_from_residue(i_residue)]
            except (IndexError, KeyError):    # skip residue i
                print "EXITING: COULD NOT FIND AA_TYPE OF RESIDUE", i_residue
                sys.exit(1)
            assigned_peaks_list = get_peaks_from_xeasy(i_resid, clean_resid_assignments_dict)     # get all assigned peaks of this residue from xeasy format
            for assigned_peak in assigned_peaks_list:
                spectrum_type = assigned_peak[5]
                Cname = assigned_peak[1]
                Creson = assigned_peak[2]    # use the TOCSY resonance by default
                Hname = assigned_peak[3]
                Hreson = assigned_peak[4]
                all_possible_assignments_list.append([aa_type, Cname, Hname, 10000000.0, Hreson, Creson, spectrum_type])  # set the probability to 1 cause it has been already assigned in iteration 1
                print "DEBUG: i_residue=", i_residue, "aa_type=", aa_type, "all_possible_assignments_list=", all_possible_assignments_list
            written_nucleiNames_list = []
            if i_resid in clean_resid_assignments_dict.keys():
                written_nucleiNames_list = [line[3]  for line in clean_resid_assignments_dict[i_resid]]  # list of nuclei that have been written already, to avoid duplication
            if len(all_possible_assignments_list) > 0:  # IF WE HAVE FOUND POSSIBLE NOESY PEAK ASSIGNMENTS OF THE MISSING CARBONS
                written_nucleiNames_list.extend(select_best_NOESY_peak_combination(i_residue, all_possible_assignments_list, partly_assigned_carbons_list, atom_index, xeasy_fout, sparky_fout,
                                       protein_alignment_list, absolute_matches_alignment_list, args, NOESY_residue_NHresonances_dict,
                                       matched_NOESY_residue_assignments_dict, TOCSY_assigned_residue_list, iter1_allowed_atoms_dict, iteration=3))


for i_aa_resid in revordered_residue_keys:
    i_residue = i_aa_resid[0]+str(i_aa_resid[1])
    aa_type = aa1to3_dict[i_aa_resid[0]]
    if not i_residue in residues_with_written_prevAssigned_peaks:
        write_peaks_assigned_in_previous_iteration(i_residue, clean_resid_assignments_dict, atom_index, xeasy_fout, 
                                               NOESY_residue_NHresonances_dict, residues_with_written_NH_list)


for residue in absolute_matches_alignment_list:
    if not residue in checked_residues_set and residue != '-':
        print "WARNING: Residue ", residue, " was not checked!!!"


write_NH_of_residues_without_CH(absolute_matches_alignment_list, NOESY_residue_NHresonances_dict, TOCSY_residue_NHresonances_dict,
                                    residues_with_written_NH_list, atom_index, xeasy_fout, sparky_fout, protein_alignment_list, args, 
                                    Cterm_residue_set, iteration=3)

xeasy_fout.close()
sparky_fout.close()


clean_resid_assignments_dict = clean_xeasy_file(args.ABSOLUTE_MATCHES_FILE+".NOESY.iter3.xeasy")  # dict with all the assignments from each residue from iteration 1
for resid, peak_list in clean_resid_assignments_dict.items():
    peak_list.sort(key=itemgetter(3))
    clean_resid_assignments_dict[resid] = peak_list

print "DEBUG: ITERATION 3 point 1 clean_resid_assignments_dict=", clean_resid_assignments_dict

xeasy_out = open(args.ABSOLUTE_MATCHES_FILE+".NOESY.iter3.xeasy", 'w')
atom_index = 0  # reinitialize atom index
for resid in clean_resid_assignments_dict.keys():
    for peak in clean_resid_assignments_dict[resid]:
        atom_index += 1
        peak[0] = atom_index
        line = "\t".join([str(p) for p in peak])
        xeasy_out.write(line + "\n")
xeasy_out.close()



clean_resid_assignments_dict = clean_xeasy_file(args.ABSOLUTE_MATCHES_FILE+".NOESY.iter3.xeasy")  
for resid, peak_list in clean_resid_assignments_dict.items():
    peak_list.sort(key=itemgetter(3))
    clean_resid_assignments_dict[resid] = peak_list

print "DEBUG: ITERATION 3 point 2 clean_resid_assignments_dict=", clean_resid_assignments_dict

xeasy_out = open(args.ABSOLUTE_MATCHES_FILE+".NOESY.cleaned.iter3.xeasy", 'w')
atom_index = 0  # reinitialize atom index
for resid in clean_resid_assignments_dict.keys():
    for peak in clean_resid_assignments_dict[resid]:
        atom_index += 1
        peak[0] = atom_index
        line = "\t".join([str(p) for p in peak])
        xeasy_out.write(line + "\n")
xeasy_out.close()


print "DEBUG: args.user_TOCSY_fname=", args.user_TOCSY_fname, "args.user_NOESY_fname=", args.user_NOESY_fname
if args.user_TOCSY_fname and args.user_NOESY_fname:
    print "STARTING PROOF-READING NOESY CHEMICAL SHIFT ASSIGNMENTS"
    for resid in clean_resid_assignments_dict:
        if not resid in proofread_clean_resid_assignments_dict.keys():
            proofread_clean_resid_assignments_dict[resid] = clean_resid_assignments_dict[resid]
            continue
        for new_assignment in clean_resid_assignments_dict[resid]:
            new_reson = new_assignment[1]
            old_reson_list = [a[1] for a in proofread_clean_resid_assignments_dict[resid]]
            if not new_reson in old_reson_list:     # if this resonance was not found (newly assigned), add it to the dictionary
                proofread_clean_resid_assignments_dict[resid].append(new_assignment)
    print "DEBUG: ITERATION 3 point 3 clean_resid_assignments_dict=", clean_resid_assignments_dict
    
    
    
    for resid in proofread_clean_resid_assignments_dict.keys():
        print "DEBUG ITERATION3: proof reading resid=", resid
        resname = proofread_clean_resid_assignments_dict[resid][0][6]
        if resname in aa_equivalentCarbons_dict.keys(): # If this is a LEU or VAL
            Carbon_pair = aa_equivalentCarbons_dict[resname][0]
            if len([a for a in proofread_clean_resid_assignments_dict[resid] if a[3] in Carbon_pair]) > 0: # AND if the methyl carbons have been assigned
                Proton_pair = ["H"+C[1:] for C in Carbon_pair]   # protons of the equivalent Carbons
                program_assignments_list = [assignment for assignment in proofread_clean_resid_assignments_dict[resid] if assignment[3] in Carbon_pair+Proton_pair]   # keep only the resonances assigned by the program to the equivalent methyl carbons and protons
                proof_read_equivalentMethyls(program_assignments_list, resid, user_VAL_LEU_NOESY_resid_assignments_dict, user_VAL_LEU_TOCSY_resid_assignments_dict, Carbon_pair, Proton_pair)   # it doesn't matter if 1 or 2 of the peaks have been proof-read in iteration 1 & 2
            print "DEBUG ITERATION3: after proof_read_equivalentMethyls() proofread_clean_resid_assignments_dict[",resid,"]=", proofread_clean_resid_assignments_dict[resid]
        for assignment in proofread_clean_resid_assignments_dict[resid]:
            if '<' in assignment[7]:    # if this assignment has been proofread, then skip it
                continue
            print "DEBUG ITERATION3: proof reading new assignment=", assignment
            spectrum_type = re.sub(r' iter[1-9]', '', assignment[7])    # remove the "iter[1-9]" suffix to be able to distinguish the type of assignment (NEESY, NOESY+TOCSY, etc.)
            comment = None
            aa_type = assignment[6]
            nucleus = assignment[3] # the nucleus name
            if aa_type in aa_equivalentCarbons_dict.keys() and nucleus in aa_equivalentCarbons_dict[aa_type][0]:    # if this peak has been assigned by the program to LEU or VAL CD1, CD2 or CG1,CG2 respectively, skip it
                continue
            
            if aa_type in aatype_carbon_nongeminalHname_multidict.keys() and nucleus in aatype_carbon_nongeminalHname_multidict[aa_type].keys():
                C_assignment = assignment # the program-made carbon assignment
                H_assignment = None # default value for the program-made proton assignment
                print "DEBUG: resid=", resid, "proofread_clean_resid_assignments_dict[resid]=", proofread_clean_resid_assignments_dict[resid]
                for a in proofread_clean_resid_assignments_dict[resid]:
                    if ((type(aatype_carbon_nongeminalHname_multidict[aa_type][nucleus]) == list and a[3] in aatype_carbon_nongeminalHname_multidict[aa_type][nucleus]) or
                        (type(aatype_carbon_nongeminalHname_multidict[aa_type][nucleus]) == str and a[3] == aatype_carbon_nongeminalHname_multidict[aa_type][nucleus])):
                        H_assignment = a
                        break
                print "DEBUG: C_assignment=", C_assignment, "H_assignment=", H_assignment
                comment = ""
                if spectrum_type in ['NOESY', 'NOESY (C-term hanging residue)', 'NOESY + TOCSY-NOESY', 'TOCSY-NOESY', 'TOCSY-NOESY + NOESY', 'TOCSY (unmatched) + TOCSY-NOESY', 'TOCSY-NOESY + TOCSY (unmatched)']:    # if this resonance comes only from NOESY
                    user_assignment = get_Carbon_resonance(resid, nucleus, 'NOESY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict,
                         aa_equivalentCarbons_dict, protein_alignment_list, absolute_matches_alignment_list, args, get_assignment=True)
                    if user_assignment == False:    # if no user assignment was made for this Carbon, add no extra comment
                        print "WARNING: no user assignment made in NOESY for resid", resid, "and Carbon", nucleus,". Checking TOCSY..."
                        user_assignment = get_Carbon_resonance(resid, nucleus, 'TOCSY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict,
                         aa_equivalentCarbons_dict, protein_alignment_list, absolute_matches_alignment_list, args, get_assignment=True)
                        if user_assignment == False:
                            print "WARNING: no user assignment made also in TOCSY for resid", resid, "and Carbon", nucleus,". Skipping nucleus."
                            assignment[7] += " <NOT ASSIGNED>"
                            continue
                    user_Creson = user_assignment[2]
                    user_Hreson = user_assignment[4]
                    print "DEBUG: user_Creson=", user_Creson, "user_Hreson=", user_Hreson
                    if approx_equal(C_assignment[1], user_Creson, 0.4) and approx_equal(H_assignment[1], user_Hreson, 0.04): # if the user-made assignment agrees with ours, add to the comment "CORRECT"
                        print "DEBUG: CORRECT C_assignment found:", C_assignment
                        print "DEBUG: CORRECT H_assignment found:", H_assignment
                        comment = " <CORRECT>"
                    else:
                        print "DEBUG: WRONG C_assignment found:", C_assignment
                        print "DEBUG: WRONG H_assignment found:", H_assignment
                        comment = " <WRONG>"
                elif spectrum_type in ['TOCSY (average)', 'TOCSY (unmatched)']:
                    user_assignment = get_Carbon_resonance(resid, nucleus, 'TOCSY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict,
                         aa_equivalentCarbons_dict, protein_alignment_list, absolute_matches_alignment_list, args, get_assignment=True)
                    if user_assignment == False:    # if no user assignment was made for this Carbon, add no extra comment
                        print "WARNING: no user assignment made in TOCSY for resid", resid, "and Carbon", nucleus,". Checking NOESY..."
                        user_assignment = get_Carbon_resonance(resid, nucleus, 'NOESY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict,
                         aa_equivalentCarbons_dict, protein_alignment_list, absolute_matches_alignment_list, args, get_assignment=True)
                        if user_assignment == False:
                            print "WARNING: no user assignment made also in NOESY for resid", resid, "and Carbon", nucleus,". Skipping nucleus."
                            assignment[7] += " <NOT ASSIGNED>"
                            continue
                    user_Creson = user_assignment[2]
                    user_Hreson = user_assignment[4]
                    if approx_equal(C_assignment[1], user_Creson, 0.4) and approx_equal(H_assignment[1], user_Hreson, 0.04): # if the user-made assignment agrees with ours, add to the comment "CORRECT"
                        print "DEBUG: CORRECT C_assignment found:", C_assignment
                        print "DEBUG: CORRECT H_assignment found:", H_assignment
                        comment = " <CORRECT>"
                    else:
                        print "DEBUG: WRONG C_assignment found:", C_assignment
                        print "DEBUG: WRONG H_assignment found:", H_assignment
                        comment = " <WRONG>"
                        
                if comment == "":
                    comment = " <NOT ASSIGNED>"
                C_assignment[7] += comment
                H_assignment[7] += comment
                comment = None  # a trick to add no more comments to "assignment" at the end of this block of code
            
            elif nucleus[0] == 'C':   # only if the assignment line belongs to a Carbon then check if the user has made the same assignment
                comment = ""
                if spectrum_type in ['NOESY', 'NOESY (C-term hanging residue)', 'NOESY + TOCSY-NOESY', 'TOCSY-NOESY', 'TOCSY-NOESY + NOESY', 'TOCSY (unmatched) + TOCSY-NOESY', 'TOCSY-NOESY + TOCSY (unmatched)']:    # if this resonance comes only from NOESY
                    user_assignment = get_Carbon_resonance(resid, nucleus, 'NOESY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict,
                         aa_equivalentCarbons_dict, protein_alignment_list, absolute_matches_alignment_list, args)
                    if user_assignment == False:    # if no user assignment was made for this Carbon, add no extra comment
                        print "WARNING: no user assignment made in NOESY for resid", resid, "and Carbon", nucleus,". Checking TOCSY..."
                        user_assignment = get_Carbon_resonance(resid, nucleus, 'TOCSY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict,
                         aa_equivalentCarbons_dict, protein_alignment_list, absolute_matches_alignment_list, args)
                        if user_assignment == False:
                            print "WARNING: no user assignment made also in TOCSY for resid", resid, "and Carbon", nucleus,". Skipping nucleus."
                            assignment[7] += " <NOT ASSIGNED>"
                            continue
                    if approx_equal(assignment[1], user_assignment, 0.4): # if the user-made assignment agrees with ours, add to the comment "CORRECT"
                        print "DEBUG: CORRECT assignment found:", assignment
                        comment = " <CORRECT>"
                    else:
                        print "DEBUG: WRONG assignment found:", assignment
                        comment = " <WRONG>"
                elif spectrum_type in ['TOCSY (average)', 'TOCSY (unmatched)']:
                    user_assignment = get_Carbon_resonance(resid, nucleus, 'TOCSY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict,
                         aa_equivalentCarbons_dict, protein_alignment_list, absolute_matches_alignment_list, args)
                    if user_assignment == False:    # if no user assignment was made for this Carbon, add no extra comment
                        print "WARNING: no user assignment made in TOCSY for resid", resid, "and Carbon", nucleus,". Checking NOESY..."
                        user_assignment = get_Carbon_resonance(resid, nucleus, 'NOESY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict,
                         aa_equivalentCarbons_dict, protein_alignment_list, absolute_matches_alignment_list, args)
                        if user_assignment == False:
                            print "WARNING: no user assignment made also in NOESY for resid", resid, "and Carbon", nucleus,". Skipping nucleus."
                            assignment[7] += " <NOT ASSIGNED>"
                            continue
                    if approx_equal(assignment[1], user_assignment, 0.4): # if the user-made assignment agrees with ours, add to the comment "CORRECT"
                        print "DEBUG: CORRECT assignment found:", assignment
                        comment = " <CORRECT>"
                    else:
                        print "DEBUG: WRONG assignment found:", assignment
                        comment = " <WRONG>"
                
            if comment == None:      # this is a proton or a methyl Carbon or a Carbon with only one proton
                continue
            elif comment == "" and not '<' in assignment[7]:    # if no comment has been appended in previous iterations
                comment = " <NOT ASSIGNED>"
            assignment[7] += comment
    
    
    print "DEBUG: ITERATION 3 proofread_clean_resid_assignments_dict=", proofread_clean_resid_assignments_dict
    
    xeasy_out = open(args.ABSOLUTE_MATCHES_FILE+".NOESY.proofread.iter3.xeasy", 'w')
    atom_index = 0  # reinitialize atom index
    for resid in proofread_clean_resid_assignments_dict.keys():
        for peak in proofread_clean_resid_assignments_dict[resid]:
            atom_index += 1
            peak[0] = atom_index
            line = "\t".join([str(p) for p in peak])
            xeasy_out.write(line + "\n")
    xeasy_out.close()


print "DEBUG: ENTERING ITERATION 4 matched_NOESY_residue_assignments_dict=", matched_NOESY_residue_assignments_dict
print "DEBUG: ENTERING ITERATION 4 clean_resid_assignments_dict=", clean_resid_assignments_dict

xeasy_fout = open(args.ABSOLUTE_MATCHES_FILE+".NOESY.iter4.xeasy", 'w')
sparky_fout = open(args.ABSOLUTE_MATCHES_FILE+".NOESY.iter4.sparky", 'w')
atom_index = 1
residues_with_written_prevAssigned_peaks = []

for i_aa_resid in revordered_residue_keys:    # both the keys and the C,H assignments correspond to residue i-1
    i_residue = i_aa_resid[0]+str(i_aa_resid[1])
    aa_type = aa1to3_dict[i_aa_resid[0]]
    i_resid = get_resid_from_residue(i_residue)
    all_possible_assignments_list = []  # list with all the peaks for this residue (TOCSY-NOESY, NOESY) and the respective possible C-H assignments
    print "ITERATION 4: Matching residue i ", i_residue," TOCSY PEAKS to residue's i NOESY peaks..."
    checked_residues_set.add(i_residue)
    tmp_missing_carbons_list, tmp_partly_assigned_carbons_list = get_missing_carbons_from_xeasy_dict(aa_type, i_resid, clean_resid_assignments_dict, matched_NOESY_residue_assignments_dict)
    missing_carbons_list = [c for c in tmp_missing_carbons_list if not (aa_type in iter1_allowed_atoms_dict.keys() and c in iter1_allowed_atoms_dict[aa_type])]
    partly_assigned_carbons_list = [c for c in tmp_partly_assigned_carbons_list if not (aa_type in iter1_allowed_atoms_dict.keys() and c in iter1_allowed_atoms_dict[aa_type])]
    unassigned_NOESY_peaks_list = []    # ini
    if len(missing_carbons_list) != 0:   # IF THERE ARE MISSING CARBONS IN ITERATION 1, TRY TO FIND THEM IN NOESY (EXCLUDE THE ASSIGNED OR MATCHED NOESY PEAKS)
        print "Residue ", i_residue, " was not completely assigned in iteration 1. Using NOESY to assigned resonances to the rest Carbon types."
        try:    # DANGEROUS but works for the N-term residue (A200 in aLP)
            unassigned_NOESY_peaks_list = [peak for peak in matched_NOESY_residue_assignments_dict[i_residue] if peak[0]=="?" and peak[1]=="?" and peak[3]=="?"]
        except KeyError:    # if this residue is not in the alignment (not in matched_NOESY_residue_assignments_dict), then we cannot assign the partly assigned Carbons and hence we move on to the next
            print "Residue ", i_residue, " is either not in the alignment or doesn't have NOESY. Therefore it is impossible to assign partly assigned Carbon peaks."
            continue
        for NOESY_peak in unassigned_NOESY_peaks_list:  # make C-H type predictions for each unassigned NOESY peak
            NOESY_Creson = NOESY_peak[2]
            NOESY_Hreson = NOESY_peak[4]
            all_possible_assignments_list.extend(get_probabilities_from_H_C_resonpair_2Dhist(i_residue, NOESY_Hreson, NOESY_Creson, missing_carbons_list,
                                                                        partly_assigned_carbons_list, args, aa_CHpair_binProbabilityList_multidict, aa_carbon_binDensityList_multidict,
                                                                        NOESY_residue_peak_intensity_multidict, matched_i_to_iplus1_peaks_multidict, aa2pthres_dict,
                                                                        aa2ithres_dict, protein_alignment_list, absolute_matches_alignment_list, 
                                                                        aa_hydrogen_binDensityList_multidict, iteration=4))
        [assignment.append('NOESY') for assignment in all_possible_assignments_list]    # we don't need to created a new list, the original all_possible_assignments_list will be updated
        print "DEBUG: i_residue=", i_residue, "all_possible_assignments_list=", all_possible_assignments_list
        if i_resid in clean_resid_assignments_dict.keys():  # if this is a C-term residue that is analyzed for the first time in iter4
            written_nucleiNames_list = [line[3]  for line in clean_resid_assignments_dict[i_resid]]   # list of nuclei that have been written already, to avoid duplication
        else:
            written_nucleiNames_list = []
        if len(all_possible_assignments_list) > 0:  # IF WE HAVE FOUND POSSIBLE NOESY PEAK ASSIGNMENTS OF THE MISSING CARBONS
            written_nucleiNames_list.extend(select_best_NOESY_peak_combination(i_residue, all_possible_assignments_list, partly_assigned_carbons_list, atom_index, xeasy_fout, sparky_fout,
                                       protein_alignment_list, absolute_matches_alignment_list, args, NOESY_residue_NHresonances_dict,
                                       matched_NOESY_residue_assignments_dict, TOCSY_assigned_residue_list, iter1_allowed_atoms_dict, iteration=4))
        residues_with_written_prevAssigned_peaks.append(i_residue)
        write_peaks_assigned_in_previous_iteration(i_residue, clean_resid_assignments_dict, atom_index, xeasy_fout, 
                                               NOESY_residue_NHresonances_dict, residues_with_written_NH_list)
    elif len(missing_carbons_list) == 0 and len(partly_assigned_carbons_list) == 0: # NOTHING TO DO FOR THIS RESIDUE
        residues_with_written_prevAssigned_peaks.append(i_residue)
        write_peaks_assigned_in_previous_iteration(i_residue, clean_resid_assignments_dict, atom_index, xeasy_fout, 
                                               NOESY_residue_NHresonances_dict, residues_with_written_NH_list)
    else:   # IF ALL THE CARBONS WERE ASSIGNED IN ITERATION 1, SAVE THEM, BUT ALSO FIND IN NOESY AND SAVE PARTLY ASSIGNED CARBONS
        print "Residue ", i_residue, " was completely assigned in iteration 1."
        if not i_residue in matched_NOESY_residue_assignments_dict.keys():
            print "Residue ", i_residue, " is either not in the alignment or doesn't have NOESY. Therefore it is impossible to assign partly assigned Carbon peaks."
            continue
        try:
            aa_type = aa1to3_dict[get_aa_type_from_residue(i_residue)]
        except (IndexError, KeyError):    # skip residue i
            print "EXITING: COULD NOT FIND AA_TYPE OF RESIDUE", i_residue
            sys.exit(1)
        assigned_peaks_list = get_peaks_from_xeasy(i_resid, clean_resid_assignments_dict)     # get all assigned peaks of this residue from xeasy format
        for assigned_peak in assigned_peaks_list:
            spectrum_type = assigned_peak[5]
            Cname = assigned_peak[1]
            Creson = assigned_peak[2]    # use the TOCSY resonance by default
            Hname = assigned_peak[3]
            Hreson = assigned_peak[4]
            all_possible_assignments_list.append([aa_type, Cname, Hname, 10000000.0, Hreson, Creson, spectrum_type])  # set the probability to 1 cause it has been already assigned in iteration 1
            print "DEBUG: i_residue=", i_residue, "aa_type=", aa_type, "all_possible_assignments_list=", all_possible_assignments_list
        written_nucleiNames_list = []
        if i_resid in clean_resid_assignments_dict.keys():
            written_nucleiNames_list = [line[3]  for line in clean_resid_assignments_dict[i_resid]]   # list of nuclei that have been written already, to avoid duplication
        if len(all_possible_assignments_list) > 0:  # IF WE HAVE FOUND POSSIBLE NOESY PEAK ASSIGNMENTS OF THE MISSING CARBONS
            written_nucleiNames_list.extend(select_best_NOESY_peak_combination(i_residue, all_possible_assignments_list, partly_assigned_carbons_list, atom_index, xeasy_fout, sparky_fout,
                                       protein_alignment_list, absolute_matches_alignment_list, args, NOESY_residue_NHresonances_dict,
                                       matched_NOESY_residue_assignments_dict, TOCSY_assigned_residue_list, iter1_allowed_atoms_dict, iteration=4))


for i_aa_resid in revordered_residue_keys:
    i_residue = i_aa_resid[0]+str(i_aa_resid[1])
    aa_type = aa1to3_dict[i_aa_resid[0]]
    if not i_residue in residues_with_written_prevAssigned_peaks:
        write_peaks_assigned_in_previous_iteration(i_residue, clean_resid_assignments_dict, atom_index, xeasy_fout, 
                                               NOESY_residue_NHresonances_dict, residues_with_written_NH_list)


for residue in absolute_matches_alignment_list:
    if not residue in checked_residues_set and residue != '-':
        print "WARNING: Residue ", residue, " was not checked!!!"


write_NH_of_residues_without_CH(absolute_matches_alignment_list, NOESY_residue_NHresonances_dict, TOCSY_residue_NHresonances_dict,
                                    residues_with_written_NH_list, atom_index, xeasy_fout, sparky_fout, protein_alignment_list, args, 
                                    Cterm_residue_set, iteration=4)

xeasy_fout.close()
sparky_fout.close()


clean_resid_assignments_dict = clean_xeasy_file(args.ABSOLUTE_MATCHES_FILE+".NOESY.iter4.xeasy")  # dict with all the assignments from each residue from iteration 1
for resid, peak_list in clean_resid_assignments_dict.items():
    peak_list.sort(key=itemgetter(3))
    clean_resid_assignments_dict[resid] = peak_list

print "DEBUG: ITERATION 4 point 1 clean_resid_assignments_dict=", clean_resid_assignments_dict

xeasy_out = open(args.ABSOLUTE_MATCHES_FILE+".NOESY.iter4.xeasy", 'w')
atom_index = 0  # reinitialize atom index
for resid in clean_resid_assignments_dict.keys():
    for peak in clean_resid_assignments_dict[resid]:
        atom_index += 1
        peak[0] = atom_index
        line = "\t".join([str(p) for p in peak])
        xeasy_out.write(line + "\n")
xeasy_out.close()



clean_resid_assignments_dict = clean_xeasy_file(args.ABSOLUTE_MATCHES_FILE+".NOESY.iter4.xeasy")  
for resid, peak_list in clean_resid_assignments_dict.items():
    peak_list.sort(key=itemgetter(3))
    clean_resid_assignments_dict[resid] = peak_list

print "DEBUG: ITERATION 4 point 2 clean_resid_assignments_dict=", clean_resid_assignments_dict

xeasy_out = open(args.ABSOLUTE_MATCHES_FILE+".NOESY.cleaned.iter4.xeasy", 'w')
atom_index = 0  # reinitialize atom index
for resid in clean_resid_assignments_dict.keys():
    for peak in clean_resid_assignments_dict[resid]:
        atom_index += 1
        peak[0] = atom_index
        line = "\t".join([str(p) for p in peak])
        xeasy_out.write(line + "\n")
xeasy_out.close()


print "DEBUG: args.user_TOCSY_fname=", args.user_TOCSY_fname, "args.user_NOESY_fname=", args.user_NOESY_fname
if args.user_TOCSY_fname and args.user_NOESY_fname:
    print "STARTING PROOF-READING NOESY CHEMICAL SHIFT ASSIGNMENTS"
    for resid in clean_resid_assignments_dict:
        if not resid in proofread_clean_resid_assignments_dict.keys():
            proofread_clean_resid_assignments_dict[resid] = clean_resid_assignments_dict[resid]
            continue
        for new_assignment in clean_resid_assignments_dict[resid]:
            new_reson = new_assignment[1]
            old_reson_list = [a[1] for a in proofread_clean_resid_assignments_dict[resid]]
            if not new_reson in old_reson_list:     # if this resonance was not found (newly assigned), add it to the dictionary
                proofread_clean_resid_assignments_dict[resid].append(new_assignment)
    print "DEBUG: ITERATION 4 point 3 clean_resid_assignments_dict=", clean_resid_assignments_dict
    
    
    
    for resid in proofread_clean_resid_assignments_dict.keys():
        print "DEBUG ITERATION4: proof reading resid=", resid
        resname = proofread_clean_resid_assignments_dict[resid][0][6]
        if resname in aa_equivalentCarbons_dict.keys(): # If this is a LEU or VAL
            Carbon_pair = aa_equivalentCarbons_dict[resname][0]
            if len([a for a in proofread_clean_resid_assignments_dict[resid] if a[3] in Carbon_pair]) > 0: # AND if the methyl carbons have been assigned
                Proton_pair = ["H"+C[1:] for C in Carbon_pair]   # protons of the equivalent Carbons
                program_assignments_list = [assignment for assignment in proofread_clean_resid_assignments_dict[resid] if assignment[3] in Carbon_pair+Proton_pair]   # keep only the resonances assigned by the program to the equivalent methyl carbons and protons
                proof_read_equivalentMethyls(program_assignments_list, resid, user_VAL_LEU_NOESY_resid_assignments_dict, user_VAL_LEU_TOCSY_resid_assignments_dict, Carbon_pair, Proton_pair)   # it doesn't matter if 1 or 2 of the peaks have been proof-read in iteration 1 & 2
        for assignment in proofread_clean_resid_assignments_dict[resid]:
            if '<' in assignment[7]:    # if this assignment has been proofread, then skip it
                continue
            print "DEBUG ITERATION4: proof reading new assignment=", assignment
            spectrum_type = re.sub(r' iter[1-9]', '', assignment[7])    # remove the "iter[1-9]" suffix to be able to distinguish the type of assignment (NEESY, NOESY+TOCSY, etc.)
            comment = None
            aa_type = assignment[6]
            nucleus = assignment[3] # the nucleus name
            if aa_type in aa_equivalentCarbons_dict.keys() and nucleus in aa_equivalentCarbons_dict[aa_type][0]:    # if this peak has been assigned by the program to LEU or VAL CD1, CD2 or CG1,CG2 respectively, skip it
                continue
            
            if aa_type in aatype_carbon_nongeminalHname_multidict.keys() and nucleus in aatype_carbon_nongeminalHname_multidict[aa_type].keys():
                C_assignment = assignment # the program-made carbon assignment
                H_assignment = None # default value for the program-made proton assignment
                print "DEBUG: resid=", resid, "proofread_clean_resid_assignments_dict[resid]=", proofread_clean_resid_assignments_dict[resid]
                for a in proofread_clean_resid_assignments_dict[resid]:
                    if ((type(aatype_carbon_nongeminalHname_multidict[aa_type][nucleus]) == list and a[3] in aatype_carbon_nongeminalHname_multidict[aa_type][nucleus]) or
                        (type(aatype_carbon_nongeminalHname_multidict[aa_type][nucleus]) == str and a[3] == aatype_carbon_nongeminalHname_multidict[aa_type][nucleus])):
                        H_assignment = a
                        break
                print "DEBUG: C_assignment=", C_assignment, "H_assignment=", H_assignment
                comment = ""
                if spectrum_type in ['NOESY', 'NOESY (C-term hanging residue)', 'NOESY + TOCSY-NOESY', 'TOCSY-NOESY', 'TOCSY-NOESY + NOESY', 'TOCSY (unmatched) + TOCSY-NOESY', 'TOCSY-NOESY + TOCSY (unmatched)']:    # if this resonance comes only from NOESY
                    user_assignment = get_Carbon_resonance(resid, nucleus, 'NOESY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict,
                         aa_equivalentCarbons_dict, protein_alignment_list, absolute_matches_alignment_list, args, get_assignment=True)
                    if user_assignment == False:    # if no user assignment was made for this Carbon, add no extra comment
                        print "WARNING: no user assignment made in NOESY for resid", resid, "and Carbon", nucleus,". Checking TOCSY..."
                        user_assignment = get_Carbon_resonance(resid, nucleus, 'TOCSY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict,
                         aa_equivalentCarbons_dict, protein_alignment_list, absolute_matches_alignment_list, args, get_assignment=True)
                        if user_assignment == False:
                            print "WARNING: no user assignment made also in TOCSY for resid", resid, "and Carbon", nucleus,". Skipping nucleus."
                            assignment[7] += " <NOT ASSIGNED>"
                            continue
                    user_Creson = user_assignment[2]
                    user_Hreson = user_assignment[4]
                    print "DEBUG: user_Creson=", user_Creson, "user_Hreson=", user_Hreson
                    if approx_equal(C_assignment[1], user_Creson, 0.4) and approx_equal(H_assignment[1], user_Hreson, 0.04): # if the user-made assignment agrees with ours, add to the comment "CORRECT"
                        print "DEBUG: CORRECT C_assignment found:", C_assignment
                        print "DEBUG: CORRECT H_assignment found:", H_assignment
                        comment = " <CORRECT>"
                    else:
                        print "DEBUG: WRONG C_assignment found:", C_assignment
                        print "DEBUG: WRONG H_assignment found:", H_assignment
                        comment = " <WRONG>"
                elif spectrum_type in ['TOCSY (average)', 'TOCSY (unmatched)']:
                    user_assignment = get_Carbon_resonance(resid, nucleus, 'TOCSY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict,
                         aa_equivalentCarbons_dict, protein_alignment_list, absolute_matches_alignment_list, args, get_assignment=True)
                    if user_assignment == False:    # if no user assignment was made for this Carbon, add no extra comment
                        print "WARNING: no user assignment made in TOCSY for resid", resid, "and Carbon", nucleus,". Checking NOESY..."
                        user_assignment = get_Carbon_resonance(resid, nucleus, 'NOESY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict,
                         aa_equivalentCarbons_dict, protein_alignment_list, absolute_matches_alignment_list, args, get_assignment=True)
                        if user_assignment == False:
                            print "WARNING: no user assignment made also in NOESY for resid", resid, "and Carbon", nucleus,". Skipping nucleus."
                            assignment[7] += " <NOT ASSIGNED>"
                            continue
                    user_Creson = user_assignment[2]
                    user_Hreson = user_assignment[4]
                    if approx_equal(C_assignment[1], user_Creson, 0.4) and approx_equal(H_assignment[1], user_Hreson, 0.04): # if the user-made assignment agrees with ours, add to the comment "CORRECT"
                        print "DEBUG: CORRECT C_assignment found:", C_assignment
                        print "DEBUG: CORRECT H_assignment found:", H_assignment
                        comment = " <CORRECT>"
                    else:
                        print "DEBUG: WRONG C_assignment found:", C_assignment
                        print "DEBUG: WRONG H_assignment found:", H_assignment
                        comment = " <WRONG>"
                        
                if comment == "":
                    comment = " <NOT ASSIGNED>"
                C_assignment[7] += comment
                H_assignment[7] += comment
                comment = None  # a trick to add no more comments to "assignment" at the end of this block of code
            
            elif nucleus[0] == 'C':   # only if the assignment line belongs to a Carbon then check if the user has made the same assignment
                comment = ""
                if spectrum_type in ['NOESY', 'NOESY (C-term hanging residue)', 'NOESY + TOCSY-NOESY', 'TOCSY-NOESY', 'TOCSY-NOESY + NOESY', 'TOCSY (unmatched) + TOCSY-NOESY', 'TOCSY-NOESY + TOCSY (unmatched)']:    # if this resonance comes only from NOESY
                    user_assignment = get_Carbon_resonance(resid, nucleus, 'NOESY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict,
                         aa_equivalentCarbons_dict, protein_alignment_list, absolute_matches_alignment_list, args)
                    if user_assignment == False:    # if no user assignment was made for this Carbon, add no extra comment
                        print "WARNING: no user assignment made in NOESY for resid", resid, "and Carbon", nucleus,". Checking TOCSY..."
                        user_assignment = get_Carbon_resonance(resid, nucleus, 'TOCSY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict,
                         aa_equivalentCarbons_dict, protein_alignment_list, absolute_matches_alignment_list, args)
                        if user_assignment == False:
                            print "WARNING: no user assignment made also in TOCSY for resid", resid, "and Carbon", nucleus,". Skipping nucleus."
                            assignment[7] += " <NOT ASSIGNED>"
                            continue
                    if approx_equal(assignment[1], user_assignment, 0.4): # if the user-made assignment agrees with ours, add to the comment "CORRECT"
                        print "DEBUG: CORRECT assignment found:", assignment
                        comment = " <CORRECT>"
                    else:
                        print "DEBUG: WRONG assignment found:", assignment
                        comment = " <WRONG>"
                elif spectrum_type in ['TOCSY (average)', 'TOCSY (unmatched)']:
                    user_assignment = get_Carbon_resonance(resid, nucleus, 'TOCSY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict,
                         aa_equivalentCarbons_dict, protein_alignment_list, absolute_matches_alignment_list, args)
                    if user_assignment == False:    # if no user assignment was made for this Carbon, add no extra comment
                        print "WARNING: no user assignment made in TOCSY for resid", resid, "and Carbon", nucleus,". Checking NOESY..."
                        user_assignment = get_Carbon_resonance(resid, nucleus, 'NOESY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict,
                         aa_equivalentCarbons_dict, protein_alignment_list, absolute_matches_alignment_list, args)
                        if user_assignment == False:
                            print "WARNING: no user assignment made also in NOESY for resid", resid, "and Carbon", nucleus,". Skipping nucleus."
                            assignment[7] += " <NOT ASSIGNED>"
                            continue
                    if approx_equal(assignment[1], user_assignment, 0.4): # if the user-made assignment agrees with ours, add to the comment "CORRECT"
                        print "DEBUG: CORRECT assignment found:", assignment
                        comment = " <CORRECT>"
                    else:
                        print "DEBUG: WRONG assignment found:", assignment
                        comment = " <WRONG>"
                
            if comment == None:      # this is a proton or a methyl Carbon or a Carbon with only one proton
                continue
            elif comment == "" and not '<' in assignment[7]:    # if no comment has been appended in previous iterations
                comment = " <NOT ASSIGNED>"
            assignment[7] += comment
    
    
    print "DEBUG: ITERATION 4 proofread_clean_resid_assignments_dict=", proofread_clean_resid_assignments_dict
    
    xeasy_out = open(args.ABSOLUTE_MATCHES_FILE+".NOESY.proofread.iter4.xeasy", 'w')
    atom_index = 0  # reinitialize atom index
    for resid in proofread_clean_resid_assignments_dict.keys():
        for peak in proofread_clean_resid_assignments_dict[resid]:
            atom_index += 1
            peak[0] = atom_index
            line = "\t".join([str(p) for p in peak])
            xeasy_out.write(line + "\n")
    xeasy_out.close()


print "## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ITERATION 5 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##"
print "DEBUG: ENTERING ITERATION 5 matched_NOESY_residue_assignments_dict=", matched_NOESY_residue_assignments_dict
print "DEBUG: ENTERING ITERATION 5 clean_resid_assignments_dict=", clean_resid_assignments_dict

allowed_aa_atoms_dict["PHE"].extend(["HD1", "CD1"])
allowed_aa_atoms_dict["TYR"].extend(["HD1", "HE1", "CD1", "CE1"])


xeasy_fout = open(args.ABSOLUTE_MATCHES_FILE+".NOESY.iter5.xeasy", 'w')
sparky_fout = open(args.ABSOLUTE_MATCHES_FILE+".NOESY.iter5.sparky", 'w')
atom_index = 1
residues_with_written_prevAssigned_peaks = []


for i_aa_resid in revordered_residue_keys:    # both the keys and the C,H assignments correspond to residue i-1
    i_residue = i_aa_resid[0]+str(i_aa_resid[1])
    aa_type = aa1to3_dict[i_aa_resid[0]]
    i_resid = get_resid_from_residue(i_residue)
    if not aa_type in ['TYR', 'PHE', 'LYS', 'ARG']:
       continue
    all_possible_assignments_list = []  # list with all the peaks for this residue (TOCSY-NOESY, NOESY) and the respective possible C-H assignments
    print "ITERATION 5: Matching residue i ", i_residue," TOCSY PEAKS to residue's i NOESY peaks..."
    checked_residues_set.add(i_residue)
    tmp_missing_carbons_list, tmp_partly_assigned_carbons_list = get_missing_carbons_from_xeasy_dict(aa_type, i_resid, clean_resid_assignments_dict, matched_NOESY_residue_assignments_dict)
    missing_carbons_list = [c for c in tmp_missing_carbons_list if not (aa_type in iter1_allowed_atoms_dict.keys() and c in iter1_allowed_atoms_dict[aa_type])]
    partly_assigned_carbons_list = [c for c in tmp_partly_assigned_carbons_list if not (aa_type in iter1_allowed_atoms_dict.keys() and c in iter1_allowed_atoms_dict[aa_type])]
    unassigned_NOESY_peaks_list = []    # ini
    if len(missing_carbons_list) != 0:   # IF THERE ARE MISSING CARBONS IN ITERATION 1 & 2, TRY TO FIND THEM IN NOESY (EXCLUDE THE ASSIGNED OR MATCHED NOESY PEAKS)
        print "Residue ", i_residue, " was not completely assigned in iteration 1 & 2. Using NOESY to assigned resonances to the rest Carbon types."
        try:    # DANGEROUS but works for the N-term residue (A200 in aLP)
            unassigned_NOESY_peaks_list = [peak for peak in matched_NOESY_residue_assignments_dict[i_residue] if peak[0]=="?" and peak[1]=="?" and peak[3]=="?"]
        except KeyError:    # if this residue is not in the alignment (not in matched_NOESY_residue_assignments_dict), then we cannot assign the partly assigned Carbons and hence we move on to the next
            print "Residue ", i_residue, " is either not in the alignment or doesn't have NOESY. Therefore it is impossible to assign partly assigned Carbon peaks."
            continue
        for NOESY_peak in unassigned_NOESY_peaks_list:  # make C-H type predictions for each unassigned NOESY peak
            NOESY_Creson = NOESY_peak[2]
            NOESY_Hreson = NOESY_peak[4]
            print "DEBUG: aa_type",aa_type,"NOESY_Hreson",NOESY_Hreson, "NOESY_Creson", NOESY_Creson, "missing_carbons_list", missing_carbons_list
            all_possible_assignments_list.extend(get_probabilities_from_H_C_resonpair_2Dhist(i_residue, NOESY_Hreson, NOESY_Creson, missing_carbons_list,
                                                            partly_assigned_carbons_list, args, aa_CHpair_binProbabilityList_multidict, aa_carbon_binDensityList_multidict,
                                                            NOESY_residue_peak_intensity_multidict, matched_i_to_iplus1_peaks_multidict, aa2pthres_dict,
                                                            aa2ithres_dict, protein_alignment_list, absolute_matches_alignment_list, 
                                                            aa_hydrogen_binDensityList_multidict, iteration=5))
        [assignment.append('NOESY') for assignment in all_possible_assignments_list]    # we don't need to created a new list, the original all_possible_assignments_list will be updated
        print "DEBUG: i_residue=", i_residue, "all_possible_assignments_list=", all_possible_assignments_list
        written_nucleiNames_list = []
        if i_resid in clean_resid_assignments_dict.keys():
            written_nucleiNames_list = [line[3]  for line in clean_resid_assignments_dict[i_resid]] # list of nuclei that have been written already, to avoid duplication
        if len(all_possible_assignments_list) > 0:  # IF WE HAVE FOUND POSSIBLE NOESY PEAK ASSIGNMENTS OF THE MISSING CARBONS
            written_nucleiNames_list.extend(select_best_NOESY_peak_combination(i_residue, all_possible_assignments_list, partly_assigned_carbons_list, atom_index, xeasy_fout, sparky_fout,
                                       protein_alignment_list, absolute_matches_alignment_list, args, NOESY_residue_NHresonances_dict,
                                       matched_NOESY_residue_assignments_dict, TOCSY_assigned_residue_list, iter1_allowed_atoms_dict, iteration=5))
        residues_with_written_prevAssigned_peaks.append(i_residue)
        write_peaks_assigned_in_previous_iteration(i_residue, clean_resid_assignments_dict, atom_index, xeasy_fout, 
                                               NOESY_residue_NHresonances_dict, residues_with_written_NH_list)
    elif len(missing_carbons_list) == 0 and len(partly_assigned_carbons_list) == 0: # NOTHING TO DO FOR THIS RESIDUE
        residues_with_written_prevAssigned_peaks.append(i_residue)
        write_peaks_assigned_in_previous_iteration(i_residue, clean_resid_assignments_dict, atom_index, xeasy_fout, 
                                               NOESY_residue_NHresonances_dict, residues_with_written_NH_list)
    else:   # IF ALL THE CARBONS WERE ASSIGNED IN PREVIOUS ITERATIONS, SAVE THEM, BUT ALSO FIND IN NOESY AND SAVE PARTLY ASSIGNED CARBONS
        print "Residue ", i_residue, " was completely assigned in iteration 1 & 2."
        if not i_residue in matched_NOESY_residue_assignments_dict.keys():
            print "Residue ", i_residue, " is either not in the alignment or doesn't have NOESY. Therefore it is impossible to assign partly assigned Carbon peaks."
            continue
        try:
            aa_type = aa1to3_dict[get_aa_type_from_residue(i_residue)]
        except (IndexError, KeyError):    # skip residue i
            print "EXITING: COULD NOT FIND AA_TYPE OF RESIDUE", i_residue
            sys.exit(1)
        assigned_peaks_list = get_peaks_from_xeasy(i_resid, clean_resid_assignments_dict)     # get all assigned peaks of this residue from xeasy format
        for assigned_peak in assigned_peaks_list:
            spectrum_type = assigned_peak[5]
            Cname = assigned_peak[1]
            Creson = assigned_peak[2]    # use the TOCSY resonance by default
            Hname = assigned_peak[3]
            Hreson = assigned_peak[4]
            all_possible_assignments_list.append([aa_type, Cname, Hname, 10000000.0, Hreson, Creson, spectrum_type])  # set the probability to 1 cause it has been already assigned in iteration 1 & 2
            print "DEBUG: i_residue=", i_residue, "aa_type=", aa_type, "all_possible_assignments_list=", all_possible_assignments_list
        written_nucleiNames_list = []
        if i_resid in clean_resid_assignments_dict.keys():
            written_nucleiNames_list = [line[3]  for line in clean_resid_assignments_dict[i_resid]]  # list of nuclei that have been written already, to avoid duplication
        if len(all_possible_assignments_list) > 0:  # IF WE HAVE FOUND POSSIBLE NOESY PEAK ASSIGNMENTS OF THE MISSING CARBONS
            written_nucleiNames_list.extend(select_best_NOESY_peak_combination(i_residue, all_possible_assignments_list, partly_assigned_carbons_list, atom_index, xeasy_fout, sparky_fout,
                                       protein_alignment_list, absolute_matches_alignment_list, args, NOESY_residue_NHresonances_dict,
                                       matched_NOESY_residue_assignments_dict, TOCSY_assigned_residue_list, iter1_allowed_atoms_dict, iteration=5))
    

for i_aa_resid in revordered_residue_keys:
    i_residue = i_aa_resid[0]+str(i_aa_resid[1])
    aa_type = aa1to3_dict[i_aa_resid[0]]
    if not i_residue in residues_with_written_prevAssigned_peaks:
        write_peaks_assigned_in_previous_iteration(i_residue, clean_resid_assignments_dict, atom_index, xeasy_fout, 
                                               NOESY_residue_NHresonances_dict, residues_with_written_NH_list)


for residue in absolute_matches_alignment_list:
    if not residue in checked_residues_set and residue != '-':
        print "WARNING: Residue ", residue, " was not checked!!!"


write_NH_of_residues_without_CH(absolute_matches_alignment_list, NOESY_residue_NHresonances_dict, TOCSY_residue_NHresonances_dict,
                                    residues_with_written_NH_list, atom_index, xeasy_fout, sparky_fout, protein_alignment_list, args, 
                                    Cterm_residue_set, iteration=5)

xeasy_fout.close()
sparky_fout.close()


clean_resid_assignments_dict = clean_xeasy_file(args.ABSOLUTE_MATCHES_FILE+".NOESY.iter5.xeasy")  # dict with all the assignments from each residue from iteration 1 & 2
for resid, peak_list in clean_resid_assignments_dict.items():
    peak_list.sort(key=itemgetter(3))
    clean_resid_assignments_dict[resid] = peak_list

print "DEBUG: ITERATION 5 point 1 clean_resid_assignments_dict=", clean_resid_assignments_dict

xeasy_out = open(args.ABSOLUTE_MATCHES_FILE+".NOESY.iter5.xeasy", 'w')
atom_index = 0  # reinitialize atom index
for resid in clean_resid_assignments_dict.keys():
    for peak in clean_resid_assignments_dict[resid]:
        atom_index += 1
        peak[0] = atom_index
        line = "\t".join([str(p) for p in peak])
        xeasy_out.write(line + "\n")
xeasy_out.close()



clean_resid_assignments_dict = clean_xeasy_file(args.ABSOLUTE_MATCHES_FILE+".NOESY.iter5.xeasy")  
for resid, peak_list in clean_resid_assignments_dict.items():
    peak_list.sort(key=itemgetter(3))
    clean_resid_assignments_dict[resid] = peak_list

print "DEBUG: ITERATION 5 point 2 clean_resid_assignments_dict=", clean_resid_assignments_dict

xeasy_out = open(args.ABSOLUTE_MATCHES_FILE+".NOESY.cleaned.iter5.xeasy", 'w')
atom_index = 0  # reinitialize atom index
for resid in clean_resid_assignments_dict.keys():
    for peak in clean_resid_assignments_dict[resid]:
        atom_index += 1
        peak[0] = atom_index
        line = "\t".join([str(p) for p in peak])
        xeasy_out.write(line + "\n")
xeasy_out.close()


print "DEBUG: args.user_TOCSY_fname=", args.user_TOCSY_fname, "args.user_NOESY_fname=", args.user_NOESY_fname
if args.user_TOCSY_fname and args.user_NOESY_fname:
    print "STARTING PROOF-READING NOESY CHEMICAL SHIFT ASSIGNMENTS"
    aa_equivalentCarbons_dict["TYR"] = [["CD1"], ["CE1"]]
    aa_equivalentCarbons_dict["PHE"] = [["CD1"], ["CE1"]]
    for resid in clean_resid_assignments_dict:
        
        if not resid in proofread_clean_resid_assignments_dict.keys():
            proofread_clean_resid_assignments_dict[resid] = clean_resid_assignments_dict[resid]
            continue
        for new_assignment in clean_resid_assignments_dict[resid]:
            new_reson = new_assignment[1]
            old_reson_list = [a[1] for a in proofread_clean_resid_assignments_dict[resid]]
            if not new_reson in old_reson_list:     # if this resonance was not found (newly assigned), add it to the dictionary
                proofread_clean_resid_assignments_dict[resid].append(new_assignment)
    print "DEBUG: ITERATION 5 clean_resid_assignments_dict=", clean_resid_assignments_dict
    
    user_TYR_PHE_NOESY_resid_assignments_dict, user_TYR_PHE_TOCSY_resid_assignments_dict = create_equivalentCarbons_assignments_dict(user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict, protein_alignment_list,
                                              absolute_matches_alignment_list, aa_equivalentCarbons_dict, args)
    print "DEBUG: user_TYR_PHE_NOESY_resid_assignments_dict=", user_TYR_PHE_NOESY_resid_assignments_dict
    print "DEBUG: user_TYR_PHE_TOCSY_resid_assignments_dict=", user_TYR_PHE_TOCSY_resid_assignments_dict
    
    
    for resid in proofread_clean_resid_assignments_dict.keys():
        resname = proofread_clean_resid_assignments_dict[resid][0][6]
        if resname in aa_equivalentCarbons_dict.keys(): # If this is a TYR or PHE or LEU or VAL
            if len(aa_equivalentCarbons_dict[resname]) > 1:   # proof-read TYR and PHE and TYR and PHE
                for Carbon_pair in aa_equivalentCarbons_dict[resname]:
                    if len([a for a in proofread_clean_resid_assignments_dict[resid] if a[3] in Carbon_pair]) > 0: # AND if the methyl carbons have been assigned
                        Proton_pair = ["H"+C[1:] for C in Carbon_pair]   # protons of the equivalent Carbons
                        program_assignments_list = [assignment for assignment in proofread_clean_resid_assignments_dict[resid] if assignment[3] in Carbon_pair+Proton_pair]   # keep only the resonances assigned by the program to the equivalent methyl carbons
                        proof_read_equivalentMethyls(program_assignments_list, resid, user_TYR_PHE_NOESY_resid_assignments_dict, user_TYR_PHE_TOCSY_resid_assignments_dict, Carbon_pair, Proton_pair)   # it doesn't matter if 1 or 2 of the peaks have been proof-read in iteration 1 & 2
            else:   # proof-read LEU and VAL
                Carbon_pair = aa_equivalentCarbons_dict[resname][0]
                if len([a for a in proofread_clean_resid_assignments_dict[resid] if a[3] in Carbon_pair]) > 0: # AND if the methyl carbons have been assigned
                    Proton_pair = ["H"+C[1:] for C in Carbon_pair]   # protons of the equivalent Carbons
                    program_assignments_list = [assignment for assignment in proofread_clean_resid_assignments_dict[resid] if assignment[3] in Carbon_pair+Proton_pair]   # keep only the resonances assigned by the program to the equivalent methyl carbons and protons
                    proof_read_equivalentMethyls(program_assignments_list, resid, user_VAL_LEU_NOESY_resid_assignments_dict, user_VAL_LEU_TOCSY_resid_assignments_dict, Carbon_pair, Proton_pair)   # it doesn't matter if 1 or 2 of the peaks have been proof-read in iteration 1 & 2
        for assignment in proofread_clean_resid_assignments_dict[resid]:
            if '<' in assignment[7]:    # if this assignment has been proofread, then skip it
                continue
            print "DEBUG ITERATION5: proof reading new assignment=", assignment
            spectrum_type = re.sub(r' iter[1-9]', '', assignment[7])    # remove the "iter[1-9]" suffix to be able to distinguish the type of assignment (NEESY, NOESY+TOCSY, etc.)
            comment = None
            aa_type = assignment[6]
            nucleus = assignment[3] # the nucleus name
            if aa_type in aa_equivalentCarbons_dict.keys() and nucleus in aa_equivalentCarbons_dict[aa_type]:    # if this peak has been assigned by the program to LEU or VAL CD1, CD2 or CG1,CG2 respectively, skip it
                continue
            
            if aa_type in aatype_carbon_nongeminalHname_multidict.keys() and nucleus in aatype_carbon_nongeminalHname_multidict[aa_type].keys():
                C_assignment = assignment # the program-made carbon assignment
                H_assignment = None # default value for the program-made proton assignment
                print "DEBUG: resid=", resid, "proofread_clean_resid_assignments_dict[resid]=", proofread_clean_resid_assignments_dict[resid]
                for a in proofread_clean_resid_assignments_dict[resid]:
                    if ((type(aatype_carbon_nongeminalHname_multidict[aa_type][nucleus]) == list and a[3] in aatype_carbon_nongeminalHname_multidict[aa_type][nucleus]) or
                        (type(aatype_carbon_nongeminalHname_multidict[aa_type][nucleus]) == str and a[3] == aatype_carbon_nongeminalHname_multidict[aa_type][nucleus])):
                        H_assignment = a
                        break
                print "DEBUG: C_assignment=", C_assignment, "H_assignment=", H_assignment
                comment = ""
                if spectrum_type in ['NOESY', 'NOESY (C-term hanging residue)', 'NOESY + TOCSY-NOESY', 'TOCSY-NOESY', 'TOCSY-NOESY + NOESY', 'TOCSY (unmatched) + TOCSY-NOESY', 'TOCSY-NOESY + TOCSY (unmatched)']:    # if this resonance comes only from NOESY
                    user_assignment = get_Carbon_resonance(resid, nucleus, 'NOESY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict,
                         aa_equivalentCarbons_dict, protein_alignment_list, absolute_matches_alignment_list, args, get_assignment=True)
                    if user_assignment == False:    # if no user assignment was made for this Carbon, add no extra comment
                        print "WARNING: no user assignment made in NOESY for resid", resid, "and Carbon", nucleus,". Checking TOCSY..."
                        user_assignment = get_Carbon_resonance(resid, nucleus, 'TOCSY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict,
                         aa_equivalentCarbons_dict, protein_alignment_list, absolute_matches_alignment_list, args, get_assignment=True)
                        if user_assignment == False:
                            print "WARNING: no user assignment made also in TOCSY for resid", resid, "and Carbon", nucleus,". Skipping nucleus."
                            assignment[7] += " <NOT ASSIGNED>"
                            continue
                    user_Creson = user_assignment[2]
                    user_Hreson = user_assignment[4]
                    print "DEBUG: user_Creson=", user_Creson, "user_Hreson=", user_Hreson
                    if approx_equal(C_assignment[1], user_Creson, 0.4) and approx_equal(H_assignment[1], user_Hreson, 0.04): # if the user-made assignment agrees with ours, add to the comment "CORRECT"
                        print "DEBUG: CORRECT C_assignment found:", C_assignment
                        print "DEBUG: CORRECT H_assignment found:", H_assignment
                        comment = " <CORRECT>"
                    else:
                        print "DEBUG: WRONG C_assignment found:", C_assignment
                        print "DEBUG: WRONG H_assignment found:", H_assignment
                        comment = " <WRONG>"
                elif spectrum_type in ['TOCSY (average)', 'TOCSY (unmatched)']:
                    user_assignment = get_Carbon_resonance(resid, nucleus, 'TOCSY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict,
                         aa_equivalentCarbons_dict, protein_alignment_list, absolute_matches_alignment_list, args, get_assignment=True)
                    if user_assignment == False:    # if no user assignment was made for this Carbon, add no extra comment
                        print "WARNING: no user assignment made in TOCSY for resid", resid, "and Carbon", nucleus,". Checking NOESY..."
                        user_assignment = get_Carbon_resonance(resid, nucleus, 'NOESY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict,
                         aa_equivalentCarbons_dict, protein_alignment_list, absolute_matches_alignment_list, args, get_assignment=True)
                        if user_assignment == False:
                            print "WARNING: no user assignment made also in NOESY for resid", resid, "and Carbon", nucleus,". Skipping nucleus."
                            assignment[7] += " <NOT ASSIGNED>"
                            continue
                    user_Creson = user_assignment[2]
                    user_Hreson = user_assignment[4]
                    if approx_equal(C_assignment[1], user_Creson, 0.4) and approx_equal(H_assignment[1], user_Hreson, 0.04): # if the user-made assignment agrees with ours, add to the comment "CORRECT"
                        print "DEBUG: CORRECT C_assignment found:", C_assignment
                        print "DEBUG: CORRECT H_assignment found:", H_assignment
                        comment = " <CORRECT>"
                    else:
                        print "DEBUG: WRONG C_assignment found:", C_assignment
                        print "DEBUG: WRONG H_assignment found:", H_assignment
                        comment = " <WRONG>"
                        
                if comment == "":
                    comment = " <NOT ASSIGNED>"
                C_assignment[7] += comment
                H_assignment[7] += comment
                comment = None  # a trick to add no more comments to "assignment" at the end of this block of code
            
            elif nucleus[0] == 'C':   # only if the assignment line belongs to a Carbon then check if the user has made the same assignment
                comment = ""
                if spectrum_type in ['NOESY', 'NOESY (C-term hanging residue)', 'NOESY + TOCSY-NOESY', 'TOCSY-NOESY', 'TOCSY-NOESY + NOESY', 'TOCSY (unmatched) + TOCSY-NOESY', 'TOCSY-NOESY + TOCSY (unmatched)']:    # if this resonance comes only from NOESY
                    user_assignment = get_Carbon_resonance(resid, nucleus, 'NOESY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict,
                         aa_equivalentCarbons_dict, protein_alignment_list, absolute_matches_alignment_list, args)
                    if user_assignment == False:    # if no user assignment was made for this Carbon, add no extra comment
                        print "WARNING: no user assignment made in NOESY for resid", resid, "and Carbon", nucleus,". Checking TOCSY..."
                        user_assignment = get_Carbon_resonance(resid, nucleus, 'TOCSY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict,
                         aa_equivalentCarbons_dict, protein_alignment_list, absolute_matches_alignment_list, args)
                        if user_assignment == False:
                            print "WARNING: no user assignment made also in TOCSY for resid", resid, "and Carbon", nucleus,". Skipping nucleus."
                            assignment[7] += " <NOT ASSIGNED>"
                            continue
                    if approx_equal(assignment[1], user_assignment, 0.4): # if the user-made assignment agrees with ours, add to the comment "CORRECT"
                        print "DEBUG: CORRECT assignment found:", assignment
                        comment = " <CORRECT>"
                    else:
                        print "DEBUG: WRONG assignment found:", assignment
                        comment = " <WRONG>"
                elif spectrum_type in ['TOCSY (average)', 'TOCSY (unmatched)']:
                    user_assignment = get_Carbon_resonance(resid, nucleus, 'TOCSY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict,
                         aa_equivalentCarbons_dict, protein_alignment_list, absolute_matches_alignment_list, args)
                    if user_assignment == False:    # if no user assignment was made for this Carbon, add no extra comment
                        print "WARNING: no user assignment made in TOCSY for resid", resid, "and Carbon", nucleus,". Checking NOESY..."
                        user_assignment = get_Carbon_resonance(resid, nucleus, 'NOESY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict,
                         aa_equivalentCarbons_dict, protein_alignment_list, absolute_matches_alignment_list, args)
                        if user_assignment == False:
                            print "WARNING: no user assignment made also in NOESY for resid", resid, "and Carbon", nucleus,". Skipping nucleus."
                            assignment[7] += " <NOT ASSIGNED>"
                            continue
                    if approx_equal(assignment[1], user_assignment, 0.4): # if the user-made assignment agrees with ours, add to the comment "CORRECT"
                        print "DEBUG: CORRECT assignment found:", assignment
                        comment = " <CORRECT>"
                    else:
                        print "DEBUG: WRONG assignment found:", assignment
                        comment = " <WRONG>"
                
            if comment == None:      # this is a proton or a methyl Carbon or a Carbon with only one proton
                continue
            elif comment == "" and not '<' in assignment[7]:    # if no comment has been appended in previous iterations
                comment = " <NOT ASSIGNED>"
            assignment[7] += comment
    
    
    print "DEBUG: ITERATION 5 proofread_clean_resid_assignments_dict=", proofread_clean_resid_assignments_dict
    
    xeasy_out = open(args.ABSOLUTE_MATCHES_FILE+".NOESY.proofread.iter5.xeasy", 'w')
    atom_index = 0  # reinitialize atom index
    for resid in proofread_clean_resid_assignments_dict.keys():
        for peak in proofread_clean_resid_assignments_dict[resid]:
            atom_index += 1
            peak[0] = atom_index
            line = "\t".join([str(p) for p in peak])
            xeasy_out.write(line + "\n")
    xeasy_out.close()


print "## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FINAL ITERATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##"
print "DEBUG: ENTERING FINAL ITERATION matched_NOESY_residue_assignments_dict=", matched_NOESY_residue_assignments_dict
print "DEBUG: ENTERING FINAL ITERATION clean_resid_assignments_dict=", clean_resid_assignments_dict
    
aatype_Hname_Qname_multidict = tree()
aatype_Hname_Qname_multidict['ALA']['HB'] = 'QB'
aatype_Hname_Qname_multidict['ILE']['HG2'] = 'QG2'
aatype_Hname_Qname_multidict['ILE']['HD1'] = 'QD1'
aatype_Hname_Qname_multidict['LEU']['HD1'] = 'QD1'
aatype_Hname_Qname_multidict['LEU']['HD2'] = 'QD2'
aatype_Hname_Qname_multidict['THR']['HG2'] = 'QG2'
aatype_Hname_Qname_multidict['VAL']['HG1'] = 'QG1'
aatype_Hname_Qname_multidict['VAL']['HG2'] = 'QG2'
aatype_Hname_Qname_multidict['MET']['HE'] = 'QE'


sorted_resids = clean_resid_assignments_dict.keys()
sorted_resids.sort()

for resid in sorted_resids:
    for assignment in clean_resid_assignments_dict[resid]:
        resname = assignment[6]
        atom_name = assignment[3]
        if resname in aatype_Hname_Qname_multidict.keys() and atom_name in aatype_Hname_Qname_multidict[resname].keys():
            print "DEBUG FINAL ITERATION: renaming", assignment[3], "to", aatype_Hname_Qname_multidict[resname][atom_name]
            assignment[3] = aatype_Hname_Qname_multidict[resname][atom_name]

print "DEBUG: FINAL ITERATION clean_resid_assignments_dict=", clean_resid_assignments_dict
if args.OUT_fname:
    out_fname = args.OUT_fname
else:
    out_fname = "4DNOESY_assignedall"
xeasy_out = open(out_fname + ".xeasy", 'w')
atom_index = 0  # reinitialize atom index
for resid in sorted_resids:
    for peak in clean_resid_assignments_dict[resid]:
        atom_index += 1
        peak[0] = atom_index
        line = "\t".join([str(p) for p in peak])
        xeasy_out.write(line + "\n")
        if peak[6] in ["TYR", "PHE"] and peak[3] in ["CD1", "CE1", "HD1", "HE1"]:   # duplicate aromatic CD1-HD1 and CE1-HE1
            peak[3] = peak[3][:-1]+'2'
            atom_index += 1
            peak[0] = atom_index
            line = "\t".join([str(p) for p in peak])
            xeasy_out.write(line + "\n")
xeasy_out.close()


if args.user_TOCSY_fname and args.user_NOESY_fname:
    
    for resid in sorted_resids:
        for assignment in proofread_clean_resid_assignments_dict[resid]:
            resname = assignment[6]
            atom_name = assignment[3]
            if resname in aatype_Hname_Qname_multidict.keys() and atom_name in aatype_Hname_Qname_multidict[resname].keys():
                print "DEBUG FINAL ITERATION: renaming", assignment[3], "to", aatype_Hname_Qname_multidict[resname][atom_name]
                assignment[3] = aatype_Hname_Qname_multidict[resname][atom_name]
    
    print "DEBUG: FINAL ITERATION proofread_clean_resid_assignments_dict=", proofread_clean_resid_assignments_dict
    if args.OUT_fname:
        out_fname = args.OUT_fname
    else:
        out_fname = "4DNOESY_assignedall"
    xeasy_out = open(out_fname + ".proofread.xeasy", 'w')
    atom_index = 0  # reinitialize atom index
    for resid in sorted_resids:
        for peak in proofread_clean_resid_assignments_dict[resid]:
            atom_index += 1
            peak[0] = atom_index
            line = "\t".join([str(p) for p in peak])
            xeasy_out.write(line + "\n")
            if peak[6] in ["TYR", "PHE"] and peak[3] in ["CD1", "CE1", "HD1", "HE1"]:   # duplicate aromatic CD1-HD1 and CE1-HE1
                peak[3] = peak[3][:-1]+'2'
                atom_index += 1
                peak[0] = atom_index
                line = "\t".join([str(p) for p in peak])
                xeasy_out.write(line + "\n")
    xeasy_out.close()




annotate_NOESY_file(args.NOESY_fname, absolute_RIGmatches_alignment_list, absolute_matches_alignment_list,
                    matched_NOESY_residue_assignments_dict, NOESY_residue_peak_intensity_multidict, out_fname="only4DNOESY_assignedall.sparky")