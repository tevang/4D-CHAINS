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
import re
import sys
from argparse import SUPPRESS, ArgumentParser

import copy


# Allow CE-HE of MET in NOESY assignments (must not be allowed in TOCSY assignments)
from operator import itemgetter

from lib.alignment import Alignment
from lib.csa_for_noesy import read_NHassigned_spectrum_file, match_TOCSY_to_HCNH_lines, get_aa_type_from_residue, \
    match_HCNH_i_iplus1_peaks, update_HCNH_peaks, get_probabilities_from_H_C_resonpair_2Dhist, \
    select_best_HCNH_peak_combination, write_TOCSY_HCNH_matched_peaks, write_unmatched_TOCSY_peaks, \
    write_NH_of_residues_without_CH, clean_xeasy_file, is_valid_residue, create_equivalentCarbons_assignments_dict, \
    proof_read_equivalentMethyls, get_Carbon_resonance, get_missing_carbons_from_xeasy_dict, \
    write_peaks_assigned_in_previous_iteration, get_peaks_from_xeasy, remove_iplus1_labels, write_flanked_residues, \
    annotate_HCNH_file
from lib.global_func import tree, approx_equal, get_resid_from_residue
from lib.global_vars import allowed_aa_atoms_dict, aa1to3_dict, aatype_carbon_nondegenerateHlist_mdict, \
    aatype_carbon_nongeminalHname_mdict, exclude_from_i_iplus_matching_dict
from lib.probhist import ProbHist_Loader

allowed_aa_atoms_dict["MET"].extend(["HE", "CE"])

## Parse command line arguments
def cmdlineparse():
    parser = ArgumentParser(description="command line arguments",
        epilog="EXAMPLE: cs_assignment_NOESY.py -nhmap results_summary.3mers_round1_rst.chainlinkers.patch -rstart 200 -atocsy results_summary.3mers_round1_rst.chainlinkers.patch.sparky -noesy noesyHCNH.21.4.2016num.list -probprod -2dhist -wcthres 100 -bcthres 25 -percentile 0.9")
    parser.add_argument("-nhmap", dest="ABSOLUTE_MATCHES_FILE", required=False,
                        help="file with absolute matches from previous run",
                        metavar="<absolute matches file>")
    parser.add_argument("-rstart", dest="FIRST_RESIDUE_NUMBER", required=False, default=1, type=int, 
                        help="The number of the first residue in the protein sequence inside the NHmap_FILE. (default: %(default)s)",
                        metavar="<first residue number>")
    parser.add_argument("-atocsy", dest="TOCSY_fname", required=True,
                        help="4D TOCSY (HCTOCSYNH) file with atom assignments",
                        metavar="<4D TOCSY input file with atom assignments>")
    parser.add_argument("-noesy", dest="NOESY_FILE", required=True,
                        help="the 4D NOESY file (*num.list) produced by 4D_assignment_parallel.py",
                        metavar="<4D NOESY input file>")
    parser.add_argument("-o", dest="OUT_fname", required=False, default=None,
                        help="output file name (assigned 4D TOCSY in sparky format). (default: like -nhmap but with different extension)",
                        metavar="<output file>")
    parser.add_argument("-wH", dest="H_weight", required=False, type=float, default=0.1,
                        help="weight in (0.0-1.0] to apply on aa type prediction from aliphatic H resonances. (default: %(default)s)",
                        metavar="<H weight>")
    parser.add_argument("-wC", dest="C_weight", required=False, type=float, default=1.0,
                        help="weight in (0.0, 1.0] to apply on aa type prediction from aliphatic C resonances. (default: %(default)s)",
                        metavar="<C weight>")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s version 0.9')
    parser.add_argument("-probprod", dest="PROB_PRODUCT", required=False, action='store_true', default=False,
                        help="select the best C-H type assignment combination based on the product of probabilities of the individual \
                             C-H assignments. (default: %(default)s)")
    parser.add_argument("-2dhist", dest="USE_2D_HISTOGRAMS", required=False, action='store_true', default=False,
                        help="use 2D BMRB histograms for aa type prediction. (default: %(default)s)")
    parser.add_argument("-wcthres1", dest="WITHIN_CGROUP_THRESHOLD_ITER1", required=False, default=0.0, type=float,
                        help="The first Carbon type must have probability greater this ratio from the second one of the SAME CARBON GROUP in order \
                             to be kept in the final NOESY assignments. (default: %(default)s)")
    parser.add_argument("-bcthres1", dest="BETWEEN_CGROUP_THRESHOLD_ITER1", required=False, default=10.0, type=float,
                        help="The same Carbon type of a Carbon group must have probability greater this ratio from the probability of the \
                             SAME CARBON TYPE in any other Carbon group in order to be kept in the final NOESY assignments. (default: %(default)s)")
    parser.add_argument("-wcthres2", dest="WITHIN_CGROUP_THRESHOLD_ITER2", required=False, default=0.0, type=float,
                        help="same as -wcthres1, but for the 2nd iteration. (default: %(default)s)")
    parser.add_argument("-bcthres2", dest="BETWEEN_CGROUP_THRESHOLD_ITER2", required=False, default=10.0, type=float,
                        help="same as -bcthres1, but for the 2nd iteration. (default: %(default)s)")
    parser.add_argument("-wcthres3", dest="WITHIN_CGROUP_THRESHOLD_ITER3", required=False, default=0.0, type=float,
                        help="same as -wcthres1, but for the 3rd iteration. (default: %(default)s)")
    parser.add_argument("-bcthres3", dest="BETWEEN_CGROUP_THRESHOLD_ITER3", required=False, default=0.0, type=float,
                        help="same as -bcthres1, but for the 3rd iteration. (default: %(default)s)")
    parser.add_argument("-wcthres4", dest="WITHIN_CGROUP_THRESHOLD_ITER4", required=False, default=0.0, type=float,
                        help="same as -wcthres1, but for the 4th iteration. (default: %(default)s)")
    parser.add_argument("-bcthres4", dest="BETWEEN_CGROUP_THRESHOLD_ITER4", required=False, default=0.0, type=float,
                        help="same as -bcthres1, but for the 4th iteration. (default: %(default)s)")
    parser.add_argument("-wcthres5", dest="WITHIN_CGROUP_THRESHOLD_ITER5", required=False, default=0.0, type=float,
                        help="same as -wcthres1, but for the 5th iteration. (default: %(default)s)")
    parser.add_argument("-bcthres5", dest="BETWEEN_CGROUP_THRESHOLD_ITER5", required=False, default=0.0, type=float,
                        help="same as -bcthres1, but for the 5th iteration. (default: %(default)s)")
    parser.add_argument("-percentile1", dest="PERCENTILE_ITER1", required=False, default=0.85, type=float,
                        help="""This arguments is used as a theshold to discard low probability C-H type predictions in iteration 1.
                        Can be 0.9, 0.85, 0.8, 0.75, or 0.0 (no filtering). A 0.9 percentile of value X means that the top 10 percent of the
                        probability density values were above X. (default: %(default)s)""")
    parser.add_argument("-percentile2", dest="PERCENTILE_ITER2", required=False, default=0.9, type=float,
                        help="""This arguments is used as a theshold to discard low probability C-H type predictions in iteration 2.
                        Can be 0.9, 0.85, 0.8, 0.75, or 0.0 (no filtering). A 0.9 percentile of value X means that the top 10 percent of the
                        probability density values were above X. (default: %(default)s)""")
    parser.add_argument("-percentile3", dest="PERCENTILE_ITER3", required=False, default=0.8, type=float,
                        help="""This arguments is used as a theshold to discard low probability C-H type predictions in iteration 3.
                        Can be 0.9, 0.85, 0.8, 0.75, or 0.0 (no filtering). A 0.9 percentile of value X means that the top 10 percent of the
                        probability density values were above X. (default: %(default)s)""")
    parser.add_argument("-percentile4", dest="PERCENTILE_ITER4", required=False, default=0.8, type=float,
                        help="""This arguments is used as a theshold to discard low probability C-H type predictions in iteration 4.
                        Can be 0.9, 0.85, 0.8, 0.75, or 0.0 (no filtering). A 0.9 percentile of value X means that the top 10 percent of the
                        probability density values were above X. (default: %(default)s)""")
    parser.add_argument("-percentile5", dest="PERCENTILE_ITER5", required=False, default=0.8, type=float,
                        help="""This arguments is used as a theshold to discard low probability C-H type predictions in iteration 5.
                        Can be 0.9, 0.85, 0.8, 0.75, or 0.0 (no filtering). A 0.9 percentile of value X means that the top 10 percent of the
                        probability density values were above X. (default: %(default)s)""")
    parser.add_argument("-usertocsy", dest="user_TOCSY_FILE", required=False,
                        help="4D TOCSY (HCTOCSYNH) file in Sparky format with atom assignments made by the user. (default: %(default)s)",
                        metavar="<Sparky 4D TOCSY input file with user-made atom assignments>")
    parser.add_argument("-usernoesy", dest="user_NOESY_FILE", required=False,
                        help="4D NOESY (HCTOCSYNH) file in Sparky format with atom assignments made by the user. (default: %(default)s)",
                        metavar="<Sparky 4D NOESY input file with user-made atom assignments>")
    parser.add_argument("-int1", dest="USE_INTENSITIES_ITERATION1", required=False, action='store_true', default=False,
                        help="use peak intensities in iteration 1. (default: %(default)s)")
    parser.add_argument("-int2", dest="USE_INTENSITIES_ITERATION2", required=False, action='store_true', default=False,
                        help="use peak intensities in iteration 2. (default: %(default)s)")
    parser.add_argument("-int3", dest="USE_INTENSITIES_ITERATION3", required=False, action='store_true', default=False,
                        help="use peak intensities in iteration 3. (default: %(default)s)")
    parser.add_argument("-int4", dest="USE_INTENSITIES_ITERATION4", required=False, action='store_true', default=False,
                        help="use peak intensities in iteration 4. (default: %(default)s)")
    parser.add_argument("-int5", dest="USE_INTENSITIES_ITERATION5", required=False, action='store_true', default=False,
                        help="use peak intensities in iteration 5. (default: %(default)s)")
    parser.add_argument("-ithres1", dest="INTENSITY_THRESHOLD_ITERATION1", required=False, type=float, default=0.1,
                        help="use relative peak intensity threshold (default: %(default)s)")
    parser.add_argument("-ithres2", dest="INTENSITY_THRESHOLD_ITERATION2", required=False, type=float, default=0.0,
                        help="use relative peak intensity threshold (default: %(default)s)")
    parser.add_argument("-ithres3", dest="INTENSITY_THRESHOLD_ITERATION3", required=False, type=float, default=0.1,
                        help="use relative peak intensity threshold (default: %(default)s)")
    parser.add_argument("-ithres4", dest="INTENSITY_THRESHOLD_ITERATION4", required=False, type=float, default=0.0,
                        help="use relative peak intensity threshold (default: %(default)s)")
    parser.add_argument("-ithres5", dest="INTENSITY_THRESHOLD_ITERATION5", required=False, type=float, default=0.0,
                        help="use relative peak intensity threshold (default: %(default)s)")
    parser.add_argument("-itrans1", dest="INTENSITY_TRANSFORM_TYPE_ITER1", required=False, type=int, default=4,
                        help="transform the intensities in iteration1. Can be 1: do nothing; 2: (default: %(default)s)")
    parser.add_argument("-itrans2", dest="INTENSITY_TRANSFORM_TYPE_ITER2", required=False, type=int, default=2,
                        help="transform the intensities in iteration2. Can be 1: do nothing; 2: (default: %(default)s)")
    parser.add_argument("-itrans3", dest="INTENSITY_TRANSFORM_TYPE_ITER3", required=False, type=int, default=3,
                        help="transform the intensities in iteration3. Can be 1: do nothing; 2: (default: %(default)s)")
    parser.add_argument("-itrans4", dest="INTENSITY_TRANSFORM_TYPE_ITER4", required=False, type=int, default=2,
                        help="transform the intensities in iteration4. Can be 1: do nothing; 2: (default: %(default)s)")
    parser.add_argument("-itrans5", dest="INTENSITY_TRANSFORM_TYPE_ITER5", required=False, type=int, default=1,
                        help="transform the intensities in iteration5. Can be 1: do nothing; 2: (default: %(default)s)")
    parser.add_argument("-probmodel", dest="PROBABILITY_MODEL", required=False, type=int, default=2,
                        help="If '1' the probability of each peak will be given by [wH*1Dhist(H)+wC*1Dhist(C)]/(wH+wC). If '2' then by \
                             1Dhist(H)*1Dhist(C). (default: %(default)s)")
    parser.add_argument("-probmode", dest="CONSENSUS_CGROUP_PROB_MODE", required=False, default=5, type=int,
                        help="""The way to calculate the total score of a set of chemical shift assignment (default: %(default)s).
                    The following values control how to calculate the consensus probability of each C-group. The total score will be the
                    product of this consensus C-group probabilities.
                    1: average;
                    2: sqrt(prob1*prob2)    ; geometric mean
                    3: (prob2-prob1)/(log(prob2)-log(prob1))    ; logarithmic mean
                    4: (prob1*prob2)/((prob1+prob2)/2.) ;composite average 
                    5: prob1 + prob2    ;
                     . (default: %(default)s)
                        """, metavar="<way to calculate the cs assignment score>")
    
    # OPTIONAL ARGUMENTS
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
    
    # Alternatively you can define one ithres and pthres for all Methyls. If args.ithresMethyl and args.pthresMethyl != None then the above 10
    # values will be overriden!
    parser.add_argument("-ithresMethyl", dest="ithresMethyl", type=float, default=0.1, help=SUPPRESS) # intensity threshold for all Methyls
    parser.add_argument("-pthresMethyl", dest="pthresMethyl", type=float, default=0.8, help=SUPPRESS) # intensity threshold for all Methyls
    
    
    args=parser.parse_args()
    return args


args = cmdlineparse()
print("Input argument values:")
for arg in vars(args):
    print(arg, "=", getattr(args, arg))

if args.ithresMethyl:
    args.ithresA, args.ithresT, args.ithresV, args.ithresI, args.ithresL = args.ithresMethyl, args.ithresMethyl, args.ithresMethyl, args.ithresMethyl, args.ithresMethyl
if args.pthresMethyl:
    args.pthresA, args.pthresT, args.pthresV, args.pthresI, args.pthresL = args.pthresMethyl, args.pthresMethyl, args.pthresMethyl, args.pthresMethyl, args.pthresMethyl
# Make file paths absolute (not real, since you work with symlinks)
args.ABSOLUTE_MATCHES_FILE = os.path.abspath(args.ABSOLUTE_MATCHES_FILE)
print("DEBUG: args.ABSOLUTE_MATCHES_FILE=", args.ABSOLUTE_MATCHES_FILE)
args.TOCSY_fname = os.path.abspath(args.TOCSY_fname)
args.NOESY_FILE = os.path.abspath(args.NOESY_FILE)

##########################################################################################################
##                                         LOADING FILES                                                ##
##########################################################################################################

# LOAD PROBABILITY HISTOGRAMS
histload = ProbHist_Loader()
histload.load_1Dhistograms()
histload.load_2Dhistograms()

# Create all alignment files
ali = Alignment(args.ABSOLUTE_MATCHES_FILE, args.FIRST_RESIDUE_NUMBER, spectrum_combo="TOCSY-HCNH")
ali.create_absolute_AAIGmatches_alignment()
ali.create_absolute_matches_alignment()  # TODO: not sure if it should be "HCNH-HCNH"!

# TODO: test if the class above works as the code below and if yes remove the code below.
# ## READ FILE WITH ABSOLUTE MATCHES
# # absolute_AAIGmatches_alignment_list contains the absolutely matched root index groups (RIG), namely the names given to each N,HN peak in the HSQC.
# # It does not necessarily contain real residue names!
# protein_alignment_list, absolute_AAIGmatches_alignment_list  = read_NHmap(args.ABSOLUTE_MATCHES_FILE)
# # Remove terminal 'N/A' otherwise they will cause you trouble later
# NA_indices = []
# if protein_alignment_list[0] == 'N/A' or absolute_AAIGmatches_alignment_list[0] == 'N/A':
#     NA_indices.append(0)
# if protein_alignment_list[-1] == 'N/A' or absolute_AAIGmatches_alignment_list[-1] == 'N/A':
#     NA_indices.append(len(protein_alignment_list)-1)
# protein_alignment_list = [protein_alignment_list[i] for i in range(len(protein_alignment_list)) if not i in NA_indices]
# absolute_AAIGmatches_alignment_list = [absolute_AAIGmatches_alignment_list[i] for i in range(len(absolute_AAIGmatches_alignment_list)) if not i in NA_indices]
#
#
# resid = args.FIRST_RESIDUE_NUMBER        # start counting from the 1st residue number specified by the user
# absolute_matches_alignment_list = ['-']*len(absolute_AAIGmatches_alignment_list)     # contains the absolutely matched residue names, namely real residue names not AAIG signatures!
# for position in range(len(protein_alignment_list)):
#     resname = protein_alignment_list[position]
#     RIG = absolute_AAIGmatches_alignment_list[position]
#     if RIG == '-':
#         absolute_matches_alignment_list[position] = '-'
#     elif RIG == 'N/A':
#         absolute_matches_alignment_list[position] = 'N/A'
#     else:
#         residue = resname + str(resid+1)
#         absolute_matches_alignment_list[position] = residue
#     resid += 1

print("DEBUG: ali.protein_alignment_list=", ali.protein_alignment_list)
print("DEBUG: ali.absolute_AAIGmatches_alignment_list=", ali.absolute_AAIGmatches_alignment_list)
print("DEBUG: ali.absolute_matches_alignment_list=", ali.absolute_matches_alignment_list)


########################################################## END OF FUNCTION DEFINITIONS #########################################################    
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ITERATION 1 only for methyl groups ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

iter1_allowed_atoms_dict = {}   # resname->list of allowed carbons in iteration1
iter1_allowed_atoms_dict["LEU"] = ["CD1", "CD2", "CG"]
iter1_allowed_atoms_dict["VAL"] = ["CG1", "CG2"]
iter1_allowed_atoms_dict["ILE"] = ["CD1", "CG2"]
iter1_allowed_atoms_dict["THR"] = ["CG2"]
iter1_allowed_atoms_dict["ALA"] = ["CB"]
# Also methyl is Met CE-QE which must not be included in iteration 1 & 3

# Create dictionary of intensity threshold for each aa type
aa2ithres_dict = {}
if type(args.ithresA) == float:    aa2ithres_dict["ALA"] = args.ithresA
if type(args.ithresT) == float:    aa2ithres_dict["THR"] = args.ithresT
if type(args.ithresV) == float:    aa2ithres_dict["VAL"] = args.ithresV
if type(args.ithresI) == float:    aa2ithres_dict["ILE"] = args.ithresI
if type(args.ithresL) == float:    aa2ithres_dict["LEU"] = args.ithresL
# Create dictionary of percentile threshold for each aa type
aa2pthres_dict = {}
print("DEBUG: args.pthresA=", args.pthresA)
if type(args.pthresA) == float:    aa2pthres_dict["ALA"] = args.pthresA
if type(args.pthresT) == float:    aa2pthres_dict["THR"] = args.pthresT
if type(args.pthresV) == float:    aa2pthres_dict["VAL"] = args.pthresV
if type(args.pthresI) == float:    aa2pthres_dict["ILE"] = args.pthresI
if type(args.pthresL) == float:    aa2pthres_dict["LEU"] = args.pthresL


# LOAD ASSIGNED TOCSY CONTENTS
print("DEBUG: loading TOCSY file", args.TOCSY_fname)
TOCSY_residue_assignments_dict, \
TOCSY_residue_NHresonances_dict, \
i_to_iminus1_dict, \
patched_residues_list, \
TOCSY_residue_peak_intensity_mdict = \
    read_NHassigned_spectrum_file(args.TOCSY_fname, "TOCSY", ali)
print("DEBUG: patched_residues_list=", patched_residues_list)
print("DEBUG: i_to_iminus1_dict=", i_to_iminus1_dict)
print("DEBUG: TOCSY_residue_assignments_dict=", TOCSY_residue_assignments_dict)
print("DEBUG: TOCSY_residue_NHresonances_dict=", TOCSY_residue_NHresonances_dict)
for k,v in list(TOCSY_residue_assignments_dict.items()):
    print(k, "-->", v)
# Add the N-terminus residue in ali.absolute_matches_alignment_list, which does not have N-H
# try:    # if the N-H of the 2nd residue have been assigned
#     print "DEBUG: patching the N-terminus residue ", i_to_iminus1_dict[ali.absolute_matches_alignment_list[1]], " to ali.absolute_matches_alignment_list"
#     ali.absolute_matches_alignment_list[0] = i_to_iminus1_dict[ali.absolute_matches_alignment_list[1]]
# except IndexError:
#     pass
## TEMPORARILY DEACTIVATE THE UPDATE OF ali.absolute_matches_alignment_list TO SEE IF ALL THE C-TERMINAL RESIDUES ARE FOUND
# for i in i_to_iminus1_dict.keys():
#     if i in ali.absolute_matches_alignment_list:
#         i_index = ali.absolute_matches_alignment_list.index(i)
#         if i_index > 0 and ali.absolute_matches_alignment_list[i_index-1] in ['-', 'N/A']:
#             ali.absolute_matches_alignment_list[i_index-1] = i_to_iminus1_dict[i]
# print "DEBUG: after patching ali.absolute_matches_alignment_list=", ali.absolute_matches_alignment_list


# LOAD NOESY FILE CONTENTS
NOESY_residue_assignments_dict, \
NOESY_residue_NHresonances_dict, \
original_NOESY_peaks_dict, \
NOESY_residue_peak_intensity_mdict = read_NHassigned_spectrum_file(args.NOESY_FILE, 'HCNH', ali)

print("DEBUG: NOESY_residue_assignments_dict=", NOESY_residue_assignments_dict)
print("DEBUG: NOESY_residue_NHresonances_dict=", NOESY_residue_NHresonances_dict)
print("DEBUG: original_NOESY_peaks_dict=", original_NOESY_peaks_dict)
print("DEBUG: NOESY_residue_peak_intensity_mdict=", NOESY_residue_peak_intensity_mdict)
for k,v in list(NOESY_residue_assignments_dict.items()):
    print(k, "-->", v)

# MATCH ASSIGNED TOCSY PEAKS TO NOESY PEAKS
matched_NOESY_residue_assignments_dict, residue_unmatched_TOCSY_peaks_dict = match_TOCSY_to_HCNH_lines(TOCSY_residue_assignments_dict,
                                                        NOESY_residue_assignments_dict, ali.protein_alignment_list, ali.absolute_matches_alignment_list, args)
print("DEBUG iteration1 point1: matched_NOESY_residue_assignments_dict=", matched_NOESY_residue_assignments_dict)
for k,v in list(matched_NOESY_residue_assignments_dict.items()):
    print(k, "-->", v)
print("DEBUG: residue_unmatched_TOCSY_peaks_dict=", residue_unmatched_TOCSY_peaks_dict)

print("DEBUG: point 2 TOCSY_residue_assignments_dict=", TOCSY_residue_assignments_dict)

# # WRITE STATISTICS
# for residue in matched_NOESY_residue_assignments_dict.keys():
#     CD_peaks_list = [(p[2],p[4]) for p in matched_NOESY_residue_assignments_dict[residue] if p[0]==residue and p[1] in ['CD1', 'CD2']]
#     print "DEBUG: CD_peaks_list=", CD_peaks_list
#     for peak in NOESY_residue_peak_intensity_mdict[residue].keys():
#         print "DEBUG: peak=", peak
#         if peak in CD_peaks_list:
#             print "STATISTICS: Residue", residue," CD ", peak, "intensity",NOESY_residue_peak_intensity_mdict[residue][peak]
#     
# print "END OF STATISTICS"
# sys.exit(1)

# This is just to see which residues have been checked, not necessarily assigned 
checked_residues_set = set()
# This is to remember the C-terminal residues
Cterm_residue_set = set()


## FIND MISSING CARBON ASSIGNMENTS AND TRY TO PREDICT THEM FROM NOESY
# ATTENTION: recall that in TOCSY_residue_assignments_dict both the keys and the C,H assignments correspond to residue i-1.
# But in matched_noesy_residue_assignments_dict the keys are residue i and the values are peaks from residues i (mostly the strongest), residue i-1 and maybe i+1 or
# other neighboring residues.
# So in order to match the NOESY peaks of residue i we need primarily the values of TOCSY_residue_assignments_dict keys i respectively.
xeasy_fout = open(args.ABSOLUTE_MATCHES_FILE+".NOESY.iter1.xeasy", 'w')
sparky_fout = open(args.ABSOLUTE_MATCHES_FILE+".NOESY.iter1.sparky", 'w')
atom_index = 1
TOCSY_assigned_residue_list = [] # list of the residues that have their peaks already assigned using TOCSY & NOESY or only TOCSY
residues_with_written_NH_list = [] # list of the residues for which the N & HN resonances have been written in the xeasy file
# start the peak assignment from the last residue to the first
revordered_residue_keys_iter1 = [(k[0], int(k[1:])) for k in list(TOCSY_residue_assignments_dict.keys())]
revordered_residue_keys_iter1.sort(key=itemgetter(1), reverse=True)  # reverse ordered keys, namely they start from the last residue in the protein sequence
matched_i_to_iplus1_peaks_mdict = tree()    # residue -> i peak -> matched i+1 peak
residues_with_written_TOCSY_peaks = []

selected_revordered_residue_keys_iter1 = []
for aa_type  in ["ALA", "THR", "VAL", "ILE", "LEU"]:
    for i_aa_resid in revordered_residue_keys_iter1:
        i_residue = i_aa_resid[0]+str(i_aa_resid[1])
        aa = aa1to3_dict[get_aa_type_from_residue(i_residue)]
        if aa == aa_type:
            selected_revordered_residue_keys_iter1.append(i_aa_resid)
print("DEBUG: selected_revordered_residue_keys_iter1=", selected_revordered_residue_keys_iter1)
# selected_revordered_residue_keys_iter1 = revordered_residue_keys_iter1


# # Measure i->i+1 statistics
# revordered_residue_keys_iter1 = [(k[0], int(k[1:])) for k in matched_NOESY_residue_assignments_dict.keys() if k in ali.absolute_matches_alignment_list]
# revordered_residue_keys_iter1.sort(key=itemgetter(1))
# N_total = 0
# N_found = 0
# for i in range(len(revordered_residue_keys_iter1)-1):
#     if revordered_residue_keys_iter1[i+1][1] != revordered_residue_keys_iter1[i][1]+1:
#         continue
#     i_res = revordered_residue_keys_iter1[i][0]+str(revordered_residue_keys_iter1[i][1])
#     iplus1_res = revordered_residue_keys_iter1[i+1][0]+str(revordered_residue_keys_iter1[i+1][1])
#     N_total += len([p for p in matched_NOESY_residue_assignments_dict[i_res] if p[0]==i_res])
#     N_found += len([p for p in matched_NOESY_residue_assignments_dict[iplus1_res] if p[0]==i_res])
# print "total=",N_total
# print "i->i+1 matched=", N_found



for i_aa_resid in selected_revordered_residue_keys_iter1:    # both the keys and the C,H assignments correspond to residue i-1
    i_residue = i_aa_resid[0]+str(i_aa_resid[1])
    # ONLY FOR ITERATION 1
    aa_type = aa1to3_dict[get_aa_type_from_residue(i_residue)]
    if not aa_type in list(iter1_allowed_atoms_dict.keys()):
        continue
    if aa_type in list(aa2ithres_dict.keys()):
        istart = int(10 * aa2ithres_dict[aa_type])
    else:
        istart = 5
    print("DEBUG: aa_type=", aa_type, "istart=", istart)
    for I_THRESHOLD in reversed(list(range(1, istart+1, 1))):
        I_THRESHOLD *= 0.1
        # if not i_residue in ['L394']:
        #    continue
        # FIRST MATCH THE TOCSY ASSIGNMENTS OF RESIDUE i TO NOESY PEAKS
        #i_residue = get_i_residue_from_iminus1_residue(iminus1_residue)
        all_possible_assignments_list = []  # list with all the peaks for this residue (TOCSY-HCNH, NOESY) and the respective possible C-H assignments
        if i_residue in list(TOCSY_residue_assignments_dict.keys()) and len(TOCSY_residue_assignments_dict[i_residue]) > 0:   # if no TOCSY assignments have been made for residue i, skip NOESY assignment of residue i
            print("ITERATION 1: Matching residue i ", i_residue," TOCSY PEAKS to residue's i NOESY peaks..., I_THRESHOLD=", I_THRESHOLD)
            checked_residues_set.add(i_residue)
            all_allowed_carbons_list = [x for x in allowed_aa_atoms_dict[aa_type] if x[0] == "C"]   # all allowed carbons for this aa type
            print("DEBUG: all_allowed_carbons_list=", all_allowed_carbons_list)
            print("DEBUG: i_residue=", i_residue, "TOCSY_residue_assignments_dict[i_residue]=", TOCSY_residue_assignments_dict[i_residue])
            tmp_assigned_carbons_list = [x[1] for x in TOCSY_residue_assignments_dict[i_residue]]   # these are the assignments of residue i
            # Remove Carbon names that belong to CH2 and had been found only once in TOCSY (namely 1 C-H peak is missing and must be found in NOESY)
            assigned_carbons_list = []  # list of carbons which have ALL their protons assigned in TOCSY
            partly_assigned_carbons_list = []   # list of methylene carbons for which only one proton frequency was found in TOCSY
            for Cname in tmp_assigned_carbons_list:
                if Cname in list(aatype_carbon_nondegenerateHlist_mdict[aa_type].keys()) and tmp_assigned_carbons_list.count(Cname) == 2:
                    assigned_carbons_list.append(Cname)
                if Cname in list(aatype_carbon_nondegenerateHlist_mdict[aa_type].keys()) and tmp_assigned_carbons_list.count(Cname) == 1:
                    partly_assigned_carbons_list.append(Cname)
                elif not Cname in list(aatype_carbon_nondegenerateHlist_mdict[aa_type].keys()):
                    assigned_carbons_list.append(Cname)
            print("DEBUG: assigned_carbons_list=", assigned_carbons_list)
            print("DEBUG: partly_assigned_carbons_list=", partly_assigned_carbons_list)
            tmp_missing_carbons_list = []
            for carbon in all_allowed_carbons_list:     # keep only the carbon types of this aa type that were not assigned to peaks
                if not carbon in assigned_carbons_list:
                    tmp_missing_carbons_list.append(carbon)
            print("DEBUG: dict key i_residue", i_residue, "tmp_missing_carbons_list=", tmp_missing_carbons_list)
            matched_iplus1, \
            matched_NOESY_residue_assignments_dict, \
            matched_i_to_iplus1_peaks_mdict = \
                match_HCNH_i_iplus1_peaks(i_residue,
                                          ali.absolute_matches_alignment_list,
                                          matched_NOESY_residue_assignments_dict,
                                          matched_i_to_iplus1_peaks_mdict)  # add labels to the i->i+1 matched peaks. If i_residue is a C-term residue or presends a gap, skip it from this iteration
            if matched_iplus1 == False:
                continue
            else:   # otherwise proceed to assigned but not from the extreme carbons in exclude_from_i_iplus_matching_dict
                missing_carbons_list = [c for c in tmp_missing_carbons_list if c in iter1_allowed_atoms_dict[aa_type]]
                print("DEBUG: dict key i_residue", i_residue, "missing_carbons_list=", missing_carbons_list)
                # ONLY FOR ITERATION 1
                if len(missing_carbons_list) == 0:
                    matched_NOESY_residue_assignments_dict = update_HCNH_peaks(i_residue, matched_NOESY_residue_assignments_dict, ali.protein_alignment_list, ali.absolute_matches_alignment_list, args, sparky_lines_list=[])   # remove "i+1" and "i-1" labels
                    continue
            unassigned_NOESY_peaks_list = []
            if len(missing_carbons_list) != 0:   # IF THERE ARE MISSING CARBONS IN TOCSY, TRY TO FIND THEM IN NOESY (EXCLUDE THE ASSIGNED OR MATCHED NOESY PEAKS)
                print("Residue ", i_residue, " was not completely assigned in TOCSY. Using NOESY to assigned resonances to the rest Carbon types.")
                try:    # DANGEROUS but works for the N-term residue (A200 in aLP)
                    unassigned_NOESY_peaks_list = [peak for peak in matched_NOESY_residue_assignments_dict[i_residue] if peak[0]=="i+1" and peak[1]=="?" and peak[3]=="?"]
                except KeyError:
                    continue
                print("DEBUG: unassigned_NOESY_peaks_list=", unassigned_NOESY_peaks_list)
                for NOESY_peak in unassigned_NOESY_peaks_list:  # make C-H type predictions for each unassigned NOESY peak
                    NOESY_Creson = NOESY_peak[2]
                    NOESY_Hreson = NOESY_peak[4]
                    all_possible_assignments_list.extend( get_probabilities_from_H_C_resonpair_2Dhist(i_residue,
                                                                                                      NOESY_Hreson,
                                                                                                      NOESY_Creson,
                                                                                                      missing_carbons_list,
                                                                                                      partly_assigned_carbons_list,
                                                                                                      args,
                                                                                                      NOESY_residue_peak_intensity_mdict,
                                                                                                      matched_i_to_iplus1_peaks_mdict,
                                                                                                      aa2pthres_dict,
                                                                                                      aa2ithres_dict,
                                                                                                      ali.protein_alignment_list,
                                                                                                      ali.absolute_matches_alignment_list,
                                                                                                      histload,
                                                                                                      iteration=1,
                                                                                                      I_THRESHOLD=I_THRESHOLD)
                                                          )
                [assignment.append('HCNH') for assignment in all_possible_assignments_list]    # we don't need to created a new list, the original all_possible_assignments_list will be updated
                print("DEBUG: i_residue=", i_residue, "all_possible_assignments_list=", all_possible_assignments_list)
                written_nucleiNames_list = []   # list of nuclei that have been written already, to avoid duplication
                if len(all_possible_assignments_list) > 0:  # IF WE HAVE FOUND POSSIBLE NOESY PEAK ASSIGNMENTS OF THE MISSING CARBONS
                    # WRITE THE BEST NOESY PEAK COMBINATION TO COMPLETE THE ASSIGNMENT OF THE MISSING CARBONS
                    written_nucleiNames_list, \
                    matched_NOESY_residue_assignments_dict, \
                    TOCSY_assigned_residue_list = \
                        select_best_HCNH_peak_combination(i_residue,
                                                          all_possible_assignments_list,
                                                          partly_assigned_carbons_list,
                                                          atom_index,
                                                          xeasy_fout,
                                                          sparky_fout,
                                                          ali.protein_alignment_list,
                                                          ali.absolute_matches_alignment_list,
                                                          args,
                                                          NOESY_residue_NHresonances_dict,
                                                          matched_NOESY_residue_assignments_dict,
                                                          TOCSY_assigned_residue_list,
                                                          iter1_allowed_atoms_dict,
                                                          iteration=1)
                # WRITE ALSO THE TOCSY-HCNH MATCHED PEAKS OF THIS RESIDUE
                residues_with_written_TOCSY_peaks.append(i_residue)
                if write_TOCSY_HCNH_matched_peaks(i_residue, atom_index, matched_NOESY_residue_assignments_dict, xeasy_fout, sparky_fout,
                                    ali.protein_alignment_list, ali.absolute_matches_alignment_list, args, TOCSY_assigned_residue_list,
                                    original_NOESY_peaks_dict, residue_unmatched_TOCSY_peaks_dict, NOESY_residue_NHresonances_dict,
                                    residues_with_written_NH_list, written_nucleiNames_list) == False:
                    continue
            else:   # IF ALL THE CARBONS WERE ASSIGNED IN TOCSY, FIND & SAVE PARTLY ASSIGNED CARBONS IN NOESY,
                    # SAVE THE MATCHED NOESY RESONANCES, AND IF NOT MATCHED, SAVE THE TOCSY RESONANCES
                print("Residue ", i_residue, " was completely assigned in TOCSY.")
                try:
                    aa_type = aa1to3_dict[get_aa_type_from_residue(i_residue)]
                except (IndexError, KeyError):    # skip residue i
                    print("EXITING: COULD NOT FIND AA_TYPE OF RESIDUE", i_residue)
                    sys.exit(1)
                for TOCSY_assignment in TOCSY_residue_assignments_dict[i_residue]:
                    source_spectrum = 'TOCSY (unmatched)'
                    Cname = TOCSY_assignment[1]
                    Creson = TOCSY_assignment[2]    # use the TOCSY resonance by default
                    Hname = TOCSY_assignment[3]
                    Hreson = TOCSY_assignment[4]
                    if i_residue in list(matched_NOESY_residue_assignments_dict.keys()):  # first check if this residue has NOESY peaks
                        matched_NOESY_peak = [peak for peak in matched_NOESY_residue_assignments_dict[i_residue] if peak[0]==i_residue and peak[1]==Cname and peak[3]==Hname]
                        if len(matched_NOESY_peak) == 1:    # if this TOCSY peak has been matched with a NOESY peak, use the NOESY resonances (more accurate)
                            Creson = matched_NOESY_peak[0][2]
                            Hreson = matched_NOESY_peak[0][4]
                            source_spectrum = 'TOCSY-HCNH'
                        elif len(matched_NOESY_peak) > 1:
                            print("ERROR: residue's ", i_residue, "TOCSY peak ", TOCSY_assignment, ", has been matched with multiple NOESY peaks", matched_NOESY_peak)
                            sys.exit(1)
                    all_possible_assignments_list.append([aa_type, Cname, Hname, 10000000.0, Hreson, Creson, source_spectrum])  # set the probability to 1 cause we don't know it
                    print("DEBUG: i_residue=", i_residue, "aa_type=", aa_type, "all_possible_assignments_list=", all_possible_assignments_list)
                written_nucleiNames_list = []   # list of nuclei that have been written already, to avoid duplication
                if len(all_possible_assignments_list) > 0:  # IF WE HAVE FOUND POSSIBLE NOESY PEAK ASSIGNMENTS OF THE MISSING CARBONS
                    # WRITE THE BEST NOESY PEAK COMBINATION TO COMPLETE THE ASSIGNMENT OF THE MISSING CARBONS
                    written_nucleiNames_list, \
                    matched_NOESY_residue_assignments_dict, \
                    TOCSY_assigned_residue_list = \
                        select_best_HCNH_peak_combination(i_residue,
                                                           all_possible_assignments_list,
                                                           partly_assigned_carbons_list,
                                                           atom_index,
                                                           xeasy_fout,
                                                           sparky_fout,
                                                           ali.protein_alignment_list,
                                                           ali.absolute_matches_alignment_list,
                                                           args,
                                                           NOESY_residue_NHresonances_dict,
                                                           matched_NOESY_residue_assignments_dict,
                                                           TOCSY_assigned_residue_list,
                                                           iter1_allowed_atoms_dict,
                                                           iteration=1)

# WRITE ALSO THE TOCSY-HCNH MATCHED PEAKS OF ALL RESIDUES THAT WERE EXCLUDED FROM THIS ITERATION
for i_aa_resid in revordered_residue_keys_iter1:
    i_residue = i_aa_resid[0]+str(i_aa_resid[1])
    aa_type = aa1to3_dict[i_aa_resid[0]]
    if not i_residue in residues_with_written_TOCSY_peaks:
        if write_TOCSY_HCNH_matched_peaks(i_residue, atom_index, matched_NOESY_residue_assignments_dict, xeasy_fout, sparky_fout,
                                    ali.protein_alignment_list, ali.absolute_matches_alignment_list, args, TOCSY_assigned_residue_list,
                                    original_NOESY_peaks_dict, residue_unmatched_TOCSY_peaks_dict, NOESY_residue_NHresonances_dict,
                                    residues_with_written_NH_list) == False:
            continue

print("DEBUG iteration1 point2: matched_NOESY_residue_assignments_dict=", matched_NOESY_residue_assignments_dict)

# Now find if there are any residue in the alignment that haven't been checked:
for residue in ali.absolute_matches_alignment_list:
    if not residue in checked_residues_set and residue != '-':
        print("WARNING: Residue ", residue, " was not checked!!!")
        

# ##
# ##  COPY FROM TOCSY ASSIGNMENTS (TOCSY_residue_assignments_dict) OF ALL THE RESONANCES THAT COULD NOT BE MATCHED WITH THE NOESY ASSIGNMENTS (TOCSY_residue_assignments_dict)
# ##
write_unmatched_TOCSY_peaks(matched_NOESY_residue_assignments_dict, TOCSY_residue_assignments_dict, atom_index, residues_with_written_NH_list,
                                xeasy_fout, sparky_fout, ali.protein_alignment_list, ali.absolute_matches_alignment_list, args,
                                NOESY_residue_NHresonances_dict, TOCSY_residue_NHresonances_dict, iteration=1)


# ##
# ##  WRITE THE HN & N RESONANCES TO XEASY ONLY FOR THOSE RESIDUES THAT NO C,H ASSIGNMENTS COULD BE MADE BUT THEIR HN & H RESONANCES ARE KNOWN FROM TOCSY OR NOESY
# ##
write_NH_of_residues_without_CH(ali.absolute_matches_alignment_list, NOESY_residue_NHresonances_dict, TOCSY_residue_NHresonances_dict,
                                    residues_with_written_NH_list, atom_index, xeasy_fout, sparky_fout, ali.protein_alignment_list, args,
                                    Cterm_residue_set, iteration=1)

xeasy_fout.close()
sparky_fout.close()

##
##  REWRITE THE XEASY FILE BY REMOVING DUPLICATES OR AVERAGING, AND CORRECTING THE ORDER
##

## CLEAN THE CONTENTS OF XEASY FILE
# clean_resid_assignments_dict contains the final, non-redundant assignments written to xeasy file. E.g.
# 1: [[2336, 56.022, 0.2, 'CA', 1, '#', 'MET', 'TOCSY (unmatched) (iter1)'], [2338, 34.786, 0.2, 'CB', 1, '#', 'MET', 'TOCSY (unmatched) (iter1)'],
# [2340, 30.946, 0.2, 'CG', 1, '#', 'MET', 'TOCSY (unmatched) (iter1)'], [2337, 4.194, 0.02, 'HA', 1, '#', 'MET', 'TOCSY (unmatched) (iter1)'],
# [2339, 2.14, 0.02, 'QB', 1, '#', 'MET', 'TOCSY (unmatched) (iter1)'], [2341, 2.487, 0.02, 'QG', 1, '#', 'MET', 'TOCSY (unmatched) (iter1)']]
clean_resid_assignments_dict = clean_xeasy_file(args.ABSOLUTE_MATCHES_FILE+".NOESY.iter1.xeasy")  
# sort the peak assignments of each resid
for resid, peak_list in list(clean_resid_assignments_dict.items()):
    peak_list.sort(key=itemgetter(3))
    clean_resid_assignments_dict[resid] = peak_list

print("DEBUG: ITERATION 1 point 1 clean_resid_assignments_dict=", clean_resid_assignments_dict)
## WRITE THE CLEAN XEASY FILE
xeasy_out = open(args.ABSOLUTE_MATCHES_FILE+".NOESY.cleaned.iter1.xeasy", 'w')
atom_index = 0  # reinitialize atom index
for resid in list(clean_resid_assignments_dict.keys()):
    for peak in clean_resid_assignments_dict[resid]:
        atom_index += 1
        peak[0] = atom_index
        line = "\t".join([str(p) for p in peak])
        xeasy_out.write(line + "\n")
xeasy_out.close()


##
##  IF PROVIDED, LOAD THE USER-ANNOTATED TOCSY AND NOESY FILES FOR PROOF-READING
##
print("DEBUG: args.user_TOCSY_FILE=", args.user_TOCSY_FILE, "args.user_NOESY_FILE=", args.user_NOESY_FILE)
if args.user_TOCSY_FILE and args.user_NOESY_FILE:
    print("STARTING PROOF-READING NOESY CHEMICAL SHIFT ASSIGNMENTS")
    aa_equivalentCarbons_dict = {}
    aa_equivalentCarbons_dict["LEU"] = [["CD1", "CD2"]]
    aa_equivalentCarbons_dict["VAL"] = [["CG1", "CG2"]]
    
    # LOAD USER  ASSIGNED TOCSY CONTENTS
    user_TOCSY_residue_assignments_dict, \
    user_TOCSY_residue_NHresonances_dict, \
    user_i_to_iminus1_dict, \
    user_patched_residues_list, \
    user_TOCSY_residue_peak_intensity_mdict = \
        read_NHassigned_spectrum_file(args.user_TOCSY_FILE, "TOCSY", ali)
    # print "DEBUG: i_to_iminus1_dict=", i_to_iminus1_dict
    print("DEBUG: user_TOCSY_residue_assignments_dict=", user_TOCSY_residue_assignments_dict)
    print("DEBUG: user_TOCSY_residue_NHresonances_dict=", user_TOCSY_residue_NHresonances_dict)
    for k,v in list(user_TOCSY_residue_assignments_dict.items()):
        print(k, "-->", v)
        
    # LOAD USER ASSIGNED NOESY FILE CONTENTS
    user_NOESY_residue_assignments_dict, \
    user_NOESY_residue_NHresonances_dict, \
    user_original_NOESY_peaks_dict, \
    user_NOESY_residue_peak_instensity_dict = \
        read_NHassigned_spectrum_file(args.user_NOESY_FILE, 'HCNH', ali)
    print("DEBUG: user_NOESY_residue_assignments_dict=", user_NOESY_residue_assignments_dict)
    print("DEBUG: user_NOESY_residue_NHresonances_dict=", user_NOESY_residue_NHresonances_dict)
    print("DEBUG: user_original_NOESY_peaks_dict=", user_original_NOESY_peaks_dict)
    # rename TYR and PHE CD-QD>CD1-HD1 and CE-QE->CE1-HE1
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
    for k,v in list(user_NOESY_residue_assignments_dict.items()):
        print(k, "-->", v)
    
    # convert the keys from residues to resids to be compatible with clean_resid_assignments_dict
    user_NOESY_resid_assignments_dict = {}  # same as user_NOESY_resid_assignments_dict, but with keys resids instead of residues
    for k,v in list(user_NOESY_residue_assignments_dict.items()):
        if not is_valid_residue(k, ali.protein_alignment_list, args, spectrum_combo='TOCSY-HCNH'):
            print("WARNING: omitting non-valid user assigned residue ", k, "!")
            continue
        user_NOESY_resid_assignments_dict[int(k[1:])] = v
    user_TOCSY_resid_assignments_dict = {}  # same as user_TOCSY_resid_assignments_dict, but with keys resids instead of residues
    for k,v in list(user_TOCSY_residue_assignments_dict.items()):
        if not is_valid_residue(k, ali.protein_alignment_list, args, spectrum_combo='TOCSY-HCNH'):
            print("WARNING: omitting non-valid user assigned residue ", k, "!")
            continue
        user_TOCSY_resid_assignments_dict[int(k[1:])] = v
    
    print("DEBUG: user_NOESY_resid_assignments_dict=", user_NOESY_resid_assignments_dict)
    print("DEBUG: ITERATION 1 point 2 clean_resid_assignments_dict=", clean_resid_assignments_dict)
    
    # CREATE ONLY ASSINGED VAL AND LEU EQUIVALENT CARBONS DICTIONARIES
    user_VAL_LEU_NOESY_resid_assignments_dict, user_VAL_LEU_TOCSY_resid_assignments_dict = create_equivalentCarbons_assignments_dict(user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict, ali.protein_alignment_list,
                                            ali.absolute_matches_alignment_list, aa_equivalentCarbons_dict, args)
    print("DEBUG: user_VAL_LEU_NOESY_resid_assignments_dict=", user_VAL_LEU_NOESY_resid_assignments_dict)
    print("DEBUG: user_VAL_LEU_TOCSY_resid_assignments_dict=", user_VAL_LEU_TOCSY_resid_assignments_dict)
    
    ##
    ##  COMPARE YOUR CARBON ASSIGNMENTS WITH THE USER ASSIGNMENTS
    ##
    ## Remember that the comments can be the following:
    # NOESY
    # NOESY (C-term hanging residue)
    # NOESY (flanked residue)
    # NOESY + TOCSY-HCNH
    # TOCSY (average)
    # TOCSY-HCNH
    # TOCSY-HCNH + NOESY
    # TOCSY-HCNH + TOCSY (unmatched)
    # TOCSY (unmatched)
    # TOCSY (unmatched) + TOCSY-HCNH
    # Recall that if we have 2 peaks in the same C-group but only one comes from NOESY ('TOCSY (unmatched) + TOCSY-HCNH', or
    # 'TOCSY-HCNH + TOCSY (unmatched)'), then we save only the NOESY Carbon resonance. 'TOCSY (unmatched)' -> only one resonance available
    # only in TOCSY, 'TOCSY (average)' -> 2 resonances available only from TOCSY and they were averaged
    
    proofread_clean_resid_assignments_dict = copy.deepcopy(clean_resid_assignments_dict)   # same as clean_resid_assignments_dict but with comments CORRECT or WRONG at the carbon lines
    for resid in list(proofread_clean_resid_assignments_dict.keys()):
        print("DEBUG ITERATION1: proof reading resid=", resid)
        resname = proofread_clean_resid_assignments_dict[resid][0][6]
        if resname in list(aa_equivalentCarbons_dict.keys()): # If this is a LEU or VAL
            Carbon_pair = aa_equivalentCarbons_dict[resname][0]
            if len([a for a in proofread_clean_resid_assignments_dict[resid] if a[3] in Carbon_pair]) > 0: # AND if the methyl carbons have been assigned
                Proton_pair = ["H"+C[1:] for C in Carbon_pair]   # protons of the equivalent Carbons
                # get only the methyl Carbon assignments (LEU --> CD1,CD2 or VAL --> CG1,CG2) made by the program, to compare them with the user made assignments
                program_assignments_list = [assignment for assignment in proofread_clean_resid_assignments_dict[resid] if assignment[3] in Carbon_pair+Proton_pair]   # keep only the resonances assigned by the program to the equivalent methyl carbons and protons
                # proof_read_equivalentCarbons(program_assignments_list, resid, user_VAL_LEU_NOESY_resid_assignments_dict, user_VAL_LEU_TOCSY_resid_assignments_dict, Carbon_pair)   # OBSOLETE!!! it doesn't matter if 1 or 2 of the peaks have been proof-read in iteration 1 & 2
                proof_read_equivalentMethyls(program_assignments_list, resid, user_VAL_LEU_NOESY_resid_assignments_dict, user_VAL_LEU_TOCSY_resid_assignments_dict, Carbon_pair, Proton_pair)   # it doesn't matter if 1 or 2 of the peaks have been proof-read in iteration 1 & 2
        # THEN PROCEED NORMALLY
        for assignment in proofread_clean_resid_assignments_dict[resid]:
            print("DEBUG ITERATION1: proof reading assignment=", assignment)
            source_spectrum = re.sub(r' iter[1-9]$', '', assignment[7])    # remove the "iter[1-9]" suffix to be able to distinguish the type of assignment (NOESY, NOESY+TOCSY, etc.)
            comment = None
            aa_type = assignment[6]
            nucleus = assignment[3] # the nucleus name
            if aa_type in list(aa_equivalentCarbons_dict.keys()) and nucleus in aa_equivalentCarbons_dict[aa_type][0]:   # if this peak has been assigned by the program to LEU or VAL CD1, CD2 or CG1,CG2 respectively, skip it
                continue
            
            # IF THIS CARBON IS A METHYL OR HAS ONLY ONE PROTON, CONSIDE ALSO THE PROTON RESONANCE FOR PROOFREADING
            if aa_type in list(aatype_carbon_nongeminalHname_mdict.keys()) and nucleus in list(aatype_carbon_nongeminalHname_mdict[aa_type].keys()):
                # get the proton assignment made by the program
                C_assignment = assignment # the program-made carbon assignment
                H_assignment = None # default value for the program-made proton assignment
                print("DEBUG: resid=", resid, "proofread_clean_resid_assignments_dict[resid]=", proofread_clean_resid_assignments_dict[resid])
                for a in proofread_clean_resid_assignments_dict[resid]:
                    if ((type(aatype_carbon_nongeminalHname_mdict[aa_type][nucleus]) == list and a[3] in aatype_carbon_nongeminalHname_mdict[aa_type][nucleus]) or
                        (type(aatype_carbon_nongeminalHname_mdict[aa_type][nucleus]) == str and a[3] == aatype_carbon_nongeminalHname_mdict[aa_type][nucleus])):
                        H_assignment = a
                        break
                print("DEBUG: C_assignment=", C_assignment, "H_assignment=", H_assignment)
                comment = ""
                if source_spectrum in ['HCNH', 'NOESY (C-term hanging residue)', 'NOESY (flanked residue)', 'NOESY + TOCSY-HCNH', 'TOCSY-HCNH', 'TOCSY-HCNH + NOESY', 'TOCSY (unmatched) + TOCSY-HCNH', 'TOCSY-HCNH + TOCSY (unmatched)']:    # if this resonance comes only from NOESY
                    user_assignment = get_Carbon_resonance(resid, nucleus, 'HCNH', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict, get_assignment=True)
                    if user_assignment == False:    # if no user assignment was made for this Carbon, add no extra comment
                        print("WARNING: no user assignment made in NOESY for resid", resid, "and Carbon", nucleus,". Checking TOCSY...")
                        user_assignment = get_Carbon_resonance(resid, nucleus, 'TOCSY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict, get_assignment=True)
                        if user_assignment == False:
                            print("WARNING: no user assignment made also in TOCSY for resid", resid, "and Carbon", nucleus,". Skipping nucleus.")
                            assignment[7] += " <NOT ASSIGNED>"
                            #proofread_clean_resid_assignments_dict[resid].append(assignment)
                            continue
                    user_Creson = user_assignment[2]
                    user_Hreson = user_assignment[4]
                    print("DEBUG: user_Creson=", user_Creson, "user_Hreson=", user_Hreson)
                    if approx_equal(C_assignment[1], user_Creson, 0.4) and approx_equal(H_assignment[1], user_Hreson, 0.04): # if the user-made assignment agrees with ours, add to the comment "CORRECT"
                        print("DEBUG: CORRECT C_assignment found:", C_assignment)
                        print("DEBUG: CORRECT H_assignment found:", H_assignment)
                        comment = " <CORRECT>"
                    else:
                        print("DEBUG: WRONG C_assignment found:", C_assignment)
                        print("DEBUG: WRONG H_assignment found:", H_assignment)
                        comment = " <WRONG>"
                elif source_spectrum in ['TOCSY (average)', 'TOCSY (unmatched)']:
                    user_assignment = get_Carbon_resonance(resid, nucleus, 'TOCSY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict, get_assignment=True)
                    if user_assignment == False:    # if no user assignment was made for this Carbon, add no extra comment
                        print("WARNING: no user assignment made in TOCSY for resid", resid, "and Carbon", nucleus,". Checking NOESY...")
                        user_assignment = get_Carbon_resonance(resid, nucleus, 'HCNH', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict, get_assignment=True)
                        if user_assignment == False:
                            print("WARNING: no user assignment made also in NOESY for resid", resid, "and Carbon", nucleus,". Skipping nucleus.")
                            assignment[7] += " <NOT ASSIGNED>"
                            #proofread_clean_resid_assignments_dict[resid].append(assignment)
                            continue
                    user_Creson = user_assignment[2]
                    user_Hreson = user_assignment[4]
                    if approx_equal(C_assignment[1], user_Creson, 0.4) and approx_equal(H_assignment[1], user_Hreson, 0.04): # if the user-made assignment agrees with ours, add to the comment "CORRECT"
                        print("DEBUG: CORRECT C_assignment found:", C_assignment)
                        print("DEBUG: CORRECT H_assignment found:", H_assignment)
                        comment = " <CORRECT>"
                    else:
                        print("DEBUG: WRONG C_assignment found:", C_assignment)
                        print("DEBUG: WRONG H_assignment found:", H_assignment)
                        comment = " <WRONG>"
                        
                if comment == "":
                    comment = " <NOT ASSIGNED>"
                C_assignment[7] += comment
                H_assignment[7] += comment
                comment = None  # a trick to add no more comments to "assignment" at the end of this block of code
            
            elif nucleus[0] == 'C':   # only if the assignment line belongs to a Carbon then check if the user has made the same assignment
                comment = ""
                if source_spectrum in ['HCNH', 'NOESY (C-term hanging residue)', 'NOESY (flanked residue)', 'NOESY + TOCSY-HCNH', 'TOCSY-HCNH', 'TOCSY-HCNH + NOESY', 'TOCSY (unmatched) + TOCSY-HCNH', 'TOCSY-HCNH + TOCSY (unmatched)']:    # if this resonance comes only from NOESY
                    user_assignment = get_Carbon_resonance(resid, nucleus, 'HCNH', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict)
                    if user_assignment == False:    # if no user assignment was made for this Carbon, add no extra comment
                        print("WARNING: no user assignment made in NOESY for resid", resid, "and Carbon", nucleus,". Checking TOCSY...")
                        user_assignment = get_Carbon_resonance(resid, nucleus, 'TOCSY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict)
                        if user_assignment == False:
                            print("WARNING: no user assignment made also in TOCSY for resid", resid, "and Carbon", nucleus,". Skipping nucleus.")
                            assignment[7] += " <NOT ASSIGNED>"
                            #proofread_clean_resid_assignments_dict[resid].append(assignment)
                            continue
                    if approx_equal(assignment[1], user_assignment, 0.4): # if the user-made assignment agrees with ours, add to the comment "CORRECT"
                        print("DEBUG: CORRECT assignment found:", assignment)
                        comment = " <CORRECT>"
                    else:
                        print("DEBUG: WRONG assignment found:", assignment)
                        comment = " <WRONG>"
                elif source_spectrum in ['TOCSY (average)', 'TOCSY (unmatched)']:
                    user_assignment = get_Carbon_resonance(resid, nucleus, 'TOCSY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict)
                    if user_assignment == False:    # if no user assignment was made for this Carbon, add no extra comment
                        print("WARNING: no user assignment made in TOCSY for resid", resid, "and Carbon", nucleus,". Checking NOESY...")
                        user_assignment = get_Carbon_resonance(resid, nucleus, 'HCNH', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict)
                        if user_assignment == False:
                            print("WARNING: no user assignment made also in NOESY for resid", resid, "and Carbon", nucleus,". Skipping nucleus.")
                            assignment[7] += " <NOT ASSIGNED>"
                            #proofread_clean_resid_assignments_dict[resid].append(assignment)
                            continue
                    if approx_equal(assignment[1], user_assignment, 0.4): # if the user-made assignment agrees with ours, add to the comment "CORRECT"
                        print("DEBUG: CORRECT assignment found:", assignment)
                        comment = " <CORRECT>"
                    else:
                        print("DEBUG: WRONG assignment found:", assignment)
                        comment = " <WRONG>"
                
            if comment == None:     # this is a proton or a methyl Carbon or a Carbon with only one proton
                # assignment[7] += iteration_comment
                continue
            elif comment == "":
                comment = " <NOT ASSIGNED>"
            assignment[7] += comment
    
    
    print("DEBUG: proofread_clean_resid_assignments_dict=", proofread_clean_resid_assignments_dict)
    
    ## WRITE THE PROOF-READ XEASY FILE
    xeasy_out = open(args.ABSOLUTE_MATCHES_FILE+".NOESY.proofread.iter1.xeasy", 'w')
    atom_index = 0  # reinitialize atom index
    for resid in list(proofread_clean_resid_assignments_dict.keys()):
        for peak in proofread_clean_resid_assignments_dict[resid]:
            atom_index += 1
            peak[0] = atom_index
            line = "\t".join([str(p) for p in peak])
            xeasy_out.write(line + "\n")
    xeasy_out.close()


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ITERATION 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
print("DEBUG: ENTERING ITERATION 2 matched_NOESY_residue_assignments_dict=", matched_NOESY_residue_assignments_dict)
print("DEBUG: ENTERING ITERATION 2 clean_resid_assignments_dict=", clean_resid_assignments_dict)

## FIND MISSING CARBON ASSIGNMENTS AND TRY TO PREDICT THEM FROM NOESY
# ATTENTION: recall that in TOCSY_residue_assignments_dict both the keys and the C,H assignments correspond to residue i-1.
# But in matched_noesy_residue_assignments_dict the keys are residue i and the values are peaks from residues i (mostly the strongest), residue i-1 and maybe i+1 or
# other neighboring residues.
# So in order to match the NOESY peaks of residue i we need primarily the values of TOCSY_residue_assignments_dict keys i respectively.
xeasy_fout = open(args.ABSOLUTE_MATCHES_FILE+".NOESY.iter2.xeasy", 'w')
sparky_fout = open(args.ABSOLUTE_MATCHES_FILE+".NOESY.iter2.sparky", 'w')
atom_index = 1
# start the peak assignment from the last residue to the first
revordered_residue_keys_iter2 = [(k[0], int(k[1:])) for k in list(TOCSY_residue_assignments_dict.keys())]
revordered_residue_keys_iter2.sort(key=itemgetter(1), reverse=True)  # reverse ordered keys, namely they start from the last residue in the protein sequence
matched_i_to_iplus1_peaks_mdict = tree()    # residue -> i peak -> matched i+1 peak
residues_with_written_prevAssigned_peaks = []

for i_aa_resid in revordered_residue_keys_iter2:    # both the keys and the C,H assignments correspond to residue i-1
    i_residue = i_aa_resid[0]+str(i_aa_resid[1])
    # if not i_residue in ['L215']:
    #    continue
    # FIRST MATCH THE TOCSY ASSIGNMENTS OF RESIDUE i TO NOESY PEAKS
    #i_residue = get_i_residue_from_iminus1_residue(iminus1_residue)
    all_possible_assignments_list = []  # list with all the peaks for this residue (TOCSY-HCNH, NOESY) and the respective possible C-H assignments
    if i_residue in list(TOCSY_residue_assignments_dict.keys()) and len(TOCSY_residue_assignments_dict[i_residue]) > 0:   # if no TOCSY assignments have been made for residue i, skip NOESY assignment of residue i
        print("ITERATION 2: Matching residue i ", i_residue," TOCSY PEAKS to residue's i NOESY peaks...")
        i_resid = get_resid_from_residue(i_residue)
        checked_residues_set.add(i_residue)
        aa_type = aa1to3_dict[get_aa_type_from_residue(i_residue)]
        # all_allowed_carbons_list = [x for x in allowed_aa_atoms_dict[aa_type] if x[0] == "C" and not (aa_type in iter1_allowed_atoms_dict.keys() and x in iter1_allowed_atoms_dict[aa_type])]   # all allowed carbons for this aa type (exclude methyls)
        all_allowed_carbons_list = [x for x in allowed_aa_atoms_dict[aa_type] if x[0] == "C" and not (aa_type in list(iter1_allowed_atoms_dict.keys()) and x in iter1_allowed_atoms_dict[aa_type])]   # all allowed carbons for this aa type (exclude methyls)
        print("DEBUG: all_allowed_carbons_list=", all_allowed_carbons_list)
        print("DEBUG: i_residue=", i_residue, "TOCSY_residue_assignments_dict[i_residue]=", TOCSY_residue_assignments_dict[i_residue])
        matched_iplus1, \
        matched_NOESY_residue_assignments_dict, \
        matched_i_to_iplus1_peaks_mdict = \
            match_HCNH_i_iplus1_peaks(i_residue,
                                      ali.absolute_matches_alignment_list,
                                      matched_NOESY_residue_assignments_dict,
                                      matched_i_to_iplus1_peaks_mdict)  # add labels to the i->i+1 matched peaks. If i_residue is a C-term residue or presends a gap, skip it from this iteration
        if matched_iplus1 == False:
            continue
        else:   # otherwise proceed to assigned but not from the extreme carbons in exclude_from_i_iplus_matching_dict
            assigned_carbons_list = []  # list of carbons which have ALL their protons assigned in TOCSY
            partly_assigned_carbons_list = []   # list of methylene carbons for which only one proton frequency was found in TOCSY
            missing_carbons_list = []
            tmp_missing_carbons_list, tmp_partly_assigned_carbons_list = get_missing_carbons_from_xeasy_dict(aa_type, i_resid, clean_resid_assignments_dict, matched_NOESY_residue_assignments_dict, allowed_aa_atoms_dict)
            print("DEBUG: tmp_missing_carbons_list=", tmp_missing_carbons_list, "tmp_partly_assigned_carbons_list=", tmp_partly_assigned_carbons_list)
            for carbon in tmp_missing_carbons_list:
                if aa_type in list(iter1_allowed_atoms_dict.keys()) and carbon in iter1_allowed_atoms_dict[aa_type]:     # to exclude methyls from this iteration
                    continue
                if not aa_type in list(exclude_from_i_iplus_matching_dict.keys()) or (not carbon in exclude_from_i_iplus_matching_dict[aa_type]):
                    missing_carbons_list.append(carbon)
            for carbon in tmp_partly_assigned_carbons_list:
                if not aa_type in list(exclude_from_i_iplus_matching_dict.keys()) or (not carbon in exclude_from_i_iplus_matching_dict[aa_type]):
                    partly_assigned_carbons_list.append(carbon)
        print("DEBUG: dict key i_residue", i_residue, "missing_carbons_list=", missing_carbons_list)
        unassigned_NOESY_peaks_list = []    # ini
        if len(missing_carbons_list) != 0:   # IF THERE ARE MISSING CARBONS IN TOCSY, TRY TO FIND THEM IN NOESY (EXCLUDE THE ASSIGNED OR MATCHED NOESY PEAKS)
            print("Residue ", i_residue, " was not completely assigned in ITERATION 1. Using NOESY to assigned resonances to the rest Carbon types.")
            try:    # DANGEROUS but works for the N-term residue (A200 in aLP)
                unassigned_NOESY_peaks_list = [peak for peak in matched_NOESY_residue_assignments_dict[i_residue] if peak[0]=="i+1" and peak[1]=="?" and peak[3]=="?"]
            except KeyError:
                continue
            print("DEBUG: unassigned_NOESY_peaks_list=", unassigned_NOESY_peaks_list)
            for NOESY_peak in unassigned_NOESY_peaks_list:  # make C-H type predictions for each unassigned NOESY peak
                NOESY_Creson = NOESY_peak[2]
                NOESY_Hreson = NOESY_peak[4]
                all_possible_assignments_list.extend( get_probabilities_from_H_C_resonpair_2Dhist(i_residue,
                                                                                                  NOESY_Hreson,
                                                                                                  NOESY_Creson,
                                                                                                  missing_carbons_list,
                                                                                                  partly_assigned_carbons_list,
                                                                                                  args,
                                                                                                  NOESY_residue_peak_intensity_mdict,
                                                                                                  matched_i_to_iplus1_peaks_mdict,
                                                                                                  aa2pthres_dict,
                                                                                                  aa2ithres_dict,
                                                                                                  ali.protein_alignment_list,
                                                                                                  ali.absolute_matches_alignment_list,
                                                                                                  histload,
                                                                                                  iteration=2)
                                                    )
            [assignment.append('HCNH') for assignment in all_possible_assignments_list]    # we don't need to created a new list, the original all_possible_assignments_list will be updated
            print("DEBUG: i_residue=", i_residue, "all_possible_assignments_list=", all_possible_assignments_list)
            written_nucleiNames_list = []
            if i_resid in list(clean_resid_assignments_dict.keys()):
                written_nucleiNames_list = [line[3]  for line in clean_resid_assignments_dict[i_resid]]   # list of nuclei that have been written already, to avoid duplication
            if len(all_possible_assignments_list) > 0:  # IF WE HAVE FOUND POSSIBLE NOESY PEAK ASSIGNMENTS OF THE MISSING CARBONS
                # WRITE THE BEST NOESY PEAK COMBINATION TO COMPLETE THE ASSIGNMENT OF THE MISSING CARBONS
                extra_written_nucleiNames_list, \
                matched_NOESY_residue_assignments_dict, \
                TOCSY_assigned_residue_list = \
                    select_best_HCNH_peak_combination(i_residue,
                                                       all_possible_assignments_list,
                                                       partly_assigned_carbons_list,
                                                       atom_index,
                                                       xeasy_fout,
                                                       sparky_fout,
                                                       ali.protein_alignment_list,
                                                       ali.absolute_matches_alignment_list,
                                                       args,
                                                       NOESY_residue_NHresonances_dict,
                                                       matched_NOESY_residue_assignments_dict,
                                                       TOCSY_assigned_residue_list,
                                                       iter1_allowed_atoms_dict,
                                                       iteration=2)
                written_nucleiNames_list.extend(extra_written_nucleiNames_list)
            # WRITE ALSO THE TOCSY-HCNH MATCHED PEAKS OF THIS RESIDUE
            residues_with_written_prevAssigned_peaks.append(i_residue)
            # if write_TOCSY_NOESY_matched_peaks(i_residue, written_nucleiNames_list) == False:
            #     continue
            # WRITE ALSO THE PEAKS ASSIGNED IN ITERATION 1
            write_peaks_assigned_in_previous_iteration(i_residue, clean_resid_assignments_dict, atom_index, xeasy_fout, 
                                               NOESY_residue_NHresonances_dict, residues_with_written_NH_list)
        elif len(missing_carbons_list) == 0 and len(partly_assigned_carbons_list) == 0: # NOTHING TO DO FOR THIS RESIDUE
            residues_with_written_prevAssigned_peaks.append(i_residue)
            # WRITE ALSO THE PEAKS ASSIGNED IN ITERATION 1 & 2
            write_peaks_assigned_in_previous_iteration(i_residue, clean_resid_assignments_dict, atom_index, xeasy_fout, 
                                               NOESY_residue_NHresonances_dict, residues_with_written_NH_list)
        else:   # IF ALL THE CARBONS WERE ASSIGNED IN TOCSY, FIND & SAVE PARTLY ASSIGNED CARBONS IN NOESY, SAVE THE MATCHED NOESY RESONANCES, AND IF NOT MATCHED, SAVE THE TOCSY RESONANCES
            print("Residue ", i_residue, " was completely assigned in ITERATION 1.")
            try:
                aa_type = aa1to3_dict[get_aa_type_from_residue(i_residue)]
            except (IndexError, KeyError):    # skip residue i
                print("EXITING: COULD NOT FIND AA_TYPE OF RESIDUE", i_residue)
                sys.exit(1)
            
            assigned_peaks_list = get_peaks_from_xeasy(i_resid, clean_resid_assignments_dict)     # get all assigned peaks of this residue from xeasy format
            for assigned_peak in assigned_peaks_list:
                source_spectrum = assigned_peak[5]
                Cname = assigned_peak[1]
                Creson = assigned_peak[2]    # use the TOCSY resonance by default
                Hname = assigned_peak[3]
                Hreson = assigned_peak[4]
                all_possible_assignments_list.append([aa_type, Cname, Hname, 10000000.0, Hreson, Creson, source_spectrum])  # set the probability to 1 cause it has been already assigned in iteration 1
                print("DEBUG: i_residue=", i_residue, "aa_type=", aa_type, "all_possible_assignments_list=", all_possible_assignments_list)
            written_nucleiNames_list = []
            if i_resid in list(clean_resid_assignments_dict.keys()):
                written_nucleiNames_list = [line[3]  for line in clean_resid_assignments_dict[i_resid]]   # list of nuclei that have been written already, to avoid duplication
            if len(all_possible_assignments_list) > 0:  # IF WE HAVE FOUND POSSIBLE NOESY PEAK ASSIGNMENTS OF THE MISSING CARBONS
                # WRITE THE BEST NOESY PEAK COMBINATION TO COMPLETE THE ASSIGNMENT OF THE MISSING CARBONS
                written_nucleiNames_list, \
                matched_NOESY_residue_assignments_dict, \
                TOCSY_assigned_residue_list = \
                    select_best_HCNH_peak_combination(i_residue,
                                                       all_possible_assignments_list,
                                                       partly_assigned_carbons_list,
                                                       atom_index,
                                                       xeasy_fout,
                                                       sparky_fout,
                                                       ali.protein_alignment_list,
                                                       ali.absolute_matches_alignment_list,
                                                       args,
                                                       NOESY_residue_NHresonances_dict,
                                                       matched_NOESY_residue_assignments_dict,
                                                       TOCSY_assigned_residue_list,
                                                       iter1_allowed_atoms_dict,
                                                       iteration=2)
        matched_NOESY_residue_assignments_dict = \
            remove_iplus1_labels(i_residue,
                                matched_NOESY_residue_assignments_dict,
                                ali.protein_alignment_list,
                                ali.absolute_matches_alignment_list,
                                args)     # remove all the 'i+1', 'i-1' labels from this residue and residue i+1 (if applicable)
        
print("DEBUG iteration2 point2: matched_NOESY_residue_assignments_dict=", matched_NOESY_residue_assignments_dict)

# WRITE ALSO THE PEAKS ASSIGNED IN PREVIOUS ITERATIONS OF ALL RESIDUES THAT WERE EXCLUDED FROM THIS ITERATION
for i_aa_resid in set(revordered_residue_keys_iter1+revordered_residue_keys_iter2):
    i_residue = i_aa_resid[0]+str(i_aa_resid[1])
    aa_type = aa1to3_dict[i_aa_resid[0]]
    if not i_residue in residues_with_written_prevAssigned_peaks:
        write_peaks_assigned_in_previous_iteration(i_residue, clean_resid_assignments_dict, atom_index, xeasy_fout, 
                                               NOESY_residue_NHresonances_dict, residues_with_written_NH_list)

# ##
# ## NOW WRITE TO XEASY ONLY FILE THE TOCSY PEAKS OF RESIDUES THAT DON'T HAVE TOCSY AT ALL (PATCHED N-TERM RESIDUES)
# ## I.e. if a residue i [R304] doesn't have TOCSY (we cannot get the resonances of residue i [G303]) but both residues i & i+1 have NOESY peaks,
# ## then use only the NOESY peaks to do the assignment for residue i. In the code, the above means that:
# ## 'G303' in TOCSY_residue_assignments_dict.keys() = False
# ## 'R304' in TOCSY_residue_assignments_dict.keys() = True
# ## 'G303' in matched_NOESY_residue_assignments_dict.keys() = True
# ## 'R304' in matched_NOESY_residue_assignments_dict.keys() = True
# ##
# ## 'N334' in TOCSY_residue_assignments_dict.keys() = False
# ## 'A335' in TOCSY_residue_assignments_dict.keys() = True
# ## 'N334' in matched_NOESY_residue_assignments_dict.keys() = True
# ## 'A335' in matched_NOESY_residue_assignments_dict.keys() = True, but 'A335' shouldn't have NOESY! It has because we named in it the root file
# as 'A335'. To Avoid this error add an extra condition was added that requires absence of gap in the alignment at the position of 'A335'

# print "DEBUG: before entering write_flanked_residues(), ali.absolute_matches_alignment_list=", ali.absolute_matches_alignment_list
revordered_flanked_residues = write_flanked_residues(ali.absolute_matches_alignment_list,
                                                     matched_NOESY_residue_assignments_dict,
                                                     TOCSY_residue_assignments_dict,
                                                     atom_index,
                                                     xeasy_fout,
                                                     sparky_fout,
                                                     patched_residues_list,
                                                     clean_resid_assignments_dict,
                                                     ali.protein_alignment_list,
                                                     args,
                                                     NOESY_residue_peak_intensity_mdict,
                                                     aa2pthres_dict,
                                                     aa2ithres_dict,
                                                     checked_residues_set,
                                                     matched_i_to_iplus1_peaks_mdict,
                                                     NOESY_residue_NHresonances_dict,
                                                     residues_with_written_NH_list,
                                                     residues_with_written_prevAssigned_peaks,
                                                     revordered_residue_keys_iter2,
                                                     iter1_allowed_atoms_dict,
                                                     allowed_aa_atoms_dict,
                                                     iteration=2)


# Now find if there are any residue in the alignment that haven't been checked:
for residue in ali.absolute_matches_alignment_list:
    if not residue in checked_residues_set and residue != '-':
        print("WARNING: Residue ", residue, " was not checked!!!")
        
# ##
# ##  COPY FROM TOCSY ASSIGNMENTS (TOCSY_residue_assignments_dict) OF ALL THE RESONANCES THAT COULD NOT BE MATCHED WITH THE NOESY ASSIGNMENTS (TOCSY_residue_assignments_dict)
# ##
# write_unmatched_TOCSY_peaks(iteration=2)


# ##
# ##  WRITE THE HN & N RESONANCES TO XEASY ONLY FOR THOSE RESIDUES THAT NO C,H ASSIGNMENTS COULD BE MADE BUT THEIR HN & H RESONANCES ARE KNOWN FROM TOCSY OR NOESY
# ##
write_NH_of_residues_without_CH(ali.absolute_matches_alignment_list, NOESY_residue_NHresonances_dict, TOCSY_residue_NHresonances_dict,
                                    residues_with_written_NH_list, atom_index, xeasy_fout, sparky_fout, ali.protein_alignment_list, args,
                                    Cterm_residue_set, iteration=2)

xeasy_fout.close()
sparky_fout.close()

##
##  REWRITE THE XEASY FILE BY REMOVING DUPLICATES OR AVERAGING, AND CORRECTING THE ORDER
##

## CLEAN THE CONTENTS OF XEASY FILE
# clean_resid_assignments_dict contains the final, non-redundant assignments written to xeasy file. E.g.
# 1: [[2336, 56.022, 0.2, 'CA', 1, '#', 'MET', 'TOCSY (unmatched) (iter2)'], [2338, 34.786, 0.2, 'CB', 1, '#', 'MET', 'TOCSY (unmatched) (iter2)'],
# [2340, 30.946, 0.2, 'CG', 1, '#', 'MET', 'TOCSY (unmatched) (iter2)'], [2337, 4.194, 0.02, 'HA', 1, '#', 'MET', 'TOCSY (unmatched) (iter2)'],
# [2339, 2.14, 0.02, 'QB', 1, '#', 'MET', 'TOCSY (unmatched) (iter2)'], [2341, 2.487, 0.02, 'QG', 1, '#', 'MET', 'TOCSY (unmatched) (iter2)']]
clean_resid_assignments_dict = clean_xeasy_file(args.ABSOLUTE_MATCHES_FILE+".NOESY.iter2.xeasy")  
# sort the peak assignments of each resid
for resid, peak_list in list(clean_resid_assignments_dict.items()):
    peak_list.sort(key=itemgetter(3))
    clean_resid_assignments_dict[resid] = peak_list

print("DEBUG: ITERATION 2 point 1 clean_resid_assignments_dict=", clean_resid_assignments_dict)
## WRITE THE CLEAN XEASY FILE
xeasy_out = open(args.ABSOLUTE_MATCHES_FILE+".NOESY.cleaned.iter2.xeasy", 'w')
atom_index = 0  # reinitialize atom index
for resid in list(clean_resid_assignments_dict.keys()):
    for peak in clean_resid_assignments_dict[resid]:
        atom_index += 1
        peak[0] = atom_index
        line = "\t".join([str(p) for p in peak])
        xeasy_out.write(line + "\n")
xeasy_out.close()


##
##  IF PROVIDED, LOAD THE USER-ANNOTATED TOCSY AND NOESY FILES FOR PROOF-READING
##
print("DEBUG: args.user_TOCSY_FILE=", args.user_TOCSY_FILE, "args.user_NOESY_FILE=", args.user_NOESY_FILE)
if args.user_TOCSY_FILE and args.user_NOESY_FILE:
    print("STARTING PROOF-READING NOESY CHEMICAL SHIFT ASSIGNMENTS")
    aa_equivalentCarbons_dict = {}
    aa_equivalentCarbons_dict["LEU"] = [["CD1", "CD2"]]
    aa_equivalentCarbons_dict["VAL"] = [["CG1", "CG2"]]
    
    # LOAD USER  ASSIGNED TOCSY CONTENTS
    user_TOCSY_residue_assignments_dict, \
    user_TOCSY_residue_NHresonances_dict, \
    user_i_to_iminus1_dict, \
    user_patched_residues_list, \
    user_TOCSY_residue_peak_intensity_mdict = \
        read_NHassigned_spectrum_file(args.user_TOCSY_FILE, "TOCSY", ali)
    # print "DEBUG: i_to_iminus1_dict=", i_to_iminus1_dict
    print("DEBUG: user_TOCSY_residue_assignments_dict=", user_TOCSY_residue_assignments_dict)
    print("DEBUG: user_TOCSY_residue_NHresonances_dict=", user_TOCSY_residue_NHresonances_dict)
    for k,v in list(user_TOCSY_residue_assignments_dict.items()):
        print(k, "-->", v)
        
    # LOAD USER ASSIGNED NOESY FILE CONTENTS
    user_NOESY_residue_assignments_dict, \
    user_NOESY_residue_NHresonances_dict, \
    user_original_NOESY_peaks_dict, \
    user_NOESY_residue_peak_instensity_dict = \
        read_NHassigned_spectrum_file(args.user_NOESY_FILE, 'HCNH', ali)
    print("DEBUG: user_NOESY_residue_assignments_dict=", user_NOESY_residue_assignments_dict)
    print("DEBUG: user_NOESY_residue_NHresonances_dict=", user_NOESY_residue_NHresonances_dict)
    print("DEBUG: user_original_NOESY_peaks_dict=", user_original_NOESY_peaks_dict)
    # rename TYR and PHE CD-QD>CD1-HD1 and CE-QE->CE1-HE1
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
    for k,v in list(user_NOESY_residue_assignments_dict.items()):
        print(k, "-->", v)
    
    # convert the keys from residues to resids to be compatible with clean_resid_assignments_dict
    user_NOESY_resid_assignments_dict = {}  # same as user_NOESY_resid_assignments_dict, but with keys resids instead of residues
    for k,v in list(user_NOESY_residue_assignments_dict.items()):
        if not is_valid_residue(k, ali.protein_alignment_list, args, spectrum_combo='TOCSY-HCNH'):
            print("WARNING: omitting non-valid user assigned residue ", k, "!")
            continue
        user_NOESY_resid_assignments_dict[int(k[1:])] = v
    user_TOCSY_resid_assignments_dict = {}  # same as user_TOCSY_resid_assignments_dict, but with keys resids instead of residues
    for k,v in list(user_TOCSY_residue_assignments_dict.items()):
        if not is_valid_residue(k, ali.protein_alignment_list, args, spectrum_combo='TOCSY-HCNH'):
            print("WARNING: omitting non-valid user assigned residue ", k, "!")
            continue
        user_TOCSY_resid_assignments_dict[int(k[1:])] = v
    
    print("DEBUG: user_NOESY_resid_assignments_dict=", user_NOESY_resid_assignments_dict)
    print("DEBUG: ITERATION 2 point 2 clean_resid_assignments_dict=", clean_resid_assignments_dict)
    
    # CREATE ONLY ASSINGED VAL AND LEU EQUIVALENT CARBONS DICTIONARIES
    user_VAL_LEU_NOESY_resid_assignments_dict, user_VAL_LEU_TOCSY_resid_assignments_dict = create_equivalentCarbons_assignments_dict(user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict, ali.protein_alignment_list,
                                                        ali.absolute_matches_alignment_list, aa_equivalentCarbons_dict, args)
    print("DEBUG: user_VAL_LEU_NOESY_resid_assignments_dict=", user_VAL_LEU_NOESY_resid_assignments_dict)
    print("DEBUG: user_VAL_LEU_TOCSY_resid_assignments_dict=", user_VAL_LEU_TOCSY_resid_assignments_dict)
    
    ##
    ##  COMPARE YOUR CARBON ASSIGNMENTS WITH THE USER ASSIGNMENTS
    ##
    ## Remember that the comments can be the following:
    # NOESY
    # NOESY (C-term hanging residue)
    # NOESY (flanked residue)
    # NOESY + TOCSY-HCNH
    # TOCSY (average)
    # TOCSY-HCNH
    # TOCSY-HCNH + NOESY
    # TOCSY-HCNH + TOCSY (unmatched)
    # TOCSY (unmatched)
    # TOCSY (unmatched) + TOCSY-HCNH
    # Recall that if we have 2 peaks in the same C-group but only one comes from NOESY ('TOCSY (unmatched) + TOCSY-HCNH', or
    # 'TOCSY-HCNH + TOCSY (unmatched)'), then we save only the NOESY Carbon resonance. 'TOCSY (unmatched)' -> only one resonance available
    # only in TOCSY, 'TOCSY (average)' -> 2 resonances available only from TOCSY and they were averaged
    
    proofread_clean_resid_assignments_dict = copy.deepcopy(clean_resid_assignments_dict)   # same as clean_resid_assignments_dict but with comments CORRECT or WRONG at the carbon lines
    for resid in list(proofread_clean_resid_assignments_dict.keys()):
        print("DEBUG ITERATION2: proof reading resid=", resid)
        resname = proofread_clean_resid_assignments_dict[resid][0][6]
        if resname in list(aa_equivalentCarbons_dict.keys()): # If this is a LEU or VAL
            Carbon_pair = aa_equivalentCarbons_dict[resname][0]
            if len([a for a in proofread_clean_resid_assignments_dict[resid] if a[3] in Carbon_pair]) > 0: # AND if the methyl carbons have been assigned
                Proton_pair = ["H"+C[1:] for C in Carbon_pair]   # protons of the equivalent Carbons
                # get only the methyl Carbon assignments (LEU --> CD1,CD2 or VAL --> CG1,CG2) made by the program, to compare them with the user made assignments
                program_assignments_list = [assignment for assignment in proofread_clean_resid_assignments_dict[resid] if assignment[3] in Carbon_pair+Proton_pair]   # keep only the resonances assigned by the program to the equivalent methyl carbons and protons
                # proof_read_equivalentCarbons(program_assignments_list, resid, user_VAL_LEU_NOESY_resid_assignments_dict, user_VAL_LEU_TOCSY_resid_assignments_dict, Carbon_pair)   # OBSOLETE!!! it doesn't matter if 1 or 2 of the peaks have been proof-read in iteration 1 & 2
                proof_read_equivalentMethyls(program_assignments_list, resid, user_VAL_LEU_NOESY_resid_assignments_dict, user_VAL_LEU_TOCSY_resid_assignments_dict, Carbon_pair, Proton_pair)   # it doesn't matter if 1 or 2 of the peaks have been proof-read in iteration 1 & 2
        # THEN PROCEED NORMALLY
        for assignment in proofread_clean_resid_assignments_dict[resid]:
            print("DEBUG ITERATION2: proof reading assignment=", assignment)
            source_spectrum = re.sub(r' iter[1-9]$', '', assignment[7])    # remove the "iter[1-9]" suffix to be able to distinguish the type of assignment (NOESY, NOESY+TOCSY, etc.)
            comment = None
            aa_type = assignment[6]
            nucleus = assignment[3] # the nucleus name
            if aa_type in list(aa_equivalentCarbons_dict.keys()) and nucleus in aa_equivalentCarbons_dict[aa_type][0]:   # if this peak has been assigned by the program to LEU or VAL CD1, CD2 or CG1,CG2 respectively, skip it
                continue
            
            # IF THIS CARBON IS A METHYL OR HAS ONLY ONE PROTON, CONSIDE ALSO THE PROTON RESONANCE FOR PROOFREADING
            if aa_type in list(aatype_carbon_nongeminalHname_mdict.keys()) and nucleus in list(aatype_carbon_nongeminalHname_mdict[aa_type].keys()):
                # get the proton assignment made by the program
                C_assignment = assignment # the program-made carbon assignment
                H_assignment = None # default value for the program-made proton assignment
                print("DEBUG: resid=", resid, "proofread_clean_resid_assignments_dict[resid]=", proofread_clean_resid_assignments_dict[resid])
                for a in proofread_clean_resid_assignments_dict[resid]:
                    if ((type(aatype_carbon_nongeminalHname_mdict[aa_type][nucleus]) == list and a[3] in aatype_carbon_nongeminalHname_mdict[aa_type][nucleus]) or
                        (type(aatype_carbon_nongeminalHname_mdict[aa_type][nucleus]) == str and a[3] == aatype_carbon_nongeminalHname_mdict[aa_type][nucleus])):
                        H_assignment = a
                        break
                print("DEBUG: C_assignment=", C_assignment, "H_assignment=", H_assignment)
                comment = ""
                if source_spectrum in ['HCNH', 'NOESY (C-term hanging residue)', 'NOESY (flanked residue)', 'NOESY + TOCSY-HCNH', 'TOCSY-HCNH', 'TOCSY-HCNH + NOESY', 'TOCSY (unmatched) + TOCSY-HCNH', 'TOCSY-HCNH + TOCSY (unmatched)']:    # if this resonance comes only from NOESY
                    user_assignment = get_Carbon_resonance(resid, nucleus, 'HCNH', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict, get_assignment=True)
                    if user_assignment == False:    # if no user assignment was made for this Carbon, add no extra comment
                        print("WARNING: no user assignment made in NOESY for resid", resid, "and Carbon", nucleus,". Checking TOCSY...")
                        user_assignment = get_Carbon_resonance(resid, nucleus, 'TOCSY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict, get_assignment=True)
                        if user_assignment == False:
                            print("WARNING: no user assignment made also in TOCSY for resid", resid, "and Carbon", nucleus,". Skipping nucleus.")
                            assignment[7] += " <NOT ASSIGNED>"
                            #proofread_clean_resid_assignments_dict[resid].append(assignment)
                            continue
                    user_Creson = user_assignment[2]
                    user_Hreson = user_assignment[4]
                    print("DEBUG: user_Creson=", user_Creson, "user_Hreson=", user_Hreson)
                    if approx_equal(C_assignment[1], user_Creson, 0.4) and approx_equal(H_assignment[1], user_Hreson, 0.04): # if the user-made assignment agrees with ours, add to the comment "CORRECT"
                        print("DEBUG: CORRECT C_assignment found:", C_assignment)
                        print("DEBUG: CORRECT H_assignment found:", H_assignment)
                        comment = " <CORRECT>"
                    else:
                        print("DEBUG: WRONG C_assignment found:", C_assignment)
                        print("DEBUG: WRONG H_assignment found:", H_assignment)
                        comment = " <WRONG>"
                elif source_spectrum in ['TOCSY (average)', 'TOCSY (unmatched)']:
                    user_assignment = get_Carbon_resonance(resid, nucleus, 'TOCSY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict, get_assignment=True)
                    if user_assignment == False:    # if no user assignment was made for this Carbon, add no extra comment
                        print("WARNING: no user assignment made in TOCSY for resid", resid, "and Carbon", nucleus,". Checking NOESY...")
                        user_assignment = get_Carbon_resonance(resid, nucleus, 'HCNH', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict, get_assignment=True)
                        if user_assignment == False:
                            print("WARNING: no user assignment made also in NOESY for resid", resid, "and Carbon", nucleus,". Skipping nucleus.")
                            assignment[7] += " <NOT ASSIGNED>"
                            #proofread_clean_resid_assignments_dict[resid].append(assignment)
                            continue
                    user_Creson = user_assignment[2]
                    user_Hreson = user_assignment[4]
                    if approx_equal(C_assignment[1], user_Creson, 0.4) and approx_equal(H_assignment[1], user_Hreson, 0.04): # if the user-made assignment agrees with ours, add to the comment "CORRECT"
                        print("DEBUG: CORRECT C_assignment found:", C_assignment)
                        print("DEBUG: CORRECT H_assignment found:", H_assignment)
                        comment = " <CORRECT>"
                    else:
                        print("DEBUG: WRONG C_assignment found:", C_assignment)
                        print("DEBUG: WRONG H_assignment found:", H_assignment)
                        comment = " <WRONG>"
                        
                if comment == "":
                    comment = " <NOT ASSIGNED>"
                C_assignment[7] += comment
                H_assignment[7] += comment
                comment = None  # a trick to add no more comments to "assignment" at the end of this block of code
            
            elif nucleus[0] == 'C':   # only if the assignment line belongs to a Carbon then check if the user has made the same assignment
                comment = ""
                if source_spectrum in ['HCNH', 'NOESY (C-term hanging residue)', 'NOESY (flanked residue)', 'NOESY + TOCSY-HCNH', 'TOCSY-HCNH', 'TOCSY-HCNH + NOESY', 'TOCSY (unmatched) + TOCSY-HCNH', 'TOCSY-HCNH + TOCSY (unmatched)']:    # if this resonance comes only from NOESY
                    user_assignment = get_Carbon_resonance(resid, nucleus, 'HCNH', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict)
                    if user_assignment == False:    # if no user assignment was made for this Carbon, add no extra comment
                        print("WARNING: no user assignment made in NOESY for resid", resid, "and Carbon", nucleus,". Checking TOCSY...")
                        user_assignment = get_Carbon_resonance(resid, nucleus, 'TOCSY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict)
                        if user_assignment == False:
                            print("WARNING: no user assignment made also in TOCSY for resid", resid, "and Carbon", nucleus,". Skipping nucleus.")
                            assignment[7] += " <NOT ASSIGNED>"
                            #proofread_clean_resid_assignments_dict[resid].append(assignment)
                            continue
                    if approx_equal(assignment[1], user_assignment, 0.4): # if the user-made assignment agrees with ours, add to the comment "CORRECT"
                        print("DEBUG: CORRECT assignment found:", assignment)
                        comment = " <CORRECT>"
                    else:
                        print("DEBUG: WRONG assignment found:", assignment)
                        comment = " <WRONG>"
                elif source_spectrum in ['TOCSY (average)', 'TOCSY (unmatched)']:
                    user_assignment = get_Carbon_resonance(resid, nucleus, 'TOCSY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict)
                    if user_assignment == False:    # if no user assignment was made for this Carbon, add no extra comment
                        print("WARNING: no user assignment made in TOCSY for resid", resid, "and Carbon", nucleus,". Checking NOESY...")
                        user_assignment = get_Carbon_resonance(resid, nucleus, 'HCNH', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict)
                        if user_assignment == False:
                            print("WARNING: no user assignment made also in NOESY for resid", resid, "and Carbon", nucleus,". Skipping nucleus.")
                            assignment[7] += " <NOT ASSIGNED>"
                            #proofread_clean_resid_assignments_dict[resid].append(assignment)
                            continue
                    if approx_equal(assignment[1], user_assignment, 0.4): # if the user-made assignment agrees with ours, add to the comment "CORRECT"
                        print("DEBUG: CORRECT assignment found:", assignment)
                        comment = " <CORRECT>"
                    else:
                        print("DEBUG: WRONG assignment found:", assignment)
                        comment = " <WRONG>"
                
            if comment == None:     # this is a proton or a methyl Carbon or a Carbon with only one proton
                # assignment[7] += iteration_comment
                continue
            elif comment == "":
                comment = " <NOT ASSIGNED>"
            assignment[7] += comment
    
    
    print("DEBUG: proofread_clean_resid_assignments_dict=", proofread_clean_resid_assignments_dict)
    
    ## WRITE THE PROOF-READ XEASY FILE
    xeasy_out = open(args.ABSOLUTE_MATCHES_FILE+".NOESY.proofread.iter2.xeasy", 'w')
    atom_index = 0  # reinitialize atom index
    for resid in list(proofread_clean_resid_assignments_dict.keys()):
        for peak in proofread_clean_resid_assignments_dict[resid]:
            atom_index += 1
            peak[0] = atom_index
            line = "\t".join([str(p) for p in peak])
            xeasy_out.write(line + "\n")
    xeasy_out.close()


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ITERATION 3 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
print("DEBUG: ENTERING ITERATION 3 matched_NOESY_residue_assignments_dict=", matched_NOESY_residue_assignments_dict)
print("DEBUG: ENTERING ITERATION 3 clean_resid_assignments_dict=", clean_resid_assignments_dict)

xeasy_fout = open(args.ABSOLUTE_MATCHES_FILE+".NOESY.iter3.xeasy", 'w')
sparky_fout = open(args.ABSOLUTE_MATCHES_FILE+".NOESY.iter3.sparky", 'w')
atom_index = 1
revordered_residue_keys_iter3 = [(k[0], int(k[1:])) for k in list(matched_NOESY_residue_assignments_dict.keys()) if k in ali.absolute_matches_alignment_list]
revordered_residue_keys_iter3.sort(key=itemgetter(1), reverse=True)  # reverse ordered keys, namely they start from the last residue in the protein sequence
residues_with_written_prevAssigned_peaks = []

selected_revordered_residue_keys_iter3 = []
for aa_type  in ["ALA", "THR", "VAL", "ILE", "LEU"]:
    for i_aa_resid in revordered_residue_keys_iter3:
        i_residue = i_aa_resid[0]+str(i_aa_resid[1])
        aa = aa1to3_dict[get_aa_type_from_residue(i_residue)]
        if aa == aa_type:
            selected_revordered_residue_keys_iter3.append(i_aa_resid)
print("DEBUG: selected_revordered_residue_keys_iter3=", selected_revordered_residue_keys_iter3)
# selected_revordered_residue_keys_iter3 = revordered_residue_keys_iter3


## THIS TIME ITERATE OVER ALL RESIDUE THAT HAVE NOESY AND ARE IN THE ALIGNMENT (INCLUDING FLANKED AND C-TERMINAL)
for i_aa_resid in selected_revordered_residue_keys_iter3:    # both the keys and the C,H assignments correspond to residue i-1
    if aa_type in list(aa2ithres_dict.keys()):
        istart = int(10 * aa2ithres_dict[aa_type])
    else:
        istart = 5
    for I_THRESHOLD in reversed(list(range(1, istart+1, 1))):
        I_THRESHOLD *= 0.1
        i_residue = i_aa_resid[0]+str(i_aa_resid[1])
        aa_type = aa1to3_dict[i_aa_resid[0]]
        i_resid = get_resid_from_residue(i_residue)
        if not aa_type in list(iter1_allowed_atoms_dict.keys()):
            continue
        # if not i_residue in ['I206']:
        #    continue
        # FIRST MATCH THE TOCSY ASSIGNMENTS OF RESIDUE i TO NOESY PEAKS
        #i_residue = get_i_residue_from_iminus1_residue(iminus1_residue)
        all_possible_assignments_list = []  # list with all the peaks for this residue (TOCSY-HCNH, NOESY) and the respective possible C-H assignments
        # TEMP DEACTIVATE: if i_residue in TOCSY_residue_assignments_dict.keys() and len(TOCSY_residue_assignments_dict[i_residue]) > 0:   # if no TOCSY assignments have been made for residue i, skip NOESY assignment of residue i
        print("ITERATION 3: Matching residue i ", i_residue," TOCSY PEAKS to residue's i NOESY peaks... I_THRESHOLD=", I_THRESHOLD)
        checked_residues_set.add(i_residue)
        tmp_missing_carbons_list, tmp_partly_assigned_carbons_list = get_missing_carbons_from_xeasy_dict(aa_type, i_resid, clean_resid_assignments_dict, matched_NOESY_residue_assignments_dict, allowed_aa_atoms_dict)
        # Keep only the allowed methyl carbons
        missing_carbons_list = [c for c in tmp_missing_carbons_list if c in iter1_allowed_atoms_dict[aa_type]]
        partly_assigned_carbons_list = [c for c in tmp_partly_assigned_carbons_list if c in iter1_allowed_atoms_dict[aa_type]]
        print("DEBUG: missing_carbons_list=", missing_carbons_list, "partly_assigned_carbons_list=", partly_assigned_carbons_list)
        unassigned_NOESY_peaks_list = []
        if len(missing_carbons_list) != 0:   # IF THERE ARE MISSING CARBONS IN ITERATION 1, TRY TO FIND THEM IN NOESY (EXCLUDE THE ASSIGNED OR MATCHED NOESY PEAKS)
            print("Residue ", i_residue, " was not completely assigned in iteration 1. Using NOESY to assigned resonances to the rest Carbon types.")
            try:    # DANGEROUS but works for the N-term residue (A200 in aLP)
                unassigned_NOESY_peaks_list = [peak for peak in matched_NOESY_residue_assignments_dict[i_residue] if peak[0]=="?" and peak[1]=="?" and peak[3]=="?"]
            except KeyError:    # if this residue is not in the alignment (not in matched_NOESY_residue_assignments_dict), then we cannot assign the partly assigned Carbons and hence we move on to the next
                print("Residue ", i_residue, " is either not in the alignment or doesn't have NOESY. Therefore it is impossible to assign partly assigned Carbon peaks.")
                continue
            for NOESY_peak in unassigned_NOESY_peaks_list:  # make C-H type predictions for each unassigned NOESY peak
                NOESY_Creson = NOESY_peak[2]
                NOESY_Hreson = NOESY_peak[4]
                all_possible_assignments_list.extend( get_probabilities_from_H_C_resonpair_2Dhist(i_residue,
                                                                                                  NOESY_Hreson,
                                                                                                  NOESY_Creson,
                                                                                                  missing_carbons_list,
                                                                                                  partly_assigned_carbons_list,
                                                                                                  args,
                                                                                                  NOESY_residue_peak_intensity_mdict,
                                                                                                  matched_i_to_iplus1_peaks_mdict,
                                                                                                  aa2pthres_dict,
                                                                                                  aa2ithres_dict,
                                                                                                  ali.protein_alignment_list,
                                                                                                  ali.absolute_matches_alignment_list,
                                                                                                  histload,
                                                                                                  iteration=3,
                                                                                                  I_THRESHOLD=I_THRESHOLD)
                                                    )
            [assignment.append('HCNH') for assignment in all_possible_assignments_list]    # we don't need to created a new list, the original all_possible_assignments_list will be updated
            print("DEBUG: i_residue=", i_residue, "all_possible_assignments_list=", all_possible_assignments_list)
            if i_resid in list(clean_resid_assignments_dict.keys()):  # if this is a C-term residue that is analyzed for the first time in iter3
                written_nucleiNames_list = [line[3]  for line in clean_resid_assignments_dict[i_resid]]   # list of nuclei that have been written already, to avoid duplication
            else:
                written_nucleiNames_list = []
            if len(all_possible_assignments_list) > 0:  # IF WE HAVE FOUND POSSIBLE NOESY PEAK ASSIGNMENTS OF THE MISSING CARBONS
                # WRITE THE BEST NOESY PEAK COMBINATION TO COMPLETE THE ASSIGNMENT OF THE MISSING CARBONS
                extra_written_nucleiNames_list, \
                matched_NOESY_residue_assignments_dict, \
                TOCSY_assigned_residue_list = \
                    select_best_HCNH_peak_combination(i_residue,
                                                       all_possible_assignments_list,
                                                       partly_assigned_carbons_list,
                                                       atom_index,
                                                       xeasy_fout,
                                                       sparky_fout,
                                                       ali.protein_alignment_list,
                                                       ali.absolute_matches_alignment_list,
                                                       args,
                                                       NOESY_residue_NHresonances_dict,
                                                       matched_NOESY_residue_assignments_dict,
                                                       TOCSY_assigned_residue_list,
                                                       iter1_allowed_atoms_dict,
                                                       iteration=3)
                written_nucleiNames_list.extend(extra_written_nucleiNames_list)
            residues_with_written_prevAssigned_peaks.append(i_residue)
            # WRITE ALSO THE PEAKS ASSIGNED IN ITERATION 1 & 2
            write_peaks_assigned_in_previous_iteration(i_residue, clean_resid_assignments_dict, atom_index, xeasy_fout, 
                                               NOESY_residue_NHresonances_dict, residues_with_written_NH_list)
        elif len(missing_carbons_list) == 0 and len(partly_assigned_carbons_list) == 0: # NOTHING TO DO FOR THIS RESIDUE
            residues_with_written_prevAssigned_peaks.append(i_residue)
            # WRITE ALSO THE PEAKS ASSIGNED IN ITERATION 1 & 2
            write_peaks_assigned_in_previous_iteration(i_residue, clean_resid_assignments_dict, atom_index, xeasy_fout, 
                                               NOESY_residue_NHresonances_dict, residues_with_written_NH_list)
        else:   # IF ALL THE CARBONS WERE ASSIGNED IN ITERATION 1, SAVE THEM, BUT ALSO FIND IN NOESY AND SAVE PARTLY ASSIGNED CARBONS
            print("Residue ", i_residue, " was completely assigned in iteration 1.")
            if not i_residue in list(matched_NOESY_residue_assignments_dict.keys()):
                print("Residue ", i_residue, " is either not in the alignment or doesn't have NOESY. Therefore it is impossible to assign partly assigned Carbon peaks.")
                continue
            try:
                aa_type = aa1to3_dict[get_aa_type_from_residue(i_residue)]
            except (IndexError, KeyError):    # skip residue i
                print("EXITING: COULD NOT FIND AA_TYPE OF RESIDUE", i_residue)
                sys.exit(1)
            assigned_peaks_list = get_peaks_from_xeasy(i_resid, clean_resid_assignments_dict)     # get all assigned peaks of this residue from xeasy format
            for assigned_peak in assigned_peaks_list:
                source_spectrum = assigned_peak[5]
                Cname = assigned_peak[1]
                Creson = assigned_peak[2]    # use the TOCSY resonance by default
                Hname = assigned_peak[3]
                Hreson = assigned_peak[4]
                all_possible_assignments_list.append([aa_type, Cname, Hname, 10000000.0, Hreson, Creson, source_spectrum])  # set the probability to 1 cause it has been already assigned in iteration 1
                print("DEBUG: i_residue=", i_residue, "aa_type=", aa_type, "all_possible_assignments_list=", all_possible_assignments_list)
            written_nucleiNames_list = []
            if i_resid in list(clean_resid_assignments_dict.keys()):
                written_nucleiNames_list = [line[3]  for line in clean_resid_assignments_dict[i_resid]]  # list of nuclei that have been written already, to avoid duplication
            if len(all_possible_assignments_list) > 0:  # IF WE HAVE FOUND POSSIBLE NOESY PEAK ASSIGNMENTS OF THE MISSING CARBONS
                # WRITE THE BEST NOESY PEAK COMBINATION TO COMPLETE THE ASSIGNMENT OF THE MISSING CARBONS
                extra_written_nucleiNames_list, \
                matched_NOESY_residue_assignments_dict, \
                TOCSY_assigned_residue_list = \
                    select_best_HCNH_peak_combination(i_residue,
                                                       all_possible_assignments_list,
                                                       partly_assigned_carbons_list,
                                                       atom_index,
                                                       xeasy_fout,
                                                       sparky_fout,
                                                       ali.protein_alignment_list,
                                                       ali.absolute_matches_alignment_list,
                                                       args,
                                                       NOESY_residue_NHresonances_dict,
                                                       matched_NOESY_residue_assignments_dict,
                                                       TOCSY_assigned_residue_list,
                                                       iter1_allowed_atoms_dict,
                                                       iteration=3)
                written_nucleiNames_list.extend(extra_written_nucleiNames_list)


# WRITE ALSO THE PEAKS ASSIGNED IN PREVIOUS ITERATIONS OF ALL RESIDUES THAT WERE EXCLUDED FROM THIS ITERATION
revordered_flanked_residue_keys = [(r[0], r[1:]) for r in revordered_flanked_residues]  # convert flanked residues to keys
for i_aa_resid in set(revordered_residue_keys_iter1 + revordered_residue_keys_iter2 + selected_revordered_residue_keys_iter3 + revordered_flanked_residue_keys):  # the same residues were assigned in both ITER1 and 2
    i_residue = i_aa_resid[0]+str(i_aa_resid[1])
    aa_type = aa1to3_dict[i_aa_resid[0]]
    if not i_residue in residues_with_written_prevAssigned_peaks:
        write_peaks_assigned_in_previous_iteration(i_residue, clean_resid_assignments_dict, atom_index, xeasy_fout, 
                                               NOESY_residue_NHresonances_dict, residues_with_written_NH_list)


# Now find if there are any residue in the alignment that haven't been checked:
for residue in ali.absolute_matches_alignment_list:
    if not residue in checked_residues_set and residue != '-':
        print("WARNING: Residue ", residue, " was not checked!!!")
        

# ##
# ##  COPY FROM TOCSY ASSIGNMENTS (TOCSY_residue_assignments_dict) OF ALL THE RESONANCES THAT COULD NOT BE MATCHED WITH THE NOESY ASSIGNMENTS (TOCSY_residue_assignments_dict)
# ##
# write_unmatched_TOCSY_peaks(iteration=3)


# ##
# ##  WRITE THE HN & N RESONANCES TO XEASY ONLY FOR THOSE RESIDUES THAT NO C,H ASSIGNMENTS COULD BE MADE BUT THEIR HN & H RESONANCES ARE KNOWN FROM TOCSY OR NOESY
# ##
write_NH_of_residues_without_CH(ali.absolute_matches_alignment_list, NOESY_residue_NHresonances_dict, TOCSY_residue_NHresonances_dict,
                                    residues_with_written_NH_list, atom_index, xeasy_fout, sparky_fout, ali.protein_alignment_list, args,
                                    Cterm_residue_set, iteration=3)

xeasy_fout.close()
sparky_fout.close()

##
##  REWRITE THE XEASY FILE BY REMOVING DUPLICATES OR AVERAGING, AND CORRECTING THE ORDER
##
print("DEBUG: ITERATION 3 point 0 clean_resid_assignments_dict=", clean_resid_assignments_dict)
## CLEAN THE CONTENTS OF XEASY FILE
clean_resid_assignments_dict = clean_xeasy_file(args.ABSOLUTE_MATCHES_FILE+".NOESY.iter3.xeasy")  # dict with all the assignments from each residue from iteration 1
# sort the peak assignments of each resid
for resid, peak_list in list(clean_resid_assignments_dict.items()):
    peak_list.sort(key=itemgetter(3))
    clean_resid_assignments_dict[resid] = peak_list

print("DEBUG: ITERATION 3 point 1 clean_resid_assignments_dict=", clean_resid_assignments_dict)

## WRITE THE CLEAN XEASY FILE
xeasy_out = open(args.ABSOLUTE_MATCHES_FILE+".NOESY.iter3.xeasy", 'w')
atom_index = 0  # reinitialize atom index
for resid in list(clean_resid_assignments_dict.keys()):
    for peak in clean_resid_assignments_dict[resid]:
        atom_index += 1
        peak[0] = atom_index
        line = "\t".join([str(p) for p in peak])
        xeasy_out.write(line + "\n")
xeasy_out.close()


##
##  REWRITE THE XEASY FILE BY REMOVING DUPLICATES OR AVERAGING, AND CORRECTING THE ORDER
##

## CLEAN THE CONTENTS OF XEASY FILE
# clean_resid_assignments_dict contains the final, non-redundant assignments written to xeasy file. E.g.
# 1: [[2336, 56.022, 0.2, 'CA', 1, '#', 'MET', 'TOCSY (unmatched) (iter1)'], [2338, 34.786, 0.2, 'CB', 1, '#', 'MET', 'TOCSY (unmatched) (iter1)'],
# [2340, 30.946, 0.2, 'CG', 1, '#', 'MET', 'TOCSY (unmatched) (iter1)'], [2337, 4.194, 0.02, 'HA', 1, '#', 'MET', 'TOCSY (unmatched) (iter1)'],
# [2339, 2.14, 0.02, 'QB', 1, '#', 'MET', 'TOCSY (unmatched) (iter1)'], [2341, 2.487, 0.02, 'QG', 1, '#', 'MET', 'TOCSY (unmatched) (iter1)']]
clean_resid_assignments_dict = clean_xeasy_file(args.ABSOLUTE_MATCHES_FILE+".NOESY.iter3.xeasy")  
# sort the peak assignments of each resid
for resid, peak_list in list(clean_resid_assignments_dict.items()):
    peak_list.sort(key=itemgetter(3))
    clean_resid_assignments_dict[resid] = peak_list

print("DEBUG: ITERATION 3 point 2 clean_resid_assignments_dict=", clean_resid_assignments_dict)

## WRITE THE CLEAN XEASY FILE
xeasy_out = open(args.ABSOLUTE_MATCHES_FILE+".NOESY.cleaned.iter3.xeasy", 'w')
atom_index = 0  # reinitialize atom index
for resid in list(clean_resid_assignments_dict.keys()):
    for peak in clean_resid_assignments_dict[resid]:
        atom_index += 1
        peak[0] = atom_index
        line = "\t".join([str(p) for p in peak])
        xeasy_out.write(line + "\n")
xeasy_out.close()


##
##  IF PROVIDED, LOAD THE USER-ANNOTATED TOCSY AND NOESY FILES FOR PROOF-READING
##
print("DEBUG: args.user_TOCSY_FILE=", args.user_TOCSY_FILE, "args.user_NOESY_FILE=", args.user_NOESY_FILE)
if args.user_TOCSY_FILE and args.user_NOESY_FILE:
    print("STARTING PROOF-READING NOESY CHEMICAL SHIFT ASSIGNMENTS")
    # THIS IS ITERATION 3, SO THERE ARE SOME EXTRA LINES THAT DON'T EXIST WHICH NEED TO BE ADDED
    for resid in clean_resid_assignments_dict:
        if not resid in list(proofread_clean_resid_assignments_dict.keys()):
            proofread_clean_resid_assignments_dict[resid] = clean_resid_assignments_dict[resid]
            continue
        for new_assignment in clean_resid_assignments_dict[resid]:
            new_reson = new_assignment[1]
            old_reson_list = [a[1] for a in proofread_clean_resid_assignments_dict[resid]]
            if not new_reson in old_reson_list:     # if this resonance was not found (newly assigned), add it to the dictionary
                proofread_clean_resid_assignments_dict[resid].append(new_assignment)
    print("DEBUG: ITERATION 3 point 3 clean_resid_assignments_dict=", clean_resid_assignments_dict)
    
    ##
    ##  COMPARE YOUR CARBON ASSIGNMENTS WITH THE USER ASSIGNMENTS
    ##
    ## Remember that the comments can be the following:
    # NOESY
    # NOESY (C-term hanging residue)
    # NOESY (flanked residue)
    # NOESY + TOCSY-HCNH
    # TOCSY (average)
    # TOCSY-HCNH
    # TOCSY-HCNH + NOESY
    # TOCSY-HCNH + TOCSY (unmatched)
    # TOCSY (unmatched)
    # TOCSY (unmatched) + TOCSY-HCNH
    # Recall that if we have 2 peaks in the same C-group but only one comes from NOESY ('TOCSY (unmatched) + TOCSY-HCNH', or
    # 'TOCSY-HCNH + TOCSY (unmatched)'), then we save only the NOESY Carbon resonance. 'TOCSY (unmatched)' -> only one resonance available
    # only in TOCSY, 'TOCSY (average)' -> 2 resonances available only from TOCSY and they were averaged
    
    
    for resid in list(proofread_clean_resid_assignments_dict.keys()):
        print("DEBUG ITERATION3: proof reading resid=", resid)
        resname = proofread_clean_resid_assignments_dict[resid][0][6]
        if resname in list(aa_equivalentCarbons_dict.keys()): # If this is a LEU or VAL
            Carbon_pair = aa_equivalentCarbons_dict[resname][0]
            if len([a for a in proofread_clean_resid_assignments_dict[resid] if a[3] in Carbon_pair]) > 0: # AND if the methyl carbons have been assigned
                Proton_pair = ["H"+C[1:] for C in Carbon_pair]   # protons of the equivalent Carbons
                # get only the methyl Carbon assignments (LEU --> CD1,CD2 or VAL --> CG1,CG2) made by the program, to compare them with the user made assignments
                print("DEBUG: proofread_clean_resid_assignments_dict[resid]=", proofread_clean_resid_assignments_dict[resid])
                program_assignments_list = [assignment for assignment in proofread_clean_resid_assignments_dict[resid] if assignment[3] in Carbon_pair+Proton_pair]   # keep only the resonances assigned by the program to the equivalent methyl carbons and protons
                # proof_read_equivalentCarbons(program_assignments_list, resid, user_VAL_LEU_NOESY_resid_assignments_dict, user_VAL_LEU_TOCSY_resid_assignments_dict, Carbon_pair)   # OBSOLETE!!! it doesn't matter if 1 or 2 of the peaks have been proof-read in iteration 1 & 2
                proof_read_equivalentMethyls(program_assignments_list, resid, user_VAL_LEU_NOESY_resid_assignments_dict, user_VAL_LEU_TOCSY_resid_assignments_dict, Carbon_pair, Proton_pair)   # it doesn't matter if 1 or 2 of the peaks have been proof-read in iteration 1 & 2
            print("DEBUG ITERATION3: after proof_read_equivalentMethyls() proofread_clean_resid_assignments_dict[",resid,"]=", proofread_clean_resid_assignments_dict[resid])
        # THEN PROCEED NORMALLY
        for assignment in proofread_clean_resid_assignments_dict[resid]:
            if '<' in assignment[7]:    # if this assignment has been proofread, then skip it
                continue
            print("DEBUG ITERATION3: proof reading new assignment=", assignment)
            spectrum_combo = re.sub(r' iter[1-9]', '', assignment[7])    # remove the "iter[1-9]" suffix to be able to distinguish the type of assignment (NEESY, NOESY+TOCSY, etc.)
            comment = None
            aa_type = assignment[6]
            nucleus = assignment[3] # the nucleus name
            if aa_type in list(aa_equivalentCarbons_dict.keys()) and nucleus in aa_equivalentCarbons_dict[aa_type][0]:    # if this peak has been assigned by the program to LEU or VAL CD1, CD2 or CG1,CG2 respectively, skip it
                continue
            
            # IF THIS CARBON IS A METHYL OR HAS ONLY ONE PROTON, CONSIDE ALSO THE PROTON RESONANCE FOR PROOFREADING
            if aa_type in list(aatype_carbon_nongeminalHname_mdict.keys()) and nucleus in list(aatype_carbon_nongeminalHname_mdict[aa_type].keys()):
                # get the proton assignment made by the program
                C_assignment = assignment # the program-made carbon assignment
                H_assignment = None # default value for the program-made proton assignment
                print("DEBUG: resid=", resid, "proofread_clean_resid_assignments_dict[resid]=", proofread_clean_resid_assignments_dict[resid])
                for a in proofread_clean_resid_assignments_dict[resid]:
                    if ((type(aatype_carbon_nongeminalHname_mdict[aa_type][nucleus]) == list and a[3] in aatype_carbon_nongeminalHname_mdict[aa_type][nucleus]) or
                        (type(aatype_carbon_nongeminalHname_mdict[aa_type][nucleus]) == str and a[3] == aatype_carbon_nongeminalHname_mdict[aa_type][nucleus])):
                        H_assignment = a
                        break
                print("DEBUG: C_assignment=", C_assignment, "H_assignment=", H_assignment)
                comment = ""
                if spectrum_combo in ['HCNH', 'NOESY (C-term hanging residue)', 'NOESY (flanked residue)', 'NOESY + TOCSY-HCNH', 'TOCSY-HCNH', 'TOCSY-HCNH + NOESY', 'TOCSY (unmatched) + TOCSY-HCNH', 'TOCSY-HCNH + TOCSY (unmatched)']:    # if this resonance comes only from NOESY
                    user_assignment = get_Carbon_resonance(resid, nucleus, 'HCNH', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict, get_assignment=True)
                    if user_assignment == False:    # if no user assignment was made for this Carbon, add no extra comment
                        print("WARNING: no user assignment made in NOESY for resid", resid, "and Carbon", nucleus,". Checking TOCSY...")
                        user_assignment = get_Carbon_resonance(resid, nucleus, 'TOCSY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict, get_assignment=True)
                        if user_assignment == False:
                            print("WARNING: no user assignment made also in TOCSY for resid", resid, "and Carbon", nucleus,". Skipping nucleus.")
                            assignment[7] += " <NOT ASSIGNED>"
                            #proofread_clean_resid_assignments_dict[resid].append(assignment)
                            continue
                    user_Creson = user_assignment[2]
                    user_Hreson = user_assignment[4]
                    print("DEBUG: user_Creson=", user_Creson, "user_Hreson=", user_Hreson)
                    if approx_equal(C_assignment[1], user_Creson, 0.4) and approx_equal(H_assignment[1], user_Hreson, 0.04): # if the user-made assignment agrees with ours, add to the comment "CORRECT"
                        print("DEBUG: CORRECT C_assignment found:", C_assignment)
                        print("DEBUG: CORRECT H_assignment found:", H_assignment)
                        comment = " <CORRECT>"
                    else:
                        print("DEBUG: WRONG C_assignment found:", C_assignment)
                        print("DEBUG: WRONG H_assignment found:", H_assignment)
                        comment = " <WRONG>"
                elif spectrum_combo in ['TOCSY (average)', 'TOCSY (unmatched)']:
                    user_assignment = get_Carbon_resonance(resid, nucleus, 'TOCSY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict, get_assignment=True)
                    if user_assignment == False:    # if no user assignment was made for this Carbon, add no extra comment
                        print("WARNING: no user assignment made in TOCSY for resid", resid, "and Carbon", nucleus,". Checking NOESY...")
                        user_assignment = get_Carbon_resonance(resid, nucleus, 'HCNH', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict, get_assignment=True)
                        if user_assignment == False:
                            print("WARNING: no user assignment made also in NOESY for resid", resid, "and Carbon", nucleus,". Skipping nucleus.")
                            assignment[7] += " <NOT ASSIGNED>"
                            #proofread_clean_resid_assignments_dict[resid].append(assignment)
                            continue
                    user_Creson = user_assignment[2]
                    user_Hreson = user_assignment[4]
                    if approx_equal(C_assignment[1], user_Creson, 0.4) and approx_equal(H_assignment[1], user_Hreson, 0.04): # if the user-made assignment agrees with ours, add to the comment "CORRECT"
                        print("DEBUG: CORRECT C_assignment found:", C_assignment)
                        print("DEBUG: CORRECT H_assignment found:", H_assignment)
                        comment = " <CORRECT>"
                    else:
                        print("DEBUG: WRONG C_assignment found:", C_assignment)
                        print("DEBUG: WRONG H_assignment found:", H_assignment)
                        comment = " <WRONG>"
                        
                if comment == "":
                    comment = " <NOT ASSIGNED>"
                C_assignment[7] += comment
                H_assignment[7] += comment
                comment = None  # a trick to add no more comments to "assignment" at the end of this block of code
            
            elif nucleus[0] == 'C':   # only if the assignment line belongs to a Carbon then check if the user has made the same assignment
                comment = ""
                if spectrum_combo in ['HCNH', 'NOESY (C-term hanging residue)', 'NOESY (flanked residue)', 'NOESY + TOCSY-HCNH', 'TOCSY-HCNH', 'TOCSY-HCNH + NOESY', 'TOCSY (unmatched) + TOCSY-HCNH', 'TOCSY-HCNH + TOCSY (unmatched)']:    # if this resonance comes only from NOESY
                    user_assignment = get_Carbon_resonance(resid, nucleus, 'HCNH', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict)
                    if user_assignment == False:    # if no user assignment was made for this Carbon, add no extra comment
                        print("WARNING: no user assignment made in NOESY for resid", resid, "and Carbon", nucleus,". Checking TOCSY...")
                        user_assignment = get_Carbon_resonance(resid, nucleus, 'TOCSY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict)
                        if user_assignment == False:
                            print("WARNING: no user assignment made also in TOCSY for resid", resid, "and Carbon", nucleus,". Skipping nucleus.")
                            assignment[7] += " <NOT ASSIGNED>"
                            # proofread_clean_resid_assignments_dict[resid].append(assignment)
                            continue
                    if approx_equal(assignment[1], user_assignment, 0.4): # if the user-made assignment agrees with ours, add to the comment "CORRECT"
                        print("DEBUG: CORRECT assignment found:", assignment)
                        comment = " <CORRECT>"
                    else:
                        print("DEBUG: WRONG assignment found:", assignment)
                        comment = " <WRONG>"
                elif spectrum_combo in ['TOCSY (average)', 'TOCSY (unmatched)']:
                    user_assignment = get_Carbon_resonance(resid, nucleus, 'TOCSY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict)
                    if user_assignment == False:    # if no user assignment was made for this Carbon, add no extra comment
                        print("WARNING: no user assignment made in TOCSY for resid", resid, "and Carbon", nucleus,". Checking NOESY...")
                        user_assignment = get_Carbon_resonance(resid, nucleus, 'HCNH', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict)
                        if user_assignment == False:
                            print("WARNING: no user assignment made also in NOESY for resid", resid, "and Carbon", nucleus,". Skipping nucleus.")
                            assignment[7] += " <NOT ASSIGNED>"
                            # proofread_clean_resid_assignments_dict[resid].append(assignment)
                            continue
                    if approx_equal(assignment[1], user_assignment, 0.4): # if the user-made assignment agrees with ours, add to the comment "CORRECT"
                        print("DEBUG: CORRECT assignment found:", assignment)
                        comment = " <CORRECT>"
                    else:
                        print("DEBUG: WRONG assignment found:", assignment)
                        comment = " <WRONG>"
                
            if comment == None:      # this is a proton or a methyl Carbon or a Carbon with only one proton
                continue
            elif comment == "" and not '<' in assignment[7]:    # if no comment has been appended in previous iterations
                comment = " <NOT ASSIGNED>"
            assignment[7] += comment
            # proofread_clean_resid_assignments_dict[resid].append(assignment)
    
    
    print("DEBUG: ITERATION 3 proofread_clean_resid_assignments_dict=", proofread_clean_resid_assignments_dict)
    
    ## WRITE THE PROOF-READ XEASY FILE
    xeasy_out = open(args.ABSOLUTE_MATCHES_FILE+".NOESY.proofread.iter3.xeasy", 'w')
    atom_index = 0  # reinitialize atom index
    for resid in list(proofread_clean_resid_assignments_dict.keys()):
        for peak in proofread_clean_resid_assignments_dict[resid]:
            atom_index += 1
            peak[0] = atom_index
            line = "\t".join([str(p) for p in peak])
            xeasy_out.write(line + "\n")
    xeasy_out.close()


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ITERATION 4 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
print("DEBUG: ENTERING ITERATION 4 matched_NOESY_residue_assignments_dict=", matched_NOESY_residue_assignments_dict)
print("DEBUG: ENTERING ITERATION 4 clean_resid_assignments_dict=", clean_resid_assignments_dict)

xeasy_fout = open(args.ABSOLUTE_MATCHES_FILE+".NOESY.iter4.xeasy", 'w')
sparky_fout = open(args.ABSOLUTE_MATCHES_FILE+".NOESY.iter4.sparky", 'w')
atom_index = 1
residues_with_written_prevAssigned_peaks = []

## THIS TIME ITERATE OVER ALL RESIDUE THAT HAVE NOESY AND ARE IN THE ALIGNMENT (INCLUDING FLANKED AND C-TERMINAL)
for i_aa_resid in revordered_residue_keys_iter3:    # both the keys and the C,H assignments correspond to residue i-1
    i_residue = i_aa_resid[0]+str(i_aa_resid[1])
    aa_type = aa1to3_dict[i_aa_resid[0]]
    i_resid = get_resid_from_residue(i_residue)
    # if not i_residue in ['I206']:
    #    continue
    # FIRST MATCH THE TOCSY ASSIGNMENTS OF RESIDUE i TO NOESY PEAKS
    #i_residue = get_i_residue_from_iminus1_residue(iminus1_residue)
    all_possible_assignments_list = []  # list with all the peaks for this residue (TOCSY-HCNH, NOESY) and the respective possible C-H assignments
    # TEMP DEACTIVATE: if i_residue in TOCSY_residue_assignments_dict.keys() and len(TOCSY_residue_assignments_dict[i_residue]) > 0:   # if no TOCSY assignments have been made for residue i, skip NOESY assignment of residue i
    print("ITERATION 4: Matching residue i ", i_residue," TOCSY PEAKS to residue's i NOESY peaks...")
    checked_residues_set.add(i_residue)
    tmp_missing_carbons_list, tmp_partly_assigned_carbons_list = get_missing_carbons_from_xeasy_dict(aa_type, i_resid, clean_resid_assignments_dict, matched_NOESY_residue_assignments_dict, allowed_aa_atoms_dict)
    missing_carbons_list = [c for c in tmp_missing_carbons_list if not (aa_type in list(iter1_allowed_atoms_dict.keys()) and c in iter1_allowed_atoms_dict[aa_type])]
    partly_assigned_carbons_list = [c for c in tmp_partly_assigned_carbons_list if not (aa_type in list(iter1_allowed_atoms_dict.keys()) and c in iter1_allowed_atoms_dict[aa_type])]
    unassigned_NOESY_peaks_list = []    # ini
    if len(missing_carbons_list) != 0:   # IF THERE ARE MISSING CARBONS IN ITERATION 1, TRY TO FIND THEM IN NOESY (EXCLUDE THE ASSIGNED OR MATCHED NOESY PEAKS)
        print("Residue ", i_residue, " was not completely assigned in iteration 1. Using NOESY to assigned resonances to the rest Carbon types.")
        try:    # DANGEROUS but works for the N-term residue (A200 in aLP)
            unassigned_NOESY_peaks_list = [peak for peak in matched_NOESY_residue_assignments_dict[i_residue] if peak[0]=="?" and peak[1]=="?" and peak[3]=="?"]
        except KeyError:    # if this residue is not in the alignment (not in matched_NOESY_residue_assignments_dict), then we cannot assign the partly assigned Carbons and hence we move on to the next
            print("Residue ", i_residue, " is either not in the alignment or doesn't have NOESY. Therefore it is impossible to assign partly assigned Carbon peaks.")
            continue
        for NOESY_peak in unassigned_NOESY_peaks_list:  # make C-H type predictions for each unassigned NOESY peak
            NOESY_Creson = NOESY_peak[2]
            NOESY_Hreson = NOESY_peak[4]
            all_possible_assignments_list.extend( get_probabilities_from_H_C_resonpair_2Dhist(i_residue,
                                                                                              NOESY_Hreson,
                                                                                              NOESY_Creson,
                                                                                              missing_carbons_list,
                                                                                              partly_assigned_carbons_list,
                                                                                              args,
                                                                                              NOESY_residue_peak_intensity_mdict,
                                                                                              matched_i_to_iplus1_peaks_mdict,
                                                                                              aa2pthres_dict,
                                                                                              aa2ithres_dict,
                                                                                              ali.protein_alignment_list,
                                                                                              ali.absolute_matches_alignment_list,
                                                                                              histload,
                                                                                              iteration=4)
                                                )
        [assignment.append('HCNH') for assignment in all_possible_assignments_list]    # we don't need to created a new list, the original all_possible_assignments_list will be updated
        print("DEBUG: i_residue=", i_residue, "all_possible_assignments_list=", all_possible_assignments_list)
        if i_resid in list(clean_resid_assignments_dict.keys()):  # if this is a C-term residue that is analyzed for the first time in iter4
            written_nucleiNames_list = [line[3]  for line in clean_resid_assignments_dict[i_resid]]   # list of nuclei that have been written already, to avoid duplication
        else:
            written_nucleiNames_list = []
        if len(all_possible_assignments_list) > 0:  # IF WE HAVE FOUND POSSIBLE NOESY PEAK ASSIGNMENTS OF THE MISSING CARBONS
            # WRITE THE BEST NOESY PEAK COMBINATION TO COMPLETE THE ASSIGNMENT OF THE MISSING CARBONS
            extra_written_nucleiNames_list, \
            matched_NOESY_residue_assignments_dict, \
            TOCSY_assigned_residue_list = \
                select_best_HCNH_peak_combination(i_residue,
                                                   all_possible_assignments_list,
                                                   partly_assigned_carbons_list,
                                                   atom_index,
                                                   xeasy_fout,
                                                   sparky_fout,
                                                   ali.protein_alignment_list,
                                                   ali.absolute_matches_alignment_list,
                                                   args,
                                                   NOESY_residue_NHresonances_dict,
                                                   matched_NOESY_residue_assignments_dict,
                                                   TOCSY_assigned_residue_list,
                                                   iter1_allowed_atoms_dict,
                                                   iteration=4)
            written_nucleiNames_list.extend(extra_written_nucleiNames_list)
        residues_with_written_prevAssigned_peaks.append(i_residue)
        # WRITE ALSO THE PEAKS ASSIGNED IN ITERATION 1 & 2
        write_peaks_assigned_in_previous_iteration(i_residue, clean_resid_assignments_dict, atom_index, xeasy_fout, 
                                               NOESY_residue_NHresonances_dict, residues_with_written_NH_list)
    elif len(missing_carbons_list) == 0 and len(partly_assigned_carbons_list) == 0: # NOTHING TO DO FOR THIS RESIDUE
        residues_with_written_prevAssigned_peaks.append(i_residue)
        # WRITE ALSO THE PEAKS ASSIGNED IN ITERATION 1 & 2
        write_peaks_assigned_in_previous_iteration(i_residue, clean_resid_assignments_dict, atom_index, xeasy_fout, 
                                               NOESY_residue_NHresonances_dict, residues_with_written_NH_list)
    else:   # IF ALL THE CARBONS WERE ASSIGNED IN ITERATION 1, SAVE THEM, BUT ALSO FIND IN NOESY AND SAVE PARTLY ASSIGNED CARBONS
        print("Residue ", i_residue, " was completely assigned in iteration 1.")
        if not i_residue in list(matched_NOESY_residue_assignments_dict.keys()):
            print("Residue ", i_residue, " is either not in the alignment or doesn't have NOESY. Therefore it is impossible to assign partly assigned Carbon peaks.")
            continue
        try:
            aa_type = aa1to3_dict[get_aa_type_from_residue(i_residue)]
        except (IndexError, KeyError):    # skip residue i
            print("EXITING: COULD NOT FIND AA_TYPE OF RESIDUE", i_residue)
            sys.exit(1)
        assigned_peaks_list = get_peaks_from_xeasy(i_resid, clean_resid_assignments_dict)     # get all assigned peaks of this residue from xeasy format
        for assigned_peak in assigned_peaks_list:
            spectrum_combo = assigned_peak[5]
            Cname = assigned_peak[1]
            Creson = assigned_peak[2]    # use the TOCSY resonance by default
            Hname = assigned_peak[3]
            Hreson = assigned_peak[4]
            all_possible_assignments_list.append([aa_type, Cname, Hname, 10000000.0, Hreson, Creson, spectrum_combo])  # set the probability to 1 cause it has been already assigned in iteration 1
            print("DEBUG: i_residue=", i_residue, "aa_type=", aa_type, "all_possible_assignments_list=", all_possible_assignments_list)
        written_nucleiNames_list = []
        if i_resid in list(clean_resid_assignments_dict.keys()):
            written_nucleiNames_list = [line[3]  for line in clean_resid_assignments_dict[i_resid]]   # list of nuclei that have been written already, to avoid duplication
        if len(all_possible_assignments_list) > 0:  # IF WE HAVE FOUND POSSIBLE NOESY PEAK ASSIGNMENTS OF THE MISSING CARBONS
            # WRITE THE BEST NOESY PEAK COMBINATION TO COMPLETE THE ASSIGNMENT OF THE MISSING CARBONS
            extra_written_nucleiNames_list, \
            matched_NOESY_residue_assignments_dict, \
            TOCSY_assigned_residue_list = \
                select_best_HCNH_peak_combination(i_residue,
                                                   all_possible_assignments_list,
                                                   partly_assigned_carbons_list,
                                                   atom_index,
                                                   xeasy_fout,
                                                   sparky_fout,
                                                   ali.protein_alignment_list,
                                                   ali.absolute_matches_alignment_list,
                                                   args,
                                                   NOESY_residue_NHresonances_dict,
                                                   matched_NOESY_residue_assignments_dict,
                                                   TOCSY_assigned_residue_list,
                                                   iter1_allowed_atoms_dict,
                                                   iteration=4)
            written_nucleiNames_list.extend(extra_written_nucleiNames_list)

# WRITE ALSO THE PEAKS ASSIGNED IN PREVIOUS ITERATIONS OF ALL RESIDUES THAT WERE EXCLUDED FROM THIS ITERATION
for i_aa_resid in set(revordered_residue_keys_iter1 + revordered_residue_keys_iter2 + revordered_residue_keys_iter3 + revordered_flanked_residue_keys):   # there is no revordered_residue_keys_iter4
    i_residue = i_aa_resid[0]+str(i_aa_resid[1])
    aa_type = aa1to3_dict[i_aa_resid[0]]
    if not i_residue in residues_with_written_prevAssigned_peaks:
        write_peaks_assigned_in_previous_iteration(i_residue, clean_resid_assignments_dict, atom_index, xeasy_fout, 
                                               NOESY_residue_NHresonances_dict, residues_with_written_NH_list)

# # print "DEBUG: before entering write_flanked_residues(), ali.absolute_matches_alignment_list=", ali.absolute_matches_alignment_list
# write_flanked_residues(iteration=4)

# Now find if there are any residue in the alignment that haven't been checked:
for residue in ali.absolute_matches_alignment_list:
    if not residue in checked_residues_set and residue != '-':
        print("WARNING: Residue ", residue, " was not checked!!!")
        

# ##
# ##  COPY FROM TOCSY ASSIGNMENTS (TOCSY_residue_assignments_dict) OF ALL THE RESONANCES THAT COULD NOT BE MATCHED WITH THE NOESY ASSIGNMENTS (TOCSY_residue_assignments_dict)
# ##
# write_unmatched_TOCSY_peaks(iteration=4)


# ##
# ##  WRITE THE HN & N RESONANCES TO XEASY ONLY FOR THOSE RESIDUES THAT NO C,H ASSIGNMENTS COULD BE MADE BUT THEIR HN & H RESONANCES ARE KNOWN FROM TOCSY OR NOESY
# ##
write_NH_of_residues_without_CH(ali.absolute_matches_alignment_list, NOESY_residue_NHresonances_dict, TOCSY_residue_NHresonances_dict,
                                    residues_with_written_NH_list, atom_index, xeasy_fout, sparky_fout, ali.protein_alignment_list, args,
                                    Cterm_residue_set, iteration=4)

xeasy_fout.close()
sparky_fout.close()

##
##  REWRITE THE XEASY FILE BY REMOVING DUPLICATES OR AVERAGING, AND CORRECTING THE ORDER
##

## CLEAN THE CONTENTS OF XEASY FILE
clean_resid_assignments_dict = clean_xeasy_file(args.ABSOLUTE_MATCHES_FILE+".NOESY.iter4.xeasy")  # dict with all the assignments from each residue from iteration 1
# sort the peak assignments of each resid
for resid, peak_list in list(clean_resid_assignments_dict.items()):
    peak_list.sort(key=itemgetter(3))
    clean_resid_assignments_dict[resid] = peak_list

print("DEBUG: ITERATION 4 point 1 clean_resid_assignments_dict=", clean_resid_assignments_dict)

## WRITE THE CLEAN XEASY FILE
xeasy_out = open(args.ABSOLUTE_MATCHES_FILE+".NOESY.iter4.xeasy", 'w')
atom_index = 0  # reinitialize atom index
for resid in list(clean_resid_assignments_dict.keys()):
    for peak in clean_resid_assignments_dict[resid]:
        atom_index += 1
        peak[0] = atom_index
        line = "\t".join([str(p) for p in peak])
        xeasy_out.write(line + "\n")
xeasy_out.close()


##
##  REWRITE THE XEASY FILE BY REMOVING DUPLICATES OR AVERAGING, AND CORRECTING THE ORDER
##

## CLEAN THE CONTENTS OF XEASY FILE
# clean_resid_assignments_dict contains the final, non-redundant assignments written to xeasy file. E.g.
# 1: [[2336, 56.022, 0.2, 'CA', 1, '#', 'MET', 'TOCSY (unmatched) (iter1)'], [2338, 34.786, 0.2, 'CB', 1, '#', 'MET', 'TOCSY (unmatched) (iter1)'],
# [2340, 30.946, 0.2, 'CG', 1, '#', 'MET', 'TOCSY (unmatched) (iter1)'], [2337, 4.194, 0.02, 'HA', 1, '#', 'MET', 'TOCSY (unmatched) (iter1)'],
# [2339, 2.14, 0.02, 'QB', 1, '#', 'MET', 'TOCSY (unmatched) (iter1)'], [2341, 2.487, 0.02, 'QG', 1, '#', 'MET', 'TOCSY (unmatched) (iter1)']]
clean_resid_assignments_dict = clean_xeasy_file(args.ABSOLUTE_MATCHES_FILE+".NOESY.iter4.xeasy")  
# sort the peak assignments of each resid
for resid, peak_list in list(clean_resid_assignments_dict.items()):
    peak_list.sort(key=itemgetter(3))
    clean_resid_assignments_dict[resid] = peak_list

print("DEBUG: ITERATION 4 point 2 clean_resid_assignments_dict=", clean_resid_assignments_dict)

## WRITE THE CLEAN XEASY FILE
xeasy_out = open(args.ABSOLUTE_MATCHES_FILE+".NOESY.cleaned.iter4.xeasy", 'w')
atom_index = 0  # reinitialize atom index
for resid in list(clean_resid_assignments_dict.keys()):
    for peak in clean_resid_assignments_dict[resid]:
        atom_index += 1
        peak[0] = atom_index
        line = "\t".join([str(p) for p in peak])
        xeasy_out.write(line + "\n")
xeasy_out.close()


##
##  IF PROVIDED, LOAD THE USER-ANNOTATED TOCSY AND NOESY FILES FOR PROOF-READING
##
print("DEBUG: args.user_TOCSY_FILE=", args.user_TOCSY_FILE, "args.user_NOESY_FILE=", args.user_NOESY_FILE)
if args.user_TOCSY_FILE and args.user_NOESY_FILE:
    print("STARTING PROOF-READING NOESY CHEMICAL SHIFT ASSIGNMENTS")
    # THIS IS ITERATION 4, SO THERE ARE SOME EXTRA LINES THAT DON'T EXIST WHICH NEED TO BE ADDED
    for resid in clean_resid_assignments_dict:
        if not resid in list(proofread_clean_resid_assignments_dict.keys()):
            proofread_clean_resid_assignments_dict[resid] = clean_resid_assignments_dict[resid]
            continue
        for new_assignment in clean_resid_assignments_dict[resid]:
            new_reson = new_assignment[1]
            old_reson_list = [a[1] for a in proofread_clean_resid_assignments_dict[resid]]
            if not new_reson in old_reson_list:     # if this resonance was not found (newly assigned), add it to the dictionary
                proofread_clean_resid_assignments_dict[resid].append(new_assignment)
    print("DEBUG: ITERATION 4 point 3 clean_resid_assignments_dict=", clean_resid_assignments_dict)
    
    ##
    ##  COMPARE YOUR CARBON ASSIGNMENTS WITH THE USER ASSIGNMENTS
    ##
    ## Remember that the comments can be the following:
    # NOESY
    # NOESY (C-term hanging residue)
    # NOESY (flanked residue)
    # NOESY + TOCSY-HCNH
    # TOCSY (average)
    # TOCSY-HCNH
    # TOCSY-HCNH + NOESY
    # TOCSY-HCNH + TOCSY (unmatched)
    # TOCSY (unmatched)
    # TOCSY (unmatched) + TOCSY-HCNH
    # Recall that if we have 2 peaks in the same C-group but only one comes from NOESY ('TOCSY (unmatched) + TOCSY-HCNH', or
    # 'TOCSY-HCNH + TOCSY (unmatched)'), then we save only the NOESY Carbon resonance. 'TOCSY (unmatched)' -> only one resonance available
    # only in TOCSY, 'TOCSY (average)' -> 2 resonances available only from TOCSY and they were averaged
    
    
    for resid in list(proofread_clean_resid_assignments_dict.keys()):
        # FIRST TAKE CARE OF LEU AND VAL EQUIVALENT METHYL CARBONS
        print("DEBUG ITERATION4: proof reading resid=", resid)
        resname = proofread_clean_resid_assignments_dict[resid][0][6]
        if resname in list(aa_equivalentCarbons_dict.keys()): # If this is a LEU or VAL
            Carbon_pair = aa_equivalentCarbons_dict[resname][0]
            if len([a for a in proofread_clean_resid_assignments_dict[resid] if a[3] in Carbon_pair]) > 0: # AND if the methyl carbons have been assigned
                Proton_pair = ["H"+C[1:] for C in Carbon_pair]   # protons of the equivalent Carbons
                # get only the methyl Carbon assignments (LEU --> CD1,CD2 or VAL --> CG1,CG2) made by the program, to compare them with the user made assignments
                program_assignments_list = [assignment for assignment in proofread_clean_resid_assignments_dict[resid] if assignment[3] in Carbon_pair+Proton_pair]   # keep only the resonances assigned by the program to the equivalent methyl carbons and protons
                # proof_read_equivalentCarbons(program_assignments_list, resid, user_VAL_LEU_NOESY_resid_assignments_dict, user_VAL_LEU_TOCSY_resid_assignments_dict, Carbon_pair)   # OBSOLETE!!! it doesn't matter if 1 or 2 of the peaks have been proof-read in iteration 1 & 2
                proof_read_equivalentMethyls(program_assignments_list, resid, user_VAL_LEU_NOESY_resid_assignments_dict, user_VAL_LEU_TOCSY_resid_assignments_dict, Carbon_pair, Proton_pair)   # it doesn't matter if 1 or 2 of the peaks have been proof-read in iteration 1 & 2
        # THEN PROCEED NORMALLY
        for assignment in proofread_clean_resid_assignments_dict[resid]:
            if '<' in assignment[7]:    # if this assignment has been proofread, then skip it
                continue
            print("DEBUG ITERATION4: proof reading new assignment=", assignment)
            spectrum_combo = re.sub(r' iter[1-9]', '', assignment[7])    # remove the "iter[1-9]" suffix to be able to distinguish the type of assignment (NEESY, NOESY+TOCSY, etc.)
            comment = None
            aa_type = assignment[6]
            nucleus = assignment[3] # the nucleus name
            if aa_type in list(aa_equivalentCarbons_dict.keys()) and nucleus in aa_equivalentCarbons_dict[aa_type][0]:    # if this peak has been assigned by the program to LEU or VAL CD1, CD2 or CG1,CG2 respectively, skip it
                continue
            
            # IF THIS CARBON IS A METHYL OR HAS ONLY ONE PROTON, CONSIDE ALSO THE PROTON RESONANCE FOR PROOFREADING
            if aa_type in list(aatype_carbon_nongeminalHname_mdict.keys()) and nucleus in list(aatype_carbon_nongeminalHname_mdict[aa_type].keys()):
                # get the proton assignment made by the program
                C_assignment = assignment # the program-made carbon assignment
                H_assignment = None # default value for the program-made proton assignment
                print("DEBUG: resid=", resid, "proofread_clean_resid_assignments_dict[resid]=", proofread_clean_resid_assignments_dict[resid])
                for a in proofread_clean_resid_assignments_dict[resid]:
                    if ((type(aatype_carbon_nongeminalHname_mdict[aa_type][nucleus]) == list and a[3] in aatype_carbon_nongeminalHname_mdict[aa_type][nucleus]) or
                        (type(aatype_carbon_nongeminalHname_mdict[aa_type][nucleus]) == str and a[3] == aatype_carbon_nongeminalHname_mdict[aa_type][nucleus])):
                        H_assignment = a
                        break
                print("DEBUG: C_assignment=", C_assignment, "H_assignment=", H_assignment)
                comment = ""
                if spectrum_combo in ['HCNH', 'NOESY (C-term hanging residue)', 'NOESY (flanked residue)', 'NOESY + TOCSY-HCNH', 'TOCSY-HCNH', 'TOCSY-HCNH + NOESY', 'TOCSY (unmatched) + TOCSY-HCNH', 'TOCSY-HCNH + TOCSY (unmatched)']:    # if this resonance comes only from NOESY
                    user_assignment = get_Carbon_resonance(resid, nucleus, 'HCNH', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict, get_assignment=True)
                    if user_assignment == False:    # if no user assignment was made for this Carbon, add no extra comment
                        print("WARNING: no user assignment made in NOESY for resid", resid, "and Carbon", nucleus,". Checking TOCSY...")
                        user_assignment = get_Carbon_resonance(resid, nucleus, 'TOCSY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict, get_assignment=True)
                        if user_assignment == False:
                            print("WARNING: no user assignment made also in TOCSY for resid", resid, "and Carbon", nucleus,". Skipping nucleus.")
                            assignment[7] += " <NOT ASSIGNED>"
                            #proofread_clean_resid_assignments_dict[resid].append(assignment)
                            continue
                    user_Creson = user_assignment[2]
                    user_Hreson = user_assignment[4]
                    print("DEBUG: user_Creson=", user_Creson, "user_Hreson=", user_Hreson)
                    if approx_equal(C_assignment[1], user_Creson, 0.4) and approx_equal(H_assignment[1], user_Hreson, 0.04): # if the user-made assignment agrees with ours, add to the comment "CORRECT"
                        print("DEBUG: CORRECT C_assignment found:", C_assignment)
                        print("DEBUG: CORRECT H_assignment found:", H_assignment)
                        comment = " <CORRECT>"
                    else:
                        print("DEBUG: WRONG C_assignment found:", C_assignment)
                        print("DEBUG: WRONG H_assignment found:", H_assignment)
                        comment = " <WRONG>"
                elif spectrum_combo in ['TOCSY (average)', 'TOCSY (unmatched)']:
                    user_assignment = get_Carbon_resonance(resid, nucleus, 'TOCSY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict, get_assignment=True)
                    if user_assignment == False:    # if no user assignment was made for this Carbon, add no extra comment
                        print("WARNING: no user assignment made in TOCSY for resid", resid, "and Carbon", nucleus,". Checking NOESY...")
                        user_assignment = get_Carbon_resonance(resid, nucleus, 'HCNH', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict, get_assignment=True)
                        if user_assignment == False:
                            print("WARNING: no user assignment made also in NOESY for resid", resid, "and Carbon", nucleus,". Skipping nucleus.")
                            assignment[7] += " <NOT ASSIGNED>"
                            #proofread_clean_resid_assignments_dict[resid].append(assignment)
                            continue
                    user_Creson = user_assignment[2]
                    user_Hreson = user_assignment[4]
                    if approx_equal(C_assignment[1], user_Creson, 0.4) and approx_equal(H_assignment[1], user_Hreson, 0.04): # if the user-made assignment agrees with ours, add to the comment "CORRECT"
                        print("DEBUG: CORRECT C_assignment found:", C_assignment)
                        print("DEBUG: CORRECT H_assignment found:", H_assignment)
                        comment = " <CORRECT>"
                    else:
                        print("DEBUG: WRONG C_assignment found:", C_assignment)
                        print("DEBUG: WRONG H_assignment found:", H_assignment)
                        comment = " <WRONG>"
                        
                if comment == "":
                    comment = " <NOT ASSIGNED>"
                C_assignment[7] += comment
                H_assignment[7] += comment
                comment = None  # a trick to add no more comments to "assignment" at the end of this block of code
            
            elif nucleus[0] == 'C':   # only if the assignment line belongs to a Carbon then check if the user has made the same assignment
                comment = ""
                if spectrum_combo in ['HCNH', 'NOESY (C-term hanging residue)', 'NOESY (flanked residue)', 'NOESY + TOCSY-HCNH', 'TOCSY-HCNH', 'TOCSY-HCNH + NOESY', 'TOCSY (unmatched) + TOCSY-HCNH', 'TOCSY-HCNH + TOCSY (unmatched)']:    # if this resonance comes only from NOESY
                    user_assignment = get_Carbon_resonance(resid, nucleus, 'HCNH', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict)
                    if user_assignment == False:    # if no user assignment was made for this Carbon, add no extra comment
                        print("WARNING: no user assignment made in NOESY for resid", resid, "and Carbon", nucleus,". Checking TOCSY...")
                        user_assignment = get_Carbon_resonance(resid, nucleus, 'TOCSY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict)
                        if user_assignment == False:
                            print("WARNING: no user assignment made also in TOCSY for resid", resid, "and Carbon", nucleus,". Skipping nucleus.")
                            assignment[7] += " <NOT ASSIGNED>"
                            # proofread_clean_resid_assignments_dict[resid].append(assignment)
                            continue
                    if approx_equal(assignment[1], user_assignment, 0.4): # if the user-made assignment agrees with ours, add to the comment "CORRECT"
                        print("DEBUG: CORRECT assignment found:", assignment)
                        comment = " <CORRECT>"
                    else:
                        print("DEBUG: WRONG assignment found:", assignment)
                        comment = " <WRONG>"
                elif spectrum_combo in ['TOCSY (average)', 'TOCSY (unmatched)']:
                    user_assignment = get_Carbon_resonance(resid, nucleus, 'TOCSY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict)
                    if user_assignment == False:    # if no user assignment was made for this Carbon, add no extra comment
                        print("WARNING: no user assignment made in TOCSY for resid", resid, "and Carbon", nucleus,". Checking NOESY...")
                        user_assignment = get_Carbon_resonance(resid, nucleus, 'HCNH', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict)
                        if user_assignment == False:
                            print("WARNING: no user assignment made also in NOESY for resid", resid, "and Carbon", nucleus,". Skipping nucleus.")
                            assignment[7] += " <NOT ASSIGNED>"
                            # proofread_clean_resid_assignments_dict[resid].append(assignment)
                            continue
                    if approx_equal(assignment[1], user_assignment, 0.4): # if the user-made assignment agrees with ours, add to the comment "CORRECT"
                        print("DEBUG: CORRECT assignment found:", assignment)
                        comment = " <CORRECT>"
                    else:
                        print("DEBUG: WRONG assignment found:", assignment)
                        comment = " <WRONG>"
                
            if comment == None:      # this is a proton or a methyl Carbon or a Carbon with only one proton
                continue
            elif comment == "" and not '<' in assignment[7]:    # if no comment has been appended in previous iterations
                comment = " <NOT ASSIGNED>"
            assignment[7] += comment
    
    
    print("DEBUG: ITERATION 4 proofread_clean_resid_assignments_dict=", proofread_clean_resid_assignments_dict)
    
    ## WRITE THE PROOF-READ XEASY FILE
    xeasy_out = open(args.ABSOLUTE_MATCHES_FILE+".NOESY.proofread.iter4.xeasy", 'w')
    atom_index = 0  # reinitialize atom index
    for resid in list(proofread_clean_resid_assignments_dict.keys()):
        for peak in proofread_clean_resid_assignments_dict[resid]:
            atom_index += 1
            peak[0] = atom_index
            line = "\t".join([str(p) for p in peak])
            xeasy_out.write(line + "\n")
    xeasy_out.close()


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ITERATION 5 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## In iteration 5 only PHE-->CD-QD and TYR-->CD-QD,CE-QE are assigned
print("## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ITERATION 5 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##")
print("DEBUG: ENTERING ITERATION 5 matched_NOESY_residue_assignments_dict=", matched_NOESY_residue_assignments_dict)
print("DEBUG: ENTERING ITERATION 5 clean_resid_assignments_dict=", clean_resid_assignments_dict)

allowed_aa_atoms_dict["PHE"].extend(["HD1", "CD1"])
allowed_aa_atoms_dict["TYR"].extend(["HD1", "HE1", "CD1", "CE1"])


xeasy_fout = open(args.ABSOLUTE_MATCHES_FILE+".NOESY.iter5.xeasy", 'w')
sparky_fout = open(args.ABSOLUTE_MATCHES_FILE+".NOESY.iter5.sparky", 'w')
atom_index = 1
residues_with_written_prevAssigned_peaks = []


## ITERATE OVER ALL RESIDUE THAT HAVE NOESY AND ARE IN THE ALIGNMENT (INCLUDING FLANKED AND C-TERMINAL)
for i_aa_resid in revordered_residue_keys_iter3:    # both the keys and the C,H assignments correspond to residue i-1
    i_residue = i_aa_resid[0]+str(i_aa_resid[1])
    aa_type = aa1to3_dict[i_aa_resid[0]]
    i_resid = get_resid_from_residue(i_residue)
    if not aa_type in ['TYR', 'PHE', 'LYS', 'ARG']:
        continue
    # FIRST MATCH THE TOCSY ASSIGNMENTS OF RESIDUE i TO NOESY PEAKS
    #i_residue = get_i_residue_from_iminus1_residue(iminus1_residue)
    all_possible_assignments_list = []  # list with all the peaks for this residue (TOCSY-HCNH, NOESY) and the respective possible C-H assignments
    # TEMP DEACTIVATE if i_residue in TOCSY_residue_assignments_dict.keys() and len(TOCSY_residue_assignments_dict[i_residue]) > 0:   # if no TOCSY assignments have been made for residue i, skip NOESY assignment of residue i
    print("ITERATION 5: Matching residue i ", i_residue," TOCSY PEAKS to residue's i NOESY peaks...")
    checked_residues_set.add(i_residue)
    tmp_missing_carbons_list, tmp_partly_assigned_carbons_list = get_missing_carbons_from_xeasy_dict(aa_type, i_resid, clean_resid_assignments_dict, matched_NOESY_residue_assignments_dict, allowed_aa_atoms_dict)
    missing_carbons_list = [c for c in tmp_missing_carbons_list if not (aa_type in list(iter1_allowed_atoms_dict.keys()) and c in iter1_allowed_atoms_dict[aa_type])]
    partly_assigned_carbons_list = [c for c in tmp_partly_assigned_carbons_list if not (aa_type in list(iter1_allowed_atoms_dict.keys()) and c in iter1_allowed_atoms_dict[aa_type])]
    unassigned_NOESY_peaks_list = []    # ini
    if len(missing_carbons_list) != 0:   # IF THERE ARE MISSING CARBONS IN ITERATION 1 & 2, TRY TO FIND THEM IN NOESY (EXCLUDE THE ASSIGNED OR MATCHED NOESY PEAKS)
        print("Residue ", i_residue, " was not completely assigned in iteration 1 & 2. Using NOESY to assigned resonances to the rest Carbon types.")
        try:    # DANGEROUS but works for the N-term residue (A200 in aLP)
            unassigned_NOESY_peaks_list = [peak for peak in matched_NOESY_residue_assignments_dict[i_residue] if peak[0]=="?" and peak[1]=="?" and peak[3]=="?"]
        except KeyError:    # if this residue is not in the alignment (not in matched_NOESY_residue_assignments_dict), then we cannot assign the partly assigned Carbons and hence we move on to the next
            print("Residue ", i_residue, " is either not in the alignment or doesn't have NOESY. Therefore it is impossible to assign partly assigned Carbon peaks.")
            continue
        for NOESY_peak in unassigned_NOESY_peaks_list:  # make C-H type predictions for each unassigned NOESY peak
            NOESY_Creson = NOESY_peak[2]
            NOESY_Hreson = NOESY_peak[4]
            print("DEBUG: aa_type",aa_type,"NOESY_Hreson",NOESY_Hreson, "NOESY_Creson", NOESY_Creson, "missing_carbons_list", missing_carbons_list)
            all_possible_assignments_list.extend( get_probabilities_from_H_C_resonpair_2Dhist(i_residue,
                                                                                              NOESY_Hreson,
                                                                                              NOESY_Creson,
                                                                                              missing_carbons_list,
                                                                                              partly_assigned_carbons_list,
                                                                                              args,
                                                                                              NOESY_residue_peak_intensity_mdict,
                                                                                              matched_i_to_iplus1_peaks_mdict,
                                                                                              aa2pthres_dict,
                                                                                              aa2ithres_dict,
                                                                                              ali.protein_alignment_list,
                                                                                              ali.absolute_matches_alignment_list,
                                                                                              histload,
                                                                                              iteration=5)
                                                )
        [assignment.append('HCNH') for assignment in all_possible_assignments_list]    # we don't need to created a new list, the original all_possible_assignments_list will be updated
        print("DEBUG: i_residue=", i_residue, "all_possible_assignments_list=", all_possible_assignments_list)
        written_nucleiNames_list = []
        if i_resid in list(clean_resid_assignments_dict.keys()):
            written_nucleiNames_list = [line[3]  for line in clean_resid_assignments_dict[i_resid]] # list of nuclei that have been written already, to avoid duplication
        if len(all_possible_assignments_list) > 0:  # IF WE HAVE FOUND POSSIBLE NOESY PEAK ASSIGNMENTS OF THE MISSING CARBONS
            # WRITE THE BEST NOESY PEAK COMBINATION TO COMPLETE THE ASSIGNMENT OF THE MISSING CARBONS
            extra_written_nucleiNames_list, \
            matched_NOESY_residue_assignments_dict, \
            TOCSY_assigned_residue_list = \
                select_best_HCNH_peak_combination(i_residue,
                                                   all_possible_assignments_list,
                                                   partly_assigned_carbons_list,
                                                   atom_index,
                                                   xeasy_fout,
                                                   sparky_fout,
                                                   ali.protein_alignment_list,
                                                   ali.absolute_matches_alignment_list,
                                                   args,
                                                   NOESY_residue_NHresonances_dict,
                                                   matched_NOESY_residue_assignments_dict,
                                                   TOCSY_assigned_residue_list,
                                                   iter1_allowed_atoms_dict,
                                                   iteration=5)
            written_nucleiNames_list.extend(extra_written_nucleiNames_list)
        residues_with_written_prevAssigned_peaks.append(i_residue)
        # WRITE ALSO THE PEAKS ASSIGNED IN ITERATION 1 & 2 & 3 & 4
        write_peaks_assigned_in_previous_iteration(i_residue, clean_resid_assignments_dict, atom_index, xeasy_fout, 
                                               NOESY_residue_NHresonances_dict, residues_with_written_NH_list)
    elif len(missing_carbons_list) == 0 and len(partly_assigned_carbons_list) == 0: # NOTHING TO DO FOR THIS RESIDUE
        residues_with_written_prevAssigned_peaks.append(i_residue)
        # WRITE ALSO THE PEAKS ASSIGNED IN ITERATION 1 & 2
        write_peaks_assigned_in_previous_iteration(i_residue, clean_resid_assignments_dict, atom_index, xeasy_fout, 
                                               NOESY_residue_NHresonances_dict, residues_with_written_NH_list)
    else:   # IF ALL THE CARBONS WERE ASSIGNED IN PREVIOUS ITERATIONS, SAVE THEM, BUT ALSO FIND IN NOESY AND SAVE PARTLY ASSIGNED CARBONS
        print("Residue ", i_residue, " was completely assigned in iteration 1 & 2.")
        if not i_residue in list(matched_NOESY_residue_assignments_dict.keys()):
            print("Residue ", i_residue, " is either not in the alignment or doesn't have NOESY. Therefore it is impossible to assign partly assigned Carbon peaks.")
            continue
        try:
            aa_type = aa1to3_dict[get_aa_type_from_residue(i_residue)]
        except (IndexError, KeyError):    # skip residue i
            print("EXITING: COULD NOT FIND AA_TYPE OF RESIDUE", i_residue)
            sys.exit(1)
        assigned_peaks_list = get_peaks_from_xeasy(i_resid, clean_resid_assignments_dict)     # get all assigned peaks of this residue from xeasy format
        for assigned_peak in assigned_peaks_list:
            spectrum_combo = assigned_peak[5]
            Cname = assigned_peak[1]
            Creson = assigned_peak[2]    # use the TOCSY resonance by default
            Hname = assigned_peak[3]
            Hreson = assigned_peak[4]
            all_possible_assignments_list.append([aa_type, Cname, Hname, 10000000.0, Hreson, Creson, spectrum_combo])  # set the probability to 1 cause it has been already assigned in iteration 1 & 2
            print("DEBUG: i_residue=", i_residue, "aa_type=", aa_type, "all_possible_assignments_list=", all_possible_assignments_list)
        written_nucleiNames_list = []
        if i_resid in list(clean_resid_assignments_dict.keys()):
            written_nucleiNames_list = [line[3]  for line in clean_resid_assignments_dict[i_resid]]  # list of nuclei that have been written already, to avoid duplication
        if len(all_possible_assignments_list) > 0:  # IF WE HAVE FOUND POSSIBLE NOESY PEAK ASSIGNMENTS OF THE MISSING CARBONS
            # WRITE THE BEST NOESY PEAK COMBINATION TO COMPLETE THE ASSIGNMENT OF THE MISSING CARBONS
            extra_written_nucleiNames_list, \
            matched_NOESY_residue_assignments_dict, \
            TOCSY_assigned_residue_list = \
                select_best_HCNH_peak_combination(i_residue,
                                                   all_possible_assignments_list,
                                                   partly_assigned_carbons_list,
                                                   atom_index,
                                                   xeasy_fout,
                                                   sparky_fout,
                                                   ali.protein_alignment_list,
                                                   ali.absolute_matches_alignment_list,
                                                   args,
                                                   NOESY_residue_NHresonances_dict,
                                                   matched_NOESY_residue_assignments_dict,
                                                   TOCSY_assigned_residue_list,
                                                   iter1_allowed_atoms_dict,
                                                   iteration=5)
            written_nucleiNames_list.extend(extra_written_nucleiNames_list)

# WRITE ALSO THE PEAKS ASSIGNED IN PREVIOUS ITERATIONS OF ALL RESIDUES THAT WERE EXCLUDED FROM THIS ITERATION
for i_aa_resid in set(revordered_residue_keys_iter1+revordered_residue_keys_iter2+revordered_residue_keys_iter3 + revordered_flanked_residue_keys):   # there is no revordered_residue_keys_iter5
    i_residue = i_aa_resid[0]+str(i_aa_resid[1])
    aa_type = aa1to3_dict[i_aa_resid[0]]
    if not i_residue in residues_with_written_prevAssigned_peaks:
        write_peaks_assigned_in_previous_iteration(i_residue, clean_resid_assignments_dict, atom_index, xeasy_fout, 
                                               NOESY_residue_NHresonances_dict, residues_with_written_NH_list)


# Now find if there are any residue in the alignment that haven't been checked:
for residue in ali.absolute_matches_alignment_list:
    if not residue in checked_residues_set and residue != '-':
        print("WARNING: Residue ", residue, " was not checked!!!")
        

# ##
# ##  COPY FROM TOCSY ASSIGNMENTS (TOCSY_residue_assignments_dict) OF ALL THE RESONANCES THAT COULD NOT BE MATCHED WITH THE NOESY ASSIGNMENTS (TOCSY_residue_assignments_dict)
# ##
# write_unmatched_TOCSY_peaks(iteration=5)


# ##
# ##  WRITE THE HN & N RESONANCES TO XEASY ONLY FOR THOSE RESIDUES THAT NO C,H ASSIGNMENTS COULD BE MADE BUT THEIR HN & H RESONANCES ARE KNOWN FROM TOCSY OR NOESY
# ##
write_NH_of_residues_without_CH(ali.absolute_matches_alignment_list, NOESY_residue_NHresonances_dict, TOCSY_residue_NHresonances_dict,
                                    residues_with_written_NH_list, atom_index, xeasy_fout, sparky_fout, ali.protein_alignment_list, args,
                                    Cterm_residue_set, iteration=5)

xeasy_fout.close()
sparky_fout.close()

##
##  REWRITE THE XEASY FILE BY REMOVING DUPLICATES OR AVERAGING, AND CORRECTING THE ORDER
##

## CLEAN THE CONTENTS OF XEASY FILE
clean_resid_assignments_dict = clean_xeasy_file(args.ABSOLUTE_MATCHES_FILE+".NOESY.iter5.xeasy")  # dict with all the assignments from each residue from iteration 1 & 2
# sort the peak assignments of each resid
for resid, peak_list in list(clean_resid_assignments_dict.items()):
    peak_list.sort(key=itemgetter(3))
    clean_resid_assignments_dict[resid] = peak_list

print("DEBUG: ITERATION 5 point 1 clean_resid_assignments_dict=", clean_resid_assignments_dict)

## WRITE THE CLEAN XEASY FILE
xeasy_out = open(args.ABSOLUTE_MATCHES_FILE+".NOESY.iter5.xeasy", 'w')
atom_index = 0  # reinitialize atom index
for resid in list(clean_resid_assignments_dict.keys()):
    for peak in clean_resid_assignments_dict[resid]:
        atom_index += 1
        peak[0] = atom_index
        line = "\t".join([str(p) for p in peak])
        xeasy_out.write(line + "\n")
xeasy_out.close()


##
##  REWRITE THE XEASY FILE BY REMOVING DUPLICATES OR AVERAGING, AND CORRECTING THE ORDER
##

## CLEAN THE CONTENTS OF XEASY FILE
# clean_resid_assignments_dict contains the final, non-redundant assignments written to xeasy file. E.g.
# 1: [[2336, 56.022, 0.2, 'CA', 1, '#', 'MET', 'TOCSY (unmatched) (iter1)'], [2338, 34.786, 0.2, 'CB', 1, '#', 'MET', 'TOCSY (unmatched) (iter1)'],
# [2340, 30.946, 0.2, 'CG', 1, '#', 'MET', 'TOCSY (unmatched) (iter1)'], [2337, 4.194, 0.02, 'HA', 1, '#', 'MET', 'TOCSY (unmatched) (iter1)'],
# [2339, 2.14, 0.02, 'QB', 1, '#', 'MET', 'TOCSY (unmatched) (iter1)'], [2341, 2.487, 0.02, 'QG', 1, '#', 'MET', 'TOCSY (unmatched) (iter1)']]
clean_resid_assignments_dict = clean_xeasy_file(args.ABSOLUTE_MATCHES_FILE+".NOESY.iter5.xeasy")  
# sort the peak assignments of each resid
for resid, peak_list in list(clean_resid_assignments_dict.items()):
    peak_list.sort(key=itemgetter(3))
    clean_resid_assignments_dict[resid] = peak_list

print("DEBUG: ITERATION 5 point 2 clean_resid_assignments_dict=", clean_resid_assignments_dict)

## WRITE THE CLEAN XEASY FILE
xeasy_out = open(args.ABSOLUTE_MATCHES_FILE+".NOESY.cleaned.iter5.xeasy", 'w')
atom_index = 0  # reinitialize atom index
for resid in list(clean_resid_assignments_dict.keys()):
    for peak in clean_resid_assignments_dict[resid]:
        atom_index += 1
        peak[0] = atom_index
        line = "\t".join([str(p) for p in peak])
        xeasy_out.write(line + "\n")
xeasy_out.close()


##
##  IF PROVIDED, LOAD THE USER-ANNOTATED TOCSY AND NOESY FILES FOR PROOF-READING
##
print("DEBUG: args.user_TOCSY_FILE=", args.user_TOCSY_FILE, "args.user_NOESY_FILE=", args.user_NOESY_FILE)
if args.user_TOCSY_FILE and args.user_NOESY_FILE:
    print("STARTING PROOF-READING NOESY CHEMICAL SHIFT ASSIGNMENTS")
    aa_equivalentCarbons_dict["TYR"] = [["CD1"], ["CE1"]]
    aa_equivalentCarbons_dict["PHE"] = [["CD1"], ["CE1"]]
    # THIS IS ITERATION 5, SO THERE ARE SOME EXTRA LINES THAT DON'T EXIST WHICH NEED TO BE ADDED
    for resid in clean_resid_assignments_dict:
        # # RENAME CD1-HD1 or CD2-HD2 --> CD-QD and CE1-HE1 or CE2-HE2 --> CE-QE if only one CD[12]-HD[12] or CE[12]-HE[12] has been assigned
        # # by the user in TYR of PHE residues
        # if clean_resid_assignments_dict[resid][0][6] in ["TYR", "PHE"]:
        #     CD_indices = [i for i,a in enumerate(clean_resid_assignments_dict[resid]) if a[3] in ["CD1", "CD2"]]
        #     if len(CD_indices) == 1:
        #         # clean_resid_assignments_dict[resid][CD_indices[0]][3] = "CD"
        #         clean_resid_assignments_dict[resid].append(clean_resid_assignments_dict[resid][CD_indices[0]])
        #     HD_indices = [i for i,a in enumerate(clean_resid_assignments_dict[resid]) if a[3] in ["HD1", "HD2"]]
        #     if len(HD_indices) == 1:
        #         clean_resid_assignments_dict[resid][HD_indices[0]][3] = "QD"
        #     
        #     CE_indices = [i for i,a in enumerate(clean_resid_assignments_dict[resid]) if a[3] in ["CE1", "CE2"]]
        #     if len(CE_indices) == 1:
        #         clean_resid_assignments_dict[resid][CE_indices[0]][3] = "CE"
        #     HE_indices = [i for i,a in enumerate(clean_resid_assignments_dict[resid]) if a[3] in ["HE1", "HE2"]]
        #     if len(HE_indices) == 1:
        #         clean_resid_assignments_dict[resid][HE_indices[0]][3] = "QE"
        
        if not resid in list(proofread_clean_resid_assignments_dict.keys()):
            proofread_clean_resid_assignments_dict[resid] = clean_resid_assignments_dict[resid]
            continue
        for new_assignment in clean_resid_assignments_dict[resid]:
            new_reson = new_assignment[1]
            old_reson_list = [a[1] for a in proofread_clean_resid_assignments_dict[resid]]
            if not new_reson in old_reson_list:     # if this resonance was not found (newly assigned), add it to the dictionary
                proofread_clean_resid_assignments_dict[resid].append(new_assignment)
    print("DEBUG: ITERATION 5 clean_resid_assignments_dict=", clean_resid_assignments_dict)
    
    # APPEND TO THE ASSINGED VAL AND LEU EQUIVALENT CARBONS DICTIONARIES THE TYR AND PHE RECORDS
    user_TYR_PHE_NOESY_resid_assignments_dict, user_TYR_PHE_TOCSY_resid_assignments_dict = create_equivalentCarbons_assignments_dict(user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict, ali.protein_alignment_list,
                                            ali.absolute_matches_alignment_list, aa_equivalentCarbons_dict, args)
    print("DEBUG: user_TYR_PHE_NOESY_resid_assignments_dict=", user_TYR_PHE_NOESY_resid_assignments_dict)
    print("DEBUG: user_TYR_PHE_TOCSY_resid_assignments_dict=", user_TYR_PHE_TOCSY_resid_assignments_dict)
    
    ##
    ##  COMPARE YOUR CARBON ASSIGNMENTS WITH THE USER ASSIGNMENTS
    ##
    ## Remember that the comments can be the following:
    # NOESY
    # NOESY (C-term hanging residue)
    # NOESY (flanked residue)
    # NOESY + TOCSY-HCNH
    # TOCSY (average)
    # TOCSY-HCNH
    # TOCSY-HCNH + NOESY
    # TOCSY-HCNH + TOCSY (unmatched)
    # TOCSY (unmatched)
    # TOCSY (unmatched) + TOCSY-HCNH
    # Recall that if we have 2 peaks in the same C-group but only one comes from NOESY ('TOCSY (unmatched) + TOCSY-HCNH', or
    # 'TOCSY-HCNH + TOCSY (unmatched)'), then we save only the NOESY Carbon resonance. 'TOCSY (unmatched)' -> only one resonance available
    # only in TOCSY, 'TOCSY (average)' -> 2 resonances available only from TOCSY and they were averaged
    
    for resid in list(proofread_clean_resid_assignments_dict.keys()):
        # FIRST TAKE CARE OF TYR AND PHE AND LEU AND VAL EQUIVALENT METHYL CARBONS
        resname = proofread_clean_resid_assignments_dict[resid][0][6]
        if resname in list(aa_equivalentCarbons_dict.keys()): # If this is a TYR or PHE or LEU or VAL
            if len(aa_equivalentCarbons_dict[resname]) > 1:   # proof-read TYR and PHE and TYR and PHE
                for Carbon_pair in aa_equivalentCarbons_dict[resname]:
                    if len([a for a in proofread_clean_resid_assignments_dict[resid] if a[3] in Carbon_pair]) > 0: # AND if the methyl carbons have been assigned
                        Proton_pair = ["H"+C[1:] for C in Carbon_pair]   # protons of the equivalent Carbons
                        # get only the methyl Carbon assignments (LEU --> CD1,CD2 or VAL --> CG1,CG2) made by the program, to compare them with the user made assignments
                        program_assignments_list = [assignment for assignment in proofread_clean_resid_assignments_dict[resid] if assignment[3] in Carbon_pair+Proton_pair]   # keep only the resonances assigned by the program to the equivalent methyl carbons
                        proof_read_equivalentMethyls(program_assignments_list, resid, user_TYR_PHE_NOESY_resid_assignments_dict, user_TYR_PHE_TOCSY_resid_assignments_dict, Carbon_pair, Proton_pair)   # it doesn't matter if 1 or 2 of the peaks have been proof-read in iteration 1 & 2
            else:   # proof-read LEU and VAL
                Carbon_pair = aa_equivalentCarbons_dict[resname][0]
                if len([a for a in proofread_clean_resid_assignments_dict[resid] if a[3] in Carbon_pair]) > 0: # AND if the methyl carbons have been assigned
                    Proton_pair = ["H"+C[1:] for C in Carbon_pair]   # protons of the equivalent Carbons
                    # get only the methyl Carbon assignments (LEU --> CD1,CD2 or VAL --> CG1,CG2) made by the program, to compare them with the user made assignments
                    program_assignments_list = [assignment for assignment in proofread_clean_resid_assignments_dict[resid] if assignment[3] in Carbon_pair+Proton_pair]   # keep only the resonances assigned by the program to the equivalent methyl carbons and protons
                    # proof_read_equivalentCarbons(program_assignments_list, resid, user_VAL_LEU_NOESY_resid_assignments_dict, user_VAL_LEU_TOCSY_resid_assignments_dict, Carbon_pair)   # OBSOLETE!!! it doesn't matter if 1 or 2 of the peaks have been proof-read in iteration 1 & 2
                    proof_read_equivalentMethyls(program_assignments_list, resid, user_VAL_LEU_NOESY_resid_assignments_dict, user_VAL_LEU_TOCSY_resid_assignments_dict, Carbon_pair, Proton_pair)   # it doesn't matter if 1 or 2 of the peaks have been proof-read in iteration 1 & 2
        # THEN PROCEED NORMALLY
        for assignment in proofread_clean_resid_assignments_dict[resid]:
            if '<' in assignment[7]:    # if this assignment has been proofread, then skip it
                continue
            print("DEBUG ITERATION5: proof reading new assignment=", assignment)
            spectrum_combo = re.sub(r' iter[1-9]', '', assignment[7])    # remove the "iter[1-9]" suffix to be able to distinguish the type of assignment (NEESY, NOESY+TOCSY, etc.)
            comment = None
            aa_type = assignment[6]
            nucleus = assignment[3] # the nucleus name
            if aa_type in list(aa_equivalentCarbons_dict.keys()) and nucleus in aa_equivalentCarbons_dict[aa_type]:    # if this peak has been assigned by the program to LEU or VAL CD1, CD2 or CG1,CG2 respectively, skip it
                continue
            
            # IF THIS CARBON IS A METHYL OR HAS ONLY ONE PROTON, CONSIDE ALSO THE PROTON RESONANCE FOR PROOFREADING
            if aa_type in list(aatype_carbon_nongeminalHname_mdict.keys()) and nucleus in list(aatype_carbon_nongeminalHname_mdict[aa_type].keys()):
                # get the proton assignment made by the program
                C_assignment = assignment # the program-made carbon assignment
                H_assignment = None # default value for the program-made proton assignment
                print("DEBUG: resid=", resid, "proofread_clean_resid_assignments_dict[resid]=", proofread_clean_resid_assignments_dict[resid])
                for a in proofread_clean_resid_assignments_dict[resid]:
                    if ((type(aatype_carbon_nongeminalHname_mdict[aa_type][nucleus]) == list and a[3] in aatype_carbon_nongeminalHname_mdict[aa_type][nucleus]) or
                        (type(aatype_carbon_nongeminalHname_mdict[aa_type][nucleus]) == str and a[3] == aatype_carbon_nongeminalHname_mdict[aa_type][nucleus])):
                        H_assignment = a
                        break
                print("DEBUG: C_assignment=", C_assignment, "H_assignment=", H_assignment)
                comment = ""
                if spectrum_combo in ['HCNH', 'NOESY (C-term hanging residue)', 'NOESY (flanked residue)', 'NOESY + TOCSY-HCNH', 'TOCSY-HCNH', 'TOCSY-HCNH + NOESY', 'TOCSY (unmatched) + TOCSY-HCNH', 'TOCSY-HCNH + TOCSY (unmatched)']:    # if this resonance comes only from NOESY
                    user_assignment = get_Carbon_resonance(resid, nucleus, 'HCNH', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict, get_assignment=True)
                    if user_assignment == False:    # if no user assignment was made for this Carbon, add no extra comment
                        print("WARNING: no user assignment made in NOESY for resid", resid, "and Carbon", nucleus,". Checking TOCSY...")
                        user_assignment = get_Carbon_resonance(resid, nucleus, 'TOCSY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict, get_assignment=True)
                        if user_assignment == False:
                            print("WARNING: no user assignment made also in TOCSY for resid", resid, "and Carbon", nucleus,". Skipping nucleus.")
                            assignment[7] += " <NOT ASSIGNED>"
                            #proofread_clean_resid_assignments_dict[resid].append(assignment)
                            continue
                    user_Creson = user_assignment[2]
                    user_Hreson = user_assignment[4]
                    print("DEBUG: user_Creson=", user_Creson, "user_Hreson=", user_Hreson)
                    if approx_equal(C_assignment[1], user_Creson, 0.4) and approx_equal(H_assignment[1], user_Hreson, 0.04): # if the user-made assignment agrees with ours, add to the comment "CORRECT"
                        print("DEBUG: CORRECT C_assignment found:", C_assignment)
                        print("DEBUG: CORRECT H_assignment found:", H_assignment)
                        comment = " <CORRECT>"
                    else:
                        print("DEBUG: WRONG C_assignment found:", C_assignment)
                        print("DEBUG: WRONG H_assignment found:", H_assignment)
                        comment = " <WRONG>"
                elif spectrum_combo in ['TOCSY (average)', 'TOCSY (unmatched)']:
                    user_assignment = get_Carbon_resonance(resid, nucleus, 'TOCSY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict, get_assignment=True)
                    if user_assignment == False:    # if no user assignment was made for this Carbon, add no extra comment
                        print("WARNING: no user assignment made in TOCSY for resid", resid, "and Carbon", nucleus,". Checking NOESY...")
                        user_assignment = get_Carbon_resonance(resid, nucleus, 'HCNH', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict, get_assignment=True)
                        if user_assignment == False:
                            print("WARNING: no user assignment made also in NOESY for resid", resid, "and Carbon", nucleus,". Skipping nucleus.")
                            assignment[7] += " <NOT ASSIGNED>"
                            #proofread_clean_resid_assignments_dict[resid].append(assignment)
                            continue
                    user_Creson = user_assignment[2]
                    user_Hreson = user_assignment[4]
                    if approx_equal(C_assignment[1], user_Creson, 0.4) and approx_equal(H_assignment[1], user_Hreson, 0.04): # if the user-made assignment agrees with ours, add to the comment "CORRECT"
                        print("DEBUG: CORRECT C_assignment found:", C_assignment)
                        print("DEBUG: CORRECT H_assignment found:", H_assignment)
                        comment = " <CORRECT>"
                    else:
                        print("DEBUG: WRONG C_assignment found:", C_assignment)
                        print("DEBUG: WRONG H_assignment found:", H_assignment)
                        comment = " <WRONG>"
                        
                if comment == "":
                    comment = " <NOT ASSIGNED>"
                C_assignment[7] += comment
                H_assignment[7] += comment
                comment = None  # a trick to add no more comments to "assignment" at the end of this block of code
            
            elif nucleus[0] == 'C':   # only if the assignment line belongs to a Carbon then check if the user has made the same assignment
                comment = ""
                if spectrum_combo in ['HCNH', 'NOESY (C-term hanging residue)', 'NOESY (flanked residue)', 'NOESY + TOCSY-HCNH', 'TOCSY-HCNH', 'TOCSY-HCNH + NOESY', 'TOCSY (unmatched) + TOCSY-HCNH', 'TOCSY-HCNH + TOCSY (unmatched)']:    # if this resonance comes only from NOESY
                    user_assignment = get_Carbon_resonance(resid, nucleus, 'HCNH', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict)
                    if user_assignment == False:    # if no user assignment was made for this Carbon, add no extra comment
                        print("WARNING: no user assignment made in NOESY for resid", resid, "and Carbon", nucleus,". Checking TOCSY...")
                        user_assignment = get_Carbon_resonance(resid, nucleus, 'TOCSY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict)
                        if user_assignment == False:
                            print("WARNING: no user assignment made also in TOCSY for resid", resid, "and Carbon", nucleus,". Skipping nucleus.")
                            assignment[7] += " <NOT ASSIGNED>"
                            # proofread_clean_resid_assignments_dict[resid].append(assignment)
                            continue
                    if approx_equal(assignment[1], user_assignment, 0.4): # if the user-made assignment agrees with ours, add to the comment "CORRECT"
                        print("DEBUG: CORRECT assignment found:", assignment)
                        comment = " <CORRECT>"
                    else:
                        print("DEBUG: WRONG assignment found:", assignment)
                        comment = " <WRONG>"
                elif spectrum_combo in ['TOCSY (average)', 'TOCSY (unmatched)']:
                    user_assignment = get_Carbon_resonance(resid, nucleus, 'TOCSY', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict)
                    if user_assignment == False:    # if no user assignment was made for this Carbon, add no extra comment
                        print("WARNING: no user assignment made in TOCSY for resid", resid, "and Carbon", nucleus,". Checking NOESY...")
                        user_assignment = get_Carbon_resonance(resid, nucleus, 'HCNH', user_NOESY_resid_assignments_dict, user_TOCSY_resid_assignments_dict)
                        if user_assignment == False:
                            print("WARNING: no user assignment made also in NOESY for resid", resid, "and Carbon", nucleus,". Skipping nucleus.")
                            assignment[7] += " <NOT ASSIGNED>"
                            # proofread_clean_resid_assignments_dict[resid].append(assignment)
                            continue
                    if approx_equal(assignment[1], user_assignment, 0.4): # if the user-made assignment agrees with ours, add to the comment "CORRECT"
                        print("DEBUG: CORRECT assignment found:", assignment)
                        comment = " <CORRECT>"
                    else:
                        print("DEBUG: WRONG assignment found:", assignment)
                        comment = " <WRONG>"
                
            if comment == None:      # this is a proton or a methyl Carbon or a Carbon with only one proton
                continue
            elif comment == "" and not '<' in assignment[7]:    # if no comment has been appended in previous iterations
                comment = " <NOT ASSIGNED>"
            assignment[7] += comment
    
    
    print("DEBUG: ITERATION 5 proofread_clean_resid_assignments_dict=", proofread_clean_resid_assignments_dict)
    
    ## WRITE THE PROOF-READ XEASY FILE
    xeasy_out = open(args.ABSOLUTE_MATCHES_FILE+".NOESY.proofread.iter5.xeasy", 'w')
    atom_index = 0  # reinitialize atom index
    for resid in list(proofread_clean_resid_assignments_dict.keys()):
        for peak in proofread_clean_resid_assignments_dict[resid]:
            atom_index += 1
            peak[0] = atom_index
            line = "\t".join([str(p) for p in peak])
            xeasy_out.write(line + "\n")
    xeasy_out.close()



## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FINAL ITERATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
# Simply renames all protons to the format that Rosetta understands. Also duplicates CD1-HD1 and CE1-HE1 peaks of TYR and PHE.
print("## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FINAL ITERATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##")
print("DEBUG: ENTERING FINAL ITERATION matched_NOESY_residue_assignments_dict=", matched_NOESY_residue_assignments_dict)
print("DEBUG: ENTERING FINAL ITERATION clean_resid_assignments_dict=", clean_resid_assignments_dict)
    
# Protons that cannot be grouped together
aatype_Hname_Qname_mdict = tree()
aatype_Hname_Qname_mdict['ALA']['HB'] = 'QB'
aatype_Hname_Qname_mdict['ILE']['HG2'] = 'QG2'
aatype_Hname_Qname_mdict['ILE']['HD1'] = 'QD1'
aatype_Hname_Qname_mdict['LEU']['HD1'] = 'QD1'
aatype_Hname_Qname_mdict['LEU']['HD2'] = 'QD2'
aatype_Hname_Qname_mdict['THR']['HG2'] = 'QG2'
aatype_Hname_Qname_mdict['VAL']['HG1'] = 'QG1'
aatype_Hname_Qname_mdict['VAL']['HG2'] = 'QG2'
aatype_Hname_Qname_mdict['MET']['HE'] = 'QE'


sorted_resids = list(clean_resid_assignments_dict.keys())
sorted_resids.sort()

for resid in sorted_resids:
    for assignment in clean_resid_assignments_dict[resid]:
        resname = assignment[6]
        atom_name = assignment[3]
        if resname in list(aatype_Hname_Qname_mdict.keys()) and atom_name in list(aatype_Hname_Qname_mdict[resname].keys()):
            print("DEBUG FINAL ITERATION: renaming", assignment[3], "to", aatype_Hname_Qname_mdict[resname][atom_name])
            assignment[3] = aatype_Hname_Qname_mdict[resname][atom_name]

print("DEBUG: FINAL ITERATION clean_resid_assignments_dict=", clean_resid_assignments_dict)
## WRITE THE CLEANED XEASY FILE
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


if args.user_TOCSY_FILE and args.user_NOESY_FILE:
    
    for resid in sorted_resids:
        for assignment in proofread_clean_resid_assignments_dict[resid]:
            resname = assignment[6]
            atom_name = assignment[3]
            if resname in list(aatype_Hname_Qname_mdict.keys()) and atom_name in list(aatype_Hname_Qname_mdict[resname].keys()):
                print("DEBUG FINAL ITERATION: renaming", assignment[3], "to", aatype_Hname_Qname_mdict[resname][atom_name])
                assignment[3] = aatype_Hname_Qname_mdict[resname][atom_name]
    
    print("DEBUG: FINAL ITERATION proofread_clean_resid_assignments_dict=", proofread_clean_resid_assignments_dict)
    ## WRITE THE PROOF-READ XEASY FILE
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



##
## Write in NOESY Sparky format the i-1 and i+1 peaks, too. Write also unassigned peaks with just the root index ID.
##

annotate_HCNH_file(args.NOESY_FILE, ali.absolute_AAIGmatches_alignment_list, ali.absolute_matches_alignment_list,
                    matched_NOESY_residue_assignments_dict, NOESY_residue_peak_intensity_mdict, out_fname="4DNOESY_assignedall.sparky")
