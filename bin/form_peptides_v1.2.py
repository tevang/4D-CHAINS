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
import sys
import traceback
from collections import OrderedDict
from operator import itemgetter

import numpy as np
from argparse import RawDescriptionHelpFormatter, ArgumentParser

from cluster import HierarchicalClustering
from lib.aa_prediction import get_aatypes_from_all_H_C_resonpairs, get_aatypes_from_H_C_resonpair_2Dhist, \
    get_aatypes_from_H_C_resonpair
from lib.bayes.statistics import Probability
from lib.global_func import bcolors, download_CS_histograms, ColorPrint
from lib.global_vars import aatype_carbon_nongeminalHname_mdict, aa_carbonBondedGeminalHydrogensDict_dict
from lib.hsqc_spectrum import HSQC_spectrum
from lib.probhist import ProbHist_Loader
from lib.tocsy_connectivities import TOCSY_Connectivities
from lib.trees.chain import Chain
from lib.trees.peptide import Peptide

CHAINS_BIN_DIR = os.path.dirname(os.path.realpath(__file__))


## Set global variables
HOME_DIR = os.path.dirname(os.path.realpath(__file__))


## Parse command line arguments
def cmdlineparse():
    parser = ArgumentParser(description="command line arguments", formatter_class=RawDescriptionHelpFormatter,
                            epilog="EXAMPLE: 4D_assignment_tree-based.py -hsqc TDhsqc.list -tocsy TDtocsy.list -noesy "
                                   "TDnoesy.list -tolH 0.01 -tolC 0.1 -rtolH 0.01 -rtolN 0.1 -mcutoff 1.0 -acutoff 1.0 "
                                   "-wH 0.5 -wC 1")
    parser.add_argument("-fasta", dest="FASTA_FILE", required=True,
                        help="teplate sequence file in fasta format", metavar="<template sequence file>")
    parser.add_argument("-hsqc", dest="HSQC_FILE", required=True,
                        help="2D N-H HSQC root spectrum", metavar="<2D N-H HSQC input file>")
    parser.add_argument("-tocsy", dest="TOCSY_fname", required=True,
                        help="4D TOCSY (HCTOCSYNH) file", metavar="<4D TOCSY input file>")
    parser.add_argument("-noesy", dest="NOESY_FILE", required=True,
                        help="4D NOESY (HCNOENH) file", metavar="<4D NOESY input file>")
    parser.add_argument("-tolH", dest="tolH", required=False, type=float, default=0.04,
                        help="tolerance for the proton resonance when matching NOESY peaks to TOCSY peaks. (default %(default)s)",
                        metavar="<proton tolerance>")
    parser.add_argument("-tolC", dest="tolC", required=False, type=float, default=0.4,
                        help="tolerance for the carbon resonance when matching NOESY peaks to TOCSY peaks. (default %(default)s)",
                        metavar="<carbon tolerance>")
    #parser.add_argument("-rtolH", dest="rtolH", required=False, type=float, default=0.02, help="tolerance for the proton resonance when matching peaks between the root spectum and TOCSY/NOESY", metavar="<proton tolerance>")
    #parser.add_argument("-rtolN", dest="rtolN", required=False, type=float, default=0.2, help="tolerance for the nitrogen resonance when matching peaks between the root spectum and TOCSY/NOESY", metavar="<carbon tolerance>")
    parser.add_argument("-update", dest="DOWNLOAD_CS_HISTOGRAMS", required=False, action='store_true',
                        help="download the latest chemical shift histgrams from BMRB")
    parser.add_argument("-wH", dest="H_weight", required=False, type=float, default=0.1,
                        help="weight in (0.0-1.0] to apply on aa type prediction from aliphatic H resonances. (default %(default)s)",
                        metavar="<H weight>")
    parser.add_argument("-wC", dest="C_weight", required=False, type=float, default=1.0,
                        help="weight in (0.0, 1.0] to apply on aa type prediction from aliphatic C resonances. (default %(default)s)",
                        metavar="<C weight>")
    #parser.add_argument("-examples", dest='EXAMPLES', action='store_true', required=False, help="Usage and Examples")
    parser.add_argument("-mcutoff", dest="RESONANCE_MATCH_CUTOFF", required=False, type=float, default=0.8,
                        help="number 0.0-1.0 saying how much of TOCSY resonances should match in the NOESY in order to consider it a possible match. (default %(default)s)",
                        metavar="<resonance match cutoff>")
    parser.add_argument("-zmcutoff", dest="ZSCORE_RESONANCE_MATCH_CUTOFF", required=False, type=float, default=0.0,
                        help="a real number specifying the lower Z-score for a resonance match to be retained for chain building. (default %(default)s)",
                        metavar="<Z-score resonance match cutoff>")
    parser.add_argument("-mratio", dest="MULTI_RATIO", required=False, type=float, default=1.0,
                        help="number 0.0-1.0 specifying the distance from the maximum multi score a possible "
                             "connectivity must have in order to be used.")
    parser.add_argument("-acutoff", dest="ASSIGNMENT_CUTOFF", required=False, type=float, default=None,
                        help="number 0.0-1.0 saying how much of TOCSY resonances should match with a particular aa type in order to be considered, e.g. \
                        if TOCSY resonances are 5, 0.8 will mean that at least 4 resonances must match. (default %(default)s)",
                        metavar="<aa type assignment cutoff>")
    parser.add_argument("-zacutoff", dest="ZSCORE_ASSIGNMENT_CUTOFF", required=False, type=float, default=-1.0,
                        help="a real number specifying the lower Z-score for an aa type prediction to be considered as valid. (default %(default)s)",
                        metavar="<Z-score aa type assignment cutoff>")
    #parser.add_argument("-moreaa", dest="ALLOW_SINGLE_CH_PAIRS", required=False, action='store_true', help="allow aa type predictions from single C-H resonance pairs. Because it increases memory consumption and lowers findelity, it should be used in the last rounds of assignment.")
    parser.add_argument("-maxlen", dest="MAX_PEPTIDE_LENGTH", required=False, type=int, default=8,
                        help="the maximum peptide length (high values require more memery; default: 8). (default %(default)s)",
                        metavar="<max peptide length>")
    parser.add_argument("-minlen", dest="MIN_PEPTIDE_LENGTH", required=False, type=int, default=3,
                        help="the minimum peptide length (low values increase noise; default: 2). (default %(default)s)",
                        metavar="<min peptide length>")
    parser.add_argument("-resoncut", dest="JUST_CARBON_MATCH_CUTOFF", required=False, type=int, default=3,
                        help="the minimum number of C-H resonance pairs for a TOCSY index group to start predicting aa types from both C-H or C \
                             only resonances. (default %(default)s)",
                        metavar="<Carbon prediction cutoff>")
    
    #parser.add_argument("-confile", dest="CONNECTIVITIES_FILE", required=False, help="connectivities file; if specified connectivity calculation from input files will be skipped", metavar="<connectivities file>")
    parser.add_argument("-poolconfile", dest="POOL_CONNECTIVITIES_FILE", required=False, default=None,
                        help="connectivities pool file; necessary if -confile specified in order to calculate correct probabilities", metavar="<connectivities pool file>")
    parser.add_argument("-allconfile", dest="COMPLETE_CONNECTIVITIES_FILE", required=False, default=None,
                        help="all connectivities file; necessary if -confile specified in order to calculate correct probabilities", metavar="<all connectivities file>")
    #parser.add_argument("-aafile", dest="AA_TYPES_FILE", required=False, help="amino acid assignment file; if specified amino acid assignment calculation from input files will be skipped", metavar="<amino acid assignment file>")
    parser.add_argument("-poolaafile", dest="POOL_AA_TYPES_FILE", required=False, default=None,
                        help="pool amino acid assignment file; if specified amino acid assignment calculation from input files will be skipped", metavar="<pool amino acid assignment file>")
    parser.add_argument("-allaafile", dest="COMPLETE_AA_TYPES_FILE", required=False, default=None,
                        help="all amino acid assignment file; necessary if -aafile specified in order to calculate correct probabilities", metavar="<all amino acid assignment file>")
    parser.add_argument("-chainfile", dest="NON_REDUNDANT_CHAINS_FILE", required=False, default=None,
                        help="non-redundant chains file; if specified, calculation of chains from connectivities will be omitted", metavar="<non-redundant chains file>")
    parser.add_argument("-zmin", dest="MIN_NUM_OF_PREDICTIONS", required=False, type=int, default=4,
                        help="minimum number of aa type predictions required to apply the Z-score cutoff. The fewer the predictions \
                        the more inaccurate is Z-score. (default %(default)s)", metavar="<minimum number of aa types prediction for Z-score filtering>")
    parser.add_argument("-transform", dest="TRANSFORM_TYPE", required=False, type=str, default='None',
                        help="type of mathematical transform to apply on the probabilities P[TAAIG(i)|aatype(i-1)]. Allowed values \
                             are: \"None\", \"log\", \"log10\", \"boxcox_pearson\", \"boxcox_mle\". (default %(default)s)",
                        metavar="<type of mathematical transform>")
    parser.add_argument("-resume", dest="RESUME", required=False, action='store_true', default=False,
                        help="resume previous run by loading all peptide sequences saved in tmp_peptide_folder/. By default tmp_peptide_folder/ \
                             will be cleaned and new peptide sequences will be writen inside it. (default %(default)s)")
    parser.add_argument("-probprod", dest="PROB_PRODUCT", required=False, action='store_true', default=False,
                        help="select the best C-H type assignment combination based on the product of probabilities of the individual \
                             C-H assignments. (default %(default)s)")
    parser.add_argument("-probmodel", dest="PROBABILITY_MODEL", required=False, type=int, default=2,
                        help="If '1' the probability of each peak will be given by [wH*1Dhist(H)+wC*1Dhist(C)]/(wH+wC). If '2' then \
                             by 1Dhist(H)*1Dhist(C). (default %(default)s)")
    parser.add_argument("-cgrpprob", dest="CONSENSUS_CGROUP_PROB_MODE", required=False, default=2, type=int,
                        help="""The way to calculate the total score of a set of chemical shift assignment (default %(default)s). Can be:
                        0: just multiply all the probabilities of the individual peaks
                        The following values control how to calculate the consensus probability of each C-group. The total score will be the
                        product of this consensus C-group probabilities.
                        1: average;
                        2: sqrt(prob1*prob2)    ; geometric mean
                        3: (prob2-prob1)/(log(prob2)-log(prob1))    ; logarithmic mean
                        4: (prob1*prob2)/((prob1+prob2)/2.) ;composite average 
                        5: prob1*prob2    ; product of probabilities
                        When the C-group contains only one peak, the respective probability will be prob1 in all cases. (default %(default)s)
                            """, metavar="<way to calculate the cs assignment score>")
    parser.add_argument("-2dhist", dest="USE_2D_HISTOGRAMS", required=False, action='store_true', default=False,
                        help="use 2D BMRB histograms for aa type prediction. (default %(default)s)")
    parser.add_argument("-log", dest="LOG_TRANSFORM", required=False, action='store_true', default=False,
                        help="convert aa type prediction probabilities to logarithmic scale and then calculate Z-scores, if the \
                             min(probability)/max(probability) > 1000. (default %(default)s)")
    parser.add_argument("-skipchains", dest="SKIP_CHAINS", required=False, action='store_true', default=False,
                        help="Generate connectivies and amino acid type predictions, but skip chain and peptide formation. (default %(default)s)")
    parser.add_argument("-delpred", dest="DELETE_AA_TYPE_PREDICTIONS", required=False, action='store_true', default=False,
                        help="delete aa type predictions with probabilities which are 1000, 10000 or 100000 times lower than the highest, if \
                             the highest is >10e-10, >10e-20, <=10e-20, respectively. (default %(default)s)")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s version 1.0')
    args=parser.parse_args()
    return args


##
##  CHECK INPUT FILES FORMAT FOR CORRECTNESS
##

def check_input_files_format(args):
    """
    Method to check input files' format for correctness.

    :return:
    """

    # Check NOESY file
    query_contents = []
    HAS_INTENSITIES = False
    with open(args.NOESY_FILE, 'r') as f:
        contents = f.readlines()
        for line in contents: # CLEANING THE SPECTRUM FROM IRRELEVANT LINES
            word_list = line.split()
            try:
                if not word_list[0] == '?-?-?-?':
                    continue
                float(word_list[1]); float(word_list[2]); float(word_list[3]); float(word_list[4]); # checking if the line contains numbers (resonances)
                if len(word_list) == 6:
                    float(word_list[5])
                    HAS_INTENSITIES = True
                query_contents.append(line)
            except (IndexError, ValueError):
                print(bcolors.FAIL + "DEBUG check_input_files_format: this NOESY line will be discarded: "+ line + bcolors.ENDC)
                continue
    # Now check if all the lines contain intensity
    if HAS_INTENSITIES:
        for word_list in query_contents:
            if len(word_list) < 6:
                print(bcolors.FAIL + "ERROR: Please check your input NOESY file. Not all lines have intensity!!!" + bcolors.ENDC)
                sys.exit(1)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##              CALCULATE BAYESIAN STATISTICS                        #
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

def split_oversized_Cgroups(possible_aatype_prob_C_H_resonpair_TAAIG_list_list):
    """
        FUNCTION to split Cgroups with more that 2 peaks into separate Cgroups. Currently works for maximum Cgroup size 3 only.
        This happens when two peaks have the same Carbon resonance but different proton resonance. So because C-grouping works
        with Carbon resonances only, all 3 peaks are included in the same group
    """
    
    aatypes_set = set([p[0] for p in possible_aatype_prob_C_H_resonpair_TAAIG_list_list])  # must be only one for CS assignment scripts
    corrected_possible_aatype_prob_C_H_resonpair_TAAIG_list_list = []
    for aatype in aatypes_set:
        aaspec_possible_aatype_prob_C_H_resonpair_TAAIG_list_list = [p for p in possible_aatype_prob_C_H_resonpair_TAAIG_list_list if p[0]==aatype] # aatype-specific possible_aatype_prob_C_H_resonpair_TAAIG_list_list
        clustID2size_dict = {}
        clustID2Creson_dict = {}
        for a in aaspec_possible_aatype_prob_C_H_resonpair_TAAIG_list_list:
            clustID = a[7]
            Creson = a[6]
            Hreson = a[5]
            if not clustID in list(clustID2Creson_dict.keys()):
                clustID2Creson_dict[clustID] = set([(Hreson, Creson)])
            else:
                clustID2Creson_dict[clustID].add((Hreson, Creson))
        
        for clustID in list(clustID2Creson_dict.keys()):
            clustID2size_dict[clustID] = len(clustID2Creson_dict[clustID])
        
        # print "DEBUG: clustID2size_dict=", clustID2size_dict
        maxClustID = max(clustID2size_dict.keys())
        for clustID in list(clustID2size_dict.keys()):
            if clustID2size_dict[clustID] == 3:
                Creson_list = [a[6] for a in aaspec_possible_aatype_prob_C_H_resonpair_TAAIG_list_list if a[7]==clustID]
                dist12 = abs(Creson_list[0] - Creson_list[1])
                dist13 = abs(Creson_list[0] - Creson_list[2])
                dist23 = abs(Creson_list[1] - Creson_list[2])
                if dist12 < dist13 and dist12 < dist23: # change the Cgroup of peak3
                    Coutlier = Creson_list[2]
                    for i in range(len(aaspec_possible_aatype_prob_C_H_resonpair_TAAIG_list_list)):
                        a = aaspec_possible_aatype_prob_C_H_resonpair_TAAIG_list_list[i]
                        if a[6] == Coutlier:
                            # print "DEBUG: changing Cgroup of ", a, "to ", maxClustID + 1
                            a[7] = maxClustID + 1
                            maxClustID += 1 # update the highest Cgroup number
                elif dist13 < dist12 and dist13 < dist23:   # change the Cgroup of peak2
                    Coutlier = Creson_list[1]
                    for i in range(len(aaspec_possible_aatype_prob_C_H_resonpair_TAAIG_list_list)):
                        a = aaspec_possible_aatype_prob_C_H_resonpair_TAAIG_list_list[i]
                        if a[6] == Coutlier:
                            # print "DEBUG: changing Cgroup of ", a, "to ", maxClustID + 1
                            a[7] = maxClustID + 1
                            maxClustID += 1 # update the highest Cgroup number
                elif dist23 < dist12 and dist23 < dist13:   # change the Cgroup of peak1
                    Coutlier = Creson_list[0]
                    for i in range(len(aaspec_possible_aatype_prob_C_H_resonpair_TAAIG_list_list)):
                        a = aaspec_possible_aatype_prob_C_H_resonpair_TAAIG_list_list[i]
                        if a[6] == Coutlier:
                            # print "DEBUG: changing Cgroup of ", a, "to ", maxClustID + 1
                            a[7] = maxClustID + 1
                            maxClustID += 1 # update the highest Cgroup number
                    
            elif clustID2size_dict[clustID] > 3:
                print("ERROR: found Cgroup with more than 3 peaks!!!")
                print("DEBUG: aaspec_possible_aatype_prob_C_H_resonpair_TAAIG_list_list=", aaspec_possible_aatype_prob_C_H_resonpair_TAAIG_list_list)
                sys.exit(1)
        corrected_possible_aatype_prob_C_H_resonpair_TAAIG_list_list.extend(aaspec_possible_aatype_prob_C_H_resonpair_TAAIG_list_list)
    
    return corrected_possible_aatype_prob_C_H_resonpair_TAAIG_list_list


def group_carbons(possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list, iteration=None):
    """
        This is a different function from the group_carbons() used in the CS assignment scripts, because there we know the aa type,
        but here we group the carbons for every possible aa type.
        ARGUMENTS:
        possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list:    list of lists of the form
                (aa type, Carbon name, Hydrogen name, weighted average probability, TOCSY_reson_index, H resonance, C resonance, carbon group)
        RETURNS:
        The same as the input but with the correct carbon groups (last element)
    """

    ## In Iterations 1 & 3 we do only the methyls therefore we must not group the carbons together
    if iteration in [1,3]:  # assigned a different C-group to each peak
        HCreson_set = set([(p[5], p[6]) for p in possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list])
        clustID = 0
        for HCreson in HCreson_set:
            clustID += 1
            for peak in possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list:
                if peak[5]==HCreson[0] and peak[6]==HCreson[1]:
                    peak[7] = clustID
        return possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list
    
    print("DEBUG group_carbons: possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list=", possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list)
    aa_types_set = set([x[0] for x in possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list])
    updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list = []    # same as possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list but with the cluster IDs for each aa type
    for aa_type in aa_types_set:
        singleAAtype_prob_C_H_resonpair_TOCSYindex_list_list = [x for x in possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list if x[0]==aa_type]
        if not aa_type in list(aa_carbonBondedGeminalHydrogensDict_dict.keys()):  # this aa does not have geminal protons, skip clustering
            clustID = 0
            for index in range(len(singleAAtype_prob_C_H_resonpair_TOCSYindex_list_list)):
                clustID += 1
                singleAAtype_prob_C_H_resonpair_TOCSYindex_list_list[index][7] = clustID
            # save the updated list entries with the correct cluster IDs
            updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list.extend(singleAAtype_prob_C_H_resonpair_TOCSYindex_list_list)
            continue    # continue with the next possible aa type
        
        new_AAIG_carbonGroupsTuple_dict = {}
        print("DEBUG: before singleAAtype_prob_C_H_resonpair_TOCSYindex_list_list=", singleAAtype_prob_C_H_resonpair_TOCSYindex_list_list)
        carbonCS_list = [l[6] for l in singleAAtype_prob_C_H_resonpair_TOCSYindex_list_list]
        print("DEBUG: carbonCS_list=", carbonCS_list)
        cl = HierarchicalClustering(carbonCS_list , lambda x,y: float(abs(x-y)))
        
        CONTINUE_CLUSTERING = True
        cutoff = 0.200001       # <== CHANGE ME (arbitrary starting value)
        while CONTINUE_CLUSTERING:  # continue clustering by lowering the cutoff until no more than 2 carbons resonances are within each cluster
            CONTINUE_CLUSTERING = False
            cluster_list = cl.getlevel(cutoff)    # <== CHANGE ME (the tolerance for carbon)
            print("DEBUG: cluster_list=", cluster_list)
            for cluster in cluster_list:
                try:
                    if len(set(cluster)) > 2:       # set because in carbonCS_list may be contained the same carbons with different atom type assignemnts
                        CONTINUE_CLUSTERING = True
                        cutoff -= 0.02
                        break
                except TypeError:   # in case only one nucleus prediction is available (and only one peak)
                    CONTINUE_CLUSTERING = False
            
        for index in range(len(singleAAtype_prob_C_H_resonpair_TOCSYindex_list_list)):
            carbon_resonance = singleAAtype_prob_C_H_resonpair_TOCSYindex_list_list[index][6]
            for clustID, clust in enumerate(cluster_list):
                print("DEBUG: clustID=", clustID, "clust=", clust, "carbon_resonance=", carbon_resonance)
                if len(cluster_list) == 1 and carbon_resonance == clust:    # to avoid "TypeError: argument of type 'float' is not iterable"
                    print("DEBUG: assigning cluster", clustID+1)
                    singleAAtype_prob_C_H_resonpair_TOCSYindex_list_list[index][7] = clustID + 1  # because enumeration starts from 0
                elif carbon_resonance in clust:
                    print("DEBUG: assigning cluster", clustID+1)
                    singleAAtype_prob_C_H_resonpair_TOCSYindex_list_list[index][7] = clustID + 1  # because enumeration starts from 0
        # save the updated list entries with the correct cluster IDs
        updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list.extend(singleAAtype_prob_C_H_resonpair_TOCSYindex_list_list)
        print("DEBUG group_carbons: appended carbon clusters with ", singleAAtype_prob_C_H_resonpair_TOCSYindex_list_list)
    
    
    ## IN the extreme scenario where a Cgroup has >2 peaks, keep the 2 peaks with closest Carbon resonance and the change the Cgroup of the others
    ## (Currently works only for 3 peaks in the same Cgroup).
    print("DEBUG: point 1 updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list=", updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list)
    updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list = \
        split_oversized_Cgroups(updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list)
    
    print("DEBUG: point 2 updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list=", updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list)
    if len(possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list[0]) == 9:    # if this is NOESY assignment (there is an 9th columns with the spectrum source)
        # For every Cname that is marked as 'TOCSY (unmatched)', search for another one Cname like this which is marked as 'TOCSY-HCNH',
        # and change its clusterID. This should happen because one Creson comes from NOESY and the other from TOCSY, therefore they could not
        # be grouped together with the current cutoff. This needs to be modified when I will search the second missing C-H peak in NOESY...
        matched_peaks = [x for x in possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list if x[8]=='TOCSY-HCNH']
        matched_Cnames = [x[1] for x in matched_peaks]
        matched_clustIDs = [x[7] for x in matched_peaks]
        for Cname, clustID in zip(matched_Cnames, matched_clustIDs):
            for index in range(len(updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list)):
                peak = updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list[index]
                if peak[8] in ['TOCSY-HCNH', 'TOCSY (unmatched)'] and peak[1]==Cname:
                    updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list[index][7] = clustID # change the cluster ID
    
    # Make sure that the same peaks (but with different predictions) belong to the same cluster ID
    print("DEBUG: point 3 updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list=", updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list)
    for aa_type in aa_types_set:
        CH_reson_tuple_set = set([(p[5],p[6]) for p in updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list if p[0]==aa_type])    # all unique C-H resonances
        for CH_reson_tuple in CH_reson_tuple_set:
            Creson = CH_reson_tuple[1]
            Hreson = CH_reson_tuple[0]
            # get all the assigned cluster IDs of this C-H peak (normally they should be one)
            clusterIDs_set = set([p[7] for p in updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list if p[5]==Hreson and p[6]==Creson and p[0]==aa_type])
            print("DEBUG: clusterIDs_set=", clusterIDs_set)
            common_clusterID = min(clusterIDs_set) # give them the same cluster ID (by default the smallest)
            for i in range(len(updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list)):
                if updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list[i][0]==aa_type and updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list[i][5]==Hreson and updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list[i][6]==Creson:
                    updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list[i][7] = common_clusterID    # reinitialize the cluster ID
    
    print("DEBUG: point 4 updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list=", updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list)
    # When 2 peaks are grouped they cannot be predicted to be a methyl (e.g. protein MS6282 L88)
    reupdated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list = []
    for aa_type in aa_types_set:
        clusterIDs_set = set([p[7] for p in updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list if p[0]==aa_type])   # all the cluster IDs
        for clustID in clusterIDs_set:
            CH_reson_tuple_set = set([(p[5],p[6]) for p in updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list if p[0]==aa_type and p[7]==clustID])    # get all peaks of this cluster
            print("DEBUG: CH_reson_tuple_set=", CH_reson_tuple_set)
            predictions_list = [p for p in updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list if p[0]==aa_type and p[7]==clustID]    # all C-H predictions of this cluster
            if len(CH_reson_tuple_set) == 2:    # if there are 2 peaks in this C-group, check if they have been predicted to be geminal
                for pred in predictions_list:
                    Cname = pred[1]
                    print("DEBUG: there are 2 peaks in this C-group, Cname=", Cname)
                    if not Cname in list(aatype_carbon_nongeminalHname_mdict[aa_type].keys()) or (aa_type in ["LEU", "VAL"] and Cname in ["CD1", "CD2", "CG1", "CG2"]):    # if this is a Carbon with geminal protons, keep the prediction, otherwise discard it
                        print("DEBUG: saving prediction:", pred)
                        reupdated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list.append(pred)
            else:   # if the C-group contains only one peak, then it can be also a methylene, so save the prediction
                print("DEBUG: there is only 1 peak in this C-group")
                for pred in predictions_list:
                    print("DEBUG: saving prediction:", pred)
                    reupdated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list.append(pred)
    
    print("DEBUG group_carbons: reupdated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list=", reupdated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list)
    return reupdated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list


##################################################### END OF FUNCTION DEFINITIONS #####################################################

if __name__ == "__main__":
    
    try:

        ############################################ SANITY CHECKS ############################################
        args = cmdlineparse()
        print("Input argument values:")
        for arg in vars(args):
            print(arg, "=", getattr(args, arg))
        if args.H_weight <= 0 or args.H_weight > 1:
            print("ERROR: -wH must be a number greater than 0 and lower or equal to 1!")
            sys.exit(1)
        if args.C_weight <= 0 or args.C_weight > 1:
            print("ERROR: -wC must be a number greater than 0 and lower or equal to 1!")
            sys.exit(1)
        if (args.ASSIGNMENT_CUTOFF != None and args.POOL_AA_TYPES_FILE != None) or (args.ASSIGNMENT_CUTOFF != None and args.COMPLETE_AA_TYPES_FILE != None):
            print("ERROR: Usage of -poolaafile or -allaafile is incompatible with -acutoff argument! Either remove -poolaafile and -allaafile to make new amino acid \
            type predictions from scratch, or remove -acutoff to use the amino acid prediction files you provided as input.")
            sys.exit(1)
        if args.ASSIGNMENT_CUTOFF == None and args.POOL_AA_TYPES_FILE == None and args.COMPLETE_AA_TYPES_FILE == None:
            print("ERROR: you must either provide an aa type assignment cutoff with -acutoff argument (recommented value: 1.0), or provide a pool amino acid assignment file with \
            argument -poolaafile and file with all the aa type assignment with -allaafile argument.")
            sys.exit(1)
        
        # Launch input file check    
        check_input_files_format(args)
        
        #if args.EXAMPLES:
        #    print "USAGE: python HCTOCSYNH_to_HCNOENH_tree-based.py -tocsy <HCTOCSYNH file> -noesy <HCNOENH file> [-stdH <stdev1>] [-stdC <stdev2>] [-update]"
        #    print "Example: ./4D_assignment_tree-based.py -hsqc TDhsqc.list -tocsy TDtocsy.list -noesy TDnoesy.list -stdH 0.03 -stdC 0.3 -rstdH 0.02 -rstdC 0.2 -cutoff 1.0"
        #    sys.exit(0)
        if args.DOWNLOAD_CS_HISTOGRAMS:
            download_CS_histograms()
        ############################################## END OF SANITY CHECKS ###################################


        # GRADUALLY REDUCE TOLERANCES TO FIND MORE NON-OVERLAPPING GROUPS IN HSQC SPECTRUM
        HSQC_spec = HSQC_spectrum(args.HSQC_FILE, args.FASTA_FILE)
        
        ##
        ## COPY AAIGs FROM HSQC SPECTRUM (HSQC-TOCSY matching & HSQC-NOESY matching)
        ##
        ## TOCSY_lines contains the AAIG of i, the H,C resonances of residue i-1 and the N,HN resonances of residue i
        TOCSY_lines, TOCSY_sidechain_resonances_list = \
            HSQC_spec.copy_AAIGnames_from_HSQC_spectrum(spectrum4D_fname=args.TOCSY_fname, spectrum4D_type="TOCSY")
        print("DEBUG: TOCSY_sidechain_resonances_list=", TOCSY_sidechain_resonances_list)
        NOESY_lines, _ = \
            HSQC_spec.copy_AAIGnames_from_HSQC_spectrum(spectrum4D_fname=args.NOESY_FILE, spectrum4D_type="HCNH",
                                                        TOCSY_sidechain_resonances_list=TOCSY_sidechain_resonances_list)
        ## Load TOCSY file
        TOCSY_contents = [l.split() for l in TOCSY_lines]
        TOCSY_contents.sort(key=itemgetter(0))  # sort TOCSY contents by the random AAIG
        
        #print "DEBUG: TOCSY_contents=",TOCSY_contents
        #print "DEBUG: NOESY_lines =",NOESY_lines
        # DEPRECATED
        #with open(args.NOESY_FILE, 'r') as f:
        #    NOESY_lines=f.readlines()
        
        ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
        ##                                          1st ROUND                                         ##
        ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
        
        ##
        ## GET POSSIBLE CONNECTIVITIES AND PROCESS THEM
        ##

        con = TOCSY_Connectivities(tolC=args.tolC, tolH=args.tolH)     # <== CHANGE ME
        con.load_TOCSY_spectrum(HSQC_spec.numlist_files['TOCSY'])   # load the labeled TOCSY
        con.load_HCNH_spectrum(HSQC_spec.numlist_files['HCNH'])  # load the labeled HCNH NOESY
        con.write_spectrum_info()   # write statistics about TOCSY and NOESY spectra
        con.__scale_intensities__()  # scale the intensities wrt the whole spectrum maximum

        def create_connectivity_dictionaries(con,
                                             args,
                                             clean_native_peaks=False,
                                             calc_intersection=True,
                                             uniquify=True,
                                             selected_i_iminus1_dict=None):

            con.match_TOCSY_to_NOESY_Peaks(tolH=args.tolH,
                                           tolC=args.tolC,
                                           calc_intersection=calc_intersection,
                                           uniquify=uniquify,
                                           selected_i_iminus1_dict=selected_i_iminus1_dict,
                                           clean_native_peaks=clean_native_peaks)
            i_iminus1_complete_dict = con.get_TOCSY_NOESY_connectivities_dict()  # all possible connectivities
            # i_iminus1_complete_dict contains all possible connectivities of every TOCSY AAIG; has keys the
                                # TOCSY AAIGs and values lists of quintuplets (tuples) \
                                # consisting of (NOESY AAIG i-1, occupancy, numOfResonances, intersection, multi score)

            #print "DEBUG:",numOfResonances, AAIG_signature, sorted_NOESYAAIG_occupancy_numOfResonances_tuple__list
            # Iterate over all NOESY aa indices and, if applicable, keep those that exceed the multi threshold
            i_iminus1_dict = {}
            for i_AAIG_name, connectivities in list(i_iminus1_complete_dict.items()):
                filtered_connectivities = []
                max_multi = np.max([float(q[4]) for q in connectivities])  # max multi for this AAIG i
                for quintuplet in connectivities:
                    multi = float(quintuplet[4])
                    if max_multi/multi <= args.MULTI_RATIO and multi > 10e-6:   # <== CHANGE ME
                        filtered_connectivities.append(quintuplet)
                i_iminus1_dict[i_AAIG_name] = filtered_connectivities
            
            return i_iminus1_dict, i_iminus1_complete_dict
        
        
        if not args.POOL_CONNECTIVITIES_FILE:
            print("Calculating Connectivities...")
            i_iminus1_dict, \
            i_iminus1_complete_dict = create_connectivity_dictionaries(con, args)
            # IMPORTANT: Load again the TOCSY and NOESY spectra to create a new TOCSY_Connectivities() object,
            # IMPORTANT: otherwise the program gets stuck and consumes all the memory!!!
            con = TOCSY_Connectivities(tolC=args.tolC, tolH=args.tolH)  # <== CHANGE ME
            con.load_TOCSY_spectrum(HSQC_spec.numlist_files['TOCSY'])
            con.load_HCNH_spectrum(HSQC_spec.numlist_files['HCNH'])
            con.__scale_intensities__()  # scale the intensities wrt the whole spectrum maximum
            i_iminus1_nonativepeaks_dict, \
            i_iminus1_complete_nonativepeaks_dict = create_connectivity_dictionaries(con, args, clean_native_peaks=True)
        else:
            print("Loading Connectivities from file", args.POOL_CONNECTIVITIES_FILE)
            i_iminus1_pool_dict = {} # same dictionary with the pool of possible connectivities remained after filtering using absolute consensus matches of a previous run
            i_iminus1_complete_dict = {} # the same dictionary but with all possible connectivities, including those below args.RESONANCE_MATCH_CUTOFF
        
            ## I THINK THIS IS NOT NECESSARY ##
            #with open(args.POOL_CONNECTIVITIES_FILE, 'r') as f:
            #    connectivities_file_contents = f.readlines()
            #for line in connectivities_file_contents[1:]:
            #    word_list = line.split()
            #    TOCSY_AAIG = word_list[0]
            #    i_iminus1_dict[TOCSY_AAIG] = []
            #    values_string = ''.join(word_list[1:])
            #    elements_string = re.sub('[\(\)\'\"]', '',  values_string ).split("),(")[0]
            #    elements_list = elements_string.split(",")
            #    for aa, occupancy, TOCSY_resonnum in zip(elements_list[0::3], elements_list[1::3], elements_list[2::3]):
            #        #print aa, occupancy, TOCSY_resonnum
            #        i_iminus1_dict[TOCSY_AAIG].append((aa, int(occupancy), int(TOCSY_resonnum)))
            ###################################
            
            # Now load the pool connectivities file to calculate correct probabilities
            i_iminus1_pool_dict = con.load_connectivities_from_file(args.POOL_CONNECTIVITIES_FILE)
            
            # Now load the complete connectivities file to calculate correct probabilities
            i_iminus1_complete_dict = con.load_connectivities_from_file(args.COMPLETE_CONNECTIVITIES_FILE)
            
            # FIND NEW CONNECTIVITIES FOR THOSE AAIGs THAT WERE NOT CORRECTED
            # FIRST QUICK RECONNAISSANCE ROUND: find missing connectivities
            ColorPrint("FIRST QUICK RECONNAISSANCE ROUND: find missing connectivities.", "BOLDGREEN")
            new_i_iminus1_dict, new_i_iminus1_complete_dict = create_connectivity_dictionaries(con,
                                                                                               args,
                                                                                               calc_intersection=False,
                                                                                               uniquify=False)    # new_i_iminus1_dict is useless here!
            selected_i_iminus1_dict = {}    # dict with the missing connectivities to be calculated later
            for TOCSYAAIG in list(new_i_iminus1_complete_dict.keys()):
                if TOCSYAAIG not in list(i_iminus1_complete_dict.keys()):  # if there was not connectivity for this TAAIG add it to the dictionaries
                    selected_i_iminus1_dict[TOCSYAAIG] = new_i_iminus1_complete_dict[TOCSYAAIG]
                elif len(new_i_iminus1_complete_dict[TOCSYAAIG]) > len(i_iminus1_complete_dict[TOCSYAAIG]): # if there were found extra connectivities using the new tolerances
                    if not TOCSYAAIG in list(i_iminus1_pool_dict.keys()):  # probably all NOESY AAIGs which had connectivities have been mapped
                                                                        # to the sequence
                        # TODO: Maybe it would make more sense to add only the new connenctivites that did not exist in
                        # TODO: i_iminus1_complete_dict to selected_i_iminus1_dict.
                        selected_i_iminus1_dict[TOCSYAAIG] = new_i_iminus1_complete_dict[TOCSYAAIG]  # update the connectivities
                    elif len(i_iminus1_pool_dict[TOCSYAAIG]) == 1 and len(i_iminus1_complete_dict[TOCSYAAIG]) > 1:
                        are_all_unique = True   # do all Tindices in the complete dict have a single connectivity in the pool
                        for tmp_quintuplet in i_iminus1_complete_dict[TOCSYAAIG]:
                            tmp_TOCSYAAIG = tmp_quintuplet[0]
                            if tmp_TOCSYAAIG in list(i_iminus1_pool_dict.keys()) and len(i_iminus1_pool_dict[tmp_TOCSYAAIG]) > 1:
                                are_all_unique = False
                                break
                        if are_all_unique == False: # this means that i_iminus1_pool_dict[AAIG_signature] has been modified and not cleaned from Tindices used elsewhere
                            continue    # leave the connectivity of this AAIG_signature unchanged
                        selected_i_iminus1_dict[TOCSYAAIG] = new_i_iminus1_complete_dict[TOCSYAAIG]
                    elif len(i_iminus1_pool_dict[TOCSYAAIG]) == 1 and len(i_iminus1_complete_dict[TOCSYAAIG]) == 1:
                        selected_i_iminus1_dict[TOCSYAAIG] = new_i_iminus1_complete_dict[TOCSYAAIG]   # update the connectivities

            ColorPrint("SECOND ROUND: calculate intersections of missing connectivities and add them.", "BOLDGREEN")
            con.clean_connectivities()  # VERY IMPORTANT!!!
            new_i_iminus1_dict, \
            new_i_iminus1_complete_dict = create_connectivity_dictionaries(con,
                                                                           args,
                                                                           calc_intersection=True,
                                                                           uniquify=True,
                                                                           selected_i_iminus1_dict=selected_i_iminus1_dict)    # new_i_iminus1_dict is useless here!
            for TOCSYAAIG in list(new_i_iminus1_complete_dict.keys()):
                if TOCSYAAIG not in list(i_iminus1_complete_dict.keys()):  # if there was not connectivity for this TAAIG add it to the dictionaries
                    i_iminus1_complete_dict[TOCSYAAIG] = new_i_iminus1_complete_dict[TOCSYAAIG]   # update the connectivities
                    i_iminus1_pool_dict[TOCSYAAIG] = new_i_iminus1_complete_dict[TOCSYAAIG]
                elif len(new_i_iminus1_complete_dict[TOCSYAAIG]) > len(i_iminus1_complete_dict[TOCSYAAIG]): # if there were found extra connectivities using the new tolerances
                    if not TOCSYAAIG in list(i_iminus1_pool_dict.keys()):  # probably all NOESY AAIGs which had connectivities have been mapped
                        # to the sequence
                        # TODO: Maybe it would make more sense to add only the new connenctivites that did not exist in
                        # TODO: i_iminus1_complete_dict to i_iminus1_pool_dict and i_iminus1_complete_dict.
                        i_iminus1_complete_dict[TOCSYAAIG] = new_i_iminus1_complete_dict[TOCSYAAIG]  # update the connectivities
                        i_iminus1_pool_dict[TOCSYAAIG] = new_i_iminus1_complete_dict[TOCSYAAIG]
                    elif len(i_iminus1_pool_dict[TOCSYAAIG]) == 1 and len(i_iminus1_complete_dict[TOCSYAAIG]) > 1:
                        are_all_unique = True   # do all Tindices in the complete dict have a single connectivity in the pool
                        for tmp_quintuplet in i_iminus1_complete_dict[TOCSYAAIG]:
                            tmp_TOCSYAAIG = tmp_quintuplet[0]
                            if tmp_TOCSYAAIG in list(i_iminus1_pool_dict.keys()) and len(i_iminus1_pool_dict[tmp_TOCSYAAIG]) > 1:
                                are_all_unique = False
                                break
                        if are_all_unique == False: # this means that i_iminus1_pool_dict[AAIG_signature] has been modified and not cleaned from Tindices used elsewhere
                            continue    # leave the connectivity of this AAIG_signature unchanged
                        i_iminus1_complete_dict[TOCSYAAIG] = new_i_iminus1_complete_dict[TOCSYAAIG]   # update the connectivities
                        i_iminus1_pool_dict[TOCSYAAIG] = new_i_iminus1_complete_dict[TOCSYAAIG]
                    elif len(i_iminus1_pool_dict[TOCSYAAIG]) == 1 and len(i_iminus1_complete_dict[TOCSYAAIG]) == 1:
                        i_iminus1_complete_dict[TOCSYAAIG] = new_i_iminus1_complete_dict[TOCSYAAIG]   # update the connectivities
                        i_iminus1_pool_dict[TOCSYAAIG] = new_i_iminus1_complete_dict[TOCSYAAIG]

            # CALCULATING CONNECTIVITIES FROM THE POOL OF CONNECTIVITIES USING THE args.MULTI_RATIO
            print("Finding Connectivities within the -mratio from the Pool of Connectivities...")
            # connectivities_mdict is a multidimensional dictionary with structure:
            # TOCSY AAIG => NOESY AAIG => [number of matched TOCSY w2,w3 resonances in NOESY for that NOESY AAIG,
            # total number of TOCSY w2,w3 resonances for that TOCSY AAIG]
            i_iminus1_dict = {} # a dictionary containing all possible connectivities of every TOCSY AAIG above 'multi' threshold;
                                # has keys the TOCSY aa indices and values lists of quintuplets (tuples) \
                                # consisting of (NOESYAAIG, occupancy, numOfResonances)
            for TOCSYAAIG in list(i_iminus1_pool_dict.keys()):
                # Iterate over all NOESY aa indices and, if applicable, keep those that match in all resonance pairs w2,w3
                max_multi = np.max([float(q[4]) for q in i_iminus1_pool_dict[TOCSYAAIG]])
                for quintuplet in i_iminus1_pool_dict[TOCSYAAIG]:
                    multi = float(quintuplet[4])
                    if max_multi/multi <= args.MULTI_RATIO and multi > 10e-6:     # adjustable threshold
                        try:
                            i_iminus1_dict[TOCSYAAIG].append(quintuplet)
                        except KeyError:
                            i_iminus1_dict[TOCSYAAIG] = [(quintuplet)]

        # NOT SURE IF THIS IS ANYWHERE USED!
        ## Convert occupancies to probabilities
        i_iminus1_normProbabilities_dict = {}   # dict of the form TAAIG(i)->[(TAAIG(i-1), probability), (...), ...];
                                                # eventually will contain only the matches above the cutoff
        AAIG_weightSum_dict = {}
        for TOCSY_AAIG in list(i_iminus1_complete_dict.keys()):
            weight_sum = 0
            for quintuplet in i_iminus1_complete_dict[TOCSY_AAIG]:  # use all possible matches to calculate the weight_sum
                weight = quintuplet[4]   # multi score
                weight_sum += weight
            AAIG_weightSum_dict[TOCSY_AAIG] = weight_sum
        
        for TOCSY_AAIG in list(i_iminus1_dict.keys()):     # use only the matches above the cutoff to save theirs probabilities
            #print "DEBUG: TOCSY_AAIG=", TOCSY_AAIG, "AAIG_weightSum_dict[TOCSY_AAIG] =",AAIG_weightSum_dict[TOCSY_AAIG]
            for quintuplet in i_iminus1_dict[TOCSY_AAIG]:
                iminus1_AAIG_name = quintuplet[0]
                weigh = quintuplet[4]   # the multi score
                duplet = (iminus1_AAIG_name, weight/AAIG_weightSum_dict[TOCSY_AAIG])           # recover the probability by dividing by the weight_sum
                try:
                    i_iminus1_normProbabilities_dict[TOCSY_AAIG].append(duplet)
                except KeyError:
                    i_iminus1_normProbabilities_dict[TOCSY_AAIG] = [duplet]
        
         
        #for k,v in i_iminus1_normProbabilities_dict.items():
        #    print k,v
        #    prob_sum = 0
        #    for duplet in v:
        #        prob_sum += duplet[1]
        #    print prob_sum
        #
        #sys.exit(1)
        
        ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
        ##              FORM CHAINS FROM CONNECTIVITIES               ##
        ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
        
        if not args.NON_REDUNDANT_CHAINS_FILE:

            with open("connectivities_mratio_"+str(args.MULTI_RATIO), 'w') as f:
                f.write("i\tpossible i-1\n")
                for k,v in list(i_iminus1_dict.items()):
                    print(k,v)
                    f.write(k+"\t"+', '.join(str(e) for e in v)+"\n")    
            
            with open("connectivities_all", 'w') as f:
                f.write("i\tpossible i-1\n")
                for k,v in list(i_iminus1_complete_dict.items()):
                    #print k,v
                    f.write(k+"\t"+', '.join(str(e) for e in v)+"\n")
            
            try:    # write connectivities without native peaks only if the respective dictionaries exist
                with open("connectivities_nonativepeaks_mratio_"+str(args.MULTI_RATIO), 'w') as f:
                    f.write("i\tpossible i-1\n")
                    for k,v in list(i_iminus1_nonativepeaks_dict.items()):
                        print(k,v)
                        f.write(k+"\t"+', '.join(str(e) for e in v)+"\n")
                with open("connectivities_all_nonativepeaks", 'w') as f:
                    f.write("i\tpossible i-1\n")
                    for k,v in list(i_iminus1_complete_nonativepeaks_dict.items()):
                        #print k,v
                        f.write(k+"\t"+', '.join(str(e) for e in v)+"\n")
            except NameError:
                pass
            
            if not args.SKIP_CHAINS:
                ColorPrint("Forming chains.", "BOLDGREEN")
                chain = Chain(i_iminus1_dict=i_iminus1_dict,
                              i_iminus1_complete_dict=i_iminus1_complete_dict)
                chain.form_chains(MAX_PEPTIDE_LENGTH=args.MAX_PEPTIDE_LENGTH,
                                  MIN_PEPTIDE_LENGTH=args.MIN_PEPTIDE_LENGTH)
                ColorPrint("Saving all the chains to chains.list file ...", "BOLDBLUE")
                chain.write_chains_file(chain.all_chainScore_set, "chains.list")
                ColorPrint("Removing redundant chains ...", "BOLDBLUE")
                all_chainScore_list = chain.remove_redundant_chains()
                chain.write_chains_file(all_chainScore_list, "chains.non-redundant.list")

                ##
                ##  Remove redudancy from the chain list
                ### Bash script to help in validating the redundancy removal
                #while read line;
                #    do chain=$(perl -pi -e "s/.* chain = \((.*)\)/\1/" <<< $line);
                #    score=$(awk '{print $4}' <<< $line);
                #    echo "---> redundant chain = $chain";
                #    grep "Score = $score chain = ($chain" chains.non-redundant.list | grep "$chain";
                #done < chains.redundant.list > redundant_lines.txt

        else:   # otherwise load the specified non-redudant chains file
            all_chainScore_list = Chain.load_chains_file(args.NON_REDUNDANT_CHAINS_FILE)
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        #   AUTOMATIC ASSIGNMENT OF AMINO ACID TYPES TO TOCSY RESONANCES    #
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        
        if not args.POOL_AA_TYPES_FILE and not args.COMPLETE_AA_TYPES_FILE:
            # Load 1D and 2D probability histograms and export hist object for probability calculations
            histload = ProbHist_Loader()
            histload.load_1Dhistograms()

            if args.USE_2D_HISTOGRAMS == True:
                ## IF ASKED IN THE COMMAND LINE, LOAD 2D BMRB HISTOGRAMS, TOO
                histload.load_2Dhistograms()
            
            
            # sys.exit(0) # TEMPTORARILY DEACTIVATE AA TYPE PREDICTION.
            iAAIG_iminus1aaTypesProbTupleList_dict = OrderedDict() # ordereddict with keys the AAIG of residue i and values lists of tuples of the form (residue i-1 matching aa type, average probability)
            Num_of_TOCSY_resonances = 0     # save here the number of TOCSY resonances for a particular AAIG
            previous_TOCSY_AAIG = None
            possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list = []
            carbon_groups_list = []  # list of tuples of the form (AAIG, H resonance, C resonance, N resonance, HN resonance, carbon group)
            AAIG_carbonGroupsTuple_dict = {} # dict of the form: AAIG -> list of the form (AAIG, H resonance, C resonance, N resonance, HN resonance, carbon group)
            counter = 0
            previous_AAIG = ""
            sorted_TOCSY_contents = sorted(TOCSY_contents, key=itemgetter(0))
            # print("sorted_TOCSY_contents=", sorted_TOCSY_contents)
            #sys.exit(1)
            #sorted_TOCSY_contents= [
            #    ['A6', '4.148', '61.783', '131.690', '8.959'], ['A6', '1.829', '33.205', '131.701', '8.960'], ['A6', '1.038', '58.474', '131.692', '8.961'],
            #    ['A6', '0.756', '57.876', '131.676', '8.961']
            #    ]
            for TOCSY_words_list in sorted_TOCSY_contents:
                try:
                    # ignore 1st column
                    TOCSY_AAIG=TOCSY_words_list[0]    # residue i
                    print("DEBUG: TOCSY_AAIG=",TOCSY_AAIG,"previous_TOCSY_AAIG=",previous_TOCSY_AAIG)
                    if previous_TOCSY_AAIG != None and TOCSY_AAIG != previous_TOCSY_AAIG: # if we look at a different AAIG in TOCSY, print the matches and occupancies
                                                                                                   # of the previous TOCSY AAIG and clear matchingNOESYAAIG_occupancy_dict
                        
                        print("Assigning possible aa types to the aa upstream of AAIG",previous_TOCSY_AAIG)
                        #print "DEBUG: switched from AAIG", previous_TOCSY_AAIG, "to", TOCSY_AAIG
                        #print "DEBUG: possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list=",possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list
                        #print "DEBUG: Num_of_TOCSY_resonances=",Num_of_TOCSY_resonances
                        # GROUP THE CARBON RESONANCE OF EACH TAAIG
                        new_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list = group_carbons(possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list)
                        iAAIG_iminus1aaTypesProbTupleList_dict[previous_TOCSY_AAIG] = get_aatypes_from_all_H_C_resonpairs(new_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list, Num_of_TOCSY_resonances, args)
                        #print "DEBUG: previous_TOCSY_AAIG=",previous_TOCSY_AAIG,"iAAIG_iminus1aaTypesProbTupleList_dict[previous_TOCSY_AAIG]=",iAAIG_iminus1aaTypesProbTupleList_dict[previous_TOCSY_AAIG]
                        possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list = [] # list of tuples of the form (aa type, probability, H_resonance, C_resonance, TOCSY_reson_index) containing all matching sets of H,C resonances
                        AAIG_carbonGroupsTuple_dict[TOCSY_AAIG] = carbon_groups_list
                        #print "DEBUG: AAIG_carbonGroupsTuple_dict=", AAIG_carbonGroupsTuple_dict
                        carbon_groups_list = []  # list of tuples of the form (TOCSY_AAIG, H resonance, C resonance, N resonance, HN, resonance, carbon group)
                        
                        Num_of_TOCSY_resonances = 0
                        #counter += 1
                        #if counter == 10:
                        #    break
                    
                    TOCSY_H_resonance=float(TOCSY_words_list[1]) # aliphatic H resonance of residue i-1 
                    TOCSY_C_resonance=float(TOCSY_words_list[2]) # aliphatic C (Ca,Cb,Cc,Cg,Ce,etc.) resonance of residue i-1; this C is covalently bonded to the above H
                    TOCSY_N_resonance=float(TOCSY_words_list[3])
                    TOCSY_HN_resonance=float(TOCSY_words_list[4])
                    Num_of_TOCSY_resonances += 1
                    #print "DEBUG: TOCSY_C_resonance=",TOCSY_C_resonance,"TOCSY_H_resonance=",TOCSY_H_resonance, "Num_of_TOCSY_resonances", Num_of_TOCSY_resonances
                    # valid_matches_list: list of lists of the form [aa type, Carbon name, Hydrogen name, weighted average probability, TOCSY_reson_index, H resonance, C resonance, carbon group]
                    if args.USE_2D_HISTOGRAMS == True:
                        valid_matches_list = get_aatypes_from_H_C_resonpair_2Dhist(TOCSY_H_resonance,
                                                                                   TOCSY_C_resonance,
                                                                                   Num_of_TOCSY_resonances,
                                                                                   histload)
                    else:
                        valid_matches_list = get_aatypes_from_H_C_resonpair(TOCSY_H_resonance,
                                                                            TOCSY_C_resonance,
                                                                            Num_of_TOCSY_resonances,
                                                                            histload,
                                                                            args)
                    #if TOCSY_AAIG == "Y130":
                    #print "DEBUG: valid_matches_list=",valid_matches_list
                    possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list.extend(valid_matches_list)
                    carbon_groups_list.append([TOCSY_AAIG, TOCSY_H_resonance, TOCSY_C_resonance, TOCSY_N_resonance, TOCSY_HN_resonance, None])
                
                    previous_TOCSY_AAIG = TOCSY_AAIG
                except (ValueError, IndexError):
                    print("WARNING: the 3rd and 4th elements of the following TOCSY file line are not numbers:")
                    print("TOCSY file line:", TOCSY_words_list)
                    continue
            # THIS IS FOR THE LAST TOCSY AAIG
            #if not counter == 10:
            print("Assigning possible aa types to the last aa, which is upstream of AAIG",previous_TOCSY_AAIG)
            #print "DEBUG: possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list=",possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list
            # GROUP THE CARBON RESONANCE OF EACH TAAIG
            new_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list = group_carbons(possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list)
            iAAIG_iminus1aaTypesProbTupleList_dict[previous_TOCSY_AAIG] = get_aatypes_from_all_H_C_resonpairs(new_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list, Num_of_TOCSY_resonances, args)
            #print "DEBUG: previous_TOCSY_AAIG=",previous_TOCSY_AAIG,"iAAIG_iminus1aaTypesProbTupleList_dict[previous_TOCSY_AAIG]=",iAAIG_iminus1aaTypesProbTupleList_dict[previous_TOCSY_AAIG]

            ## OBSOLETE!
            # ##### NEW WAY: CALCULATE Z-SCORE FOR EACH TOCSY AAIG INDIVIDUAL #####
            # ## Calculate Z-scores from ALL aa type prediction probabilites and keep only those predictions above the cutoff
            # iAAIG_iminus1aaTypesZscoreTupleList_dict = OrderedDict() # ordereddict with keys the AAIG of residue i and values lists of tuples of the form (residue i-1 matching aa type, Z-score)
            # iAAIG_iminus1aaTypesCutoffProbTupleList_dict = OrderedDict() # contains only the tuples (residue i-1 matching aa type, probability) with Z-score above the cutoff
            # for iAAIG, iminus1aaTypesProbTuple_list in list(iAAIG_iminus1aaTypesProbTupleList_dict.items()):
            #     # DECIDE THE PROBABILITY THRESHOLD TO DISCARD AA TYPE PREDICTIONS BELOW IT
            #     try:
            #         max_prob = np.max([duplet[1] for duplet in iminus1aaTypesProbTuple_list])   # highest aa type probability
            #     except ValueError:  # if no aa type predictions exist for this TAAIG group
            #         max_prob = 0
            #     if max_prob > 10e-10:   # Dedice the probability threshold
            #         prob_threshold = 1000
            #     elif max_prob > 10e-20:
            #         prob_threshold = 10000
            #     elif max_prob <= 10e-20:
            #         prob_threshold = 100000
            #     print("DEBUG: max_prob", max_prob, "prob_threshold", prob_threshold, "args.LOG_TRANSFORM=", args.LOG_TRANSFORM)
            #     aatype_list, prob_list = [], []
            #     for duplet in iminus1aaTypesProbTuple_list:
            #         if args.DELETE_AA_TYPE_PREDICTIONS == True:     # delete low probability aa type predictions
            #             try:
            #                 if max_prob/float(duplet[1]) < prob_threshold:
            #                     aatype_list.append(duplet[0])
            #                     prob_list.append(duplet[1])  # unsorted list of probabilities
            #             except ZeroDivisionError:
            #                 print("DEBUG: ZeroDivisionError ", max_prob, "/", float(duplet[1]))
            #         elif args.DELETE_AA_TYPE_PREDICTIONS == False:
            #             aatype_list.append(duplet[0])
            #             prob_list.append(duplet[1])  # unsorted list of probabilities
            #     print("DEBUG: prob_list=", prob_list)
            #     if len(prob_list) == 1:
            #         zscore_array = np.array([10.00000])     # default Z-score in the case of just one prediction
            #     elif len(prob_list) > 2:
            #         if args.LOG_TRANSFORM == True:
            #             try:
            #                 ratio = float(np.max(prob_list))/np.min(prob_list)    # if the min probability is at least 3 orders of magnitude smaller, convert them to logarithmic scale
            #             except ZeroDivisionError:
            #                 print("DEBUG: ZeroDivisionError ", np.max(prob_list), "/", np.min(prob_list))
            #                 ratio = -1.0
            #             if ratio > 1000:
            #                 zscore_array = zscore(np.log(prob_list))  # convert probabilities the logarithms and then to Z-scores
            #             else:
            #                 zscore_array = zscore(prob_list)  # convert probabilities to Z-scores
            #         elif args.LOG_TRANSFORM == False:
            #             zscore_array = zscore(prob_list)  # convert probabilities to Z-scores
            #     elif len(prob_list) == 2:   # if only 2 predictions exist, set the 2nd to 0.000 to allow you to use high -zcutoff values (<= 0)
            #         zscore_array = np.array([1.00000, 0.00000])
            #     else:   # if no aa type prediction was made create an empty array
            #         zscore_array = np.array([])
            #     print("DEBUG2: saving zscore_array = ", zscore_array)
            #     iminus1aaTypesZscoreTuple_list = []
            #     iminus1aaTypesCutoffProbTuple_list = []
            #     for aatype, Zscore, prob in zip(aatype_list, zscore_array, prob_list):
            #         ## ONLY IF THE Z-SCORE OF THE AA TYPE PREDICTION IS GREATER THAN THE CUTOFF AND THE AA TYPE PREDICTIONS ARE AT LEAST args.MIN_NUM_OF_PREDICTIONS, THEN INCLUDE IT IN THE LIST
            #         if Zscore > args.ZSCORE_ASSIGNMENT_CUTOFF and len(zscore_array) >= args.MIN_NUM_OF_PREDICTIONS:
            #             aatypeProb_tuple = (aatype, prob)
            #             aatypeZscore_tuple = (aatype, Zscore)  # replace the probability with the respective Z-score
            #             iminus1aaTypesZscoreTuple_list.append(aatypeZscore_tuple)
            #             iminus1aaTypesCutoffProbTuple_list.append(aatypeProb_tuple) # but also save the aatype probabilities that satisfy the Z-score cutoff
            #         ## IF THE AA TYPE PREDICTIONS ARE LESS THAN args.MIN_NUM_OF_PREDICTIONS, THEN INCLUDED ALL OF THEM IN THE LIST
            #         elif len(zscore_array) < args.MIN_NUM_OF_PREDICTIONS:
            #             aatypeProb_tuple = (aatype, prob)
            #             aatypeZscore_tuple = (aatype, Zscore)  # replace the probability with the respective Z-score
            #             iminus1aaTypesZscoreTuple_list.append(aatypeZscore_tuple)
            #             iminus1aaTypesCutoffProbTuple_list.append(aatypeProb_tuple) # but also save the aatype probabilities that satisfy the Z-score cutoff
            #     iAAIG_iminus1aaTypesZscoreTupleList_dict[iAAIG] = iminus1aaTypesZscoreTuple_list
            #     iAAIG_iminus1aaTypesCutoffProbTupleList_dict[iAAIG] = iminus1aaTypesCutoffProbTuple_list

            # Use ALL aa-type predictions (we don't have a pool yet)
            iAAIG_iminus1aaTypesCutoffProbTupleList_dict, \
            iAAIG_iminus1aaTypesZscoreTupleList_dict = \
                Probability.get_pruned_prob_dicts_by_zscore(iAAIG_iminus1aaTypesProbTupleList_dict,
                                        ZSCORE_ASSIGNMENT_CUTOFF=args.ZSCORE_ASSIGNMENT_CUTOFF,
                                        DELETE_AA_TYPE_PREDICTIONS=args.DELETE_AA_TYPE_PREDICTIONS,
                                        LOG_TRANSFORM=args.LOG_TRANSFORM,
                                        MIN_NUM_OF_PREDICTIONS=args.MIN_NUM_OF_PREDICTIONS)

            ## OBSOLETE!
            # # Save ALL amino acid type predictions with the respective probabilities
            # #with open("amino_acid_type_prediction_probabilities.mcutoff"+str(args.RESONANCE_MATCH_CUTOFF)+".acutoff"+str(args.ASSIGNMENT_CUTOFF)+".zacutoff"+str(args.ZSCORE_ASSIGNMENT_CUTOFF)
            #           #+".wH"+str(args.H_weight)+".wC"+str(args.C_weight)+".tolH"+str(args.tolH)+".tolC"+str(args.tolC)+".rtolH"+str(args.rtolH)+".rtolN"+str(args.rtolN), 'w') as f:
            # with open("amino_acid_type_prediction_probabilities", 'w') as f:
            #     f.write("i AAIG\tpossible i-1 aa types\n")
            #     print("Probabilities:")
            #     for i_AAIG in list(iAAIG_iminus1aaTypesProbTupleList_dict.keys()):
            #         print(i_AAIG,"---> (i-1) aa type",sorted(iAAIG_iminus1aaTypesProbTupleList_dict[i_AAIG], key=itemgetter(1), reverse=True))
            #         f.write(i_AAIG + " ---> (i-1) aa type " + str(sorted(iAAIG_iminus1aaTypesProbTupleList_dict[i_AAIG], key=itemgetter(1), reverse=True)) + "\n")
            Probability.save_aa_type_predictions_to_file(iAAIG_iminus1aaTypesProbTupleList_dict,
                                                         "amino_acid_type_prediction_probabilities",
                                                         spectrum_type="TOCSY")


            ## OBSOLETE!
            # # Save amino acid type predictions that are above the Z-Score cutoff, along with the respective probabilities
            # #with open("amino_acid_type_prediction_probabilities.mcutoff"+str(args.RESONANCE_MATCH_CUTOFF)+".acutoff"+str(args.ASSIGNMENT_CUTOFF)+".zacutoff"+str(args.ZSCORE_ASSIGNMENT_CUTOFF)
            #           #+".wH"+str(args.H_weight)+".wC"+str(args.C_weight)+".tolH"+str(args.tolH)+".tolC"+str(args.tolC)+".rtolH"+str(args.rtolH)+".rtolN"+str(args.rtolN), 'w') as f:
            # with open("amino_acid_type_prediction_probabilities_above_cutoff.zacutoff"+str(args.ZSCORE_ASSIGNMENT_CUTOFF), 'w') as f:
            #     f.write("i AAIG\tpossible i-1 aa types\n")
            #     print("Probabilities:")
            #     for i_AAIG in list(iAAIG_iminus1aaTypesCutoffProbTupleList_dict.keys()):
            #         print(i_AAIG,"---> (i-1) aa type",sorted(iAAIG_iminus1aaTypesCutoffProbTupleList_dict[i_AAIG], key=itemgetter(1), reverse=True))
            #         f.write(i_AAIG + " ---> (i-1) aa type " + str(sorted(iAAIG_iminus1aaTypesCutoffProbTupleList_dict[i_AAIG], key=itemgetter(1), reverse=True)) + "\n")

            Probability.save_aa_type_predictions_to_file(iAAIG_iminus1aaTypesCutoffProbTupleList_dict,
                                "amino_acid_type_prediction_probabilities_above_cutoff.zacutoff"+str(args.ZSCORE_ASSIGNMENT_CUTOFF),
                                spectrum_type="TOCSY")

            ## OBSOLETE!
            # # Save amino acid type predictions with the respective Z-scores
            # #with open("amino_acid_type_prediction_Z-scores.mcutoff"+str(args.RESONANCE_MATCH_CUTOFF)+".acutoff"+str(args.ASSIGNMENT_CUTOFF)+".zacutoff"+str(args.ZSCORE_ASSIGNMENT_CUTOFF)+
            #           #".wH"+str(args.H_weight)+".wC"+str(args.C_weight)+".tolH"+str(args.tolH)+".tolC"+str(args.tolC)+".rtolH"+str(args.rtolH)+".rtolN"+str(args.rtolN), 'w') as f:
            # with open("amino_acid_type_prediction_Z-scores.zacutoff"+str(args.ZSCORE_ASSIGNMENT_CUTOFF), 'w') as f:
            #     f.write("i AAIG\tpossible i-1 aa types\n")
            #     # print("Z-scores:")
            #     for i_AAIG in list(iAAIG_iminus1aaTypesZscoreTupleList_dict.keys()):
            #         # print(i_AAIG,"---> (i-1) aa type",sorted(iAAIG_iminus1aaTypesZscoreTupleList_dict[i_AAIG], key=itemgetter(1), reverse=True))
            #         f.write(i_AAIG + " ---> (i-1) aa type " + str(sorted(iAAIG_iminus1aaTypesZscoreTupleList_dict[i_AAIG], key=itemgetter(1), reverse=True)) + "\n")
            Probability.save_aa_type_predictions_to_file(iAAIG_iminus1aaTypesZscoreTupleList_dict,
                                                        "amino_acid_type_prediction_Z-scores.zacutoff" + str(
                                                         args.ZSCORE_ASSIGNMENT_CUTOFF),
                                                         spectrum_type="TOCSY")

        else:
            ##########################################################################################################################################################################
            ##  LOAD AMINO ACID TYPE PREDICTIONS FOR FASTER DEBUGGING
            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            print("")    # change line after the previous progress bar
            print("Loading amino acid prediction files ...")
            # minimum_Zscore = 10000000

            ## OBSOLETE!
            # iAAIG_iminus1aaTypesProbPoolTupleList_dict = OrderedDict() # ordereddict with keys the AAIG of residue i and values lists of tuples of the form
            #                                                             # (residue i-1 matching aa type, Z-score)
            # with open(args.POOL_AA_TYPES_FILE, 'r') as f:
            #     pool_aa_type_file_contents = f.readlines()
            #     for line in pool_aa_type_file_contents[1:]:
            #         #print "DEBUG: line=",line
            #         word_list = re.sub('---> \(i-1\) aa type', '', line).split()
            #         key = word_list[0]
            #         iAAIG_iminus1aaTypesProbPoolTupleList_dict[key] = []
            #         values_string = ''.join(word_list[1:])
            #         elements_string = re.sub('[\(\)\[\]\'\"]', '',  values_string ).split("),(")[0]
            #         elements_list = elements_string.split(",")
            #         for aa, Zscore in zip(elements_list[0::2], elements_list[1::2]):
            #             #print "aa=",aa,"probability=",probability
            #             duplet = (aa, float(Zscore))
            #             iAAIG_iminus1aaTypesProbPoolTupleList_dict[key].append(duplet)
            iAAIG_iminus1aaTypesProbPoolTupleList_dict = Probability.load_aatype_probs_from_file(fname=args.POOL_AA_TYPES_FILE)

            # ## Now Calculate Z-scores from the pool of aa type prediction probabilites and keep only those predictions above
            # ## the cutoff
            ## OBSOLETE!
            # iAAIG_iminus1aaTypesZscoreTupleList_dict = OrderedDict() # ordereddict with keys the AAIG of residue i and values lists of tuples of the form (residue i-1 matching aa type, Z-score)
            # iAAIG_iminus1aaTypesCutoffProbTupleList_dict = OrderedDict() # contains only the tuples (residue i-1 matching aa type, probability) with Z-score above the cutoff
            # for iAAIG, iminus1aaTypesProbTuple_list in list(iAAIG_iminus1aaTypesProbPoolTupleList_dict.items()):
            #     # DECIDE THE PROBABILITY THRESHOLD TO DISCARD AA TYPE PREDICTIONS BELOW IT
            #     try:
            #         max_prob = np.max([duplet[1] for duplet in iminus1aaTypesProbTuple_list])   # highest aa type probability
            #     except ValueError:  # if no aa type predictions exist for this TAAIG group
            #         max_prob = 0
            #     if max_prob > 10e-10:   # Decide the probability threshold
            #         prob_threshold = 1000
            #     elif max_prob > 10e-20:
            #         prob_threshold = 10000
            #     elif max_prob <= 10e-20:
            #         prob_threshold = 100000
            #     print("DEBUG: max_prob", max_prob, "prob_threshold", prob_threshold, "args.LOG_TRANSFORM=", args.LOG_TRANSFORM)
            #     aatype_list, prob_list = [], []
            #     for duplet in iminus1aaTypesProbTuple_list:
            #         if args.DELETE_AA_TYPE_PREDICTIONS == True:     # delete low probability aa type predictions
            #             try:
            #                 if max_prob/float(duplet[1]) < prob_threshold:
            #                     aatype_list.append(duplet[0])
            #                     prob_list.append(duplet[1])  # unsorted list of probabilities
            #             except ZeroDivisionError:
            #                 print("DEBUG: ZeroDivisionError ", max_prob, "/", float(duplet[1]))
            #         elif args.DELETE_AA_TYPE_PREDICTIONS == False:
            #             aatype_list.append(duplet[0])
            #             prob_list.append(duplet[1])  # unsorted list of probabilities
            #     print("DEBUG: prob_list=", prob_list)
            #     if len(prob_list) == 1:
            #         zscore_array = np.array([10.00000])     # default Z-score in the case of just one prediction
            #     elif len(prob_list) > 1:
            #         if args.LOG_TRANSFORM == True:
            #             try:
            #                 ratio = float(np.max(prob_list))/np.min(prob_list)    # if the min probability is at least 3 orders of magnitude smaller, convert them to logarithmic scale
            #             except ZeroDivisionError:
            #                 print("DEBUG: ZeroDivisionError ", np.max(prob_list), "/", np.min(prob_list))
            #                 ratio = -1.0
            #             if ratio > 1000:
            #                 zscore_array = zscore(np.log(prob_list))  # convert probabilities to logarithms and then to Z-scores
            #             else:
            #                 zscore_array = zscore(prob_list)  # convert probabilities to Z-scores
            #         elif args.LOG_TRANSFORM == False:
            #             zscore_array = zscore(prob_list)  # convert probabilities to Z-scores
            #     else:   # if no aa type prediction was made create an empty array
            #         zscore_array = np.array([])
            #     # print "DEBUG1: saving zscore_array = ", zscore_array
            #     iminus1aaTypesZscoreTuple_list = []
            #     iminus1aaTypesCutoffProbTuple_list = []
            #     for aatype, Zscore, prob in zip(aatype_list, zscore_array, prob_list):
            #         if Zscore < minimum_Zscore:
            #             minimum_Zscore = Zscore
            #         ## ONLY IF THE Z-SCORE OF THE AA TYPE PREDICTION IS GREATER THAN THE CUTOFF AND THE AA TYPE PREDICTIONS ARE AT LEAST args.MIN_NUM_OF_PREDICTIONS, THEN INCLUDE IT IN THE LIST
            #         if Zscore > args.ZSCORE_ASSIGNMENT_CUTOFF and len(zscore_array) >= args.MIN_NUM_OF_PREDICTIONS:
            #             aatypeProb_tuple = (aatype, prob)
            #             aatypeZscore_tuple = (aatype, Zscore)  # replace the probability with the respective Z-score
            #             iminus1aaTypesZscoreTuple_list.append(aatypeZscore_tuple)
            #             iminus1aaTypesCutoffProbTuple_list.append(aatypeProb_tuple) # but also save the aatype probabilities that satisfy the Z-score cutoff
            #         ## IF THE AA TYPE PREDICTIONS ARE LESS THAN args.MIN_NUM_OF_PREDICTIONS, THEN INCLUDED ALL OF THEM IN THE LIST
            #         elif len(zscore_array) < args.MIN_NUM_OF_PREDICTIONS:
            #             aatypeProb_tuple = (aatype, prob)
            #             aatypeZscore_tuple = (aatype, Zscore)  # replace the probability with the respective Z-score
            #             iminus1aaTypesZscoreTuple_list.append(aatypeZscore_tuple)
            #             iminus1aaTypesCutoffProbTuple_list.append(aatypeProb_tuple) # but also save the aatype probabilities that satisfy the Z-score cutoff
            #     iAAIG_iminus1aaTypesZscoreTupleList_dict[iAAIG] = iminus1aaTypesZscoreTuple_list
            #     iAAIG_iminus1aaTypesCutoffProbTupleList_dict[iAAIG] = iminus1aaTypesCutoffProbTuple_list

            # Use pool of aa type prediction probabilities to calculate z-scores
            iAAIG_iminus1aaTypesCutoffProbTupleList_dict, \
            iAAIG_iminus1aaTypesZscoreTupleList_dict = \
                Probability.get_pruned_prob_dicts_by_zscore(iAAIG_iminus1aaTypesProbPoolTupleList_dict,
                                                            ZSCORE_ASSIGNMENT_CUTOFF=args.ZSCORE_ASSIGNMENT_CUTOFF,
                                                            DELETE_AA_TYPE_PREDICTIONS=args.DELETE_AA_TYPE_PREDICTIONS,
                                                            LOG_TRANSFORM=args.LOG_TRANSFORM,
                                                            MIN_NUM_OF_PREDICTIONS=args.MIN_NUM_OF_PREDICTIONS)

            ## OBSOLETE!
            # iAAIG_iminus1aaTypesProbTupleList_dict = OrderedDict()   # ordereddict with keys the AAIG of residue i and values lists of tuples of the form
            #                                                             # (residue i-1 matching aa type, average probability)
            # with open(args.COMPLETE_AA_TYPES_FILE, 'r') as f:
            #     complete_aa_type_file_contents = f.readlines()
            #     for line in complete_aa_type_file_contents[1:]:
            #         #print "DEBUG: line=",line
            #         word_list = re.sub('---> \(i-1\) aa type', '', line).split()
            #         key = word_list[0]
            #         iAAIG_iminus1aaTypesProbTupleList_dict[key] = []
            #         values_string = ''.join(word_list[1:])
            #         elements_string = re.sub('[\(\)\[\]\'\"]', '',  values_string ).split("),(")[0]
            #         elements_list = elements_string.split(",")
            #         for aa, hist_prob in zip(elements_list[0::2], elements_list[1::2]):
            #             #print "aa=",aa,"probability=",probability
            #             duplet = (aa, float(hist_prob))
            #             iAAIG_iminus1aaTypesProbTupleList_dict[key].append(duplet)

            iAAIG_iminus1aaTypesProbTupleList_dict = Probability.load_aatype_probs_from_file(fname=args.COMPLETE_AA_TYPES_FILE)

            #for i_AAIG in iAAIG_iminus1aaTypesZscoreTupleList_dict.keys():
            #    print i_AAIG,"---> (i-1) aa type",sorted(iAAIG_iminus1aaTypesZscoreTupleList_dict[i_AAIG], key=itemgetter(1), reverse=True)
            ##########################################################################################################################################################################

            ## OBSOLETE!
            # with open("amino_acid_type_prediction_Z-scores.zacutoff"+str(args.ZSCORE_ASSIGNMENT_CUTOFF), 'w') as f:
            #     f.write("i AAIG\tpossible i-1 aa types\n")
            #     print("Z-scores:")
            #     for i_AAIG in list(iAAIG_iminus1aaTypesZscoreTupleList_dict.keys()):
            #         print(i_AAIG, "---> (i-1) aa type",sorted(iAAIG_iminus1aaTypesZscoreTupleList_dict[i_AAIG], key=itemgetter(1), reverse=True))
            #         f.write(i_AAIG + " ---> (i-1) aa type " + str(sorted(iAAIG_iminus1aaTypesZscoreTupleList_dict[i_AAIG], key=itemgetter(1), reverse=True)) + "\n")

            Probability.save_aa_type_predictions_to_file(iAAIG_iminus1aaTypesZscoreTupleList_dict,
                                        fname="amino_acid_type_prediction_Z-scores.zacutoff"+str(args.ZSCORE_ASSIGNMENT_CUTOFF),
                                                         spectrum_type="TOCSY")


        #sys.exit(1)

        ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        ##              CALCULATE BAYESIAN STATISTICS                        #
        ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

        prob = Probability(fasta_file=args.FASTA_FILE)
        prob.calc_TOCSY_aa_type_conditional_probs(iAAIG_iminus1aaTypesProbTupleList_dict,
                                                         transform_type=args.TRANSFORM_TYPE)
        Probability.save_aa_type_predictions_to_file(prob.iAAIG_iminus1aaType_Condprob_mdict,
            fname="amino_acid_type_prediction_conditional_probabilities",
            spectrum_type="TOCSY")  # Save ALL amino acid type predictions with the respective conditional probabilities

        ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        ##       BUILD TREES OF POSSIBLE PEPTIDE SEQUENCES                   #
        ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        
        if not args.SKIP_CHAINS:

            # NEW WAY with a Peptide object
            # Eventually the Conditional AA-type Probabilities will be used to derive the peptide scores.
            pept = Peptide(total_chain_number=len(all_chainScore_list),
                 min_peptide_length=args.MIN_PEPTIDE_LENGTH,
                 AAIG_aaType_Condprob_mdict=prob.iAAIG_iminus1aaType_Condprob_mdict,
                 AAIG_aaTypesCutoffProbTupleList_dict=iAAIG_iminus1aaTypesCutoffProbTupleList_dict,
                 AAIG_aaTypesTransformedProbTupleList_dict=prob.iAAIG_iminus1aaTypesTransformedProbTupleList_dict,
                 aatype_P_dict=prob.aatype_P_dict,
                 connectivity_dict=i_iminus1_dict)
            pept.build_all_peptide_trees(all_chainScore_list, resume=args.RESUME, is_parallel=True)
            # Repeat peptide building for a second time to create whatever was left out accindentally.
            pept.build_all_peptide_trees(all_chainScore_list, resume=True, is_parallel=True)
            # Finally extract the peptide from the pickle files and write them to text files
            Peptide.write_all_peptides_to_files(all_chainScore_list, max_peptide_length=args.MAX_PEPTIDE_LENGTH)
    
    except:
        type, value, tb = sys.exc_info()
        lines = traceback.format_exception(type, value, tb)
        print((''.join(lines)))
        raise