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



import sys, re, os, cPickle, traceback, shutil, bz2, math
from scoop import futures, shared
import numpy as np
from operator import itemgetter
from ordereddict import OrderedDict
from ete3 import Tree
import ftplib
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from scipy.stats.mstats import zscore
from scipy import stats, sqrt
import collections
import gc
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
from open_func import *


HOME_DIR = os.path.dirname(os.path.realpath(__file__))


def cmdlineparse():
    parser = ArgumentParser(description="command line arguments", formatter_class=RawDescriptionHelpFormatter,
                            epilog="EXAMPLE: 4D_assignment_tree-based.py -root TDhsqc.list -tocsy TDtocsy.list -noesy TDnoesy.list -tolH 0.01 -tolC 0.1 -rtolH 0.01 -rtolN 0.1 -mcutoff 1.0 -acutoff 1.0 -wH 0.5 -wC 1")
    parser.add_argument("-tseq", dest="template_sequence_file", required=True,
                        help="teplate sequence file in fasta format", metavar="<template sequence file>")
    parser.add_argument("-root", dest="ROOT_fname", required=True,
                        help="2D N-H HSQC root spectrum", metavar="<2D N-H HSQC input file>")
    parser.add_argument("-tocsy", dest="TOCSY_fname", required=True,
                        help="4D TOCSY (HCTOCSYNH) file", metavar="<4D TOCSY input file>")
    parser.add_argument("-noesy", dest="NOESY_fname", required=True,
                        help="4D NOESY (HCNOENH) file", metavar="<4D NOESY input file>")
    parser.add_argument("-tolH", dest="tolH", required=False, type=float, default=0.04,
                        help="tolerance for the proton resonance when matching NOESY peaks to TOCSY peaks. (default %(default)s)",
                        metavar="<proton tolerance>")
    parser.add_argument("-tolC", dest="tolC", required=False, type=float, default=0.4,
                        help="tolerance for the carbon resonance when matching NOESY peaks to TOCSY peaks. (default %(default)s)",
                        metavar="<carbon tolerance>")
    parser.add_argument("-update", dest="DOWNLOAD_CS_HISTOGRAMS", required=False, action='store_true',
                        help="download the latest chemical shift histgrams from BMRB")
    parser.add_argument("-wH", dest="H_weight", required=False, type=float, default=0.1,
                        help="weight in (0.0-1.0] to apply on aa type prediction from aliphatic H resonances. (default %(default)s)",
                        metavar="<H weight>")
    parser.add_argument("-wC", dest="C_weight", required=False, type=float, default=1.0,
                        help="weight in (0.0, 1.0] to apply on aa type prediction from aliphatic C resonances. (default %(default)s)",
                        metavar="<C weight>")
    parser.add_argument("-mcutoff", dest="RESONANCE_MATCH_CUTOFF", required=False, type=float, default=0.8,
                        help="number 0.0-1.0 saying how much of TOCSY resonances should match in the NOESY in order to consider it a possible match. (default %(default)s)",
                        metavar="<resonance match cutoff>")
    parser.add_argument("-zmcutoff", dest="ZSCORE_RESONANCE_MATCH_CUTOFF", required=False, type=float, default=0.0,
                        help="a real number specifying the lower Z-score for a resonance match to be retained for chain building. (default %(default)s)",
                        metavar="<Z-score resonance match cutoff>")
    parser.add_argument("-acutoff", dest="ASSIGNMENT_CUTOFF", required=False, type=float, default=None,
                        help="number 0.0-1.0 saying how much of TOCSY resonances should match with a particular aa type in order to be considered, e.g. \
                        if TOCSY resonances are 5, 0.8 will mean that at least 4 resonances must match. (default %(default)s)",
                        metavar="<aa type assignment cutoff>")
    parser.add_argument("-zacutoff", dest="ZSCORE_ASSIGNMENT_CUTOFF", required=False, type=float, default=-1.0,
                        help="a real number specifying the lower Z-score for an aa type prediction to be considered as valid. (default %(default)s)",
                        metavar="<Z-score aa type assignment cutoff>")
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
    
    parser.add_argument("-poolconfile", dest="POOL_CONNECTIVITIES_FILE", required=False, default=None,
                        help="connectivities pool file; necessary if -confile specified in order to calculate correct probabilities", metavar="<connectivities pool file>")
    parser.add_argument("-allconfile", dest="COMPLETE_CONNECTIVITIES_FILE", required=False, default=None,
                        help="all connectivities file; necessary if -confile specified in order to calculate correct probabilities", metavar="<all connectivities file>")
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
                        help="type of mathematical transform to apply on the probabilities P[Tindex(i)|aatype(i-1)]. Allowed values are: \"None\", \"log\", \"log10\", \"boxcox_pearson\", \"boxcox_mle\"",
                        metavar="<type of mathematical transform>")
    parser.add_argument("-resume", dest="RESUME", required=False, action='store_true', default=False,
                        help="resume previous run by loading all peptide sequences saved in tmp_peptide_folder/. By default tmp_peptide_folder/ will be cleaned and new peptide sequences will be writen inside it.")
    parser.add_argument("-probprod", dest="PROB_PRODUCT", required=False, action='store_true', default=False,
                        help="select the best C-H type assignment combination based on the product of probabilities of the individual C-H assignments")
    parser.add_argument("-probmodel", dest="PROBABILITY_MODEL", required=False, type=int, default=2,
                        help="If '1' the probability of each peak will be given by [wH*1Dhist(H)+wC*1Dhist(C)]/(wH+wC). If '2' then by 1Dhist(H)*1Dhist(C)")
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
                        When the C-group contains only one peak, the respective probability will be prob1 in all cases.
                            """, metavar="<way to calculate the cs assignment score>")
    parser.add_argument("-2dhist", dest="USE_2D_HISTOGRAMS", required=False, action='store_true', default=False,
                        help="use 2D BMRB histograms for aa type prediction")
    parser.add_argument("-log", dest="LOG", required=False, action='store_true', default=False,
                        help="convert aa type prediction probabilities to logarithmic scale and then calculate Z-scores, if the min(probability)/max(probability) > 1000")
    parser.add_argument("-skipchains", dest="SKIP_CHAINS", required=False, action='store_true', default=False,
                        help="Generate connectivies and amino acid type predictions, but skip chain and peptide formation.")
    parser.add_argument("-delpred", dest="DELETE_AA_TYPE_PREDICTIONS", required=False, action='store_true', default=False,
                        help="delete aa type predictions with probabilities which are 1000, 10000 or 100000 times lower than the highest, if the highest is >10e-10, >10e-20, <=10e-20, respectively.")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s version 1.0')
    args=parser.parse_args()
    return args



def check_input_files_format():
    global args
    
    query_contents = []
    HAS_INTENSITIES = False
    with open(args.NOESY_fname, 'r') as f:
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
                print "DEBUG check_input_files_format: this NOESY line will be discarded:", line
                continue
    if HAS_INTENSITIES:
        for word_list in query_contents:
            if len(word_list) < 6:
                print "ERROR: Please check your input NOESY file. Not all lines have intensity!!!"
                sys.exit(1)



def transform(probability_array):
    """
        FUNCTION to do a mathematical transform on a value
    """
    if args.TRANSFORM_TYPE == "None":
        return probability_array
    elif args.TRANSFORM_TYPE == "log":
        return np.log(probability_array)
    elif args.TRANSFORM_TYPE == "log10":
        return np.log10(probability_array)
    elif args.TRANSFORM_TYPE == "boxcox_pearson":
        lmax_pearsonr = stats.boxcox_normmax(probability_array)
        prob_pearson = stats.boxcox(probability_array, lmbda=lmax_pearsonr)
        return prob_pearson
    elif args.TRANSFORM_TYPE == "boxcox_mle":
        prob_mle, lmax_mle = stats.boxcox(probability_array)
        return prob_mle


def chunkIt(seq, num):
    """
        split a list into a specified number of approximately equal sublists.
    """
    avg = len(seq) / float(num)
    out = []
    last = 0.0
  
    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg
  
    return out


def remove_redundant_chains(chainScore_list):
    
    progbar = shared.getConst('PROGBAR')
    all_chainScore_list = shared.getConst('ALL_CHAINSCORE_LIST')
    chainScores2remove_set = set()
    chain_num = len(chainScore_list)
    for i, chainScore1 in enumerate(chainScore_list):
        chain1 = chainScore1[0:-1]
        score1 = chainScore1[-1]
        progbar.set_progress((float(i)/chain_num))
        for chainScore2 in all_chainScore_list:
            chain2 = chainScore2[0:-1]
            score2 = chainScore2[-1]
            if is_sublist(chain1, chain2) and score1 == score2: # this is the correct one
                chainScores2remove_set.add(tuple(chainScore1))   # make first the list immutable to be able to insert it into a set
    
    return chainScores2remove_set


def save_peptides_to_file(peptide_files_list, part):
    
    progbar = shared.getConst('PROGBAR')
    args = shared.getConst('ARGS')
    all_chainScore_list = shared.getConst('NEW_ALL_CHAINSCORE_LIST')
    f = open("peptides."+str(args.MAX_PEPTIDE_LENGTH)+"mers.list.chunk"+part, 'w')
    index = 0
    peptide_file_index = 0
    peptide_file_num = len(peptide_files_list)
    for peptide_file in peptide_files_list:
        mo = re.search('chainIndex_([0-9]+).pickle.bz2', peptide_file)
        if mo:
            peptide_file_index += 1
            progbar.set_progress((float(peptide_file_index)/peptide_file_num))
            chainIndex = int(mo.group(1))
            try:
                with bz2.BZ2File('tmp_peptide_folder/' + peptide_file, "rb") as pickled_file:
                    pickle_object = cPickle.load(pickled_file)
            except EOFError:
                print "ERROR: file tmp_peptide_folder/" + peptide_file, " seems to be corrupted perhaps due to abrupt termination of the program before the file was properly saved.\
                Please run again the program without the -resume option."
                sys.exit(1)
            peptideScoreChainindexList_list = pickle_object[0]  # a list of lists containing the possible peptide sequences and the overall score as the last element
            peptideProbList_list = pickle_object[1] # a list of lists containing the probability that each amino acid of the peptide is correctly predicted from the chain
            for peptideScoreChainindexList, peptideProbList in zip(peptideScoreChainindexList_list, peptideProbList_list):
                chain = all_chainScore_list[chainIndex]
                peptide = peptideScoreChainindexList[0:-3]
                Cterm = peptideScoreChainindexList[-3]
                peptide_score = peptideScoreChainindexList[-2]
                chainIndex2 = peptideScoreChainindexList[-1]
                if chainIndex != chainIndex2:
                    print "ERROR: the chainIndex of file ", peptide_file, " should be ", chainIndex2, ", not ", chainIndex
                    sys.exit(1)
                if peptide_score == 0.0:    # do not save peptides containing aa types that are not present in the protein
                    continue
                f.write("Chain probability = "+str(chain[-1])+" peptide P_"+str(len(peptide))+"_"+str(index+1) + " = "+''.join(peptide)+" "+'-'.join(chain[0:-1])+" "+','.join(map(str, peptideProbList[:-1]))+"\n")
                index += 1
    f.close()


def split_oversized_Cgroups(possible_aatype_prob_C_H_resonpair_Tindex_list_list):
    """
        FUNCTION to split Cgroups with more that 2 peaks into separate Cgroups. Currently works for maximum Cgroup size 3 only. This happens when
        two peaks have the same Carbon resonance but different proton resonance. So because C-grouping works with Carbon resonances only, all 3 peaks
        are included in the same group
    """
    
    aatypes_set = set([p[0] for p in possible_aatype_prob_C_H_resonpair_Tindex_list_list])  # must be only one for CS assignment scripts
    corrected_possible_aatype_prob_C_H_resonpair_Tindex_list_list = []
    for aatype in aatypes_set:
        aaspec_possible_aatype_prob_C_H_resonpair_Tindex_list_list = [p for p in possible_aatype_prob_C_H_resonpair_Tindex_list_list if p[0]==aatype] # aatype-specific possible_aatype_prob_C_H_resonpair_Tindex_list_list
        clustID2size_dict = {}
        clustID2Creson_dict = {}
        for a in aaspec_possible_aatype_prob_C_H_resonpair_Tindex_list_list:
            clustID = a[7]
            Creson = a[6]
            Hreson = a[5]
            if not clustID in clustID2Creson_dict.keys():
                clustID2Creson_dict[clustID] = set([(Hreson, Creson)])
            else:
                clustID2Creson_dict[clustID].add((Hreson, Creson))
        
        for clustID in clustID2Creson_dict.keys():
            clustID2size_dict[clustID] = len(clustID2Creson_dict[clustID])
        
        maxClustID = max(clustID2size_dict.keys())
        for clustID in clustID2size_dict.keys():
            if clustID2size_dict[clustID] == 3:
                Creson_list = [a[6] for a in aaspec_possible_aatype_prob_C_H_resonpair_Tindex_list_list if a[7]==clustID]
                dist12 = abs(Creson_list[0] - Creson_list[1])
                dist13 = abs(Creson_list[0] - Creson_list[2])
                dist23 = abs(Creson_list[1] - Creson_list[2])
                if dist12 < dist13 and dist12 < dist23: # change the Cgroup of peak3
                    Coutlier = Creson_list[2]
                    for i in range(len(aaspec_possible_aatype_prob_C_H_resonpair_Tindex_list_list)):
                        a = aaspec_possible_aatype_prob_C_H_resonpair_Tindex_list_list[i]
                        if a[6] == Coutlier:
                            a[7] = maxClustID + 1
                            maxClustID += 1 # update the highest Cgroup number
                elif dist13 < dist12 and dist13 < dist23:   # change the Cgroup of peak2
                    Coutlier = Creson_list[1]
                    for i in range(len(aaspec_possible_aatype_prob_C_H_resonpair_Tindex_list_list)):
                        a = aaspec_possible_aatype_prob_C_H_resonpair_Tindex_list_list[i]
                        if a[6] == Coutlier:
                            a[7] = maxClustID + 1
                            maxClustID += 1 # update the highest Cgroup number
                elif dist23 < dist12 and dist23 < dist13:   # change the Cgroup of peak1
                    Coutlier = Creson_list[0]
                    for i in range(len(aaspec_possible_aatype_prob_C_H_resonpair_Tindex_list_list)):
                        a = aaspec_possible_aatype_prob_C_H_resonpair_Tindex_list_list[i]
                        if a[6] == Coutlier:
                            a[7] = maxClustID + 1
                            maxClustID += 1 # update the highest Cgroup number
                    
            elif clustID2size_dict[clustID] > 3:
                print "ERROR: found Cgroup with more than 3 peaks!!!"
                print "DEBUG: aaspec_possible_aatype_prob_C_H_resonpair_Tindex_list_list=", aaspec_possible_aatype_prob_C_H_resonpair_Tindex_list_list
                sys.exit(1)
        corrected_possible_aatype_prob_C_H_resonpair_Tindex_list_list.extend(aaspec_possible_aatype_prob_C_H_resonpair_Tindex_list_list)
    
    return corrected_possible_aatype_prob_C_H_resonpair_Tindex_list_list


def group_carbons(possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list, iteration=None):
    """
        This is a different function for the group_carbons() used in the CS assignment scripts, because there we know the aa type, but here we
        group the carbons for every possible aa type.
        ARGUMENTS:
        possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list:    list of lists of the form
                (aa type, Carbon name, Hydrogen name, weighted average probability, TOCSY_reson_index, H resonance, C resonance, carbon group)
        RETURNS:
        The same as the input but with the correct carbon groups (last element)
    """
    global aa_carbonBondedGeminalHydrogensDict_dict
    
    if iteration in [1,3]:  # assigned a different C-group to each peak
        HCreson_set = set([(p[5], p[6]) for p in possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list])
        clustID = 0
        for HCreson in HCreson_set:
            clustID += 1
            for peak in possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list:
                if peak[5]==HCreson[0] and peak[6]==HCreson[1]:
                    peak[7] = clustID
        return possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list
    
    print "DEBUG group_carbons: possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list=", possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list
    aa_types_set = set([x[0] for x in possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list])
    updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list = []    # same as possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list but with the cluster IDs for each aa type
    for aa_type in aa_types_set:
        singleAAtype_prob_C_H_resonpair_TOCSYindex_list_list = [x for x in possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list if x[0]==aa_type]
        if not aa_type in aa_carbonBondedGeminalHydrogensDict_dict.keys():  # this aa does not have geminal protons, skip clustering
            clustID = 0
            for index in range(len(singleAAtype_prob_C_H_resonpair_TOCSYindex_list_list)):
                clustID += 1
                singleAAtype_prob_C_H_resonpair_TOCSYindex_list_list[index][7] = clustID
            updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list.extend(singleAAtype_prob_C_H_resonpair_TOCSYindex_list_list)
            continue    # continue with the next possible aa type
        
        new_aaindex_carbonGroupsTuple_dict = {}
        print "DEBUG: before singleAAtype_prob_C_H_resonpair_TOCSYindex_list_list=", singleAAtype_prob_C_H_resonpair_TOCSYindex_list_list
        carbonCS_list = [l[6] for l in singleAAtype_prob_C_H_resonpair_TOCSYindex_list_list]
        print "DEBUG: carbonCS_list=", carbonCS_list
        cl = HierarchicalClustering(carbonCS_list , lambda x,y: float(abs(x-y)))
        
        CONTINUE_CLUSTERING = True
        cutoff = 0.200001       # <== CHANGE ME (arbitrary starting value)
        while CONTINUE_CLUSTERING:  # continue clustering by lowering the cutoff until no more than 2 carbons resonances are within each cluster
            CONTINUE_CLUSTERING = False
            cluster_list = cl.getlevel(cutoff)    # <== CHANGE ME (the tolerance for carbon)
            print "DEBUG: cluster_list=", cluster_list
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
                print "DEBUG: clustID=", clustID, "clust=", clust, "carbon_resonance=", carbon_resonance
                if len(cluster_list) == 1 and carbon_resonance == clust:    # to avoid "TypeError: argument of type 'float' is not iterable"
                    print "DEBUG: assigning cluster", clustID+1
                    singleAAtype_prob_C_H_resonpair_TOCSYindex_list_list[index][7] = clustID + 1  # because enumeration starts from 0
                elif carbon_resonance in clust:
                    print "DEBUG: assigning cluster", clustID+1
                    singleAAtype_prob_C_H_resonpair_TOCSYindex_list_list[index][7] = clustID + 1  # because enumeration starts from 0
        updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list.extend(singleAAtype_prob_C_H_resonpair_TOCSYindex_list_list)
        print "DEBUG group_carbons: appended carbon clusters with ", singleAAtype_prob_C_H_resonpair_TOCSYindex_list_list
    
    
    print  "DEBUG: point 1 updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list=", updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list
    updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list = split_oversized_Cgroups(updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list)
    
    print  "DEBUG: point 2 updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list=", updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list
    if len(possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list[0]) == 9:    # if this is NOESY assignment (there is an 9th columns with the spectrum source)
        matched_peaks = [x for x in possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list if x[8]=='TOCSY-NOESY']
        matched_Cnames = [x[1] for x in matched_peaks]
        matched_clustIDs = [x[7] for x in matched_peaks]
        for Cname, clustID in zip(matched_Cnames, matched_clustIDs):
            for index in range(len(updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list)):
                peak = updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list[index]
                if peak[8] in ['TOCSY-NOESY', 'TOCSY (unmatched)'] and peak[1]==Cname:
                    updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list[index][7] = clustID # change the cluster ID
    
    print  "DEBUG: point 3 updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list=", updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list
    for aa_type in aa_types_set:
        CH_reson_tuple_set = set([(p[5],p[6]) for p in updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list if p[0]==aa_type])    # all unique C-H resonances
        for CH_reson_tuple in CH_reson_tuple_set:
            Creson = CH_reson_tuple[1]
            Hreson = CH_reson_tuple[0]
            clusterIDs_set = set([p[7] for p in updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list if p[5]==Hreson and p[6]==Creson and p[0]==aa_type])
            print "DEBUG: clusterIDs_set=", clusterIDs_set
            common_clusterID = min(clusterIDs_set) # give them the same cluster ID (by default the smallest)
            for i in range(len(updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list)):
                if updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list[i][0]==aa_type and updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list[i][5]==Hreson and updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list[i][6]==Creson:
                    updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list[i][7] = common_clusterID    # reinitialize the cluster ID
    
    print  "DEBUG: point 4 updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list=", updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list
    reupdated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list = []
    for aa_type in aa_types_set:
        clusterIDs_set = set([p[7] for p in updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list if p[0]==aa_type])   # all the cluster IDs
        for clustID in clusterIDs_set:
            CH_reson_tuple_set = set([(p[5],p[6]) for p in updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list if p[0]==aa_type and p[7]==clustID])    # get all peaks of this cluster
            print "DEBUG: CH_reson_tuple_set=", CH_reson_tuple_set
            predictions_list = [p for p in updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list if p[0]==aa_type and p[7]==clustID]    # all C-H predictions of this cluster
            if len(CH_reson_tuple_set) == 2:    # if there are 2 peaks in this C-group, check if they have been predicted to be geminal
                for pred in predictions_list:
                    Cname = pred[1]
                    print "DEBUG: there are 2 peaks in this C-group, Cname=", Cname
                    if not Cname in aatype_carbon_nongeminalHname_multidict[aa_type].keys() or (aa_type in ["LEU", "VAL"] and Cname in ["CD1", "CD2", "CG1", "CG2"]):    # if this is a Carbon with geminal protons, keep the prediction, otherwise discard it
                        print "DEBUG: saving prediction:", pred
                        reupdated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list.append(pred)
            else:   # if the C-group contains only one peak, then it can be also a methylene, so save the prediction
                print "DEBUG: there is only 1 peak in this C-group"
                for pred in predictions_list:
                    print "DEBUG: saving prediction:", pred
                    reupdated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list.append(pred)
    
    print "DEBUG group_carbons: reupdated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list=", reupdated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list
    return reupdated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list



if __name__ == "__main__":
    
    try:
        
        args = cmdlineparse()
        if args.H_weight <= 0 or args.H_weight > 1:
            print "ERROR: -wH must be a number greater than 0 and lower or equal to 1!"
            sys.exit(1)
        if args.C_weight <= 0 or args.C_weight > 1:
            print "ERROR: -wC must be a number greater than 0 and lower or equal to 1!"
            sys.exit(1)
        if (args.ASSIGNMENT_CUTOFF != None and args.POOL_AA_TYPES_FILE != None) or (args.ASSIGNMENT_CUTOFF != None and args.COMPLETE_AA_TYPES_FILE != None):
            print "ERROR: Usage of -poolaafile or -allaafile is incompatible with -acutoff argument! Either remove -poolaafile and -allaafile to make new amino acid \
            type predictions from scratch, or remove -acutoff to use the amino acid prediction files you provided as input."
            sys.exit(1)
        if args.ASSIGNMENT_CUTOFF == None and args.POOL_AA_TYPES_FILE == None and args.COMPLETE_AA_TYPES_FILE == None:
            print "ERROR: you must either provide an aa type assignment cutoff with -acutoff argument (recommented value: 1.0), or provide a pool amino acid assignment file with \
            argument -poolaafile and file with all the aa type assignment with -allaafile argument."
            sys.exit(1)
        
        check_input_files_format()
        
        if args.DOWNLOAD_CS_HISTOGRAMS:
            download_CS_histograms()
        
        
        with open(args.ROOT_fname, 'r') as f:
            root_contents=f.readlines()
        remaining_root_contents = []    # list of Root spectrum lines with overlapping groups
        nonoverlapping_root_contents_and_tolerances = []    # list wiht the lines of the Root spectrum that contain non-overlapping groups and the associated rtoH, rtolN
        for rtolH, rtolN in zip([0.02, 0.01, 0.01, 0.005], [0.2, 0.1, 0.05, 0.015]):
            nonoverlapping_root_contents_and_tolerances.extend(find_nonoverlapping_groups_in_root(root_contents, rtolH, rtolN))
            for triplet in nonoverlapping_root_contents_and_tolerances:
                line = triplet[0]
                try:
                    root_contents.remove(line)
                except ValueError:
                    continue
        
        
        
        
        print "DEBUG: final nonoverlapping_root_contents_and_tolerances=", nonoverlapping_root_contents_and_tolerances
        TOCSY_lines, TOCSY_sidechain_resonances_list = copy_aaindices_from_root_spectrum_2(nonoverlapping_root_contents_and_tolerances, args.TOCSY_fname, "TOCSY")
        NOESY_lines = copy_aaindices_from_root_spectrum_2(nonoverlapping_root_contents_and_tolerances, args.NOESY_fname, "NOESY", TOCSY_sidechain_resonances_list)
        
        TOCSY_contents = []
        print "DEBUG: TOCSY_lines=", TOCSY_lines
        for TOCSY_line in TOCSY_lines:
            TOCSY_words_list = filter(lambda a: a!= '' ,re.split('\s+', TOCSY_line)[:-1]) # discard the last element which is ''
            try:
                if len(TOCSY_words_list)==5 and float(TOCSY_words_list[1]) and float(TOCSY_words_list[2]) and float(TOCSY_words_list[3]) and float(TOCSY_words_list[4]):
                    TOCSY_contents.append(TOCSY_words_list)
                elif len(TOCSY_words_list)==6 and float(TOCSY_words_list[1]) and float(TOCSY_words_list[2]) and float(TOCSY_words_list[3]) and float(TOCSY_words_list[4]) and float(TOCSY_words_list[5]):
                    TOCSY_contents.append(TOCSY_words_list)
            except ValueError:
                print "Discarding the following line from TOCSY:"
                TOCSY_line
                continue
        TOCSY_contents.sort(key=itemgetter(0))  # sort TOCSY lines by the random aa index (2nd column)
        
        
        
        
        def create_connectivity_dictionaries():
            connectivities_multidict = get_possible_connectivities(TOCSY_contents,NOESY_lines, args.tolH, args.tolC)
            i_iminus1_dict = {} # a dictionary containing all possible connectivities of every TOCSY aa index; has keys the TOCSY aa indices and values lists of triplets (tuples) \
            i_iminus1_complete_dict = {} # the same dictionary but with all possible connectivities, including those below args.RESONANCE_MATCH_CUTOFF
            for TOCSYaaindex in connectivities_multidict.keys():
                NOESYaaindex_list = connectivities_multidict[TOCSYaaindex].keys()
                if TOCSYaaindex in NOESYaaindex_list:
                    NOESYaaindex_list.remove(TOCSYaaindex)  # remove residue i from the list of possible residues (i-1)
                NOESYaaindex_occupancy_numOfResonances_tuple__list = [] 
                for NOESYaaindex in NOESYaaindex_list:
                    occupancy = connectivities_multidict[TOCSYaaindex][NOESYaaindex][0]
                    numOfResonances = connectivities_multidict[TOCSYaaindex][NOESYaaindex][1]
                    NOESYaaindex_occupancy_numOfResonances_tuple__list.append((NOESYaaindex, occupancy, numOfResonances))
                sorted_NOESYaaindex_occupancy_numOfResonances_tuple__list = sorted(NOESYaaindex_occupancy_numOfResonances_tuple__list, key=itemgetter(1, 2), reverse=True)
                for triplet in sorted_NOESYaaindex_occupancy_numOfResonances_tuple__list:
                    if float(triplet[1])/numOfResonances >= float(args.RESONANCE_MATCH_CUTOFF):     # adjustable cutoff
                        try:
                            i_iminus1_dict[TOCSYaaindex].append(triplet)
                        except KeyError:
                            i_iminus1_dict[TOCSYaaindex] = [(triplet)]
                    try:
                        i_iminus1_complete_dict[TOCSYaaindex].append(triplet)
                    except KeyError:
                        i_iminus1_complete_dict[TOCSYaaindex] = [(triplet)]
            
            new_i_iminus1_dict = {}
            for i in i_iminus1_dict.keys():
                prob_list = [float(triplet[1])/triplet[2] for triplet in i_iminus1_dict[i]]
                if len(prob_list) == 1:     
                    zscore_array = np.array([10.00000])     # default Z-score in the case of just one connectivity
                elif len(set(prob_list)) == 1:  # if all connectivities have the same value keep them all
                    zscore_array = np.array([10.00000]*len(prob_list))
                elif len(prob_list) > 2:
                    zscore_array = zscore(prob_list)  # convert probabilities to Z-scores
                elif len(prob_list) == 2:   # if only 2 connectivities exist, set the 2nd to 0.000 to allow you to use high -zcutoff values (<= 0)
                    zscore_array = np.array([1.00000, 0.00000])
                elif len(prob_list) == 0:   # if no connectivity was made create an empty array
                    zscore_array = np.array([])
                
                new_i_iminus1_dict[i] = []  # add i to the new dictionary
                for triplet, Zscore in zip(i_iminus1_dict[i], zscore_array):
                    if Zscore >= args.ZSCORE_RESONANCE_MATCH_CUTOFF or approx_equal(Zscore, args.ZSCORE_RESONANCE_MATCH_CUTOFF):
                        new_i_iminus1_dict[i].append(triplet)
            
            del i_iminus1_dict
            i_iminus1_dict = new_i_iminus1_dict
            del new_i_iminus1_dict
            
            return i_iminus1_dict, i_iminus1_complete_dict
        
        
        if not args.POOL_CONNECTIVITIES_FILE:
            print "Calculating Connectivities..."
            i_iminus1_dict, i_iminus1_complete_dict = create_connectivity_dictionaries()
            
        else:
            print "Loading Connectivities from file", args.POOL_CONNECTIVITIES_FILE
            i_iminus1_pool_dict = {} # same dictionary with the pool of possible connectivities remained after filtering using absolute consensus matches of a previous run
            i_iminus1_complete_dict = {} # the same dictionary but with all possible connectivities, including those below args.RESONANCE_MATCH_CUTOFF
        
            
            with open(args.POOL_CONNECTIVITIES_FILE, 'r') as f:
                pool_connectivities_file_contents = f.readlines()
            for line in pool_connectivities_file_contents[1:]:
                word_list = line.split()
                TOCSY_aaindex = word_list[0]
                i_iminus1_pool_dict[TOCSY_aaindex] = []
                values_string = ''.join(word_list[1:])
                elements_string = re.sub('[\(\)\'\"]', '',  values_string ).split("),(")[0]
                elements_list = elements_string.split(",")
                for aa, occupancy, TOCSY_resonnum in zip(elements_list[0::3], elements_list[1::3], elements_list[2::3]):
                    i_iminus1_pool_dict[TOCSY_aaindex].append((aa, int(occupancy), int(TOCSY_resonnum)))
            
            with open(args.COMPLETE_CONNECTIVITIES_FILE, 'r') as f:
                complete_connectivities_file_contents = f.readlines()
            for line in complete_connectivities_file_contents[1:]:
                word_list = line.split()
                TOCSY_aaindex = word_list[0]
                i_iminus1_complete_dict[TOCSY_aaindex] = []
                values_string = ''.join(word_list[1:])
                elements_string = re.sub('[\(\)\'\"]', '',  values_string ).split("),(")[0]
                elements_list = elements_string.split(",")
                for aa, occupancy, TOCSY_resonnum in zip(elements_list[0::3], elements_list[1::3], elements_list[2::3]):
                    i_iminus1_complete_dict[TOCSY_aaindex].append((aa, int(occupancy), int(TOCSY_resonnum)))
            
            if '-tolH' in sys.argv and '-tolC' in sys.argv:
                print "Calculating Connectivities using H tolerance", args.tolH, "and C torelance", args.tolC
                new_i_iminus1_dict, new_i_iminus1_complete_dict = create_connectivity_dictionaries()
                for TOCSYaaindex in new_i_iminus1_complete_dict.keys():
                    if TOCSYaaindex not in i_iminus1_complete_dict.keys():  # if there was not connectivity for this Tindex add it to the dictionaries               
                        i_iminus1_complete_dict[TOCSYaaindex] = new_i_iminus1_complete_dict[TOCSYaaindex]   # update the connectivities
                        i_iminus1_pool_dict[TOCSYaaindex] = new_i_iminus1_complete_dict[TOCSYaaindex]
                    elif len(new_i_iminus1_complete_dict[TOCSYaaindex]) > len(i_iminus1_complete_dict[TOCSYaaindex]): # if there were found extra connectivities using the new tolerances
                        if len(i_iminus1_pool_dict[TOCSYaaindex]) == 1 and len(i_iminus1_complete_dict[TOCSYaaindex]) > 1:
                            are_all_unique = True   # do all Tindices in the complete dict have a single connectivity in the pool
                            for tmp_triplet in i_iminus1_complete_dict[TOCSYaaindex]:
                                tmp_TOCSYaaindex = tmp_triplet[0]
                                if tmp_TOCSYaaindex in i_iminus1_pool_dict.keys() and len(i_iminus1_pool_dict[tmp_TOCSYaaindex]) > 1:
                                    are_all_unique = False
                                    break
                            if are_all_unique == False: # this means that i_iminus1_pool_dict[TOCSYaaindex] has been modified and not cleaned from Tindices used elsewhere
                                continue    # leave the connectivity of this TOCSYaaindex unchanged
                            i_iminus1_complete_dict[TOCSYaaindex] = new_i_iminus1_complete_dict[TOCSYaaindex]   # update the connectivities
                            i_iminus1_pool_dict[TOCSYaaindex] = new_i_iminus1_complete_dict[TOCSYaaindex]
                        elif len(i_iminus1_pool_dict[TOCSYaaindex]) == 1 and len(i_iminus1_complete_dict[TOCSYaaindex]) == 1:
                            i_iminus1_complete_dict[TOCSYaaindex] = new_i_iminus1_complete_dict[TOCSYaaindex]   # update the connectivities
                            i_iminus1_pool_dict[TOCSYaaindex] = new_i_iminus1_complete_dict[TOCSYaaindex]
            
            print "Finding Connectivities above the -mcutoff from the Pool of Connectivities..."
            i_iminus1_dict = {} # a dictionary containing all possible connectivities of every TOCSY aa index above args.RESONANCE_MATCH_CUTOFF; has keys the TOCSY aa indices and values lists of triplets (tuples) \
            for TOCSYaaindex in i_iminus1_pool_dict.keys():
                for triplet in i_iminus1_pool_dict[TOCSYaaindex]:
                    numOfResonances = triplet[2]
                    if float(triplet[1])/numOfResonances >= float(args.RESONANCE_MATCH_CUTOFF):     # adjustable cutoff
                        try:
                            i_iminus1_dict[TOCSYaaindex].append(triplet)
                        except KeyError:
                            i_iminus1_dict[TOCSYaaindex] = [(triplet)]
        
            new_i_iminus1_dict = {}
            for i in i_iminus1_dict.keys():
                prob_list = [float(triplet[1])/triplet[2] for triplet in i_iminus1_dict[i]]
                if len(prob_list) == 1:     
                    zscore_array = np.array([10.00000])     # default Z-score in the case of just one connectivity
                elif len(set(prob_list)) == 1:  # if all connectivities have the same value keep them all
                    zscore_array = np.array([10.00000]*len(prob_list))
                elif len(prob_list) > 2:
                    zscore_array = zscore(prob_list)  # convert probabilities to Z-scores
                elif len(prob_list) == 2:   # if only 2 connectivities exist, set the 2nd to 0.000 to allow you to use high -zcutoff values (<= 0)
                    zscore_array = np.array([1.00000, 0.00000])
                elif len(prob_list) == 0:   # if no connectivity was made create an empty array
                    zscore_array = np.array([])
                new_i_iminus1_dict[i] = []  # add i to the new dictionary
                for triplet, Zscore in zip(i_iminus1_dict[i], zscore_array):
                    if Zscore >= args.ZSCORE_RESONANCE_MATCH_CUTOFF or approx_equal(Zscore, args.ZSCORE_RESONANCE_MATCH_CUTOFF):
                        new_i_iminus1_dict[i].append(triplet)
            
            del i_iminus1_dict
            i_iminus1_dict = new_i_iminus1_dict
            del new_i_iminus1_dict
            
            
        i_iminus1_normProbabilities_dict = {}   # dict of the form Tindex(i)->[(Tindex(i-1), probability), (...), ...]; eventually will contain only the matches above the cutoff
        aaindex_weightSum_dict = {}
        for TOCSY_aaindex in i_iminus1_complete_dict.keys():
            weight_sum = 0
            for triplet in i_iminus1_complete_dict[TOCSY_aaindex]:  # use all possible matches to calculate the weight_sum
                occupancy = triplet[1]
                TOCSY_resonnum = triplet[2]
                weight = float(occupancy)/TOCSY_resonnum   # occupancy/total number of TOCSY resonances
                weight_sum += weight
            aaindex_weightSum_dict[TOCSY_aaindex] = weight_sum
        
        for TOCSY_aaindex in i_iminus1_dict.keys():     # use only the matches above the cutoff to save theirs probabilities
            for triplet in i_iminus1_dict[TOCSY_aaindex]:
                aa = triplet[0]
                occupancy = triplet[1]
                TOCSY_resonnum = triplet[2]
                weight = float(occupancy)/TOCSY_resonnum   # occupancy/total number of TOCSY resonances
                duplet = (aa, weight/aaindex_weightSum_dict[TOCSY_aaindex])           # recover the probability by dividing by the weight_sum
                try:
                    i_iminus1_normProbabilities_dict[TOCSY_aaindex].append(duplet)
                except KeyError:
                    i_iminus1_normProbabilities_dict[TOCSY_aaindex] = [duplet]
        
         
        
        
        if not args.NON_REDUNDANT_CHAINS_FILE:
            
            all_chainScore_set = set()    # a list of lists containing the possible connectivities and the overall score as the last element
            with open("connectivities_cutoff_"+str(args.RESONANCE_MATCH_CUTOFF)+"_zmcutoff"+str(args.ZSCORE_RESONANCE_MATCH_CUTOFF), 'w') as f:
                f.write("i\tpossible i-1\n")
                for k,v in i_iminus1_dict.items():
                    print k,v
                    f.write(k+"\t"+', '.join(str(e) for e in v)+"\n")
            
            with open("connectivities_all", 'w') as f:
                f.write("i\tpossible i-1\n")
                for k,v in i_iminus1_complete_dict.items():
                    f.write(k+"\t"+', '.join(str(e) for e in v)+"\n")
            
            
            if not args.SKIP_CHAINS:
                for i in i_iminus1_dict.keys():
                    build_Chain_Tree(i, all_chainScore_set, i_iminus1_dict, i_iminus1_normProbabilities_dict, args.MAX_PEPTIDE_LENGTH, args.MIN_PEPTIDE_LENGTH)
                
                
                print "Saving all the chains to chains.list file ..."
                def write_chains_file(all_chainScore_set, outfname):
                    """
                    ARGUMENTS:
                    all_chainScore_set: can be both a set or a list
                    """
                    all_chainScore_list = list(all_chainScore_set)
                    all_chainScore_list.sort(key=itemgetter(-1), reverse=True)
                    with open(outfname, 'w') as f:
                        for chainScore_list in all_chainScore_list:
                            chain = chainScore_list[0:-1]
                            score = chainScore_list[-1]
                            f.write("Overall Score = "+str(score)+" chain = "+str(chain)+"\n")
                
                write_chains_file(all_chainScore_set, "chains.list")
                
                
                print "Removing redundant chains ..."
                all_chainScore_list = list(all_chainScore_set)  # convert the set to a list to keep track of the sequence that each score belongs to
                del all_chainScore_set  # delete it to save memory
                progbar = ProgressBar(100)
                all_chainScore_listList_list = chunkIt(all_chainScore_list, scoop.SIZE)
                try:
                    shared.setConst(PROGBAR=progbar)
                except TypeError:
                    pass
                shared.setConst(ALL_CHAINSCORE_LIST=all_chainScore_list)
                results = list(futures.map(remove_redundant_chains, all_chainScore_listList_list))
                chainScores2remove_set = set()
                for s in results:
                    chainScores2remove_set.union(s)
                
                for chainScore in chainScores2remove_set:
                    all_chainScore_list.remove(tuple(chainScore))
                
                write_chains_file(all_chainScore_list, "chains.non-redundant.list")
        
        else:   # otherwise load the specified non-redudant chains file
            all_chainScore_list = []
            with open(args.NON_REDUNDANT_CHAINS_FILE, 'r') as f:
                for line in f:
                    word_list = line.split()
                    chain_score = float(word_list[3])
                    values_string = ''.join(word_list[6:])
                    Tindices_string = re.sub('[\(\)\'\"]', '',  values_string ).split("),(")[0]
                    chain = Tindices_string.split(",")
                    chain.append(chain_score)   # append the chain score to the end of the chain list
                    all_chainScore_list.append(chain)
        
        
        
        if not args.POOL_AA_TYPES_FILE and not args.COMPLETE_AA_TYPES_FILE:
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
            
            
            iaaindex_iminus1aaTypesProbTupleList_dict = OrderedDict() # ordereddict with keys the aa index of residue i and values lists of tuples of the form (residue i-1 matching aa type, average probability)
            Num_of_TOCSY_resonances = 0 # save here the number of TOCSY resonances for a particular aa index
            previous_TOCSY_aaindex = None
            possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list = []
            carbon_groups_list = []  # list of tuples of the form (aaindex, H resonance, C resonance, N resonance, HN resonance, carbon group)
            aaindex_carbonGroupsTuple_dict = {} # dict of the form: aaindex -> list of the form (aaindex, H resonance, C resonance, N resonance, HN resonance, carbon group)
            counter = 0
            previous_aaindex = ""
            sorted_TOCSY_contents = sorted(TOCSY_contents, key=itemgetter(0))
            print "sorted_TOCSY_contents=", sorted_TOCSY_contents
            for TOCSY_words_list in sorted_TOCSY_contents:
                try:
                    TOCSY_aaindex=TOCSY_words_list[0]    # residue i
                    print "DEBUG: TOCSY_aaindex=",TOCSY_aaindex,"previous_TOCSY_aaindex=",previous_TOCSY_aaindex
                    if previous_TOCSY_aaindex != None and TOCSY_aaindex != previous_TOCSY_aaindex: # if we look at a different aa index in TOCSY, print the matches and occupancies
                        
                        print "Assigning possible aa types to the aa upstream of aa index",previous_TOCSY_aaindex
                        new_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list = group_carbons(possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list)
                        iaaindex_iminus1aaTypesProbTupleList_dict[previous_TOCSY_aaindex] = get_aatypes_from_all_H_C_resonpairs(new_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list, Num_of_TOCSY_resonances, args)
                        possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list = [] # list of tuples of the form (aa type, probability, H_resonance, C_resonance, TOCSY_reson_index) containing all matching sets of H,C resonances
                        aaindex_carbonGroupsTuple_dict[TOCSY_aaindex] = carbon_groups_list
                        carbon_groups_list = []  # list of tuples of the form (TOCSY_aaindex, H resonance, C resonance, N resonance, HN, resonance, carbon group)
                        
                        Num_of_TOCSY_resonances = 0
                    
                    TOCSY_H_resonance=float(TOCSY_words_list[1]) # aliphatic H resonance of residue i-1 
                    TOCSY_C_resonance=float(TOCSY_words_list[2]) # aliphatic C (Ca,Cb,Cc,Cg,Ce,etc.) resonance of residue i-1; this C is covalently bonded to the above H
                    TOCSY_N_resonance=float(TOCSY_words_list[3])
                    TOCSY_HN_resonance=float(TOCSY_words_list[4])
                    Num_of_TOCSY_resonances += 1
                    if args.USE_2D_HISTOGRAMS == True:
                        valid_matches_list = get_aatypes_from_H_C_resonpair_2Dhist(TOCSY_H_resonance, TOCSY_C_resonance, Num_of_TOCSY_resonances, aa_CHpair_binProbabilityList_multidict, aa_carbon_binDensityList_multidict)
                    else:
                        valid_matches_list = get_aatypes_from_H_C_resonpair(TOCSY_H_resonance, TOCSY_C_resonance, Num_of_TOCSY_resonances, aa_carbon_binDensityList_multidict, aa_hydrogen_binDensityList_multidict, args)
                    possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list.extend(valid_matches_list)
                    carbon_groups_list.append([TOCSY_aaindex, TOCSY_H_resonance, TOCSY_C_resonance, TOCSY_N_resonance, TOCSY_HN_resonance, None])
                
                    previous_TOCSY_aaindex = TOCSY_aaindex
                except (ValueError, IndexError):
                    print "WARNING: the 3rd and 4th elements of the following TOCSY file line are not numbers:"
                    print "TOCSY file line:", TOCSY_words_list
                    continue
            print "Assigning possible aa types to the last aa, which is upstream of aa index",previous_TOCSY_aaindex
            new_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list = group_carbons(possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list)
            iaaindex_iminus1aaTypesProbTupleList_dict[previous_TOCSY_aaindex] = get_aatypes_from_all_H_C_resonpairs(new_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list, Num_of_TOCSY_resonances, args)
            
            iaaindex_iminus1aaTypesZscoreTupleList_dict = OrderedDict() # ordereddict with keys the aa index of residue i and values lists of tuples of the form (residue i-1 matching aa type, Z-score)
            iaaindex_iminus1aaTypesCutoffProbTupleList_dict = OrderedDict() # contains only the tuples (residue i-1 matching aa type, probability) with Z-score above the cutoff
            for iaaindex, iminus1aaTypesProbTuple_list in iaaindex_iminus1aaTypesProbTupleList_dict.items():
                try:
                    max_prob = np.max([duplet[1] for duplet in iminus1aaTypesProbTuple_list])   # highest aa type probability
                except ValueError:  # if no aa type predictions exist for this Tindex group
                    max_prob = 0
                if max_prob > 10e-10:   # Dedice the probability threshold
                    prob_threshold = 1000
                elif max_prob > 10e-20:
                    prob_threshold = 10000
                elif max_prob <= 10e-20:
                    prob_threshold = 100000
                print "DEBUG: max_prob", max_prob, "prob_threshold", prob_threshold, "args.LOG=", args.LOG
                aatype_list, prob_list = [], []
                for duplet in iminus1aaTypesProbTuple_list:
                    if args.DELETE_AA_TYPE_PREDICTIONS == True:     # delete low probability aa type predictions
                        try:
                            if max_prob/float(duplet[1]) < prob_threshold:
                                aatype_list.append(duplet[0])
                                prob_list.append(duplet[1])  # unsorted list of probabilities
                        except ZeroDivisionError:
                            print "DEBUG: ZeroDivisionError ", max_prob, "/", float(duplet[1])
                    elif args.DELETE_AA_TYPE_PREDICTIONS == False:
                        aatype_list.append(duplet[0])
                        prob_list.append(duplet[1])  # unsorted list of probabilities
                print "DEBUG: prob_list=", prob_list
                if len(prob_list) == 1:     
                    zscore_array = np.array([10.00000])     # default Z-score in the case of just one prediction
                elif len(prob_list) > 2:
                    if args.LOG == True:
                        try:
                            ratio = float(np.max(prob_list))/np.min(prob_list)    # if the min probability is at least 3 orders of magnitude smaller, convert them to logarithmic scale
                        except ZeroDivisionError:
                            print "DEBUG: ZeroDivisionError ", np.max(prob_list), "/", np.min(prob_list) 
                            ratio = -1.0
                        if ratio > 1000:
                            zscore_array = zscore(np.log(prob_list))  # convert probabilities the logarithms and then to Z-scores
                        else:
                            zscore_array = zscore(prob_list)  # convert probabilities to Z-scores
                    elif args.LOG == False:
                        zscore_array = zscore(prob_list)  # convert probabilities to Z-scores
                elif len(prob_list) == 2:   # if only 2 predictions exist, set the 2nd to 0.000 to allow you to use high -zcutoff values (<= 0)
                    zscore_array = np.array([1.00000, 0.00000]) 
                else:   # if no aa type prediction was made create an empty array
                    zscore_array = np.array([])
                print "DEBUG2: saving zscore_array = ", zscore_array
                iminus1aaTypesZscoreTuple_list = []
                iminus1aaTypesCutoffProbTuple_list = []
                for aatype, Zscore, prob in zip(aatype_list, zscore_array, prob_list):
                    if Zscore > args.ZSCORE_ASSIGNMENT_CUTOFF and len(zscore_array) >= args.MIN_NUM_OF_PREDICTIONS:
                        aatypeProb_tuple = (aatype, prob)
                        aatypeZscore_tuple = (aatype, Zscore)  # replace the probability with the respective Z-score
                        iminus1aaTypesZscoreTuple_list.append(aatypeZscore_tuple)
                        iminus1aaTypesCutoffProbTuple_list.append(aatypeProb_tuple) # but also save the aatype probabilities that satisfy the Z-score cutoff
                    elif len(zscore_array) < args.MIN_NUM_OF_PREDICTIONS:
                        aatypeProb_tuple = (aatype, prob)
                        aatypeZscore_tuple = (aatype, Zscore)  # replace the probability with the respective Z-score
                        iminus1aaTypesZscoreTuple_list.append(aatypeZscore_tuple)
                        iminus1aaTypesCutoffProbTuple_list.append(aatypeProb_tuple) # but also save the aatype probabilities that satisfy the Z-score cutoff
                iaaindex_iminus1aaTypesZscoreTupleList_dict[iaaindex] = iminus1aaTypesZscoreTuple_list
                iaaindex_iminus1aaTypesCutoffProbTupleList_dict[iaaindex] = iminus1aaTypesCutoffProbTuple_list
            
            
            with open("amino_acid_type_prediction_probabilities", 'w') as f:
                f.write("i aa index\tpossible i-1 aa types\n")
                print "Probabilities:"
                for i_aaindex in iaaindex_iminus1aaTypesProbTupleList_dict.keys():
                    print i_aaindex,"---> (i-1) aa type",sorted(iaaindex_iminus1aaTypesProbTupleList_dict[i_aaindex], key=itemgetter(1), reverse=True)
                    f.write(i_aaindex + " ---> (i-1) aa type " + str(sorted(iaaindex_iminus1aaTypesProbTupleList_dict[i_aaindex], key=itemgetter(1), reverse=True)) + "\n")
            
            with open("amino_acid_type_prediction_probabilities_above_cutoff.zacutoff"+str(args.ZSCORE_ASSIGNMENT_CUTOFF), 'w') as f:
                f.write("i aa index\tpossible i-1 aa types\n")
                print "Probabilities:"
                for i_aaindex in iaaindex_iminus1aaTypesCutoffProbTupleList_dict.keys():
                    print i_aaindex,"---> (i-1) aa type",sorted(iaaindex_iminus1aaTypesCutoffProbTupleList_dict[i_aaindex], key=itemgetter(1), reverse=True)
                    f.write(i_aaindex + " ---> (i-1) aa type " + str(sorted(iaaindex_iminus1aaTypesCutoffProbTupleList_dict[i_aaindex], key=itemgetter(1), reverse=True)) + "\n")
            
            with open("amino_acid_type_prediction_Z-scores.zacutoff"+str(args.ZSCORE_ASSIGNMENT_CUTOFF), 'w') as f:
                f.write("i aa index\tpossible i-1 aa types\n")
                print "Z-scores:"
                for i_aaindex in iaaindex_iminus1aaTypesZscoreTupleList_dict.keys():
                    print i_aaindex,"---> (i-1) aa type",sorted(iaaindex_iminus1aaTypesZscoreTupleList_dict[i_aaindex], key=itemgetter(1), reverse=True)
                    f.write(i_aaindex + " ---> (i-1) aa type " + str(sorted(iaaindex_iminus1aaTypesZscoreTupleList_dict[i_aaindex], key=itemgetter(1), reverse=True)) + "\n")
        else:
            print ""    # change line after the previous progress bar
            print "Loading amino acid prediction files ..."
            minimum_Zscore = 10000000
            iaaindex_iminus1aaTypesProbPoolTupleList_dict = OrderedDict() # ordereddict with keys the aa index of residue i and values lists of tuples of the form
            with open(args.POOL_AA_TYPES_FILE, 'r') as f:
                pool_aa_type_file_contents = f.readlines()
                for line in pool_aa_type_file_contents[1:]:
                    word_list = re.sub('---> \(i-1\) aa type', '', line).split()
                    key = word_list[0]
                    iaaindex_iminus1aaTypesProbPoolTupleList_dict[key] = []
                    values_string = ''.join(word_list[1:])
                    elements_string = re.sub('[\(\)\[\]\'\"]', '',  values_string ).split("),(")[0]
                    elements_list = elements_string.split(",")
                    for aa, Zscore in zip(elements_list[0::2], elements_list[1::2]):
                        duplet = (aa, float(Zscore))
                        iaaindex_iminus1aaTypesProbPoolTupleList_dict[key].append(duplet)
            
            iaaindex_iminus1aaTypesZscoreTupleList_dict = OrderedDict() # ordereddict with keys the aa index of residue i and values lists of tuples of the form (residue i-1 matching aa type, Z-score)
            iaaindex_iminus1aaTypesCutoffProbTupleList_dict = OrderedDict() # contains only the tuples (residue i-1 matching aa type, probability) with Z-score above the cutoff
            for iaaindex, iminus1aaTypesProbTuple_list in iaaindex_iminus1aaTypesProbPoolTupleList_dict.items():
                try:
                    max_prob = np.max([duplet[1] for duplet in iminus1aaTypesProbTuple_list])   # highest aa type probability
                except ValueError:  # if no aa type predictions exist for this Tindex group
                    max_prob = 0
                if max_prob > 10e-10:   # Dedice the probability threshold
                    prob_threshold = 1000
                elif max_prob > 10e-20:
                    prob_threshold = 10000
                elif max_prob <= 10e-20:
                    prob_threshold = 100000
                print "DEBUG: max_prob", max_prob, "prob_threshold", prob_threshold, "args.LOG=", args.LOG
                aatype_list, prob_list = [], []
                for duplet in iminus1aaTypesProbTuple_list:
                    if args.DELETE_AA_TYPE_PREDICTIONS == True:     # delete low probability aa type predictions
                        try:
                            if max_prob/float(duplet[1]) < prob_threshold:
                                aatype_list.append(duplet[0])
                                prob_list.append(duplet[1])  # unsorted list of probabilities
                        except ZeroDivisionError:
                            print "DEBUG: ZeroDivisionError ", max_prob, "/", float(duplet[1])
                    elif args.DELETE_AA_TYPE_PREDICTIONS == False:
                        aatype_list.append(duplet[0])
                        prob_list.append(duplet[1])  # unsorted list of probabilities
                print "DEBUG: prob_list=", prob_list
                if len(prob_list) == 1:     
                    zscore_array = np.array([10.00000])     # default Z-score in the case of just one prediction
                elif len(prob_list) > 1:
                    if args.LOG == True:
                        try:
                            ratio = float(np.max(prob_list))/np.min(prob_list)    # if the min probability is at least 3 orders of magnitude smaller, convert them to logarithmic scale
                        except ZeroDivisionError:
                            print "DEBUG: ZeroDivisionError ", np.max(prob_list), "/", np.min(prob_list) 
                            ratio = -1.0
                        if ratio > 1000:
                            zscore_array = zscore(np.log(prob_list))  # convert probabilities the logarithms and then to Z-scores
                        else:
                            zscore_array = zscore(prob_list)  # convert probabilities to Z-scores
                    elif args.LOG == False:
                        zscore_array = zscore(prob_list)  # convert probabilities to Z-scores
                else:   # if no aa type prediction was made create an empty array
                    zscore_array = np.array([])
                iminus1aaTypesZscoreTuple_list = []
                iminus1aaTypesCutoffProbTuple_list = []
                for aatype, Zscore, prob in zip(aatype_list, zscore_array, prob_list):
                    if Zscore < minimum_Zscore:
                        minimum_Zscore = Zscore
                    if Zscore > args.ZSCORE_ASSIGNMENT_CUTOFF and len(zscore_array) >= args.MIN_NUM_OF_PREDICTIONS:
                        aatypeProb_tuple = (aatype, prob)
                        aatypeZscore_tuple = (aatype, Zscore)  # replace the probability with the respective Z-score
                        iminus1aaTypesZscoreTuple_list.append(aatypeZscore_tuple)
                        iminus1aaTypesCutoffProbTuple_list.append(aatypeProb_tuple) # but also save the aatype probabilities that satisfy the Z-score cutoff
                    elif len(zscore_array) < args.MIN_NUM_OF_PREDICTIONS:
                        aatypeProb_tuple = (aatype, prob)
                        aatypeZscore_tuple = (aatype, Zscore)  # replace the probability with the respective Z-score
                        iminus1aaTypesZscoreTuple_list.append(aatypeZscore_tuple)
                        iminus1aaTypesCutoffProbTuple_list.append(aatypeProb_tuple) # but also save the aatype probabilities that satisfy the Z-score cutoff
                iaaindex_iminus1aaTypesZscoreTupleList_dict[iaaindex] = iminus1aaTypesZscoreTuple_list
                iaaindex_iminus1aaTypesCutoffProbTupleList_dict[iaaindex] = iminus1aaTypesCutoffProbTuple_list
            
            
            iaaindex_iminus1aaTypesProbTupleList_dict = OrderedDict()   # ordereddict with keys the aa index of residue i and values lists of tuples of the form
            with open(args.COMPLETE_AA_TYPES_FILE, 'r') as f:
                complete_aa_type_file_contents = f.readlines()
                for line in complete_aa_type_file_contents[1:]:
                    word_list = re.sub('---> \(i-1\) aa type', '', line).split()
                    key = word_list[0]
                    iaaindex_iminus1aaTypesProbTupleList_dict[key] = []
                    values_string = ''.join(word_list[1:])
                    elements_string = re.sub('[\(\)\[\]\'\"]', '',  values_string ).split("),(")[0]
                    elements_list = elements_string.split(",")
                    for aa, hist_prob in zip(elements_list[0::2], elements_list[1::2]):
                        duplet = (aa, float(hist_prob))
                        iaaindex_iminus1aaTypesProbTupleList_dict[key].append(duplet)
            
            
            with open("amino_acid_type_prediction_Z-scores.zacutoff"+str(args.ZSCORE_ASSIGNMENT_CUTOFF), 'w') as f:
                f.write("i aa index\tpossible i-1 aa types\n")
                print "Z-scores:"
                for i_aaindex in iaaindex_iminus1aaTypesZscoreTupleList_dict.keys():
                    print i_aaindex,"---> (i-1) aa type",sorted(iaaindex_iminus1aaTypesZscoreTupleList_dict[i_aaindex], key=itemgetter(1), reverse=True)
                    f.write(i_aaindex + " ---> (i-1) aa type " + str(sorted(iaaindex_iminus1aaTypesZscoreTupleList_dict[i_aaindex], key=itemgetter(1), reverse=True)) + "\n")
        
        
        template_sequence_string = ""
        with open(args.template_sequence_file, 'r') as f:
            for line in f:
                if not re.match("^>", line):
                    template_sequence_string += re.sub(r"[^A-Z]", "", line)
        template_sequence_list = list(template_sequence_string)
        
        aatype_P_dict = {}    # dict with the P[aatype(i-1)] for every aatype(i-1): the percentage of each aatype in the template sequence
        for aa_type in aatype_maxH_C_pairs_dict.keys():
            aatype_P_dict[aa_type] = template_sequence_list.count(aa3to1_dict[aa_type]) / float(len(template_sequence_list))
        
        all_prob_list = []
        for iaaindex, iminus1aaTypesProbTuple_list in iaaindex_iminus1aaTypesProbTupleList_dict.items():
            for duplet in iminus1aaTypesProbTuple_list:
                prob = duplet[1]
                all_prob_list.append(prob)
        
        all_prob_array = np.array(all_prob_list)
        transformed_all_prob_array = transform(all_prob_array)
        
        iaaindex_iminus1aaTypesTransformedProbTupleList_dict = OrderedDict()
        array_index = 0
        for iaaindex, iminus1aaTypesProbTuple_list in iaaindex_iminus1aaTypesProbTupleList_dict.items():
            iaaindex_iminus1aaTypesTransformedProbTupleList_dict[iaaindex] = []
            for duplet in iminus1aaTypesProbTuple_list:
                aatype = duplet[0]
                trans_prob = transformed_all_prob_array[array_index]
                iaaindex_iminus1aaTypesTransformedProbTupleList_dict[iaaindex].append((aatype, trans_prob))
                array_index += 1
        
        def calculate_Paaiminus1_Tindexi(Tindexi, aaiminus1):
            """
                Calculate the conditional probability P(AA^{u} | AAIG^{CS}), which quantifies how probable is to get an amino-acid of type u given a set of
                chemical shifts AAIG^{CS}.
                ARGS:
                Tindexi:    the AAIG or TOCSY Index Group aligned to position i
                aaiminus1:  the aa type if the residue in position i-1
            """
            global aatype_P_dict, iaaindex_iminus1aaTypesTransformedProbTupleList_dict
            
            PTindexi = 0    # P[Tindex(i)] = Sum{P[Tindex(i)|aatype(i-1)]}
            for duplet in iaaindex_iminus1aaTypesTransformedProbTupleList_dict[Tindexi]:
                aatype = duplet[0]
                PTindexi += duplet[1]
                if aatype == aaiminus1:
                    PTindexi_aaiminus1 = duplet[1]  # this is P[Tindex(i)|aatype(i-1)]
            
            Paaiminus1 = aatype_P_dict[aaiminus1]
            Paaiminus1_Tindexi = PTindexi_aaiminus1 * Paaiminus1 / PTindexi
            return Paaiminus1_Tindexi
        
        
        print "DEBUG: iaaindex_iminus1aaTypesProbTupleList_dict=", iaaindex_iminus1aaTypesProbTupleList_dict
        iaaindex_iminus1aaType_Condprob_multidict = tree() # multidict of the form aa index of residue i --> residue i-1 matching aa type --> conditional probability
        for i_aaindex in iaaindex_iminus1aaTypesProbTupleList_dict.keys():
            aatypeCondprob_list = []    # list of the form [(aa type, conditional probability), (aa type, conditional probability), ...]
            for aa_type,prob in sorted(iaaindex_iminus1aaTypesProbTupleList_dict[i_aaindex], key=itemgetter(1), reverse=True):
                condprob = calculate_Paaiminus1_Tindexi(i_aaindex, aa_type)
                iaaindex_iminus1aaType_Condprob_multidict[i_aaindex][aa_type] = condprob
        
            with open("amino_acid_type_prediction_conditional_probabilities", 'w') as f:
                f.write("i aa index\tpossible i-1 aa types\n")
                print "Conditional Probabilities:"
                for i_aaindex in iaaindex_iminus1aaType_Condprob_multidict.keys():
                    aa_type_condprob_list = []  # list of the form [(aa type, conditional probability), (aa type, conditional probability), ...]
                    for aa_type in iaaindex_iminus1aaType_Condprob_multidict[i_aaindex].keys():
                        aa_type_condprob_list.append( (aa_type, iaaindex_iminus1aaType_Condprob_multidict[i_aaindex][aa_type]) )
                    print i_aaindex,"---> (i-1) aa type",sorted(aa_type_condprob_list, key=itemgetter(1), reverse=True)
                    f.write(i_aaindex + " ---> (i-1) aa type " + str(sorted(aa_type_condprob_list, key=itemgetter(1), reverse=True)) + "\n")
            
        
        
        if not args.SKIP_CHAINS:
            if os.path.exists('tmp_peptide_folder/') == True and args.RESUME == False:   # if the folder exists and this is not a resumed run, remove it and make a new one
                shutil.rmtree('tmp_peptide_folder/', ignore_errors=True)
                os.makedirs('tmp_peptide_folder/')
            elif os.path.exists('tmp_peptide_folder/') == False and args.RESUME == False: # if the folder does not exists, make a new one
                os.makedirs('tmp_peptide_folder/')
            elif os.path.exists('tmp_peptide_folder/') == False and args.RESUME == True:
                print "ERROR: folder tmp_peptide_folder/ with peptide sequence files does not exist! The process cannot be resumed!"
            fnames=os.listdir('tmp_peptide_folder/')
            fpattern = re.compile('chainIndex_[0-9]+.pickle.bz2')
            peptide_files_list = filter(fpattern.search, fnames)
            chainIndices2skip_set = set()
            for peptide_file in peptide_files_list:
                mo = re.search('chainIndex_([0-9]+).pickle.bz2', peptide_file)
                if mo:
                    chainIndex = int(mo.group(1))
                    chainIndices2skip_set.add(chainIndex)
            
            total_chain_number = len(all_chainScore_list)
            remaining_chainIndex_list = []
            remaining_chainScore_list = []
            shared.setConst(TOTAL_CHAIN_NUMBER=total_chain_number)
            shared.setConst(ARGS=args)
            shared.setConst(IAAINDEX_IMINUS1AATYPE_CONDPROB_MULTIDICT = iaaindex_iminus1aaType_Condprob_multidict)
            shared.setConst(IAAINDEX_IMINUS1AATYPESCUTOFFPROBTUPLELIST_DICT=iaaindex_iminus1aaTypesCutoffProbTupleList_dict)
            shared.setConst(AA3TO1_DICT=aa3to1_dict)
            shared.setConst(IAAINDEX_IMINUS1AATYPESTRANSFORMEDPROBTUPLELIST_DICT=iaaindex_iminus1aaTypesTransformedProbTupleList_dict)
            shared.setConst(AATYPE_P_DICT=aatype_P_dict)
            shared.setConst(I_IMINUS1_DICT=i_iminus1_dict)
            for chainIndex, chainScore in enumerate(all_chainScore_list):
                if not chainIndex in chainIndices2skip_set:
                    remaining_chainIndex_list.append(chainIndex)
                    remaining_chainScore_list.append(chainScore)
            
            results = list(futures.map(build_Peptide_Tree, remaining_chainIndex_list, remaining_chainScore_list))   # build peptide tree from chain
            
            progbar = ProgressBar(100)
            peptide_file_num = 0
            fnames=os.listdir('tmp_peptide_folder/')
            fpattern = re.compile('chainIndex_[0-9]+.pickle.bz2')
            peptide_files_list = filter(fpattern.search, fnames)
            peptide_file_num = len(peptide_files_list)
            for chainIndex, chainScore in enumerate(all_chainScore_list):
                if not 'chainIndex_'+str(chainIndex)+'.pickle.bz2' in peptide_files_list:
                    build_Peptide_Tree(chainIndex, chainScore)
            
            print "Saving all the peptide sequences to peptides.list file ..."
            peptide_filesList_list = chunkIt(peptide_files_list, scoop.SIZE)
            parts_list = [str(part) for part in range(1, scoop.SIZE +1)]
            try:
                shared.setConst(ARGS=args)
            except TypeError:
                pass
            try:
                shared.setConst(PROGBAR=progbar)
            except TypeError:
                pass
            shared.setConst(NEW_ALL_CHAINSCORE_LIST=all_chainScore_list)
            results = list(futures.map(save_peptides_to_file, peptide_filesList_list, parts_list))
            fasta_filehandler = open("peptides."+str(args.MAX_PEPTIDE_LENGTH)+"mers.fasta", 'w')
            f = open("peptides."+str(args.MAX_PEPTIDE_LENGTH)+"mers.list", 'w')
            index = 0
            for part in parts_list:
                with open("peptides."+str(args.MAX_PEPTIDE_LENGTH)+"mers.list.chunk"+part, 'r') as fin:
                    for line in fin:
                        mo = re.search("^(.* peptide )P_[0-9]+_[0-9]+ = ([A-Z]+)( .*)$", line)
                        if mo:
                            line_part1 = mo.group(1)
                            peptide_seq = mo.group(2)
                            line_part2 = mo.group(3)
                            peptide_name = "P_" + str(len(peptide_seq)) + "_" + str(index)
                            f.write(line_part1 + peptide_name + " = " + peptide_seq + line_part2 + "\n")
                            fasta_filehandler.write(">" + peptide_name + "\n" + peptide_seq + "\n")
                            index += 1
                os.remove("peptides."+str(args.MAX_PEPTIDE_LENGTH)+"mers.list.chunk"+part)
            shutil.rmtree('tmp_peptide_folder/', ignore_errors=True)
    
    except:
        type, value, tb = sys.exc_info()
        lines = traceback.format_exception(type, value, tb)
        print(''.join(lines))
        raise