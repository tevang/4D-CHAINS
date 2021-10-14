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

# Import 4D-CHAINS libraries
CHAINS_BIN_DIR = os.path.dirname(os.path.realpath(__file__))
from lib.spectrum_processing import *
from lib.aa_prediction import *
from lib.trees.peptide import *
from lib.trees.chain import *


## Set global variables
HOME_DIR = os.path.dirname(os.path.realpath(__file__))


## Parse command line arguments
def cmdlineparse():
    parser = ArgumentParser(description="command line arguments", formatter_class=RawDescriptionHelpFormatter,
                            epilog="EXAMPLE: 4D_assignment_tree-based.py -hsqc TDhsqc.list -tocsy TDtocsy.list -noesy TDnoesy.list -tolH 0.01 -tolC 0.1 -rtolH 0.01 -rtolN 0.1 -mcutoff 1.0 -acutoff 1.0 -wH 0.5 -wC 1")
    parser.add_argument("-fasta", dest="template_sequence_file", required=True,
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
    parser.add_argument("-moffset", dest="MULTI_OFFSET", required=False, type=float, default=0.15,
                        help="number 0.0-1.0 specifying the distance from the maximum multi score a possible connectivity must have"
                             "in order to be used.")
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
    parser.add_argument("-log", dest="LOG", required=False, action='store_true', default=False,
                        help="convert aa type prediction probabilities to logarithmic scale and then calculate Z-scores, if the \
                             min(probability)/max(probability) > 1000. (default %(default)s)")
    parser.add_argument("-skipchains", dest="SKIP_CHAINS", required=False, action='store_true', default=False,
                        help="Generate connectivies and amino acid type predictions, but skip chain and peptide formation. (default %(default)s)")
    parser.add_argument("-delpred", dest="DELETE_AA_TYPE_PREDICTIONS", required=False, action='store_true', default=False,
                        help="delete aa type predictions with probabilities which are 1000, 10000 or 100000 times lower than the highest, if \
                             the highest is >10e-10, >10e-20, <=10e-20, respectively. (default %(default)s)")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s version 1.2')
    args=parser.parse_args()
    return args


##
##  CHECK INPUT FILES FORMAT FOR CORRECTNESS
##

def check_input_files_format():
    global args
    
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
                # print "WARNING: Discarding", spectrum_combo, "line:", line
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
            if Chain.is_sublist(chain1, chain2) and score1 == score2: # this is the correct one
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
            #print "DEBUG: loading tmp_peptide_folder/" + peptide_file
            try:
                with bz2.BZ2File('tmp_peptide_folder/' + peptide_file, "rb") as pickled_file:
                    pickle_object = pickle.load(pickled_file)
            except EOFError:
                print("ERROR: file tmp_peptide_folder/" + peptide_file, " seems to be corrupted perhaps due to abrupt termination of the program before the file was properly saved.\
                Please run again the program without the -resume option.")
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
                    print("ERROR: the chainIndex of file ", peptide_file, " should be ", chainIndex2, ", not ", chainIndex)
                    sys.exit(1)
                if peptide_score == 0.0:    # do not save peptides containing aa types that are not present in the protein
                    continue
                #print "DEBUG: chain=", chain, " peptide = ", peptide, "index=", index, "peptideProbList=", peptideProbList
                f.write("Chain probability = "+str(chain[-1])+" peptide P_"+str(len(peptide))+"_"+str(index+1) + " = "+''.join(peptide)+" "+'-'.join(chain[0:-1])+" "+','.join(map(str, peptideProbList[:-1]))+"\n")
                index += 1
    f.close()


def split_oversized_Cgroups(possible_aatype_prob_C_H_resonpair_TAAIG_list_list):
    """
        FUNCTION to split Cgroups with more that 2 peaks into separate Cgroups. Currently works for maximum Cgroup size 3 only. This happens when
        two peaks have the same Carbon resonance but different proton resonance. So because C-grouping works with Carbon resonances only, all 3 peaks
        are included in the same group
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
        This is a different function for the group_carbons() used in the CS assignment scripts, because there we know the aa type, but here we
        group the carbons for every possible aa type.
        ARGUMENTS:
        possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list:    list of lists of the form
                (aa type, Carbon name, Hydrogen name, weighted average probability, TOCSY_reson_index, H resonance, C resonance, carbon group)
        RETURNS:
        The same as the input but with the correct carbon groups (last element)
    """
    global aa_carbonBondedGeminalHydrogensDict_dict
    
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
        
        new_aaindex_carbonGroupsTuple_dict = {}
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
    updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list = split_oversized_Cgroups(updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list)
    
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



if __name__ == "__main__":
    
    try:
        
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
        check_input_files_format()
        
        #if args.EXAMPLES:
        #    print "USAGE: python HCTOCSYNH_to_HCNOENH_tree-based.py -tocsy <HCTOCSYNH file> -noesy <HCNOENH file> [-stdH <stdev1>] [-stdC <stdev2>] [-update]"
        #    print "Example: ./4D_assignment_tree-based.py -hsqc TDhsqc.list -tocsy TDtocsy.list -noesy TDnoesy.list -stdH 0.03 -stdC 0.3 -rstdH 0.02 -rstdC 0.2 -cutoff 1.0"
        #    sys.exit(0)
        if args.DOWNLOAD_CS_HISTOGRAMS:
            download_CS_histograms()
        
        
        # GRADUALLY REDUCE TOLERANCES TO FIND MORE NON-OVERLAPPING GROUPS IN HSQC SPECTRUM
        with open(args.HSQC_FILE, 'r') as f:
            root_contents = []
            root_CS_set = set() # to remove replicate lines
            for line in f:
                try:
                    words = line.split()
                    CS_values = (float(words[1]), float(words[2]))
                    if CS_values in root_CS_set:
                        print(bcolors.WARNING + "WARNING: discarding replicate HSQC line: "+line + bcolors.ENDC)
                        continue
                    root_contents.append(line)
                    root_CS_set.add(CS_values)
                except (IndexError, ValueError):
                    print(bcolors.WARNING + "WARNING: Discarding invalid HSQC line: " + line + bcolors.ENDC)
        remaining_root_contents = []    # list of Root spectrum lines with overlapping groups
        nonoverlapping_root_contents_and_tolerances = []    # list wiht the lines of the Root spectrum that contain non-overlapping groups and the associated rtoH, rtolN
        for rtolH, rtolN in zip([0.02, 0.01, 0.01, 0.005], [0.2, 0.1, 0.05, 0.015]):
        #for rtolH, rtolN in zip([0.02, 0.01, 0.01], [0.2, 0.1, 0.05]):
            #print "DEBUG: len(root_contents)=", len(root_contents)
            # FIND NON-OVERLAPPING RESONANCE GROUPS IN THE HSQC SPECTRUM
            nonoverlapping_root_contents_and_tolerances.extend(find_nonoverlapping_groups_in_root(root_contents, rtolH, rtolN))
            for triplet in nonoverlapping_root_contents_and_tolerances:
                line = triplet[0]
                try:
                    root_contents.remove(line)
                except ValueError:
                    continue
        
        #sys.exit(1)
        ## AT THE END GROUP WHICHEVER PEAKS IS LEFT TO THE GROUP IN THE HSQC SPECTRUM THAT HAS THE NEAREST H,N RESONANCE DISTANCE
        #nonoverlapping_root_contents_and_tolerances.extend(find_nonoverlapping_groups_in_root(root_contents, rtolH, rtolN, True))
        #for triplet in nonoverlapping_root_contents_and_tolerances:
        #    line = triplet[0]
        #    try:
        #        root_contents.remove(line)
        #    except ValueError:
        #        continue
        
        #print "DEBUG: len(nonoverlapping_root_contents_and_tolerances)=", len(nonoverlapping_root_contents_and_tolerances)
        #for triplet in nonoverlapping_root_contents_and_tolerances:
        #    print triplet
        
        
        ##
        ## COPY AAIGs FROM HSQC SPECTRUM
        ##
        ## TOCSY_lines contains the AAIG of i, the H,C resonances of residue i-1 and the N,HN resonances of residue i
        print("DEBUG: final nonoverlapping_root_contents_and_tolerances=", nonoverlapping_root_contents_and_tolerances)
        TOCSY_lines, TOCSY_sidechain_resonances_list = copy_aaindices_from_root_spectrum_2(nonoverlapping_root_contents_and_tolerances, args.TOCSY_fname, "TOCSY")
        NOESY_lines = copy_aaindices_from_root_spectrum_2(nonoverlapping_root_contents_and_tolerances, args.NOESY_FILE, "NOESY", TOCSY_sidechain_resonances_list)
        
        ## Load TOCSY file
        TOCSY_contents = []
        print("DEBUG: TOCSY_lines=", TOCSY_lines)
        for TOCSY_line in TOCSY_lines:
            TOCSY_words_list = [a for a in re.split('\s+', TOCSY_line)[:-1] if a!= ''] # discard the last element which is ''
            try:
                [float(w) for w in TOCSY_words_list[1:]];   # works with any number of numbers (intensities or not)
                TOCSY_contents.append(TOCSY_words_list)
            except ValueError:
                print("Discarding the following line from TOCSY:")
                TOCSY_line
                continue
        TOCSY_contents.sort(key=itemgetter(0))  # sort TOCSY lines by the random AAIG (2nd column)
        
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
        
        def create_connectivity_dictionaries(clean_native_peaks=False):
            connectivities_mdict = get_possible_connectivities(TOCSY_contents, NOESY_lines, args.tolH, args.tolC,
                                                                   clean_native_peaks=clean_native_peaks)
            # connectivities_mdict is a multidimensional dictionary with structure:
            # TOCSY AAIG => NOESY AAIG => [number of matched TOCSY w2,w3 resonances in NOESY for that NOESY AAIG, total number of TOCSY w2,w3 resonances for that TOCSY AAIG]
            i_iminus1_dict = {} # a dictionary containing all possible connectivities of every TOCSY AAIG; has keys the TOCSY aa indices and values lists of triplets (tuples) \
                                # consisting of (NOESYaaindex, occupancy, numOfResonances)
            i_iminus1_complete_dict = {} # the same dictionary but with all possible connectivities, including those below args.RESONANCE_MATCH_CUTOFF
            for TOCSYaaindex in list(connectivities_mdict.keys()):
                NOESYaaindex_list = list(connectivities_mdict[TOCSYaaindex].keys())
                if TOCSYaaindex in NOESYaaindex_list:
                    NOESYaaindex_list.remove(TOCSYaaindex)  # remove residue i from the list of possible residues (i-1)
                # Put the NOESYaaindex, occupancy and total number of resonances for that particular NOESYaaindex into a tuple for convenience
                NOESYaaindex_occupancy_numOfResonances_tuple__list = [] 
                for NOESYaaindex in NOESYaaindex_list:
                    occupancy = connectivities_mdict[TOCSYaaindex][NOESYaaindex][0]
                    numOfResonances = connectivities_mdict[TOCSYaaindex][NOESYaaindex][1]
                    NOESYaaindex_occupancy_numOfResonances_tuple__list.append((NOESYaaindex, occupancy, numOfResonances))
                sorted_NOESYaaindex_occupancy_numOfResonances_tuple__list = sorted(NOESYaaindex_occupancy_numOfResonances_tuple__list,
                                                                                   key=itemgetter(1, 2), reverse=True)
                #print "DEBUG:",numOfResonances, TOCSYaaindex, sorted_NOESYaaindex_occupancy_numOfResonances_tuple__list
                # Iterate over all NOESY aa indices and, if applicable, keep those that match in all resonance pairs w2,w3
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
            
            # FILTER i_iminus1_dict ACCORDING TO THE Z-SCORE CUTOFF
            new_i_iminus1_dict = {}
            for i in list(i_iminus1_dict.keys()):
                prob_list = [float(triplet[1])/triplet[2] for triplet in i_iminus1_dict[i]]
                if len(prob_list) == 1:     
                    zscore_array = np.array([10.00000])     # default Z-score in the case of just one connectivity
                elif len(set(prob_list)) == 1:  # if all connectivities have the same value keep them all
                    #print "DEBUG: all connectivities of ", i," have the same value:", prob_list
                    zscore_array = np.array([10.00000]*len(prob_list))
                    #print "DEBUG: zscore_array=", zscore_array
                elif len(prob_list) > 2:
                    zscore_array = zscore(prob_list)  # convert probabilities to Z-scores
                elif len(prob_list) == 2:   # if only 2 connectivities exist, set the 2nd to 0.000 to allow you to use high -zcutoff values (<= 0)
                    zscore_array = np.array([1.00000, 0.00000])
                elif len(prob_list) == 0:   # if no connectivity was made create an empty array
                    zscore_array = np.array([])
                
                new_i_iminus1_dict[i] = []  # add i to the new dictionary
                for triplet, Zscore in zip(i_iminus1_dict[i], zscore_array):
                    ## ONLY IF THE Z-SCORE OF THE CONNECTIVITY IS GREATER THAN THE CUTOFF
                    if Zscore >= args.ZSCORE_RESONANCE_MATCH_CUTOFF or approx_equal(Zscore, args.ZSCORE_RESONANCE_MATCH_CUTOFF):
                        new_i_iminus1_dict[i].append(triplet)
            
            del i_iminus1_dict
            i_iminus1_dict = new_i_iminus1_dict
            del new_i_iminus1_dict
            
            return i_iminus1_dict, i_iminus1_complete_dict
        
        
        if not args.POOL_CONNECTIVITIES_FILE:
            print("Calculating Connectivities...")
            i_iminus1_dict, i_iminus1_complete_dict = create_connectivity_dictionaries()
            i_iminus1_nonativepeaks_dict, i_iminus1_complete_nonativepeaks_dict = create_connectivity_dictionaries(clean_native_peaks=True)
        else:
            print("Loading Connectivities from file", args.POOL_CONNECTIVITIES_FILE)
            i_iminus1_pool_dict = {} # same dictionary with the pool of possible connectivities remained after filtering using absolute consensus matches of a previous run
            i_iminus1_complete_dict = {} # the same dictionary but with all possible connectivities, including those below args.RESONANCE_MATCH_CUTOFF
        
            ## I THINK THIS IS NOT NECESSARY ##
            #with open(args.POOL_CONNECTIVITIES_FILE, 'r') as f:
            #    connectivities_file_contents = f.readlines()
            #for line in connectivities_file_contents[1:]:
            #    word_list = line.split()
            #    TOCSY_aaindex = word_list[0]
            #    i_iminus1_dict[TOCSY_aaindex] = []
            #    values_string = ''.join(word_list[1:])
            #    elements_string = re.sub('[\(\)\'\"]', '',  values_string ).split("),(")[0]
            #    elements_list = elements_string.split(",")
            #    for aa, occupancy, TOCSY_resonnum in zip(elements_list[0::3], elements_list[1::3], elements_list[2::3]):
            #        #print aa, occupancy, TOCSY_resonnum
            #        i_iminus1_dict[TOCSY_aaindex].append((aa, int(occupancy), int(TOCSY_resonnum)))
            ###################################
            
            # Now load the pool connectivities file to calculate correct probabilities
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
                    #print aa, occupancy, TOCSY_resonnum
                    i_iminus1_pool_dict[TOCSY_aaindex].append((aa, int(occupancy), int(TOCSY_resonnum)))
            
            # Now load the complete connectivities file to calculate correct probabilities
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
                    #print aa, occupancy, TOCSY_resonnum
                    i_iminus1_complete_dict[TOCSY_aaindex].append((aa, int(occupancy), int(TOCSY_resonnum)))

            
            # IF THE USER HAS SPECIFIED TOLERANCES FOR MATCHING, USE THEM TO FIND NEW CONNECTIVITIES FOR THOSE AAIGs THAT
            # WERE NOT CORRECTED
            if '-tolH' in sys.argv and '-tolC' in sys.argv:
                print("Calculating Connectivities using H tolerance", args.tolH, "and C torelance", args.tolC)
                new_i_iminus1_dict, new_i_iminus1_complete_dict = create_connectivity_dictionaries()
                for TOCSYaaindex in list(new_i_iminus1_complete_dict.keys()):
                    if TOCSYaaindex not in list(i_iminus1_complete_dict.keys()):  # if there was not connectivity for this TAAIG add it to the dictionaries               
                        i_iminus1_complete_dict[TOCSYaaindex] = new_i_iminus1_complete_dict[TOCSYaaindex]   # update the connectivities
                        i_iminus1_pool_dict[TOCSYaaindex] = new_i_iminus1_complete_dict[TOCSYaaindex]
                    elif len(new_i_iminus1_complete_dict[TOCSYaaindex]) > len(i_iminus1_complete_dict[TOCSYaaindex]): # if there were found extra connectivities using the new tolerances
                        if len(i_iminus1_pool_dict[TOCSYaaindex]) == 1 and len(i_iminus1_complete_dict[TOCSYaaindex]) > 1:
                            are_all_unique = True   # do all Tindices in the complete dict have a single connectivity in the pool
                            for tmp_triplet in i_iminus1_complete_dict[TOCSYaaindex]:
                                tmp_TOCSYaaindex = tmp_triplet[0]
                                if tmp_TOCSYaaindex in list(i_iminus1_pool_dict.keys()) and len(i_iminus1_pool_dict[tmp_TOCSYaaindex]) > 1:
                                    are_all_unique = False
                                    break
                            if are_all_unique == False: # this means that i_iminus1_pool_dict[TOCSYaaindex] has been modified and not cleaned from Tindices used elsewhere
                                continue    # leave the connectivity of this TOCSYaaindex unchanged
                            i_iminus1_complete_dict[TOCSYaaindex] = new_i_iminus1_complete_dict[TOCSYaaindex]   # update the connectivities
                            i_iminus1_pool_dict[TOCSYaaindex] = new_i_iminus1_complete_dict[TOCSYaaindex]
                        elif len(i_iminus1_pool_dict[TOCSYaaindex]) == 1 and len(i_iminus1_complete_dict[TOCSYaaindex]) == 1:
                            i_iminus1_complete_dict[TOCSYaaindex] = new_i_iminus1_complete_dict[TOCSYaaindex]   # update the connectivities
                            i_iminus1_pool_dict[TOCSYaaindex] = new_i_iminus1_complete_dict[TOCSYaaindex]


            # CALCULATING CONNECTIVITIES FROM THE POOL OF CONNECTIVITIES USING THE args.RESONANCE_MATCH_CUTOFF
            print("Finding Connectivities above the -mcutoff from the Pool of Connectivities...")
            # connectivities_mdict is a multidimensional dictionary with structure:
            # TOCSY AAIG => NOESY AAIG => [number of matched TOCSY w2,w3 resonances in NOESY for that NOESY AAIG, total number of TOCSY w2,w3 resonances for that TOCSY AAIG]
            i_iminus1_dict = {} # a dictionary containing all possible connectivities of every TOCSY AAIG above args.RESONANCE_MATCH_CUTOFF; has keys the TOCSY aa indices and values lists of triplets (tuples) \
                                # consisting of (NOESYaaindex, occupancy, numOfResonances)
            for TOCSYaaindex in list(i_iminus1_pool_dict.keys()):
                # Iterate over all NOESY aa indices and, if applicable, keep those that match in all resonance pairs w2,w3
                for triplet in i_iminus1_pool_dict[TOCSYaaindex]:
                    numOfResonances = triplet[2]
                    if float(triplet[1])/numOfResonances >= float(args.RESONANCE_MATCH_CUTOFF):     # adjustable cutoff
                        try:
                            i_iminus1_dict[TOCSYaaindex].append(triplet)
                        except KeyError:
                            i_iminus1_dict[TOCSYaaindex] = [(triplet)]
        
            # FILTER i_iminus1_dict ACCORDING TO THE Z-SCORE CUTOFF
            new_i_iminus1_dict = {}
            for i in list(i_iminus1_dict.keys()):
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
                    ## ONLY IF THE Z-SCORE OF THE CONNECTIVITY IS GREATER THAN THE CUTOFF
                    if Zscore >= args.ZSCORE_RESONANCE_MATCH_CUTOFF or approx_equal(Zscore, args.ZSCORE_RESONANCE_MATCH_CUTOFF):
                        new_i_iminus1_dict[i].append(triplet)
            
            del i_iminus1_dict
            i_iminus1_dict = new_i_iminus1_dict
            del new_i_iminus1_dict
            
            
        ## Convert occupancies to probabilities
        i_iminus1_normProbabilities_dict = {}   # dict of the form TAAIG(i)->[(TAAIG(i-1), probability), (...), ...]; eventually will contain only the matches above the cutoff
        aaindex_weightSum_dict = {}
        for TOCSY_aaindex in list(i_iminus1_complete_dict.keys()):
            weight_sum = 0
            for triplet in i_iminus1_complete_dict[TOCSY_aaindex]:  # use all possible matches to calculate the weight_sum
                occupancy = triplet[1]
                TOCSY_resonnum = triplet[2]
                weight = float(occupancy)/TOCSY_resonnum   # occupancy/total number of TOCSY resonances
                weight_sum += weight
            aaindex_weightSum_dict[TOCSY_aaindex] = weight_sum
        
        for TOCSY_aaindex in list(i_iminus1_dict.keys()):     # use only the matches above the cutoff to save theirs probabilities
            #print "DEBUG: TOCSY_aaindex=", TOCSY_aaindex, "aaindex_weightSum_dict[TOCSY_aaindex] =",aaindex_weightSum_dict[TOCSY_aaindex]
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
            
            all_chainScore_set = set()    # a list of lists containing the possible connectivities and the overall score as the last element
            with open("connectivities_cutoff_"+str(args.RESONANCE_MATCH_CUTOFF)+"_zmcutoff"+str(args.ZSCORE_RESONANCE_MATCH_CUTOFF), 'w') as f:
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
                with open("connectivities_nonativepeaks_cutoff_" + str(args.RESONANCE_MATCH_CUTOFF)+"_zmcutoff" +
                          str(args.ZSCORE_RESONANCE_MATCH_CUTOFF), 'w') as f:
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
                for i in list(i_iminus1_dict.keys()):
                    build_Chain_Tree(i, all_chainScore_set, i_iminus1_dict, i_iminus1_normProbabilities_dict, args.MAX_PEPTIDE_LENGTH, args.MIN_PEPTIDE_LENGTH)
                
                ##
                ## Write chains to files
                ##
                
                print("Saving all the chains to chains.list file ...")
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
                            #print "Overall Score = ",score," chain = ", chain
                            f.write("Overall Score = "+str(score)+" chain = "+str(chain)+"\n")
                
                write_chains_file(all_chainScore_set, "chains.list")
                
                ##
                ##  Remove redudancy from the chain list
                ### Bash script to help in validating the redundancy removal
                #while read line;
                #    do chain=$(perl -pi -e "s/.* chain = \((.*)\)/\1/" <<< $line);
                #    score=$(awk '{print $4}' <<< $line);
                #    echo "---> redundant chain = $chain";
                #    grep "Score = $score chain = ($chain" chains.non-redundant.list | grep "$chain";
                #done < chains.redundant.list > redundant_lines.txt
                
                print("Removing redundant chains ...")
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
                # results is a list of sets with chains to be removed
                chainScores2remove_set = set()
                for s in results:
                    chainScores2remove_set.union(s)
                
                # remove subchains to keep only the longest, unique, non-overlapping ones
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
        
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        #   AUTOMATIC ASSIGNMENT OF AMINO ACID TYPES TO TOCSY RESONANCES    #
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        
        if not args.POOL_AA_TYPES_FILE and not args.COMPLETE_AA_TYPES_FILE:
            ## LOAD 1D BMRB HISTOGRAMS
            print("\nLoading BMRB chemical shift histograms...")
            aa_carbon_binDensityList_mdict = tree()   # multidict with amino acid name -> Carbon atom type -> [array of bin limits, array of probability density]
            aa_hydrogen_binDensityList_mdict = tree()   # multidict with amino acid name -> Hydrogen atom type -> [array of bin limits, array of probability density]
            fnames = os.listdir(CHAINS_BIN_DIR+"/../databases/histograms/")
            fpattern = re.compile("[A-Z]{3}_[A-Z0-9]+_hist.txt$")
            hist_files_list = list(filter(fpattern.search, fnames))
            for hist_file in hist_files_list:
                aa = hist_file.split("_")[0]
                atom = hist_file.split("_")[1]
                #try:
                if atom in allowed_aa_atoms_dict[aa]:   # load histograms of the allowed atom types only
                    bin_list = []
                    density_list = []
                    with open(CHAINS_BIN_DIR+"/../databases/histograms/"+hist_file, 'r') as f:
                        for line in f:
                            word_list = line.split()
                            bin_list.append(float(word_list[0]))
                            density_list.append(float(word_list[1]))
                            # print "DEBUG: saving bin ", word_list[0], "and density", word_list[1]
                    bin_array = np.array(bin_list)
                    density_array = np.array(density_list)/sum(density_list)
                    # use C & H resonances to assign aa type to the (i-1) residue 
                    if atom[0] == "C":
                        aa_carbon_binDensityList_mdict[aa][atom] = [bin_array, density_array]
                    elif atom[0] == "H":
                        aa_hydrogen_binDensityList_mdict[aa][atom] = [bin_array, density_array]
                #except KeyError:   # active this exception if use N, HN resonance to assign aa type to the ith residue
                #    if aa == "PRO": # Prolines cannot be detected in root spectrum (N-H HSQC) because they don't have "HN" atom
                #        continue
            
            
            if args.USE_2D_HISTOGRAMS == True:
                ## IF ASKED IN THE COMMAND LINE, LOAD 2D BMRB HISTOGRAMS, TOO
                print("\nLoading 2D BMRB chemical shift histograms...")
                aa_CHpair_binProbabilityList_mdict = tree()   # multidict with amino acid name -> Carbon atom type -> [array of bin limits, array of probability density]
                fnames = os.listdir(CHAINS_BIN_DIR+"/../databases/histograms/")
                fpattern = re.compile("[A-Z]{3}_[A-Z0-9-]+_correlated_2Dhist.smoothed.txt$")
                hist_files_list = list(filter(fpattern.search, fnames))
                for hist_file in hist_files_list:
                    #print "Loading histogram", hist_file 
                    aa = hist_file.split("_")[0]
                    CH_pair = hist_file.split("_")[1]
                    #try:
                    if CH_pair in list(aa_CHpair_2Dhist_mdict[aa].keys()):   # load histograms of the allowed atom types only
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
                        # use C & H resonances to assign aa type to the (i-1) residue
                        aa_CHpair_binProbabilityList_mdict[aa][CH_pair] = [x_bin_array, y_bin_array, probability_array]
                    #except KeyError:   # active this exception if use N, HN resonance to assign aa type to the ith residue
                    #    if aa == "PRO": # Prolines cannot be detected in root spectrum (N-H HSQC) because they don't have "HN" atom
                    #        continue
            
            
            # sys.exit(0) # TEMPTORARILY DEACTIVATE AA TYPE PREDICTION.
            iaaindex_iminus1aaTypesProbTupleList_dict = OrderedDict() # ordereddict with keys the AAIG of residue i and values lists of tuples of the form (residue i-1 matching aa type, average probability)
            Num_of_TOCSY_resonances = 0 # save here the number of TOCSY resonances for a particular AAIG
            previous_TOCSY_aaindex = None
            possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list = []
            carbon_groups_list = []  # list of tuples of the form (aaindex, H resonance, C resonance, N resonance, HN resonance, carbon group)
            aaindex_carbonGroupsTuple_dict = {} # dict of the form: aaindex -> list of the form (aaindex, H resonance, C resonance, N resonance, HN resonance, carbon group)
            counter = 0
            previous_aaindex = ""
            sorted_TOCSY_contents = sorted(TOCSY_contents, key=itemgetter(0))
            print("sorted_TOCSY_contents=", sorted_TOCSY_contents)
            #sys.exit(1)
            #sorted_TOCSY_contents= [
            #    ['A6', '4.148', '61.783', '131.690', '8.959'], ['A6', '1.829', '33.205', '131.701', '8.960'], ['A6', '1.038', '58.474', '131.692', '8.961'],
            #    ['A6', '0.756', '57.876', '131.676', '8.961']
            #    ]
            for TOCSY_words_list in sorted_TOCSY_contents:
                try:
                    # ignore 1st column
                    TOCSY_aaindex=TOCSY_words_list[0]    # residue i
                    print("DEBUG: TOCSY_aaindex=",TOCSY_aaindex,"previous_TOCSY_aaindex=",previous_TOCSY_aaindex)
                    if previous_TOCSY_aaindex != None and TOCSY_aaindex != previous_TOCSY_aaindex: # if we look at a different AAIG in TOCSY, print the matches and occupancies
                                                                                                   # of the previous TOCSY AAIG and clear matchingNOESYaaindex_occupancy_dict
                        
                        print("Assigning possible aa types to the aa upstream of AAIG",previous_TOCSY_aaindex)
                        #print "DEBUG: switched from AAIG", previous_TOCSY_aaindex, "to", TOCSY_aaindex
                        #print "DEBUG: possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list=",possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list
                        #print "DEBUG: Num_of_TOCSY_resonances=",Num_of_TOCSY_resonances
                        # GROUP THE CARBON RESONANCE OF EACH TAAIG
                        new_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list = group_carbons(possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list)
                        iaaindex_iminus1aaTypesProbTupleList_dict[previous_TOCSY_aaindex] = get_aatypes_from_all_H_C_resonpairs(new_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list, Num_of_TOCSY_resonances, args)
                        #print "DEBUG: previous_TOCSY_aaindex=",previous_TOCSY_aaindex,"iaaindex_iminus1aaTypesProbTupleList_dict[previous_TOCSY_aaindex]=",iaaindex_iminus1aaTypesProbTupleList_dict[previous_TOCSY_aaindex]
                        possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list = [] # list of tuples of the form (aa type, probability, H_resonance, C_resonance, TOCSY_reson_index) containing all matching sets of H,C resonances
                        aaindex_carbonGroupsTuple_dict[TOCSY_aaindex] = carbon_groups_list
                        #print "DEBUG: aaindex_carbonGroupsTuple_dict=", aaindex_carbonGroupsTuple_dict
                        carbon_groups_list = []  # list of tuples of the form (TOCSY_aaindex, H resonance, C resonance, N resonance, HN, resonance, carbon group)
                        
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
                        valid_matches_list = get_aatypes_from_H_C_resonpair_2Dhist(TOCSY_H_resonance, TOCSY_C_resonance, Num_of_TOCSY_resonances, aa_CHpair_binProbabilityList_mdict, aa_carbon_binDensityList_mdict)
                    else:
                        valid_matches_list = get_aatypes_from_H_C_resonpair(TOCSY_H_resonance, TOCSY_C_resonance, Num_of_TOCSY_resonances, aa_carbon_binDensityList_mdict, aa_hydrogen_binDensityList_mdict, args)
                    #if TOCSY_aaindex == "Y130":
                    #print "DEBUG: valid_matches_list=",valid_matches_list
                    possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list.extend(valid_matches_list)
                    carbon_groups_list.append([TOCSY_aaindex, TOCSY_H_resonance, TOCSY_C_resonance, TOCSY_N_resonance, TOCSY_HN_resonance, None])
                
                    previous_TOCSY_aaindex = TOCSY_aaindex
                except (ValueError, IndexError):
                    print("WARNING: the 3rd and 4th elements of the following TOCSY file line are not numbers:")
                    print("TOCSY file line:", TOCSY_words_list)
                    continue
            # THIS IS FOR THE LAST TOCSY INDEX
            #if not counter == 10:
            print("Assigning possible aa types to the last aa, which is upstream of AAIG",previous_TOCSY_aaindex)
            #print "DEBUG: possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list=",possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list
            # GROUP THE CARBON RESONANCE OF EACH TAAIG
            new_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list = group_carbons(possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list)
            iaaindex_iminus1aaTypesProbTupleList_dict[previous_TOCSY_aaindex] = get_aatypes_from_all_H_C_resonpairs(new_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list, Num_of_TOCSY_resonances, args)
            #print "DEBUG: previous_TOCSY_aaindex=",previous_TOCSY_aaindex,"iaaindex_iminus1aaTypesProbTupleList_dict[previous_TOCSY_aaindex]=",iaaindex_iminus1aaTypesProbTupleList_dict[previous_TOCSY_aaindex]
            
            ##### NEW WAY: CALCULATE Z-SCORE FOR EACH TOCSY INDEX INDIVIDUAL #####
            ## Calculate Z-scores from aa type prediction probabilites and keep only those predictions above the cutoff
            iaaindex_iminus1aaTypesZscoreTupleList_dict = OrderedDict() # ordereddict with keys the AAIG of residue i and values lists of tuples of the form (residue i-1 matching aa type, Z-score)
            iaaindex_iminus1aaTypesCutoffProbTupleList_dict = OrderedDict() # contains only the tuples (residue i-1 matching aa type, probability) with Z-score above the cutoff
            for iaaindex, iminus1aaTypesProbTuple_list in list(iaaindex_iminus1aaTypesProbTupleList_dict.items()):
                # DECIDE THE PROBABILITY THRESHOLD TO DISCARD AA TYPE PREDICTIONS BELOW IT
                try:
                    max_prob = np.max([duplet[1] for duplet in iminus1aaTypesProbTuple_list])   # highest aa type probability
                except ValueError:  # if no aa type predictions exist for this TAAIG group
                    max_prob = 0
                if max_prob > 10e-10:   # Dedice the probability threshold
                    prob_threshold = 1000
                elif max_prob > 10e-20:
                    prob_threshold = 10000
                elif max_prob <= 10e-20:
                    prob_threshold = 100000
                print("DEBUG: max_prob", max_prob, "prob_threshold", prob_threshold, "args.LOG=", args.LOG)
                aatype_list, prob_list = [], []
                for duplet in iminus1aaTypesProbTuple_list:
                    if args.DELETE_AA_TYPE_PREDICTIONS == True:     # delete low probability aa type predictions
                        try:
                            if max_prob/float(duplet[1]) < prob_threshold:
                                aatype_list.append(duplet[0])
                                prob_list.append(duplet[1])  # unsorted list of probabilities
                        except ZeroDivisionError:
                            print("DEBUG: ZeroDivisionError ", max_prob, "/", float(duplet[1]))
                    elif args.DELETE_AA_TYPE_PREDICTIONS == False:
                        aatype_list.append(duplet[0])
                        prob_list.append(duplet[1])  # unsorted list of probabilities
                print("DEBUG: prob_list=", prob_list)
                if len(prob_list) == 1:     
                    zscore_array = np.array([10.00000])     # default Z-score in the case of just one prediction
                elif len(prob_list) > 2:
                    if args.LOG == True:
                        try:
                            ratio = float(np.max(prob_list))/np.min(prob_list)    # if the min probability is at least 3 orders of magnitude smaller, convert them to logarithmic scale
                        except ZeroDivisionError:
                            print("DEBUG: ZeroDivisionError ", np.max(prob_list), "/", np.min(prob_list)) 
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
                print("DEBUG2: saving zscore_array = ", zscore_array)
                iminus1aaTypesZscoreTuple_list = []
                iminus1aaTypesCutoffProbTuple_list = []
                for aatype, Zscore, prob in zip(aatype_list, zscore_array, prob_list):
                    ## ONLY IF THE Z-SCORE OF THE AA TYPE PREDICTION IS GREATER THAN THE CUTOFF AND THE AA TYPE PREDICTIONS ARE AT LEAST args.MIN_NUM_OF_PREDICTIONS, THEN INCLUDE IT IN THE LIST
                    if Zscore > args.ZSCORE_ASSIGNMENT_CUTOFF and len(zscore_array) >= args.MIN_NUM_OF_PREDICTIONS:
                        aatypeProb_tuple = (aatype, prob)
                        aatypeZscore_tuple = (aatype, Zscore)  # replace the probability with the respective Z-score
                        iminus1aaTypesZscoreTuple_list.append(aatypeZscore_tuple)
                        iminus1aaTypesCutoffProbTuple_list.append(aatypeProb_tuple) # but also save the aatype probabilities that satisfy the Z-score cutoff
                    ## IF THE AA TYPE PREDICTIONS ARE LESS THAN args.MIN_NUM_OF_PREDICTIONS, THEN INCLUDED ALL OF THEM IN THE LIST
                    elif len(zscore_array) < args.MIN_NUM_OF_PREDICTIONS:
                        aatypeProb_tuple = (aatype, prob)
                        aatypeZscore_tuple = (aatype, Zscore)  # replace the probability with the respective Z-score
                        iminus1aaTypesZscoreTuple_list.append(aatypeZscore_tuple)
                        iminus1aaTypesCutoffProbTuple_list.append(aatypeProb_tuple) # but also save the aatype probabilities that satisfy the Z-score cutoff
                iaaindex_iminus1aaTypesZscoreTupleList_dict[iaaindex] = iminus1aaTypesZscoreTuple_list
                iaaindex_iminus1aaTypesCutoffProbTupleList_dict[iaaindex] = iminus1aaTypesCutoffProbTuple_list
            
            
            # Save ALL amino acid type predictions with the respective probabilities
            #with open("amino_acid_type_prediction_probabilities.mcutoff"+str(args.RESONANCE_MATCH_CUTOFF)+".acutoff"+str(args.ASSIGNMENT_CUTOFF)+".zacutoff"+str(args.ZSCORE_ASSIGNMENT_CUTOFF)
                      #+".wH"+str(args.H_weight)+".wC"+str(args.C_weight)+".tolH"+str(args.tolH)+".tolC"+str(args.tolC)+".rtolH"+str(args.rtolH)+".rtolN"+str(args.rtolN), 'w') as f:
            with open("amino_acid_type_prediction_probabilities", 'w') as f:
                f.write("i AAIG\tpossible i-1 aa types\n")
                print("Probabilities:")
                for i_aaindex in list(iaaindex_iminus1aaTypesProbTupleList_dict.keys()):
                    print(i_aaindex,"---> (i-1) aa type",sorted(iaaindex_iminus1aaTypesProbTupleList_dict[i_aaindex], key=itemgetter(1), reverse=True))
                    f.write(i_aaindex + " ---> (i-1) aa type " + str(sorted(iaaindex_iminus1aaTypesProbTupleList_dict[i_aaindex], key=itemgetter(1), reverse=True)) + "\n")
            
            # Save amino acid type predictions that are above the Z-Score cutoff, along with the respective probabilities
            #with open("amino_acid_type_prediction_probabilities.mcutoff"+str(args.RESONANCE_MATCH_CUTOFF)+".acutoff"+str(args.ASSIGNMENT_CUTOFF)+".zacutoff"+str(args.ZSCORE_ASSIGNMENT_CUTOFF)
                      #+".wH"+str(args.H_weight)+".wC"+str(args.C_weight)+".tolH"+str(args.tolH)+".tolC"+str(args.tolC)+".rtolH"+str(args.rtolH)+".rtolN"+str(args.rtolN), 'w') as f:
            with open("amino_acid_type_prediction_probabilities_above_cutoff.zacutoff"+str(args.ZSCORE_ASSIGNMENT_CUTOFF), 'w') as f:
                f.write("i AAIG\tpossible i-1 aa types\n")
                print("Probabilities:")
                for i_aaindex in list(iaaindex_iminus1aaTypesCutoffProbTupleList_dict.keys()):
                    print(i_aaindex,"---> (i-1) aa type",sorted(iaaindex_iminus1aaTypesCutoffProbTupleList_dict[i_aaindex], key=itemgetter(1), reverse=True))
                    f.write(i_aaindex + " ---> (i-1) aa type " + str(sorted(iaaindex_iminus1aaTypesCutoffProbTupleList_dict[i_aaindex], key=itemgetter(1), reverse=True)) + "\n")
            
            # Save amino acid type predictions with the respective Z-scores
            #with open("amino_acid_type_prediction_Z-scores.mcutoff"+str(args.RESONANCE_MATCH_CUTOFF)+".acutoff"+str(args.ASSIGNMENT_CUTOFF)+".zacutoff"+str(args.ZSCORE_ASSIGNMENT_CUTOFF)+
                      #".wH"+str(args.H_weight)+".wC"+str(args.C_weight)+".tolH"+str(args.tolH)+".tolC"+str(args.tolC)+".rtolH"+str(args.rtolH)+".rtolN"+str(args.rtolN), 'w') as f:
            with open("amino_acid_type_prediction_Z-scores.zacutoff"+str(args.ZSCORE_ASSIGNMENT_CUTOFF), 'w') as f:
                f.write("i AAIG\tpossible i-1 aa types\n")
                print("Z-scores:")
                for i_aaindex in list(iaaindex_iminus1aaTypesZscoreTupleList_dict.keys()):
                    print(i_aaindex,"---> (i-1) aa type",sorted(iaaindex_iminus1aaTypesZscoreTupleList_dict[i_aaindex], key=itemgetter(1), reverse=True))
                    f.write(i_aaindex + " ---> (i-1) aa type " + str(sorted(iaaindex_iminus1aaTypesZscoreTupleList_dict[i_aaindex], key=itemgetter(1), reverse=True)) + "\n")
        else:
            ##########################################################################################################################################################################
            ##  LOAD AMINO ACID TYPE PREDICTIONS FOR FASTER DEBUGGING
            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            print("")    # change line after the previous progress bar
            print("Loading amino acid prediction files ...")
            minimum_Zscore = 10000000
            iaaindex_iminus1aaTypesProbPoolTupleList_dict = OrderedDict() # ordereddict with keys the AAIG of residue i and values lists of tuples of the form
                                                                        # (residue i-1 matching aa type, Z-score)
            with open(args.POOL_AA_TYPES_FILE, 'r') as f:
                pool_aa_type_file_contents = f.readlines()
                for line in pool_aa_type_file_contents[1:]:
                    #print "DEBUG: line=",line
                    word_list = re.sub('---> \(i-1\) aa type', '', line).split()
                    key = word_list[0]
                    iaaindex_iminus1aaTypesProbPoolTupleList_dict[key] = []
                    values_string = ''.join(word_list[1:])
                    elements_string = re.sub('[\(\)\[\]\'\"]', '',  values_string ).split("),(")[0]
                    elements_list = elements_string.split(",")
                    for aa, Zscore in zip(elements_list[0::2], elements_list[1::2]):
                        #print "aa=",aa,"probability=",probability
                        duplet = (aa, float(Zscore))
                        iaaindex_iminus1aaTypesProbPoolTupleList_dict[key].append(duplet)
            
            ## Now Calculate Z-scores from the pool of aa type prediction probabilites and keep only those predictions above the cutoff
            iaaindex_iminus1aaTypesZscoreTupleList_dict = OrderedDict() # ordereddict with keys the AAIG of residue i and values lists of tuples of the form (residue i-1 matching aa type, Z-score)
            iaaindex_iminus1aaTypesCutoffProbTupleList_dict = OrderedDict() # contains only the tuples (residue i-1 matching aa type, probability) with Z-score above the cutoff
            for iaaindex, iminus1aaTypesProbTuple_list in list(iaaindex_iminus1aaTypesProbPoolTupleList_dict.items()):
                # DECIDE THE PROBABILITY THRESHOLD TO DISCARD AA TYPE PREDICTIONS BELOW IT
                try:
                    max_prob = np.max([duplet[1] for duplet in iminus1aaTypesProbTuple_list])   # highest aa type probability
                except ValueError:  # if no aa type predictions exist for this TAAIG group
                    max_prob = 0
                if max_prob > 10e-10:   # Dedice the probability threshold
                    prob_threshold = 1000
                elif max_prob > 10e-20:
                    prob_threshold = 10000
                elif max_prob <= 10e-20:
                    prob_threshold = 100000
                print("DEBUG: max_prob", max_prob, "prob_threshold", prob_threshold, "args.LOG=", args.LOG)
                aatype_list, prob_list = [], []
                for duplet in iminus1aaTypesProbTuple_list:
                    if args.DELETE_AA_TYPE_PREDICTIONS == True:     # delete low probability aa type predictions
                        try:
                            if max_prob/float(duplet[1]) < prob_threshold:
                                aatype_list.append(duplet[0])
                                prob_list.append(duplet[1])  # unsorted list of probabilities
                        except ZeroDivisionError:
                            print("DEBUG: ZeroDivisionError ", max_prob, "/", float(duplet[1]))
                    elif args.DELETE_AA_TYPE_PREDICTIONS == False:
                        aatype_list.append(duplet[0])
                        prob_list.append(duplet[1])  # unsorted list of probabilities
                print("DEBUG: prob_list=", prob_list)
                if len(prob_list) == 1:     
                    zscore_array = np.array([10.00000])     # default Z-score in the case of just one prediction
                elif len(prob_list) > 1:
                    if args.LOG == True:
                        try:
                            ratio = float(np.max(prob_list))/np.min(prob_list)    # if the min probability is at least 3 orders of magnitude smaller, convert them to logarithmic scale
                        except ZeroDivisionError:
                            print("DEBUG: ZeroDivisionError ", np.max(prob_list), "/", np.min(prob_list)) 
                            ratio = -1.0
                        if ratio > 1000:
                            zscore_array = zscore(np.log(prob_list))  # convert probabilities the logarithms and then to Z-scores
                        else:
                            zscore_array = zscore(prob_list)  # convert probabilities to Z-scores
                    elif args.LOG == False:
                        zscore_array = zscore(prob_list)  # convert probabilities to Z-scores
                else:   # if no aa type prediction was made create an empty array
                    zscore_array = np.array([])
                # print "DEBUG1: saving zscore_array = ", zscore_array
                iminus1aaTypesZscoreTuple_list = []
                iminus1aaTypesCutoffProbTuple_list = []
                for aatype, Zscore, prob in zip(aatype_list, zscore_array, prob_list):
                    if Zscore < minimum_Zscore:
                        minimum_Zscore = Zscore
                    ## ONLY IF THE Z-SCORE OF THE AA TYPE PREDICTION IS GREATER THAN THE CUTOFF AND THE AA TYPE PREDICTIONS ARE AT LEAST args.MIN_NUM_OF_PREDICTIONS, THEN INCLUDE IT IN THE LIST
                    if Zscore > args.ZSCORE_ASSIGNMENT_CUTOFF and len(zscore_array) >= args.MIN_NUM_OF_PREDICTIONS:
                        aatypeProb_tuple = (aatype, prob)
                        aatypeZscore_tuple = (aatype, Zscore)  # replace the probability with the respective Z-score
                        iminus1aaTypesZscoreTuple_list.append(aatypeZscore_tuple)
                        iminus1aaTypesCutoffProbTuple_list.append(aatypeProb_tuple) # but also save the aatype probabilities that satisfy the Z-score cutoff
                    ## IF THE AA TYPE PREDICTIONS ARE LESS THAN args.MIN_NUM_OF_PREDICTIONS, THEN INCLUDED ALL OF THEM IN THE LIST
                    elif len(zscore_array) < args.MIN_NUM_OF_PREDICTIONS:
                        aatypeProb_tuple = (aatype, prob)
                        aatypeZscore_tuple = (aatype, Zscore)  # replace the probability with the respective Z-score
                        iminus1aaTypesZscoreTuple_list.append(aatypeZscore_tuple)
                        iminus1aaTypesCutoffProbTuple_list.append(aatypeProb_tuple) # but also save the aatype probabilities that satisfy the Z-score cutoff
                iaaindex_iminus1aaTypesZscoreTupleList_dict[iaaindex] = iminus1aaTypesZscoreTuple_list
                iaaindex_iminus1aaTypesCutoffProbTupleList_dict[iaaindex] = iminus1aaTypesCutoffProbTuple_list
            
            
            iaaindex_iminus1aaTypesProbTupleList_dict = OrderedDict()   # ordereddict with keys the AAIG of residue i and values lists of tuples of the form
                                                                        # (residue i-1 matching aa type, average probability)
            with open(args.COMPLETE_AA_TYPES_FILE, 'r') as f:
                complete_aa_type_file_contents = f.readlines()
                for line in complete_aa_type_file_contents[1:]:
                    #print "DEBUG: line=",line
                    word_list = re.sub('---> \(i-1\) aa type', '', line).split()
                    key = word_list[0]
                    iaaindex_iminus1aaTypesProbTupleList_dict[key] = []
                    values_string = ''.join(word_list[1:])
                    elements_string = re.sub('[\(\)\[\]\'\"]', '',  values_string ).split("),(")[0]
                    elements_list = elements_string.split(",")
                    for aa, hist_prob in zip(elements_list[0::2], elements_list[1::2]):
                        #print "aa=",aa,"probability=",probability
                        duplet = (aa, float(hist_prob))
                        iaaindex_iminus1aaTypesProbTupleList_dict[key].append(duplet)
            
            #for i_aaindex in iaaindex_iminus1aaTypesZscoreTupleList_dict.keys():
            #    print i_aaindex,"---> (i-1) aa type",sorted(iaaindex_iminus1aaTypesZscoreTupleList_dict[i_aaindex], key=itemgetter(1), reverse=True)
            ##########################################################################################################################################################################
            
            with open("amino_acid_type_prediction_Z-scores.zacutoff"+str(args.ZSCORE_ASSIGNMENT_CUTOFF), 'w') as f:
                f.write("i AAIG\tpossible i-1 aa types\n")
                print("Z-scores:")
                for i_aaindex in list(iaaindex_iminus1aaTypesZscoreTupleList_dict.keys()):
                    print(i_aaindex,"---> (i-1) aa type",sorted(iaaindex_iminus1aaTypesZscoreTupleList_dict[i_aaindex], key=itemgetter(1), reverse=True))
                    f.write(i_aaindex + " ---> (i-1) aa type " + str(sorted(iaaindex_iminus1aaTypesZscoreTupleList_dict[i_aaindex], key=itemgetter(1), reverse=True)) + "\n")
        
        #sys.exit(1)
        ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        ##              CALCULATE BAYESIAN STATISTICS                        #
        ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        
        ## CALCULATE THE PROBABILITY P[AATYPE(I-1)]
        ############################### NEW WAY #####################################
        template_sequence_string = ""
        with open(args.template_sequence_file, 'r') as f:
            for line in f:
                if not re.match("^>", line):
                    template_sequence_string += re.sub(r"[^A-Z]", "", line)
        template_sequence_list = list(template_sequence_string)
        
        aatype_P_dict = {}    # dict with the P[aatype(i-1)] for every aatype(i-1): the percentage of each aatype in the template sequence
        for aa_type in list(aatype_maxH_C_pairs_dict.keys()):
            aatype_P_dict[aa_type] = template_sequence_list.count(aa3to1_dict[aa_type]) / float(len(template_sequence_list))
        
        ## DO A MATHEMATIC TRANSFORM ON P
        # iaaindex_iminus1aaTypesProbTupleList_dict is an ordereddict with keys the AAIG of residue i and values lists of tuples of the form
        # (residue i-1 matching aa type, average probability)
        all_prob_list = []
        for iaaindex, iminus1aaTypesProbTuple_list in list(iaaindex_iminus1aaTypesProbTupleList_dict.items()):
            for duplet in iminus1aaTypesProbTuple_list:
                prob = duplet[1]
                all_prob_list.append(prob)
        
        all_prob_array = np.array(all_prob_list)
        transformed_all_prob_array = transform(all_prob_array)
        
        iaaindex_iminus1aaTypesTransformedProbTupleList_dict = OrderedDict()
        array_index = 0
        for iaaindex, iminus1aaTypesProbTuple_list in list(iaaindex_iminus1aaTypesProbTupleList_dict.items()):
            iaaindex_iminus1aaTypesTransformedProbTupleList_dict[iaaindex] = []
            for duplet in iminus1aaTypesProbTuple_list:
                aatype = duplet[0]
                trans_prob = transformed_all_prob_array[array_index]
                iaaindex_iminus1aaTypesTransformedProbTupleList_dict[iaaindex].append((aatype, trans_prob))
                array_index += 1
        
        def calculate_Paaiminus1_TAAIGi(TAAIGi, aaiminus1):
            """
                Calculate the conditional probability P(AA^{u} | AAIG^{CS}), which quantifies how probable is to get an amino-acid of type u given a set of
                chemical shifts AAIG^{CS}.
                ARGS:
                TAAIGi:    the AAIG or TOCSY Index Group aligned to position i
                aaiminus1:  the aa type if the residue in position i-1
            """
            # aatype_P_dict:                                  probability P[aatype(i-1)], dict with the P[aatype(i-1)] for every aatype(i-1)
            global aatype_P_dict, iaaindex_iminus1aaTypesTransformedProbTupleList_dict
            # aatype_P_dict = shared.getConst('AATYPE_P_DICT')
            # iaaindex_iminus1aaTypesTransformedProbTupleList_dict = shared.getConst('IAAINDEX_IMINUS1AATYPESTRANSFORMEDPROBTUPLELIST_DICT')
            
            PTAAIGi = 0    # P[TAAIG(i)] = Sum{P[TAAIG(i)|aatype(i-1)]}
            for duplet in iaaindex_iminus1aaTypesTransformedProbTupleList_dict[TAAIGi]:
                #print "DEBUG: duplet=", duplet
                aatype = duplet[0]
                PTAAIGi += duplet[1]
                if aatype == aaiminus1:
                    PTAAIGi_aaiminus1 = duplet[1]  # this is P[TAAIG(i)|aatype(i-1)]
            
            #print "DEBUG: aaiminus1=",aaiminus1
            #print "DEBUG: aatype_P_dict=", aatype_P_dict
            Paaiminus1 = aatype_P_dict[aaiminus1]
            #print "DEBUG: PTAAIGi_aaiminus1=",PTAAIGi_aaiminus1
            #print "DEBUG: Paaiminus1=",Paaiminus1
            #print "DEBUG: PTAAIGi=", PTAAIGi
            Paaiminus1_TAAIGi = PTAAIGi_aaiminus1 * Paaiminus1 / PTAAIGi
            #print "DEBUG: Paaiminus1_TAAIGi=", Paaiminus1_TAAIGi
            return Paaiminus1_TAAIGi
        
        
        # Calculate and save the Conditional Probabilities for every aa type prediction
        print("DEBUG: iaaindex_iminus1aaTypesProbTupleList_dict=", iaaindex_iminus1aaTypesProbTupleList_dict)
        iaaindex_iminus1aaType_Condprob_mdict = tree() # multidict of the form AAIG of residue i --> residue i-1 matching aa type --> conditional probability
        for i_aaindex in list(iaaindex_iminus1aaTypesProbTupleList_dict.keys()):
            aatypeCondprob_list = []    # list of the form [(aa type, conditional probability), (aa type, conditional probability), ...]
            for aa_type,prob in sorted(iaaindex_iminus1aaTypesProbTupleList_dict[i_aaindex], key=itemgetter(1), reverse=True):
                condprob = calculate_Paaiminus1_TAAIGi(i_aaindex, aa_type)
                iaaindex_iminus1aaType_Condprob_mdict[i_aaindex][aa_type] = condprob
        
        # Save ALL amino acid type predictions with the respective conditional probabilities
            with open("amino_acid_type_prediction_conditional_probabilities", 'w') as f:
                f.write("i AAIG\tpossible i-1 aa types\n")
                print("Conditional Probabilities:")
                for i_aaindex in list(iaaindex_iminus1aaType_Condprob_mdict.keys()):
                    aa_type_condprob_list = []  # list of the form [(aa type, conditional probability), (aa type, conditional probability), ...]
                    for aa_type in list(iaaindex_iminus1aaType_Condprob_mdict[i_aaindex].keys()):
                        aa_type_condprob_list.append( (aa_type, iaaindex_iminus1aaType_Condprob_mdict[i_aaindex][aa_type]) )
                    print(i_aaindex,"---> (i-1) aa type",sorted(aa_type_condprob_list, key=itemgetter(1), reverse=True))
                    f.write(i_aaindex + " ---> (i-1) aa type " + str(sorted(aa_type_condprob_list, key=itemgetter(1), reverse=True)) + "\n")
            
        
        ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        ##       BUILD TREES OF POSSIBLE PEPTIDE SEQUENCES                   #
        ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        
        if not args.SKIP_CHAINS:
            if os.path.exists('tmp_peptide_folder/') == True and args.RESUME == False:   # if the folder exists and this is not a resumed run, remove it and make a new one
                shutil.rmtree('tmp_peptide_folder/', ignore_errors=True)
                os.makedirs('tmp_peptide_folder/')
            elif os.path.exists('tmp_peptide_folder/') == False and args.RESUME == False: # if the folder does not exists, make a new one
                os.makedirs('tmp_peptide_folder/')
            elif os.path.exists('tmp_peptide_folder/') == False and args.RESUME == True:
                print("ERROR: folder tmp_peptide_folder/ with peptide sequence files does not exist! The process cannot be resumed!")
            # LOAD ALL PREVIOUSLY WRITEN PEPTIDE FILES
            fnames=os.listdir('tmp_peptide_folder/')
            fpattern = re.compile('chainIndex_[0-9]+.pickle.bz2')
            peptide_files_list = list(filter(fpattern.search, fnames))
            chainIndices2skip_set = set()
            for peptide_file in peptide_files_list:
                mo = re.search('chainIndex_([0-9]+).pickle.bz2', peptide_file)
                if mo:
                    chainIndex = int(mo.group(1))
                    chainIndices2skip_set.add(chainIndex)
            
            total_chain_number = len(all_chainScore_list)
            remaining_chainIndex_list = []
            remaining_chainScore_list = []
            # SET SHARED VARIABLES BETWEEN THREADS
            shared.setConst(TOTAL_CHAIN_NUMBER=total_chain_number)
            shared.setConst(ARGS=args)
            #print "DEBUG: = iaaindex_iminus1aaTypesCutoffProbTupleList_dict=", iaaindex_iminus1aaTypesCutoffProbTupleList_dict
            shared.setConst(IAAINDEX_IMINUS1AATYPE_CONDPROB_MULTIDICT = iaaindex_iminus1aaType_Condprob_mdict)
            shared.setConst(IAAINDEX_IMINUS1AATYPESCUTOFFPROBTUPLELIST_DICT=iaaindex_iminus1aaTypesCutoffProbTupleList_dict)
            shared.setConst(AA3TO1_DICT=aa3to1_dict)
            shared.setConst(IAAINDEX_IMINUS1AATYPESTRANSFORMEDPROBTUPLELIST_DICT=iaaindex_iminus1aaTypesTransformedProbTupleList_dict)
            #shared.setConst(IAAINDEX_IMINUS1AATYPESPAATYPETUPLELIST_DICT=iaaindex_iminus1aaTypesPaatypeTupleList_dict)
            shared.setConst(AATYPE_P_DICT=aatype_P_dict)
            shared.setConst(I_IMINUS1_DICT=i_iminus1_dict)
            # OMMIT BUILDING PEPTIDE FOR THE CHAINS THAT ALREADY EXIST IN tmp_peptide_folder/
            for chainIndex, chainScore in enumerate(all_chainScore_list):
                if not chainIndex in chainIndices2skip_set:
                    remaining_chainIndex_list.append(chainIndex)
                    remaining_chainScore_list.append(chainScore)
            
            results = list(futures.map(build_Peptide_Tree, remaining_chainIndex_list, remaining_chainScore_list))   # build peptide tree from chain
            # results will be empty!!
            #for chainIndex, chainScore in zip(remaining_chainIndex_list, remaining_chainScore_list):
            #    build_Peptide_Tree(remaining_chainIndex_list, remaining_chainScore_list)
            
            ##
            ## Write peptides to files
            ##
            # LOAD ALL WRITEN PEPTIDE FILES
            progbar = ProgressBar(100)
            peptide_file_num = 0
            #while peptide_file_num < total_chain_number:    # CREATE MISSING PEPTIDE TREES
            fnames=os.listdir('tmp_peptide_folder/')
            fpattern = re.compile('chainIndex_[0-9]+.pickle.bz2')
            peptide_files_list = list(filter(fpattern.search, fnames))
            peptide_file_num = len(peptide_files_list)
            for chainIndex, chainScore in enumerate(all_chainScore_list):
                if not 'chainIndex_'+str(chainIndex)+'.pickle.bz2' in peptide_files_list:
                    build_Peptide_Tree(chainIndex, chainScore)
            
            print("Saving all the peptide sequences to peptides.list file ...")
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
            # results will be an empty list
            # Now concatenate all parts of peptides.list and peptides.fasta
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
                            #f.write("Chain probability = "+str(chain[-1])+" peptide P_"+str(len(peptide))+"_"+str(index+1) + " = "+''.join(peptide)+" "+'-'.join(chain[0:-1])+" "+','.join(map(str, peptideProbList[:-1]))+"\n")
                            #fasta_filehandler.write(">P_"+str(len(peptide))+"_"+str(index+1)+"\n"+''.join(peptide)+"\n")
                            f.write(line_part1 + peptide_name + " = " + peptide_seq + line_part2 + "\n")
                            fasta_filehandler.write(">" + peptide_name + "\n" + peptide_seq + "\n")
                            index += 1
                os.remove("peptides."+str(args.MAX_PEPTIDE_LENGTH)+"mers.list.chunk"+part)
            shutil.rmtree('tmp_peptide_folder/', ignore_errors=True)
    
    except:
        type, value, tb = sys.exc_info()
        lines = traceback.format_exception(type, value, tb)
        print((''.join(lines)))
        raise