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

import sys, re, os, pickle, traceback, shutil, bz2, math
from scoop import futures, shared
import numpy as np
from operator import itemgetter
from collections import OrderedDict
from ete3 import Tree
import ftplib
from argparse import ArgumentParser
from scipy.stats.mstats import zscore
from scipy import stats, sqrt
import collections
import gc
from cluster import HierarchicalClustering
from .global_vars import *
from .global_func import *
from .csa_for_tocsy import *
from .probhist import *

################################# AUTOMATIC ASSIGNMENT OF AMINO ACID TYPES TO TOCSY RESONANCES FUNCTION DEFINITIONS ####################################################

def get_aatypes_from_H_C_resonpair_2Dhist(Hreson,
                                          Creson,
                                          TOCSY_reson_index,
                                          histload):
    """
        FUNCTION to find all possible aa types and C-H types and the respective probabilities from a pair of aliphatic H,C resonances.
        If the 2D-hist probability is 0.0 then the 1D-hist carbon probability is used instead.
        ARGUMENTS:
        TOCSY_reson_index:  an index indicating the position of this H-C resonance pair in the TOCSY group it belongs to
        RETURNS:
        matches_list: list of lists of the form [aa type, Carbon name, Hydrogen name, C-H 2D-hist probability,
                                                TOCSY_reson_index, H resonance, C resonance, carbon group]
    """
    
    #print "DEBUG: Hreson=",Hreson,"Creson=",Creson
    matches_list = [] # list of tuples of the form (aa type, Carbon name, Hydrogen name, overall probability, )
    histcalc = ProbHist_Calculator()
    for aa in list(histload.aa_CHpair_binProbabilityList_mdict.keys()):  # iterate over all amino acids
        carbon_matches_list = []    # list of tuples of the form (aa type,carbon,probability); here aa should be the same in all the tuples of the list!
        
        for carbon in list(histload.aa_carbon_binDensityList_mdict[aa].keys()):    # iterate over all carbons of the current amino acid
            #print "DEBUG: carbon=", carbon
            probability = histcalc.get_C_probability_from_1Dhistogram(aa, carbon, Creson, histload)
            #print "DEBUG: aa=",aa,"carbon=",carbon,"probability=",probability
            if probability > 0.0:
                carbon_matches_list.append((aa,carbon,probability))
        
        for CH_pair in list(histload.aa_CHpair_binProbabilityList_mdict[aa].keys()):
            probability = histcalc.get_CH_probability_from_2Dhistogram(CH_pair, aa, Hreson, Creson, histload)
            # Keep only the amino acids for which were found both the carbon and the respective covalently bonded hydrogen
            # The probability of each prediction is given by the 2D histogram
            if probability > 0.0:
                Cname, Hname = CH_pair.split("-")
                matches_list.append([aa, Cname, Hname, probability, TOCSY_reson_index, Hreson, Creson, None])  # add also the index of this C-H pair in the TOCSY group it
            elif probability == 0.0:    # if the probability of the 2D hist at this point is 0,
                                        # use only the Carbon 1D histogram to get the probability
                Cname, Hname = CH_pair.split("-")
                #print "DEBUG: Cname=", Cname, "carbon_matches_list=", carbon_matches_list
                carbon_match = [match for match in carbon_matches_list if match[1]==Cname]
                if len(carbon_match) == 0:  # if the hist(C) of this Cname was zero, don't save it 
                    continue
                #print "DEBUG: len(carbon_match) is ", len(carbon_match)," and it should be 1 !"
                #print "DEBUG: carbon_match=", carbon_match
                weighted_average_probability = -1 * carbon_match[0][2] #  but first make it negative to distiguish it from weighted average probabilities
                #print "DEBUG: only the Carbon 1D histogram used:", [aa, Cname, Hname, weighted_average_probability, TOCSY_reson_index, Hreson, Creson, None]
                matches_list.append([aa, Cname, Hname, weighted_average_probability, TOCSY_reson_index, Hreson, Creson, None]) 
            
    #print "DEBUG: returning matches_list=", matches_list
    return matches_list


def get_aatypes_from_H_C_resonpair(Hreson,
                                   Creson,
                                   TOCSY_reson_index,
                                   hist,
                                   args):
    """
        FUNCTION to find all possible amino acid types and the respective probabilities from a pair of aliphatic H,C resonances.
        ARGUMENTS:
        TOCSY_reson_index:  an index indicating the position of this H-C resonance pair in the TOCSY group it belongs to
        RETURNS:
        valid_matches_list: list of lists of the form [aa type, Carbon name, Hydrogen name, weighted average probability, TOCSY_reson_index, H resonance, C resonance, carbon group]
    """
    global aa_carbonBondedHydrogensDict_dict
    
    def is_valid_hydrogen(aa_type,
                          Cname,
                          Hname,
                          C_bondedH_dict):
        """
            FUNCTION to check whether a given hydrogen is covalently bonded to the given carbon of a given aa type
            ARGUMENTS:
            C_bondedH_dict:     dictonary with keys the carbon names of a given amino acid and values the respective covalently hydrogen names
        """
        
        if Hname in C_bondedH_dict[Cname]:
            return True
        
        return False
    
    # print "DEBUG: Hreson=",Hreson,"Creson=",Creson
    hydrogen_matches_list = []   # list of tuples of the form (aa type,hydrogen,probability)
    carbon_matches_list = []    # list of tuples of the form (aa type,carbon,probability)
    for aa in list(hist.aa_hydrogen_binDensityList_mdict.keys()):  # iterate over all amino acids
        
        for hydrogen in list(hist.aa_hydrogen_binDensityList_mdict[aa].keys()):    # iterate over all hydrogens of the current amino acid
            bin_array, density_array = hist.aa_hydrogen_binDensityList_mdict[aa][hydrogen]
            probability = hist.get_probability_from_histogram(Hreson, bin_array, density_array)
            if probability > 0.0:
                hydrogen_matches_list.append((aa,hydrogen,probability))
            elif probability == 0.0:    # consider all cases where hist(H)=0
                hydrogen_matches_list.append((aa,hydrogen, 0.0))
        
        for carbon in list(hist.aa_carbon_binDensityList_mdict[aa].keys()):    # iterate over all carbons of the current amino acid
            bin_array, density_array = hist.aa_carbon_binDensityList_mdict[aa][carbon]
            probability = hist.get_probability_from_histogram(Creson, bin_array, density_array)
            #print "DEBUG: aa=",aa,"carbon=",carbon,"probability=",probability
            if probability > 0.0:
                carbon_matches_list.append((aa,carbon,probability))
    
    # Keep only the amino acids for which were found both the carbon and the respective covalently bonded hydrogen
    # The probability of each prediction will be (Carbon_probability + Hydrogen_probability) /2
    matches_list = [] # list of tuples of the form (aa type, Carbon name, Hydrogen name, overall probability, )
    for hydrogen_match in hydrogen_matches_list:    # iterate over all hydrogen matches found in the previous step
        for carbon_match in carbon_matches_list:    # iterate over all carbons matches found in the previous step
            if hydrogen_match[0] == carbon_match[0]:    # if they belong to the same amino acid
                aa_type = hydrogen_match[0]
                Hname = hydrogen_match[1]
                Cname = carbon_match[1]
                if hydrogen_match[2] > 0.0:
                    if args.PROBABILITY_MODEL == 1: # treat the H and C probabilities as dependent events
                        # print "DEBUG: treating the H and C probabilities as dependent events."
                        if aa == "LEU":
                            weighted_average_probability = (1.0 * hydrogen_match[2] + 0.1 * carbon_match[2])/float((1.0 + 0.1))
                        else:
                            weighted_average_probability = (args.H_weight * hydrogen_match[2] + args.C_weight * carbon_match[2])/float((args.H_weight + args.C_weight))
                    elif args.PROBABILITY_MODEL == 2:
                        # print  "DEBUG: treating the H and C probabilities as independent events."
                        weighted_average_probability = hydrogen_match[2] * carbon_match[2]
                elif hydrogen_match[2] == 0.0:  # if hist(H)=0, consider only hist(C) but make it negative to distiguish it from weighted average probabilities
                    weighted_average_probability = -1 * carbon_match[2]
                matches_list.append((aa_type, Cname, Hname, weighted_average_probability))
    
    # TODO: calculate Z-scores from the histogram probabilities of the remaining unmatched carbons and ignore the hydrogens for those C with Z-score > 0.
    valid_matches_list = []    # the same list but contains only the most probable C-H matching pair from each aminoacid. 
    previous_aatype = None
    for quartet in matches_list:
        #print "DEBUG: quartet=",quartet
        aa_type = quartet[0]
        Cname = quartet[1]
        Hname = quartet[2]
        if is_valid_hydrogen(aa_type, Cname, Hname, aa_carbonBondedHydrogensDict_dict[aa_type]) == False: # if this carbon is not covalently bonded to this hydrogen, skip it
            continue
        #OBSOLETE # if aa_type == previous_aatype:  # elegant way to keep only the (aatype, Cname, Hname) with the highest probability, e.g. we may have:
        #OBSOLETE # ('VAL', 'CB', 'HB', 1.4736221632773356e-05), ('VAL', 'CG2', 'HG2', 0.01135282415016002), ('VAL', 'CG1', 'HG1', 0.01433564978569348) but we keep only the 3rd
        #OBSOLETE # othewise we would keep 'CG2_HG2', 'CG1_HG1', 'CB_HB' for just one C-H resonance pair!
        #OBSOLETE #    continue    
        valid_matches_list.append([aa_type, Cname, Hname, quartet[3], TOCSY_reson_index, Hreson, Creson, None])  # add also the index of this H-C pair in the TOCSY group it
        # belongs to, but also the H and C resonances for future usage
        previous_aatype = aa_type
    
    #print "DEBUG: matches_list=", matches_list
    #print "DEBUG: valid_matches_list=", valid_matches_list
    return valid_matches_list


def select_correct_H_C_resonpairs(aatypeResonpairMatchesTuple_list):
    """
        FUNCTION  to match the TOCSY resonance pairs to as much as possible C-H resonance pairs of a PARTICULAR aatype. I may need to used a Genetic Algorith to do that
        if the current implementation proves to be error prone.
        
        ARGUMETS:
        aatypeResonpairMatchesTuple_list:   list of lists of the form [aa_type, C, H, weighted average probability, selected TOCSY reson index, H resonance, C resonance, clusterID] where aa_type is the same in the whole list (e.g. "PRO")
        
        RETURNS:
        correct_C_H_resonpair_matches_list: list of list of the form [aa_type, C, H, weighted average probability, selected TOCSY reson index], where aa_type is the same in the whole list (e.g. "PRO"). This list contains
                                            the final, uniq nucleus type assignment for the current amino acid. E.g.
        [['MET', 'CA', 'HA', 0.01993300618212495, 2, 3.843, 55.13, 1], ['MET', 'CB', 'HB2', 0.012634347058708313, 4, 1.266, 31.504, 2],
        ['MET', 'CG', 'HG2', 0.02281033340649995, 1, 1.508, 31.481, 2], ['MET', 'CB', 'HB3', 0.009737867955406978, 3, 1.911, 31.223, 3],
        ['MET', 'CG', 'HG3', 0.01607381733870664, 5, 0.403, 31.186, 3]]
    """
    global aa_carbonBondedHydrogensDict_dict
    
    def do_carbons_match(match, correct_C_H_resonpair_matches_list):
        """
            FUNCTION to check if the carbon resonances of geminal protons match. E.g. If he have added 
        """
        Cname = match[1]
        Creson = match[6]
        for correct_match in correct_C_H_resonpair_matches_list:
            # if the carbon name is already in the correct matches but its resonance is different, return false
            #print "DEBUG: Cname=",Cname, "correct_match[1]=", correct_match[1], "Creson=", Creson, "correct_match[6]=", correct_match[6]
            if Cname == correct_match[1] and approx_equal(Creson, correct_match[6], 0.3) == False:
                #print "DEBUG: conflict found, do_carbons_match() returning False!"
                return False
        
        return True # otherwise return true
        
        
    #print "DEBUG: entered function select_correct_H_C_resonpairs() with aatypeResonpairMatchesTuple_list=",aatypeResonpairMatchesTuple_list
    aa_type = aatypeResonpairMatchesTuple_list[0][0]    # recall that we have only one aa type in this script
    correct_C_H_resonpair_matches_list = [] # list of lists, see the Function definition
    used_TOCSY_reson_indices_list = []  # list of assigned peaks (TOCSY reson indices)
    used_C_H_pairs_list = []    # list of assigned C and H nucleus types
    
    # FIRST TAKE CARE OF CARBON GROUPING AND THEIR PROTONS (e.g. CB, CB2, CB3)
    clustIDs_list = [x[7] for x in aatypeResonpairMatchesTuple_list]    # list of cluster ID of each possible assignment
    clustIDs_set = set(clustIDs_list)
    # FIRST OF ALL SORT THE CLUSTERS ACCORDING TO THE TOTAL PROBABILITY OF THEIR CARBONS
    clustID_highestTotProb_dict = {}
    for clustID in clustIDs_set:
        TresonIndex_set = set([x[4] for x in aatypeResonpairMatchesTuple_list if x[7]==clustID])    # set with all different peaks (TOCSY reson indices) of the current cluster ID
        # NOW USE THE Treson indices WITH MULTIPLE TYPES OF CARBONS
        #print "DEBUG: clustIDs_list.count(clustID) =", clustIDs_list.count(clustID), " len(TresonIndex_set)=", len(TresonIndex_set)
        #if clustIDs_list.count(clustID) > 1 and len(TresonIndex_set) > 1:    # if there are 2 or more possible assignments and at least 2 different peaks within the same cluster ID
        nucleusAssignmentsforCurrentClustID_list = [x for x in aatypeResonpairMatchesTuple_list if x[7]==clustID]    # save all the members of this cluster ID
        # Recall that nucleusAssignmentsforCurrentClustID_list is already sorted by lines 1,2,3, BUT WE MUST ...
        nucleusAssignmentsforCurrentClustID_list.sort(key=itemgetter(3), reverse=True)  # ... SORT BY THE PREDICTION PROBABILITY
        
        # Make a list of the carbon types in each Treson index of this cluster and its maximum probability
        C_totalProb_dict = {}   # the total probability of each C type in this cluster. Use this to decide the assignment
        carbon_set = set([x[1] for x in nucleusAssignmentsforCurrentClustID_list])
        TresonIndex_set = set([x[4] for x in nucleusAssignmentsforCurrentClustID_list])
        for TRI in TresonIndex_set:
            nucleusAssignmentsforTRI_list = [x for x in nucleusAssignmentsforCurrentClustID_list if x[4]]
            nucleusAssignmentsforTRI_list.sort(key=itemgetter(3), reverse=True)  # sort by the prediction probability
            C_prob_dict = {}
            for prediction_list in nucleusAssignmentsforTRI_list:
                C = prediction_list[1]
                if C not in list(C_prob_dict.keys()):
                    C_prob_dict[C] = prediction_list[3]
            #print "DEBUG: TRI = ", TRI
            for C in carbon_set:
                if C not in list(C_prob_dict.keys()): # if this C type was not in possible predictions if this TRI, set its total probability to 0 (it will be completely excluded!)
                    C_totalProb_dict[C] = 0
                    continue
                try:
                    C_totalProb_dict[C] *= C_prob_dict[C]
                except KeyError:
                    C_totalProb_dict[C] = C_prob_dict[C]
        
        #print "DEBUG: C_totalProb_dict=", C_totalProb_dict
        sorted_C_totalProb_list = sorted(list(C_totalProb_dict.items()), key=itemgetter(1), reverse=True)   # sorted C_totalProb_list, e.g. [(CB, 0.9), (CG, 0.7)]
        
        clustID_highestTotProb_dict[clustID] = sorted_C_totalProb_list[0][1]    # save the total probability of the first element only, since the list was sorted
    # AND THIS IS THE LIST OF CLUSTER IDS SORTED BY THEIR HIGHEST TOTAL PROBABILITY
    #print "DEBUG: clustID_highestTotProb_dict = ", clustID_highestTotProb_dict
    sorted_clustIDs_list = [x[0] for x in sorted(list(clustID_highestTotProb_dict.items()), key=itemgetter(1), reverse=True)]
    #print "DEBUG: sorted_clustIDs_list=", sorted_clustIDs_list
    for clustID in clustIDs_set:
        if clustID not in sorted_clustIDs_list:
            sorted_clustIDs_list.append(clustID)
    #print "DEBUG: after addition of missing clusters sorted_clustIDs_list=", sorted_clustIDs_list
    
    
    CONTINUE_1ST_GROUPING = True # if no new assignments are saved, set to False to stop the iterations
    while CONTINUE_1ST_GROUPING:
        CONTINUE_1ST_GROUPING = False
        for clustID in sorted_clustIDs_list:
            TresonIndex_set = set([x[4] for x in aatypeResonpairMatchesTuple_list if x[7]==clustID])    # set with all different peaks (TOCSY reson indices) of the current cluster ID
            if clustIDs_list.count(clustID) > 1 and len(TresonIndex_set) > 1:    # if there are 2 or more possible assignments and at least 2 different peaks within the same cluster ID
                nucleusAssignmentsforCurrentClustID_list = [x for x in aatypeResonpairMatchesTuple_list if x[7]==clustID]    # save all the members of this cluster ID
                # Recall that nucleusAssignmentsforCurrentClustID_list is already sorted by lines 1,2,3
                # BUT FIRST CHECK IF ANY OF THE Treson INDICES CONTAINS A SINGLE TYPE OF CARBON. IN THAT CASE, THAT CARBON HAS PRIORITY IN THE ASSIGNMENT.
                for TresonIndex in TresonIndex_set:
                    carbon_list = [assignment[1] for assignment in nucleusAssignmentsforCurrentClustID_list if assignment[4]==TresonIndex]  # save all types of C in this Treson index
                    carbon_set = set(carbon_list)
                    if len(carbon_set) == 1:    # CONDITION: if there is only a single type of carbon within the Treson index 
                        assignment_list = [assignment for assignment in nucleusAssignmentsforCurrentClustID_list if assignment[4]==TresonIndex]  # save all types of C in this Treson index
                        assignment_list.sort(key=itemgetter(3), reverse=True)  # sort by the prediction probability
                        assignment = assignment_list[0] # USE ONLY THE ASSIGNMENT WITH THE HIGHEST PREDICTION PROBABILITY
                        C = assignment[1]
                        CARBON_OF_THIS_CLUSTER = C  # all assignments of this clustID must contain the this type of Carbon!!!
                        H = assignment[2]
                        # CONDITIONS:
                        # 1) TOCSY reson index has not been used again
                        # 2) this H must be covalently bonded to this C
                        # 3) this C,H pair has not been used again
                        if not TresonIndex in used_TOCSY_reson_indices_list and H in aa_carbonBondedHydrogensDict_dict[aa_type][C] and not C+"_"+H in used_C_H_pairs_list:
                            #print "1st GROUPING EXAMPLE: nucleusAssignmentsforCurrentClustID_list=", nucleusAssignmentsforCurrentClustID_list
                            correct_C_H_resonpair_matches_list.append(assignment)
                            used_TOCSY_reson_indices_list.append(TresonIndex)
                            used_C_H_pairs_list.append(C+"_"+H)
                            CONTINUE_1ST_GROUPING = True # since a new assignment has been saved, iterate the cycle once more
                            # CLEAN aatypeResonpairMatchesTuple_list from other possible assignments of the same clustID that do not contain this C type
                            for x in aatypeResonpairMatchesTuple_list:
                                # CONDITIONS:
                                # 1) the possible assignment must belong to the current clusterID
                                # 2) if the Carbon type in the assignment is not the CARBON_OF_THIS_CLUSTER, remove the assignment
                                if x[7] == clustID and x[1] != CARBON_OF_THIS_CLUSTER:
                                    aatypeResonpairMatchesTuple_list.remove(x)
        #print "DEBUG: after 1st grouping correct_C_H_resonpair_matches_list=", correct_C_H_resonpair_matches_list
    
    for clustID in sorted_clustIDs_list:
        TresonIndex_set = set([x[4] for x in aatypeResonpairMatchesTuple_list if x[7]==clustID])    # set with all different peaks (TOCSY reson indices) of the current cluster ID
        # NOW USE THE Treson indices WITH MULTIPLE TYPES OF CARBONS
        if clustIDs_list.count(clustID) > 1 and len(TresonIndex_set) > 1:    # if there are 2 or more possible assignments and at least 2 different peaks within the same cluster ID
            nucleusAssignmentsforCurrentClustID_list = [x for x in aatypeResonpairMatchesTuple_list if x[7]==clustID]    # save all the members of this cluster ID
            # Recall that nucleusAssignmentsforCurrentClustID_list is already sorted by lines 1,2,3, BUT WE MUST ...
            nucleusAssignmentsforCurrentClustID_list.sort(key=itemgetter(3), reverse=True)  # ... SORT BY THE PREDICTION PROBABILITY
            
            # Make a list of the carbon types in each Treson index of this cluster and its maximum probability
            C_totalProb_dict = {}   # the total probability of each C type in this cluster. Use this to decide the assignment
            carbon_set = set([x[1] for x in nucleusAssignmentsforCurrentClustID_list])
            TresonIndex_set = set([x[4] for x in nucleusAssignmentsforCurrentClustID_list])
            for TRI in TresonIndex_set:
                nucleusAssignmentsforTRI_list = [x for x in nucleusAssignmentsforCurrentClustID_list if x[4]]
                nucleusAssignmentsforTRI_list.sort(key=itemgetter(3), reverse=True)  # sort by the prediction probability
                C_prob_dict = {}
                for prediction_list in nucleusAssignmentsforTRI_list:
                    C = prediction_list[1]
                    if C not in list(C_prob_dict.keys()):
                        C_prob_dict[C] = prediction_list[3]
                
                for C in carbon_set:
                    if C not in list(C_prob_dict.keys()): # if this C type was not in possible predictions if this TRI, set its total probability to 0 (it will be completely excluded!)
                        C_totalProb_dict[C] = 0
                        continue
                    try:
                        C_totalProb_dict[C] *= C_prob_dict[C]
                    except KeyError:
                        C_totalProb_dict[C] = C_prob_dict[C]
            
            sorted_C_totalProb_list = sorted(list(C_totalProb_dict.items()), key=itemgetter(1), reverse=True)   # sorted C_totalProb_list, e.g. [(CB, 0.9), (CG, 0.7)]
            #print "DEBUG: clustID=", clustID, "sorted_C_totalProb_list=", sorted_C_totalProb_list
            for C_prob in sorted_C_totalProb_list:
                CARBON_OF_THIS_CLUSTER = C_prob[0] # all assignments of this clustID must contain this type of Carbon!!!
                
                for assignment in nucleusAssignmentsforCurrentClustID_list: # assignmend are sorted by probability; assignment example: ['MET', 'CG', 'HG3', 0.02342287175549055, 4, 1.266, 31.504, 2]
                    C = assignment[1]
                    if C != CARBON_OF_THIS_CLUSTER: # if this possible assignment contains a different C type, skip it
                        continue
                    H = assignment[2]
                    TresonIndex = assignment[4]
                    # CONDITIONS:
                    # 1) TOCSY reson index has not been used again
                    # 2) this H must be covalently bonded to this C
                    # 3) this C,H pair has not been used again
                    if not TresonIndex in used_TOCSY_reson_indices_list and H in aa_carbonBondedHydrogensDict_dict[aa_type][C] and not C+"_"+H in used_C_H_pairs_list:
                        #print "2nd GROUPING EXAMPLE: nucleusAssignmentsforCurrentClustID_list=", nucleusAssignmentsforCurrentClustID_list
                        correct_C_H_resonpair_matches_list.append(assignment)
                        used_TOCSY_reson_indices_list.append(TresonIndex)
                        used_C_H_pairs_list.append(C+"_"+H)
                        # CLEAN aatypeResonpairMatchesTuple_list from other possible assignments of the same clustID that do not contain this C type
                        for x in aatypeResonpairMatchesTuple_list:
                            # CONDITIONS:
                            # 1) the possible assignment must belong to the current clusterID
                            # 2) if the Carbon type in the assignment is not the CARBON_OF_THIS_CLUSTER, remove the assignment
                            if x[7] == clustID and x[1] != CARBON_OF_THIS_CLUSTER:
                                aatypeResonpairMatchesTuple_list.remove(x)
        #print "DEBUG: after 2nd grouping correct_C_H_resonpair_matches_list=", correct_C_H_resonpair_matches_list
    
    #print "DEBUG: after grouping all clustIDs correct_C_H_resonpair_matches_list=", correct_C_H_resonpair_matches_list
    
    
    # NOW TAKE CARE OF UNGROUPED CARBONS, SPECIFICALLY CLUSTERS WITH ONLY ONE Treson index AND ONLY ONE POSSIBLE ASSIGNMENT
    for clustID in sorted_clustIDs_list:
        TresonIndex_set = set([x[4] for x in aatypeResonpairMatchesTuple_list if x[7]==clustID])    # set with all different peaks (TOCSY reson indices) of the current cluster ID
        # NOW USE THE Treson indices WITH MULTIPLE TYPES OF CARBONS
        if clustIDs_list.count(clustID) == 1 and len(TresonIndex_set) == 1:    # if there is only 1 possible assignment and only 1 peak within the same cluster ID
            nucleusAssignmentsforCurrentClustID_list = [x for x in aatypeResonpairMatchesTuple_list if x[7]==clustID]    # save the single member of this cluster ID
            # SINCE IT IS ONLY ONE POSSIBLE ASSIGNMENT, WE DON'T NEED TO SORT nucleusAssignmentsforCurrentClustID_list
            assignment = nucleusAssignmentsforCurrentClustID_list[0]
            C = assignment[1]
            CARBON_OF_THIS_CLUSTER = C  # all assignments of this clustID must contain the this type of Carbon!!!
            H = assignment[2]
            TresonIndex = assignment[4]
            # CONDITIONS:
            # 1) TOCSY reson index has not been used again
            # 2) this H must be covalently bonded to this C
            # 3) this C,H pair has not been used again
            if not TresonIndex in used_TOCSY_reson_indices_list and H in aa_carbonBondedHydrogensDict_dict[aa_type][C] and not C+"_"+H in used_C_H_pairs_list:
                #print "3rd GROUPING EXAMPLE: nucleusAssignmentsforCurrentClustID_list=", nucleusAssignmentsforCurrentClustID_list
                correct_C_H_resonpair_matches_list.append(assignment)
                used_TOCSY_reson_indices_list.append(TresonIndex)
                used_C_H_pairs_list.append(C+"_"+H)
                # CLEAN aatypeResonpairMatchesTuple_list from other possible assignments of the same clustID that do not contain this C type
                for x in aatypeResonpairMatchesTuple_list:
                    # CONDITIONS:
                    # 1) the possible assignment must belong to the current clusterID
                    # 2) if the Carbon type in the assignment is not the CARBON_OF_THIS_CLUSTER, remove the assignment
                    if x[7] == clustID and x[1] != CARBON_OF_THIS_CLUSTER:
                        aatypeResonpairMatchesTuple_list.remove(x)
    #print "DEBUG: after clusters with a single TRI and single assignment correct_C_H_resonpair_matches_list=", correct_C_H_resonpair_matches_list
    
    # THEN, TAKE CARE OF UNGROUPED CARBONS, SPECIFICALLY CLUSTERS WITH ONLY ONE Treson index BUT MULTIPLE POSSIBLE ASSIGNMENTS
    for clustID in sorted_clustIDs_list:
        TresonIndex_set = set([x[4] for x in aatypeResonpairMatchesTuple_list if x[7]==clustID])    # set with all different peaks (TOCSY reson indices) of the current cluster ID
        # NOW USE THE Treson indices WITH MULTIPLE TYPES OF CARBONS
        if clustIDs_list.count(clustID) > 1 and len(TresonIndex_set) == 1:    # if there is only 1 possible assignment and only 1 peak within the same cluster ID
            nucleusAssignmentsforCurrentClustID_list = [x for x in aatypeResonpairMatchesTuple_list if x[7]==clustID]    # save the single member of this cluster ID
            TresonIndex = TresonIndex_set.pop() # substract the single TresonIndex of this clustID
            assignment_list = [assignment for assignment in nucleusAssignmentsforCurrentClustID_list if assignment[4]==TresonIndex]  # save all types of C in this Treson index
            assignment_list.sort(key=itemgetter(3), reverse=True)  # SORT BY THE PREDICTION PROBABILITY
            assignment = assignment_list[0]                        # USE ONLY THE ASSIGNMENT WITH THE HIGHEST PREDICTION PROBABILITY
            C = assignment[1]
            CARBON_OF_THIS_CLUSTER = C  # all assignments of this clustID must contain the this type of Carbon!!!
            H = assignment[2]
            TresonIndex = assignment[4]
            # CONDITIONS:
            # 1) TOCSY reson index has not been used again
            # 2) this H must be covalently bonded to this C
            # 3) this C,H pair has not been used again
            if not TresonIndex in used_TOCSY_reson_indices_list and H in aa_carbonBondedHydrogensDict_dict[aa_type][C] and not C+"_"+H in used_C_H_pairs_list:
                #print "4th GROUPING EXAMPLE: nucleusAssignmentsforCurrentClustID_list=", nucleusAssignmentsforCurrentClustID_list
                correct_C_H_resonpair_matches_list.append(assignment)
                used_TOCSY_reson_indices_list.append(TresonIndex)
                used_C_H_pairs_list.append(C+"_"+H)
                # CLEAN aatypeResonpairMatchesTuple_list from other possible assignments of the same clustID that do not contain this C type
                for x in aatypeResonpairMatchesTuple_list:
                    # CONDITIONS:
                    # 1) the possible assignment must belong to the current clusterID
                    # 2) if the Carbon type in the assignment is not the CARBON_OF_THIS_CLUSTER, remove the assignment
                    if x[7] == clustID and x[1] != CARBON_OF_THIS_CLUSTER:
                        aatypeResonpairMatchesTuple_list.remove(x)
    #print "DEBUG: after clusters with a single TRI but multiple assignments correct_C_H_resonpair_matches_list=", correct_C_H_resonpair_matches_list
    
    # FINALLY TAKE CARE OF WHATEVER IS LEFT IN aatypeResonpairMatchesTuple_list. DON'T KNOW IF THIS PART IS NECESSARY ANY MORE!
    C_H_occupancy_dict = {} # dict with keys C_H pairs and values the number of TOCSY peaks of the currect group this C_H pair matched
    C_H_matches_dict = {}
    for C, H_list in list(aa_carbonBondedHydrogensDict_dict[aa_type].items()):
        for H in H_list:
            #print "DEBUG:"+C+"_"+H
            if C+"_"+H in used_C_H_pairs_list:  # if this C,H pair has been used again, skip it
                continue
            matchesTuples_list = [t for t in aatypeResonpairMatchesTuple_list if t[1]==C and t[2]==H]    # keep all the tuples matching the current C-H pair
            C_H_occupancy_dict[C+"_"+H] = len(matchesTuples_list)
            C_H_matches_dict[C+"_"+H] = matchesTuples_list
    
    #print "DEBUG: C_H_occupancy_dict = ", C_H_occupancy_dict
    #print "DEBUG: C_H_matches_dict = ", C_H_matches_dict
    sorted_C_H_occupancyTuple_list = sorted(list(C_H_occupancy_dict.items()), key=itemgetter(1))  # the C_H_occupancyTuple tuples sorted by occupancy, so ambiguous assignments (higher occupancy) are left for the end
    #print "DEBUG: sorted_C_H_occupancyTuple_list=", sorted_C_H_occupancyTuple_list
    for C_H_occupancy_tuple in sorted_C_H_occupancyTuple_list:
        C_H_pair = C_H_occupancy_tuple[0]
        C = C_H_pair.split("_")[0]
        H = C_H_pair.split("_")[1]
        for match in C_H_matches_dict[C_H_pair]:   # we don't need to sort here the alternative C-H resonances according to the weighted average prob,
                                # because it has been done already at "possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list.sort(key=itemgetter(0,1,2,3), reverse=True)"
            #print "DEBUG: match = ", match
            #print "DEBUG: correct_C_H_resonpair_matches_list=", correct_C_H_resonpair_matches_list
            TOCSY_reson_index = match[4]
            # CONDITIONS:
            # 1) This Treson index has not been assigned yet
            # 2) This condition was useful before the incorporation of carbon grouping. I.e. it checked whether HB2-CB and HB3-CB had the same carbon resonance (+- 0.3)
            if not TOCSY_reson_index in used_TOCSY_reson_indices_list and do_carbons_match(match, correct_C_H_resonpair_matches_list):
                #print "5th GROUPING EXAMPLE: nucleusAssignmentsforCurrentClustID_list=", nucleusAssignmentsforCurrentClustID_list
                correct_C_H_resonpair_matches_list.append(match)
                used_TOCSY_reson_indices_list.append(TOCSY_reson_index)
                break   # move to the next C-H resonance pair
    
    #print "DEBUG: final correct_C_H_resonpair_matches_list=",correct_C_H_resonpair_matches_list
    return correct_C_H_resonpair_matches_list


def populate_leaves_for_correct_H_C_resonpairs2(Assignment_Tree,
                                                possible_aatype_prob_C_H_resonpair_TOCSYindex_assignID_list_list,
                                                TRI):
    """
        FUNCTION that adds new branches to the leaves of the Tree.
        ARGUMENTS:
        Assignment_Tree:    The Tree structure with connectivities
        RETURNS:
        (Assignment_Tree, BOOLEAN):    A tuple with elements the input Tree structure with new branches (if applicable), and a BOOLEAN value which is True if the function added
                                       new leaves to the Tree, or False otherwise
    """
    nucleusAssignmentsforCurrentTRI_list = [x for x in possible_aatype_prob_C_H_resonpair_TOCSYindex_assignID_list_list if x[4]==TRI] # save all the members of this cluster ID
    #TresonIndex_set = set([x[4] for x in possible_aatype_prob_C_H_resonpair_TOCSYindex_assignID_list_list if x[7]==clustID])    # set with all different peaks (TOCSY reson indices) of the current cluster ID
    assignments_list = []
    #if len(TresonIndex_set) == 1:
    assignments_list = nucleusAssignmentsforCurrentTRI_list
    number_of_new_leaves = 0
    # ATTENTION: never use Assignment_Tree.iter_leaf_names(), it doesn't return the names in the order
    # ATTENTION: corresponding to Assignment_Tree.get_leaves()!!!
    for leaf in Assignment_Tree.get_leaves():
        try:
            for assignment in assignments_list:
                _Ctype = assignment[1]
                _Htype = assignment[2]
                _Prob = assignment[3]
                _TRI = assignment[4]
                _H = assignment[5]
                _C = assignment[6]
                _clustID = assignment[7]
                assignID = assignment[8]
                ancestors_list = [ancestor.Ctype + "_" + ancestor.Htype for ancestor in leaf.get_ancestors()]
                if _Ctype + "_" + _Htype in ancestors_list:   # if this C-H type is already in the Tree, skip it
                    continue
                new_child = leaf.add_child(name=assignID) # add a new brach to the current TOCSY add index (leaf) with length the respective probability and return it
                new_child.add_features(Ctype=_Ctype, Htype=_Htype, Prob=_Prob, TRI=_TRI, H=_H, C=_C, clustID=_clustID)
                number_of_new_leaves += 1
        except KeyError:
            continue
    
    # print Assignment_Tree.get_ascii(show_internal=True, compact=False, attributes=["name", "Ctype", "Htype", "TRI"])
    # print Assignment_Tree.get_ascii(show_internal=True, compact=False)
    if number_of_new_leaves > 0:
        return (Assignment_Tree, True)
    else:
        return (Assignment_Tree, False)


def get_combined_weighted_and_presence_probability(AA_type,
                                                   aatype_probsum_dict,
                                                   aatype_CHnucleiType_presenceProbSum_mdict,
                                                   args):
    """
    FUNCTION to calculate the probability of an amino acid type AA_type given an assignment, as in eq. 6 of RESCUE2 paper.
    But instead of multiplying all C & H probabilities, we multiply the weighted averages of C-H pair probabilities.
    
    ARGUMENTS:
    aatype_probsum_dict:    in this case it is a dictionary with the product of the weighted average probabilities of all C-H pairs assigned.
    aatype_CHnucleiType_presenceProbSum_mdict:  multidimensiona dict of the form aa type --> C_H pair --> weighted average presence probabilities as in eq. 6 of RESCUE2 paper.
                                                    It misses the absence probabilities, which must be computed afterwards.
    """
    global aa_carbonBondedHydrogensDict_dict, Prob_CS
    
    carbon_bondexHydrogensList_dict = aa_carbonBondedHydrogensDict_dict[AA_type]
    
    for Cname in list(carbon_bondexHydrogensList_dict.keys()):
        for Hname in carbon_bondexHydrogensList_dict[Cname]:
            #print "DEBUG: checking ", Cname+"_"+Hname, "of aa type", AA_type
            if not Cname+"_"+Hname in list(aatype_CHnucleiType_presenceProbSum_mdict[AA_type].keys()):
                # if there is no weighted presence probability for the current C-H pair, save the respective absence probability
                aatype_CHnucleiType_presenceProbSum_mdict[AA_type][Cname+"_"+Hname] = 1 - (args.H_weight * Prob_CS[AA_type][Hname] + args.C_weight * Prob_CS[AA_type][Cname])/float((args.H_weight + args.C_weight))
                #print "DEBUG: Assigning absence to ", Cname+"_"+Hname, ":1 - *(",args.H_weight,"*", Prob_CS[AA_type][Hname], "+", args.C_weight, "*", Prob_CS[AA_type][Cname], ")/float((args.H_weight + args.C_weight)) =", aatype_CHnucleiType_presenceProbSum_mdict[AA_type][Cname+"_"+Hname]
    
    #print "DEBUG: aatype_CHnucleiType_presenceProbSum_mdict[AA_type]=", aatype_CHnucleiType_presenceProbSum_mdict[AA_type].items()
    for weighted_presence_probabilty in list(aatype_CHnucleiType_presenceProbSum_mdict[AA_type].values()):
        try:
            weighted_presence_probabilities_product *= weighted_presence_probabilty
        except NameError:
            weighted_presence_probabilities_product = weighted_presence_probabilty
    
    #print "DEBUG: final probability=", weighted_presence_probabilities_product * aatype_probsum_dict[AA_type]
    return weighted_presence_probabilities_product * aatype_probsum_dict[AA_type]


def consider_only_Carbon_matches(possible_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list,
                                 only_missing=True,
                                 remove_Cbased_predictions=False):
    """
        FUNCTION to check if there are C-H resonance pairs in the currect TAAIG group for every aa type, and to allow C-H type
        prediction only from C resonances when one of the following conditions is met:
                1) The H resonance is negative.
                2) No aa type prediction can be made for a particular TAAIG group.
        
        ARGUMENTS:
        possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list:   list of lists of the form (aa type, Carbon name, Hydrogen name, weighted average probability, TOCSY_reson_index, H resonance, C resonance, carbon group)
                                                                        This list includes also negative probabilities wherever only the Carbon resonance match was considered.
        possible_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list:   list of tuples of the form (aa type, Carbon name, Hydrogen name, weighted average probability, TOCSY_reson_index, H resonance, C resonance)
        only_missing:   if True only the C hist probabilities will be considered in C-H resonance cases where no C-H type prediction could be made for a particular aa. If False,
                        then only the C hist probabilities will be considered in every C-H resonance case. E.g. If we have the following predictions:
                        ('THR', 'CB', 'HB', -0.007423257864062593, 1, 4.239, 61.662), ('THR', 'CA', 'HA', 0.022600159959100728, 1, 4.239, 61.662)
                        ('THR', 'CB', 'HB', -0.020221544866552674, 2, 1.634, 68.57), ('THR', 'CA', 'HA', -0.0014439440299790283, 2, 1.634, 68.57)
                        ('THR', 'CG2', 'HG2', 0.001264928291792422, 3, 0.836, 26.318)
                        If only_missing=True, the allowed predictions will become:
                        ('THR', 'CA', 'HA', 0.022600159959100728, 1, 4.239, 61.662)
                        ('THR', 'CB', 'HB', -0.020221544866552674, 2, 1.634, 68.57), ('THR', 'CA', 'HA', -0.0014439440299790283, 2, 1.634, 68.57)
                        ('THR', 'CG2', 'HG2', 0.001264928291792422, 3, 0.836, 26.318)
                        But if only_missing=False, they will remain the same:
                        ('THR', 'CB', 'HB', -0.007423257864062593, 1, 4.239, 61.662), ('THR', 'CA', 'HA', 0.022600159959100728, 1, 4.239, 61.662)
                        ('THR', 'CB', 'HB', -0.020221544866552674, 2, 1.634, 68.57), ('THR', 'CA', 'HA', -0.0014439440299790283, 2, 1.634, 68.57)
                        ('THR', 'CG2', 'HG2', 0.001264928291792422, 3, 0.836, 26.318)
    """
    #print "DEBUG: possible_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list=", possible_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list
    possible_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list.sort(key=itemgetter(0,4), reverse=True)  # sort by aa type and Treson index
    previous_aatype = None
    aatype_rest7_dict = {}  # dict of the form: aa type -> (Carbon name, Hydrogen name, weighted average probability, TOCSY_reson_index, H resonance, C resonance)
    for group8 in possible_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list:
        aa_type = group8[0]
        if aa_type != previous_aatype and previous_aatype != None:
            aatype_rest7_dict[aa_type] = [(group8[1:])]
        elif aa_type == previous_aatype or previous_aatype == None:
            try:
                aatype_rest7_dict[aa_type].append((group8[1:]))
            except KeyError:
                aatype_rest7_dict[aa_type] = [(group8[1:])]
            
        previous_aatype = aa_type    
    
    new_possible_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list = [] # possible_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list after filtering
    for aa_type in list(aatype_rest7_dict.keys()):
        # print("DEBUG: aa_type=", aa_type)
        previous_TresonIndex = None
        C_H_type_prediction_list = []   # list with the valid 
        for group7 in aatype_rest7_dict[aa_type]:
            # print("DEBUG: aa_type=", aa_type, "group7", group7)
            if remove_Cbased_predictions == True and group7[2] < 0:
                continue
            if only_missing == False:   # if instructed add both C- and C-H-based aa type predictions
                new_possible_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list.append( (aa_type, group7[0], group7[1], abs(group7[2]), group7[3], group7[4], group7[5], group7[6]) )
            
            elif only_missing == True:
                # Recall the TresonIndex numbers are sorted in aatype_rest7_dict as a result of sorting possible_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list
                # print("DEBUG: previous_TresonIndex=", previous_TresonIndex, "group7[3]=", group7[3])
                if previous_TresonIndex != group7[3] and previous_TresonIndex != None:
                    # do the appropriate changes and save the C-H type predictions
                    positive_C_H_type_prediction_list = [g7 for g7 in C_H_type_prediction_list if g7[3] > 0]    # list of the C-H type predictions with positive weighted average probabilities
                    # print("DEBUG: positive_C_H_type_prediction_list=", positive_C_H_type_prediction_list)
                    for C_H_type_prediction in C_H_type_prediction_list:
                        # ADD ALL POSITIVE PREDICTIONS AND THOSE NEGATIVE PREDICTIONS WITH PROBABILITY ABOVE THE AVERAGE PROBABILITY
                        # OF ALL POSITIVE PREDICTIONS, ONLY IF THE H RESONACNE IS WITHIN THE LOW OCCUPANCY REGIONS
                        if C_H_type_prediction[3] > 0 or ( (C_H_type_prediction[3] < 0 and is_H_in_low_occupancy_region(C_H_type_prediction[5]) ) and ( len(positive_C_H_type_prediction_list) == 0 or abs(C_H_type_prediction[3]) >= np.average([x[3] for x in positive_C_H_type_prediction_list]) ) ):
                            lst = list(C_H_type_prediction)
                            if (C_H_type_prediction[3] < 0 and is_H_in_low_occupancy_region(C_H_type_prediction[5])):
                                lst[5] = str(lst[5]) + "C"   # MARK THE H RESONANCE TO KNOW THAT IT WAS NOT USED FOR AA TYPE PREDICTION (marking C resonance will raise error)
                            lst[3] = abs(C_H_type_prediction[3])
                            C_H_type_prediction = tuple(lst)
                            new_possible_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list.append(C_H_type_prediction)
                    
                    # now reinitialize C_H_type_prediction_list
                    C_H_type_prediction_list = [ (aa_type, group7[0], group7[1], group7[2], group7[3], group7[4], group7[5], group7[6]) ]
                
                elif previous_TresonIndex == group7[3] or previous_TresonIndex == None:
                    C_H_type_prediction_list.append( (aa_type, group7[0], group7[1], group7[2], group7[3], group7[4], group7[5], group7[6]) )
        
            previous_TresonIndex = group7[3]
        
        if only_missing == True:  # NOW PROCESS THE LAST TresonIndex
            # do the appropriate changes and save the C-H type predictions
            positive_C_H_type_prediction_list = [g7 for g7 in C_H_type_prediction_list if g7[3] > 0]    # list of the C-H type predictions with positive weighted average probabilities
            #print "DEBUG: positive_C_H_type_prediction_list=", positive_C_H_type_prediction_list
            for C_H_type_prediction in C_H_type_prediction_list:
                # ADD ALL POSITIVE PREDICTIONS AND THOSE NEGATIVE PREDICTIONS WITH PROBABILITY ABOVE THE AVERAGE PROBABILITY OF ALL POSITIVE PREDICTIONS, ONLY IF THE H RESONACNE
                # IS NEGATIVE
                if C_H_type_prediction[3] > 0 or ( (C_H_type_prediction[3] < 0 and is_H_in_low_occupancy_region(C_H_type_prediction[5])) and ( len(positive_C_H_type_prediction_list) == 0 or abs(C_H_type_prediction[3]) >= np.average([x[3] for x in positive_C_H_type_prediction_list]) ) ):
                    lst = list(C_H_type_prediction)
                    if (C_H_type_prediction[3] < 0 and is_H_in_low_occupancy_region(C_H_type_prediction[5])):
                        lst[5] = str(lst[5]) + "C"   # MARK THE H RESONANCE TO KNOW THAT IT WAS NOT USED FOR AA TYPE PREDICTION
                                                     # (marking C resonance will raise error)
                    lst[3] = abs(C_H_type_prediction[3])
                    C_H_type_prediction = tuple(lst)
                    new_possible_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list.append(C_H_type_prediction)
        
    #print "DEBUG: new_possible_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list="
    #for t in new_possible_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list:
    #    print t
    
    return new_possible_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list


def get_aatypes_from_all_H_C_resonpairs(raw_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list,
                                        Num_of_TOCSY_resonances,
                                        args):
    """
        FUNCTION to return the matched amino acid types for a particular TOCSY AAIG by using all its H,C resonance pairs.
        
        ARGUMENTS:
        raw_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list:   list of lists of the form (aa type, Carbon name, Hydrogen name, weighted average probability, TOCSY_reson_index, H resonance, C resonance, carbon group)
                                                                        This list includes also negative probabilities wherever only the Carbon resonance match was considered.
        
        RETURN:
        matching_aaTypesProbTuple_list:     list of tuples of the form (aa type, average probability)
    """
    global aatype_maxH_C_pairs_dict, Prob_CS
    
    #print "DEBUG: raw_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list=", raw_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list
    # If the last arguments in True, it will remove all C-based predictions
    possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list = consider_only_Carbon_matches(raw_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list, True, False)
    #print "DEBUG get_aatypes_from_all_H_C_resonpairs: possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list=",possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list
    # The following sorting ensures that the alternative C-H pairs are sorted according to their weighted average probability 
    possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list.sort(key=itemgetter(0,1,2,3), reverse=True)  # sort by every column except the last 4
    
    ## The following block keeps only the matches for a particular aatype
    TRI_set = set([x[4] for x in possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list])
    TRI_num = len(TRI_set)  # the number of different peaks available for this spin system
    previous_aatype = None
    aatypeResonpairMatchesTuple_list = []
    selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list = []
    aatype_predictionScore_dict = {}     # dict aa type -> total prediction score of this aa type.
                                         # Used only if args.PROB_PRODUCT == True (TRUE by default)!
    for group8 in possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list:
        #print "DEBUG: group8=", group8
        aa_type = group8[0]
        weighted_average_probability = group8[3]    # = [wH * hist(H) + wC * hist(C)]/2
        if aa_type != previous_aatype and previous_aatype != None:
            #print "DEBUG: aa_type=",aa_type,"previous_aatype=",previous_aatype,"aatypeResonpairMatchesTuple_list=",aatypeResonpairMatchesTuple_list
            if args.PROB_PRODUCT == True:
                new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list, \
                aatype2save, \
                aatype_predictionScore = select_correct_H_C_resonpairs2(aatypeResonpairMatchesTuple_list,
                                                                          args,
                                                                          get_aa_type=True)
                aatype_predictionScore_dict[aatype2save] = aatype_predictionScore
            else:
                new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list = select_correct_H_C_resonpairs(aatypeResonpairMatchesTuple_list)
            # CHECK IF ALL THE PEAKS HAVE BEEN ASSIGNED TO C-H TYPES FOR THE CURRENT AA TYPE. IF NOT DO NOT SAVE
            # THIS AA TYPE AS POSSIBLE PREDICTION FOR THE CURRENT TOCSY index group
            if len(new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list) >= (args.ASSIGNMENT_CUTOFF * TRI_num):
                selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list.extend(new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list)
            else:
                pass
                #print "DEBUG: the assigned C-H pairs for aa type", previous_aatype, "were below the requested threshold (",args.ASSIGNMENT_CUTOFF * TRI_num,")"
                #print "DEBUG: new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list=", new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list
            #print "DEBUG: selected correct H-C resonpairs!"
            aatypeResonpairMatchesTuple_list = []
        
        aatypeResonpairMatchesTuple_list.append(group8)
        previous_aatype = aa_type
    #print "DEBUG get_aatypes_from_all_H_C_resonpairs: aa_type=",aa_type,"previous_aatype=",previous_aatype,"aatypeResonpairMatchesTuple_list=",aatypeResonpairMatchesTuple_list
    if len(selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list) > 0:
    # SELECT THE OPTIMUM C-H COMBINATIONS FOR THE LAST aa_type
        if args.PROB_PRODUCT == True:
            new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list, aatype2save, aatype_predictionScore = select_correct_H_C_resonpairs2(aatypeResonpairMatchesTuple_list,
                                                                                                                                      args,
                                                                                                                                      get_aa_type=True)  # do the last aa_type
            aatype_predictionScore_dict[aatype2save] = aatype_predictionScore
        else:
            new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list = select_correct_H_C_resonpairs(aatypeResonpairMatchesTuple_list)  # do the last aa_type
        # CHECK IF ALL THE PEAKS HAVE BEEN ASSIGNED TO C-H TYPES FOR THE CURRENT AA TYPE. IF NOT DO NOT SAVE THIS AA TYPE AS POSSIBLE PREDICTION FOR THE CURRENT TOCSY index group
        if len(new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list) >= (args.ASSIGNMENT_CUTOFF * TRI_num):
            selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list.extend(new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list)
        else:
            pass
            #print "DEBUG: the assigned C-H pairs for aa type", previous_aatype, "were below the requested threshold (",args.ASSIGNMENT_CUTOFF * TRI_num,")"
            #print "DEBUG: new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list=", new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list
        #print "DEBUG: selected correct H-C resonpairs!"
    # IF NO C-H TYPE ASSIGNMENTS WERE MADE, REPEAT THE PROCEDURE INCLUDING ONLY CARBON RESONANCES
    elif len(selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list) == 0:
        print(bcolors.BOLDBLUE + "NO C-H TYPE ASSIGNMENTS WERE MADE, REPEAT THE PROCEDURE INCLUDING ONLY CARBON RESONANCES" + bcolors.ENDBOLD)

        possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list = consider_only_Carbon_matches(raw_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list, False, False)
        # The following sorting ensures that the alternative C-H pairs are sorted according to their weighted average probability 
        possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list.sort(key=itemgetter(0,1,2,3), reverse=True)  # sort by every column except the last 3
        ## The following block keeps only the matches for a particular aatype
        previous_aatype = None
        aatypeResonpairMatchesTuple_list = []
        selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list = []
        #print "DEBUG: 2032 possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list=", possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list
        for group8 in possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list:
            #print "DEBUG: group8=", group8
            aa_type = group8[0]
            weighted_average_probability = group8[3]    # = [wH * hist(H) + wC * hist(C)]/2
            if aa_type != previous_aatype and previous_aatype != None:
                #print "DEBUG: aa_type=",aa_type,"previous_aatype=",previous_aatype,"aatypeResonpairMatchesTuple_list=",aatypeResonpairMatchesTuple_list
                if args.PROB_PRODUCT == True:
                    new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list, aatype2save, aatype_predictionScore = select_correct_H_C_resonpairs2(aatypeResonpairMatchesTuple_list,
                                                                                                                                              args,
                                                                                                                                              get_aa_type=True)
                    aatype_predictionScore_dict[aatype2save] = aatype_predictionScore
                else:
                    new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list = select_correct_H_C_resonpairs(aatypeResonpairMatchesTuple_list)
                # CHECK IF ALL THE PEAKS HAVE BEEN ASSIGNED TO C-H TYPES FOR THE CURRENT AA TYPE. IF NOT DO NOT SAVE THIS AA TYPE AS POSSIBLE PREDICTION FOR THE CURRENT TOCSY index group
                if len(new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list) >= (args.ASSIGNMENT_CUTOFF * TRI_num):
                    selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list.extend(new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list)
                else:
                    pass
                    #print "DEBUG: the assigned C-H pairs for aa type", previous_aatype, "were below the requested threshold (",args.ASSIGNMENT_CUTOFF * TRI_num,")"
                    #print "DEBUG: new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list=", new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list
                #print "DEBUG: selected correct H-C resonpairs!"
                aatypeResonpairMatchesTuple_list = []
            
            aatypeResonpairMatchesTuple_list.append(group8)
            previous_aatype = aa_type
        #print "DEBUG: aa_type=",aa_type,"previous_aatype=",previous_aatype,"aatypeResonpairMatchesTuple_list=",aatypeResonpairMatchesTuple_list
        if args.PROB_PRODUCT == True:
            new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list, aatype2save, aatype_predictionScore = select_correct_H_C_resonpairs2(aatypeResonpairMatchesTuple_list,
                                                                                                                                      args,
                                                                                                                                      get_aa_type=True)  # do the last aa_type
            aatype_predictionScore_dict[aatype2save] = aatype_predictionScore
        else:
            new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list = select_correct_H_C_resonpairs(aatypeResonpairMatchesTuple_list)  # do the last aa_type
        # CHECK IF ALL THE PEAKS HAVE BEEN ASSIGNED TO C-H TYPES FOR THE CURRENT AA TYPE. IF NOT DO NOT SAVE THIS AA TYPE AS POSSIBLE PREDICTION FOR THE CURRENT TOCSY index group
        if len(new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list) >= (args.ASSIGNMENT_CUTOFF * TRI_num):
            selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list.extend(new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list)
        else:
            pass
            #print "DEBUG: the assigned C-H pairs for aa type", previous_aatype, "were below the requested threshold (",args.ASSIGNMENT_CUTOFF * TRI_num,")"
            #print "DEBUG: new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list=", new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list
        #print "DEBUG: selected correct H-C resonpairs!"
    
    #print "DEBUG: selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list=", selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list
    aatype_occupancy_dict = {}  # dict with key the aa types and values how many times they were found during the H,C resonance analysis
    aatype_probsum_dict = {}    # dictionary with the aa type --> sum of probabilities of individual matched H,C resonance pairs
    aatype_CHnucleiType_presenceProbSum_mdict = tree()    # multidimensional dict: aa type --> C,H nuclei pair type --> probability of presence or absence (eq 6 in RESCUE2 paper)
    aatype_CHnucleiType_TresonIndex_mdict = tree()    # multidimensional dict: aa type --> C,H nuclei pair type --> TOCSY_reson_index
    aa_type_set = set()
    aatype_validC_H_pairsList_dict = {}
    #print "DEBUG: selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list=", selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list
    for group8 in selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list:
        #print "DEBUG: group8=",group8
        aa_type = group8[0]
        Cname = group8[1]
        Hname = group8[2]
        weighted_average_probability = group8[3]
        TOCSY_reson_index = group8[4]
        #TOCSY_reson_index = group5[4]  # we don't need this
        #if (aa_type, Cname, Hname) == previous_aatype_C_H_tuple:  # elegant way to keep only the (aatype, Cname, Hname) with the highest weighted_average_probability, e.g. we may have:
        #    # ('CYS', 'CA', 'HA', 0.012993371769822381), ('CYS', 'CA', 'HA', 0.002004997884084288), ('CYS', 'CA', 'HA', 0.0018243946693470656) but we keep only the 1st 
        #    continue
        try:
            if (Cname+"_"+Hname in aatype_validC_H_pairsList_dict[aa_type]):  # if we have saved twice this aatype,C,H combination print ERROR & exit
                print("ERROR: ",Cname+"_"+Hname,"already found !!!")
                sys.exit(1)
            aatype_validC_H_pairsList_dict[aa_type].append(Cname+"_"+Hname)
            aa_type_set.add(aa_type)
            aatype_occupancy_dict[aa_type] += 1
            ## THERE 2 DIFFENT WAYS TO CALCULATE THE OVERALL WEIGHTED AVERAGE PROBABILITY
            #aatype_probsum_dict[aa_type] += weighted_average_probability    # 1st WAY
            aatype_probsum_dict[aa_type] *= weighted_average_probability    # 2nd WAY
            aatype_CHnucleiType_presenceProbSum_mdict[aa_type][Cname+"_"+Hname] = (args.H_weight * Prob_CS[aa_type][Hname] + args.C_weight * Prob_CS[aa_type][Cname])/float((args.H_weight + args.C_weight))  # 2nd WAY
            aatype_CHnucleiType_TresonIndex_mdict[aa_type][Cname+"_"+Hname] = TOCSY_reson_index
        except KeyError:
            aatype_validC_H_pairsList_dict[aa_type] = [Cname+"_"+Hname]
            #if args.ALLOW_SINGLE_CH_PAIRS:
            aa_type_set.add(aa_type)
            aatype_occupancy_dict[aa_type] = 1
            aatype_probsum_dict[aa_type] = weighted_average_probability     # COMMON IN 1st & 2nd WAY
            aatype_CHnucleiType_presenceProbSum_mdict[aa_type][Cname+"_"+Hname] = (args.H_weight * Prob_CS[aa_type][Hname] + args.C_weight * Prob_CS[aa_type][Cname])/float((args.H_weight + args.C_weight))  # 2nd WAY
            aatype_CHnucleiType_TresonIndex_mdict[aa_type][Cname+"_"+Hname] = TOCSY_reson_index
            
        previous_aatype_C_H_tuple = (aa_type, Cname, Hname)
    
    matching_aaTypesProbTuple_list = []
    #print "DEBUG: Num_of_TOCSY_resonances=", Num_of_TOCSY_resonances
    #print "DEBUG: args.ASSIGNMENT_CUTOFF=", args.ASSIGNMENT_CUTOFF
    #print "DEBUG: aatype_occupancy_dict=",aatype_occupancy_dict
    #print "DEBUG: aatype_validC_H_pairsList_dict=",aatype_validC_H_pairsList_dict
    #print "DEBUG: aa_type_set=", aa_type_set
    for AA_type in aa_type_set:
        if aatype_occupancy_dict[AA_type] > aatype_maxH_C_pairs_dict[AA_type]:    # if we found more C-H pairs than the maximum number for that aa type, ommit it
            #print "DEBUG: ommiting AA_type=",AA_type,"due to higher occupancy (",aatype_occupancy_dict[AA_type],") than the maximum number of allowed C-H pairs (",aatype_maxH_C_pairs_dict[AA_type],")!"
            continue
        # consider this aa type as a match ONLY IF A PORTION OF THE TOCSY RESONANCES MATCH AND if ...
        # The latter excludes cases were we had 5 TOCSY resonances, we found 4 mathcing Val, but Val has only 4 C-H pairs in maximum
        # But when we have 1 TOCSY resonance, we found matching Cys, but Cys has 3 aatype_maxH_C_pairs_dict
        #print "DEBUG: AA_type=",AA_type,"aatype_occupancy_dict[AA_type]=", aatype_occupancy_dict[AA_type], ">= args.ASSIGNMENT_CUTOFF * Num_of_TOCSY_resonances=",args.ASSIGNMENT_CUTOFF * Num_of_TOCSY_resonances
        #print "DEBUG: aatype_maxH_C_pairs_dict[AA_type]=", aatype_maxH_C_pairs_dict[AA_type], ">= Num_of_TOCSY_resonances=", Num_of_TOCSY_resonances
        if (aatype_occupancy_dict[AA_type] >= (args.ASSIGNMENT_CUTOFF * Num_of_TOCSY_resonances)) and (aatype_maxH_C_pairs_dict[AA_type] >= Num_of_TOCSY_resonances):   
            try:
                #print "DEBUG: Saving AA_type",AA_type," as valid prediction."
                ## THERE ARE A FEW ALTERNATIVE SCORING SCHEMES HERE:
                if args.PROBABILITY_MODEL == 1:
                    if args.PROB_PRODUCT == True:
                        score = aatype_predictionScore_dict[AA_type]
                    else:
                        score = get_combined_weighted_and_presence_probability(AA_type, aatype_probsum_dict, aatype_CHnucleiType_presenceProbSum_mdict)
                elif args.PROBABILITY_MODEL == 2:
                    if args.PROB_PRODUCT == True:
                        score = aatype_predictionScore_dict[AA_type]
                    else:
                        score = aatype_probsum_dict[AA_type]
                #score = aatype_probsum_dict[AA_type]/float(Num_of_TOCSY_resonances) # divide by the total num of TOCSY resonance pairs; 1st WAY
                # CYS GET WORSE & Z-SCORES MORE NEGATIVE # score = aatype_probsum_dict[AA_type]/float(aatype_occupancy_dict[AA_type]) # divide by the prob sum by the occupancy of that aatype
                # BAD! # score = float(aatype_maxH_C_pairs_dict[AA_type]) * aatype_probsum_dict[AA_type]/(aatype_occupancy_dict[AA_type]**2)
                matching_aaTypesProbTuple_list.append( (AA_type, score) )
            except TypeError:
                print(aatype_probsum_dict[aa_type], Num_of_TOCSY_resonances)
                sys.exit(1)
    
    
    # FINALLY IF NO AA TYPE PREDICTIONS WERE MADE, REPEAT THE PROCEDURE INCLUDING ONLY CARBON RESONANCES
    if len(matching_aaTypesProbTuple_list) == 0:
        print("NO AA TYPE PREDICTIONS WERE MADE, THUS ACTIVATING CARBON-BASED PREDICTIONS!")
        possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list = consider_only_Carbon_matches(raw_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list, False, False)
        # The following sorting ensures that the alternative C-H pairs are sorted according to their weighted average probability 
        possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list.sort(key=itemgetter(0,1,2,3), reverse=True)  # sort by every column except the last 3
        ## The following block keeps only the matches for a particular aatype
        previous_aatype = None
        aatypeResonpairMatchesTuple_list = []
        selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list = []
        for group8 in possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list:
            #print "DEBUG: group8=", group8
            aa_type = group8[0]
            weighted_average_probability = group8[3]    # = [wH * hist(H) + wC * hist(C)]/2
            if aa_type != previous_aatype and previous_aatype != None:
                #print "DEBUG: aa_type=",aa_type,"previous_aatype=",previous_aatype,"aatypeResonpairMatchesTuple_list=",aatypeResonpairMatchesTuple_list
                if args.PROB_PRODUCT == True:
                    # selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list.extend( select_correct_H_C_resonpairs2(aatypeResonpairMatchesTuple_list) )
                    new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list, aatype2save, aatype_predictionScore = select_correct_H_C_resonpairs2(aatypeResonpairMatchesTuple_list,
                                                                                                                                              args,
                                                                                                                                              get_aa_type=True)  # do the last aa_type
                    selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list.extend(new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list)
                    aatype_predictionScore_dict[aatype2save] = aatype_predictionScore
                else:
                    selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list.extend( select_correct_H_C_resonpairs(aatypeResonpairMatchesTuple_list) )
                aatypeResonpairMatchesTuple_list = []
            
            aatypeResonpairMatchesTuple_list.append(group8)
            previous_aatype = aa_type
        #print "DEBUG: aa_type=",aa_type,"previous_aatype=",previous_aatype,"aatypeResonpairMatchesTuple_list=",aatypeResonpairMatchesTuple_list
        if args.PROB_PRODUCT == True:
            # selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list.extend( select_correct_H_C_resonpairs2(aatypeResonpairMatchesTuple_list) )  # do the last aa_type
            new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list, aatype2save, aatype_predictionScore = select_correct_H_C_resonpairs2(aatypeResonpairMatchesTuple_list,
                                                                                                                                      args,
                                                                                                                                      get_aa_type=True)  # do the last aa_type
            selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list.extend(new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list)
            aatype_predictionScore_dict[aatype2save] = aatype_predictionScore
        else:
            selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list.extend( select_correct_H_C_resonpairs(aatypeResonpairMatchesTuple_list) )  # do the last aa_type
    
        #print "DEBUG: selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list=", selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list
        aatype_occupancy_dict = {}  # dict with key the aa types and values how many times they were found during the H,C resonance analysis
        aatype_probsum_dict = {}    # dictionary with the aa type --> sum of probabilities of individual matched H,C resonance pairs
        aatype_CHnucleiType_presenceProbSum_mdict = tree()    # multidimensional dict: aa type --> C,H nuclei pair type --> probability of presence or absence (eq 6 in RESCUE2 paper)
        aatype_CHnucleiType_TresonIndex_mdict = tree()    # multidimensional dict: aa type --> C,H nuclei pair type --> TOCSY_reson_index
        aa_type_set = set()
        aatype_validC_H_pairsList_dict = {}
        print("DEBUG: selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list=", selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list)
        for group8 in selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list:
            #print "DEBUG: group8=",group8
            aa_type = group8[0]
            Cname = group8[1]
            Hname = group8[2]
            weighted_average_probability = group8[3]
            TOCSY_reson_index = group8[4]
            #if (aa_type, Cname, Hname) == previous_aatype_C_H_tuple:  # elegant way to keep only the (aatype, Cname, Hname) with the highest weighted_average_probability, e.g. we may have:
            #    # ('CYS', 'CA', 'HA', 0.012993371769822381), ('CYS', 'CA', 'HA', 0.002004997884084288), ('CYS', 'CA', 'HA', 0.0018243946693470656) but we keep only the 1st 
            #    continue
            try:
                if (Cname+"_"+Hname in aatype_validC_H_pairsList_dict[aa_type]):  # if we have saved twice this aatype,C,H combination print ERROR & exit
                    print("ERROR: ",Cname+"_"+Hname,"already found !!!")
                    sys.exit(1)
                aatype_validC_H_pairsList_dict[aa_type].append(Cname+"_"+Hname)
                aa_type_set.add(aa_type)
                aatype_occupancy_dict[aa_type] += 1
                ## THERE 2 DIFFENT WAYS TO CALCULATE THE OVERALL WEIGHTED AVERAGE PROBABILITY
                #aatype_probsum_dict[aa_type] += weighted_average_probability    # 1st WAY
                aatype_probsum_dict[aa_type] *= weighted_average_probability    # 2nd WAY
                aatype_CHnucleiType_presenceProbSum_mdict[aa_type][Cname+"_"+Hname] = (args.H_weight * Prob_CS[aa_type][Hname] + args.C_weight * Prob_CS[aa_type][Cname])/float((args.H_weight + args.C_weight))  # 2nd WAY
                aatype_CHnucleiType_TresonIndex_mdict[aa_type][Cname+"_"+Hname] = TOCSY_reson_index
            except KeyError:
                aatype_validC_H_pairsList_dict[aa_type] = [Cname+"_"+Hname]
                #if args.ALLOW_SINGLE_CH_PAIRS:
                aa_type_set.add(aa_type)
                aatype_occupancy_dict[aa_type] = 1
                aatype_probsum_dict[aa_type] = weighted_average_probability     # COMMON IN 1st & 2nd WAY
                aatype_CHnucleiType_presenceProbSum_mdict[aa_type][Cname+"_"+Hname] = (args.H_weight * Prob_CS[aa_type][Hname] + args.C_weight * Prob_CS[aa_type][Cname])/float((args.H_weight + args.C_weight))  # 2nd WAY
                aatype_CHnucleiType_TresonIndex_mdict[aa_type][Cname+"_"+Hname] = TOCSY_reson_index
                
            previous_aatype_C_H_tuple = (aa_type, Cname, Hname)
        
        matching_aaTypesProbTuple_list = []
        #print "DEBUG: Num_of_TOCSY_resonances=", Num_of_TOCSY_resonances
        #print "DEBUG: args.ASSIGNMENT_CUTOFF=", args.ASSIGNMENT_CUTOFF
        #print "DEBUG: aatype_occupancy_dict=",aatype_occupancy_dict
        #print "DEBUG: aatype_validC_H_pairsList_dict=",aatype_validC_H_pairsList_dict
        #print "DEBUG: aa_type_set=", aa_type_set
        for AA_type in aa_type_set:
            if aatype_occupancy_dict[AA_type] > aatype_maxH_C_pairs_dict[AA_type]:    # if we found more C-H pairs than the maximum number for that aa type, ommit it
                #print "DEBUG: ommiting AA_type=",AA_type,"due to higher occupancy (",aatype_occupancy_dict[AA_type],") than the maximum number of allowed C-H pairs (",aatype_maxH_C_pairs_dict[AA_type],")!"
                continue
            # consider this aa type as a match ONLY IF A PORTION OF THE TOCSY RESONANCES MATCH AND if ...
            # The latter excludes cases were we had 5 TOCSY resonances, we found 4 mathcing Val, but Val has only 4 C-H pairs in maximum
            # But when we have 1 TOCSY resonance, we found matching Cys, but Cys has 3 aatype_maxH_C_pairs_dict
            #print "DEBUG: AA_type=",AA_type,"aatype_occupancy_dict[AA_type]=", aatype_occupancy_dict[AA_type], ">= args.ASSIGNMENT_CUTOFF * Num_of_TOCSY_resonances=",args.ASSIGNMENT_CUTOFF * Num_of_TOCSY_resonances
            #print "DEBUG: aatype_maxH_C_pairs_dict[AA_type]=", aatype_maxH_C_pairs_dict[AA_type], ">= Num_of_TOCSY_resonances=", Num_of_TOCSY_resonances
            if (aatype_occupancy_dict[AA_type] >= (args.ASSIGNMENT_CUTOFF * Num_of_TOCSY_resonances)) and (aatype_maxH_C_pairs_dict[AA_type] >= Num_of_TOCSY_resonances):   
                try:
                    #print "DEBUG: Saving AA_type",AA_type," as valid prediction."
                    ## THERE ARE A FEW ALTERNATIVE SCORING SCHEMES HERE:
                    if args.PROBABILITY_MODEL == 1:
                        if args.PROB_PRODUCT == True:
                            score = aatype_predictionScore_dict[AA_type]
                        else:
                            score = get_combined_weighted_and_presence_probability(AA_type, aatype_probsum_dict, aatype_CHnucleiType_presenceProbSum_mdict)
                    elif args.PROBABILITY_MODEL == 2:
                        if args.PROB_PRODUCT == True:
                            score = aatype_predictionScore_dict[AA_type]
                        else:
                            score = aatype_probsum_dict[AA_type]
                    #score = aatype_probsum_dict[AA_type]/float(Num_of_TOCSY_resonances) # divide by the total num of TOCSY resonance pairs; 1st WAY
                    # CYS GET WORSE & Z-SCORES MORE NEGATIVE # score = aatype_probsum_dict[AA_type]/float(aatype_occupancy_dict[AA_type]) # divide by the prob sum by the occupancy of that aatype
                    # BAD! # score = float(aatype_maxH_C_pairs_dict[AA_type]) * aatype_probsum_dict[AA_type]/(aatype_occupancy_dict[AA_type]**2)
                    matching_aaTypesProbTuple_list.append( (AA_type, score) )
                except TypeError:
                    print(aatype_probsum_dict[aa_type], Num_of_TOCSY_resonances)
                    sys.exit(1)
    
    
    return matching_aaTypesProbTuple_list
########################################################## END OF FUNCTION DEFINITIONS #########################################################