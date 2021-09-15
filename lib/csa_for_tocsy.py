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

import sys, os

# Import 4D-CHAINS libraries
sys.path.append( os.path.dirname( os.path.dirname( os.path.abspath(__file__) ) ) )  # import the top level package directory
from lib.csa_for_noesy import *

from ete3 import Tree
import gc
from operator import itemgetter


def get_protein_sequence(FASTA_fname):

    protein_sequence_string = ""
    with open(FASTA_fname, 'r') as f:
        for line in f:
            if not re.match("^>", line):
                protein_sequence_string += re.sub(r"[^A-Z]", "", line)
    protein_sequence_list = list(protein_sequence_string)
    return protein_sequence_list


################################################# AUTOMATIC ASSIGNMENT OF AMINO ACID TYPES TO TOCSY RESONANCES FUNCTION DEFINITIONS ####################################################

def get_aatypes_from_H_C_resonpair(aa,
                                   Hreson,
                                   Creson,
                                   TOCSY_reson_index,
                                   histload,
                                   args):
    """
        FUNCTION to find the 1D-hist probabilities of all possible aliphatic H,C types, given an aa type and a pair of H,C resonances.
        
        ARGUMENTS:
        aa:                 the assigned aa type for the current C-H resonanse pair
        TOCSY_reson_index:  an index indicating the position of this H-C resonance pair in the TOCSY group it belongs to
        
        RETURNS:
        valid_matches_list:   list of lists of the form (aa type, Carbon name, Hydrogen name, weighted average probability, TOCSY_reson_index, H resonance, C resonance, carbon group)
    """
    
    global aa_carbonBondedHydrogensDict_dict
    
    def is_valid_hydrogen(aa_type, C_name, H_name, C_bondedH_dict):
        """
            FUNCTION to check whether a given hydrogen is covalently bonded to the given carbon of a given aa type
            ARGUMENTS:
            C_bondedH_dict:     dictonary with keys the carbon names of a given amino acid and values the respective covalently hydrogen names
        """
        
        if H_name in C_bondedH_dict[C_name]:
            return True
        
        return False
    
    print("DEBUG: Hreson=",Hreson,"Creson=",Creson)
    hydrogen_matches_list = []   # list of tuples of the form (aa type,hydrogen,probability)
    carbon_matches_list = []    # list of tuples of the form (aa type,carbon,probability)
    histcalc = ProbHist_Calculator()

    for hydrogen in list(histload.aa_hydrogen_binDensityList_mdict[aa].keys()):    # iterate over all hydrogens of the current amino acid
        #print "DEBUG: aa=", aa, "hydrogen=", hydrogen, "Hreson=", Hreson, "bin_array=", bin_array, "density_array=", density_array
        probability = histcalc.get_H_probability_from_1Dhistogram(aa=aa, Hname=hydrogen, Hreson=Hreson, histload=histload)
        #print "DEBUG: probability=", probability
        if probability > 0.0:
            hydrogen_matches_list.append((aa,hydrogen,probability))
    #print "DEBUG: hydrogen_matches_list=", hydrogen_matches_list
    
    for carbon in list(histload.aa_carbon_binDensityList_mdict[aa].keys()):    # iterate over all carbons of the current amino acid
        probability = histcalc.get_C_probability_from_1Dhistogram(aa=aa, Cname=carbon, Creson=Creson, histload=histload)
        #print "DEBUG: aa=",aa,"carbon=",carbon, "Creson=", Creson, "probability=",probability
        if probability > 0.0:
            carbon_matches_list.append((aa,carbon,probability))
    
    # Keep only the amino acids for which were found both the carbon and the respective covalently bonded hydrogen
    # The probability of each prediction will be (Carbon_probability + Hydrogen_probability) /2
    matches_list = [] # list of tuples of the form (aa type, Carbon name, Hydrogen name, overall probability, )
    for hydrogen_match in hydrogen_matches_list:    # iterate over all hydrogen matches found in the previous step
        for carbon_match in carbon_matches_list:    # iterate over all carbons matches found in the previous step
            if hydrogen_match[0] == carbon_match[0]:    # if they belong to the same amino acid
                aa_type = hydrogen_match[0]
                H_name = hydrogen_match[1]
                C_name = carbon_match[1]
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
                matches_list.append((aa_type, C_name, H_name, weighted_average_probability))
    
    # TODO: calculate Z-scores from the histogram probabilities of the remaining unmatched carbons and ignore the hydrogens for those C with Z-score > 0.
    valid_matches_list = []    # the same list but contains only the most probable C-H matching pair from each aminoacid. 
    previous_aatype = None
    for quartet in matches_list:
        print("DEBUG: quartet=",quartet)
        aa_type = quartet[0]
        C_name = quartet[1]
        H_name = quartet[2]
        if is_valid_hydrogen(aa_type, C_name, H_name, aa_carbonBondedHydrogensDict_dict[aa_type]) == False: # if this carbon is not covalently bonded to this hydrogen, skip it
            continue
        #OBSOLETE # if aa_type == previous_aatype:  # elegant way to keep only the (aatype, C_name, H_name) with the highest probability, e.g. we may have:
        #OBSOLETE # ('VAL', 'CB', 'HB', 1.4736221632773356e-05), ('VAL', 'CG2', 'HG2', 0.01135282415016002), ('VAL', 'CG1', 'HG1', 0.01433564978569348) but we keep only the 3rd
        #OBSOLETE # othewise we would keep 'CG2_HG2', 'CG1_HG1', 'CB_HB' for just one C-H resonance pair!
        #OBSOLETE #    continue    
        valid_matches_list.append([aa_type, C_name, H_name, quartet[3], TOCSY_reson_index, Hreson, Creson, None])  # add also the index of this H-C pair in the TOCSY group it
        # belongs to, but also the H and C resonances for future usage
        previous_aatype = aa_type
    
    print("DEBUG: matches_list=", matches_list)
    print("DEBUG: valid_matches_list=", valid_matches_list)
    return valid_matches_list


def get_all_assignments_from_H_C_resonpair_2Dhist(aa,
                                                  Hreson,
                                                  Creson,
                                                  TOCSY_reson_index,
                                                  histload,
                                                  useonlyCarbon=False):
    """
        FUNCTION to find the 2D-hist probabilities of all possible aliphatic H,C types, given an aa type and a pair of H,C resonances.
        
        ARGUMENTS:
        TOCSY_reson_index:  an index indicating the position of this H-C resonance pair in the TOCSY group it belongs to
        useonlyCarbon:      use only Carbon 1D histograms when 2Dhist probability is 0. However, this generated problems like 3 Carbons in the same C-group
                            in cs_assignment.py and cs_assignment_CAB.py therefore this option was deactivated.
        RETURNS:
        matches_list: list of lists of the form [aa type, Carbon name, Hydrogen name, 2D-hist probability, TOCSY_reson_index, H resonance, C resonance, None]
    """

    #print "DEBUG: Hreson=",Hreson,"Creson=",Creson
    matches_list = [] # list of tuples of the form (aa type, Carbon name, Hydrogen name, overall probability, )
    carbon_matches_list = []    # list of tuples of the form (aa type,carbon,probability); here aa should be the same in all the tuples of the list!
    histcalc = ProbHist_Calculator()
    for carbon in list(histload.aa_carbon_binDensityList_mdict[aa].keys()):    # iterate over all carbons of the current amino acid
        probability = histcalc.get_C_probability_from_1Dhistogram(aa, carbon, Creson, histload)
        # print "DEBUG: aa=",aa,"carbon=",carbon,"probability=",probability
        if probability > 0.0:
            carbon_matches_list.append((aa,carbon,probability))
        
    for CH_pair in list(histload.aa_CHpair_binProbabilityList_mdict[aa].keys()):
        # print "DEBUG: checking 2D hist probability for aa=",aa, "CH_pair=", CH_pair, "Hreson=", Hreson, "Creson", Creson
        probability = histcalc.get_CH_probability_from_2Dhistogram(CH_pair, aa, Hreson, Creson, histload)
        # print "DEBUG: 2D hist probability for %s is %f" %(CH_pair, probability)
        # Keep only the amino acids for which were found both the carbon and the respective covalently bonded hydrogen
        # The probability of each prediction is given by the 2D histogram
        if probability > 0.0:
            C_name, H_name = CH_pair.split("-")
            matches_list.append([aa, C_name, H_name, probability, TOCSY_reson_index, Hreson, Creson, None])  # add also the index of this C-H pair in the TOCSY group it
        elif probability == 0.0 and useonlyCarbon == True:    # if the probability of the 2D hist at this point is 0, use only the Carbon 1D histogram to get the probability
            print("WARNING: 2D hist probability is 0. Using Carbon 1D histogram to get the probability!")
            C_name, H_name = CH_pair.split("-")
            #print "DEBUG: C_name=", C_name, "carbon_matches_list=", carbon_matches_list
            carbon_match = [match for match in carbon_matches_list if match[1]==C_name]
            if len(carbon_match) == 0:  # if the hist(C) of this C_name was zero, don't save it 
                continue
            #print "DEBUG: len(carbon_match) is ", len(carbon_match)," and it should be 1 !"
            #print "DEBUG: carbon_match=", carbon_match
            weighted_average_probability = -1 * carbon_match[0][2] #  but first make it negative to distiguish it from weighted average probabilities
            print("DEBUG: only the Carbon 1D histogram used:", [aa, C_name, H_name, weighted_average_probability, TOCSY_reson_index, Hreson, Creson, None])
            matches_list.append([aa, C_name, H_name, weighted_average_probability, TOCSY_reson_index, Hreson, Creson, None]) 
            
    print("DEBUG: returning matches_list=", matches_list)
    return matches_list


def select_correct_H_C_resonpairs(aatypeResonpairMatchesTuple_list):
    """
        FUNCTION  to match the TOCSY resonance pairs to as much as possible C-H resonance pairs of a PARTICULAR aatype. I may need to used a Genetic Algorith to do that
        if the current implementation proves to be error prone.
        
        ARGUMETS:
        aatypeResonpairMatchesTuple_list:   list of lists of the form [aa_type, C, H, weighted average probability, selected TOCSY reson index, H resonance, C resonance, ClusterID] where aa_type is the same in the whole list (e.g. "PRO")
        
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
        C_name = match[1]
        Creson = match[6]
        for correct_match in correct_C_H_resonpair_matches_list:
            # if the carbon name is already in the correct matches but its resonance is different, return false
            print("DEBUG: C_name=",C_name, "correct_match[1]=", correct_match[1], "Creson=", Creson, "correct_match[6]=", correct_match[6])
            if C_name == correct_match[1] and approx_equal(Creson, correct_match[6], 0.3) == False:
                print("DEBUG: conflict found, do_carbons_match() returning False!")
                return False
        
        return True # otherwise return true
        
        
    print("DEBUG: entered function select_correct_H_C_resonpairs() with aatypeResonpairMatchesTuple_list=",aatypeResonpairMatchesTuple_list)
    aa_type = aatypeResonpairMatchesTuple_list[0][0]    # recall that we have only one aa type in this script
    correct_C_H_resonpair_matches_list = [] # list of lists, see the Function definition
    used_TOCSY_reson_indices_list = []  # list of assigned peaks (TOCSY reson indices)
    used_C_H_pairs_list = []    # list of assigned C and H nucleus types
    
    # FIRST TAKE CARE OF CARBON GROUPING AND THEIR PROTONS (e.g. CB, CB2, CB3)
    clustIDs_list = [x[7] for x in aatypeResonpairMatchesTuple_list]    # list of Cluster ID of each possible assignment
    clustIDs_set = set(clustIDs_list)
    # FIRST OF ALL SORT THE CLUSTERS ACCORDING TO THE TOTAL PROBABILITY OF THEIR CARBONS
    clustID_highestTotProb_dict = {}
    for clustID in clustIDs_set:
        TresonIndex_set = set([x[4] for x in aatypeResonpairMatchesTuple_list if x[7]==clustID])    # set with all different peaks (TOCSY reson indices) of the current Cluster ID
        # NOW USE THE Treson indices WITH MULTIPLE TYPES OF CARBONS
        print("DEBUG: clustIDs_list.count(clustID) =", clustIDs_list.count(clustID), " len(TresonIndex_set)=", len(TresonIndex_set))
        #if clustIDs_list.count(clustID) > 1 and len(TresonIndex_set) > 1:    # if there are 2 or more possible assignments and at least 2 different peaks within the same Cluster ID
        nucleusAssignmentsforCurrentClustID_list = [x for x in aatypeResonpairMatchesTuple_list if x[7]==clustID]    # save all the members of this Cluster ID
        # Recall that nucleusAssignmentsforCurrentClustID_list is already sorted by lines 1,2,3, BUT WE MUST ...
        nucleusAssignmentsforCurrentClustID_list.sort(key=itemgetter(3), reverse=True)  # ... SORT BY THE PREDICTION PROBABILITY
        
        # Make a list of the carbon types in each Treson index of this Cluster and its maximum probability
        C_totalProb_dict = {}   # the total probability of each C type in this Cluster. Use this to decide the assignment
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
            print("DEBUG: TRI = ", TRI)
            for C in carbon_set:
                if C not in list(C_prob_dict.keys()): # if this C type was not in possible predictions if this TRI, set its total probability to 0 (it will be completely excluded!)
                    C_totalProb_dict[C] = 0
                    continue
                try:
                    C_totalProb_dict[C] *= C_prob_dict[C]
                except KeyError:
                    C_totalProb_dict[C] = C_prob_dict[C]
        
        print("DEBUG: C_totalProb_dict=", C_totalProb_dict)
        sorted_C_totalProb_list = sorted(list(C_totalProb_dict.items()), key=itemgetter(1), reverse=True)   # sorted C_totalProb_list, e.g. [(CB, 0.9), (CG, 0.7)]
        
        clustID_highestTotProb_dict[clustID] = sorted_C_totalProb_list[0][1]    # save the total probability of the first element only, since the list was sorted
    # AND THIS IS THE LIST OF CLUSTER IDS SORTED BY THEIR HIGHEST TOTAL PROBABILITY
    print("DEBUG: clustID_highestTotProb_dict = ", clustID_highestTotProb_dict)
    sorted_clustIDs_list = [x[0] for x in sorted(list(clustID_highestTotProb_dict.items()), key=itemgetter(1), reverse=True)]
    print("DEBUG: sorted_clustIDs_list=", sorted_clustIDs_list)
    for clustID in clustIDs_set:
        if clustID not in sorted_clustIDs_list:
            sorted_clustIDs_list.append(clustID)
    print("DEBUG: after addition of missing Clusters sorted_clustIDs_list=", sorted_clustIDs_list)
    
    
    CONTINUE_1ST_GROUPING = True # if no new assignments are saved, set to False to stop the iterations
    while CONTINUE_1ST_GROUPING:
        CONTINUE_1ST_GROUPING = False
        for clustID in sorted_clustIDs_list:
            TresonIndex_set = set([x[4] for x in aatypeResonpairMatchesTuple_list if x[7]==clustID])    # set with all different peaks (TOCSY reson indices) of the current Cluster ID
            if clustIDs_list.count(clustID) > 1 and len(TresonIndex_set) > 1:    # if there are 2 or more possible assignments and at least 2 different peaks within the same Cluster ID
                nucleusAssignmentsforCurrentClustID_list = [x for x in aatypeResonpairMatchesTuple_list if x[7]==clustID]    # save all the members of this Cluster ID
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
                            print("1st GROUPING EXAMPLE: nucleusAssignmentsforCurrentClustID_list=", nucleusAssignmentsforCurrentClustID_list)
                            correct_C_H_resonpair_matches_list.append(assignment)
                            used_TOCSY_reson_indices_list.append(TresonIndex)
                            used_C_H_pairs_list.append(C+"_"+H)
                            CONTINUE_1ST_GROUPING = True # since a new assignment has been saved, iterate the cycle once more
                            # CLEAN aatypeResonpairMatchesTuple_list from other possible assignments of the same clustID that do not contain this C type
                            for x in aatypeResonpairMatchesTuple_list:
                                # CONDITIONS:
                                # 1) the possible assignment must belong to the current ClusterID
                                # 2) if the Carbon type in the assignment is not the CARBON_OF_THIS_CLUSTER, remove the assignment
                                if x[7] == clustID and x[1] != CARBON_OF_THIS_CLUSTER:
                                    aatypeResonpairMatchesTuple_list.remove(x)
        print("DEBUG: after 1st grouping correct_C_H_resonpair_matches_list=", correct_C_H_resonpair_matches_list)
    
    for clustID in sorted_clustIDs_list:
        TresonIndex_set = set([x[4] for x in aatypeResonpairMatchesTuple_list if x[7]==clustID])    # set with all different peaks (TOCSY reson indices) of the current Cluster ID
        # NOW USE THE Treson indices WITH MULTIPLE TYPES OF CARBONS
        if clustIDs_list.count(clustID) > 1 and len(TresonIndex_set) > 1:    # if there are 2 or more possible assignments and at least 2 different peaks within the same Cluster ID
            nucleusAssignmentsforCurrentClustID_list = [x for x in aatypeResonpairMatchesTuple_list if x[7]==clustID]    # save all the members of this Cluster ID
            # Recall that nucleusAssignmentsforCurrentClustID_list is already sorted by lines 1,2,3, BUT WE MUST ...
            nucleusAssignmentsforCurrentClustID_list.sort(key=itemgetter(3), reverse=True)  # ... SORT BY THE PREDICTION PROBABILITY
            
            # Make a list of the carbon types in each Treson index of this Cluster and its maximum probability
            C_totalProb_dict = {}   # the total probability of each C type in this Cluster. Use this to decide the assignment
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
            print("DEBUG: clustID=", clustID, "sorted_C_totalProb_list=", sorted_C_totalProb_list)
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
                        print("2nd GROUPING EXAMPLE: nucleusAssignmentsforCurrentClustID_list=", nucleusAssignmentsforCurrentClustID_list)
                        correct_C_H_resonpair_matches_list.append(assignment)
                        used_TOCSY_reson_indices_list.append(TresonIndex)
                        used_C_H_pairs_list.append(C+"_"+H)
                        # CLEAN aatypeResonpairMatchesTuple_list from other possible assignments of the same clustID that do not contain this C type
                        for x in aatypeResonpairMatchesTuple_list:
                            # CONDITIONS:
                            # 1) the possible assignment must belong to the current ClusterID
                            # 2) if the Carbon type in the assignment is not the CARBON_OF_THIS_CLUSTER, remove the assignment
                            if x[7] == clustID and x[1] != CARBON_OF_THIS_CLUSTER:
                                aatypeResonpairMatchesTuple_list.remove(x)
        print("DEBUG: after 2nd grouping correct_C_H_resonpair_matches_list=", correct_C_H_resonpair_matches_list)
    
    print("DEBUG: after grouping all clustIDs correct_C_H_resonpair_matches_list=", correct_C_H_resonpair_matches_list)
    
    
    # NOW TAKE CARE OF UNGROUPED CARBONS, SPECIFICALLY CLUSTERS WITH ONLY ONE Treson index AND ONLY ONE POSSIBLE ASSIGNMENT
    for clustID in sorted_clustIDs_list:
        TresonIndex_set = set([x[4] for x in aatypeResonpairMatchesTuple_list if x[7]==clustID])    # set with all different peaks (TOCSY reson indices) of the current Cluster ID
        # NOW USE THE Treson indices WITH MULTIPLE TYPES OF CARBONS
        if clustIDs_list.count(clustID) == 1 and len(TresonIndex_set) == 1:    # if there is only 1 possible assignment and only 1 peak within the same Cluster ID
            nucleusAssignmentsforCurrentClustID_list = [x for x in aatypeResonpairMatchesTuple_list if x[7]==clustID]    # save the single member of this Cluster ID
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
                print("3rd GROUPING EXAMPLE: nucleusAssignmentsforCurrentClustID_list=", nucleusAssignmentsforCurrentClustID_list)
                correct_C_H_resonpair_matches_list.append(assignment)
                used_TOCSY_reson_indices_list.append(TresonIndex)
                used_C_H_pairs_list.append(C+"_"+H)
                # CLEAN aatypeResonpairMatchesTuple_list from other possible assignments of the same clustID that do not contain this C type
                for x in aatypeResonpairMatchesTuple_list:
                    # CONDITIONS:
                    # 1) the possible assignment must belong to the current ClusterID
                    # 2) if the Carbon type in the assignment is not the CARBON_OF_THIS_CLUSTER, remove the assignment
                    if x[7] == clustID and x[1] != CARBON_OF_THIS_CLUSTER:
                        aatypeResonpairMatchesTuple_list.remove(x)
    print("DEBUG: after Clusters with a single TRI and single assignment correct_C_H_resonpair_matches_list=", correct_C_H_resonpair_matches_list)
    
    # THEN, TAKE CARE OF UNGROUPED CARBONS, SPECIFICALLY CLUSTERS WITH ONLY ONE Treson index BUT MULTIPLE POSSIBLE ASSIGNMENTS
    for clustID in sorted_clustIDs_list:
        TresonIndex_set = set([x[4] for x in aatypeResonpairMatchesTuple_list if x[7]==clustID])    # set with all different peaks (TOCSY reson indices) of the current Cluster ID
        # NOW USE THE Treson indices WITH MULTIPLE TYPES OF CARBONS
        if clustIDs_list.count(clustID) > 1 and len(TresonIndex_set) == 1:    # if there is only 1 possible assignment and only 1 peak within the same Cluster ID
            nucleusAssignmentsforCurrentClustID_list = [x for x in aatypeResonpairMatchesTuple_list if x[7]==clustID]    # save the single member of this Cluster ID
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
                print("4th GROUPING EXAMPLE: nucleusAssignmentsforCurrentClustID_list=", nucleusAssignmentsforCurrentClustID_list)
                correct_C_H_resonpair_matches_list.append(assignment)
                used_TOCSY_reson_indices_list.append(TresonIndex)
                used_C_H_pairs_list.append(C+"_"+H)
                # CLEAN aatypeResonpairMatchesTuple_list from other possible assignments of the same clustID that do not contain this C type
                for x in aatypeResonpairMatchesTuple_list:
                    # CONDITIONS:
                    # 1) the possible assignment must belong to the current ClusterID
                    # 2) if the Carbon type in the assignment is not the CARBON_OF_THIS_CLUSTER, remove the assignment
                    if x[7] == clustID and x[1] != CARBON_OF_THIS_CLUSTER:
                        aatypeResonpairMatchesTuple_list.remove(x)
    print("DEBUG: after Clusters with a single TRI but multiple assignments correct_C_H_resonpair_matches_list=", correct_C_H_resonpair_matches_list)
    
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
    
    print("DEBUG: C_H_occupancy_dict = ", C_H_occupancy_dict)
    print("DEBUG: C_H_matches_dict = ", C_H_matches_dict)
    sorted_C_H_occupancyTuple_list = sorted(list(C_H_occupancy_dict.items()), key=itemgetter(1))  # the C_H_occupancyTuple tuples sorted by occupancy, so ambiguous assignments (higher occupancy) are left for the end
    print("DEBUG: sorted_C_H_occupancyTuple_list=", sorted_C_H_occupancyTuple_list)
    for C_H_occupancy_tuple in sorted_C_H_occupancyTuple_list:
        C_H_pair = C_H_occupancy_tuple[0]
        C = C_H_pair.split("_")[0]
        H = C_H_pair.split("_")[1]
        for match in C_H_matches_dict[C_H_pair]:   # we don't need to sort here the alternative C-H resonances according to the weighted average prob,
                                # because it has been done already at "possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list.sort(key=itemgetter(0,1,2,3), reverse=True)"
            print("DEBUG: match = ", match)
            print("DEBUG: correct_C_H_resonpair_matches_list=", correct_C_H_resonpair_matches_list)
            TOCSY_reson_index = match[4]
            # CONDITIONS:
            # 1) This Treson index has not been assigned yet
            # 2) This condition was useful before the incorporation of carbon grouping. I.e. it checked whether HB2-CB and HB3-CB had the same carbon resonance (+- 0.3)
            if not TOCSY_reson_index in used_TOCSY_reson_indices_list and do_carbons_match(match, correct_C_H_resonpair_matches_list):
                print("5th GROUPING EXAMPLE: nucleusAssignmentsforCurrentClustID_list=", nucleusAssignmentsforCurrentClustID_list)
                correct_C_H_resonpair_matches_list.append(match)
                used_TOCSY_reson_indices_list.append(TOCSY_reson_index)
                break   # move to the next C-H resonance pair
    
    print("DEBUG: final correct_C_H_resonpair_matches_list=",correct_C_H_resonpair_matches_list)
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
    nucleusAssignmentsforCurrentTRI_list = [x for x in possible_aatype_prob_C_H_resonpair_TOCSYindex_assignID_list_list
                                            if x[4]==TRI] # save all the members of this Cluster ID
    #TresonIndex_set = set([x[4] for x in possible_aatype_prob_C_H_resonpair_TOCSYindex_assignID_list_list if x[7]==clustID])    # set with all different peaks (TOCSY reson indices) of the current Cluster ID
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
    
    # print Assignment_Tree.get_ascii(show_internal=True, compact=False, attributes=["name", "dist", "occupancy", "numOfResonances"])
    # print Assignment_Tree.get_ascii(show_internal=True, compact=False)
    if number_of_new_leaves > 0:
        return (Assignment_Tree, True)
    else:
        return (Assignment_Tree, False)


def get_prediction_score(prob_clustID_tuple_list, pscore_mode=1):
    """
    For each C-group (aka clustID), let the individual probabilities are prob1 and prob2, this function will calculate a
    consensus probability in one of the following ways (determined by 'cons_mode' variable), and at the end return the
    product if these consensus probabilities.
    cons_mode:   1: average;
            2: sqrt(prob1*prob2)    ; geometric mean
            3: (prob2-prob1)/(log(prob2)-log(prob1))    ; logarithmic mean
            4: (prob1*prob2)/((prob1+prob2)/2.) ;composite average
            5: prob1*prob2    ; product of probabilities
    """
    
    prediction_score = 1 # product of the consensus probabilites of each C-group
    Cgroup_set = set([t[1] for t in prob_clustID_tuple_list])
    for Cgroup in Cgroup_set:
        prob_set = set([t[0] for t in prob_clustID_tuple_list if t[1]==Cgroup]) # get the peak probabilities of this C-group
        prob_list = list(prob_set)
        prediction_score *= get_Cgroup_consensus_probability(prob_list, cons_mode=pscore_mode)
    
    print("DEBUG: prediction_score=", prediction_score)
    return prediction_score


def select_correct_H_C_resonpairs2(aatypeResonpairMatchesTuple_list,
                                   args,
                                   get_aa_type=False):
    """
        NEW FUNCTION  to match the TOCSY resonance pairs to as much as possible C-H resonance pairs of a PARTICULAR aa type.
        This function selects the best C-H type assignment combination based on the product of probabilities of the
        individual C-H assignments.
        *** This function is still used for aa type prediction ***
        
        ARGUMETS:
        aatypeResonpairMatchesTuple_list:   list of lists of the form [aa_type, C, H, weighted average probability, selected TOCSY reson index, H resonance, C resonance, ClusterID] where aa_type is the same in the whole list (e.g. "PRO")
        
        RETURNS:
        correct_C_H_resonpair_matches_list: list of list of the form [aa_type, C, H, weighted average probability, selected TOCSY reson index], where aa_type is the same in the whole list (e.g. "PRO"). This list contains
                                            the final, uniq nucleus type assignment for the current amino acid. E.g.
        [['MET', 'CA', 'HA', 0.01993300618212495, 2, 3.843, 55.13, 1], ['MET', 'CB', 'HB2', 0.012634347058708313, 4, 1.266, 31.504, 2],
        ['MET', 'CG', 'HG2', 0.02281033340649995, 1, 1.508, 31.481, 2], ['MET', 'CB', 'HB3', 0.009737867955406978, 3, 1.911, 31.223, 3],
        ['MET', 'CG', 'HG3', 0.01607381733870664, 5, 0.403, 31.186, 3]]
    """ 

    global aa_carbonBondedHydrogensDict_dict

    possible_aatype_prob_C_H_resonpair_TOCSYindex_assignID_list_list = []   # same as possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list but with an assignment extra index at the end
    for i,x in enumerate(aatypeResonpairMatchesTuple_list):
        possible_aatype_prob_C_H_resonpair_TOCSYindex_assignID_list_list.append([x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], i])
    
    # print("DEBUG: possible_aatype_prob_C_H_resonpair_TOCSYindex_assignID_list_list=", possible_aatype_prob_C_H_resonpair_TOCSYindex_assignID_list_list)
    clustIDs_list = [x[7] for x in possible_aatype_prob_C_H_resonpair_TOCSYindex_assignID_list_list]    # list of Cluster ID of each possible assignment
    clustIDs_set = set(clustIDs_list)
    TRI_set = set([x[4] for x in possible_aatype_prob_C_H_resonpair_TOCSYindex_assignID_list_list])
    TRI_num = len(TRI_set)
    # print("DEBUG: TRI_set =", TRI_set)
    
    print("Building Tree starting from amino acid index...")
    aa_type = possible_aatype_prob_C_H_resonpair_TOCSYindex_assignID_list_list[0][0]
    expand_tree = True
    Assignment_Tree = Tree()
    Root = Assignment_Tree.get_tree_root()
    Root.add_features(name=aa_type, Ctype="None", Htype="None")
    level = 0
    sys.stdout.write("Expanding tree from level ")
    for TRI in TRI_set:
        sys.stdout.write(str(level)+" ")
        sys.stdout.flush()
        Assignment_Tree, expand_tree = populate_leaves_for_correct_H_C_resonpairs2(Assignment_Tree,
                                                                                   possible_aatype_prob_C_H_resonpair_TOCSYindex_assignID_list_list,
                                                                                   TRI)
        if expand_tree == False:
            break
        level += 1
        if level == TRI_num:
            break
    # if the maximum number of Cgroups for this aa_Type has been reached, discard combinations with less that Cluster_num peaks
    # WARNING: might not work for aromatic aa types like PHE where the CZ-HZ is included in aa_carbonBondedHydrogensDict_dict
    # TODO: dectivating the ERROR-check for aa type prediction is wrong, but in any case this function will be replaced
    # TODO: by NNs.
    if not get_aa_type and \
            level+1 < TRI_num and \
            level+1 < len(list(aa_carbonBondedHydrogensDict_dict[aa_type].keys())):
        print(bcolors.FAIL + "\nERROR: NOT ALL C-GROUPS WERE USED IN THE ASSIGNMENT OF " + aa_type + bcolors.ENDC)
        sys.exit(1)
      
    
    print("\nSaving chains from Tree...")
    all_chainScore_list = []
    for leaf in Assignment_Tree.get_leaves():
        chain = []
        prob_product = leaf.Prob
        prob_clustID_tuple_list = [(leaf.Prob, leaf.clustID)]
        assignID = int(leaf.name)
        chain.append(possible_aatype_prob_C_H_resonpair_TOCSYindex_assignID_list_list[assignID])
        # print("DEBUG: leaf.name=", leaf.name)
        for ancestor in leaf.get_ancestors():
            # NOTE: instead of 'if ancestor == Root: break' , I check if it's name is an int (only Root has name string).
            # print("DEBUG: ancestor.name=", ancestor.name)
            try:
                assignID = int(ancestor.name)
            except ValueError:
                continue
            chain.append(possible_aatype_prob_C_H_resonpair_TOCSYindex_assignID_list_list[assignID])
            prob_product *= ancestor.Prob
            prob_clustID_tuple_list.append((ancestor.Prob, ancestor.clustID))
        if args.CONSENSUS_CGROUP_PROB_MODE == 0:    # multiply all 2D-hist probabilities (including the 2 C-H of methylenes) the get the
                                                    # aa type score
            chain.append(prob_product)
        else:   # For each C-group with 2 peaks (methylene) calculate one probability value according to the given
                # CONSENSUS_CGROUP_PROB_MODE (2 by default, namely sqrt(prob1*prob2)). For each C-group with 1 peak keep the original
                # 2D-hist probability. Then multiply the probabilities of all C-groups to get the aa type score.
            pred_score = get_prediction_score(prob_clustID_tuple_list,
                                     pscore_mode=args.CONSENSUS_CGROUP_PROB_MODE)
            if pred_score == None:  # Carbon group contains more than 2 peaks!!!
                print(bcolors.FAIL + "ERROR: this Carbon group contains more than two peaks!!!!" + bcolors.ENDC)
                clustID = prob_clustID_tuple_list[0][1]
                print([p for p in possible_aatype_prob_C_H_resonpair_TOCSYindex_assignID_list_list if p[7]==clustID])
                sys.exit(1)
            chain.append(pred_score)
        # print("DEBUG: saving chain=", chain)
        all_chainScore_list.append(tuple(chain))
        del chain
        del ancestor
        del leaf
    del Assignment_Tree
    gc.collect()
    
    
    correct_chainScore_list = []
    print("DEBUG: all_chainScore_list=", all_chainScore_list)
    carbon_exceptions_list = ["CD1", "CD2"] # carbons that may belong to the same group but have different names
    aatype_Ctype_exceptions_dict = {}   # not in use yet
    aatype_Ctype_exceptions_dict['GLY'] = ["CA"]
    for x in all_chainScore_list:
        DISCARD = False # discard this chain 'x'
        # discard the prediction if it does not contain all the TRIs
        if len(x) < len(TRI_set)+1 :
            DISCARD = True
        # check the C-groups
        Ctype2ClustID_dict = {}
        for clustID in clustIDs_set:
            Ctype_list = list(set([a[1] for a in x[:-1] if a[7]==clustID]))    # C-type(s) of this Cluster
            if len(Ctype_list) > 1 and not Ctype_list[0] in carbon_exceptions_list and not Ctype_list[1] in carbon_exceptions_list:
                DISCARD = True
                break
        # check that each C-group contains a different C-type
        Ctype_set = set([a[1] for a in x[:-1]])
        for Ctype in Ctype_set:
            if len(set([a[7] for a in x[:-1] if a[1]==Ctype and not (a[0]=='GLY' and a[1]=='CA')])) > 1:
                DISCARD = True
                break
        # check that each C-type belongs to only one C-group
        Ctype_set = set([a[1] for a in x[:-1]])
        for Ctype in Ctype_set:
            Cgroups_list = list(set([a[-2] for a in x[:-1] if a[1]==Ctype]))    # C-groups of this Cluster
            if len(Cgroups_list) > 1:
                DISCARD = True
                break
        # check for H-type duplication
        Htype_set = set([a[2] for a in x[:-1]])
        Htype_list = [a[2] for a in x[:-1]]
        if len(Htype_set) !=  len(Htype_list):
            DISCARD = True
        
        if DISCARD:
            continue
        else:
            correct_chainScore_list.append(x)
    
    
    correct_chainScore_list.sort(key=itemgetter(len(TRI_set)), reverse=True) # sort probability by descending order; the highest probability combination must be 1st
    print("DEBUG: all possible predictions correct_chainScore_list=", correct_chainScore_list)
    
    print("\nDEBUG: Best Prediction:")
    try:
        for x in correct_chainScore_list[0][:-1]:
            if x[3] <0:
                print("DEBUG: Best Prediction contains C-based probability!!!")
            print(x)
    except IndexError:
        print("DEBUG: correct_chainScore_list=", correct_chainScore_list)
        #sys.exit(1)    # temporarily inactive
        if get_aa_type:
            return [], aa_type, 0.0
        else:
            return []
    if get_aa_type:
        return correct_chainScore_list[0][:-1], aa_type, correct_chainScore_list[0][-1]
    else:
        return correct_chainScore_list[0][:-1]


def find_peakID(Hreson, Creson, aatypeResonpairMatchesTuple_list):
    """
        FUNCTION to recover the original Peak ID.
        ARGS:
        Hreson: float
        Creson: float
        aatypeResonpairMatchesTuple_list:   e.g.
                        [
                        ['LYS', 'CG', 'HG3', 1.064222338e-08, 2.14, 34.786, 2, 3, 'TOCSY'],
                        ['LYS', 'CG', 'HG3', 2.199171449e-23, 2.487, 30.946, 3, 2, 'TOCSY'],
                        ['LYS', 'CG', 'HG2', 1.064222338e-08, 2.14, 34.786, 2, 3, 'TOCSY'],
                        ['LYS', 'CG', 'HG2', 2.199171449e-23, 2.487, 30.946, 3, 2, 'TOCSY'],
                        ['LYS', 'CE', 'HE3', 9.582673271e-06, 2.14, 34.786, 2, 3, 'TOCSY'],
                        ['LYS', 'CE', 'HE2', 9.582673271e-06, 2.14, 34.786, 2, 3, 'TOCSY'],
                        ['LYS', 'CD', 'HD3', 5.104388159e-07, 2.14, 34.786, 2, 3, 'TOCSY'],
                        ['LYS', 'CD', 'HD3', 2.878109776e-08, 2.487, 30.946, 3, 2, 'TOCSY'],
                        ['LYS', 'CD', 'HD2', 5.104388159e-07, 2.14, 34.786, 2, 3, 'TOCSY'],
                        ['LYS', 'CD', 'HD2', 2.878109776e-08, 2.487, 30.946, 3, 2, 'TOCSY'],
                        ['LYS', 'CB', 'HB3', 0.0001074583495, 2.14, 34.786, 2, 3, 'TOCSY'],
                        ['LYS', 'CB', 'HB3', 8.821068844e-06, 2.487, 30.946, 3, 2, 'TOCSY'],
                        ['LYS', 'CB', 'HB2', 0.0001074583495, 2.14, 34.786, 2, 3, 'TOCSY'],
                        ['LYS', 'CB', 'HB2', 8.821068844e-06, 2.487, 30.946, 3, 2, 'TOCSY'],
                        ['LYS', 'CA', 'HA', 0.00130609899, 4.194, 56.022, 1, 1, 'TOCSY']
                        ]
        RETURNS:
        peak_ID:    int >= 0
    """
    
    peak_ID_set = set([p[6] for p in aatypeResonpairMatchesTuple_list if p[5]==Hreson and p[6]==Creson])
    if len(peak_ID_set) > 1:
        print("ERROR: more than one peak IDs correspond to resonances ", Hreson, Creson)
        print("aatypeResonpairMatchesTuple_list=", aatypeResonpairMatchesTuple_list)
        sys.exit(1)
    elif len(peak_ID_set) == 1:
        peak_ID = list(peak_ID_set)[0]
        return peak_ID


def select_correct_H_C_resonpairs4(aatypeResonpairMatchesTuple_list, args, iteration=None):
    """
        NEW FUNCTION  to match the C-groups to as much as possible C-types of a PARTICULAR aatype. This function selects the
        best C-group assignment combination based on the product of probabilities of the individual C-group assignments. The
        difference here is that we work with C-groups, not with individual peaks! Tested only for TOCSY assignments! Does not
        support Carbon only based prediction (look at eliminate_orphan_Ctype_predictions).
        
        ARGUMETS:
        aatypeResonpairMatchesTuple_list:   list of lists of the form [aa_type, C, H, weighted average probability, peak_ID, H resonance, C resonance, ClusterID, spectrum type] where aa_type is the same in the whole list (e.g. "PRO")
        iteration:      used only in filter_for_thresholds, to eventually control the thresholds
        
        RETURNS:
        correct_C_H_resonpair_matches_list: list of list of the form [aa_type, C, H, weighted average probability, selected TOCSY reson index], where aa_type is the same in the whole list (e.g. "PRO"). This list contains
                                            the final, uniq nucleus type assignment for the current amino acid. E.g.
        [['MET', 'CA', 'HA', 0.01993300618212495, 2, 3.843, 55.13, 1], ['MET', 'CB', 'HB2', 0.012634347058708313, 4, 1.266, 31.504, 2],
        ['MET', 'CG', 'HG2', 0.02281033340649995, 1, 1.508, 31.481, 2], ['MET', 'CB', 'HB3', 0.009737867955406978, 3, 1.911, 31.223, 3],
        ['MET', 'CG', 'HG3', 0.01607381733870664, 5, 0.403, 31.186, 3]]
    """
    
    global aa_carbonBondedHydrogensDict_dict
    
    if iteration == None:   # if this is TOCSY assignments, append the spectrum type at the end of each peak and convert it to NOESY format
                            # namely from [aa_type, C, H, weighted average probability, peak_ID, H resonance, C resonance, ClusterID, spectrum type]
                            # to          [aa_type, C, H, weighted average probability, H resonance, C resonance, peak_ID, ClusterID, spectrum type] 
        aatypeResonpairMatchesTuple_list = [ [p[0], p[1], p[2], p[3], p[5], p[6], p[4], p[7], "TOCSY"]
                                             for p in aatypeResonpairMatchesTuple_list ]
    
    print("DEBUG select_correct_H_C_resonpairs4: aatypeResonpairMatchesTuple_list=", aatypeResonpairMatchesTuple_list)
    # Eliminate C-type predictions that do not contain both peaks of the same C-group
    filtered_aatypeResonpairMatchesTuple_list = eliminate_orphan_Ctype_predictions(aatypeResonpairMatchesTuple_list)
    Cname_set = set([x[1] for x in filtered_aatypeResonpairMatchesTuple_list if x[3]>0]) # exclude Carbon-only based assignments
    Cgroups_set = set([x[7] for x in filtered_aatypeResonpairMatchesTuple_list if x[3]>0]) # exclude Carbon-only based assignments
    aa_type =  filtered_aatypeResonpairMatchesTuple_list[0][0]
    possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_list_list = []   # [[aa_type, Carbon, tuple(Hlist_of_this_Cgroup), Cgroup_probability, tuple(HCresonpair_set_of_this_Cgroup), Cgroup], ...]
    for Cgroup in Cgroups_set:
        for Cname in Cname_set:
            print("DEBUG: Cgroup=", Cgroup, "Cname=", Cname)
            # matching assingments can be 2 or 4, because each peak can be e.g. CB-HB2 or CB-HB3 with equal probability. So it is fair to just average the probabilities.
            matching_assignments = list(set([tuple(x) for x in filtered_aatypeResonpairMatchesTuple_list
                                             if x[1]==Cname and x[7]==Cgroup and x[3]>0]))   # exclude Carbon based only assignments (they have negative probability)
            if len(matching_assignments) == 0:  # if there are not assignments of this C-type for this C-group, continue to the next C-type
                continue
            Hlist_of_this_Cgroup = list(set([x[2] for x in matching_assignments]))  # it doesn't matter if the order will be changed, both protons are equivalent
            Hlist_of_this_Cgroup.sort()
            HCresonpair_spectrumType_set_of_this_Cgroup = tuple(set([(str(x[5])+"_"+str(x[6]), x[8]) for x in matching_assignments]))
            Cgroup_probability = get_Cgroup_consensus_probability([x[3] for x in matching_assignments
                                                                   if x[1] == Cname and x[7] == Cgroup],
                                                                  cons_mode=args.CONSENSUS_CGROUP_PROB_MODE)
            if len(HCresonpair_spectrumType_set_of_this_Cgroup) == 1:    # if the C-group contains only one peak
                Hlist_of_this_Cgroup = [Hlist_of_this_Cgroup[0]]    # keep only one of the methylene proton names
            HCresonpair_set_of_this_Cgroup = [t[0] for t in HCresonpair_spectrumType_set_of_this_Cgroup]
            HCresonpair_set_of_this_Cgroup.sort()   # always the lowest proton resonance goes 1st
            HCresonpair_set_of_this_Cgroup = tuple(HCresonpair_set_of_this_Cgroup)  # convert the list to tuple
            spectrum_combos = tuple([t[1] for t in HCresonpair_spectrumType_set_of_this_Cgroup])
            Cgroup_assignment = [aa_type, Cname, tuple(Hlist_of_this_Cgroup), Cgroup_probability, HCresonpair_set_of_this_Cgroup, Cgroup, spectrum_combos]
            possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_list_list.append(Cgroup_assignment)
    print("DEBUG: possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_list_list=",
          possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_list_list)
    print("DEBUG: Cgroups_set=", Cgroups_set)
    if (iteration==1 and args.WITHIN_CGROUP_THRESHOLD_ITER1==0 and args.BETWEEN_CGROUP_THRESHOLD_ITER1==0) or (
        iteration==2 and args.WITHIN_CGROUP_THRESHOLD_ITER2==0 and args.BETWEEN_CGROUP_THRESHOLD_ITER2==0) or (
        iteration==3 and args.WITHIN_CGROUP_THRESHOLD_ITER3==0 and args.BETWEEN_CGROUP_THRESHOLD_ITER3==0) or (
        iteration==4 and args.WITHIN_CGROUP_THRESHOLD_ITER4==0 and args.BETWEEN_CGROUP_THRESHOLD_ITER4==0) or (
        iteration==5 and args.WITHIN_CGROUP_THRESHOLD_ITER5==0 and args.BETWEEN_CGROUP_THRESHOLD_ITER5==0) or (
        iteration == None): # or if this is TOCSY assignment (iteration=None)
        print("DEBUG: no filtering for WITHIN_CGROUP_THRESHOLD and BETWEEN_CGROUP_THRESHOLD")
        filtered_possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_list_list = possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_list_list
    elif iteration != None:
        print("DEBUG: filtering for WITHIN_CGROUP_THRESHOLD and BETWEEN_CGROUP_THRESHOLD")
        filtered_possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_list_list = filter_for_thresholds(possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_list_list, iteration, args)
    if len(filtered_possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_list_list) == 0:
        return []
    
    possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_assignID_list_list = []   # same as filtered_possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_list_list but with an assignment extra index at the end
    for i,x in enumerate(filtered_possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_list_list):
        possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_assignID_list_list.append([x[0], x[1], x[2], x[3], x[4], x[5], i, x[6]])  # x[6] is the spectrum type
    print("DEBUG: possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_assignID_list_list=", possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_assignID_list_list)
    
    Cgroups_list = [x[5] for x in possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_assignID_list_list]    # list of Cluster ID of each possible assignment
    Cgroups_set = set(Cgroups_list)
    filtered_Cnames_list = [x[1] for x in possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_assignID_list_list]    # C-names after filtering
    filtered_Cnames_set = set(filtered_Cnames_list)
    filtered_Cnames_num = len(filtered_Cnames_set)
    if filtered_Cnames_num > len(Cgroups_set):
        if aa_type == "LEU" and 'CD1' in filtered_Cnames_set and 'CD2' in filtered_Cnames_set:
            filtered_Cnames_num = len(Cgroups_set)
        elif aa_type == "VAL" and 'CG1' in filtered_Cnames_set and 'CG2' in filtered_Cnames_set:
            filtered_Cnames_num = len(Cgroups_set)
        # else:
        #     print "ERROR: more C-names than C-groups, no assignment combinations will be made!"
        #     print "Cgroups_set=", Cgroups_set, "filtered_Cnames_set=", filtered_Cnames_set
        #     sys.exit(1)
    print("DEBUG: filtered_Cnames_set=", filtered_Cnames_set, "filtered_Cnames_num=", filtered_Cnames_num)
    print("Building Tree starting from amino acid index...")
    all_chainScore_list = []
    sys.stdout.write("Expanding tree from level ")
    # possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_assignID_list_list is a list [[aa_type, Carbon, tuple(Hlist_of_this_Cgroup), Cgroup_probability, tuple(HCresonpair_set_of_this_Cgroup), Cgroup, assignID, spectrum type tuple], ...]
    if len(Cgroups_set) >= filtered_Cnames_num:
        combs = combinations(Cgroups_set, filtered_Cnames_num)
    elif len(Cgroups_set) < filtered_Cnames_num:    # if the missing Cnames are more than the avialable Cgroups, use all Cgroups!
        combs = [tuple(Cgroups_set)]
    for Cluster_combination in combs:  # after filtering each C-group is contained only once in possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_assignID_list_list
        print("\nDEBUG: Cluster_combination=", Cluster_combination)
        # To make sure that all the Cgroups and possible Cnames are used, for each Cluster combination build a tree each time starting from a different Cgroup.
        # That will produce redundant chains but will not miss any Cgroup combination.
        for permutation in permutations(Cluster_combination):
            print("DEBUG: building new tree using ClusterID permutation:", permutation)
            expand_tree = True
            Assignment_Tree = Tree()
            Root = Assignment_Tree.get_tree_root()
            Root.add_features(name=aa_type, Ctype="None", Htype="None")
            level = 0
            for Cluster in permutation:
                sys.stdout.write(str(level)+" ")
                sys.stdout.flush()
                Assignment_Tree, expand_tree = populate_leaves_for_correct_H_C_resonpairs3(Assignment_Tree, possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_assignID_list_list, Cluster)
                print("DEBUG: expand_tree=", expand_tree)
                if expand_tree == False:
                    break
                print("DEBUG: added Cgroup", Cluster, "to the Tree.")
                level += 1
                #if level == CH_pair_num:
                #    break
            # if level < filtered_Cnames_num: # DANGEROUS: discard combinations with less than Cgroups_num peaks
            #     print "ERROR: NOT ALL C-GROUPS WERE USED IN THE ASSIGNMENT OF ", aa_type
            #     sys.exit(1)
                
            print("\nSaving chains from Tree...")
            for leaf in Assignment_Tree.get_leaves():
                chain = []
                prob_product = leaf.Prob
                assignID = int(leaf.name)
                chain.append(possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_assignID_list_list[assignID])
                print("DEBUG: leaf.name=", leaf.name, "leaf.Prob=", leaf.Prob)
                for ancestor in leaf.get_ancestors():
                    # NOTE: instead of 'if ancestor == Root: break' , I check if it's name is an int (only Root has name string).
                    print("DEBUG: ancestor.name=", ancestor.name, type(ancestor.name))
                    try:
                        assignID = int(ancestor.name)
                    except ValueError:
                        continue
                    print("DEBUG: ancestor.Prob=", ancestor.Prob)
                    chain.append(possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_assignID_list_list[assignID])
                    prob_product *= ancestor.Prob
                if iteration is not None and iteration > 1:
                    if len(chain) < filtered_Cnames_num: # DANGEROUS: discard combinations with less than Cgroups_num peaks: THIS DOES NOT APPLY IN TOCSY or
                                                    # NOESY i->i+1 ASSIGNMENT BECAUSE WE ARE NEVER SURE WHETHER ALL THE i->i+1 CGROUPS CORRESPOND TO RESIDUE i
                        print("DEBUG: discarding short chain=", chain)
                        continue
                chain.append(prob_product)
                print("DEBUG: saving chain=", chain)
                all_chainScore_list.append(tuple(chain))
                del chain
                del ancestor
                del leaf
            del Assignment_Tree
            gc.collect()
    if iteration == 1:  # keep also combinations consisting of only 1 peak
        for Cluster in Cgroups_set:
            print("DEBUG: saving single length chains: ", [(a, a[3]) for a in possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_assignID_list_list if a[5]==Cluster])
            all_chainScore_list.extend([(a, a[3]) for a in possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_assignID_list_list if a[5]==Cluster])
    correct_chainScore_list = []
    print("DEBUG: all_chainScore_list=", all_chainScore_list)
    aatype_Ctype_exceptions_dict = {}   # not in use yet
    aatype_Ctype_exceptions_dict['GLY'] = ["CA"]
    for x in all_chainScore_list:
        DISCARD = False
        # # discard the prediction if it does not contain all the Clusters
        # if len(x) != missing_Cname_num+1:
        #     DISCARD = True
        # check that each C-group is assigned to only one C-type
        Ctype2ClustID_dict = {}
        for Cgroup in Cgroups_set:
            Ctype_list = list(set([a[1] for a in x[:-1] if a[5]==Cgroup]))    # C-type(s) of this Cluster in the current prediction (x)
            # if two different C-types have been assigned to the peaks of this C-group, discart the prediction
            if len(Ctype_list) > 1: # check for carbons exceptions that may belong to the same C-group but have different names
                if aa_type == "LEU" and Ctype_list[0] in ["CD1", "CD2"] and Ctype_list[1] in ["CD1", "CD2"]:
                    pass
                elif aa_type == "VAL" and Ctype_list[0] in ["CG1", "CG2"] and Ctype_list[1] in ["CG1", "CG2"]:
                    pass
                else:
                    DISCARD = True
                    break
        # check that each C-group contains a different C-type
        Ctype_set = set([a[1] for a in x[:-1]])
        assigned_Ctypes_list = []   # list with the C-types encountered already
        for Cgroup in Cgroups_set:
            Ctypes_of_this_Cgroup_list = [a[1] for a in x[:-1] if a[5]==Cgroup]
            if len(Ctypes_of_this_Cgroup_list)==1 and Ctypes_of_this_Cgroup_list[0] in assigned_Ctypes_list: # if this C-type has already been assigned to another C-group, discard the assignment
                DISCARD = True
                break
            # if there are two C-types but they are within the exceptions, save the assigned C-types
            elif len(Ctypes_of_this_Cgroup_list)==2:
                if aa_type == "LEU" and Ctypes_of_this_Cgroup_list[0] in ["CD1", "CD2"] and Ctypes_of_this_Cgroup_list[1] in ["CD1", "CD2"]:
                    assigned_Ctypes_list.append(Ctypes_of_this_Cgroup_list[0])  # save the encountered C-type
                    assigned_Ctypes_list.append(Ctypes_of_this_Cgroup_list[1])  # save the encountered C-type
                elif aa_type == "VAL" and Ctypes_of_this_Cgroup_list[0] in ["CG1", "CG2"] and Ctypes_of_this_Cgroup_list[1] in ["CG1", "CG2"]:
                    assigned_Ctypes_list.append(Ctypes_of_this_Cgroup_list[0])  # save the encountered C-type
                    assigned_Ctypes_list.append(Ctypes_of_this_Cgroup_list[1])  # save the encountered C-type
            elif len(Ctypes_of_this_Cgroup_list)==1:   # in all other cases save the assigned C-types
                assigned_Ctypes_list.append(Ctypes_of_this_Cgroup_list[0])  # save the encountered C-type
        # check that each C-type belongs to only one C-group
        Ctype_set = set([a[1] for a in x[:-1]])
        for Ctype in Ctype_set:
            Cgroups_of_this_Ctype_list = list(set([a[5] for a in x[:-1] if a[1]==Ctype]))    # C-groups of this Cluster
            if len(Cgroups_of_this_Ctype_list) > 1 and not ( (set(Cgroups_of_this_Ctype_list)==set(["CG1", "CG2"]) and aa_type!="ILE") or
                (set(Cgroups_of_this_Ctype_list)==set(["CD1", "CD2"]) and aa_type!="TRP") or
                set(Cgroups_of_this_Ctype_list)==set(["CE1", "CE2"]) ):
                DISCARD = True
                break
        # If this carbon is partially assigned and one assignment have probability 10000000.0 (TOCSY-HCNH matched peak), then make sure that
        # the second peak assigned to this C-type belongs to the same C-group. E.g.
        # x = chainScore_list= [(['GLN', 'CG', ('HG2',), 11.562660130867117, ('1.832_32.524',), 2, 1, ('HCNH',)], 11.562660130867117), 
        # (['GLN', 'CG', ('HG2',), 10000000.0, ('2.289_33.986',), 1, 0, ('TOCSY-HCNH',)], 1.0)]
        # for Ctype in Ctype_set:
        #     Cgroups_of_this_Ctype = [a[5] for a in x[:-1] if a[1]==Ctype and a[3]==10000000.0]    # C-group of this C-type
        #     if len(Cgroups_of_this_Ctype) == 0:
        #         continue
        #     elif len(Cgroups_of_this_Ctype) == 1:
        #         Cgroup_of_this_Ctype = Cgroups_of_this_Ctype[0]
        #     else:
        #         print "ERROR: found possible chain with multiple C-groups with probability 10000000.0 for the same C-type!!!"
        #         print "x=", x
        #         sys.exit(1)
        #     peaks_of_this_Ctype = [a for a in x[:-1] if a[1]==Ctype]
        #     if len(set([p for p in peaks_of_this_Ctype if p[5]==Cgroup_of_this_Ctype])) > 1:    # if multiple C-groups are found, discard it
        #         DISCARD = True
        #         break
        # # check for H-type duplication
        # Htype_set = set([a[2] for a in x[:-1]])
        # Htype_list = [a[2] for a in x[:-1]]
        # if len(Htype_set) !=  len(Htype_list):
        #     DISCARD = True
        
        if DISCARD:
            continue
        else:
            correct_chainScore_list.append(x)
    
    correct_chainScore_list.sort(key=itemgetter(-1), reverse=True) # sort probability by descending order; the highest probability
                                                                   # combination must be 1st
    print("DEBUG: all possible predictions correct_chainScore_list=", correct_chainScore_list)
    best_prediction_list = get_best_prediction(correct_chainScore_list, iteration)
    
    print("\nDEBUG: Best Prediction:")
    for x in best_prediction_list:
        print(x)
    # print "\nDEBUG: Best Prediction after filtering with percentiles:"
    # filtered_best_prediction_list = []
    # for pred in best_prediction_list:
    #     print "DEBUG: pred=",pred,pred[3], ">=", percentile_mdict[aa_type][CH_pair][args.PERCENTILE], "aa_type=", aa_type, "CH_pair=", CH_pair, "args.PERCENTILE=", args.PERCENTILE
    #     if pred[3] >= percentile_mdict[aa_type][CH_pair][args.PERCENTILE]:   # save this C-H type prediction, only if the 2D-hist probability is above the value of the designated percentile
    #         filtered_best_prediction_list.append(pred)
    # return list(filtered_best_prediction_list)
    
    if iteration != None:   # if this is NOESY assignment
        return list(best_prediction_list)
    else:   # if TOCSY assignment, convert the list to the appropriate format
        # from ['ARG', 'CG', ('HG2', 'HG3'), 0.00241622327626462, ('1.575_27.175', '2.073_27.186'), 3, 2, ('TOCSY', 'TOCSY')]
        # to ['ARG', 'CG', 'HG2', 0.00241622327626462, '1.575', '27.175', 3, 2],
        #    ['ARG', 'CG', 'HG3', 0.00241622327626462, '2.073', '27.186', 4, 2]
        best_prediction_list_TOCSY = []
        for p in list(best_prediction_list):
            aa_type = p[0]
            Cname = p[1]
            Cgroup_score = p[3]
            Cgroup = p[6]
            for Hname, CHresonpair in zip(p[2], p[4]):
                Hreson = float(CHresonpair.split("_")[0])
                Creson = float(CHresonpair.split("_")[1])
                peak_ID = find_peakID(Hreson, Creson, aatypeResonpairMatchesTuple_list)
                print("DEBUG: [ aa_type, Cname, Hname, Cgroup_score, peak_ID, Hreson, Creson, Cgroup ] =", [ aa_type, Cname, Hname, Cgroup_score, peak_ID, Hreson, Creson, Cgroup ])
                best_prediction_list_TOCSY.append( [ aa_type, Cname, Hname, Cgroup_score, peak_ID, Hreson, Creson, Cgroup ] )
        print("\nDEBUG: Best Prediction in TOCSY format:")
        for x in best_prediction_list_TOCSY:
            print(x)
        return best_prediction_list_TOCSY


def get_presenceprob_from_all_H_C_resonpairs(possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list, Num_of_TOCSY_resonances, args):
    """
        FUNCTION to return the matched amino acid types for a particular TOCSY AAIG by using all its H,C resonance pairs.
        This function is used only for TOCSY assignment.
        
        ARGUMENTS:
        possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list:   list of lists of the form (aa type, Carbon name, Hydrogen name, weighted average probability, TOCSY_reson_index, H resonance, C resonance, carbon group)
        
        RETURN:
        aatype_CHnucleiType_presenceProbSum_mdict:  multidimensional dict: aa type --> C,H nuclei pair type --> probability of presence or absence (eq 6 in RESCUE2 paper)
        aatype_CHnucleiType_TresonIndex_mdict:      multidimensional dict: aa type --> C,H nuclei pair type --> TOCSY_reson_index
        
    """
    
    from .global_vars import Prob_CS
    Probability_CS = Prob_CS
    
    # The following sorting ensures that the alternative C-H pairs are sorted according to their weighted average probability 
    possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list.sort(key=itemgetter(0,1,2,3), reverse=True)  # sort by every column except the last 4
    #print "DEBUG: possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list=",possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list
    aatype_occupancy_dict = {}  # dict with key the aa types and values how many times they were found during the H,C resonance analysis
    aatype_probsum_dict = {}    # dictionary with the aa type --> sum of probabilities of individual matched H,C resonance pairs
    aatype_CHnucleiType_presenceProbSum_mdict = tree()    # multidimensional dict: aa type --> C,H nuclei pair type --> probability of presence or absence (eq 6 in RESCUE2 paper)
    aatype_CHnucleiType_TresonIndex_mdict = tree()    # multidimensional dict: aa type --> C,H nuclei pair type --> TOCSY_reson_index
    aa_type_set = set()
    aatype_validC_H_pairsList_dict = {}
    print("DEBUG get_presenceprob_from_all_H_C_resonpairs: possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list=",possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list)
    #sys.exit(1)
    ## The following block keeps only the matches for a particular aatype
    aa_type = None
    previous_aatype = None
    aatypeResonpairMatchesTuple_list = []
    selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list = []
    for group8 in possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list:
        aa_type = group8[0]
        weighted_average_probability = group8[3]    # = [wH * hist(H) + wC * hist(C)]/2
        if aa_type != previous_aatype and previous_aatype != None:  # this condition must be never satisfied in this script because we already know the aa_type (there is only
            # one aa_type in possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list)
            print("DEBUG: aa_type=",aa_type,"previous_aatype=",previous_aatype,"aatypeResonpairMatchesTuple_list=",aatypeResonpairMatchesTuple_list)
            if args.PROB_PRODUCT == True:
                # selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list.extend( select_correct_H_C_resonpairs2(aatypeResonpairMatchesTuple_list, args, aa_carbonBondedHydrogensDict_dict) )
                selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list.extend( select_correct_H_C_resonpairs4(aatypeResonpairMatchesTuple_list, args) )
            else:
                selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list.extend( select_correct_H_C_resonpairs(aatypeResonpairMatchesTuple_list) )
            aatypeResonpairMatchesTuple_list = []
        
        aatypeResonpairMatchesTuple_list.append(group8)
        previous_aatype = aa_type
    print("DEBUG: aa_type=",aa_type,"previous_aatype=",previous_aatype,"aatypeResonpairMatchesTuple_list=",aatypeResonpairMatchesTuple_list)
    if args.PROB_PRODUCT == True:
        # selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list.extend( select_correct_H_C_resonpairs2(aatypeResonpairMatchesTuple_list, args, aa_carbonBondedHydrogensDict_dict) )  # do the last aa_type
        selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list.extend( select_correct_H_C_resonpairs4(aatypeResonpairMatchesTuple_list, args) )
    else:
        selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list.extend( select_correct_H_C_resonpairs(aatypeResonpairMatchesTuple_list) )  # do the last aa_type
    
    if len(selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list) == 0:  # if no atom assignment could be made, skip this aa! (TEMPORARY FIX)
        return {}, {}
    #print "DEBUG: selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list=", selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list
    for group8 in selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list:
        print("DEBUG: group8=",group8)
        aa_type = group8[0]
        C_name = group8[1]
        H_name = group8[2]
        weighted_average_probability = group8[3]
        TOCSY_reson_index = group8[4]
        #if (aa_type, C_name, H_name) == previous_aatype_C_H_tuple:  # elegant way to keep only the (aatype, C_name, H_name) with the highest weighted_average_probability, e.g. we may have:
        #    # ('CYS', 'CA', 'HA', 0.012993371769822381), ('CYS', 'CA', 'HA', 0.002004997884084288), ('CYS', 'CA', 'HA', 0.0018243946693470656) but we keep only the 1st 
        #    continue
        try:
            if (C_name+"_"+H_name in aatype_validC_H_pairsList_dict[aa_type]):  # if we have saved twice this aatype,C,H combination print ERROR & exit
                print("ERROR: ",C_name+"_"+H_name,"already found !!!")
                sys.exit(1)
            aatype_validC_H_pairsList_dict[aa_type].append(C_name+"_"+H_name)
            aa_type_set.add(aa_type)
            aatype_occupancy_dict[aa_type] += 1
            ## THERE 2 DIFFENT WAYS TO CALCULATE THE OVERALL WEIGHTED AVERAGE PROBABILITY
            #aatype_probsum_dict[aa_type] += weighted_average_probability    # 1st WAY
            aatype_probsum_dict[aa_type] *= weighted_average_probability    # 2nd WAY
            aatype_CHnucleiType_presenceProbSum_mdict[aa_type][C_name+"_"+H_name] = (args.H_weight * Probability_CS[aa_type][H_name] + args.C_weight * Probability_CS[aa_type][C_name])/float((args.H_weight + args.C_weight))  # 2nd WAY
            aatype_CHnucleiType_TresonIndex_mdict[aa_type][C_name+"_"+H_name] = TOCSY_reson_index
        except KeyError:
            aatype_validC_H_pairsList_dict[aa_type] = [C_name+"_"+H_name]
            #if args.ALLOW_SINGLE_CH_PAIRS:
            aa_type_set.add(aa_type)
            aatype_occupancy_dict[aa_type] = 1
            aatype_probsum_dict[aa_type] = weighted_average_probability     # COMMON IN 1st & 2nd WAY
            aatype_CHnucleiType_presenceProbSum_mdict[aa_type][C_name+"_"+H_name] = (args.H_weight * Probability_CS[aa_type][H_name] + args.C_weight * Probability_CS[aa_type][C_name])/float((args.H_weight + args.C_weight))  # 2nd WAY
            aatype_CHnucleiType_TresonIndex_mdict[aa_type][C_name+"_"+H_name] = TOCSY_reson_index
    
    #for aa_type in aatype_CHnucleiType_presenceProbSum_mdict.keys():
    #    for k,v in aatype_CHnucleiType_presenceProbSum_mdict[aa_type].items():
    #        print "DEBUG: aa_type=", aa_type, k, v
    #for aa_type in aatype_CHnucleiType_TresonIndex_mdict.keys():
    #    for k,v in aatype_CHnucleiType_TresonIndex_mdict[aa_type].items():
    #        print "DEBUG: aa_type=", aa_type, k, v
    
    return aatype_CHnucleiType_presenceProbSum_mdict, aatype_CHnucleiType_TresonIndex_mdict
