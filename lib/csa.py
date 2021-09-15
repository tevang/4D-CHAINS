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


import math


def eliminate_orphan_Ctype_predictions(aatypeResonpairMatchesTuple_list):
    """
        FUNCTION to eliminate C-type predictions that do not contain both peaks of the same C-group. E.g.
        aatypeResonpairMatchesTuple_list= [['ARG', 'CB', 'HB2', 10000000.0, 1, 2.093, 30.83, 1, 'TOCSY-HCNH'],
        ['ARG', 'CG', 'HG2', 0.0014075550018389298, 2, 1.633, 30.846, 1, 'HCNH'],
        ['ARG', 'CG', 'HG3', 0.0014075550018389298, 2, 1.633, 30.846, 1, 'HCNH'],
        ['ARG', 'CB', 'HB3', 0.06925795331510819, 2, 1.633, 30.846, 1, 'HCNH'],
        ['ARG', 'CB', 'HB2', 0.06925795331510819, 2, 1.633, 30.846, 1, 'HCNH']]
        In this case we have one C-group which contains the peaks (1.633, 30.846) and (2.093, 30.83). However, only peak (1.633, 30.846) can be CG
        according to the 2D-histograms, therefore CG must be eliminated from all possible CS assignments of this C-group:
        filtered_aatypeResonpairMatchesTuple_list= [['ARG', 'CB', 'HB2', 10000000.0, 1, 2.093, 30.83, 1, 'TOCSY-HCNH'],
        ['ARG', 'CB', 'HB3', 0.06925795331510819, 2, 1.633, 30.846, 1, 'HCNH'],
        ['ARG', 'CB', 'HB2', 0.06925795331510819, 2, 1.633, 30.846, 1, 'HCNH']]

    """
    print("DEBUG eliminate_orphan_Ctype_predictions: aatypeResonpairMatchesTuple_list=", aatypeResonpairMatchesTuple_list)
    Cname_set = set([x[1] for x in aatypeResonpairMatchesTuple_list if x[3]>0]) # exclude Carbon-only based assignments
    Cgroups_set = set([x[7] for x in aatypeResonpairMatchesTuple_list if x[3]>0]) # exclude Carbon-only based assignments
    Cgroup_peakNum_dict = {}
    filtered_aatypeResonpairMatchesTuple_list = []
    for Cgroup in Cgroups_set:
        Cgroup_peakNum = len(set([(x[5], x[6]) for x in aatypeResonpairMatchesTuple_list if x[3]>0 and x[7]==Cgroup]))  # number of peaks of this Cgroup
        for Cname in Cname_set: # keep only the Cname predictions that are found in both peaks of this Cgroup
            if len(set([(x[5], x[6]) for x in aatypeResonpairMatchesTuple_list if x[3]>0 and x[7]==Cgroup and x[1] == Cname])) == Cgroup_peakNum:
                filtered_aatypeResonpairMatchesTuple_list.extend([x for x in aatypeResonpairMatchesTuple_list if x[3]>0 and x[7]==Cgroup and x[1]==Cname])

    print("DEBUG eliminate_orphan_Ctype_predictions: filtered_aatypeResonpairMatchesTuple_list=", filtered_aatypeResonpairMatchesTuple_list)
    return filtered_aatypeResonpairMatchesTuple_list


def get_Cgroup_consensus_probability(prob_list, cons_mode=1):
    """
    ARGUMENTS:
    prob_list:  list of the peak probabilities of this C-group
    cons_mode:       controls how to calculate the consensus probability of this C-group.
                1: average;
                2: sqrt(prob1*prob2)    ; geometric mean
                3: (prob2-prob1)/(log(prob2)-log(prob1))    ; logarithmic mean
                4: (prob1*prob2)/((prob1+prob2)/2.) ; composite average
                5: prob1*prob2    ; product of probabilities
    """
    consensus_probability = None
    prob_list = list(set(prob_list))    # remove replicate values

    if len(prob_list) == 1:
        prob1 = prob_list[0]
        # print "DEBUG get_Cgroup_consensus_probability: only one probability is available prob1=", prob1, "cons_mode=", cons_mode
        consensus_probability = prob1
    elif len(prob_list) == 2:
        prob1 = prob_list[0]
        prob2 = prob_list[1]
        # print "DEBUG get_Cgroup_consensus_probability: prob1=", prob1, "prob2=", prob2, "cons_mode=", cons_mode
        if cons_mode == 1:   # arithmetic average
            consensus_probability = (prob1 + prob2)/2.
        elif cons_mode == 2: # geometric average
            # use abs() in case when only Carbon predictions were activated and hence the probability is negative
            consensus_probability = math.sqrt(abs(prob1*prob2))
        elif cons_mode == 3: # logarithmic average
            if prob1 == 0 and prob2 == 0:
                consensus_probability = 0
            elif prob1 == prob2:
                consensus_probability = prob1
            else:
                # use abs() in case when only Carbon predictions were activated and hence the probability is negative
                consensus_probability = (prob2-prob1)/(math.log10(abs(prob2))-math.log10(abs(prob1)))
        elif cons_mode == 4: # composite average
            consensus_probability = (prob1*prob2)/((prob1+prob2)/2.)
        elif cons_mode == 5: # product of probabilities
            consensus_probability = prob1*prob2
    elif len(prob_list) > 2:
        print("ERROR: this Carbon group contains more than two peaks!!!!")
        print("DEBUG: prob_list=", prob_list)
        return None
        # sys.exit(1)

    # print "DEBUG: consensus_probability=", consensus_probability
    return consensus_probability