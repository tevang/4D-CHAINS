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


import gc
import re
import sys
from collections import defaultdict, OrderedDict
from itertools import combinations, permutations
from operator import itemgetter

import numpy as np

# def get_resid_from_residue(residue):
#
#     global absolute_matches_alignment_list, args
#
#     resid = absolute_matches_alignment_list.index(residue) + args.FIRST_RESIDUE_NUMBER
#     return resid
from cluster import HierarchicalClustering
from ete3 import Tree
from lib.csa import eliminate_orphan_Ctype_predictions, get_Cgroup_consensus_probability
from lib.global_func import ColorPrint, bcolors, approx_equal_proportional, tree, remove_NH_suffix, \
    get_residue_from_AAIGsignature, approx_equal, Debuginfo, get_resid_from_residue
from lib.global_vars import aatype_carbon_degenerateH_mdict, aa3to1_dict, aa1to3_dict, \
    aatype_carbon_nongeminalHname_mdict, aatype_carbon_methylHydrogens_mdict, aa_CarbonListwithDegenerateH_dict, \
    aa_carbonBondedGeminalHydrogensDict_dict, aatype_carbon_nondegenerateHlist_mdict, aa_carbonBondedHydrogensDict_dict
from lib.probhist import ProbHist_Calculator


def get_iminus1_residue_from_i_residue(i_residue, protein_alignment_list, absolute_matches_alignment_list, args):
    """
        FUNCTION to get the i-1 residue from the i residue using the sequence of the protein.
        ARGS:
        i_residue:       e.g. I313
        RETURNS:
        iminus1_residue: e.g. T312
    """

    i_index = absolute_matches_alignment_list.index(i_residue)
    iminus1 = (i_index-1) + int(args.FIRST_RESIDUE_NUMBER) +1 # +1 for the 'N/A' at the 1st position that was removed
    iminus1_residue = protein_alignment_list[i_index-1] + str(iminus1)
    return iminus1_residue


def get_iplus1_residue_from_i_residue(i_residue, protein_alignment_list, absolute_matches_alignment_list, args):
    """
        FUNCTION to get the i+1 residue from the i residue using the sequence of the protein.
        ARGS:
        i_residue:       e.g. I313
        RETURNS:
        iplus1_residue: e.g. T314
    """

    i_index = absolute_matches_alignment_list.index(i_residue)
    iplus1 = (i_index+1) + int(args.FIRST_RESIDUE_NUMBER) +1   # +1 for the 'N/A' at the 1st position that was removed
    try:
        iplus1_residue = protein_alignment_list[i_index+1] + str(iplus1)
    except IndexError:  # if it raises IndexError here but it didn't 2 lines above, then this is the C-term residue in NHmap table ('N/A')
        if protein_alignment_list[i_index] != 'N/A':
            return 'N/A'
        else:
            print("ERROR at get_iplus1_residue_from_i_residue().")
            print(i_residue, protein_alignment_list, absolute_matches_alignment_list, args.FIRST_RESIDUE_NUMBER)
            sys.exit(1)
    return iplus1_residue


def get_i_residue_from_iplus1_residue(iplus1_residue, protein_alignment_list, absolute_matches_alignment_list, args):
    """
        FUNCTION to get the i-1 residue from the i residue using the sequence of the protein.
        ARGS:
        i_residue:       e.g. I313
        RETURNS:
        iplus1_residue: e.g. T314
    """

    iplus1_index = absolute_matches_alignment_list.index(iplus1_residue)
    i_index = iplus1_index -1
    i = (iplus1_index-1) + int(args.FIRST_RESIDUE_NUMBER) +1   # +1 for the 'N/A' at the 1st position that was removed
    i_residue = protein_alignment_list[i_index] + str(i)    # get the aa type from the alignment and prepend it to the resid
    return i_residue


def get_i_residue_from_iminus1_residue(iminus1_residue, protein_alignment_list, absolute_matches_alignment_list, args):
    """
        FUNCTION to get the i residue from the i-1 residue using the sequence of the protein.
        ARGS:
        i_residue:       e.g. T312
        RETURNS:
        iminus1_residue: e.g. I313
    """

    iminus1_index = absolute_matches_alignment_list.index(iminus1_residue)
    i = (iminus1_index+1) + int(args.FIRST_RESIDUE_NUMBER) +1 # +1 for the 'N/A' at the 1st position that was removed
    i_residue = protein_alignment_list[iminus1_index+1] + str(i)    # get the aa type from the alignment and prepend it to the resid
    return i_residue


def get_aa_type_from_residue(residue):
    """
    RETURN:
    aa_type:    in 1-letter code
    """
    aa_type = residue[0]
    return aa_type


# def get_aa_type_from_residue(residue):
#     """
#     RETURN:
#     aa_type:    in 1-letter code
#     """
#
#     global protein_alignment_list, absolute_matches_alignment_list, args
#
#     index = absolute_matches_alignment_list.index(residue)
#     aa_type = protein_alignment_list[index]
#     return aa_type


def get_aa_type_from_resid(resid, protein_alignment_list, first_resid=1):
    """
    RETURN:
    aa_type:    in 1-letter code
    """
    print("DEBUG get_aa_type_from_resid: resid=", resid, "first_resid=", first_resid)
    aa_type = aa1to3_dict[protein_alignment_list[resid -1 - first_resid]]   # -1 for the 'N/A' at the 1st position that was removed
    # aa_type = aa1to3_dict[protein_alignment_list[resid - first_resid]]   # use this if the 'N/A' at the 1st position was retained (I guess onlyNOESY)
    return aa_type


def get_iminus1_aa_type_from_i_residue(i_residue, protein_alignment_list, absolute_matches_alignment_list, i_to_iminus1_dict):

    iminus1_residue = i_to_iminus1_dict[i_residue]
    i_index = absolute_matches_alignment_list.index(i_residue)
    iminus1_index = i_index -1
    aa_type = protein_alignment_list[iminus1_index]

    return aa_type


def is_valid_residue(residue, protein_alignment_list, args, spectrum_combo):
    """
        FUNCTION to check if a residue has the format [A-Z][0-9]+ and if the number is a valid resid.
    """
    mo = re.search("([ACDEFGHIKLMNPQRSTVWY])([0-9]+)", residue)
    if not mo:
        return False
    else:
        aa_type = mo.group(1)
        resid = int(mo.group(2))

    if spectrum_combo == 'TOCSY-HCNH':
        index = resid - args.FIRST_RESIDUE_NUMBER -1 # -1 because we removed the first residue that was mapped to 'N/A' (works for cs_assignment_HCNH.py)
    elif spectrum_combo == 'HCNH-HCNH':
        index = resid - args.FIRST_RESIDUE_NUMBER   # works for cs_assignment_onlyHCNH only.
    try:
        if protein_alignment_list[index] == aa_type:    # if both the resid and the aa_type are correct
            return True
        else:
            return False
    except IndexError:  # for the case that the resid is not valid
        return False


def is_in_alignment(RIG, absolute_matches_alignment_list):
    """
        onlyHCNH function.
        FUNCTION to check if a root index group (RIG) has been placed in the alignment.
    """
    if RIG in absolute_matches_alignment_list:
        return True
    else:
        return False


def rename_lone_methylene_proton(p, aatype_carbon_degenerateH_mdict, aatype_carbon_nondegenerateHlist_mdict):
    """
        RENAME THE METHYLENE PROTONS FROM "QX"->"HG[23]".
        ARGS:
        p:  a peak in the form ['S393', 'CA', 56.639, 'HA', 4.593]
    """

    # aatype_carbon_nondegenerateHlist_mdict['LYS']['CG'] = ['HG2', 'HG3']
    if p[3][0]=="Q":
        aatype = aa1to3_dict[p[0][0]]
        if p[3] in aatype_carbon_degenerateH_mdict[aatype][p[1]]:
            p[3] = aatype_carbon_nondegenerateHlist_mdict[aatype][p[1]][0]  # rename it to the first methylene proton, e.g. LYS QG->HG2

    return tuple(p)


def read_NHassigned_spectrum_file(query_fname, spectrum_combo, ali, i_to_iminus1_dict=None):
    """
    Method to read the TOCSY or HCNH file with the assignments in SPARKY format, as produced by
    cs_assignment.py script.

    :param query_fname:
    :param spectrum_combo: "TOCSY" or "HCNH"
    :param ali:
    :param i_to_iminus1_dict:
    :return:  i_residue & iplus1_residue in this function are RIG not valid residue names!
    :return:  residue_assignments_dict:    i_residue ->   [ (iminus1_residue, Cname, Cresonance, Hname, Hresonance),
                                                            (iminus1_residue, Cname, Cresonance, Hname, Hresonance), ... ]
    :return:  residue_NHresonances_dict:   i_residue -> (average_Nresonance, average_Hresonance, stdev_N_resonance, stdev_Hreson)
    """

    print("Loading "+spectrum_combo+" file "+query_fname)
    def normalize_intensities():
        """
            FUNCTION to normalize the intensities; each residue must have peaks with max intensity = 1.0.
        """
        for residue in list(residue_peak_intensity_mdict.keys()):
            max_intensity = max([ residue_peak_intensity_mdict[residue][peak] for peak in residue_peak_intensity_mdict[residue] ])
            for peak in residue_peak_intensity_mdict[residue]:
                # print "DEBUG: residue", residue, "peak", peak, "intensity normalization:", residue_peak_intensity_mdict[residue][peak], "/", max_intensity
                residue_peak_intensity_mdict[residue][peak] /= max_intensity    # normalize by dividing by the max intensity value for this residue
                # print "DEBUG: intensity after normalization:", residue_peak_intensity_mdict[residue][peak]

    with open(query_fname, 'r') as f:
        tmp_query_contents=f.readlines()    # contents of original query_fname (4D TOCSY or 4D HCNH) in 5 column format (name H C N HN)
    query_contents = []
    query_CS_set = set() # store the numbers here to discard replicate lines
    if spectrum_combo == "TOCSY":
        for line in tmp_query_contents: # CLEANING THE SPECTRUM FROM IRRELEVANT LINES
            word_list = line.split()
            try:
                float(word_list[1]); float(word_list[2]); float(word_list[3]); float(word_list[4]); # checking if the line contains numbers (resonances)
                if len(word_list)==6:   # if there is an extra column chech if it is intensity number
                    float(word_list[5])
                    CS_values = (float(word_list[1]), float(word_list[2]), float(word_list[3]),
                                 float(word_list[4]), float(word_list[5]))
                else:
                    CS_values = (float(word_list[1]), float(word_list[2]), float(word_list[3]),
                                float(word_list[4]))
                # Keep only only assigned lines
                if not '?' in word_list[0][0:7] and not CS_values in query_CS_set:
                    #print "DEBUG: appending line:", line
                    query_contents.append(line)
                    query_CS_set.add(CS_values)
                else:
                    ColorPrint("WARNING: Discarding TOCSY line: " + line, "WARNING")
            except (IndexError, ValueError):
                print(bcolors.WARNING + "WARNING: Discarding TOCSY line: " + line + bcolors.ENDC)
        # Read the patched N-terminal residues written as a comment assigned TOCSY file (in sparky format)
        patched_residues_list = []
        for line in reversed(tmp_query_contents):   # ATTENTION: read the original file contents, not query_contents!!!
            if line[:19] == '# PATCHED RESIDUES:':
                patched_residues_list = line[19:].split()
                break
    elif spectrum_combo == "HCNH":
        for line in tmp_query_contents: # CLEANING THE SPECTRUM FROM IRRELEVANT LINES
            word_list = line.split()
            try:
                float(word_list[1]); float(word_list[2]); float(word_list[3]); float(word_list[4]); # checking if the line contains numbers (resonances)
                if len(word_list)==6:   # if there is an extra column check if it is intensity number
                    float(word_list[5])
                    CS_values = (float(word_list[1]), float(word_list[2]), float(word_list[3]),
                                 float(word_list[4]), float(word_list[5]))
                else:
                    CS_values = (float(word_list[1]), float(word_list[2]), float(word_list[3]),
                                float(word_list[4]))
                # if word_list[0][0:7] != "?-?-?-?":
                if not CS_values in query_CS_set:
                    print("DEBUG: appending HCNH line:", line)
                    query_contents.append(line)
                    query_CS_set.add(CS_values)
                else:
                    ColorPrint("WARNING: Discarding replicate HCNH line: " + line, "WARNING")
            except (IndexError, ValueError):
                ColorPrint("WARNING: Discarding invalid HCNH line:" + line, "WARNING")

    # remove duplicate lines from query_fname (4D TOCSY or 4D HCNH)
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
                        #print "DEBUG: will remove line",qline2
        except (ValueError, IndexError):
            #print "WARNING: the 2nd and 3rd elements of the following root file line are not numbers:"
            #print "Root file line:", root_line
            continue
    # now remove the duplicate lines
    for qline in lines2remove_set:
        query_contents.remove(qline)

    # populate a dictionary with keys the resids and values the list of the respective C-H & N-H resonances from the spectrum file (TOCSY or HCNH)
    query_lineLists_list = []   # list of the lines of query frame in list form not in string
    for qline in query_contents:
        qline_list = qline.split()
        query_lineLists_list.append(qline_list)
    sorted_query_lineLists_list = sorted(query_lineLists_list, key=itemgetter(0))   # sort the spectrum lines by the assignment
    residue_assignments_dict = {}    # i_residue -> [ (i_residue, Cname, Cresonance, Hname, Hresonance), (i_residue, Cname, Cresonance, Hname, Hresonance), ... ]
    residue_NHresonances_dict = {}   # iplus1_residue -> [(iplus1_residue, Nresonance, Hresonance), (iplus1_residue, Nresonance, Hresonance), ...]
    original_HCNH_peaks_dict = {}   # residue -> [Hreson, Creson, Nreson, HNreson]
    residue_peak_intensity_mdict = tree()    # peak (e.g. (23.45, 0,823)) -> intensity
    if spectrum_combo == "TOCSY":
        i_to_iminus1_dict = {}  # dictionary to be populated from TOCSY assignment file
    for qline_list in sorted_query_lineLists_list:
        print("DEBUG: ------------------------->")
        print("DEBUG: reading line:", qline_list)
        if qline_list[0] == '?-?-?-?': # if this is an non-labeled peak (?-?-?-?), skipt it
            print("WARNING: the following line will be omitted because it is not assigned:", qline_list)
            continue
        elif not (qline_list[0].endswith('N-H') or qline_list[0].endswith('NX-HX')):
            print("WARNING: the following line will be omitted because the amide is not of the backbone's:", qline_list)
            continue
        components = remove_NH_suffix(qline_list[0]).split('-')    # this will give something like ['A297HA', 'CA', 'A298']
        #print "DEBUG: iplus1_residue=", iplus1_residue
        i_Cname = components[1]
        if spectrum_combo == "TOCSY":
            mo = re.search('([A-Za-z0-9][0-9]+)([HQM][ABGDEHZ][1-3]{0,2})', components[0])    # read only lines with atom type assignments
            if mo:
                i_residue = mo.group(1)
                i_Hname = mo.group(2)
            else:
                print("ERROR: cannot find Hname and TAAIG in string:", components[0], ". Please correct the following line in TOCSY:", qline_list)
                sys.exit(1)
        elif spectrum_combo == 'HCNH':
            if components[0] == '?' and components[1] == '?':
                i_Hname = "?"
                i_residue = "?"
            else:   # only if this is a user-provided annotated HCNH file, it will contain the i residue label (apart from the i-1)
                mo = re.search('([A-Za-z0-9][0-9]+)([HQM][ABGDEHZ][1-3]{0,2})', components[0])     # read only lines with atom type assignments
                if mo:
                    i_residue = mo.group(1)
                    i_Hname = mo.group(2)
                else:
                    print("WARNING: the following line will be omitted because Hname and TAAIG could not be found:", components[0])
                    # sys.exit(1)
                    continue

        # Rename VAL QG1-->HG1, QG2-->HG2 and LEU QD1-->HD1, QD2-->HD2 for compatibility with the code and the histogram file names
        if i_residue[0] == 'L' and i_Hname == "QD1":
            i_Hname = "HD1"
        elif i_residue[0] == 'L' and i_Hname == "QD2":
            i_Hname = "HD2"
        if i_residue[0] == 'V' and i_Hname == "QG1":
            i_Hname = "HG1"
        elif i_residue[0] == 'V' and i_Hname == "QG2":
            i_Hname = "HG2"

        if components[2] == '':
            iplus1_residue = i_residue
        else:   # if this is a HCNH *num.list file or user-provided annotated HCNH file,
                # it will contain the i residue label (apart from the i-1)
            iplus1_residue = get_residue_from_AAIGsignature(components[2],
                                                            ali.absolute_AAIGmatches_alignment_list,
                                                            ali.absolute_matches_alignment_list)  # rename the AAIG to the residue actual name for consistency with TOCSY

        print("DEBUG: i_residue=", i_residue, "iplus1_residue=", iplus1_residue,
              "spectrum_combo=", spectrum_combo, "components=", components)
        #print "DEBUG: point 0, iplus1_residue=", iplus1_residue
        if spectrum_combo == "TOCSY" and i_residue not in list(residue_assignments_dict.keys()) and not iplus1_residue in list(residue_NHresonances_dict.keys()):
            print("DEBUG read_assigned_TOCSY_file: point 1")
            residue_assignments_dict[i_residue] = []
            residue_NHresonances_dict[iplus1_residue] = []
        elif spectrum_combo == 'HCNH' and \
                iplus1_residue not in list(residue_assignments_dict.keys()) and\
                iplus1_residue not in list(residue_NHresonances_dict.keys()):
            print("DEBUG: point 1, iplus1_residue=", iplus1_residue)
            # if is_in_alignment(iplus1_residue, absolute_matches_alignment_list) == True:   # if this is in the alignment
            residue_assignments_dict[iplus1_residue] = []
            original_HCNH_peaks_dict[iplus1_residue] = []
            NresonIndex = 0
            residue_NHresonances_dict[iplus1_residue] = []
        if spectrum_combo == "TOCSY":
            print("DEBUG read_assigned_TOCSY_file: point 2")
            try:
                residue_assignments_dict[i_residue].append( rename_lone_methylene_proton([i_residue, i_Cname, float(qline_list[2]), i_Hname, float(qline_list[1])], aatype_carbon_degenerateH_mdict, aatype_carbon_nondegenerateHlist_mdict) )
            except KeyError:
                residue_assignments_dict[i_residue] = [ rename_lone_methylene_proton([i_residue, i_Cname, float(qline_list[2]), i_Hname, float(qline_list[1])], aatype_carbon_degenerateH_mdict, aatype_carbon_nondegenerateHlist_mdict) ]
            try:    # if the last column is the intensity, save it
                residue_peak_intensity_mdict[i_residue][(float(qline_list[2]), float(qline_list[1]))] = float(qline_list[5])
            except IndexError:  # if the intensity column is missing, add intensity 1.0 to all peaks
                residue_peak_intensity_mdict[i_residue][(float(qline_list[2]), float(qline_list[1]))] = 1.0
                pass
            i_to_iminus1_dict[iplus1_residue] = i_residue
        elif spectrum_combo == 'HCNH':
            print("DEBUG: i_residue=", i_residue, "iplus1_residue=", iplus1_residue,"i_Hname=", i_Hname,"i_Cname=", i_Cname)
            # if is_in_alignment(iplus1_residue, absolute_matches_alignment_list) == False:   # if this is not in the alignment
            #     print "DEBUG: iplus1_residue", iplus1_residue, "is not in alignment!"
            #     print "WARNING: the following line will be omitted because",iplus1_residue,"is not in alignment:", absolute_matches_alignment_list
            #     continue
            print("DEBUG: point 2, iplus1_residue=", iplus1_residue)
            residue_assignments_dict[iplus1_residue].append([i_residue, i_Cname, float(qline_list[2]), i_Hname, float(qline_list[1])])
            try:   # if the last column is the intensity, save it
                residue_peak_intensity_mdict[iplus1_residue][(float(qline_list[2]), float(qline_list[1]))] = float(qline_list[5])
            except IndexError:  # if the intensity column is missing, add intensity 1.0 to all peaks
                residue_peak_intensity_mdict[iplus1_residue][(float(qline_list[2]), float(qline_list[1]))] = 1.0
            NresonIndex += 1
            original_HCNH_peaks_dict[iplus1_residue].append([NresonIndex, float(qline_list[1]), float(qline_list[2]), float(qline_list[3]), float(qline_list[4])])
        try:
            residue_NHresonances_dict[iplus1_residue].append((iplus1_residue, float(qline_list[3]), float(qline_list[4])))
        except KeyError:
            residue_NHresonances_dict[iplus1_residue] = [ (iplus1_residue, float(qline_list[3]), float(qline_list[4])) ]

    # finaly average the N and HN resonances for each residue
    # after that residue_NHresonances_dict: residue -> (average_Nresonance, average_Hresonance)
    # print "DEBUG: point X residue_NHresonances_dict=", residue_NHresonances_dict
    for residue in list(residue_NHresonances_dict.keys()):
        N_resonance_list = []
        Hreson_list = []
        for triplet in residue_NHresonances_dict[residue]:
            N_resonance_list.append(triplet[1])
            Hreson_list.append(triplet[2])
        average_N_resonance = round(np.average(N_resonance_list), 3)
        average_Hreson = round(np.average(Hreson_list), 3)
        stdev_N_resonance = np.std(N_resonance_list)
        stdev_Hreson = np.std(Hreson_list)
        residue_NHresonances_dict[residue] = (average_N_resonance, average_Hreson, stdev_N_resonance, stdev_Hreson)

    # Normalize the intensities; each residue must have peaks with max intensity = 1.0
    if len(list(residue_peak_intensity_mdict.keys())) > 0:
        normalize_intensities()
        for residue in list(residue_peak_intensity_mdict.keys()):
            print("DEBUG: residue=", residue, "normalized intensities=", residue_peak_intensity_mdict[residue])

    # print "DEBUG read_assigned_TOCSY_file: residue_peak_intensity_mdict.keys()=", residue_peak_intensity_mdict.keys()
    if spectrum_combo == "TOCSY":
        # print "DEBUG read_assigned_TOCSY_file: returning", residue_assignments_dict, residue_NHresonances_dict, i_to_iminus1_dict, patched_residues_list, residue_peak_intensity_mdict
        return residue_assignments_dict, \
               residue_NHresonances_dict, \
               i_to_iminus1_dict, \
               patched_residues_list, \
               residue_peak_intensity_mdict
    elif spectrum_combo == 'HCNH':
        print("DEBUG read_assigned_HCNH_file: returning",
              residue_assignments_dict,
              residue_NHresonances_dict,
              original_HCNH_peaks_dict,
              residue_peak_intensity_mdict)
        return residue_assignments_dict, \
               residue_NHresonances_dict, \
               original_HCNH_peaks_dict, \
               residue_peak_intensity_mdict


def create_equivalentCarbons_assignments_dict(user_HCNH_resid_assignments_dict,
                                              user_TOCSY_resid_assignments_dict,
                                              protein_alignment_list,
                                              absolute_matches_alignment_list,
                                              aa_equivalentCarbons_dict,
                                              args):
    """

    :param user_HCNH_resid_assignments_dict:
    :param user_TOCSY_resid_assignments_dict:
    :param protein_alignment_list:
    :param absolute_matches_alignment_list:
    :param aa_equivalentCarbons_dict:
    :param args:
    :return user_special_HCNH_resid_assignments_dict: same as user_HCNH_resid_assignments_dict but only for VAL and LEU,
                                                      and at the end has an extra element True/False to indicate if
                                                      this resonance has been assigned
    :return user_special_TOCSY_resid_assignments_dict: same as user_TOCSY_resid_assignments_dict but only for VAL
                                                      and LEU, and at the end has an extra element True/False to
                                                      indicate if this resonance has been assigned
    """
    # Find the user assignments of all VAL and LEU and save them in a separate dict with an extra True/False indicator
    user_special_HCNH_resid_assignments_dict = {}  # same as user_HCNH_resid_assignments_dict but only for VAL and LEU, and at the end has an extra element True/False to indicate if this resonance has been assigned
    for resid in list(user_HCNH_resid_assignments_dict.keys()):
        aa_type = get_aa_type_from_resid(resid, protein_alignment_list, first_resid=args.FIRST_RESIDUE_NUMBER)
        print("DEBUG: protein_alignment_list=", protein_alignment_list)
        print("DEBUG: resid=%i aa_type=%s" % (resid, aa_type))
        residue = aa3to1_dict[aa_type] + str(resid)
        if aa_type in list(aa_equivalentCarbons_dict.keys()):
            all_equivalentCarbons_list = []
            [all_equivalentCarbons_list.extend(C_pair) for C_pair in aa_equivalentCarbons_dict[aa_type]]
            print("DEBUG HCNH: resid=%i , all_equivalentCarbons_list=%s" % (resid, all_equivalentCarbons_list))
            for equivalentCarbon in all_equivalentCarbons_list:
                user_assignment = get_Carbon_resonance(resid,
                                                       equivalentCarbon,
                                                       "HCNH",
                                                       user_HCNH_resid_assignments_dict,
                                                       user_TOCSY_resid_assignments_dict,
                                                       get_assignment=True)
                if user_assignment != False:
                    try:
                        user_special_HCNH_resid_assignments_dict[resid].append(user_assignment)
                    except KeyError:
                        user_special_HCNH_resid_assignments_dict[resid] = [user_assignment]
            # # # for assignment in user_HCNH_resid_assignments_dict[resid]:
            # # #     if assignment[1] in all_equivalentCarbons_list and assignment[0] == residue:    # sometimes it contains peaks of the i-1
            # # #         try:
            # # #             user_special_HCNH_resid_assignments_dict[resid].append([assignment[0], assignment[1], assignment[2], assignment[3], assignment[4], False])
            # # #         except:
            # # #             user_special_HCNH_resid_assignments_dict[resid] = [[assignment[0], assignment[1], assignment[2], assignment[3], assignment[4], False]]

    user_special_TOCSY_resid_assignments_dict = {}  # same as user_TOCSY_resid_assignments_dict but only for VAL and LEU, and at the end has an extra element True/False to indicate if this resonance has been assigned
    for resid in list(user_TOCSY_resid_assignments_dict.keys()):
        aa_type = get_aa_type_from_resid(resid, protein_alignment_list, first_resid=args.FIRST_RESIDUE_NUMBER)
        residue = aa3to1_dict[aa_type] + str(resid)
        if aa_type in list(aa_equivalentCarbons_dict.keys()):
            all_equivalentCarbons_list = []
            [all_equivalentCarbons_list.extend(C_pair) for C_pair in aa_equivalentCarbons_dict[aa_type]]
            print("DEBUG TOCSY: resid=%i , all_equivalentCarbons_list=%s" % (resid, all_equivalentCarbons_list))
            for equivalentCarbon in all_equivalentCarbons_list:
                user_assignment = get_Carbon_resonance(resid, equivalentCarbon, "TOCSY", user_HCNH_resid_assignments_dict, user_TOCSY_resid_assignments_dict,
                                                       get_assignment=True)
                if user_assignment != False:
                    try:
                        user_special_TOCSY_resid_assignments_dict[resid].append(user_assignment)
                    except KeyError:
                        user_special_TOCSY_resid_assignments_dict[resid] = [user_assignment]
            # # # for assignment in user_TOCSY_resid_assignments_dict[resid]:
            # # #     if assignment[1] in all_equivalentCarbons_list and assignment[0] == residue:    # sometimes it contains peaks of the i-1
            # # #         try:
            # # #             user_special_TOCSY_resid_assignments_dict[resid].append([assignment[0], assignment[1], assignment[2], assignment[3], assignment[4], False])
            # # #         except:
            # # #             user_special_TOCSY_resid_assignments_dict[resid] = [[assignment[0], assignment[1], assignment[2], assignment[3], assignment[4], False]]

    print("DEBUG: user_special_HCNH_resid_assignments_dict=", user_special_HCNH_resid_assignments_dict)
    print("DEBUG: user_special_TOCSY_resid_assignments_dict=", user_special_TOCSY_resid_assignments_dict)
    return user_special_HCNH_resid_assignments_dict, user_special_TOCSY_resid_assignments_dict


def get_Carbon_resonance(resid,
                         Cname,
                         spectrum_combo,
                         user_HCNH_resid_assignments_dict,
                         user_TOCSY_resid_assignments_dict,
                         get_assignment=False):
    """
        FUNCTION to get the user-assigned resonance for a specific Carbon and residue. If there are two peaks
        for the same Carbon (CH2), it returns the average resonance. If get_assignment=True it returns the
        whole peak assignments, e.g. ['E128', 'CA', 60.684, 'HA', 3.834]
    """
    print("DEBUG get_Carbon_resonance: resid=", resid)

    if spectrum_combo == "HCNH":
        Creson_list = []
        Hreson_list = []
        Cassignment_list = []
        # check in the HCNH peaks of i residue
        if resid in user_HCNH_resid_assignments_dict.keys():
            for assignment in user_HCNH_resid_assignments_dict[resid]:
                if assignment[0] != '?' and int(assignment[0][1:]) == resid and assignment[1] == Cname:
                    Creson_list.append(assignment[2])
                    Hreson_list.append(assignment[4])
                    Cassignment_list.append(assignment)
        # check in the HCNH peaks of i-1 residue
        if resid-1 in user_HCNH_resid_assignments_dict.keys():
            for assignment in user_HCNH_resid_assignments_dict[resid-1]:
                if assignment[0] != '?' and int(assignment[0][1:]) == resid and assignment[1] == Cname:
                    Creson_list.append(assignment[2])
                    Hreson_list.append(assignment[4])
                    Cassignment_list.append(assignment)
        # check in the HCNH peaks of i+1 residue
        if resid+1 in user_HCNH_resid_assignments_dict.keys():
            for assignment in user_HCNH_resid_assignments_dict[resid+1]:
                print("DEBUG: assignment=", assignment)
                if assignment[0] != '?' and int(assignment[0][1:]) == resid and assignment[1] == Cname:
                    Creson_list.append(assignment[2])
                    Hreson_list.append(assignment[4])
                    Cassignment_list.append(assignment)

        if len(Creson_list) == 0: # if not found return False
            return False
        else:   # otherwise return the average resonance of as many times as Cname was found in the user HCNH assignments
            if get_assignment == False:
                return round(np.average(Creson_list), 3)
            elif get_assignment == True:
                first_assignment = Cassignment_list[0]
                return [first_assignment[0], first_assignment[1], round(np.average(Creson_list), 3), first_assignment[3], round(np.average(Hreson_list), 3), False]

    elif spectrum_combo == "TOCSY":
        Creson_list = []
        Hreson_list = []
        Cassignment_list = []
        # check in the TOCSY peaks of i residue
        if resid in list(user_TOCSY_resid_assignments_dict.keys()):
            for assignment in user_TOCSY_resid_assignments_dict[resid]:
                if assignment[0] != '?' and int(assignment[0][1:]) == resid and assignment[1] == Cname:
                    Creson_list.append(assignment[2])
                    Hreson_list.append(assignment[4])
                    Cassignment_list.append(assignment)

        if len(Creson_list) == 0: # if not found return False
            return False
        else:   # otherwise return the average resonance
            if get_assignment == False:
                return round(np.average(Creson_list), 3)
            elif get_assignment == True:
                first_assignment = Cassignment_list[0]
                return [first_assignment[0], first_assignment[1], round(np.average(Creson_list), 3), first_assignment[3], round(np.average(Hreson_list), 3), False]


def are_peaks_equal(program_Cassignment, program_Hassignment, user_assignment, ignoreH=False):
    """
        FUNCTION to compare the C and H resonances of 2 peaks in the form of assignments made by the program
        and the user.
        E.g.
        program_Cassignment= [1094, 23.351, 0.2, 'CD1', 135, '#', 'LEU', 'TOCSY-HCNH']
        program_Hassignment= [1101, 0.823, 0.02, 'HD1', 135, '#', 'LEU', 'TOCSY-HCNH']
        user_assignment= ['L135', 'CD2', 23.346, 'HD2', 0.823, False]
    """

    Htol = 0.04
    Ctol = 0.4

    program_Creson = program_Cassignment[1]
    program_Hreson = program_Hassignment[1]
    user_Creson = user_assignment[2]
    user_Hreson = user_assignment[4]

    if not ignoreH:
        return approx_equal(program_Creson, user_Creson, Ctol) and approx_equal(program_Hreson, user_Hreson, Htol)
    else:
        return approx_equal(program_Creson, user_Creson, Ctol)


def compare_equivalentMethyl_peaks(program_assignments_list, user_assignments_list, ignoreH=False):
    """
    Method to compare a list of program assigned resonances for the same type of Carbon
    (may be one of the equivalent methyl Carbons of LEU or VAL) with a list of user assigned resonances.
    It updates user_assignments_list and its parent dicts
    user_VAL_LEU_HCNH_resid_assignments_dict/user_VAL_LEU_TOCSY_resid_assignments_dict.
    # IF 1 carbon is curated but I find 2 that do not match, print <WRONG> for both
    # IF 1 carbon is curated but I find 2 one of which matches, print <CORRECT> and <WRONG>
                (I don't check which is closer to the user assigned resonance)

    Input Examples:
    program_assignments_list = [[855, 23.825, 0.2, 'CD1', 97, '#', 'LEU', 'TOCSY-HCNH'],
                                [856, 26.265, 0.2, 'CD2', 97, '#', 'LEU', 'TOCSY-HCNH']
                                plus the 2 protons]
    user_assignments_list = [['L97', 'CD1', 26.265, 'QD1', 0.967, False],
                            ['L97', 'CD2', 23.825, 'QD2', 0.819, False]]

    ARGS:
    ignoreH:    proof-read only the Carbons of the methyls. The same label CORRECT/WRONG is added to the protons, too.
    """
    print("DEBUG compare_equivalentMethyl_peaks: program_assignments_list=", program_assignments_list)
    print("DEBUG compare_equivalentMethyl_peaks: user_assignments_list=", user_assignments_list)

    program_Cassignments_list = [a for a in program_assignments_list if a[3][0]=='C']
    program_Hassignments_list = []
    Carbon_pair = [a[3] for a in program_Cassignments_list]
    Proton_pair = ['H'+C[1:] for C in Carbon_pair]
    for H in Proton_pair:
        for a in program_assignments_list:
            if a[3] == H:
                program_Hassignments_list.append(a)

    print("DEBUG compare_equivalentMethyl_peaks: program_Cassignments_list=", program_Cassignments_list)
    print("DEBUG compare_equivalentMethyl_peaks: program_Hassignments_list=", program_Hassignments_list)

    if len(program_Cassignments_list) == 1 and len(user_assignments_list) == 1:
        # program_assignment1 = program_assignments_list[0]   ;# if you do this then program_assignment1 will not be a reference of program_assignments_list[0]
        # user_assignment1 = user_assignments_list[0]
        if are_peaks_equal(program_Cassignments_list[0], program_Hassignments_list[0], user_assignments_list[0], ignoreH=ignoreH) == True:
            user_assignments_list[0][5] = True     # this must update automatically user_VAL_LEU_HCNH_resid_assignments_dict/user_VAL_LEU_TOCSY_resid_assignments_dict
            if " <CORRECT>" not in program_Cassignments_list[0][7]: program_Cassignments_list[0][7] += " <CORRECT>" # + iteration_comment
            if " <CORRECT>" not in program_Hassignments_list[0][7]: program_Hassignments_list[0][7] += " <CORRECT>" # + iteration_comment
        else:
            if " <WRONG>" not in program_Cassignments_list[0][7]: program_Cassignments_list[0][7] += " <WRONG>" # + iteration_comment
            if " <WRONG>" not in program_Hassignments_list[0][7]: program_Hassignments_list[0][7] += " <WRONG>" # + iteration_comment

    elif len(program_Cassignments_list) == 1 and len(user_assignments_list) == 2:
        if are_peaks_equal(program_Cassignments_list[0], program_Hassignments_list[0], user_assignments_list[0], ignoreH=ignoreH) == True:
            user_assignments_list[0][5] = True     # this must update automatically user_VAL_LEU_HCNH_resid_assignments_dict/user_VAL_LEU_TOCSY_resid_assignments_dict
            if " <CORRECT>" not in program_Cassignments_list[0][7]: program_Cassignments_list[0][7] += " <CORRECT>" # + iteration_comment
            if " <CORRECT>" not in program_Hassignments_list[0][7]: program_Hassignments_list[0][7] += " <CORRECT>" # + iteration_comment
        elif are_peaks_equal(program_Cassignments_list[0], program_Hassignments_list[0], user_assignments_list[1], ignoreH=ignoreH) == True:
            user_assignments_list[1][5] = True     # this must update automatically user_VAL_LEU_HCNH_resid_assignments_dict/user_VAL_LEU_TOCSY_resid_assignments_dict
            if " <CORRECT>" not in program_Cassignments_list[0][7]: program_Cassignments_list[0][7] += " <CORRECT>" # + iteration_comment
            if " <CORRECT>" not in program_Hassignments_list[0][7]: program_Hassignments_list[0][7] += " <CORRECT>" # + iteration_comment
        else:
            if " <WRONG>" not in program_Cassignments_list[0][7]: program_Cassignments_list[0][7] += " <WRONG>" # + iteration_comment
            if " <WRONG>" not in program_Hassignments_list[0][7]: program_Hassignments_list[0][7] += " <WRONG>" # + iteration_comment

    elif len(program_Cassignments_list) == 2 and len(user_assignments_list) == 1:
        # program_assignment1 = program_Cassignments_list[0]
        # program_assignment2 = program_Cassignments_list[1]
        # user_assignment1 = user_assignments_list[0]
        if are_peaks_equal(program_Cassignments_list[0], program_Hassignments_list[0], user_assignments_list[0], ignoreH=ignoreH) == True:
            user_assignments_list[0][5] = True     # this must update automatically user_VAL_LEU_HCNH_resid_assignments_dict/user_VAL_LEU_TOCSY_resid_assignments_dict
            if " <CORRECT>" not in program_Cassignments_list[0][7]: program_Cassignments_list[0][7] += " <CORRECT>" # + iteration_comment
            if " <CORRECT>" not in program_Hassignments_list[0][7]: program_Hassignments_list[0][7] += " <CORRECT>" # + iteration_comment
        else:
            if " <CHECK MANUALLY>" not in program_Cassignments_list[0][7]: program_Cassignments_list[0][7] += " <CHECK MANUALLY>" # + iteration_comment
            if " <CHECK MANUALLY>" not in program_Hassignments_list[0][7]: program_Hassignments_list[0][7] += " <CHECK MANUALLY>" # + iteration_comment

        if are_peaks_equal(program_Cassignments_list[1], program_Hassignments_list[1], user_assignments_list[0], ignoreH=ignoreH) == True and user_assignments_list[0][5] == False:
            user_assignments_list[0][5] = True     # this must update automatically user_VAL_LEU_HCNH_resid_assignments_dict/user_VAL_LEU_TOCSY_resid_assignments_dict
            if " <CORRECT>" not in program_Cassignments_list[1][7]: program_Cassignments_list[1][7] += " <CORRECT>" # + iteration_comment
            if " <CORRECT>" not in program_Hassignments_list[1][7]: program_Hassignments_list[1][7] += " <CORRECT>" # + iteration_comment
        elif user_assignments_list[0][5] == True:
            if " <CHECK MANUALLY>" not in program_Cassignments_list[1][7]: program_Cassignments_list[1][7] += " <CHECK MANUALLY>" # + iteration_comment
            if " <CHECK MANUALLY>" not in program_Hassignments_list[1][7]: program_Hassignments_list[1][7] += " <CHECK MANUALLY>" # + iteration_comment

    elif len(program_Cassignments_list) == 2 and len(user_assignments_list) == 2:
        # CASES:
        # 1. p1==u1, p1==u2, p2==u1, p2==u2
        # 2. ... many cases ... uncosidered here!!!


        # If a program-assigned peak is equal (within the tolerances) with both user-assigned peaks
        # (e.g. protein MS6282 and residue V65 CG1-QG1 & CG2-QG2), then select the closest peak when
        # you do the proofreading
        if ( are_peaks_equal(program_Cassignments_list[0], program_Hassignments_list[0], user_assignments_list[0], ignoreH=ignoreH) == True and
            are_peaks_equal(program_Cassignments_list[0], program_Hassignments_list[0], user_assignments_list[1], ignoreH=ignoreH) == True):
            dist00 = np.sqrt((program_Hassignments_list[0][1]-user_assignments_list[0][4])**2 + (((program_Cassignments_list[0][1]-user_assignments_list[0][2])/6)**2))
            dist01 = np.sqrt((program_Hassignments_list[0][1]-user_assignments_list[1][4])**2 + (((program_Cassignments_list[0][1]-user_assignments_list[1][2])/6)**2))
            if dist00 <= dist01:
                user_assignments_list[0][5] = True     # this must update automatically user_VAL_LEU_HCNH_resid_assignments_dict/user_VAL_LEU_TOCSY_resid_assignments_dict
            elif dist00 > dist01:
                user_assignments_list[1][5] = True     # this must update automatically user_VAL_LEU_HCNH_resid_assignments_dict/user_VAL_LEU_TOCSY_resid_assignments_dict
            if " <CORRECT>" not in program_Cassignments_list[0][7]: program_Cassignments_list[0][7] += " <CORRECT>" # + iteration_comment
            if " <CORRECT>" not in program_Hassignments_list[0][7]: program_Hassignments_list[0][7] += " <CORRECT>" # + iteration_comment
        elif are_peaks_equal(program_Cassignments_list[0], program_Hassignments_list[0], user_assignments_list[0], ignoreH=ignoreH) == True:
            user_assignments_list[0][5] = True     # this must update automatically user_VAL_LEU_HCNH_resid_assignments_dict/user_VAL_LEU_TOCSY_resid_assignments_dict
            if " <CORRECT>" not in program_Cassignments_list[0][7]: program_Cassignments_list[0][7] += " <CORRECT>" # + iteration_comment
            if " <CORRECT>" not in program_Hassignments_list[0][7]: program_Hassignments_list[0][7] += " <CORRECT>" # + iteration_comment
        elif are_peaks_equal(program_Cassignments_list[0], program_Hassignments_list[0], user_assignments_list[1], ignoreH=ignoreH) == True:
            user_assignments_list[1][5] = True     # this must update automatically user_VAL_LEU_HCNH_resid_assignments_dict/user_VAL_LEU_TOCSY_resid_assignments_dict
            if " <CORRECT>" not in program_Cassignments_list[0][7]: program_Cassignments_list[0][7] += " <CORRECT>" # + iteration_comment
            if " <CORRECT>" not in program_Hassignments_list[0][7]: program_Hassignments_list[0][7] += " <CORRECT>" # + iteration_comment
        else:
            if " <WRONG>" not in program_Cassignments_list[0][7]: program_Cassignments_list[0][7] += " <WRONG>" # + iteration_comment
            if " <WRONG>" not in program_Hassignments_list[0][7]: program_Hassignments_list[0][7] += " <WRONG>" # + iteration_comment


        if ( user_assignments_list[0][5] == False and are_peaks_equal(program_Cassignments_list[1], program_Hassignments_list[1], user_assignments_list[0], ignoreH=ignoreH) == True and
            user_assignments_list[1][5] == False and are_peaks_equal(program_Cassignments_list[1], program_Hassignments_list[1], user_assignments_list[1], ignoreH=ignoreH) == True ):
            dist10 = np.sqrt((program_Hassignments_list[1][1]-user_assignments_list[0][4])**2 + (((program_Cassignments_list[1][1]-user_assignments_list[0][2])/6)**2))
            dist11 = np.sqrt((program_Hassignments_list[1][1]-user_assignments_list[1][4])**2 + (((program_Cassignments_list[1][1]-user_assignments_list[1][2])/6)**2))
            if dist10 <= dist11:
                user_assignments_list[0][5] = True     # this must update automatically user_VAL_LEU_HCNH_resid_assignments_dict/user_VAL_LEU_TOCSY_resid_assignments_dict
            elif dist10 > dist11:
                user_assignments_list[1][5] = True     # this must update automatically user_VAL_LEU_HCNH_resid_assignments_dict/user_VAL_LEU_TOCSY_resid_assignments_dict
            if " <CORRECT>" not in program_Cassignments_list[1][7]: program_Cassignments_list[1][7] += " <CORRECT>" # + iteration_comment
            if " <CORRECT>" not in program_Hassignments_list[1][7]: program_Hassignments_list[1][7] += " <CORRECT>" # + iteration_comment
        elif user_assignments_list[0][5] == False and are_peaks_equal(program_Cassignments_list[1], program_Hassignments_list[1], user_assignments_list[0], ignoreH=ignoreH) == True:
            user_assignments_list[0][5] = True     # this must update automatically user_VAL_LEU_HCNH_resid_assignments_dict/user_VAL_LEU_TOCSY_resid_assignments_dict
            if " <CORRECT>" not in program_Cassignments_list[1][7]: program_Cassignments_list[1][7] += " <CORRECT>" # + iteration_comment
            if " <CORRECT>" not in program_Hassignments_list[1][7]: program_Hassignments_list[1][7] += " <CORRECT>" # + iteration_comment
        elif user_assignments_list[1][5] == False and are_peaks_equal(program_Cassignments_list[1], program_Hassignments_list[1], user_assignments_list[1], ignoreH=ignoreH) == True:
            user_assignments_list[1][5] = True     # this must update automatically user_VAL_LEU_HCNH_resid_assignments_dict/user_VAL_LEU_TOCSY_resid_assignments_dict
            if " <CORRECT>" not in program_Cassignments_list[1][7]: program_Cassignments_list[1][7] += " <CORRECT>" # + iteration_comment
            if " <CORRECT>" not in program_Hassignments_list[1][7]: program_Hassignments_list[1][7] += " <CORRECT>" # + iteration_comment
        else:
            if " <WRONG>" not in program_Cassignments_list[1][7]: program_Cassignments_list[1][7] += " <WRONG>" # + iteration_comment
            if " <WRONG>" not in program_Hassignments_list[1][7]: program_Hassignments_list[1][7] += " <WRONG>" # + iteration_comment

    print("DEBUG compare_Cresonances: program_Cassignments_list=", program_Cassignments_list)
    print("DEBUG compare_Cresonances: program_Hassignments_list=", program_Hassignments_list)
    print("DEBUG compare_Cresonances: program_assignments_list=", program_assignments_list)
    print("DEBUG compare_Cresonances: user_assignments_list=", user_assignments_list)
    return user_assignments_list

def proof_read_equivalentMethyls(program_assignments_list,
                                 resid,
                                 user_special_HCNH_resid_assignments_dict,
                                 user_special_TOCSY_resid_assignments_dict,
                                 Carbon_pair,
                                 Proton_pair):
    """
    Method to compare the program- and use-assigned LEU and VAL methyl carbons and to proof read them.
    Currently works only in cases where both program assignments come from HCNH of both come from TOCSY.
    It does not work with hybrid HCNH/TOCSY assignments.

    ARGS:
    Carbon_pair:        list of the 2 equivalent Carbons, e.g. ['CG1', 'CG2'] in VAL, etc.
    user_special_HCNH_resid_assignments_dict:  in iteration1 & 2 it is user_VAL_LEU_HCNH_resid_assignments_dict
    user_special_TOCSY_resid_assignments_dict:  in iteration1 & 2 it is user_VAL_LEU_TOCSY_resid_assignments_dict
    program_assignments_list:     list of peaks assigned by the program to the equivalent methyl carbons of all
                                  LEU (CD1,CD2) and VAL (CG1,CG2). E.g.
                                  [['L217', 'CD1', 25.857, 'QD1', 0.855, False], ['L217', 'CD2', 23.09, 'QD2', 0.821, False]]
    """

    print("DEBUG proof_read_equivalentMethyls: program_assignments_list=", program_assignments_list)
    print("DEBUG proof_read_equivalentMethyls: user_special_HCNH_resid_assignments_dict=", user_special_HCNH_resid_assignments_dict)
    print("DEBUG proof_read_equivalentMethyls: user_special_TOCSY_resid_assignments_dict=", user_special_TOCSY_resid_assignments_dict)
    print("DEBUG proof_read_equivalentMethyls: Carbon_pair=", Carbon_pair)
    print("DEBUG proof_read_equivalentMethyls: Proton_pair=", Proton_pair)
    Carbon2Proton_mdict = tree()
    Carbon2Proton_mdict["LEU"]["CD1"] = "QD1"
    Carbon2Proton_mdict["LEU"]["CD2"] = "QD2"
    Carbon2Proton_mdict["LEU"]["CG1"] = "QG1"
    Carbon2Proton_mdict["LEU"]["CG2"] = "QG2"
    Carbon2Proton_mdict["TYR"]["CD1"] = "QD1"
    Carbon2Proton_mdict["PHE"]["CD1"] = "QD1"
    Carbon2Proton_mdict["PHE"]["CE1"] = "QE1"
    aa_type = program_assignments_list[0][6]


    HCNH_peak_num = 0
    TOCSY_peak_num = 0
    for assignment in program_assignments_list:
        spectrum_combo = re.sub(r" iter[0-9]$", "", assignment[7])   # applies to the case 'HCNH iter1'
        if spectrum_combo in ['HCNH',
                             'HCNH (C-term hanging residue)',
                             'HCNH (flanked residue)',
                             'HCNH + TOCSY-HCNH',
                             'TOCSY-HCNH',
                             'TOCSY-HCNH + HCNH',
                             'TOCSY (unmatched) + TOCSY-HCNH',
                             'TOCSY-HCNH + TOCSY (unmatched)']:    # if this resonance comes only from HCNH
            HCNH_peak_num += 1
        elif spectrum_combo in ['TOCSY (average)', 'TOCSY (unmatched)']:
            TOCSY_peak_num += 1
    user_assignments_list = []    # list with the C and H assignments
    if HCNH_peak_num >= 1 and TOCSY_peak_num == 0:     # if both assigned peaks come from HCNH, proof-read them
        if resid in list(user_special_HCNH_resid_assignments_dict.keys()):   # if it has been assigned from the user
            for assignment in [a for a in user_special_HCNH_resid_assignments_dict[resid] if a[1] in Carbon_pair+Proton_pair]: # iterate over the Carbon pair and Proton pair assignments
                if int(assignment[0][1:]) == resid:
                    user_assignments_list.append(assignment)
            user_assignments_list = compare_equivalentMethyl_peaks(program_assignments_list, user_assignments_list)
    elif TOCSY_peak_num >= 1 and HCNH_peak_num == 0:   # if both assigned peaks come from TOCSY, proof-read them
        if resid in list(user_special_TOCSY_resid_assignments_dict.keys()):   # if it has been assigned from the user
            for assignment in [a for a in user_special_TOCSY_resid_assignments_dict[resid] if a[1] in Carbon_pair+Proton_pair]:
                if int(assignment[0][1:]) == resid:
                    user_assignments_list.append(assignment)
            user_assignments_list = compare_equivalentMethyl_peaks(program_assignments_list, user_assignments_list)
    elif TOCSY_peak_num == 1 and HCNH_peak_num == 1:   # if one assigned peak comes from TOCSY and the other from HCNH
        program_assignments_list[0][7] += " <CHECK MANUALLY>"
        program_assignments_list[1][7] += " <CHECK MANUALLY>"

    # UPDATE THE ASSIGNMENTS BY APPENDING <CORRECT>/<WRONG> COMMENTS AT THE END
    print("DEBUG proof_read_equivalentMethyls: program_assignments_list=", program_assignments_list)
    print("DEBUG proof_read_equivalentMethyls: user_assignments_list=", user_assignments_list)
    # compare_Cresonances(program_assignments_list, user_assignments_list) # WHY DID I WRITE THIS AGAIN HERE ???


def match_TOCSY_to_HCNH_lines(TOCSY_residue_assignments_dict,
                               HCNH_residue_assignments_dict,
                               protein_alignment_list,
                               absolute_matches_alignment_list,
                               args):

    """
    FUNCTION that matches each assignment made in TOCSY to the respective peak of the HCNH. If the lower distance of a TOCSY peak from any HCNH peak is
    greater than (0.04, 0.4) then it will be considered as garbage and will be discarded from TOCSYnum.list file. Note that the keys of the input
    dictionaries are AAIGs while the keys of the output dictionaries are the actual residues of the protein

    ARGS:
    TOCSY_residue_assignments_dict:     the TOCSY peak assignments, residue -> [ (Cname, Cresonance, Hname, Hresonance), (Cname, Cresonance, Hname, Hresonance), ... ]
    HCNH_residue_assignments_dict:     the unassigned HCNH peaks, residue -> [ ("?", "?", Cresonance, "?", Hresonance), ("?", "?", Cresonance, "?", Hresonance), ... ]

    RETURNS:
    HCNH_residue_assignments_dict:     the input HCNH_residue_assignments_dict but with the TOCSY C & H assignments that corresponded to each peak.
    """
    tolH = 0.04
    tolC = 0.4

    # ATTENTION: recall that in TOCSY_residue_assignments_dict both the keys and the C,H assignments correspond to residue i-1.
    # But in matched_HCNH_residue_assignments_dict the keys are residue i and the values are peaks from residues i (mostly the strongest), residue i-1 and maybe i+1 or
    # other neighboring residues.
    # So in order to match the HCNH peaks of residue i we need primarily the values of TOCSY_residue_assignments_dict keys i respectively.

    #print "DEBUG: point 3 TOCSY_residue_assignments_dict=", TOCSY_residue_assignments_dict
    print("DEBUG: HCNH_residue_assignments_dict=", HCNH_residue_assignments_dict)
    residue_unmatched_TOCSY_peaks_dict = {}  # residues -> unmatched TOCSY peaks
    for i_residue in list(HCNH_residue_assignments_dict.keys()):
        # STUPID PYTHON, if I do: TOCSY_assigned_peaks_list = TOCSY_residue_assignments_dict[i_residue], TOCSY_assigned_peaks_list will be a reference of TOCSY_residue_assignments_dict[i_residue], not a new list!!!
        try:
            TOCSY_assigned_peaks_list = [x for x in TOCSY_residue_assignments_dict[i_residue]]  # [ (i_residue, TOCSY_Cname, HCNH_Creson, TOCSY_Hname, HCNH_Hreson), ..., (i_residue, TOCSY_Cname, HCNH_Creson, TOCSY_Hname, HCNH_Hreson) ]
        except KeyError:    # this i_residue does not have TOCSY C-H peaks
            TOCSY_assigned_peaks_list = []
        try:   # CONSIDER THE ASSIGNED i-1 TOCSY PEAKS, TOO
            iminus1_residue = get_iminus1_residue_from_i_residue(i_residue, protein_alignment_list, absolute_matches_alignment_list, args)
            if iminus1_residue in list(TOCSY_residue_assignments_dict.keys()):    # this will give C,H assignments of residue i-1. If they exist, include them in the matching TOCSY-->HCNH
                TOCSY_assigned_peaks_list.extend(TOCSY_residue_assignments_dict[iminus1_residue])    # now we have TOCSY C,H assignments for both residues i and i-1
                # so TOCSY_assigned_peaks_list = [ (iminus1_residue, TOCSY_Cname, HCNH_Creson, TOCSY_Hname, HCNH_Hreson), ..., (i_residue, TOCSY_Cname, HCNH_Creson, TOCSY_Hname, HCNH_Hreson),... ]
        except ValueError:  # there are no TOCSY C,H assignments for residue i-1, so skip it
            print("WARNING: there are no TOCSY C,H assignments for the residue downstream of ",i_residue,", so skip it")
            pass
        try:   # CONSIDER THE ASSIGNED i+1 TOCSY PEAKS, TOO
            print("DEBUG: i_residue=", i_residue)
            iplus1_residue = get_iplus1_residue_from_i_residue(i_residue, protein_alignment_list, absolute_matches_alignment_list, args)
            print("DEBUG: matching assigned TOCSY PEAKS of i+1 residue", iplus1_residue)
            if iplus1_residue in list(TOCSY_residue_assignments_dict.keys()):    # this will give C,H assignments of residue i+1. If they exist, include them in the matching TOCSY-->HCNH
                TOCSY_assigned_peaks_list.extend(TOCSY_residue_assignments_dict[iplus1_residue])    # now we have TOCSY C,H assignments for residues i, i-1 and i+1
                print("DEBUG: addinge i+1 residue TOCSY PEAK assignments to TOCSY_assigned_peaks_list=", TOCSY_assigned_peaks_list)
                # so TOCSY_assigned_peaks_list = [ (iplus1_residue, TOCSY_Cname, HCNH_Creson, TOCSY_Hname, HCNH_Hreson), ..., (i_residue, TOCSY_Cname, HCNH_Creson, TOCSY_Hname, HCNH_Hreson),... ]
        except (ValueError, IndexError):  # there are no TOCSY C,H assignments for residue i+1, so skip it
            print("WARNING: there are no TOCSY C,H assignments for the residue upstream of ",i_residue,", so skip it")
            pass
        for TOCSY_peak in TOCSY_assigned_peaks_list:
            TOCSY_residue = TOCSY_peak[0]   # either residue i-1 or i or i+1
            TOCSY_Cname = TOCSY_peak[1]
            TOCSY_Creson = TOCSY_peak[2]
            TOCSY_Hname = TOCSY_peak[3]
            TOCSY_Hreson = TOCSY_peak[4]
            if not i_residue in list(HCNH_residue_assignments_dict.keys()):    # DANGEROUS but works for the N-term (A200 in aLP that has not HCNH peaks)
                print("WARNING: residue", i_residue," does not  have any HCNH match (e.g. N-term or Proline)!")
                continue
            for Nindex in range(0, len(HCNH_residue_assignments_dict[i_residue])):   # includes peaks from residues i (strongest), i-1 and the neighboring residues
                HCNH_peak = HCNH_residue_assignments_dict[i_residue][Nindex]
                HCNH_residue = HCNH_peak[0]
                HCNH_Cname = HCNH_peak[1]
                HCNH_Creson = HCNH_peak[2]
                HCNH_Hname = HCNH_peak[3]
                HCNH_Hreson = HCNH_peak[4]
                # First ignore the garbage that are far away from any peak
                print("DEBUG: trying to match TOCSY_peak=", TOCSY_peak, "to HCNH_peak=", HCNH_peak)
                if ( (HCNH_Hreson -tolH) <= TOCSY_Hreson <= (HCNH_Hreson +tolH) ) and ( (HCNH_Creson -tolC) <= TOCSY_Creson <= (HCNH_Creson +tolC) ):
                    print("DEBUG: this HCNH peak satisfies the C & H tolerances:", HCNH_peak)
                    # if HCNH_residue == "?" and HCNH_Cname == "?" and HCNH_Hname == "?":   # if this HCNH peak does not have C,H assignments already, assign the currect C,H names to it
                    if HCNH_residue == "?":   # if this HCNH peak does not have C,H assignments already, assign the currect C,H names to it
                        print("DEBUG: this HCNH peak does not have C,H assignments already, therefore assigning the currect C,H names to it:", HCNH_peak)
                        has_been_matched = False # if this TOCSY peak has been matched already
                        for Ni in range(len(HCNH_residue_assignments_dict[i_residue])):    # first check if we already matched the current TOCSY peak to any other HCNH peak Ni
                            if HCNH_residue_assignments_dict[i_residue][Ni][0] == TOCSY_residue and HCNH_residue_assignments_dict[i_residue][Ni][1] == TOCSY_Cname and HCNH_residue_assignments_dict[i_residue][Ni][3] == TOCSY_Hname:
                                has_been_matched = True
                                Ni_Creson = HCNH_residue_assignments_dict[i_residue][Ni][2]
                                Ni_Hreson = HCNH_residue_assignments_dict[i_residue][Ni][4]
                                delta_Ni = np.sqrt((TOCSY_Hreson - Ni_Hreson)**2 + ((TOCSY_Creson - Ni_Creson)/6)**2)   # distance between the current TOCSY peak and the previously matched HCNH peak
                                delta_Nindex = np.sqrt((TOCSY_Hreson - HCNH_Hreson)**2 + ((TOCSY_Creson - HCNH_Creson)/6)**2) # distance between the current TOCSY peak and the current HCNH peak
                                if delta_Nindex < delta_Ni:  # label the Nindex HCNH peak and unlabel the Ni HCNH peak
                                    HCNH_residue_assignments_dict[i_residue][Ni][0] = "?"
                                    HCNH_residue_assignments_dict[i_residue][Ni][1] = "?"
                                    HCNH_residue_assignments_dict[i_residue][Ni][3] = "?"
                                    HCNH_residue_assignments_dict[i_residue][Nindex] = [TOCSY_residue, TOCSY_Cname, HCNH_Creson, TOCSY_Hname, HCNH_Hreson]
                                else:   # if the previously matched Ni HCNH peak is closer to the TOCSY peak than the current Nindex HCNH peak, continue
                                    continue
                        if has_been_matched == False:   # if this TOCSY peak has not been matched with any HCNH peak, match it with the current one (Nindex)
                            HCNH_residue_assignments_dict[i_residue][Nindex] = [TOCSY_residue, TOCSY_Cname, HCNH_Creson, TOCSY_Hname, HCNH_Hreson]
                    else: # if it has already been assigned C,H names, check if the currect TOCSY peak has closer C & H resonances
                        print("DEBUG: this HCNH peak has already been assigned C,H names:", HCNH_peak)
                        for tmp_TOCSY_peak in TOCSY_assigned_peaks_list:
                            tmp_TOCSY_residue = TOCSY_peak[0]
                            tmp_TOCSY_Cname = TOCSY_peak[1]
                            tmp_TOCSY_Creson = TOCSY_peak[2]
                            tmp_TOCSY_Hname = TOCSY_peak[3]
                            tmp_TOCSY_Hreson = TOCSY_peak[4]
                            if tmp_TOCSY_residue == HCNH_residue and tmp_TOCSY_Cname == HCNH_Cname and tmp_TOCSY_Hname == HCNH_Hname:   # we found the assigned TOCSY C,H peak that was previously metched to this HCNH peak
                                # if the current TOCSY peak has C & H resonances closer to this HCNH peak than the previously assigned TOCSY peak, change the C,H assignments of this HCNH peak
                                print("DEBUG: HCNH_Creson=", HCNH_Creson, "tmp_TOCSY_Creson", tmp_TOCSY_Creson)
                                print("DEBUG: HCNH_Hreson=", HCNH_Hreson, "tmp_TOCSY_Hreson=", tmp_TOCSY_Hreson)
                                assigned_TOCSY_Creson = tmp_TOCSY_Creson    # just for clarity, rename the resonances
                                assigned_TOCSY_Hreson = tmp_TOCSY_Hreson
                                delta_assigned_TOCSY = np.sqrt((assigned_TOCSY_Hreson - HCNH_Hreson)**2 + ((assigned_TOCSY_Creson - HCNH_Creson)/6)**2)
                                delta_TOCSY = np.sqrt((TOCSY_Hreson - HCNH_Hreson)**2 + ((TOCSY_Creson - HCNH_Creson)/6)**2)
                                if delta_assigned_TOCSY > delta_TOCSY:  # replace the previous C,H assignment
                                    # MEANWHILE CLEAN THE PREVIOUS HCNH PEAK TO WHICH YOU ASSIGNED THIS TOCSY PEAK
                                    for Ni in range(len(HCNH_residue_assignments_dict[i_residue])):
                                        if HCNH_residue_assignments_dict[i_residue][Ni][0] == TOCSY_residue and HCNH_residue_assignments_dict[i_residue][Ni][1] == TOCSY_Cname and HCNH_residue_assignments_dict[i_residue][Ni][3] == TOCSY_Hname:
                                            HCNH_residue_assignments_dict[i_residue][Ni][0] = "?"
                                            HCNH_residue_assignments_dict[i_residue][Ni][1] = "?"
                                            HCNH_residue_assignments_dict[i_residue][Ni][3] = "?"
                                    print("DEBUG: substituting HCNH peak", HCNH_residue_assignments_dict[i_residue][Nindex], " for ", [TOCSY_residue, TOCSY_Cname, HCNH_Creson, TOCSY_Hname, HCNH_Hreson])
                                    print("DEBUG: the updated HCNH_residue_assignments_dict[i_residue]=", HCNH_residue_assignments_dict[i_residue])
                                    HCNH_residue_assignments_dict[i_residue][Nindex] = [TOCSY_residue, TOCSY_Cname, HCNH_Creson, TOCSY_Hname, HCNH_Hreson]
                                # otherwise keep the previous AAIG group
                                # ATTENTION: THIS IS NOT ENTIRELY CORRECT, BECAUSE THE C RESONANCE MAY BE CLOSER BUT THE H NOT!
                                break

        # NOW CHECK WHICH TOCSY PEAKS COULD NOT BE MATCHED WITH ANY HCNH PEAK
        print("DEBUG: HCNH_residue_assignments_dict=", HCNH_residue_assignments_dict)
        for TOCSY_peak in TOCSY_assigned_peaks_list:
            TOCSY_residue = TOCSY_peak[0]
            TOCSY_Cname = TOCSY_peak[1]
            TOCSY_Hname = TOCSY_peak[3]
            if TOCSY_residue == i_residue:   # print only unmatched i residue peaks
                if not TOCSY_residue in list(HCNH_residue_assignments_dict.keys()):
                    print("WARNING residue", TOCSY_residue, "has no HCNH peaks (e.g. N-term or Proline)!")
                    continue
                HCNH_matched_peak_list = [peak for peak in HCNH_residue_assignments_dict[i_residue] if peak[0]==TOCSY_residue and peak[1]==TOCSY_Cname and peak[3]==TOCSY_Hname]
                if len(HCNH_matched_peak_list) == 0:
                    print("WARNING: the following peak of residue", i_residue, " could not be matched with any HCNH peak:")
                    print(TOCSY_peak)
                    try:
                        residue_unmatched_TOCSY_peaks_dict[i_residue].append(TOCSY_peak)
                    except KeyError:
                        residue_unmatched_TOCSY_peaks_dict[i_residue] = [TOCSY_peak]

    #print "DEBUG: point 4 TOCSY_residue_assignments_dict=", TOCSY_residue_assignments_dict
    return HCNH_residue_assignments_dict, residue_unmatched_TOCSY_peaks_dict


def transform_intensity(intensity, tans_mode=1):
    print("DEBUG: intensity=", intensity)
    if tans_mode == -1:
        return 1
    elif tans_mode == 1:
        return intensity
    elif tans_mode == 2:
        return 100*intensity**2
    elif tans_mode == 3:
        return 1000*intensity**3
    elif tans_mode == 4:
        return 10000*intensity**4
    elif tans_mode == 5:
        return 100000*intensity**5
    elif tans_mode == 6:
        return 1000000*intensity**6
    elif tans_mode == 7:
        return 10000000*intensity**7
    elif tans_mode == 10:
        return 10000000000*intensity**10


def get_probabilities_from_H_C_resonpair_2Dhist(residue,
                                                Hreson,
                                                Creson,
                                                missing_carbons_list,
                                                partly_assigned_carbons_list,
                                                args,
                                                HCNH_residue_peak_intensity_mdict,
                                                matched_i_to_iplus1_peaks_mdict,
                                                aa2pthres_dict,
                                                aa2ithres_dict,
                                                protein_alignment_list,
                                                absolute_matches_alignment_list,
                                                histload,
                                                iteration=1,
                                                I_THRESHOLD=None):
    """
        Method to find and the respective probabilities from a pair of aliphatic H,C resonances.

        ARGUMENTS:
        missing_carbons_list:            list of carbons which have either all or one of their protons missing from TOCSY

        RETURNS:
        matches_list: list of lists of the form [aa type, Carbon name, Hydrogen name, 2D-hist C-H probability, H resonance, C resonance]
    """
    from .global_vars import percentile_mdict
    Percentile_mdict = percentile_mdict
    histcalc = ProbHist_Calculator()

    # print "DEBUG: HCNH_residue_peak_intensity_mdict=", HCNH_residue_peak_intensity_mdict

    aa_type = aa1to3_dict[get_aa_type_from_residue(residue)]
    if iteration == 1:
        if aa_type in list(aa2pthres_dict.keys()):
            PERCENTILE = aa2pthres_dict[aa_type]
        else:
            PERCENTILE = args.PERCENTILE_ITER1
        TRANSFORMATION_TYPE = args.INTENSITY_TRANSFORM_TYPE_ITER1
        if I_THRESHOLD:
            INTENSITY_THRESHOLD = I_THRESHOLD
        else:
            INTENSITY_THRESHOLD = args.INTENSITY_THRESHOLD_ITERATION1
    elif iteration ==2:
        PERCENTILE = args.PERCENTILE_ITER2
        TRANSFORMATION_TYPE = args.INTENSITY_TRANSFORM_TYPE_ITER2
        INTENSITY_THRESHOLD = args.INTENSITY_THRESHOLD_ITERATION2
    elif iteration ==3:
        if aa_type in list(aa2pthres_dict.keys()):
            PERCENTILE = aa2pthres_dict[aa_type]
        else:
            PERCENTILE = args.PERCENTILE_ITER3
        TRANSFORMATION_TYPE = args.INTENSITY_TRANSFORM_TYPE_ITER3
        if I_THRESHOLD:
            INTENSITY_THRESHOLD = I_THRESHOLD
        else:
            INTENSITY_THRESHOLD = args.INTENSITY_THRESHOLD_ITERATION3
    elif iteration ==4:
        PERCENTILE = args.PERCENTILE_ITER4
        TRANSFORMATION_TYPE = args.INTENSITY_TRANSFORM_TYPE_ITER4
        INTENSITY_THRESHOLD = args.INTENSITY_THRESHOLD_ITERATION4
    elif iteration ==5:
        PERCENTILE = args.PERCENTILE_ITER5
        TRANSFORMATION_TYPE = args.INTENSITY_TRANSFORM_TYPE_ITER5
        INTENSITY_THRESHOLD = args.INTENSITY_THRESHOLD_ITERATION5

    print("DEBUG: HCNH_residue_peak_intensity_mdict=", HCNH_residue_peak_intensity_mdict)
    if len([HCNH_residue_peak_intensity_mdict[r][p] for r in list(HCNH_residue_peak_intensity_mdict.keys()) for p in HCNH_residue_peak_intensity_mdict[r]]) == 0:    # If the HCNH file does not contain intensity column, set them all to 1
        TRANSFORMATION_TYPE = -1    ; # this means multiply by x1

    aa = aa1to3_dict[residue[0]]

    print("DEBUG: iteration=", iteration, "PERCENTILE=", PERCENTILE, "TRANSFORMATION_TYPE=", TRANSFORMATION_TYPE)
    #print "DEBUG: Hreson=",Hreson,"Creson=",Creson
    matches_list = [] # list of tuples of the form (aa type, Carbon name, Hydrogen name, overall probability, )
    carbon_matches_list = []    # list of tuples of the form (aa type,carbon,probability); here aa should be the same in all the tuples of the list!

    # for carbon in missing_carbons_list:    # iterate over all missing carbons of this residue
    #     print "DEBUG: carbon=", carbon
    #     bin_array, density_array = histload.aa_carbon_binDensityList_mdict[aa][carbon]
    #     probability = histcalc.get_probability_from_histogram(Creson, bin_array, density_array)
    #     #print "DEBUG: aa=",aa,"carbon=",carbon,"probability=",probability
    #     if probability > 0.0:
    #         carbon_matches_list.append((aa,carbon,probability))
    all_missing_carbons_list = missing_carbons_list + partly_assigned_carbons_list  #
    for CH_pair in list(histload.aa_CHpair_binProbabilityList_mdict[aa].keys()):
        #print "DEBUG: CH_pair.split('-')[0]=", CH_pair.split('-')[0], "missing_carbons_list =", missing_carbons_list
        if CH_pair.split('-')[0] in all_missing_carbons_list:   # get 2D hist probabilities for all missing carbons (partly or not)
            probability = 0.0
            if args.USE_2D_HISTOGRAMS == True:
                #print "DEBUG: checking 2D hist probability for aa=",aa, "CH_pair=", CH_pair, "Hreson=", Hreson, "Creson", Creson
                probability = histcalc.get_CH_probability_from_2Dhistogram(CH_pair, aa, Hreson, Creson, histload)
            elif args.USE_2D_HISTOGRAMS == False:
                Cname = CH_pair.split('-')[0]
                Hname = CH_pair.split('-')[1]
                print("DEBUG: aa=", aa, "CH_pair=", CH_pair)
                probability = histcalc.get_CH_probability_from_1Dhistograms(aa=aa,
                                                                            Hname=Hname,
                                                                            Hreson=Hreson,
                                                                            Cname=Cname,
                                                                            Creson=Creson,
                                                                            histload=histload,
                                                                            PROBABILITY_MODEL=args.PROBABILITY_MODEL,
                                                                            H_weight=args.H_weight,
                                                                            C_weight=args.C_weight)

            # Keep only the amino acids for which were found both the carbon and the respective covalently bonded hydrogen
            # The probability of each prediction is given by the 2D histogram
            if probability > 0.0:
                C_name, H_name = CH_pair.split("-")
                if aa in ["LYS", "ARG"] and CH_pair in ["CD-HD2", "CD-HD3", "CE-HE2", "CE-HE3"] and probability < Percentile_mdict[aa][CH_pair][0.9]:   # strict exception for LYS and ARG
                    print("DEBUG: the following prediction did not pass the percentile threshold:", [aa, C_name, H_name, probability, Hreson, Creson], "percentile=", Percentile_mdict[aa][CH_pair][0.9])
                elif probability >= Percentile_mdict[aa][CH_pair][PERCENTILE]:   # save this C-H type prediction, only if the 2D-hist probability is above the value of the designated percentile
                    # MULTIPLY BY THE INTENSITY AFTER THE PERCENTILE FILTER
                    if (iteration==1 and args.USE_INTENSITIES_ITERATION1) or (iteration==2 and args.USE_INTENSITIES_ITERATION2) or (iteration==3 and args.USE_INTENSITIES_ITERATION3) or (iteration==4 and args.USE_INTENSITIES_ITERATION4) or (iteration==5 and args.USE_INTENSITIES_ITERATION5):
                        # check if the intensity is above the threshold only if the HCNH file has intensity column (TRANSFORMATION_TYPE!=-1)
                        if TRANSFORMATION_TYPE != -1 and HCNH_residue_peak_intensity_mdict[residue][(Creson, Hreson)] < INTENSITY_THRESHOLD:
                            print("DEBUG: The relative intensity ", HCNH_residue_peak_intensity_mdict[residue][(Creson, Hreson)], "of residue", residue, " peak ", (Creson, Hreson), "is lower than the threshold ", INTENSITY_THRESHOLD)
                            continue
                        print("DEBUG: multiplying with normalized intensity residue", residue," peak",(Creson, Hreson) , HCNH_residue_peak_intensity_mdict[residue][(Creson, Hreson)], ". Probability before:", probability)
                        if iteration <= 2:
                            try:
                                iplus1_residue = get_iplus1_residue_from_i_residue(residue, protein_alignment_list, absolute_matches_alignment_list, args)
                                if residue not in matched_i_to_iplus1_peaks_mdict.keys() or (Creson, Hreson) not in matched_i_to_iplus1_peaks_mdict[residue].keys():
                                    raise ValueError
                                iplus1_peak = matched_i_to_iplus1_peaks_mdict[residue][(Creson, Hreson)]
                                print("DEBUG: iplus1_residue=", iplus1_residue, "iplus1_peak=", iplus1_peak, "residue=", residue, "(Creson, Hreson)=", (Creson, Hreson))
                                intensity_i = HCNH_residue_peak_intensity_mdict[residue][(Creson, Hreson)]
                                print("DEBUG: intensity_i=", intensity_i)
                                intensity_iplus1 = HCNH_residue_peak_intensity_mdict[iplus1_residue][iplus1_peak]
                                transformed_intensity_i = transform_intensity(intensity_i, TRANSFORMATION_TYPE)
                                transformed_intensity_iplus1 = transform_intensity(intensity_iplus1, TRANSFORMATION_TYPE)
                                print("DEBUG: transformed_intensity_i=", transformed_intensity_i, " transformed_intensity_iplus1=", transformed_intensity_iplus1)
                                probability *=  transformed_intensity_i * transformed_intensity_iplus1 # multiply the probability by the normalized intensity of this peak and that of i+1 residue
                            except ValueError:  # this peak does not exist in the HCNH of i+1 residue
                                probability *=  transform_intensity(HCNH_residue_peak_intensity_mdict[residue][(Creson, Hreson)], TRANSFORMATION_TYPE)
                            print("DEBUG: probability after:", probability)
                        elif iteration >= 3:
                            probability *=  transform_intensity(HCNH_residue_peak_intensity_mdict[residue][(Creson, Hreson)], TRANSFORMATION_TYPE)
                            print("DEBUG: probability after:", probability)
                    matches_list.append([aa, C_name, H_name, probability, Hreson, Creson])  # add also the index of this C-H pair in the TOCSY group it
                else:   # just for debugging
                    print("DEBUG: the following prediction did not pass the percentile threshold:", [aa, C_name, H_name, probability, Hreson, Creson], "percentile=", Percentile_mdict[aa][CH_pair][PERCENTILE])
            # elif probability == 0.0:    # if the probability of the 2D hist at this point is 0, use only the Carbon 1D histogram to get the probability
            #     print "WARNING: 2D hist probability is 0. Using Carbon 1D histogram to get the probability!"
            #     C_name, H_name = CH_pair.split("-")
            #     #print "DEBUG: C_name=", C_name, "carbon_matches_list=", carbon_matches_list
            #     carbon_match = [match for match in carbon_matches_list if match[1]==C_name]
            #     if len(carbon_match) == 0:  # if the hist(C) of this C_name was zero, don't save it
            #         continue
            #     #print "DEBUG: len(carbon_match) is ", len(carbon_match)," and it should be 1 !"
            #     #print "DEBUG: carbon_match=", carbon_match
            #     Carbon_probability = -1 * carbon_match[0][2] #  but first make it negative to distiguish it from weighted average probabilities
            #     #print "DEBUG: only the Carbon 1D histogram used:", [aa, C_name, H_name, weighted_average_probability, TOCSY_reson_index, Hreson, Creson, None]
            #     matches_list.append([aa, C_name, H_name, Carbon_probability, Hreson, Creson])

    #print "DEBUG: returning matches_list=", matches_list
    return matches_list


def split_oversized_Cgroups(possible_aatype_prob_C_H_resonpair_TAAIG_list_list):
    """
        FUNCTION to split Cgroups with more that 2 peaks into separate Cgroups. Currently works for maximum Cgroup size 3 only. This happens when
        two peaks have the same Carbon resonance but different proton resonance. So because C-grouping works with Carbon resonances only, all 3 peaks
        are included in the same group
    """

    clustID2size_dict = {}
    clustID2Creson_dict = defaultdict(set)
    for a in possible_aatype_prob_C_H_resonpair_TAAIG_list_list:
        clustID = a[7]
        Creson = a[6]
        Hreson = a[5]
        clustID2Creson_dict[clustID].add((Hreson, Creson))
    # TODO: there is a very special scenario where 2 peaks may have the same Creson and Hreson but different
    # TODO: Nreson and HNreson. Here we don't work with the latter two, therefore we cannot identify that case.

    for clustID in list(clustID2Creson_dict.keys()):
        clustID2size_dict[clustID] = len(clustID2Creson_dict[clustID])

    # print "DEBUG: clustID2size_dict=", clustID2size_dict
    maxClustID = max(clustID2size_dict.keys())
    for clustID in list(clustID2size_dict.keys()):
        if clustID2size_dict[clustID] == 3:
            ColorPrint("WARNING: detected carbon cluster with 3 peaks! Trying to split it.\n%s" %
                       [a for a in possible_aatype_prob_C_H_resonpair_TAAIG_list_list if a[7]==clustID],
                       "WARNING")
            # ATTENTION: possible_aatype_prob_C_H_resonpair_TAAIG_list_list contains multiple alternative
            # assignments of the same peaks. Therefore, use clustID2Creson_dict to get the 3 Cresons safely.
            Creson_list = [t[1] for t in clustID2Creson_dict[clustID]]
            dist12 = abs(Creson_list[0] - Creson_list[1])
            dist13 = abs(Creson_list[0] - Creson_list[2])
            dist23 = abs(Creson_list[1] - Creson_list[2])
            if len(set(Creson_list)) == 2:  # special case where two peaks have exactly the same Creson
                for Creson in Creson_list:
                    if Creson_list.count(Creson) == 1:
                        Coutlier = Creson
                        break
            elif dist12 < dist13 and dist12 < dist23:   # change the Cgroup of peak3
                Coutlier = Creson_list[2]
            elif dist13 < dist12 and dist13 < dist23:   # change the Cgroup of peak2
                Coutlier = Creson_list[1]
            elif dist23 < dist12 and dist23 < dist13:   # change the Cgroup of peak1
                Coutlier = Creson_list[0]

            # The following is common in all cases (change the outlier peak's Cgroup)
            for i in range(len(possible_aatype_prob_C_H_resonpair_TAAIG_list_list)):
                a = possible_aatype_prob_C_H_resonpair_TAAIG_list_list[i]
                try:
                    if a[6] == Coutlier:
                        # print "DEBUG: changing Cgroup of ", a, "to ", maxClustID + 1
                        a[7] = maxClustID + 1
                except UnboundLocalError:
                    raise Exception(Debuginfo("FAIL: there is a new case of oversized C-groups which you did "
                                                   "not think of!", fail=True))

        elif clustID2size_dict[clustID] > 3:
            ColorPrint("ERROR: found Cgroup with more than 3 peaks!!!", "FAIL")
            ColorPrint("DEBUG: possible_aatype_prob_C_H_resonpair_TAAIG_list_list=%s" %
                  possible_aatype_prob_C_H_resonpair_TAAIG_list_list, "FAIL")
            sys.exit(1)

    print("DEBUG: after oversized carbon cluster splitting possible_aatype_prob_C_H_resonpair_TAAIG_list_list=%s" %
          possible_aatype_prob_C_H_resonpair_TAAIG_list_list)
    return possible_aatype_prob_C_H_resonpair_TAAIG_list_list


def group_carbons(possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list, iteration=None):
    """
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
    updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list = []    # same as possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list but with the Cluster IDs for each aa type
    for aa_type in aa_types_set:
        singleAAtype_prob_C_H_resonpair_TOCSYindex_list_list = [x for x in possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list if x[0]==aa_type]
        if not aa_type in list(aa_carbonBondedGeminalHydrogensDict_dict.keys()):  # this aa does not have geminal protons, skip Clustering
            clustID = 0
            for index in range(len(singleAAtype_prob_C_H_resonpair_TOCSYindex_list_list)):
                clustID += 1
                singleAAtype_prob_C_H_resonpair_TOCSYindex_list_list[index][7] = clustID
            # save the updated list entries with the correct Cluster IDs
            updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list.extend(singleAAtype_prob_C_H_resonpair_TOCSYindex_list_list)
            continue    # continue with the next possible aa type

        new_aaindex_carbonGroupsTuple_dict = {}
        #print "DEBUG: before singleAAtype_prob_C_H_resonpair_TOCSYindex_list_list=", singleAAtype_prob_C_H_resonpair_TOCSYindex_list_list
        carbonCS_list = [l[6] for l in singleAAtype_prob_C_H_resonpair_TOCSYindex_list_list]
        #print "DEBUG: carbonCS_list=", carbonCS_list
        cl = HierarchicalClustering(carbonCS_list , lambda x,y: float(abs(x-y)))

        CONTINUE_CLUSTERING = True
        if iteration == None:   # If TOCSY, set the cutoff to 3.0
            cutoff = 0.300001       # For TOCSY it must be 3.0 otherwise raises "ERROR: found Cgroup with more than 3 peaks!!!" in MS6282.
                                    # For HCNH it can be 2.0 without problem (default value).
        else:
            cutoff = 0.200001
        while CONTINUE_CLUSTERING:  # continue Clustering by lowering the cutoff until no more than 2 carbons resonances are within each Cluster
            CONTINUE_CLUSTERING = False
            Cluster_list = cl.getlevel(cutoff)    # <== CHANGE ME (the tolerance for carbon)
            #print "DEBUG: Cluster_list=", Cluster_list
            for Cluster in Cluster_list:
                try:
                    if len(set(Cluster)) > 2:
                        CONTINUE_CLUSTERING = True
                        cutoff -= 0.02
                        break
                except TypeError:   # in case only one nucleus prediction is available (and only one peak)
                    CONTINUE_CLUSTERING = False

        for index in range(len(singleAAtype_prob_C_H_resonpair_TOCSYindex_list_list)):
            carbon_resonance = singleAAtype_prob_C_H_resonpair_TOCSYindex_list_list[index][6]
            for clustID, clust in enumerate(Cluster_list):
                #print "DEBUG: clustID=", clustID, "clust=", clust, "carbon_resonance=", carbon_resonance
                if len(Cluster_list) == 1 and carbon_resonance == clust:    # to avoid "TypeError: argument of type 'float' is not iterable"
                    #print "DEBUG: assigning Cluster", clustID+1
                    singleAAtype_prob_C_H_resonpair_TOCSYindex_list_list[index][7] = clustID + 1  # because enumeration starts from 0
                elif carbon_resonance in clust:
                    #print "DEBUG: assigning Cluster", clustID+1
                    singleAAtype_prob_C_H_resonpair_TOCSYindex_list_list[index][7] = clustID + 1  # because enumeration starts from 0
        # save the updated list entries with the correct Cluster IDs
        updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list.extend(singleAAtype_prob_C_H_resonpair_TOCSYindex_list_list)
        #print "DEBUG group_carbons: appended carbon Clusters with ", singleAAtype_prob_C_H_resonpair_TOCSYindex_list_list

    ## In the extreme scenario where a Cgroup has >2 peaks, keep the 2 peaks with closest Carbon resonance and the change
    # the Cgroup of the others (Currently works only for 3 peaks in the same Cgroup).
    updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list = split_oversized_Cgroups(updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list)

    print("DEBUG: point 1 updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list=", updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list)
    if len(possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list[0]) == 9:    # if this is HCNH assignment (the 9th element is the spectrum_combo)
        # For every Cname that is marked as 'TOCSY (unmatched)', search for another one Cname like this which is marked as 'TOCSY-HCNH',
        # and change its ClusterID. This should happen because one Creson comes from HCNH and the other from TOCSY, therefore they could not
        # be grouped together with the current cutoff. This needs to be modified when I will search the second missing C-H peak in HCNH...
        matched_peaks = [x for x in possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list if x[8]=='TOCSY-HCNH']
        matched_Cnames = [x[1] for x in matched_peaks]
        matched_clustIDs = [x[7] for x in matched_peaks]
        for Cname, clustID in zip(matched_Cnames, matched_clustIDs):
            for index in range(len(updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list)):
                peak = updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list[index]
                if peak[8] in ['TOCSY-HCNH', 'TOCSY (unmatched)'] and peak[1]==Cname:
                    updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list[index][7] = clustID # change the Cluster ID

    # Make sure that the same peaks (but with different predictions) belong to the same Cluster ID
    print("DEBUG: point 2 updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list=", updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list)
    CH_reson_tuple_set = set([(p[5],p[6]) for p in updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list])    # all unique C-H resonances
    for CH_reson_tuple in CH_reson_tuple_set:
        Creson = CH_reson_tuple[1]
        Hreson = CH_reson_tuple[0]
        # get all the assigned Cluster IDs of this C-H peak (normally they should be one)
        ClusterIDs_set = set([p[7] for p in updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list if p[5]==Hreson and p[6]==Creson])
        print("DEBUG: ClusterIDs_set=", ClusterIDs_set)
        common_ClusterID = min(ClusterIDs_set) # give them the same Cluster ID (by default the smallest)
        for i in range(len(updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list)):
            if updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list[i][5]==Hreson and updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list[i][6]==Creson:
                updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list[i][7] = common_ClusterID    # reinitialize the Cluster ID

    # When 2 peaks are grouped they cannot be predicted to be a methyl (e.g. protein MS6282 L88)
    reupdated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list = []
    ClusterIDs_set = set([p[7] for p in updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list])   # all the Cluster IDs
    for clustID in ClusterIDs_set:
        CH_reson_tuple_set = set([(p[5],p[6]) for p in updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list if p[7]==clustID])    # get all peaks of this Cluster
        print("DEBUG: CH_reson_tuple_set=", CH_reson_tuple_set)
        predictions_list = [p for p in updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list if p[7]==clustID]    # all C-H predictions of this Cluster
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

    print("DEBUG group_carbons: updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list=", updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list)
    print("DEBUG group_carbons: reupdated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list=", reupdated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list)
    return reupdated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list


def is_peak_assigned(prediction_list, prediction):
    """
        FUNCTION to check whether prediction is already in prediction_list.
    """
    Hreson = prediction[5]
    Creson = prediction[6]
    for pred in prediction_list:
        if pred[5] == Hreson and pred[6] == Creson:
            return True
    return False


def select_best_peak_combination(valid_combinations):
    """
        FUNCTION to select the best peak assignment combination according to the score product. Works for combinations of every size >=2.
    """
    valid_combinations = [list(comb) for comb in valid_combinations]
    for comb in valid_combinations:
        comb.append(np.product([c[3] for c in comb]))
    N = len(comb) # number of elements in each comb list
    valid_combinations.sort(key=itemgetter(N-1), reverse=True)    # sort by last elements (score product)
    top_combination = [tuple(c) for c in valid_combinations[0][:-1]]
    return top_combination   # return the top scored peak combination, not their score product


def filter_for_thresholds(possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_list_list, iteration, args):
    """
        RETURNS:
        possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_list_list:    the same list but without the C-groups that do not satisfy the "within" and "between" C-group threshold
    """
    print("DEBUG: iteration=", iteration)
    if iteration == 1:
        WITHIN_CGROUP_THRESHOLD = args.WITHIN_CGROUP_THRESHOLD_ITER1
        BETWEEN_CGROUP_THRESHOLD = args.BETWEEN_CGROUP_THRESHOLD_ITER1
    elif iteration == 2:
        WITHIN_CGROUP_THRESHOLD = args.WITHIN_CGROUP_THRESHOLD_ITER2
        BETWEEN_CGROUP_THRESHOLD = args.BETWEEN_CGROUP_THRESHOLD_ITER2
    elif iteration == 3:
        WITHIN_CGROUP_THRESHOLD = args.WITHIN_CGROUP_THRESHOLD_ITER3
        BETWEEN_CGROUP_THRESHOLD = args.BETWEEN_CGROUP_THRESHOLD_ITER3
    elif iteration == 4:
        WITHIN_CGROUP_THRESHOLD = args.WITHIN_CGROUP_THRESHOLD_ITER4
        BETWEEN_CGROUP_THRESHOLD = args.BETWEEN_CGROUP_THRESHOLD_ITER4
    elif iteration == 5:
        WITHIN_CGROUP_THRESHOLD = args.WITHIN_CGROUP_THRESHOLD_ITER5
        BETWEEN_CGROUP_THRESHOLD = args.BETWEEN_CGROUP_THRESHOLD_ITER5


    # THE FILTERING SHOULD NOT BE APPLIED TO 'TOCSY (unmatched)' OR 'TOCSY-HCNH' PEAKS

    print("DEBUG: before filtering possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_list_list=", possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_list_list)

    # IF A PEAK IS TOCSY-HCNH AND HAS CARBON TYPE c' REMOVE ALL OTHER C-GROUP PREDICTIONS THAT WERE OF TYPE c'
    TOCSY_HCNH_peak_pairs_list = [peak_pair for peak_pair in possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_list_list if 'TOCSY-HCNH' in peak_pair[6]]
    for TOCSY_HCNH_peak_pair in TOCSY_HCNH_peak_pairs_list:
        Cname = TOCSY_HCNH_peak_pair[1]    # this Cname must be fixed to this Cgroup
        Cgroup = TOCSY_HCNH_peak_pair[5]
        for i in range(len(possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_list_list)):
            if possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_list_list[i][1] == Cname and possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_list_list[i][5] != Cgroup:
                possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_list_list[i] = [None, None, None, None, None, None, None]
    cleaned_possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_list_list = [peak_pair for peak_pair in possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_list_list if peak_pair[0] != None]
    Cgroup_set = set([peak_pair[5] for peak_pair in cleaned_possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_list_list])
    equivalentCname_dict = {"CD1": "CD2", "CG1": "CG2", "CD2": "CD1", "CG2": "CG1"}
    equivalentHname_dict = {"HD1": "HD2", "HG1": "HG2", "HD2": "HD1", "HG2": "HG1"}
    # CHECK THE WITHIN C-GROUP THRESHOLD CRITERION
    filtered1_assignments_list = []
    for Cgroup in Cgroup_set:
        print("DEBUG filter1: Cgroup=", Cgroup)
        Cgroup_assignments = [x for x in cleaned_possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_list_list if x[5]==Cgroup]
        Cgroup_assignments.sort(key=itemgetter(3), reverse=True) # sort by probability
        if len(Cgroup_assignments) == 1:
            print("DEBUG filter1: point 1")
            filtered1_assignments_list.append(tuple(Cgroup_assignments[0]))
        elif len(Cgroup_assignments) >= 2 and (     # treat specially the equivalent Carbons
            #(Cgroup_assignments[0][0] == "LEU" and Cgroup_assignments[0][1] in ['CD1', 'CD2'] and Cgroup_assignments[1][1] in ['CD1', 'CD2']) or
            (Cgroup_assignments[0][0]=="VAL" and Cgroup_assignments[0][1] in ['CG1', 'CG2'] and Cgroup_assignments[1][1] in ['CG1', 'CG2']) ):
            if len(Cgroup_assignments) == 2:    # save both equivalent Carbons as possible assignments to the same peak
                print("DEBUG filter1: 2.1")
                filtered1_assignments_list.append(tuple(Cgroup_assignments[0]))
                filtered1_assignments_list.append(tuple(Cgroup_assignments[1]))
            elif Cgroup_assignments[0][3]==Cgroup_assignments[1][3] and Cgroup_assignments[1][3]/Cgroup_assignments[2][3] > WITHIN_CGROUP_THRESHOLD:   # compare the 2nd with the 3rd because the probabilities of the 1st and the 2nd are equal
                print("DEBUG filter1: point 2.2")
                filtered1_assignments_list.append(tuple(Cgroup_assignments[1])) # Cgroup_assignments[0] adds "CD1" and Cgroup_assignments[1] "CD2" for TYR
                filtered1_assignments_list.append(tuple(Cgroup_assignments[0]))

        elif len(Cgroup_assignments) >= 2 and Cgroup_assignments[0][0] == "LEU":    # exceptionally skip WITHIN_CGROUP_THRESHOLD for LEU because it may have 2 methyls (CD1, CD2, CG)
            for a in Cgroup_assignments:    # save all Cgroup assignments
                filtered1_assignments_list.append(tuple(a))
        elif len(Cgroup_assignments) > 1 and Cgroup_assignments[0][3]/Cgroup_assignments[1][3] > WITHIN_CGROUP_THRESHOLD:
            print("DEBUG filter1: point 3")
            filtered1_assignments_list.append(tuple(Cgroup_assignments[0]))
        elif len(Cgroup_assignments) == 0:
            print("ERROR: empty C-group!!")
            sys.exit(1)
        else:
            print("DEBUG: WITHIN C-GROUP THRESHOLD CRITERION VIOLATED BY ", tuple(Cgroup_assignments[0]), "and", tuple(Cgroup_assignments[1]))

    # CHECK THE BETWEEN C-GROUP THRESHOLD CRITERION
    print("DEBUG: filtered1_assignments_list=", filtered1_assignments_list)
    filtered2_assignments_list = []
    Cname_set = set([x[1] for x in filtered1_assignments_list])
    aa_type = cleaned_possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_list_list[0][0]
    skip_LEU = False
    if aa_type == "LEU" and "CG" in Cname_set:
        skip_LEU = True
    for Cname in Cname_set:
        if skip_LEU:
            break
        Cname_assignments = [x for x in cleaned_possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_list_list if x[1]==Cname]
        if Cname_assignments[0][0] == 'ILE' and iteration in [1,3]: # treat ILE methyls specially without BETWEEN_CGROUP_THRESHOLD
            continue
        Cname_assignments.sort(key=itemgetter(3), reverse=True) # sort by probability
        print("DEBUG: Cname_assignments=", Cname_assignments)
        print("DEBUG: filtered2_assignments_list=", filtered2_assignments_list)
        if len(Cname_assignments) == 1:
            print("point 1")
            filtered2_assignments_list.append(tuple(Cname_assignments[0]))
        # If there are 2 C-groups for the same C-type, each one containing one peak, and there are matched or unmatched from TOCSY, save them without applying the filter
        elif len(Cname_assignments) == 2 and Cname_assignments[0][6][0] in ['TOCSY-HCNH','TOCSY (unmatched)'] and Cname_assignments[1][6][0] in ['TOCSY-HCNH','TOCSY (unmatched)']:
            print("point 2")
            filtered2_assignments_list.append(tuple(Cname_assignments[0]))
            filtered2_assignments_list.append(tuple(Cname_assignments[1]))
        elif len(Cname_assignments) >= 2 and (     # treat specially the equivalent Carbons
            (Cname_assignments[0][0] == "LEU" and Cname_assignments[0][1] in ['CD1', 'CD2'] and Cname_assignments[1][1] in ['CD1', 'CD2']) or
            (Cname_assignments[0][0]=="VAL" and Cname_assignments[0][1] in ['CG1', 'CG2'] and Cname_assignments[1][1] in ['CG1', 'CG2']) ):
            print("point 3")
            if len(Cname_assignments) == 2:    # save both equivalent Carbons as possible assignments to the 2 different peaks
                print("point 3.1")
                peaks_matching_Cname = [c[4] for c in filtered2_assignments_list if c[1]==Cname]
                all_filtered2_peaks = [c[4] for c in filtered2_assignments_list]
                if len(peaks_matching_Cname) == 0:  # if there are no peaks with this Cname in filtered2_assignments_list
                    if not Cname_assignments[0][4] in all_filtered2_peaks:    # if the 1st peak has not been saved already as 'Cname', then keep it
                        print("point 3.1.1")
                        filtered2_assignments_list.append(tuple(Cname_assignments[0]))
                    elif not Cname_assignments[1][4] in all_filtered2_peaks:  # otherwise keep the 2nd peak if it has not been saved already as 'Cname'
                        print("point 3.1.2")
                        filtered2_assignments_list.append(tuple(Cname_assignments[1]))
                #elif len(peaks_matching_Cname) > 0: # TODO
            elif len(Cname_assignments) > 2 and Cname_assignments[0][3] == Cname_assignments[1][3] and Cname_assignments[1][3]/Cname_assignments[2][3] > BETWEEN_CGROUP_THRESHOLD:   # compare the 2nd with the 3rd because the probabilities of the 1st and the 2nd are equal
                print("point 3.2")
                peaks_matching_Cname = [c[4] for c in filtered2_assignments_list if c[1]==Cname]
                all_filtered2_peaks = [c[4] for c in filtered2_assignments_list]
                if len(peaks_matching_Cname) == 0:  # if there are no peaks with this Cname in filtered2_assignments_list
                    if not Cname_assignments[0][4] in all_filtered2_peaks:    # if the 1st peak has not been saved already as 'Cname', then keep it
                        print("point 3.2.1")
                        filtered2_assignments_list.append(tuple(Cname_assignments[0]))
                    elif not Cname_assignments[1][4] in all_filtered2_peaks:  # otherwise keep the 2nd peak if it has not been saved already as 'Cname'
                        print("point 3.2.2")
                        filtered2_assignments_list.append(tuple(Cname_assignments[1]))
            elif len(Cname_assignments) > 2 and Cname_assignments[0][3] != Cname_assignments[1][3]:   # compare the 2nd with the 3rd because the probabilities of the 1st and the 2nd are not equal
                print("point 3.3")
                peaks_matching_Cname = [c[4] for c in filtered2_assignments_list if c[1]==Cname]
                peaks_matching_equivalentCname = [c[4] for c in filtered2_assignments_list if c[1]==equivalentCname_dict[Cname]]
                all_filtered2_peaks = [c[4] for c in filtered2_assignments_list]
                if len(peaks_matching_Cname) == 0 and len(peaks_matching_equivalentCname) == 0:  # if there are no peaks with this Cname or its equivalent Cname in filtered2_assignments_list
                    if not Cname_assignments[0][4] in all_filtered2_peaks and Cname_assignments[0][3]/Cname_assignments[2][3] > BETWEEN_CGROUP_THRESHOLD:    # if the 1st peak has not been saved already as 'Cname', and its score is more than BETWEEN_CGROUP_THRESHOLD times greater than the score of the 3rd peak, then keep it
                        print("point 3.2.1")
                        filtered2_assignments_list.append(tuple(Cname_assignments[0]))
                    if not Cname_assignments[1][4] in all_filtered2_peaks and Cname_assignments[1][3]/Cname_assignments[2][3] > BETWEEN_CGROUP_THRESHOLD:    # if the 2nd peak has not been saved already as 'Cname', and its score is more than BETWEEN_CGROUP_THRESHOLD times greater than the score of the 3rd peak, then keep it
                        print("point 3.2.2")
                        # if both the 1st and 2nd peak have the same equivalent carbon name (e.g. CD1 and CD1), change the Cname of the second before you save it (e.g. CD1->CD2)
                        if Cname_assignments[0][1] == Cname_assignments[1][1]:
                            # e.g. Cname_assignments[1] = ['LEU', 'CD1', ('HD1',), 22.75368059708115, ('0.62_24.034',), 2, ('HCNH',)]
                            aa, C, H_t, score, peak_t, Cgrp, spectrum_t = Cname_assignments[1]
                            C = equivalentCname_dict[C]
                            H_t = (equivalentHname_dict[H_t[0]],)
                            filtered2_assignments_list.append((aa, C, H_t, score, peak_t, Cgrp, spectrum_t))
                        else:
                            filtered2_assignments_list.append(tuple(Cname_assignments[1]))
                #elif len(peaks_matching_Cname) > 0: # TODO
        elif len(Cname_assignments) > 1 and Cname_assignments[0][3]/Cname_assignments[1][3] > BETWEEN_CGROUP_THRESHOLD:
            print("point 4")
            filtered2_assignments_list.append(tuple(Cname_assignments[0]))
        elif len(Cname_assignments) == 0:
            print("ERROR: empty C-group!!")
            sys.exit(1)
        else:
            print("DEBUG: BETWEEN C-GROUP THRESHOLD CRITERION VIOLATED BY ", tuple(Cname_assignments[0]), "and", tuple(Cname_assignments[1]))

    # TREAT ILE METHYLS SPECIALLY
    aa_type = cleaned_possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_list_list[0][0]
    if aa_type == "ILE" and iteration in [1,3]:
        print("point 5")
        filtered1_assignments_list.sort(key=itemgetter(3), reverse=True)    # highest score first
        Cname_set = set([x[1] for x in filtered1_assignments_list])
        if "CG2" in Cname_set and "CD1" in Cname_set:
            for Cname in ["CG2", "CD1"]:
                for peak in filtered1_assignments_list:
                    if peak[1] == Cname and not peak[4][0] in [p[4][0] for p in filtered2_assignments_list]:
                        filtered2_assignments_list.append(tuple(peak))
                        break
        elif "CG2" in Cname_set and not "CD1" in Cname_set:
            CG2_assignments = [x for x in cleaned_possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_list_list if x[1]=="CG2"]
            CG2_assignments.sort(key=itemgetter(3), reverse=True) # sort by probability
            if len(CG2_assignments) == 1:
                print("point 5.2.1")
                filtered2_assignments_list.append(tuple(CG2_assignments[0]))
            elif len(CG2_assignments) > 1 and CG2_assignments[0][3]/CG2_assignments[1][3] > BETWEEN_CGROUP_THRESHOLD:
                print("point 5.2.2")
                filtered2_assignments_list.append(tuple(CG2_assignments[0]))
        elif not "CG2" in Cname_set and "CD1" in Cname_set:
            CD1_assignments = [x for x in cleaned_possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_list_list if x[1]=="CD1"]
            CD1_assignments.sort(key=itemgetter(3), reverse=True) # sort by probability
            if len(CD1_assignments) == 1:
                print("point 5.3.1")
                filtered2_assignments_list.append(tuple(CD1_assignments[0]))
            elif len(CD1_assignments) > 1 and CD1_assignments[0][3]/CD1_assignments[1][3] > BETWEEN_CGROUP_THRESHOLD:
                print("point 5.3.2")
                filtered2_assignments_list.append(tuple(CD1_assignments[0]))

    # If CG is also in the possible Ctype assignments of LEU, then examine all permutation and select the best by the probability product.
    # Recover the original missing_carbons_list
    missing_carbons_list = list(set([a[1] for a in possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_list_list]))
    if aa_type == "LEU" and skip_LEU == True:
        print("point 6")
        filtered1_assignments_list.sort(key=itemgetter(3), reverse=True)    # highest score first
        Cname_set = set([x[1] for x in filtered1_assignments_list])
        Cgroup_set = set([x[5] for x in filtered1_assignments_list])
        if len(filtered1_assignments_list) == 1:
            print("point 6.1")
            filtered2_assignments_list.append(tuple(filtered1_assignments_list[0]))
        elif len(filtered1_assignments_list) == 2:
            print("point 6.2")
            if len(Cgroup_set) == 1: # keep the most probable assignment for only Cgroup
                print("point 6.2.1")
                filtered2_assignments_list.append(tuple(filtered1_assignments_list[0]))
            elif len(Cgroup_set) == 2 and len(Cname_set) == 1:  # if only one available Cname, then save the most probable Cgroup
                print("point 6.2.2")
                filtered2_assignments_list.append(tuple(filtered1_assignments_list[0]))
            elif len(Cgroup_set) == 2 and len(Cname_set) == 2:  # save them both
                print("point 6.2.3")
                filtered2_assignments_list.append(tuple(filtered1_assignments_list[0]))
                filtered2_assignments_list.append(tuple(filtered1_assignments_list[1]))
        elif len(filtered1_assignments_list) >= 3:
            if len(Cgroup_set) == 1: # keep the most probable assignment for only Cgroup
                print("point 6.3.1")
                filtered2_assignments_list.append(tuple(filtered1_assignments_list[0]))
            elif len(Cgroup_set) >= 2 and len(Cname_set) == 1:  # if only one available Cname, then save the most probable Cgroup
                print("point 6.3.2")
                filtered2_assignments_list.append(tuple(filtered1_assignments_list[0]))
            elif len(Cgroup_set) >= 3 and len(Cname_set) == 3:
                print("point 6.4")
                # check all ternary combinations of Cgroups. If no combination contains all 3 Cnames, then...
                found_ternary_combination = False
                found_binary_combination = False
                valid_triplets = []
                for comb in combinations(filtered1_assignments_list, 3):
                    Cgroup_list = [c[5] for c in comb]
                    Cgroup_set = set(Cgroup_list)
                    Cnames_list = [c[1] for c in comb]
                    if len(Cgroup_list) < 3 or len(Cgroup_set) < 3:    # we want only combinations of 3 different Cgroups
                        continue
                    if "CG" in Cnames_list and "CD1" in Cnames_list and "CD2" in Cnames_list:
                        print("point 6.4.1 appending triplet=", [comb[0], comb[1], comb[2]])
                        valid_triplets.append([comb[0], comb[1], comb[2]])
                        found_ternary_combination = True
                    elif "CG" in Cnames_list and ( Cnames_list.count("CD1") == 2 or Cnames_list.count("CD2") == 2 ) and "CD1" in missing_carbons_list and "CD2" in missing_carbons_list:
                    # rename the first Cgroup to its equivalent carbon (CD1->CD2, or CD2->CD1) in order to have both carbons, provided that both CD1,CD2 are missing
                        triplet = []
                        if "CD1" in Cnames_list:
                            print("point 6.4.2.1")
                            indx = Cnames_list.index("CD1")
                        elif "CD2" in Cnames_list:
                            print("point 6.4.2.2")
                            indx = Cnames_list.index("CD2")
                        for i in [0,1,2]:   # save the assignments
                            print("point 6.4.2.3")
                            if i == indx:
                                C = comb[indx][1]
                                H = comb[indx][2][0]
                                triplet.append( (comb[indx][0], equivalentCname_dict[C], (equivalentHname_dict[H],), comb[indx][3], comb[indx][4], comb[indx][5], comb[indx][6]) )
                            else:
                                triplet.append(comb[i])
                        print("point 6.4.2.4 appending triplet=", triplet)
                        valid_triplets.append(triplet)
                        found_ternary_combination = True

                # if found, save the valid duplet with the highest score product
                if found_ternary_combination == True:
                    print("point 6.4.3")
                    print("DEBUG: valid_triplets=", valid_triplets)
                    best_ternary_combination = select_best_peak_combination(valid_triplets)
                    filtered2_assignments_list.extend(best_ternary_combination)
                # ... check also all binary combinations of Cgroups. If no combination contains 2 Cnames, then ...
                # elif found_ternary_combination == False:
                print("point 6.4.4")
                valid_duplets = []
                for comb in combinations(filtered1_assignments_list, 2):
                    Cgroup_list = [c[5] for c in comb]
                    Cgroup_set = set(Cgroup_list)
                    Cnames_list = [c[1] for c in comb]
                    if len(Cgroup_list) < 2 or len(Cgroup_set) < 2:    # we want only combinations of 2 different Cgroups
                        continue
                    if "CG" in Cnames_list and ("CD1" in Cnames_list or "CD2" in Cnames_list):
                        print("point 6.4.4.1 appending duplet=", [comb[0], comb[1]])
                        valid_duplets.append([comb[0], comb[1]])
                        found_binary_combination = True
                    elif "CD1" in Cnames_list and "CD2" in Cnames_list:
                        print("point 6.4.4.2 appending duplet=", [comb[0], comb[1]])
                        valid_duplets.append([comb[0], comb[1]])
                        found_binary_combination = True
                    elif (Cnames_list == ["CD1", "CD1"] or Cnames_list == ["CD2", "CD2"]) and "CD1" in missing_carbons_list and "CD2" in missing_carbons_list:
                        # rename the first Cgroup to its equivalent carbon (e.g. CD1->CD2) in order to have both carbons, provided that both CD1, CD2 are missing
                        C = comb[0][1]
                        H = comb[0][2][0]
                        print("point 6.4.4.3 appending duplet=", ( [comb[0][0], equivalentCname_dict[C], (equivalentHname_dict[H],), comb[0][3], comb[0][4], comb[0][5], comb[0][6]], comb[1] ))
                        valid_duplets.append( ( [comb[0][0], equivalentCname_dict[C], (equivalentHname_dict[H],), comb[0][3], comb[0][4], comb[0][5], comb[0][6]], comb[1] ) )
                        found_binary_combination = True
                    else:
                        print("point 6.4.4.5")
                        continue
                # if found, save the valid duplet with the highest score product
                if found_binary_combination == True:
                    print("point 6.4.5")
                    print("DEBUG: valid_duplets=", valid_duplets)
                    best_binary_combination = select_best_peak_combination(valid_duplets)
                    filtered2_assignments_list.extend(best_binary_combination)
                # ... keep also the top scored singlet
                # elif found_binary_combination == False and found_ternary_combination == False:
                print("point 6.4.6")
                filtered2_assignments_list.append(tuple(filtered1_assignments_list[0]))
            elif len(Cgroup_set) >= 2 and len(Cname_set) >= 2:
                print("point 6.3.3")
                # check all binary combinations of Cgroups. If no combination contains 2 different Cgroups and 2 Cnames, then ...
                found_binary_combination = False
                valid_duplets = []
                for comb in combinations(filtered1_assignments_list, 2):
                    Cgroup_list = [c[5] for c in comb]
                    Cgroup_set = set(Cgroup_list)
                    Cnames_list = [c[1] for c in comb]
                    if len(Cgroup_list) < 2 or len(Cgroup_set) < 2:    # we want only combinations of 2 different Cgroups
                        continue
                    if "CG" in Cnames_list and ("CD1" in Cnames_list or "CD2" in Cnames_list):
                        print("point 6.3.3.1 appending duplet=", [comb[0], comb[1]])
                        valid_duplets.append([comb[0], comb[1]])
                        found_binary_combination = True
                    elif "CD1" in Cnames_list and "CD2" in Cnames_list:
                        print("point 6.3.3.2 appending duplet=", [comb[0], comb[1]])
                        valid_duplets.append([comb[0], comb[1]])
                        found_binary_combination = True
                    elif (Cnames_list == ["CD1", "CD1"] or Cnames_list == ["CD2", "CD2"]) and "CD1" in missing_carbons_list and "CD2" in missing_carbons_list:
                        # rename the first Cgroup to its equivalent carbon (e.g. CD1->CD2) in order to have both carbons, provided that both CD1,CD2 are missing
                        C = comb[0][1]
                        H = comb[0][2][0]
                        print("point 6.3.3.3 appending duplet=", ( [comb[0][0], equivalentCname_dict[C], (equivalentHname_dict[H],), comb[0][3], comb[0][4], comb[0][5], comb[0][6]], comb[1] ))
                        valid_duplets.append( ( [comb[0][0], equivalentCname_dict[C], (equivalentHname_dict[H],), comb[0][3], comb[0][4], comb[0][5], comb[0][6]], comb[1] ) )
                        found_binary_combination = True
                    else:
                        print("point 6.3.3.4")
                        continue
                # if found, save the valid duplet with the highest score product
                if found_binary_combination == True:
                    print("point 6.3.4")
                    print("DEBUG: valid_duplets=", valid_duplets)
                    best_binary_combination = select_best_peak_combination(valid_duplets)
                    filtered2_assignments_list.extend(best_binary_combination)
                # ...keep only the top scored singlet
                # elif found_binary_combination == False:
                print("point 6.3.5")
                filtered2_assignments_list.append(tuple(filtered1_assignments_list[0]))


    print("DEBUG: final filtered1_assignments_list=", filtered1_assignments_list)
    print("DEBUG: final filtered2_assignments_list=", filtered2_assignments_list)
    filtered1_assignments_set = set(filtered1_assignments_list)
    filtered2_assignments_set = set(filtered2_assignments_list)
    filtered_possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_list_list = list(filtered1_assignments_set.intersection(filtered2_assignments_set))
    filtered_possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_list_list.sort(key=itemgetter(3), reverse=True) # sort by probability

    print("DEBUG: after filtering filtered_possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_list_list=", filtered_possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_list_list)
    return filtered_possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_list_list


def get_best_prediction(correct_chainScore_list, iteration=None):
    """
        FUNCTION correct_chainScore_list is sorted to the last element. This method treats specially LEU that contains 3 carbons
        with only one visible proton and a methylene that may contain one or two peaks. In that case not necessary the combination
        with the most peaks wins.
    """

    if len(correct_chainScore_list) == 0:
            return []

    aa_type = correct_chainScore_list[0][0][0]
    if iteration == 1 and aa_type == "LEU":
        # In LEU we have 3 methyls, therefore keep the triplet or duplet with the highest score
        # first find the best single peak
        singlets = [comb for comb in correct_chainScore_list if len(comb)==2]
        singlets.sort(key=itemgetter(1), reverse=True)
        # then find the duplet with the highest score
        duplets = [comb for comb in correct_chainScore_list if len(comb)==3]
        duplets.sort(key=itemgetter(2), reverse=True)
        # then find the triplet with the highest score
        triplets = [comb for comb in correct_chainScore_list if len(comb)==4]
        triplets.sort(key=itemgetter(3), reverse=True)

        if len(triplets) > 0 and len(duplets) > 0:
            if triplets[0][-1] > duplets[0][-1]:    # if the score of the triplet is higher than the score of the duplet
                return triplets[0][:-1]
            elif triplets[0][-1] < duplets[0][-1]:
                return duplets[0][:-1]
        elif len(triplets) == 0 and len(duplets) > 0:
            # return duplets[0][:-1]  # TEMPORARILY IGNORE SINGLES
            if duplets[0][-1] > singlets[0][-1]:
                return duplets[0][:-1]
            elif duplets[0][-1] < singlets[0][-1]:
                return singlets[0][:-1]
        elif len(triplets) > 0 and len(duplets) == 0:
            return triplets[0][:-1]
        elif len(triplets) == 0 and len(duplets) == 0:
            return singlets[0][:-1]

    elif iteration == 1:
        # recover the missing_carbons_list
        missing_carbons_list = list(set([p[1] for chain in correct_chainScore_list for p in chain[:-1]]))
        index = 0
        best_prediction_list = []
        assigned_carbons_list = []
        assigned_peaks_list = []
        skip_chain = False
        while len(best_prediction_list) < len(missing_carbons_list) and index < len(correct_chainScore_list):    # in iter1 we work with methyls, therefore each peak corresponds to 1 carbon!
            for peak in correct_chainScore_list[index][:-1]:
                if peak[1] in assigned_carbons_list or peak[4] in assigned_peaks_list:    # rule out chain that contained already assigned carbons
                    skip_chain = True
                    break
            if skip_chain == False:
                best_prediction_list.extend(correct_chainScore_list[index][:-1])
            assigned_carbons_list = [p[1] for p in best_prediction_list]    # update the assigned_carbons_list
            assigned_peaks_list = [p[4] for p in best_prediction_list]    # update the assigned_peaks_list
            index += 1  # move to the next chain
            skip_chain = False

        return best_prediction_list
    else:
        return correct_chainScore_list[0][:-1]  # list with the assignment combination that has the highest probability


def select_correct_H_C_resonpairs3(aatypeResonpairMatchesTuple_list, iter1_allowed_atoms_dict, iteration, args):
    """
        NEW FUNCTION to match the C-groups to as much as possible C-types of a PARTICULAR aatype. This function selects the
        best C-group assignment combination based on the product of probabilities of the individual C-group assignments.
        *** The difference here is that we work with C-groups, not with individual peaks! ***

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

    print("DEBUG select_correct_H_C_resonpairs3: aatypeResonpairMatchesTuple_list=", aatypeResonpairMatchesTuple_list)
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
            matching_assignments = list(set([tuple(x) for x in filtered_aatypeResonpairMatchesTuple_list if x[1]==Cname and x[7]==Cgroup and x[3]>0]))   # exclude Carbon based only assignments (they have negative probability)
            if len(matching_assignments) == 0:  # if there are not assignments of this C-type for this C-group, continue to the next C-type
                continue
            Hlist_of_this_Cgroup = list(set([x[2] for x in matching_assignments]))  # it doesn't matter if the order will be changed, both protons are equivalent
            Hlist_of_this_Cgroup.sort()
            HCresonpair_spectrumType_set_of_this_Cgroup = tuple(set([(str(x[5])+"_"+str(x[6]), x[8]) for x in matching_assignments]))
            Cgroup_probability = get_Cgroup_consensus_probability([x[3] for x in matching_assignments
                                                                   if x[1]==Cname and x[7]==Cgroup],
                                                                  cons_mode=args.CONSENSUS_CGROUP_PROB_MODE)
            if len(HCresonpair_spectrumType_set_of_this_Cgroup) == 1:    # if the C-group contains only one peak
                Hlist_of_this_Cgroup = [Hlist_of_this_Cgroup[0]]    # keep only one of the ethyl proton names
            HCresonpair_set_of_this_Cgroup = [t[0] for t in HCresonpair_spectrumType_set_of_this_Cgroup]
            HCresonpair_set_of_this_Cgroup.sort()   # always the lowest proton resonance goes 1st
            HCresonpair_set_of_this_Cgroup = tuple(HCresonpair_set_of_this_Cgroup)  # convert the list to tuple
            spectrum_combos = tuple([t[1] for t in HCresonpair_spectrumType_set_of_this_Cgroup])
            Cgroup_assignment = [aa_type, Cname, tuple(Hlist_of_this_Cgroup), Cgroup_probability, HCresonpair_set_of_this_Cgroup, Cgroup, spectrum_combos]
            possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_list_list.append(Cgroup_assignment)
    print("DEBUG: possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_list_list=", possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_list_list)
    print("DEBUG: Cgroups_set=", Cgroups_set)
    if (iteration==1 and args.WITHIN_CGROUP_THRESHOLD_ITER1==0 and args.BETWEEN_CGROUP_THRESHOLD_ITER1==0) or (
        iteration==2 and args.WITHIN_CGROUP_THRESHOLD_ITER2==0 and args.BETWEEN_CGROUP_THRESHOLD_ITER2==0) or (
        iteration==3 and args.WITHIN_CGROUP_THRESHOLD_ITER3==0 and args.BETWEEN_CGROUP_THRESHOLD_ITER3==0) or (
        iteration==4 and args.WITHIN_CGROUP_THRESHOLD_ITER4==0 and args.BETWEEN_CGROUP_THRESHOLD_ITER4==0) or (
        iteration==5 and args.WITHIN_CGROUP_THRESHOLD_ITER5==0 and args.BETWEEN_CGROUP_THRESHOLD_ITER5==0):
        print("DEBUG: no filtering for WITHIN_CGROUP_THRESHOLD and BETWEEN_CGROUP_THRESHOLD")
        filtered_possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_list_list = possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_list_list
    else:
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
        # To make sure that all the Cgroups and possible Cnames are used, for each Cluster combination build a tree each time starting
        # from a different Cgroup. That will produce redundant chains but will not miss any Cgroup combination.
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
                #print "DEBUG: leaf.name=", leaf.name
                for ancestor in leaf.get_ancestors():
                    # NOTE: instead of 'if ancestor == Root: break' , I check if it's name is an int (only Root has name string).
                    #print "DEBUG: ancestor.name=", ancestor.name
                    try:
                        assignID = int(ancestor.name)
                    except ValueError:
                        continue
                    chain.append(possible_aatype_C_Hset_aveprob_HCresonapairset_Cgroup_assignID_list_list[assignID])
                    prob_product *= ancestor.Prob
                if iteration > 1:
                    if len(chain) < filtered_Cnames_num: # DANGEROUS: discard combinations with less than Cgroups_num peaks: THIS DOES NOT APPLY IN
                                                    # HCNH i->i+1 ASSIGNMENT BECAUSE WE ARE NEVER SURE WHETHER ALL THE i->i+1 CGROUPS CORRESPOND TO RESIDUE i
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
        assigned_Ctypes_list = []   # list with the Ctypes encountered already
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

    correct_chainScore_list.sort(key=itemgetter(-1), reverse=True) # sort probability by descending order; the highest probability combination must be 1st
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

    return list(best_prediction_list)


def populate_leaves_for_correct_H_C_resonpairs3(Assignment_Tree,
                                                possible_aatype_prob_C_H_resonpair_TOCSYindex_assignID_list_list,
                                                Cluster):
    """
        FUNCTION that adds new branches to the leaves of the Tree.
        ARGUMENTS:
        Assignment_Tree:    The Tree structure with connectivities
        RETURNS:
        (Assignment_Tree, BOOLEAN):    A tuple with elements the input Tree structure with new branches (if applicable), and a BOOLEAN value which is True if the function added
                                       new leaves to the Tree, or False otherwise
    """
    print("DEBUG: possible_aatype_prob_C_H_resonpair_TOCSYindex_assignID_list_list=",
          possible_aatype_prob_C_H_resonpair_TOCSYindex_assignID_list_list, "Cluster=", Cluster)

    nucleusAssignmentsforCurrentCluster_list = [x for x in possible_aatype_prob_C_H_resonpair_TOCSYindex_assignID_list_list if x[5]==Cluster] # save all the members of this Cluster ID
    # e.g. [['ARG', 'CB', ('HB2',), 1.79677020689e-30, ('2.648_41.314',), 4, 11], ['ARG', 'CD', ('HD2',), 7.50809746006e-05, ('2.648_41.314',), 4, 12]]
    assignments_list = []
    #if len(NresonIndex_set) == 1:
    assignments_list = nucleusAssignmentsforCurrentCluster_list
    print("DEBUG: assignments_list=", assignments_list)
    number_of_new_leaves = 0
    # ATTENTION: never use Assignment_Tree.iter_leaf_names(), it doesn't return the names in the order
    # ATTENTION: corresponding to Assignment_Tree.get_leaves()!!!
    for leaf in Assignment_Tree.get_leaves():
        print("DEBUG: leaf=", leaf.get_ascii(), "leaf.name=", leaf.name)
        try:
            for assignment in assignments_list:
                _Ctype = assignment[1]
                _Htypes = assignment[2]
                _Prob = assignment[3]
                if _Prob == None:
                    raise Exception(Debuginfo("FAIL: the following assignment is wrong, it has None weighted probability!"
                                              "\n %s" % assignment, fail=True))
                _Resonances = assignment[4]
                _clustID = assignment[5]
                assignID = assignment[6]
                ancestors_Ctype_list = [ancestor.Ctype for ancestor in leaf.get_ancestors()]    # list with all the C-types currently in the branch starting from this leaf
                ancestors_Ctype_list.append(leaf.Ctype)
                print("DEBUG: ancestors_Ctype_list=", ancestors_Ctype_list)
                if _Ctype in ancestors_Ctype_list:   # if this C-type is already in the Tree, skip it
                    print("DEBUG: _Ctype", _Ctype, "in ancestors_Ctype_list", ancestors_Ctype_list)
                    continue
                new_child = leaf.add_child(name=assignID) # add a new brach to the current TOCSY add index (leaf) with length the respective probability and return it
                print("DEBUG: adding new leaf: ", assignment)
                new_child.add_features(Ctype=_Ctype, Htypes=_Htypes, Prob=_Prob, clustID=_clustID)
                number_of_new_leaves += 1
        except KeyError:
            continue

    print(Assignment_Tree.get_ascii(show_internal=True, compact=False, attributes=["name", "Prob", "clustID"]))
    # print Assignment_Tree.get_ascii(show_internal=True, compact=False)
    if number_of_new_leaves > 0:
        return (Assignment_Tree, True)
    else:
        return (Assignment_Tree, False)


def convert_Cgroup_to_peak_format(best_prediction_list_CgroupFormat, updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list):
    """
    FUNCTION to convert this format:
    [['THR', 'CG2', ('HG2',), 0.000268993386994, ('0.728_22.147',), 6, 0], ...]
    to this format:

    ARGUMENTS:
    updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list:    the original assignment list in peak format to read the correct Nreson_indices
    """

    best_prediction_list_PeakFormat = []
    print("DEBUG: best_prediction_list_CgroupFormat=", best_prediction_list_CgroupFormat)
    for Cgroup_assignment in best_prediction_list_CgroupFormat:
        Hnames_tuple = Cgroup_assignment[2]
        CH_reson_tuple = Cgroup_assignment[4]
        spectrum_combo_tuple = Cgroup_assignment[7]
        for i in range(len(Hnames_tuple)):
            Hname = Hnames_tuple[i]
            Creson = float(CH_reson_tuple[i].split("_")[1])
            Hreson = float(CH_reson_tuple[i].split("_")[0])
            Cname = Cgroup_assignment[1]
            aa_type = Cgroup_assignment[0]
            probability = Cgroup_assignment[3]
            Cgroup = Cgroup_assignment[5]
            spectrum_combo = spectrum_combo_tuple[i]
            try:
                Nreson_index = [x[4] for x in updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list if approx_equal(Hreson, x[5], 0.00001) and approx_equal(Creson, x[6], 0.00001)][0]
            except IndexError:
                print("ERROR:")
                print("DEBUG: updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list=", updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list)
                print("DEBUG: Cname=", Cname, "Hreson=", Hreson, "Creson=", Creson, "Cgroup=", Cgroup)
                sys.exit(1)
            print("DEBUG: converting Cgroup to peak format:", [aa_type, Cname, Hname, probability, Nreson_index, Hreson, Creson, Cgroup, spectrum_combo])
            best_prediction_list_PeakFormat.append([aa_type, Cname, Hname, probability, Nreson_index, Hreson, Creson, Cgroup, spectrum_combo])

    return best_prediction_list_PeakFormat


def get_aatypes_from_all_H_C_resonpairs(all_possible_assignments_list,
                                        i_residue,
                                        partly_assigned_carbons_list,
                                        iteration,
                                        Probability_CS,
                                        HCNH_residue_NHresonances_dict,
                                        matched_HCNH_residue_assignments_dict,
                                        protein_alignment_list,
                                        absolute_matches_alignment_list,
                                        iter1_allowed_atoms_dict,
                                        args):
    """
        FUNCTION to return the matched amino acid types for a particular TOCSY AAIG by using all its H,C resonance pairs.

        ARGUMENTS:
        selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list:   list of lists of the form (aa type, Carbon name, Hydrogen name, weighted average probability, HCNH_reson_index, H resonance, C resonance, carbon group)
        iteration:      used only in select_correct_H_C_resonpairs3, to eventually control the thresholds

        RETURN:
        aatype_CHnucleiType_presenceProbSum_mdict:  multidimensional dict: aa type --> C,H nuclei pair type --> probability of presence or absence (eq 6 in RESCUE2 paper)
        aatype_CHnucleiType_NresonIndex_mdict:      multidimensional dict: aa type --> C,H nuclei pair type --> TOCSY_reson_index

    """

    aa_type = aa1to3_dict[get_aa_type_from_residue(i_residue)]
    partly_matched_HCNH_peaks_list = []
    for Cname in partly_assigned_carbons_list:
        partly_matched_HCNH_peaks_list.extend([p for p in matched_HCNH_residue_assignments_dict[i_residue] if p[0]==i_residue and p[1]==Cname])
    Creson_Hreson_tuple_list = [(p[5],p[4]) for p in all_possible_assignments_list]     # save the C & H resonances of all the peaks in all_possible_assignments_list to avoid double saving them with different proton and spectrum names
    for peak in partly_matched_HCNH_peaks_list:    # save the partly assigned TOCSY-HCNH peaks to all the assignemnts
        Cname = peak[1]
        Hname = aatype_carbon_nondegenerateHlist_mdict[aa_type][Cname][0]  # by default the one found in TOCSY methylene proton is named HB2 or HG2, etc.
        if not (peak[2], peak[4]) in Creson_Hreson_tuple_list:  # if this peaks in not already in all_possible_assignments_list, save it
            all_possible_assignments_list.append([aa_type, peak[1], Hname, 10000000.0, peak[4], peak[2], 'TOCSY-HCNH'])
            Creson_Hreson_tuple_list.append(((peak[2], peak[4])))   # update the list
    # ASSIGN HCNH RESONANCE INDICES TO THE PEAKS
    possible_aatype_prob_C_H_resonpair_HCNHindex_list_list = []
    previous_Creson, previous_Hreson = None, None
    NresonIndex = 0
    HCNHindex_ResonancesTuple_dict = {}
    HCNHindex_spectrumType_dict = {}
    all_possible_assignments_list = sorted(all_possible_assignments_list, key=itemgetter(4,5), reverse=True)
    print("DEBUG: get_aatypes_from_all_H_C_resonpairs: all_possible_assignments_list=", all_possible_assignments_list)
    for assignment in all_possible_assignments_list:
        Creson = assignment[5]
        Hreson = assignment[4]
        spectrum_combo = assignment[6]
        if previous_Hreson != Hreson or previous_Creson != Creson:
            NresonIndex += 1
        print("DEBUG: assignment=", assignment)
        possible_aatype_prob_C_H_resonpair_HCNHindex_list_list.append([assignment[0], assignment[1], assignment[2], assignment[3], NresonIndex, assignment[4], assignment[5], None, assignment[6]])
        try:
            HCNHindex_ResonancesTuple_dict[NresonIndex] = (Hreson, Creson, HCNH_residue_NHresonances_dict[i_residue][0], HCNH_residue_NHresonances_dict[i_residue][1])
        except KeyError:    # this residue has no HCNH, namely no N-H, save only the C,H resonaces
            HCNHindex_ResonancesTuple_dict[NresonIndex] = (Hreson, Creson, None, None)
        HCNHindex_spectrumType_dict[NresonIndex] = spectrum_combo
        previous_Creson = Creson
        previous_Hreson = Hreson

    # GROUP THE CARBONS
    updated_possible_aatype_prob_C_H_resonpair_HCNHindex_list_list = group_carbons(possible_aatype_prob_C_H_resonpair_HCNHindex_list_list,
                                                                                    iteration)
    print("DEBUG: updated_possible_aatype_prob_C_H_resonpair_HCNHindex_list_list=", updated_possible_aatype_prob_C_H_resonpair_HCNHindex_list_list)
    if len(updated_possible_aatype_prob_C_H_resonpair_HCNHindex_list_list) == 0:
        return {}, {}, {}, {}
    # GET THE BEST ASSIGNMENTS FOR THE MISSING CARBONS
    # each element of best_prediction_list_CgroupFormat has the format: ['THR', 'CG2', ('HG2',), 0.000268993386994, ('0.728_22.147',), 6, 0]
    best_prediction_list_CgroupFormat = select_correct_H_C_resonpairs3(updated_possible_aatype_prob_C_H_resonpair_HCNHindex_list_list, iter1_allowed_atoms_dict, iteration, args)
    print("DEBUG: selected_aatype_prob_C_H_resonpair_HCNHindex_tuple_list=", best_prediction_list_CgroupFormat)

    # Now convert the compact C-group based format of best_prediction_list_CgroupFormat to the old peak-based format
    best_prediction_list_PeakFormat = convert_Cgroup_to_peak_format(best_prediction_list_CgroupFormat, updated_possible_aatype_prob_C_H_resonpair_HCNHindex_list_list)
    # The following sorting ensures that the alternative C-H pairs are sorted according to their weighted average probability
    best_prediction_list_PeakFormat.sort(key=itemgetter(0,1,2,3), reverse=True)  # sort by every column except the last 4
    print("DEBUG: best_prediction_list_PeakFormat=",best_prediction_list_PeakFormat)
    aatype_occupancy_dict = {}  # dict with key the aa types and values how many times they were found during the H,C resonance analysis
    aatype_probsum_dict = {}    # dictionary with the aa type --> sum of probabilities of individual matched H,C resonance pairs
    aatype_CHnucleiType_presenceProbSum_mdict = tree()    # multidimensional dict: aa type --> C,H nuclei pair type --> probability of presence or absence (eq 6 in RESCUE2 paper)
    aatype_CHnucleiType_NresonIndex_mdict = tree()    # multidimensional dict: aa type --> C,H nuclei pair type --> HCNH_reson_index
    aa_type_set = set()
    aatype_validC_H_pairsList_dict = {}
    #print "best_prediction_list_PeakFormat=",best_prediction_list_PeakFormat
    #sys.exit(1)

    if len(best_prediction_list_PeakFormat) == 0:  # if no atom assignment could be made, skip this aa! (TEMPORARY FIX)
        return {}, {}, {}, {}
    print("DEBUG: best_prediction_list_PeakFormat=", best_prediction_list_PeakFormat)
    for group8 in best_prediction_list_PeakFormat:
        print("DEBUG: group8=",group8)
        aa_type = group8[0]
        C_name = group8[1]
        H_name = group8[2]
        weighted_average_probability = group8[3]
        HCNH_reson_index = group8[4]
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
            if H_name[0] == 'Q':     # TEMPORARY FIX UNTIL YOU AVERAGE THE DEGENERATE PROTON PROBABILITIES ON Probability_CS dict
                H_name = aa_carbonBondedHydrogensDict_dict[aa_type][C_name][0]
            aatype_CHnucleiType_presenceProbSum_mdict[aa_type][C_name+"_"+H_name] = (args.H_weight * Probability_CS[aa_type][H_name] + args.C_weight * Probability_CS[aa_type][C_name])/float((args.H_weight + args.C_weight))  # 2nd WAY
            aatype_CHnucleiType_NresonIndex_mdict[aa_type][C_name+"_"+H_name] = HCNH_reson_index
        except KeyError:
            aatype_validC_H_pairsList_dict[aa_type] = [C_name+"_"+H_name]
            #if args.ALLOW_SINGLE_CH_PAIRS:
            aa_type_set.add(aa_type)
            aatype_occupancy_dict[aa_type] = 1
            aatype_probsum_dict[aa_type] = weighted_average_probability     # COMMON IN 1st & 2nd WAY
            if H_name[0] == 'Q':     # TEMPORARY FIX UNTIL YOU AVERAGE THE DEGENERATE PROTON PROBABILITIES ON Probability_CS dict
                H_name = aa_carbonBondedHydrogensDict_dict[aa_type][C_name][0]
            try:
                aatype_CHnucleiType_presenceProbSum_mdict[aa_type][C_name+"_"+H_name] = (args.H_weight * Probability_CS[aa_type][H_name] + args.C_weight * Probability_CS[aa_type][C_name])/float((args.H_weight + args.C_weight))  # 2nd WAY
            except TypeError:
                print("ERROR: cannot find either H_name or C_name in Probability_CS multidict!")
                print("DEBUG: aa_type=", aa_type, "H_name=", H_name, "C_name=", C_name, "Probability_CS[aa_type]=", Probability_CS[aa_type])
                sys.exit(1)
            aatype_CHnucleiType_NresonIndex_mdict[aa_type][C_name+"_"+H_name] = HCNH_reson_index

    for aa_type in list(aatype_CHnucleiType_presenceProbSum_mdict.keys()):
       for k,v in list(aatype_CHnucleiType_presenceProbSum_mdict[aa_type].items()):
           print("DEBUG: aa_type=", aa_type, k, v)
    for aa_type in list(aatype_CHnucleiType_NresonIndex_mdict.keys()):
       for k,v in list(aatype_CHnucleiType_NresonIndex_mdict[aa_type].items()):
           print("DEBUG: aa_type=", aa_type, k, v)

    return aatype_CHnucleiType_presenceProbSum_mdict, aatype_CHnucleiType_NresonIndex_mdict, HCNHindex_ResonancesTuple_dict, HCNHindex_spectrumType_dict


def rename_protons(xeasy_lines_list, sparky_lines_list):
    """
        xeasy_lines_list.append( [atom_index, float(Carbon_resonance), 0.2, Carbon_name, int(residue[3:]), aa_type, spectrum_combo] )
        sparky_lines_list.append( [aa3to1_dict[aa_type] + residue[3:] + Hydrogen_name, Carbon_name, str(int(residue[3:])+1) + "N", float(Hydrogen_resonance), float(Carbon_resonance), float(N_resonance), float(HN_resonance)] )

    """
    print("DEBUG: before sparky_lines_list=", sparky_lines_list)

    aa_type = xeasy_lines_list[0][5]
    proton_names_list = [x[3] for x in xeasy_lines_list if x[3][0]=="H"]    # make a list with all the proton names for this residue
    carbon_names_list = [x[3] for x in xeasy_lines_list if x[3][0]=="C"]    # make a list with all the carbon names for this residue
    proton_names_set = set(proton_names_list)
    carbon_names_set = set(carbon_names_list)
    for carbon_name in carbon_names_set:
        if carbon_name in list(aatype_carbon_nondegenerateHlist_mdict[aa_type].keys()):
            bonded_protons_found_list = []
            for proton_name in aatype_carbon_nondegenerateHlist_mdict[aa_type][carbon_name]:
                if proton_name in proton_names_set:
                    bonded_protons_found_list.append(proton_name)
            # DO NOT RENAME VAL-HG1,HG2 --> VAL-QG1,QG2 and LEU-HD1,HD2-->LEU-QD1,QD2 because the code will be confused! Do it at the very end!!!
            print("DEBUG rename_protons: aa_type=", aa_type, "bonded_protons_found_list=", bonded_protons_found_list, "carbon_name=", carbon_name)
            if len(bonded_protons_found_list) == 1 and not (aa_type in ["LEU", "VAL"] and carbon_name in ["CG1", "CG2", "CD1", "CD2"]):
                proton_name = bonded_protons_found_list[0]
                # correct proton names in xeasy list
                for xeasy_line_list in xeasy_lines_list:
                    if xeasy_line_list[3] == proton_name:
                        if len(proton_name) == 3:
                            xeasy_line_list[3] = "Q" + proton_name[1]
                        else:   # in case of ILE, CG, which has 'HG12' and 'HG13'
                            xeasy_line_list[3] = "Q" + proton_name[1:3]
                # correct proton names in sparky list
                for sparky_line_list in sparky_lines_list:
                    if sparky_line_list[0][-3:] == proton_name:
                        if len(proton_name) == 3:
                            sparky_line_list[0] = sparky_line_list[0][:-3] + "Q" + sparky_line_list[0][-2]
                        else:   # in case of ILE, CG, which has 'HG12' and 'HG13'
                            sparky_line_list[0] = sparky_line_list[0][:-4] + "Q" + sparky_line_list[0][-3:-1]

    print("DEBUG: after sparky_lines_list=", sparky_lines_list)
    return xeasy_lines_list, sparky_lines_list


def update_HCNH_peaks(i_residue, matched_HCNH_residue_assignments_dict, protein_alignment_list, absolute_matches_alignment_list, args, sparky_lines_list=[]):
    """
        FUNCTION to update the HCNH peaks in matched_HCNH_residue_assignments_dict by adding labels after the assignemnt. It also reverts
        unassigned 'i+1' labels of residue i back to '?' and 'i-1' of residue i+1 back to '?'.
    """

    tolH = 0.02
    tolC = 0.2

    for line in sparky_lines_list:
        Hname = line[0].replace(i_residue, "")
        Cname = line[1]
        Hreson = line[3]
        Creson = line[4]
        #print "DEBUG: updating HCNH peaks using sparky line=", line
        ## Label peaks assigned in residue i
        for index in range(len(matched_HCNH_residue_assignments_dict[i_residue])):
            peak = matched_HCNH_residue_assignments_dict[i_residue][index]
            if peak[0] != i_residue and peak[2] == Creson and peak[4] == Hreson and peak[1] != Cname and peak[3] != Hname:    # matches also 'i+1' peak labels
                print("Updating the HCNH peaks of i_residue", i_residue, ", peak=", peak)
                matched_HCNH_residue_assignments_dict[i_residue][index][0] = i_residue
                matched_HCNH_residue_assignments_dict[i_residue][index][1] = Cname
                matched_HCNH_residue_assignments_dict[i_residue][index][3] = Hname
                break
        ## Label the i-1 peaks that match with the peaks assigned in residue i
        try:
            iminus1_residue = get_iminus1_residue_from_i_residue(i_residue, protein_alignment_list, absolute_matches_alignment_list, args)
            print("DEBUG: iminus1_residue=", iminus1_residue, "i_residue=", i_residue)
            if not iminus1_residue in ['N/A', '-'] and iminus1_residue in list(matched_HCNH_residue_assignments_dict.keys()):
                for index in range(len(matched_HCNH_residue_assignments_dict[iminus1_residue])):
                    peak = matched_HCNH_residue_assignments_dict[iminus1_residue][index]
                    print("DEBUG: iminus1_peak=", peak, "i_peak=", [i_residue, Cname, Creson, Hname, Hreson])
                    if approx_equal(peak[2], Creson, tolC) and approx_equal(peak[4], Hreson, tolH) and peak[1] != Cname and peak[3] != Hname:
                        print("Updating the HCNH peaks of iminus1_residue", iminus1_residue, ", peak=", peak)
                        matched_HCNH_residue_assignments_dict[iminus1_residue][index][0] = i_residue
                        matched_HCNH_residue_assignments_dict[iminus1_residue][index][1] = Cname
                        matched_HCNH_residue_assignments_dict[iminus1_residue][index][3] = Hname
                        break
        except ValueError:  # if the iminus_residue is not in the absolute_matches_alignment_list, continue
            #print "DEBUG: i_residue is not in the absolute_matches_alignment_list, there is not need to update labels 'i-1' in residue i+1"
            pass
        ## Label the i+1 peaks that match with the peaks assigned in residue i. If no i+1 residue in the alignment, then skip the next step
        ## that removes "i+1" and "i-1" labels.
        try:
            #iplus1_residue = absolute_matches_alignment_list[absolute_matches_alignment_list.index(i_residue)+1]
            iplus1_residue = get_iplus1_residue_from_i_residue(i_residue, protein_alignment_list, absolute_matches_alignment_list, args)
            if iplus1_residue in ['N/A', '-'] or not iplus1_residue in list(matched_HCNH_residue_assignments_dict.keys()):
                #print "DEBUG: iplus1_residue in ['N/A', '-'] ", iplus1_residue
                return matched_HCNH_residue_assignments_dict
        except ValueError:  # if the  i_residue is not in the absolute_matches_alignment_list, there is not need to update labels 'i-1' in residue i+1
            #print "DEBUG: i_residue is not in the absolute_matches_alignment_list, there is not need to update labels 'i-1' in residue i+1"
            return matched_HCNH_residue_assignments_dict
        except IndexError:  # if i_residue is the last residue
            return matched_HCNH_residue_assignments_dict
        for index in range(len(matched_HCNH_residue_assignments_dict[iplus1_residue])):
            peak = matched_HCNH_residue_assignments_dict[iplus1_residue][index]
            if peak[0] == 'i-1' and approx_equal(peak[2], Creson, tolC) and approx_equal(peak[4], Hreson, tolH) and peak[1] != Cname and peak[3] != Hname:
                print("Updating the HCNH peaks of iplus1_residue", iplus1_residue, ", peak=", peak)
                matched_HCNH_residue_assignments_dict[iplus1_residue][index][0] = i_residue
                matched_HCNH_residue_assignments_dict[iplus1_residue][index][1] = Cname
                matched_HCNH_residue_assignments_dict[iplus1_residue][index][3] = Hname
                break

    # now remove labels 'i+1' and 'i-1' from unassigned peaks in residue i and i+1, respectively
    try:
        #iplus1_residue = absolute_matches_alignment_list[absolute_matches_alignment_list.index(i_residue)+1]
        iplus1_residue = get_iplus1_residue_from_i_residue(i_residue, protein_alignment_list, absolute_matches_alignment_list, args)
        if iplus1_residue in ['N/A', '-'] or not iplus1_residue in list(matched_HCNH_residue_assignments_dict.keys()):
            #print "DEBUG: iplus1_residue in ['N/A', '-'] or not iplus1_residue in matched_HCNH_residue_assignments_dict.keys()", iplus1_residue
            return matched_HCNH_residue_assignments_dict
    except ValueError:  # if the i_residue is not in the absolute_matches_alignment_list, there is not need to revent labels 'i+1' and 'i-1' to '?'
        #print "DEBUG: i_residue is not in the absolute_matches_alignment_list, there is not need to revent labels 'i+1' and 'i-1' to '?'"
        return matched_HCNH_residue_assignments_dict
    if '-' in [i_residue, iplus1_residue] or 'N/A' in [i_residue, iplus1_residue] or not iplus1_residue in list(matched_HCNH_residue_assignments_dict.keys()):  # if there is a gap in the 1 or i+1 position, there is no need to do anything
        #print  "DEBUG: '-' in [i_residue, iplus1_residue] or 'N/A' in [i_residue, iplus1_residue] or not iplus1_residue in matched_HCNH_residue_assignments_dict.keys()"
        return matched_HCNH_residue_assignments_dict
    else:   # otherwise revent labels 'i+1' and 'i-1' to '?'
        #print "DEBUG: before reverting 'i+1' labels matched_HCNH_residue_assignments_dict[i_residue]=", matched_HCNH_residue_assignments_dict[i_residue]
        #print "DEBUG: before reverting 'i-1' labels matched_HCNH_residue_assignments_dict[iplus1_residue]=", matched_HCNH_residue_assignments_dict[iplus1_residue]
        for index in range(len(matched_HCNH_residue_assignments_dict[i_residue])):
            peak = matched_HCNH_residue_assignments_dict[i_residue][index]
            if peak[0] == 'i+1':
                print("DEBUG: Reverting to '?' the 'i+1' HCNH peaks of i_residue", i_residue, ", peak=", peak)
                matched_HCNH_residue_assignments_dict[i_residue][index][0] = '?'
        for index in range(len(matched_HCNH_residue_assignments_dict[iplus1_residue])):
            peak = matched_HCNH_residue_assignments_dict[iplus1_residue][index]
            if peak[0] == 'i-1':
                print("DEBUG: Reverting to '?' the 'i-1' HCNH peaks of iplus1_residue", iplus1_residue, ", peak=", peak)
                matched_HCNH_residue_assignments_dict[iplus1_residue][index][0] = '?'

    return matched_HCNH_residue_assignments_dict

def remove_iplus1_labels(i_residue,
                         matched_HCNH_residue_assignments_dict,
                         protein_alignment_list,
                         absolute_matches_alignment_list,
                         args):

    try:
        #iplus1_residue = absolute_matches_alignment_list[absolute_matches_alignment_list.index(i_residue)+1]
        iplus1_residue = get_iplus1_residue_from_i_residue(i_residue,
                                                           protein_alignment_list,
                                                           absolute_matches_alignment_list,
                                                           args)
        if iplus1_residue in ['N/A', '-'] or not iplus1_residue in list(matched_HCNH_residue_assignments_dict.keys()):
            #print "DEBUG: iplus1_residue in ['N/A', '-'] or not iplus1_residue in matched_HCNH_residue_assignments_dict.keys()", iplus1_residue
            return matched_HCNH_residue_assignments_dict
    except ValueError:  # if the i_residue is not in the absolute_matches_alignment_list, there is not need to revent labels 'i+1' and 'i-1' to '?'
        #print "DEBUG: i_residue is not in the absolute_matches_alignment_list, there is not need to revent labels 'i+1' and 'i-1' to '?'"
        return matched_HCNH_residue_assignments_dict
    if '-' in [i_residue, iplus1_residue] or 'N/A' in [i_residue, iplus1_residue] or not iplus1_residue in list(matched_HCNH_residue_assignments_dict.keys()):  # if there is a gap in the 1 or i+1 position, there is no need to do anything
        #print  "DEBUG: '-' in [i_residue, iplus1_residue] or 'N/A' in [i_residue, iplus1_residue] or not iplus1_residue in matched_HCNH_residue_assignments_dict.keys()"
        return matched_HCNH_residue_assignments_dict
    else:   # otherwise revent labels 'i+1' and 'i-1' to '?'
        for index in range(len(matched_HCNH_residue_assignments_dict[i_residue])):
            peak = matched_HCNH_residue_assignments_dict[i_residue][index]
            if peak[0] == 'i+1':
                print("DEBUG: Reverting to '?' the 'i+1' HCNH peaks of i_residue", i_residue, ", peak=", peak)
                matched_HCNH_residue_assignments_dict[i_residue][index][0] = '?'


        for index in range(len(matched_HCNH_residue_assignments_dict[iplus1_residue])):
            peak = matched_HCNH_residue_assignments_dict[iplus1_residue][index]
            if peak[0] == 'i-1':
                print("DEBUG: Reverting to '?' the 'i-1' HCNH peaks of iplus1_residue", iplus1_residue, ", peak=", peak)
                matched_HCNH_residue_assignments_dict[iplus1_residue][index][0] = '?'
    return matched_HCNH_residue_assignments_dict

def update_flanked_HCNH_peaks(i_residue,
                               matched_HCNH_residue_assignments_dict,
                               protein_alignment_list,
                               absolute_matches_alignment_list,
                               args):
    """
        FUNCTION to update the HCNH peaks in matched_HCNH_residue_assignments_dict of the flanked residue (see write_flanked_residues())
        by adding labels from residues i+1 (the labels of i-1 have been already added to the HCNH peaks of i through update_HCNH_peaks()).

        i_residue:  the flanked residue
    """

    tolH = 0.02
    tolC = 0.2

    iplus1_residue = get_iplus1_residue_from_i_residue(i_residue, protein_alignment_list, absolute_matches_alignment_list, args)

    for iplus1_peak in matched_HCNH_residue_assignments_dict[iplus1_residue]:
        if iplus1_peak[0] == '?':
            continue
        for index in range(len(matched_HCNH_residue_assignments_dict[i_residue])):
            i_peak = matched_HCNH_residue_assignments_dict[i_residue][index]
            if i_peak[0] == '?' and approx_equal(i_peak[2], iplus1_peak[2], tolC) and approx_equal(i_peak[4], iplus1_peak[4], tolH):
                #print "Updating the HCNH peaks of iplus1_residue", iplus1_residue, ", peak=", peak
                i_peak[0] = iplus1_peak[0]  # add residue label; matched_HCNH_residue_assignments_dict[i_residue] will be updated automatically
                i_peak[1] = iplus1_peak[1]  # add label
                i_peak[3] = iplus1_peak[3]  # add label
                break
    return matched_HCNH_residue_assignments_dict


def append_unmatched_TOCSY_peaks(residue,
                                 xeasy_lines_list,
                                 atom_index,
                                 resid,
                                 aa_type,
                                 residue_unmatched_TOCSY_peaks_dict):
    """
        FUNCTION to write the unmatched TOCSY peaks.
    """

    if residue in residue_unmatched_TOCSY_peaks_dict.keys():
        for TOCSY_peak in residue_unmatched_TOCSY_peaks_dict[residue]:
            Hreson = TOCSY_peak[4]
            Hname = TOCSY_peak[3]
            Creson = TOCSY_peak[2]
            Cname = TOCSY_peak[1]
            xeasy_lines_list.append( [atom_index, float(Hreson), 0.02, Hname, resid, aa_type, 'TOCSY (unmatched)'] )  # append the Hydrogen resonance line
            atom_index += 1
            xeasy_lines_list.append( [atom_index, float(Creson), 0.2, Cname, resid, aa_type, 'TOCSY (unmatched)'] )  # append the Carbon resonance line
            atom_index += 1

    return xeasy_lines_list, atom_index


def write_TOCSY_HCNH_matched_peaks(i_residue,
                                    atom_index,
                                    matched_HCNH_residue_assignments_dict,
                                    xeasy_fout,
                                    sparky_fout,
                                    protein_alignment_list,
                                    absolute_matches_alignment_list,
                                    args,
                                    TOCSY_assigned_residue_list,
                                    original_HCNH_peaks_dict,
                                    residue_unmatched_TOCSY_peaks_dict,
                                    HCNH_residue_NHresonances_dict,
                                    residues_with_written_NH_list,
                                    written_nucleiNames_list=[]):


    print("DEBUG write_TOCSY_HCNH_matched_peaks: i_residue=", i_residue, "written_nucleiNames_list=", written_nucleiNames_list)
    ##
    ## NOW WRITE THE TOCSY-HCNH MATCHED PEAKS. THESE HAVE ALREADY THE CORRECT PROTON NAMES SINCE THESE WERE TRANSFERED FROM TOCSY.
    ##
    TOCSY_assigned_residue_list.append(i_residue)   # write down that you have assigned the peaks of this residue already
    # resid = absolute_matches_alignment_list.index(i_residue) + int(args.FIRST_RESIDUE_NUMBER)
    resid = get_resid_from_residue(i_residue)
    aa_type = aa1to3_dict[get_aa_type_from_residue(i_residue)]
    xeasy_lines_list, sparky_lines_list = [], []    # lines to save the xeasy and sparky lines for modification
    Creson_Cname_lines_list = [] # list with the resonances and names of the printed carbons to avoid double printing (e.g. CB, CG, CD, CE)
    try:    # DANGEROUS but works for the N-term residue (A200 in aLP)
        matched_HCNH_peaks_list = [peak for peak in matched_HCNH_residue_assignments_dict[i_residue] if peak[0]==i_residue]
    except KeyError:
        return False    # function has been terminated prematurily
    for peak in matched_HCNH_peaks_list:
        Carbon_name = peak[1]
        Carbon_resonance = peak[2]
        Hydrogen_name = peak[3]
        Hydrogen_resonance = peak[4]
        N_resonance = [peak[3] for peak in original_HCNH_peaks_dict[i_residue] if peak[1]==Hydrogen_resonance and peak[2]==Carbon_resonance][0]
        HN_resonance = [peak[4] for peak in original_HCNH_peaks_dict[i_residue] if peak[1]==Hydrogen_resonance and peak[2]==Carbon_resonance][0]
        Creson_Cname_lines_list.append([float(Carbon_resonance), Carbon_name])
        if not Hydrogen_name in written_nucleiNames_list:
            xeasy_lines_list.append( [atom_index, float(Hydrogen_resonance), 0.02, Hydrogen_name, resid, aa_type, "TOCSY-HCNH"] )  # append the Hydrogen resonance line
            atom_index += 1
            try:
                #sparky_fout.write("\t%s-%s-%s-H\t\t%.3f\t%.3f\t%.3f\t%.3f\t\n" % ( aa3to1_dict[aa_type] + residue[3:] + Hydrogen_name, Carbon_name, protein_alignment_list[residue2position_dict[residue]] + str(int(residue[3:])+1) + "N",float(Hydrogen_resonance), float(Carbon_resonance), float(N_resonance), float(HN_resonance)) )
                sparky_lines_list.append( [i_residue + Hydrogen_name, Carbon_name, i_residue + "N",float(Hydrogen_resonance), float(Carbon_resonance), float(N_resonance), float(HN_resonance)] )
            except IndexError:
                #sparky_fout.write("\t%s-%s-%s-H\t\t%.3f\t%.3f\t%.3f\t%.3f\t\n" % ( aa3to1_dict[aa_type] + residue[3:] + Hydrogen_name, Carbon_name, str(int(residue[3:])+1) + "N", float(Hydrogen_resonance), float(Carbon_resonance), float(N_resonance), float(HN_resonance)) )
                sparky_lines_list.append( [i_residue + Hydrogen_name, Carbon_name, i_residue + "N", float(Hydrogen_resonance), float(Carbon_resonance), float(N_resonance), float(HN_resonance)] )
    Cnames_list = [cl[1] for cl in Creson_Cname_lines_list]
    Cnames_set = set(Cnames_list)
    for Carbon_name in Cnames_set:
        if Cnames_list.count(Carbon_name) > 1:  # if this carbon occurs in multiple C-H pairs, calculate its average resonance for the XEASY format
            Carbon_resonance = round(np.average([cl[0] for cl in Creson_Cname_lines_list if cl[1]==Carbon_name]), 3)
        else:
            Carbon_resonance = round(np.average([cl[0] for cl in Creson_Cname_lines_list if cl[1]==Carbon_name]), 3)
        #xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s\n" % (atom_index, float(Carbon_resonance), 0.2, Carbon_name, int(residue[3:]), aa_type))
        if not Carbon_name in written_nucleiNames_list:
            xeasy_lines_list.append( [atom_index, float(Carbon_resonance), 0.2, Carbon_name, resid, aa_type, "TOCSY-HCNH"] )
            atom_index += 1
    # append unmatched TOCSY peaks, if any
    if i_residue in residue_unmatched_TOCSY_peaks_dict.keys(): print("DEBUG: point 2, appending unmatched TOCSY peaks for residue", i_residue, residue_unmatched_TOCSY_peaks_dict[i_residue])
    xeasy_lines_list, atom_index = append_unmatched_TOCSY_peaks(i_residue, xeasy_lines_list, atom_index, resid, aa_type, residue_unmatched_TOCSY_peaks_dict)


    # WRITE C AND H LINES IN XEASY FORMAT
    print("DEBUG: writing TOCSY-HCNH matched peaks xeasy_lines_list=", xeasy_lines_list)
    for xeasy_line in xeasy_lines_list:
        xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s %s\n" % (xeasy_line[0], xeasy_line[1], xeasy_line[2], xeasy_line[3], xeasy_line[4], xeasy_line[5], xeasy_line[6]) )
    # WRITE C AND H LINES IN SPARKY FORMAT
    for sparky_line in sparky_lines_list:
        sparky_fout.write("\t%s-%s-%s-H\t\t%.3f\t%.3f\t%.3f\t%.3f\t\n" % (sparky_line[0], sparky_line[1], sparky_line[2], sparky_line[3], sparky_line[4], sparky_line[5], sparky_line[6]))
    matched_HCNH_residue_assignments_dict = update_HCNH_peaks(i_residue, matched_HCNH_residue_assignments_dict, protein_alignment_list, absolute_matches_alignment_list, args, sparky_lines_list) # update the labels in matched_HCNH_residue_assignments_dict

    # FINALLY WRITE N AND H LINES IN XEASY FORMAT
    if i_residue in list(HCNH_residue_NHresonances_dict.keys()): # only if this residue has N-H
        average_HN_resonance = HCNH_residue_NHresonances_dict[i_residue][1]
        stdev_HN_resonance = HCNH_residue_NHresonances_dict[i_residue][3]
        average_N_resonance = HCNH_residue_NHresonances_dict[i_residue][0]
        stdev_N_resonance = HCNH_residue_NHresonances_dict[i_residue][2]
        try:
            xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s %s\n" % (atom_index, float(average_HN_resonance), stdev_HN_resonance, "H", resid, aa_type, '')) # The N-HN of i sees the C-H of i-1
            atom_index += 1
            #print "\t"+str(atom_index)+"\t"+str(float(average_N_resonance))+"\t0.2\tN\t"+str(int(residue[3:])+1)
            xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s %s\n" % (atom_index, average_N_resonance, stdev_N_resonance, "N", resid, aa_type, ''))   # The N-HN of i sees the C-H of i-1
            atom_index += 1
        except IndexError:
            xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s %s\n" % (atom_index, float(average_HN_resonance), stdev_HN_resonance, "H", resid, aa_type, '')) # The N-HN of i sees the C-H of i-1
            atom_index += 1
            #print "\t"+str(atom_index)+"\t"+str(float(average_N_resonance))+"\t0.2\tN\t"+str(int(residue[3:])+1)
            xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s %s\n" % (atom_index, average_N_resonance, stdev_N_resonance, "N", resid, aa_type, ''))   # The N-HN of i sees the C-H of i-1
            atom_index += 1
        residues_with_written_NH_list.append(i_residue)

    return True


def select_best_HCNH_peak_combination(i_residue,
                                       all_possible_assignments_list,
                                       partly_assigned_carbons_list,
                                       atom_index,
                                       xeasy_fout,
                                       sparky_fout,
                                       protein_alignment_list,
                                       absolute_matches_alignment_list,
                                       args,
                                       HCNH_residue_NHresonances_dict,
                                       matched_HCNH_residue_assignments_dict,
                                       TOCSY_assigned_residue_list,
                                       iter1_allowed_atoms_dict,
                                       iteration=1):

    from .global_vars import Prob_CS
    Probability_CS = Prob_CS

    #print "DEBUG: i_residue=", i_residue, "all_possible_assignments_list=", all_possible_assignments_list
    ## SELECT THE BEST ASSIGNMENT COMBINATION AND SAVE IT INTO XEASY AND SPARKY FILES
    written_nucleiNames_list = []   # list with the nuclei names that will be written to files (initialized as empty, will be populated later)
    residue_residueTypesProbTupleList_dict = OrderedDict() # ordereddict with keys the AAIG of residue i and values lists of tuples of the form (residue i-1 matching aa type, average probability)
    previous_aaindex = ""
    print("Assigning HCNH chemical shifts to i_residue ", i_residue, "args=", args)
    #aa_type = residue[0:3]  # aa type in 3-letter code
    aa_type = aa1to3_dict[get_aa_type_from_residue(i_residue)]
    aatype_CHnucleiType_presenceProbSum_mdict, \
    aatype_CHnucleiType_NresonIndex_mdict, \
    HCNHindex_ResonancesTuple_dict, \
    HCNHindex_spectrumType_dict = \
        get_aatypes_from_all_H_C_resonpairs(all_possible_assignments_list,
                                            i_residue,
                                            partly_assigned_carbons_list,
                                            iteration,
                                            Probability_CS,
                                            HCNH_residue_NHresonances_dict,
                                            matched_HCNH_residue_assignments_dict,
                                            protein_alignment_list,
                                            absolute_matches_alignment_list,
                                            iter1_allowed_atoms_dict,
                                            args)
    if aatype_CHnucleiType_presenceProbSum_mdict != {} and aatype_CHnucleiType_NresonIndex_mdict != {}: # if no atom assignment could be made, skip this aa! (TEMPORARY FIX)
        ## WRITE THE CHEMICAL SHIFTS ASSIGNMENTS IN XEASY FORMAT FOR ROSETTA
        # resid = absolute_matches_alignment_list.index(i_residue) + int(args.FIRST_RESIDUE_NUMBER) # OLD WAY; returns ValueError for A335
        resid = get_resid_from_residue(i_residue)
        xeasy_lines_list, sparky_lines_list = [], []    # lines to save the xeasy and sparky lines for modification
        HN_resonances_list, N_resonances_list = [], []
        Creson_Cname_NresonIndex_lines_list = [] # list with the resonances, names and Nreson indices of the printed carbons to avoid double printing (e.g. CB, CG, CD, CEfilter
        print("DEBUG: aatype_CHnucleiType_NresonIndex_mdict[aa_type].keys() =", list(aatype_CHnucleiType_NresonIndex_mdict[aa_type].keys()))
        print("DEBUG: aatype_CHnucleiType_NresonIndex_mdict =", aatype_CHnucleiType_NresonIndex_mdict)
        for CHnucleiType, NresonIndex  in list(aatype_CHnucleiType_NresonIndex_mdict[aa_type].items()):
            ResonancesTuple = HCNHindex_ResonancesTuple_dict[NresonIndex]
            #print "DEBUG: aa_type=", aa_type, "CHnucleiType=", CHnucleiType, "NresonIndex=", NresonIndex, "ResonancesTuple=", ResonancesTuple
            Carbon_name = CHnucleiType.split('_')[0]
            Hydrogen_name = CHnucleiType.split('_')[1]
            # Rename un-resolved Geminal protons to QA, QB, QG, QD, QE.
            if len(aa_carbonBondedHydrogensDict_dict[aa_type][Carbon_name]) > 1:    # if this carbon has geminal protons
                #print "DEBUG: aa_type=", aa_type, "Carbon_name=", Carbon_name
                #print "DEBUG: aatype_CHnucleiType_NresonIndex_mdict[aa_type].keys()=", aatype_CHnucleiType_NresonIndex_mdict[aa_type].keys()
                #print "DEBUG: aa_carbonBondedHydrogensDict_dict[aa_type][Carbon_name]=", aa_carbonBondedHydrogensDict_dict[aa_type][Carbon_name]
                if [Carbon_name+"_" in x for x in list(aatype_CHnucleiType_NresonIndex_mdict[aa_type].keys())].count(True) < len(aa_carbonBondedHydrogensDict_dict[aa_type][Carbon_name]):
                    #print "DEBUG entered condition 1!"
                    if aa_type in list(aatype_carbon_methylHydrogens_mdict.keys()) and Carbon_name in list(aatype_carbon_methylHydrogens_mdict[aa_type].keys()):
                        #print "DEBUG entered condition 2!"
                        new_Hydrogen_name = aatype_carbon_methylHydrogens_mdict[aa_type][Carbon_name]
                        Hydrogen_name = new_Hydrogen_name
                    #else:
                    #    print "DEBUG entered condition 3!"
                    #    new_Hydrogen_name = 'Q'+Hydrogen_name[1]
                    #    Hydrogen_name = new_Hydrogen_name
            # CONDITIONS:
            # 1)
            if ( aa_type in list(aatype_carbon_methylHydrogens_mdict.keys()) and
                Carbon_name in list(aatype_carbon_methylHydrogens_mdict[aa_type].keys()) and
                (not aa_type in list(aa_CarbonListwithDegenerateH_dict.keys()) or
                Carbon_name in aa_CarbonListwithDegenerateH_dict[aa_type]) ):
                # DO NOT RENAME VAL-HG1,HG2 --> VAL-QG1,QG2 and LEU-HD1,HD2-->LEU-QD1,QD2 because the code will be confused! Do it at the very end!!!
                if not(aa_type in ["LEU", "VAL"] and Carbon_name in ["CG1", "CG2", "CD1", "CD2"]):
                    print("DEBUG entered condition 4! aa_type = ", aa_type, "Carbon_name=", Carbon_name)
                    new_Hydrogen_name = aatype_carbon_methylHydrogens_mdict[aa_type][Carbon_name]
                    Hydrogen_name = new_Hydrogen_name
            Carbon_resonance = ResonancesTuple[1]
            Hydrogen_resonance = ResonancesTuple[0]
            N_resonance = ResonancesTuple[2]
            HN_resonance = ResonancesTuple[3]
            #HN_resonances_list.append(float(ResonancesTuple[3]))   # UNECESSARY CAUSE WE ALREADY SAVED THEM IN HCNH_residue_NHresonances_dict
            #N_resonances_list.append(float(ResonancesTuple[2]))

            #print "\t"+str(atom_index)+"\t"+str(Carbon_resonance)+"\t0.2\t"+str(Carbon_name)+"\t"+str(residue[3:])
            Creson_Cname_NresonIndex_lines_list.append([float(Carbon_resonance), Carbon_name, NresonIndex])

            #print "\t"+str(atom_index)+"\t"+str(Hydrogen_resonance)+"\t0.02\t"+str(Hydrogen_name)+"\t"+str(residue[3:])
            #xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s\n" % (atom_index, float(Hydrogen_resonance), 0.02, Hydrogen_name, int(residue[3:]), aa_type))
            if HCNHindex_spectrumType_dict[NresonIndex] in ['TOCSY-HCNH', 'TOCSY (unmatched)', 'TOCSY (unmatched) + TOCSY-HCNH', 'TOCSY-HCNH + TOCSY (unmatched)', 'HCNH iter1 + TOCSY-HCNH', 'TOCSY-HCNH + HCNH iter1'] or 'iter'+str(iteration-1) in HCNHindex_spectrumType_dict[NresonIndex]:  # no iteration# comment at TOCSY-HCNH matched peaks
                xeasy_lines_list.append( [atom_index, float(Hydrogen_resonance), 0.02, Hydrogen_name, resid, aa_type, HCNHindex_spectrumType_dict[NresonIndex] ] )  # append the Hydrogen resonance line
            else:   # if not TOCSY-HCNH matched peak, write in which iteration it was assigned
                xeasy_lines_list.append( [atom_index, float(Hydrogen_resonance), 0.02, Hydrogen_name, resid, aa_type, HCNHindex_spectrumType_dict[NresonIndex] + " iter"+str(iteration)] )  # append the Hydrogen resonance line
                print("DEBUG: point 1 appending iter: ", [atom_index, float(Hydrogen_resonance), 0.02, Hydrogen_name, resid, aa_type, HCNHindex_spectrumType_dict[NresonIndex] + " iter"+str(iteration)])
            atom_index += 1
            if N_resonance != None and HN_resonance != None:    # if no N-H was found (e.g. N-term residue) don't write the i-i Sparky lines for this residue, only the xeasy lines
                try:
                    #sparky_fout.write("\t%s-%s-%s-H\t\t%.3f\t%.3f\t%.3f\t%.3f\t\n" % ( aa3to1_dict[aa_type] + residue[3:] + Hydrogen_name, Carbon_name, protein_alignment_list[residue2position_dict[residue]] + str(int(residue[3:])+1) + "N",float(Hydrogen_resonance), float(Carbon_resonance), float(N_resonance), float(HN_resonance)) )
                    sparky_lines_list.append( [i_residue + Hydrogen_name, Carbon_name, i_residue + "N",float(Hydrogen_resonance), float(Carbon_resonance), float(N_resonance), float(HN_resonance)] )
                except IndexError:
                    #sparky_fout.write("\t%s-%s-%s-H\t\t%.3f\t%.3f\t%.3f\t%.3f\t\n" % ( aa3to1_dict[aa_type] + residue[3:] + Hydrogen_name, Carbon_name, str(int(residue[3:])+1) + "N", float(Hydrogen_resonance), float(Carbon_resonance), float(N_resonance), float(HN_resonance)) )
                    sparky_lines_list.append( [i_residue + Hydrogen_name, Carbon_name, i_residue + "N", float(Hydrogen_resonance), float(Carbon_resonance), float(N_resonance), float(HN_resonance)] )


        Cnames_list = [cl[1] for cl in Creson_Cname_NresonIndex_lines_list]
        Cnames_set = set(Cnames_list)
        for Carbon_name in Cnames_set:
            if Cnames_list.count(Carbon_name) > 1:  # if this carbon occurs in multiple C-H pairs, calculate its average resonance for the XEASY format
                Carbon_resonance = round(np.average([cl[0] for cl in Creson_Cname_NresonIndex_lines_list if cl[1]==Carbon_name]), 3)
                NresonIndex_list = [cl[2] for cl in Creson_Cname_NresonIndex_lines_list if cl[1]==Carbon_name]
                spectrum_combo_list = [HCNHindex_spectrumType_dict[NRI] for NRI in NresonIndex_list]    # a list of the spectrum type of each instance of this C-type
                if len(set(spectrum_combo_list)) > 1:  # if this C-type has resonances that come from more that one spectrum source
                    spectrum_combo = spectrum_combo_list[0]
                    for st in spectrum_combo_list[1:]:
                        spectrum_combo += " + " + st  # append all the spectrum sources
                else:   # if this C-type has resonances that come from only one spectrum source
                    spectrum_combo = spectrum_combo_list[0]   # save only that spectrum source
            else:   # if this carbon occurs only in one C-H pair, save  that resonance in the XEASY format
                Carbon_resonance = [cl[0] for cl in Creson_Cname_NresonIndex_lines_list if cl[1]==Carbon_name][0]
                NresonIndex = [cl[2] for cl in Creson_Cname_NresonIndex_lines_list if cl[1]==Carbon_name][0]
                spectrum_combo = HCNHindex_spectrumType_dict[NresonIndex]
            #xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s\n" % (atom_index, float(Carbon_resonance), 0.2, Carbon_name, int(residue[3:]), aa_type))
            if spectrum_combo in ['TOCSY-HCNH', 'TOCSY (unmatched)', 'TOCSY (unmatched) + TOCSY-HCNH', 'TOCSY-HCNH + TOCSY (unmatched)', 'HCNH iter1 + TOCSY-HCNH'] or 'iter'+str(iteration-1) in HCNHindex_spectrumType_dict[NresonIndex]:  # no iteration# comment at TOCSY-HCNH matched peaks
                xeasy_lines_list.append( [atom_index, float(Carbon_resonance), 0.2, Carbon_name, resid, aa_type, spectrum_combo ] )
            else:   # if not TOCSY-HCNH matched peak, write in which iteration it was assigned
                xeasy_lines_list.append( [atom_index, float(Carbon_resonance), 0.2, Carbon_name, resid, aa_type, spectrum_combo + " iter"+str(iteration) ] )
                print("DEBUG: point 2 appending iter: ", [atom_index, float(Carbon_resonance), 0.2, Carbon_name, resid, aa_type, spectrum_combo + " iter"+str(iteration) ])
            atom_index += 1
        # append unmatched TOCSY peaks, if any
        # if i_residue in residue_unmatched_TOCSY_peaks_dict.keys(): print "DEBUG: point 1, appending unmatched TOCSY peaks for residue", i_residue, residue_unmatched_TOCSY_peaks_dict[i_residue]
        # xeasy_lines_list, atom_index = append_unmatched_TOCSY_peaks(i_residue, xeasy_lines_list, atom_index, resid, aa_type)

        # RENAME DEGENERATE PROTONS
        new_xeasy_lines_list, new_sparky_lines_list = rename_protons(xeasy_lines_list, sparky_lines_list)
        written_nucleiNames_list = [xline[3] for xline in new_xeasy_lines_list] # list with updated Hydrogen names, to prevent duplicate writing
        print("DEBUG: point 1 new_xeasy_lines_list=", new_xeasy_lines_list)
        # WRITE C AND H LINES IN XEASY FORMAT
        for xeasy_line in new_xeasy_lines_list:
            formated_xeasy_line = "\t%g\t%.3f\t%.3f\t%s\t%g\t# %s %s\n" % (xeasy_line[0], xeasy_line[1], xeasy_line[2], xeasy_line[3], xeasy_line[4], xeasy_line[5], xeasy_line[6])
            print("DEBUG select_best_HCNH_peak_combination: formated_xeasy_line=", formated_xeasy_line)
            xeasy_fout.write(formated_xeasy_line)
        # WRITE C AND H LINES IN SPARKY FORMAT
        print("DEBUG: point 1 new_sparky_lines_list=", new_sparky_lines_list)
        for sparky_line in new_sparky_lines_list:
            sparky_fout.write("\t%s-%s-%s-H\t\t%.3f\t%.3f\t%.3f\t%.3f\t\n" % (sparky_line[0], sparky_line[1], sparky_line[2], sparky_line[3], sparky_line[4], sparky_line[5], sparky_line[6]))
        # MODIFY matched_HCNH_residue_assignments_dict
        matched_HCNH_residue_assignments_dict = update_HCNH_peaks(i_residue,
                                                                    matched_HCNH_residue_assignments_dict,
                                                                    protein_alignment_list,
                                                                    absolute_matches_alignment_list,
                                                                    args,
                                                                    new_sparky_lines_list)    # update the labels in matched_HCNH_residue_assignments_dict
        #print "\t"+str(atom_index)+"\t"+str(float(average_HN_resonance))+"\t0.02\tH\t"+str(int(residue[3:])+1)
        TOCSY_assigned_residue_list.append(i_residue)   # write down that you have assigned the peaks of this residue already
    else:
        # remove "i+1" and "i-1" labels in case this is ITERATION 1 or 2
        matched_HCNH_residue_assignments_dict = update_HCNH_peaks(i_residue,
                                                                    matched_HCNH_residue_assignments_dict,
                                                                    protein_alignment_list,
                                                                    absolute_matches_alignment_list,
                                                                    args)    # update the labels in matched_HCNH_residue_assignments_dict

    print("DEBUG select_best_HCNH_peak_combination: returning written_nucleiNames_list=", written_nucleiNames_list)
    return written_nucleiNames_list, matched_HCNH_residue_assignments_dict, TOCSY_assigned_residue_list


def write_flanked_residues(absolute_matches_alignment_list,
                           matched_HCNH_residue_assignments_dict,
                           TOCSY_residue_assignments_dict,
                           atom_index,
                           xeasy_fout,
                           sparky_fout,
                           patched_residues_list,
                           clean_resid_assignments_dict,
                           protein_alignment_list,
                           args,
                           HCNH_residue_peak_intensity_mdict,
                           aa2pthres_dict,
                           aa2ithres_dict,
                           checked_residues_set,
                           matched_i_to_iplus1_peaks_mdict,
                           HCNH_residue_NHresonances_dict,
                           residues_with_written_NH_list,
                           residues_with_written_prevAssigned_peaks,
                           revordered_residue_keys,
                           iter1_allowed_atoms_dict,
                           allowed_aa_atoms_dict,
                           iteration=2):


    ##
    ## NOW WRITE TO XEASY ONLY FILE THE HCNH PEAKS OF RESIDUES THAT DON'T HAVE TOCSY AT ALL BUT ARE FLANKED BY RESIDUES THAT HAVE BOTH TOCSY
    ## AND HCNH. E.g. in aLP S302 - G303 - R304, G303 has only HCNH peaks but unlike R304 & S302
    ## I.e. if a residue i+1 [R304] doesn't have TOCSY (we cannot get the resonances of residue i [G303]) but both residues i & i+1 have
    ## HCNH peaks. then use only the HCNH peaks to do the assignment for residue i. In the code, the above means that:
    ## 'S302' in TOCSY_residue_assignments_dict.keys() = True
    ## 'G303' in TOCSY_residue_assignments_dict.keys() = False
    ## 'R304' in TOCSY_residue_assignments_dict.keys() = True
    ## 'S302' in matched_HCNH_residue_assignments_dict.keys() = True
    ## 'G303' in matched_HCNH_residue_assignments_dict.keys() = True
    ## 'R304' in matched_HCNH_residue_assignments_dict.keys() = True
    ##


    from .global_vars import Prob_CS
    Probability_CS = Prob_CS

    ## Find the flanked residues
    revordered_flanked_residues = []    # residues with not TOCSY but with HCNH, surrounded by residues with both TOCSY & HCNH
    TOCSY_keys = list(TOCSY_residue_assignments_dict.keys())
    HCNH_keys = list(matched_HCNH_residue_assignments_dict.keys())
    for index in reversed(list(range(1, len(absolute_matches_alignment_list)-1))):   # first and last elements are always gaps
        if (absolute_matches_alignment_list[index] != '-' and
            absolute_matches_alignment_list[index+1] != '-' and
            absolute_matches_alignment_list[index-1] != '-'):
            i_res = absolute_matches_alignment_list[index]
            iplus1_res = absolute_matches_alignment_list[index+1]
            iminus1_res = absolute_matches_alignment_list[index-1]
            if (iminus1_res in HCNH_keys and i_res in HCNH_keys and iplus1_res in HCNH_keys and
                iminus1_res in TOCSY_keys and not i_res in TOCSY_keys and iplus1_res in TOCSY_keys):
                revordered_flanked_residues.append(i_res)

    print("DEBUG: revordered_flanked_residues=", revordered_flanked_residues)
    for i_residue in revordered_flanked_residues:    # both the keys and the C,H assignments correspond to residue i-1
        print("Assigning flanked residue ", i_residue, " only from HCNH...")
        matched_HCNH_residue_assignments_dict = \
            update_flanked_HCNH_peaks(i_residue,
                                       matched_HCNH_residue_assignments_dict,
                                       protein_alignment_list,
                                       absolute_matches_alignment_list,
                                       args)   # first add labels from i+1 (if any label is available) to residue i HCNH peaks
        aa_type = aa1to3_dict[get_aa_type_from_residue(i_residue)]
        all_possible_assignments_list = []
        checked_residues_set.add(i_residue)
        all_allowed_carbons_list = [x for x in allowed_aa_atoms_dict[aa_type] if x[0] == "C"]   # ALL ALLOWED CARBONS for this aa type
        print("DEBUG: 3-letter code aa_type=", aa_type)

        matched_iplus1, \
        matched_HCNH_residue_assignments_dict, \
        matched_i_to_iplus1_peaks_mdict = \
         match_HCNH_i_iplus1_peaks(i_residue,
                                    absolute_matches_alignment_list,
                                    matched_HCNH_residue_assignments_dict,
                                    matched_i_to_iplus1_peaks_mdict)  # add labels to the i->i+1 matched peaks. If i_residue is a C-term residue or
                                                                      # presents a gap, skip it from this iteration
        if matched_iplus1 == False:

            continue

        ## NOTE: we left the flanked residues last deliberately in order to assign as many peaks as possible in residues i+1 and i-1.
        ## These peaks may also exist in the HCNH of residue i and will have labels. However we expect all the carbons of residue i to be
        ## missing and none to be partially assigned, since i has not TOCSY-HCNH matched peaks.
        tmp_assigned_carbons_list = [peak[1] for peak in matched_HCNH_residue_assignments_dict[i_residue] if peak[0]==i_residue and peak[1]!="?" and peak[3]!="?"]   # these are the HCNH assignments of residue i made in previous iterations
        # Remove Carbon names that belong to CH2 and had been found only once in TOCSY (namely 1 C-H peak is missing and must be found in HCNH)
        assigned_carbons_list = []  # list of carbons which have ALL their protons assigned in TOCSY
        partly_assigned_carbons_list = []   # list of methylene carbons for which only one proton frequency was found in TOCSY
        print("DEBUG: tmp_assigned_carbons_list=", tmp_assigned_carbons_list)
        for Cname in set(tmp_assigned_carbons_list): # now decide which carbons were completely and which partly assigned
            if Cname in list(aatype_carbon_nondegenerateHlist_mdict[aa_type].keys()) and tmp_assigned_carbons_list.count(Cname) == 2:
                assigned_carbons_list.append(Cname)
            if Cname in list(aatype_carbon_nondegenerateHlist_mdict[aa_type].keys()) and tmp_assigned_carbons_list.count(Cname) == 1:
                partly_assigned_carbons_list.append(Cname)
            elif not Cname in list(aatype_carbon_nondegenerateHlist_mdict[aa_type].keys()):
                assigned_carbons_list.append(Cname)
        print("DEBUG: assigned_carbons_list=", assigned_carbons_list)
        print("DEBUG: partly_assigned_carbons_list=", partly_assigned_carbons_list)
        missing_carbons_list = []
        for carbon in all_allowed_carbons_list:     # keep only the carbon types of this aa type that were not assigned to peaks
            if not carbon in assigned_carbons_list:
                missing_carbons_list.append(carbon)
        print("DEBUG: dict key i_residue", i_residue, "missing_carbons_list=", missing_carbons_list)

        unassigned_HCNH_peaks_list = []    # ini
        if len(missing_carbons_list) != 0:   # IF THERE ARE MISSING CARBONS, TRY TO FIND THEM (EXCLUDE THE ASSIGNED OR MATCHED HCNH PEAKS)
            # try:    # DANGEROUS but works for the N-term residue (A200 in aLP)
            if iteration in [1, 2]: # in ITER1 & 2 we assign only the "i+1" matched peaks
                unassigned_HCNH_peaks_list = [peak for peak in matched_HCNH_residue_assignments_dict[i_residue] if peak[0]=="i+1" and peak[1]=="?" and peak[3]=="?"]
            else:
                unassigned_HCNH_peaks_list = [peak for peak in matched_HCNH_residue_assignments_dict[i_residue] if peak[0]=="?" and peak[1]=="?" and peak[3]=="?"]
            # except KeyError:
            #     continue
            for HCNH_peak in unassigned_HCNH_peaks_list:
                HCNH_Creson = HCNH_peak[2]
                HCNH_Hreson = HCNH_peak[4]
                all_possible_assignments_list.extend(get_probabilities_from_H_C_resonpair_2Dhist(i_residue,
                                                                                                 HCNH_Hreson,
                                                                                                 HCNH_Creson,
                                                                                                 missing_carbons_list,
                                                                                                 partly_assigned_carbons_list,
                                                                                                 args,
                                                                                                 HCNH_residue_peak_intensity_mdict,
                                                                                                 matched_i_to_iplus1_peaks_mdict,
                                                                                                 aa2pthres_dict,
                                                                                                 aa2ithres_dict,
                                                                                                 protein_alignment_list,
                                                                                                 absolute_matches_alignment_list,
                                                                                                 iteration=iteration))
            [assignment.append('HCNH') for assignment in all_possible_assignments_list]    # we don't need to created a new list, the original all_possible_assignments_list will be updated
            print("DEBUG: i_residue=", i_residue, "all_possible_assignments_list=", all_possible_assignments_list)
            written_nucleiNames_list = []
            i_resid = get_resid_from_residue(i_residue)
            if i_resid in list(clean_resid_assignments_dict.keys()):
                written_nucleiNames_list = [line[3]  for line in clean_resid_assignments_dict[i_resid]]   # list of nuclei that have been written already, to avoid duplication

            if len(all_possible_assignments_list) > 0:
                ## SELECT THE BEST ASSIGNMENT COMBINATION AND SAVE IT INTO XEASY AND SPARKY FILES
                residue_residueTypesProbTupleList_dict = OrderedDict() # ordereddict with keys the AAIG of residue i and values lists of tuples of the form (residue i-1 matching aa type, average probability)
                previous_aaindex = ""
                print("Assigning HCNH chemical shifts to i_residue ", i_residue)
                #aa_type = residue[0:3]  # aa type in 3-letter code
                aa_type = aa1to3_dict[get_aa_type_from_residue(i_residue)]
                partly_assigned_carbons_list = []   # for iteration one there should not be any partly assigned Carbons because this residue has no TOCSY-HCNH matched peaks
                aatype_CHnucleiType_presenceProbSum_mdict, aatype_CHnucleiType_NresonIndex_mdict, HCNHindex_ResonancesTuple_dict, HCNHindex_spectrumType_dict = get_aatypes_from_all_H_C_resonpairs(all_possible_assignments_list, i_residue, partly_assigned_carbons_list, iteration, Probability_CS,
                                        HCNH_residue_NHresonances_dict, matched_HCNH_residue_assignments_dict, protein_alignment_list, absolute_matches_alignment_list, iter1_allowed_atoms_dict, args)
                if aatype_CHnucleiType_presenceProbSum_mdict != {} and aatype_CHnucleiType_NresonIndex_mdict != {}: # if no atom assignment could be made, skip this aa! (TEMPORARY FIX)
                    #print "DEBUG: previous_TOCSY_aaindex=",previous_TOCSY_aaindex,"residue_residueTypesProbTupleList_dict[previous_TOCSY_aaindex]=",residue_residueTypesProbTupleList_dict[previous_TOCSY_aaindex]
                    ## WRITE THE CHEMICAL SHIFTS ASSIGNMENTS IN XEASY FORMAT FOR ROSETTA
                    # resid = absolute_matches_alignment_list.index(i_residue) + int(args.FIRST_RESIDUE_NUMBER) # OLD WAY
                    resid = get_resid_from_residue(i_residue)
                    xeasy_lines_list, sparky_lines_list = [], []    # lines to save the xeasy and sparky lines for modification
                    HN_resonances_list, N_resonances_list = [], []
                    Creson_Cname_lines_list = [] # list with the names of the printed carbons to avoid double printing (e.g. CB, CG, CD, CE)
                    print("DEBUG: aatype_CHnucleiType_NresonIndex_mdict[aa_type].keys() =", list(aatype_CHnucleiType_NresonIndex_mdict[aa_type].keys()))
                    print("DEBUG: aatype_CHnucleiType_NresonIndex_mdict =", aatype_CHnucleiType_NresonIndex_mdict)
                    for CHnucleiType, NresonIndex  in list(aatype_CHnucleiType_NresonIndex_mdict[aa_type].items()):
                        ResonancesTuple = HCNHindex_ResonancesTuple_dict[NresonIndex]
                        #print "DEBUG: aa_type=", aa_type, "CHnucleiType=", CHnucleiType, "NresonIndex=", NresonIndex, "ResonancesTuple=", ResonancesTuple
                        Carbon_name = CHnucleiType.split('_')[0]
                        Hydrogen_name = CHnucleiType.split('_')[1]
                        # Rename un-resolved Geminal protons to QA, QB, QG, QD, QE.
                        if len(aa_carbonBondedHydrogensDict_dict[aa_type][Carbon_name]) > 1:    # if this carbon has geminal protons
                            #print "DEBUG: aa_type=", aa_type, "Carbon_name=", Carbon_name
                            #print "DEBUG: aatype_CHnucleiType_NresonIndex_mdict[aa_type].keys()=", aatype_CHnucleiType_NresonIndex_mdict[aa_type].keys()
                            #print "DEBUG: aa_carbonBondedHydrogensDict_dict[aa_type][Carbon_name]=", aa_carbonBondedHydrogensDict_dict[aa_type][Carbon_name]
                            if [Carbon_name+"_" in x for x in list(aatype_CHnucleiType_NresonIndex_mdict[aa_type].keys())].count(True) < len(aa_carbonBondedHydrogensDict_dict[aa_type][Carbon_name]):
                                #print "DEBUG entered condition 1!"
                                if aa_type in list(aatype_carbon_methylHydrogens_mdict.keys()) and Carbon_name in list(aatype_carbon_methylHydrogens_mdict[aa_type].keys()):
                                    #print "DEBUG entered condition 2!"
                                    new_Hydrogen_name = aatype_carbon_methylHydrogens_mdict[aa_type][Carbon_name]
                                    Hydrogen_name = new_Hydrogen_name
                                #else:
                                #    print "DEBUG entered condition 3!"
                                #    new_Hydrogen_name = 'Q'+Hydrogen_name[1]
                                #    Hydrogen_name = new_Hydrogen_name
                        # CONDITIONS:
                        # 1)
                        if aa_type in list(aatype_carbon_methylHydrogens_mdict.keys()) and Carbon_name in list(aatype_carbon_methylHydrogens_mdict[aa_type].keys()) and (not aa_type in list(aa_CarbonListwithDegenerateH_dict.keys()) or Carbon_name in aa_CarbonListwithDegenerateH_dict[aa_type]):
                            # DO NOT RENAME VAL-HG1,HG2 --> VAL-QG1,QG2 and LEU-HD1,HD2-->LEU-QD1,QD2 because the code will be confused! Do it at the very end!!!
                            if not(aa_type in ["LEU", "VAL"] and Carbon_name in ["CG1", "CG2", "CD1", "CD2"]):
                                #print "DEBUG entered condition 4! aa_type = ", aa_type, "Carbon_name=", Carbon_name
                                new_Hydrogen_name = aatype_carbon_methylHydrogens_mdict[aa_type][Carbon_name]
                                Hydrogen_name = new_Hydrogen_name
                        Carbon_resonance = ResonancesTuple[1]
                        Hydrogen_resonance = ResonancesTuple[0]
                        N_resonance = ResonancesTuple[2]
                        HN_resonance = ResonancesTuple[3]
                        #HN_resonances_list.append(float(ResonancesTuple[3]))   # UNECESSARY CAUSE WE ALREADY SAVED THEM IN HCNH_residue_NHresonances_dict
                        #N_resonances_list.append(float(ResonancesTuple[2]))

                        #print "\t"+str(atom_index)+"\t"+str(Carbon_resonance)+"\t0.2\t"+str(Carbon_name)+"\t"+str(residue[3:])
                        Creson_Cname_lines_list.append([float(Carbon_resonance), Carbon_name])

                        #print "\t"+str(atom_index)+"\t"+str(Hydrogen_resonance)+"\t0.02\t"+str(Hydrogen_name)+"\t"+str(residue[3:])
                        #xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s\n" % (atom_index, float(Hydrogen_resonance), 0.02, Hydrogen_name, int(residue[3:]), aa_type))
                        xeasy_lines_list.append( [atom_index, float(Hydrogen_resonance), 0.02, Hydrogen_name, resid, aa_type, "HCNH (flanked residue) iter"+str(iteration)] )  # append the Hydrogen resonance line
                        print("DEBUG: point 3 appending iter: ", [atom_index, float(Hydrogen_resonance), 0.02, Hydrogen_name, resid, aa_type, "HCNH (flanked residue) iter"+str(iteration)])
                        atom_index += 1
                        if N_resonance != None and HN_resonance != None:    # if no N-H was found (e.g. N-term residue) don't write the i-i Sparky lines for this residue, only the xeasy lines
                            try:
                                #sparky_fout.write("\t%s-%s-%s-H\t\t%.3f\t%.3f\t%.3f\t%.3f\t\n" % ( aa3to1_dict[aa_type] + residue[3:] + Hydrogen_name, Carbon_name, protein_alignment_list[residue2position_dict[residue]] + str(int(residue[3:])+1) + "N",float(Hydrogen_resonance), float(Carbon_resonance), float(N_resonance), float(HN_resonance)) )
                                sparky_lines_list.append( [i_residue + Hydrogen_name, Carbon_name, i_residue + "N",float(Hydrogen_resonance), float(Carbon_resonance), float(N_resonance), float(HN_resonance)] )
                            except IndexError:
                                #sparky_fout.write("\t%s-%s-%s-H\t\t%.3f\t%.3f\t%.3f\t%.3f\t\n" % ( aa3to1_dict[aa_type] + residue[3:] + Hydrogen_name, Carbon_name, str(int(residue[3:])+1) + "N", float(Hydrogen_resonance), float(Carbon_resonance), float(N_resonance), float(HN_resonance)) )
                                sparky_lines_list.append( [i_residue + Hydrogen_name, Carbon_name, i_residue + "N", float(Hydrogen_resonance), float(Carbon_resonance), float(N_resonance), float(HN_resonance)] )


                    Cnames_list = [cl[1] for cl in Creson_Cname_lines_list]
                    Cnames_set = set(Cnames_list)
                    for Carbon_name in Cnames_set:
                        if Cnames_list.count(Carbon_name) > 1:  # if this carbon occurs in multiple C-H pairs, calculate its average resonance for the XEASY format
                            Carbon_resonance = round(np.average([cl[0] for cl in Creson_Cname_lines_list if cl[1]==Carbon_name]), 3)
                        else:
                            Carbon_resonance = round(np.average([cl[0] for cl in Creson_Cname_lines_list if cl[1]==Carbon_name]), 3)
                        #xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s\n" % (atom_index, float(Carbon_resonance), 0.2, Carbon_name, int(residue[3:]), aa_type))
                        xeasy_lines_list.append( [atom_index, float(Carbon_resonance), 0.2, Carbon_name, resid, aa_type, "HCNH (flanked residue) iter"+str(iteration)] )
                        print("DEBUG: point 4 appending iter: ", [atom_index, float(Carbon_resonance), 0.2, Carbon_name, resid, aa_type, "HCNH (flanked residue) iter"+str(iteration)])
                        atom_index += 1

                    # RENAME DEGENERATE PROTONS
                    new_xeasy_lines_list, new_sparky_lines_list = rename_protons(xeasy_lines_list, sparky_lines_list)
                    print("DEBUG: point 2 new_xeasy_lines_list=", new_xeasy_lines_list)
                    # WRITE C AND H LINES IN XEASY FORMAT
                    for xeasy_line in new_xeasy_lines_list:
                        xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s %s\n" % (xeasy_line[0], xeasy_line[1], xeasy_line[2], xeasy_line[3], xeasy_line[4], xeasy_line[5], xeasy_line[6]) )
                    # WRITE C AND H LINES IN SPARKY FORMAT
                    for sparky_line in new_sparky_lines_list:
                        sparky_fout.write("\t%s-%s-%s-H\t\t%.3f\t%.3f\t%.3f\t%.3f\t\n" % (sparky_line[0], sparky_line[1], sparky_line[2], sparky_line[3], sparky_line[4], sparky_line[5], sparky_line[6]))
                    # WRITE N AND HN LINES IN XEASY FORMAT
                    Nreson, HNreson, Nstdev, HNstdev = HCNH_residue_NHresonances_dict[i_residue]
                    xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s %s\n" % (atom_index, HNreson, HNstdev, "H", resid, aa_type, 'HCNH (flanked residue)' ) )
                    atom_index += 1
                    xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s %s\n" % (atom_index, Nreson, Nstdev, "N", resid, aa_type, 'HCNH (flanked residue)' ) )
                    atom_index += 1
                    residues_with_written_NH_list.append(i_residue)
                    matched_HCNH_residue_assignments_dict = update_HCNH_peaks(i_residue, matched_HCNH_residue_assignments_dict, protein_alignment_list, absolute_matches_alignment_list, args, new_sparky_lines_list)
            matched_HCNH_residue_assignments_dict = \
                remove_iplus1_labels(i_residue,
                                     matched_HCNH_residue_assignments_dict,
                                     protein_alignment_list,
                                     absolute_matches_alignment_list,
                                     args)
        elif len(missing_carbons_list) == 0 and len(partly_assigned_carbons_list) == 0:   # IF NO CARBON IS MISSING, SKIP I_RESIDUE i_residue
            update_flanked_HCNH_peaks(i_residue, matched_HCNH_residue_assignments_dict, protein_alignment_list,
                               absolute_matches_alignment_list, args)
            matched_HCNH_residue_assignments_dict = \
                remove_iplus1_labels(i_residue,
                                     matched_HCNH_residue_assignments_dict,
                                     protein_alignment_list,
                                     absolute_matches_alignment_list,
                                     args)
            residues_with_written_prevAssigned_peaks.append(i_residue)
            # WRITE ALSO THE PEAKS ASSIGNED IN ITERATION 1 & 2
            write_peaks_assigned_in_previous_iteration(i_residue, clean_resid_assignments_dict, atom_index, xeasy_fout,
                                               HCNH_residue_NHresonances_dict, residues_with_written_NH_list)
            continue
        elif len(missing_carbons_list) == 0 and len(partly_assigned_carbons_list) > 0:   # IF ALL THE CARBONS WERE ASSIGNED, FIND & SAVE PARTLY ASSIGNED CARBONS, SAVE THE MATCHED HCNH RESONANCES, AND IF NOT MATCHED, SAVE THE TOCSY RESONANCES
            # TODO
            print("ERROR: missing code to find & save partly assigned carbons of flanked residues.")
            sys.exit(1)

        matched_HCNH_residue_assignments_dict = \
            remove_iplus1_labels(i_residue,
                                 matched_HCNH_residue_assignments_dict,
                                 protein_alignment_list,
                                 absolute_matches_alignment_list,
                                 args)     # remove all the 'i+1', 'i-1' labels from this residue and residue i+1 (if applicable)

        # WRITE ALSO THE PEAKS ASSIGNED IN PREVIOUS ITERATIONS OF ALL RESIDUES THAT WERE EXCLUDED FROM THIS ITERATION
        for i_aa_resid in revordered_residue_keys:
            i_residue = i_aa_resid[0]+str(i_aa_resid[1])
            aa_type = aa1to3_dict[i_aa_resid[0]]
            if not i_residue in residues_with_written_prevAssigned_peaks:
                write_peaks_assigned_in_previous_iteration(i_residue, clean_resid_assignments_dict, atom_index, xeasy_fout,
                                               HCNH_residue_NHresonances_dict, residues_with_written_NH_list)

        matched_HCNH_residue_assignments_dict = \
            remove_iplus1_labels(i_residue,
                                 matched_HCNH_residue_assignments_dict,
                                 protein_alignment_list,
                                 absolute_matches_alignment_list,
                                 args)

    return revordered_flanked_residues


def average_methylene_Carbons(peak_list):
    """
        FUNCTION to average methylene Carbon resonance for the xeasy format.
        E.g.
        peak_list = [('Q3', 'CA', 58.802, 'HA', 3.98), ('Q3', 'CB', 29.152, 'HB2', 2.265), ('Q3', 'CB', 29.159, 'HB3', 2.046),
        ('Q3', 'CG', 34.548, 'HG2', 2.265), ('Q3', 'CG', 34.563, 'HG3', 2.438)]
        xeasy_CarbonResonances = [('CB', 29.155), ('CG', 34.556), ('CA', 58.802)]
        xeasy_ProtonResonances = [('HA', 3.98), ('HB2', 2.265), ('HB3', 2.046), ('HG2', 2.265), ('HG3', 2.438)]
    """
    Carbon_list = [p[1] for p in peak_list]
    methylene_Carbons = set([c for c in Carbon_list if Carbon_list.count(c)>1])
    xeasy_CarbonResonances = []   # atom names and their resonances, but methyle resonances are averaged
    xeasy_ProtonResonances = []   # atom names and their resonances, but methyle resonances are averaged
    for Cname in methylene_Carbons:
        aveCreson = round(np.mean([p[2] for p in peak_list if p[1]==Cname]), 3)
        xeasy_CarbonResonances.append((Cname, aveCreson))
    # add the rest of the non-methylene Carbons and all protons
    for p in peak_list:
        if not p[1] in methylene_Carbons:   # add non-methylene Carbons
            xeasy_CarbonResonances.append((p[1], p[2]))
        xeasy_ProtonResonances.append((p[3], p[4])) # add all protons (methylene or not)

    return xeasy_CarbonResonances, xeasy_ProtonResonances


def write_unmatched_TOCSY_peaks(matched_HCNH_residue_assignments_dict,
                                TOCSY_residue_assignments_dict,
                                atom_index,
                                residues_with_written_NH_list,
                                xeasy_fout,
                                sparky_fout,
                                protein_alignment_list,
                                absolute_matches_alignment_list,
                                args,
                                HCNH_residue_NHresonances_dict,
                                TOCSY_residue_NHresonances_dict,
                                iteration=1):
    ##
    ##  COPY FROM TOCSY ASSIGNMENTS (TOCSY_residue_assignments_dict) ALL THE RESONANCES THAT COULD NOT BE MATCHED WITH THE HCNH ASSIGNMENTS
    ## (TOCSY_residue_assignments_dict)
    ##
    print("DEBUG write_unmatched_TOCSY_peaks: matched_HCNH_residue_assignments_dict=", matched_HCNH_residue_assignments_dict)
    for i_residue in list(TOCSY_residue_assignments_dict.keys()):
        i_aa_type = aa1to3_dict[get_aa_type_from_residue(i_residue)]
        # resid = absolute_matches_alignment_list.index(i_residue) + int(args.FIRST_RESIDUE_NUMBER) # OLD WAY; returns ValueError for A335
        resid = get_resid_from_residue(i_residue)
        if i_residue in list(HCNH_residue_NHresonances_dict.keys()):
            Nreson = HCNH_residue_NHresonances_dict[i_residue][0]
            HNreson = HCNH_residue_NHresonances_dict[i_residue][1]
            Nstdev = HCNH_residue_NHresonances_dict[i_residue][2]
            HNstdev = HCNH_residue_NHresonances_dict[i_residue][3]
        elif i_residue in list(TOCSY_residue_NHresonances_dict.keys()):
            Nreson = TOCSY_residue_NHresonances_dict[i_residue][0]
            HNreson = TOCSY_residue_NHresonances_dict[i_residue][1]
            Nstdev = TOCSY_residue_NHresonances_dict[i_residue][2]
            HNstdev = TOCSY_residue_NHresonances_dict[i_residue][3]
        else:
            HNreson = -1000
            Nreson = -1000
        if not i_residue in list(matched_HCNH_residue_assignments_dict.keys()):  # if this residue has not HCNH (e.g. Proline)
            print("Residue", i_residue, " has not HCNH peaks. Adding it's assigned TOCSY peaks:")
            print(TOCSY_residue_assignments_dict[i_residue])
            xeasy_CarbonResonances, xeasy_ProtonResonances = average_methylene_Carbons(TOCSY_residue_assignments_dict[i_residue])
            for Cname,Creson in xeasy_CarbonResonances:
                # WRITE C AND H LINES IN XEASY FORMAT
                xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s %s\n" % (atom_index, Creson, 0.2, Cname, resid, i_aa_type, 'TOCSY (unmatched)' ) )
                atom_index += 1
            for Hname,Hreson in xeasy_ProtonResonances:
                xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s %s\n" % (atom_index, Hreson, 0.02, Hname, resid, i_aa_type, 'TOCSY (unmatched)' ) )
                atom_index += 1
            # WRITE C AND H LINES IN SPARKY FORMAT
            for TOCSY_peak in TOCSY_residue_assignments_dict[i_residue]:
                if not -1000 in [HNreson, Nreson]:
                    sparky_fout.write("\t%s-%s-%s-H\t\t%.3f\t%.3f\t%.3f\t%.3f\t\n" % (i_residue+TOCSY_peak[3], TOCSY_peak[1], i_residue+"N", TOCSY_peak[4], TOCSY_peak[2], Nreson, HNreson ))
                # update_HCNH_peaks(i_residue, new_sparky_lines_list)  # these are TOCSY resonances, therefore don't update the matched_HCNH_residue_assignments_dict
            # WRITE HN AND N LINES IN XEASY FORMAT
            if not -1000 in [HNreson, Nreson]:
                # because the keys of matched_HCNH_residue_assignments_dict and HCNH_residue_NHresonances_dict are the same, there are no HCNH HN & N resonances,
                # therefore we write only TOCSY HN & N resonances
                xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s %s\n" % (atom_index, HNreson, HNstdev, "H", resid, i_aa_type, 'TOCSY (unmatched)' ) )
                atom_index += 1
                xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s %s\n" % (atom_index, Nreson, Nstdev, "N", resid, i_aa_type, 'TOCSY (unmatched)' ) )
                atom_index += 1
                residues_with_written_NH_list.append(i_residue)
        else:   # if this residue has HCNH, find and write the unmatched peaks from TOCSY
            unmatched_TOCSY_peaks = []
            for TOCSY_peak in TOCSY_residue_assignments_dict[i_residue]:
                TOCSY_residue = TOCSY_peak[0]
                TOCSY_Cname = TOCSY_peak[1]
                TOCSY_Creson = TOCSY_peak[2]
                TOCSY_Hname = TOCSY_peak[3]
                TOCSY_Hreson = TOCSY_peak[4]
                CH_pair_found = False
                for HCNH_peak in matched_HCNH_residue_assignments_dict[i_residue]:
                    HCNH_residue = TOCSY_peak[0]
                    HCNH_Cname = TOCSY_peak[1]
                    HCNH_Creson = TOCSY_peak[2]
                    HCNH_Hname = TOCSY_peak[3]
                    HCNH_Hreson = TOCSY_peak[4]
                    if TOCSY_residue == HCNH_residue and TOCSY_Cname == HCNH_Cname and TOCSY_Hname == HCNH_Hname:
                        CH_pair_found = True
                        break
                if CH_pair_found == False:
                    print("Adding missing CH resonances of residue ", TOCSY_residue, " from assigned TOCSY file:")
                    print(TOCSY_peak)
                    unmatched_TOCSY_peaks.append(TOCSY_peak)
            xeasy_CarbonResonances, xeasy_ProtonResonances = average_methylene_Carbons(unmatched_TOCSY_peaks)
            for Cname,Creson in xeasy_CarbonResonances:
                # WRITE C AND H LINES IN XEASY FORMAT
                xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s %s\n" % (atom_index, Creson, 0.2, Cname, resid, i_aa_type, 'TOCSY (unmatched)' ) )
                atom_index += 1
            for Hname,Hreson in xeasy_ProtonResonances:
                xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s %s\n" % (atom_index, Hreson, 0.02, Hname, resid, i_aa_type, 'TOCSY (unmatched)' ) )
                atom_index += 1
            # WRITE C AND H LINES IN SPARKY FORMAT
            for TOCSY_peak in unmatched_TOCSY_peaks:
                if not -1000 in [HNreson, Nreson]:
                    sparky_fout.write("\t%s-%s-%s-H\t\t%.3f\t%.3f\t%.3f\t%.3f\t\n" % (i_residue+TOCSY_peak[3], TOCSY_peak[1], i_residue+"N", TOCSY_peak[4], TOCSY_peak[2], Nreson, HNreson ))
                    # update_HCNH_peaks(i_residue, new_sparky_lines_list) # these are TOCSY resonances, therefore don't update the matched_HCNH_residue_assignments_dict
            # WRITE HN AND N LINES IN XEASY FORMAT
            if not -1000 in [HNreson, Nreson]:
                # because the keys of matched_HCNH_residue_assignments_dict and HCNH_residue_NHresonances_dict are the same, there are HCNH HN & N resonances,
                # therefore we write them into the xeasy file
                xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s\n" % (atom_index, HNreson, HNstdev, "H", resid, i_aa_type ) )
                atom_index += 1
                xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s\n" % (atom_index, Nreson, Nstdev, "N", resid, i_aa_type ) )
                atom_index += 1
                residues_with_written_NH_list.append(i_residue)


def write_NH_of_residues_without_CH(absolute_matches_alignment_list,
                                    HCNH_residue_NHresonances_dict,
                                    TOCSY_residue_NHresonances_dict,
                                    residues_with_written_NH_list,
                                    atom_index,
                                    xeasy_fout,
                                    sparky_fout,
                                    protein_alignment_list,
                                    args,
                                    Cterm_residue_set,
                                    iteration=1):


    ##
    ##  WRITE THE HN & N RESONANCES TO XEASY ONLY FOR THOSE RESIDUES THAT NO C,H ASSIGNMENTS COULD BE MADE BUT
    ##  THEIR HN & H RESONANCES ARE KNOWN FROM TOCSY OR HCNH
    ##
    for residue in HCNH_residue_NHresonances_dict.keys():
        if residue not in residues_with_written_NH_list and residue in Cterm_residue_set:
            try:
                i_aa_type = aa1to3_dict[get_aa_type_from_residue(residue)]
                resid = get_resid_from_residue(residue)
            except ValueError:  # if this residue is not  in absolute_matches_alignment_list, skip it
                continue
            HNreson = HCNH_residue_NHresonances_dict[residue][1]
            Nreson = HCNH_residue_NHresonances_dict[residue][0]
            HNstdev = HCNH_residue_NHresonances_dict[residue][2]
            Nstdev = HCNH_residue_NHresonances_dict[residue][3]
            spectrum_combo = 'HCNH (residue with no C-H assignments)'
            if residue in Cterm_residue_set:
                spectrum_combo = 'HCNH (C-term hanging residue)'
            xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s %s\n" % (atom_index, HNreson, HNstdev, "H", resid, i_aa_type, spectrum_combo ) )
            atom_index += 1
            xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s %s\n" % (atom_index, Nreson, Nstdev, "N", resid, i_aa_type, spectrum_combo ) )
            atom_index += 1
            residues_with_written_NH_list.append(residue)

    for residue in TOCSY_residue_NHresonances_dict.keys():
        if residue not in residues_with_written_NH_list and residue in Cterm_residue_set:
            try:
                i_aa_type = aa1to3_dict[get_aa_type_from_residue(residue)]
                resid = get_resid_from_residue(residue)
            except ValueError:  # if this residue is not  in absolute_matches_alignment_list, skip it
                continue
            HNreson = TOCSY_residue_NHresonances_dict[residue][1]
            Nreson = TOCSY_residue_NHresonances_dict[residue][0]
            HNstdev = TOCSY_residue_NHresonances_dict[residue][2]
            Nstdev = TOCSY_residue_NHresonances_dict[residue][3]
            spectrum_combo = 'TOCSY (residue with no C-H assignments)'
            if residue in Cterm_residue_set:
                spectrum_combo = 'TOCSY (C-term hanging residue)'
            xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s %s\n" % (atom_index, HNreson, HNstdev, "H", resid, i_aa_type, spectrum_combo) )
            atom_index += 1
            xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s %s\n" % (atom_index, Nreson, Nstdev, "N", resid, i_aa_type, spectrum_combo) )
            atom_index += 1
            residues_with_written_NH_list.append(residue)


def rename_methylene_protons(raw_resid_assignments_dict):
    """
        RENAME THE METHYLENE PROTONS SO THAT THE LOWEST RESONANCE IS NAMED ALWAYS H*2 AND THE HIGHEST H*3
    """
    for resid in list(raw_resid_assignments_dict.keys()):
        resname = raw_resid_assignments_dict[resid][0][6]
        if not resname in list(aatype_carbon_nondegenerateHlist_mdict.keys()):
            continue
        atom_names_set = set([p[3] for p in raw_resid_assignments_dict[resid]])
        carbons_set = set([p[3] for p in raw_resid_assignments_dict[resid] if p[3][0]=='C'])
        for carbon in carbons_set:
            if carbon in list(aatype_carbon_nondegenerateHlist_mdict[resname].keys()):
                try:
                    H2, H3 = aatype_carbon_nondegenerateHlist_mdict[resname][carbon]
                except ValueError as e:
                    print(e)
                    print("resid=", resid, "raw_resid_assignments_dict[resid]", raw_resid_assignments_dict[resid])
                    print("carbon=", carbon, "resname=", resname)
                    print("aatype_carbon_nondegenerateHlist_mdict[resname][carbon]=", aatype_carbon_nondegenerateHlist_mdict[resname][carbon])
                    raise ValueError
                except Exception as e:
                    print(e)
                    raise Exception
                if H2 in atom_names_set and H3 in atom_names_set:   # if both methylene protons were assigned
                    reson_set = set([p[1] for p in raw_resid_assignments_dict[resid] if p[3] in [H2, H3]])
                    if len(reson_set) == 1:
                        print("ERROR: two methylene protons with the same resonance detected (", H2, H3, ")!")
                        print("DEBUG: resid=", resid, "resname=", resname, "raw_resid_assignments_dict[resid]=", raw_resid_assignments_dict[resid])
                        raise Exception
                    else:   # if there are two different resonances, name the lowest H*2 and the highest H*1
                        reson_list = list(reson_set)
                        reson_list.sort()   # lowest resonance goes 1st
                        for p in raw_resid_assignments_dict[resid]:
                            if p[1] == reson_list[0] and p[3] in [H2,H3]:
                                p[3] = H2
                            elif p[1] == reson_list[1] and p[3] in [H2,H3]:
                                p[3] = H3

    return raw_resid_assignments_dict


def clean_xeasy_file(xeasy_fname):
    """
    Method to REWRITE THE XEASY FILE BY REMOVING DUPLICATES OR AVERAGING, AND CORRECTING THE ORDER
    :param xeasy_fname:
    :return:
    """

    def save_peak(peak, atom_index):
        peak.insert(0, atom_index)  # insert -1 in the atom index position to designate that this assignment was corrected
        try:    # save only the 1st since all of them are the same
            clean_resid_assignments_dict[resid].append(peak)
        except KeyError:
            clean_resid_assignments_dict[resid] = [peak]

    xeasy_lines_list = []
    with open(xeasy_fname, 'r') as f:
        for line in f:
            word_list = line.split()
            xeasy_lines_list.append( [int(word_list[0]), float(word_list[1]), float(word_list[2]), word_list[3], int(word_list[4]), word_list[5], word_list[6], " ".join(word_list[7:])] )

    xeasy_lines_list.sort(key=itemgetter(4, 3))
    print("DEBUG: sorted xeasy_lines_list=", xeasy_lines_list)
    raw_resid_assignments_dict = {}
    for line in xeasy_lines_list:
        resid = line[4]
        try:
            raw_resid_assignments_dict[resid].append(line)
        except KeyError:
            raw_resid_assignments_dict[resid] = [line]
    # RENAME THE METHYLENE PROTONS SO THAT THE LOWEST RESONANCE IS NAMED ALWAYS H*2 AND THE HIGHEST H*3
    resid_assignments_dict = rename_methylene_protons(raw_resid_assignments_dict)
    # print  "DEBUG: resid_assignments_dict=", resid_assignments_dict
    # REMOVE REDUNDANT DUPLICATE RESONANCES
    clean_resid_assignments_dict = {}   # same as resid_assignments_dict but the duplicates have been removed or averaged
    for resid in list(resid_assignments_dict.keys()):
        print("DEBUG: resid=", resid, "resid_assignments_dict[resid]=", resid_assignments_dict[resid])
        resname = resid_assignments_dict[resid][0][6]
        atom_names_set = set([p[3] for p in resid_assignments_dict[resid]])
        for atom_name in atom_names_set:
            # remove identical peaks (ignore the atom index)
            matched_peaks = set([tuple(p[1:]) for p in resid_assignments_dict[resid] if p[3]==atom_name])   # WARNING: omitting the atom index
            unique_matched_peaks = [list(p) for p in matched_peaks]    # convert the set back to list of lists
            print("DEBUG: resid", resid, "atom_name=", atom_name, "unique_matched_peaks=", unique_matched_peaks)
            matched_atom_indices = [p[0] for p in resid_assignments_dict[resid] if p[3]==atom_name]
            if len(unique_matched_peaks) == 1:
                reson = unique_matched_peaks[0][0]
                name = unique_matched_peaks[0][2]
                spectrum_combo = unique_matched_peaks[0][6]
                if name[0] == 'Q' and spectrum_combo == 'TOCSY-HCNH':
                    identical_proton_peak_list = [p for p in resid_assignments_dict[resid] if p[1]==reson and p[3][0]=='H' and p[3][1]==name[1] and p[7]==spectrum_combo]
                    if len(identical_proton_peak_list)==0: # save this unique resonance, otherwise avoid double writing TOCSY-HCNH methylene protons (e.g. QB & HB2)
                        save_peak(unique_matched_peaks[0], matched_atom_indices[0])
                    else:
                        print("DEBUG: not saving xeasy_line=", unique_matched_peaks[0])
                        continue
                else:   # save this unique resonance
                    print("DEBUG: saving peak:", unique_matched_peaks[0], matched_atom_indices[0], "from file:", xeasy_fname)
                    save_peak(unique_matched_peaks[0], matched_atom_indices[0])
            else:
                matched_resonances = [(p[0],p[1]) for p in unique_matched_peaks]
                if len(set(matched_resonances)) == 1:   # these are all the same (duplicates, triplicates, etc.)
                    save_peak(unique_matched_peaks[0], 0)
                else:   # in this case we have different resonances for the same nucleous
                    # CASE 1: only one of them comes from HCNH, then keep that one
                    HCNH_peaks = [p for p in unique_matched_peaks if not 'TOCSY' in p[6] or p[6]=='TOCSY-HCNH']
                    if len(HCNH_peaks) == 1:
                        print("WARNING: only one of the duplicate peaks comes from HCNH. Saving that one. unique_matched_peaks=", unique_matched_peaks)
                        save_peak(HCNH_peaks[0], -1)
                    # CASE 2: none of them comes only from HCNH (namely they are 'TOCSY')
                    elif len(HCNH_peaks) == 0:
                        print("WARNING: none of the duplicate resonances comes from HCNH. They will be averaged. unique_matched_peaks=", unique_matched_peaks)
                        ave_reson = round(np.average([p[0] for p in unique_matched_peaks]), 3)
                        ave_stdev = round(np.average([p[1] for p in unique_matched_peaks]), 3)
                        average_peak = [ave_reson, ave_stdev, atom_name, resid, '#', resname, 'TOCSY (average)']
                        save_peak(average_peak, -2)
                    # CASE 3: more that one comes from HCNH (namely none is 'TOCSY')
                    elif len(HCNH_peaks) > 1:
                        print("WARNING: none of the duplicate resonances comes from TOCSY. They will be averaged. unique_matched_peaks=", unique_matched_peaks)
                        ave_reson = round(np.average([p[0] for p in unique_matched_peaks]), 3)
                        ave_stdev = round(np.average([p[1] for p in unique_matched_peaks]), 3)
                        average_peak = [ave_reson, ave_stdev, atom_name, resid, '#', resname, 'HCNH (average)']
                        save_peak(average_peak, -3)

    return clean_resid_assignments_dict


def match_HCNH_i_iplus1_peaks(i_residue,
                               absolute_matches_alignment_list,
                               matched_HCNH_residue_assignments_dict,
                               matched_i_to_iplus1_peaks_mdict):
    """
    Method to identify the common HCNH peaks between residues i and i+1 and to update the
    matched_HCNH_residue_assignments_dict and matched_i_to_iplus1_peaks_mdict dictionaries
    every time for 2 residues only (residue i and residue i+1). These peaks in residue i are named
    'i+1" whereas in residue i+1 are named 'i-1'.

    :param i_residue:
    :param absolute_matches_alignment_list:
    :param matched_HCNH_residue_assignments_dict:
    :param matched_i_to_iplus1_peaks_mdict:
    :return False:  if the i_residue is in the C-term or after that there is a gap in the
        absolue_matches_alignment_list. Otherwise returns True.
    :return matched_HCNH_residue_assignments_dict: updated
    :return matched_i_to_iplus1_peaks_mdict:    updated
    """

    tolH = 0.02
    tolC = 0.2

    if i_residue not in matched_HCNH_residue_assignments_dict.keys():  # if this residues does not have HCNH peaks
        return False, matched_HCNH_residue_assignments_dict, matched_i_to_iplus1_peaks_mdict

    try:
        i_index = absolute_matches_alignment_list.index(i_residue)
    except ValueError:  # if i_residue is not in the alignment it means it doesn't have HCNH (?), e.g. PRO,
                        # thus skip it. Recall that this function is called in the 1st iteration
                        # only. In the next interations the residues that are not in the alignment
                        # will be assigned, e.g. PRO from the TOCSY of residue i+1.
        return False, matched_HCNH_residue_assignments_dict, matched_i_to_iplus1_peaks_mdict

    try:
        iplus1_residue = absolute_matches_alignment_list[i_index+1]
    except IndexError:  # If i_residue is the last residue
        return False, matched_HCNH_residue_assignments_dict, matched_i_to_iplus1_peaks_mdict
    if '-' in [i_residue, iplus1_residue] or 'N/A' in [i_residue, iplus1_residue] or not iplus1_residue in list(matched_HCNH_residue_assignments_dict.keys()):  # if there is a gap in the 1 or i+1 position or if this is the C-term residue (e.g. NAB E117) you cannot match peaks
        return False, matched_HCNH_residue_assignments_dict, matched_i_to_iplus1_peaks_mdict

    for i_peak_index in range(len(matched_HCNH_residue_assignments_dict[i_residue])):
        # e.g. i_peak = ['K173', 'CE', 42.178, 'QE', 3.042]
        i_peak = matched_HCNH_residue_assignments_dict[i_residue][i_peak_index]
        if i_peak[0] != '?':    # if this HCNH peak has been matched with a TOCSY peak, skip it
            continue
        for iplus1_peak_index in range(len(matched_HCNH_residue_assignments_dict[iplus1_residue])):
            iplus1_peak = matched_HCNH_residue_assignments_dict[iplus1_residue][iplus1_peak_index]
            print("DEBUG: comparing i_peak=", i_peak, " with iplus1_peak=", iplus1_peak)
            if iplus1_peak[0] != '?':    # similarly if this HCNH peak has been matched with a TOCSY peak, skip it
                continue
            elif approx_equal(i_peak[2], iplus1_peak[2], tolC) and approx_equal(i_peak[4], iplus1_peak[4], tolH):
                i_peak[0] = 'i+1'
                matched_HCNH_residue_assignments_dict[i_residue][i_peak_index] = i_peak
                iplus1_peak[0] = 'i-1'
                matched_HCNH_residue_assignments_dict[iplus1_residue][iplus1_peak_index] = iplus1_peak
                matched_i_to_iplus1_peaks_mdict[i_residue][(float(i_peak[2]), float(i_peak[4]))] = (float(iplus1_peak[2]), float(iplus1_peak[4]))
                print("DEBUG: matched i->i+1 HCNH peaks of residues", i_residue," and ", iplus1_residue, "i_peak=", i_peak, "iplus1_peak=", iplus1_peak)
    ## TODO: keep only the closest i+1 peak in matched_i_to_iplus1_peaks_mdict
    # print "DEBUG: matched_i_to_iplus1_peaks_mdict=", matched_i_to_iplus1_peaks_mdict
    return True, matched_HCNH_residue_assignments_dict, matched_i_to_iplus1_peaks_mdict


def get_missing_carbons_from_xeasy_dict(aa_type, resid, clean_resid_assignments_dict, matched_HCNH_residue_assignments_dict, allowed_aa_atoms_dict):
    """
        FUNCTION to find the Carbons of residue i that were not assigned at all in the previous iteration and those methylene Carbons for
        which only one proton has been assigned. Apart from the cleaned xeasy file from the previous iteration this function also checks if
        the missing carbons have been assigned in the matched_HCNH_residue_assignments_dict, which happens for instance in ITERATION 3 where
        a peak may have been assigned in one of the I_THRESHOLD iterations but not have been written in the xeasy file yet.

        ARGS:
        clean_resid_assignments_dict:    the contents of cleaned xeasy file in dict format

        RETURNS:
        missing_carbons_list:            list with completely unassgined carbons (each carbon is contained only once)
        partly_assigned_carbons_list:    list with partially assigned methylene carbons (does not overlap with missing_carbons_list)
    """

    if not type(resid) == int:
        print("DEBUG: ERROR resid must be int!!!")
        sys.exit(1)

    assigned_carbons_list = []  # list of carbons which have ALL their protons assigned in iteration 1
    partly_assigned_carbons_list = []   # list of methylene carbons for which only one proton frequency was found in iteration 1
    i_residue = aa3to1_dict[aa_type] + str(resid)
    HCNH_assigned_carbons_list = [x[1] for x in matched_HCNH_residue_assignments_dict[i_residue] if x[0]==i_residue]   # these are the most current assignments of residue i in HCNH dict
    all_allowed_carbons_list = [x for x in allowed_aa_atoms_dict[aa_type] if x[0] == "C"]   # all allowed carbons for this aa type
    print("DEBUG get_missing_carbons_from_xeasy_dict: HCNH_assigned_carbons_list=", HCNH_assigned_carbons_list)
    print("DEBUG get_missing_carbons_from_xeasy_dict: all_allowed_carbons_list=", all_allowed_carbons_list)
    if not resid in list(clean_resid_assignments_dict.keys()):    # If this is a C-terminal residue, it would have been skiped in iteration1
        return all_allowed_carbons_list, []     # recall C-term have HCNH but no TOCSY, therefore all Carbons are missing

    print("DEBUG get_missing_carbons_from_xeasy_dict: resid=", resid, "clean_resid_assignments_dict[resid]=", clean_resid_assignments_dict[resid])
    tmp_assigned_carbons_list = [x[3] for x in clean_resid_assignments_dict[resid] if x[3][0]=='C']   # these are the assigned Carbons of residue i in iteration 1
    for Cname in tmp_assigned_carbons_list:
        if Cname in list(aatype_carbon_nondegenerateHlist_mdict[aa_type].keys()):
            assigned_ethylene_protons = [x[3] for x in clean_resid_assignments_dict[resid] if x[3] in aatype_carbon_nondegenerateHlist_mdict[aa_type][Cname]]  # find the ethylen protons that were assigned for this Carbon in iteration 1
            if len(assigned_ethylene_protons) == 0 and HCNH_assigned_carbons_list.count(Cname) == 1:     # if only 1 was assigned, it was named QH* and hence was not found in aatype_carbon_nondegenerateHlist_mdict
                partly_assigned_carbons_list.append(Cname)
            elif len(assigned_ethylene_protons) == 2 and HCNH_assigned_carbons_list.count(Cname) == 2:
                assigned_carbons_list.append(Cname)
            elif len(assigned_ethylene_protons) == 2 and 'TOCSY (average)' in [c for c in clean_resid_assignments_dict[resid] if c[3]==Cname][0][7]:
                # 'TOCSY (average)' -> 2 resonances available only from TOCSY and they were averaged, therefore HCNH_assigned_carbons_list.count(Cname) == 0
                assigned_carbons_list.append(Cname)
        else:   # if this Carbon is not methylene, save it as assigned
            assigned_carbons_list.append(Cname)
    print("DEBUG get_missing_carbons_from_xeasy_dict: assigned_carbons_list=", assigned_carbons_list)
    print("DEBUG get_missing_carbons_from_xeasy_dict: partly_assigned_carbons_list=", partly_assigned_carbons_list)
    missing_carbons_list = []
    for carbon in all_allowed_carbons_list:     # keep only the carbon types of this aa type that were not assigned to peaks
        if not carbon in assigned_carbons_list and not carbon in partly_assigned_carbons_list:
            if not carbon in HCNH_assigned_carbons_list:   # if this carbon was not assigned during one of I_THRESHOLD iterations
                missing_carbons_list.append(carbon)
    print("DEBUG get_missing_carbons_from_xeasy_dict: missing_carbons_list=", missing_carbons_list)


    return missing_carbons_list, partly_assigned_carbons_list


def get_peaks_from_xeasy(resid, clean_resid_assignments_dict):
    """
        FUNCTION to convert xeasy format lines to peak format: [residue, Cname, Creson, atom_name, reson, spectrum_comment]
    """

    peak_list = []
    Cname_CresonHresonTupleList_dict = {}
    methylene_carbons_list = []
    non_methylene_carbons_list = []
    for line in clean_resid_assignments_dict[resid]:
        residue = aa3to1_dict[line[6]] + str(resid)
        atom_name = line[3]
        reson = line[1]
        aa_type = line[6]
        # CHECK IF THIS IS A METHYLENE CARBON
        if atom_name[0] == 'C' and atom_name in list(aatype_carbon_nondegenerateHlist_mdict[aa_type].keys()):
            methylene_carbons_list.append((atom_name, reson))
        elif atom_name[0] == 'C':
            non_methylene_carbons_list.append((atom_name, reson))
    print("DEBUG get_peaks_from_xeasy: resid=", resid, "methylene_carbons_list=", methylene_carbons_list, "non_methylene_carbons_list=", non_methylene_carbons_list)
    # SAVE METHYLENE CARBON PEAKS
    for duplet in methylene_carbons_list:
        Cname = duplet[0]
        Creson = duplet[1]
        protons_found = 0
        for line in clean_resid_assignments_dict[resid]:
            atom_name = line[3]
            reson = line[1]
            spectrum_comment = line[7]  # the spectrum type, but may also have the iteration number
            if atom_name in aatype_carbon_nondegenerateHlist_mdict[aa_type][Cname]:
                protons_found += 1
                peak_list.append([residue, Cname, Creson, atom_name, reson, spectrum_comment])
            elif atom_name in aatype_carbon_degenerateH_mdict[aa_type][Cname]:  # in this case there is onle one methylene proton assigned, so break the loop
                peak_list.append([residue, Cname, Creson, atom_name, reson, spectrum_comment])
                break
    # SAVE NON-METHYLENE CARBON PEAKS
    for duplet in non_methylene_carbons_list:
        Cname = duplet[0]
        Creson = duplet[1]
        for line in clean_resid_assignments_dict[resid]:
            atom_name = line[3]
            reson = line[1]
            spectrum_comment = line[7]  # the spectrum type, but may also have the iteration number
            if ( (aa_type in list(aatype_carbon_nongeminalHname_mdict.keys()) and Cname in list(aatype_carbon_nongeminalHname_mdict[aa_type].keys())) and
                ((type(aatype_carbon_nongeminalHname_mdict[aa_type][Cname]) == list and atom_name in aatype_carbon_nongeminalHname_mdict[aa_type][Cname]) or
                (type(aatype_carbon_nongeminalHname_mdict[aa_type][Cname]) == str and atom_name == aatype_carbon_nongeminalHname_mdict[aa_type][Cname]) or
                (aa_type=="LEU" and Cname=="CD1" and atom_name == "HD1") or
                (aa_type=="LEU" and Cname=="CD2" and atom_name == "HD2") or
                (aa_type=="VAL" and Cname=="CG1" and atom_name == "HG1") or
                (aa_type=="VAL" and Cname=="CG2" and atom_name == "HG2"))
                ):
                peak_list.append([residue, Cname, Creson, atom_name, reson, spectrum_comment])
                break

    return peak_list


def write_peaks_assigned_in_previous_iteration(i_residue,
                                               clean_resid_assignments_dict,
                                               atom_index,
                                               xeasy_fout,
                                               HCNH_residue_NHresonances_dict,
                                               residues_with_written_NH_list):

    resid = get_resid_from_residue(i_residue)
    aa_type = aa1to3_dict[get_aa_type_from_residue(i_residue)]
    if resid in list(clean_resid_assignments_dict.keys()):    # if this was not a C-term, then peaks were assigned in iteration 1, otherwise write
                                                         # only the N & HN
        xeasy_lines_list, sparky_lines_list = [], []    # lines to save the xeasy and sparky lines for modification
        for xeasy_line_list in clean_resid_assignments_dict[resid]:
            xeasy_line_list[0] = atom_index
            xeasy_lines_list.append( xeasy_line_list )
            atom_index += 1

        ## NOW SAVE IN SPARKY FORMAT THE PEAKS THAT WERE ASSIGNED IN THE PREVIOUS ITERATION
        #assigned_peaks_list = get_peaks_from_xeasy(resid, clean_resid_assignments_dict)     # get all assigned peaks of this residue from xeasy format
        #for i in range:
        #    code
        #
        #    sparky_lines_list.append( [i_residue + Hydrogen_name, Carbon_name, i_residue + "N",float(Hydrogen_resonance), float(Carbon_resonance), float(N_resonance), float(HN_resonance)] )


        #Creson_Cname_lines_list = [] # list with the resonances and names of the printed carbons to avoid double printing (e.g. CB, CG, CD, CE)
        #try:    # DANGEROUS but works for the N-term residue (A200 in aLP)
        #    matched_HCNH_peaks_list = [peak for peak in matched_HCNH_residue_assignments_dict[i_residue] if peak[0]==i_residue]
        #except KeyError:
        #    return False    # function has been terminated prematurily
        #for peak in matched_HCNH_peaks_list:
        #    Carbon_name = peak[1]
        #    Carbon_resonance = peak[2]
        #    Hydrogen_name = peak[3]
        #    Hydrogen_resonance = peak[4]
        #    N_resonance = [peak[3] for peak in original_HCNH_peaks_dict[i_residue] if peak[1]==Hydrogen_resonance and peak[2]==Carbon_resonance][0]
        #    HN_resonance = [peak[4] for peak in original_HCNH_peaks_dict[i_residue] if peak[1]==Hydrogen_resonance and peak[2]==Carbon_resonance][0]
        #    Creson_Cname_lines_list.append([float(Carbon_resonance), Carbon_name])
        #    if not Hydrogen_name in written_nucleiNames_list:
        #        xeasy_lines_list.append( [atom_index, float(Hydrogen_resonance), 0.02, Hydrogen_name, resid, aa_type, "TOCSY-HCNH"] )  # append the Hydrogen resonance line
        #        atom_index += 1
        #        try:
        #            #sparky_fout.write("\t%s-%s-%s-H\t\t%.3f\t%.3f\t%.3f\t%.3f\t\n" % ( aa3to1_dict[aa_type] + residue[3:] + Hydrogen_name, Carbon_name, protein_alignment_list[residue2position_dict[residue]] + str(int(residue[3:])+1) + "N",float(Hydrogen_resonance), float(Carbon_resonance), float(N_resonance), float(HN_resonance)) )
        #            sparky_lines_list.append( [i_residue + Hydrogen_name, Carbon_name, i_residue + "N",float(Hydrogen_resonance), float(Carbon_resonance), float(N_resonance), float(HN_resonance)] )
        #        except IndexError:
        #            #sparky_fout.write("\t%s-%s-%s-H\t\t%.3f\t%.3f\t%.3f\t%.3f\t\n" % ( aa3to1_dict[aa_type] + residue[3:] + Hydrogen_name, Carbon_name, str(int(residue[3:])+1) + "N", float(Hydrogen_resonance), float(Carbon_resonance), float(N_resonance), float(HN_resonance)) )
        #            sparky_lines_list.append( [i_residue + Hydrogen_name, Carbon_name, i_residue + "N", float(Hydrogen_resonance), float(Carbon_resonance), float(N_resonance), float(HN_resonance)] )
        #Cnames_list = [cl[1] for cl in Creson_Cname_lines_list]
        #Cnames_set = set(Cnames_list)
        #for Carbon_name in Cnames_set:
        #    if Cnames_list.count(Carbon_name) > 1:  # if this carbon occurs in multiple C-H pairs, calculate its average resonance for the XEASY format
        #        Carbon_resonance = np.average([cl[0] for cl in Creson_Cname_lines_list if cl[1]==Carbon_name])
        #    else:
        #        Carbon_resonance = np.average([cl[0] for cl in Creson_Cname_lines_list if cl[1]==Carbon_name])
        #    #xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s\n" % (atom_index, float(Carbon_resonance), 0.2, Carbon_name, int(residue[3:]), aa_type))
        #    if not Carbon_name in written_nucleiNames_list:
        #        xeasy_lines_list.append( [atom_index, float(Carbon_resonance), 0.2, Carbon_name, resid, aa_type, "TOCSY-HCNH"] )
        #        atom_index += 1
        ## append unmatched TOCSY peaks, if any
        #if i_residue in residue_unmatched_TOCSY_peaks_dict.keys(): print "DEBUG: point 2, appending unmatched TOCSY peaks for residue", i_residue, residue_unmatched_TOCSY_peaks_dict[i_residue]
        #xeasy_lines_list, atom_index = append_unmatched_TOCSY_peaks(i_residue, xeasy_lines_list, atom_index, resid, aa_type)



        # WRITE C AND H LINES IN XEASY FORMAT
        for xeasy_line in xeasy_lines_list:
            xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t%s %s %s\n" % (xeasy_line[0], xeasy_line[1], xeasy_line[2], xeasy_line[3], xeasy_line[4], xeasy_line[5], xeasy_line[6], xeasy_line[7]) )
        ## WRITE C AND H LINES IN SPARKY FORMAT
        #for sparky_line in sparky_lines_list:
        #    sparky_fout.write("\t%s-%s-%s-H\t\t%.3f\t%.3f\t%.3f\t%.3f\t\n" % (sparky_line[0], sparky_line[1], sparky_line[2], sparky_line[3], sparky_line[4], sparky_line[5], sparky_line[6]))
        #update_HCNH_peaks(i_residue, sparky_lines_list) # update the labels in matched_HCNH_residue_assignments_dict

    # FINALLY WRITE N AND H LINES IN XEASY FORMAT
    if i_residue in list(HCNH_residue_NHresonances_dict.keys()): # only if this residue has N-H
        average_HN_resonance = HCNH_residue_NHresonances_dict[i_residue][1]
        stdev_HN_resonance = HCNH_residue_NHresonances_dict[i_residue][3]
        average_N_resonance = HCNH_residue_NHresonances_dict[i_residue][0]
        stdev_N_resonance = HCNH_residue_NHresonances_dict[i_residue][2]
        try:
            xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s %s\n" % (atom_index, float(average_HN_resonance), stdev_HN_resonance, "H", resid, aa_type, '')) # The N-HN of i sees the C-H of i-1
            atom_index += 1
            #print "\t"+str(atom_index)+"\t"+str(float(average_N_resonance))+"\t0.2\tN\t"+str(int(residue[3:])+1)
            xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s %s\n" % (atom_index, average_N_resonance, stdev_N_resonance, "N", resid, aa_type, ''))   # The N-HN of i sees the C-H of i-1
            atom_index += 1
        except IndexError:
            xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s %s\n" % (atom_index, float(average_HN_resonance), stdev_HN_resonance, "H", resid, aa_type, '')) # The N-HN of i sees the C-H of i-1
            atom_index += 1
            #print "\t"+str(atom_index)+"\t"+str(float(average_N_resonance))+"\t0.2\tN\t"+str(int(residue[3:])+1)
            xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s %s\n" % (atom_index, average_N_resonance, stdev_N_resonance, "N", resid, aa_type, ''))   # The N-HN of i sees the C-H of i-1
            atom_index += 1
        residues_with_written_NH_list.append(i_residue)


def annotate_HCNH_file(query_fname,
                        absolute_AAIGmatches_alignment_list,
                        absolute_matches_alignment_list,
                        matched_HCNH_residue_assignments_dict,
                        HCNH_residue_peak_intensity_mdict,
                        out_fname="4DHCNH_assignedall.sparky"):
    """
        FUNCTION to read the HCNH file with the assignments in SPARKY format, and to add labels.

        ARGS:
        query_fname:  HCNH file name

        RETURNS:
        i_residue & iplus1_residue in this function are AAIG not valid residue names!
        residue_assignments_dict:    i_residue -> [ (iminus1_residue, Cname, Cresonance, Hname, Hresonance), (iminus1_residue, Cname, Cresonance, Hname, Hresonance), ... ]
        residue_NHresonances_dict:   i_residue -> (average_Nresonance, average_Hresonance, stdev_N_resonance, stdev_Hreson)
    """
    print("Annotating HCNH Sparky file ", query_fname)

    with open(query_fname, 'r') as f:
        tmp_query_contents=f.readlines()    # contents of original query_fname (4D TOCSY or 4D HCNH) in 5 column format (name H C N HN)
    query_contents=[]
    for line in tmp_query_contents: # CLEANING THE SPECTRUM FROM IRRELEVANT LINES
        word_list = line.split()
        try:
            float(word_list[1]); float(word_list[2]); float(word_list[3]); float(word_list[4]); # checking if the line contains numbers (resonances)
            if len(word_list)==6:   # if there is an extra column chech if it is intensity number
                float(word_list[5])
            #print "DEBUG: appending line:", line
            query_contents.append(line)
        except (IndexError, ValueError):
            print("WARNING: Discarding HCNH line:", line)

    # remove duplicate lines from query_fname (4D HCNH)
    lines2remove_set = set()
    for qline1 in query_contents:
        try:
            counter = 0
            query1_words_list = qline1.split()
            q1_w1=float(query1_words_list[1])
            q1_w2=float(query1_words_list[2])
            q1_w3=float(query1_words_list[3])
            q1_w4=float(query1_words_list[4])
            q1_w5=float(query1_words_list[5])
            if len(query1_words_list) == 5: # if there is no intensity column, set it to 1.0
                q1_w5 = 1.0
            for qline2 in query_contents:
                query2_words_list = qline2.split()
                q2_w1=float(query2_words_list[1])
                q2_w2=float(query2_words_list[2])
                q2_w3=float(query2_words_list[3])
                q2_w4=float(query2_words_list[4])
                q2_w5=float(query2_words_list[5])
                if len(query2_words_list) == 5: # if there is no intensity column, set it to 1.0
                    q2_w5 = 1.0
                if approx_equal_proportional(q1_w1, q2_w1) and approx_equal_proportional(q1_w2, q2_w2) and approx_equal_proportional(q1_w3, q2_w3) and approx_equal_proportional(q1_w4, q2_w4) and approx_equal_proportional(q1_w5, q2_w5):
                    counter += 1
                    if counter > 1:
                        lines2remove_set.add(qline2)
                        #print "DEBUG: will remove line",qline2
        except (ValueError, IndexError):
            #print "WARNING: the 2nd and 3rd elements of the following root file line are not numbers:"
            #print "Root file line:", root_line
            continue
    # now remove the duplicate lines
    for qline in lines2remove_set:
        query_contents.remove(qline)

    # populate a dictionary with keys the resids and values the list of the respective C-H & N-H resonances from the spectrum file (TOCSY or HCNH)
    query_lineLists_list = []   # list of the lines of query frame in list form not in string
    for qline in query_contents:
        qline_list = qline.split()
        query_lineLists_list.append(qline_list)
    sorted_query_lineLists_list = sorted(query_lineLists_list, key=itemgetter(0))   # sort the spectrum lines by the assignment
    # print "DEBUG: sorted_query_lineLists_list=", sorted_query_lineLists_list

    # iterate over all lines of HCNH file and find the respective entries in matched_HCNH_residue_assignments_dict to copy the label (if it exists)
    query_lineList_forwriting_list = []
    qindex = 0
    for qindex in range(len(sorted_query_lineLists_list)):
        qline = sorted_query_lineLists_list[qindex]
        label_list = qline[0].split('-')
        if label_list[2] == '?':
            query_lineList_forwriting_list.append(qline)
            continue
        Hreson = float(qline[1])
        Creson = float(qline[2])
        # the follwoing applies only to N-H AAIGs (N-H mapped), not NX-HX or sidechains.
        residue = get_residue_from_AAIGsignature(remove_NH_suffix(qline[0].replace("?-?-", "")),
                                                 absolute_AAIGmatches_alignment_list, absolute_matches_alignment_list)
        try:
            print("DEBUG: trying to label qline=", qline)
            for peak in matched_HCNH_residue_assignments_dict[residue]:
                print("DEBUG: checking peak=", peak)
                print("DEBUG: Checking if ",peak[4],"==",Hreson,peak[2],"==", Creson, peak[0],'!=?' , peak[1],'!=?' , peak[3],'!=?')
                if peak[4]==Hreson and peak[2]==Creson and peak[0]!='?' and peak[1]!='?' and peak[3]!='?': # if this is the same peak, copy the label (if it exists)
                    # first assemble the label in the Sparky format, e.g. E292HB2-CB-E292N-H,
                    # from matched_HCNH_residue_assignments_dict['E292'] -> ['E292', 'CB', 63.879, 'HB2', 4.09],
                    # or T291HA-CA-E292N-H from matched_HCNH_residue_assignments_dict['E292'] -> ['T291', 'CA', 64.875, 'HA', 3.887]
                    label = peak[0]+peak[3]+"-"+peak[1]+"-"+residue+"N-H"   # only N-H mapped, not NX-HX or sidechains
                    norm_intensity = HCNH_residue_peak_intensity_mdict[peak[0]][(float(qline[2]), float(qline[1]))]
                    if len(qline) == 6 and type(norm_intensity) == float or type(norm_intensity) == int:    # if the intensity column is missing don't enter this condition
                        print("DEBUG: saving labeled qline:", [label, qline[1], qline[2], qline[3], qline[4], qline[5], str(norm_intensity)])
                        sorted_query_lineLists_list[qindex] = [label, qline[1], qline[2], qline[3], qline[4], qline[5], str(norm_intensity)]      # copy the label to the HCNH Sparky file line
                    else:   # this condition is satisfied even in the absence of intensity column
                        if len(qline) == 6:
                            print("DEBUG: saving labeled qline:", [label, qline[1], qline[2], qline[3], qline[4], qline[5]])
                            sorted_query_lineLists_list[qindex] = [label, qline[1], qline[2], qline[3], qline[4], qline[5]]      # copy the label to the HCNH Sparky file line
                        elif len(qline) == 5:
                            print("DEBUG: saving labeled qline:", [label, qline[1], qline[2], qline[3], qline[4]])
                            sorted_query_lineLists_list[qindex] = [label, qline[1], qline[2], qline[3], qline[4]]      # copy the label to the HCNH Sparky file line
                        else:
                            print("ERROR: this line does not contain the right number of columns:", qline)
                            raise IndexError
        except KeyError:    # if this residue in not in the HCNH dictionary, keep the original line
            continue
    # now write the annotated HCNH Sparky file
    with open(out_fname, 'w') as f:
        for qline_list in sorted_query_lineLists_list:
            f.write("\t".join(qline_list) + "\n")

# def write_Cterm_residues(absolute_matches_alignment_list,
#                          matched_HCNH_residue_assignments_dict,
#                          atom_index,
#                          xeasy_fout,
#                          sparky_fout,
#                          Cterm_residue_set,
#                          protein_alignment_list,
#                          args,
#                          TOCSY_assigned_residue_list,
#                          checked_residues_set,
#                          HCNH_residue_NHresonances_dict,
#                          residue_unmatched_TOCSY_peaks_dict,
#                          iter1_allowed_atoms_dict,
#                          allowed_aa_atoms_dict,
#                          iteration=1):
#
#     """
#     ##      OBSOLETE FUNCTION!
#     ## NOW DO COMPLETE PEAK ASSIGNMENT OF THE RESIDUES THAT HAVE HCNH BUT NO TOCSY PEAKS (C-TERM RESIDUES)
#     ## I.e. if a residue i+1 doesn't have TOCSY (we cannot get the resonances of residue i) but residue i has HCNH peaks, then use only the
#     ## HCNH peaks to do the assignment for residue i.
#     ##
#     """
#     print "DOING COMPLETE PEAK ASSIGNMENT OF THE RESIDUES THAT HAVE HCNH BUT NO TOCSY PEAKS."
#     print "DEBUG: TOCSY_assigned_residue_list=", TOCSY_assigned_residue_list
#     for i_residue in matched_HCNH_residue_assignments_dict.keys():
#         all_possible_assignments_list = []
#         if not i_residue in TOCSY_assigned_residue_list and len([x for x in matched_HCNH_residue_assignments_dict[i_residue] if x[0] == i_residue])==0 and i_residue in absolute_matches_alignment_list and absolute_matches_alignment_list[absolute_matches_alignment_list.index(i_residue)+1] == '-':
#             print "Assigning C-terminal residue ", i_residue, " only from HCNH."
#             # FIRST OF ALL MATCH THE TOCSY PEAKS OF RESIDUES i-1 and i+1 (IF APPLICABLE) TO THE HCNH PEAKS OF RESIDUE i
#
#             Cterm_residue_set.add(i_residue)    # remember that this residue was C-terminal
#             checked_residues_set.add(i_residue)
#             aa_type = aa1to3_dict[get_aa_type_from_residue(i_residue)]
#             print "DEBUG: 3-letter code aa_type=", aa_type
#             missing_carbons_list = [x for x in allowed_aa_atoms_dict[aa_type] if x[0] == "C"]   # all allowed carbons for this aa type
#             partly_assigned_carbons_list = []   # this list should be empty since the C-Term residues do not have TOCSY
#             # try:    # DANGEROUS but works for the N-term residue (A200 in aLP)
#             unassigned_HCNH_peaks_list = [peak for peak in matched_HCNH_residue_assignments_dict[i_residue] if peak[0]=="?" and peak[1]=="?" and peak[3]=="?"]
#             # except KeyError:
#             #     continue
#             for HCNH_peak in unassigned_HCNH_peaks_list:
#                 HCNH_Creson = HCNH_peak[2]
#                 HCNH_Hreson = HCNH_peak[4]
#                 all_possible_assignments_list.extend(get_probabilities_from_H_C_resonpair_2Dhist(i_residue,
#                                                                                                  HCNH_Hreson,
#                                                                                                  HCNH_Creson,
#                                                                                                  missing_carbons_list,
#                                                                                                  partly_assigned_carbons_list,
#                                                                                                  args,
#                                                                                                  HCNH_residue_peak_intensity_mdict,
#                                                                                                  matched_i_to_iplus1_peaks_mdict,
#                                                                                                  aa2pthres_dict,
#                                                                                                  aa2ithres_dict,
#                                                                                                  protein_alignment_list,
#                                                                                                  absolute_matches_alignment_list,
#                                                                                                  iteration=iteration))
#             [assignment.append('HCNH') for assignment in all_possible_assignments_list]    # we don't need to created a new list, the original all_possible_assignments_list will be updated
#             print "DEBUG: i_residue=", i_residue, "all_possible_assignments_list=", all_possible_assignments_list
#
#         if len(all_possible_assignments_list) != 0:
#             ## SELECT THE BEST ASSIGNMENT COMBINATION AND SAVE IT INTO XEASY AND SPARKY FILES
#             residue_residueTypesProbTupleList_dict = OrderedDict() # ordereddict with keys the AAIG of residue i and values lists of tuples of the form (residue i-1 matching aa type, average probability)
#             previous_aaindex = ""
#             print "Assigning HCNH chemical shifts to i_residue ", i_residue
#             #aa_type = residue[0:3]  # aa type in 3-letter code
#             aa_type = aa1to3_dict[get_aa_type_from_residue(i_residue)]
#             aatype_CHnucleiType_presenceProbSum_mdict, aatype_CHnucleiType_NresonIndex_mdict, HCNHindex_ResonancesTuple_dict, HCNHindex_spectrumType_dict = get_aatypes_from_all_H_C_resonpairs(all_possible_assignments_list, i_residue, partly_assigned_carbons_list, iteration, Probability_CS,
#                                         HCNH_residue_NHresonances_dict, matched_HCNH_residue_assignments_dict, protein_alignment_list, absolute_matches_alignment_list, iter1_allowed_atoms_dict, args)
#             if aatype_CHnucleiType_presenceProbSum_mdict != {} and aatype_CHnucleiType_NresonIndex_mdict != {}: # if no atom assignment could be made, skip this aa! (TEMPORARY FIX)
#                 #print "DEBUG: previous_TOCSY_aaindex=",previous_TOCSY_aaindex,"residue_residueTypesProbTupleList_dict[previous_TOCSY_aaindex]=",residue_residueTypesProbTupleList_dict[previous_TOCSY_aaindex]
#                 ## WRITE THE CHEMICAL SHIFTS ASSIGNMENTS IN XEASY FORMAT FOR ROSETTA
#                 # resid = absolute_matches_alignment_list.index(i_residue) + int(args.FIRST_RESIDUE_NUMBER) # OLD WAY
#                 resid = get_resid_from_residue(i_residue)
#                 xeasy_lines_list, sparky_lines_list = [], []    # lines to save the xeasy and sparky lines for modification
#                 HN_resonances_list, N_resonances_list = [], []
#                 Creson_Cname_lines_list = [] # list with the names of the printed carbons to avoid double printing (e.g. CB, CG, CD, CE)
#                 print "DEBUG: aatype_CHnucleiType_NresonIndex_mdict[aa_type].keys() =", aatype_CHnucleiType_NresonIndex_mdict[aa_type].keys()
#                 print "DEBUG: aatype_CHnucleiType_NresonIndex_mdict =", aatype_CHnucleiType_NresonIndex_mdict
#                 for CHnucleiType, NresonIndex  in aatype_CHnucleiType_NresonIndex_mdict[aa_type].items():
#                     ResonancesTuple = HCNHindex_ResonancesTuple_dict[NresonIndex]
#                     #print "DEBUG: aa_type=", aa_type, "CHnucleiType=", CHnucleiType, "NresonIndex=", NresonIndex, "ResonancesTuple=", ResonancesTuple
#                     Carbon_name = CHnucleiType.split('_')[0]
#                     Hydrogen_name = CHnucleiType.split('_')[1]
#                     # Rename un-resolved Geminal protons to QA, QB, QG, QD, QE.
#                     if len(aa_carbonBondedHydrogensDict_dict[aa_type][Carbon_name]) > 1:    # if this carbon has geminal protons
#                         #print "DEBUG: aa_type=", aa_type, "Carbon_name=", Carbon_name
#                         #print "DEBUG: aatype_CHnucleiType_NresonIndex_mdict[aa_type].keys()=", aatype_CHnucleiType_NresonIndex_mdict[aa_type].keys()
#                         #print "DEBUG: aa_carbonBondedHydrogensDict_dict[aa_type][Carbon_name]=", aa_carbonBondedHydrogensDict_dict[aa_type][Carbon_name]
#                         if [Carbon_name+"_" in x for x in aatype_CHnucleiType_NresonIndex_mdict[aa_type].keys()].count(True) < len(aa_carbonBondedHydrogensDict_dict[aa_type][Carbon_name]):
#                             #print "DEBUG entered condition 1!"
#                             if aa_type in aatype_carbon_methylHydrogens_mdict.keys() and Carbon_name in aatype_carbon_methylHydrogens_mdict[aa_type].keys():
#                                 #print "DEBUG entered condition 2!"
#                                 new_Hydrogen_name = aatype_carbon_methylHydrogens_mdict[aa_type][Carbon_name]
#                                 Hydrogen_name = new_Hydrogen_name
#                             #else:
#                             #    print "DEBUG entered condition 3!"
#                             #    new_Hydrogen_name = 'Q'+Hydrogen_name[1]
#                             #    Hydrogen_name = new_Hydrogen_name
#                     # CONDITIONS:
#                     # 1)
#                     if aa_type in aatype_carbon_methylHydrogens_mdict.keys() and Carbon_name in aatype_carbon_methylHydrogens_mdict[aa_type].keys() and (not aa_type in aa_CarbonListwithDegenerateH_dict.keys() or Carbon_name in aa_CarbonListwithDegenerateH_dict[aa_type]):
#                         # DO NOT RENAME VAL-HG1,HG2 --> VAL-QG1,QG2 and LEU-HD1,HD2-->LEU-QD1,QD2 because the code will be confused! Do it at the very end!!!
#                         if not(aa_type in ["LEU", "VAL"] and Carbon_name in ["CG1", "CG2", "CD1", "CD2"]):
#                             #print "DEBUG entered condition 4! aa_type = ", aa_type, "Carbon_name=", Carbon_name
#                             new_Hydrogen_name = aatype_carbon_methylHydrogens_mdict[aa_type][Carbon_name]
#                             Hydrogen_name = new_Hydrogen_name
#                     Carbon_resonance = ResonancesTuple[1]
#                     Hydrogen_resonance = ResonancesTuple[0]
#                     N_resonance = ResonancesTuple[2]
#                     HN_resonance = ResonancesTuple[3]
#                     #HN_resonances_list.append(float(ResonancesTuple[3]))   # UNECESSARY CAUSE WE ALREADY SAVED THEM IN HCNH_residue_NHresonances_dict
#                     #N_resonances_list.append(float(ResonancesTuple[2]))
#
#                     #print "\t"+str(atom_index)+"\t"+str(Carbon_resonance)+"\t0.2\t"+str(Carbon_name)+"\t"+str(residue[3:])
#                     Creson_Cname_lines_list.append([float(Carbon_resonance), Carbon_name])
#
#                     #print "\t"+str(atom_index)+"\t"+str(Hydrogen_resonance)+"\t0.02\t"+str(Hydrogen_name)+"\t"+str(residue[3:])
#                     #xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s\n" % (atom_index, float(Hydrogen_resonance), 0.02, Hydrogen_name, int(residue[3:]), aa_type))
#                     xeasy_lines_list.append( [atom_index, float(Hydrogen_resonance), 0.02, Hydrogen_name, resid, aa_type, "HCNH (C-term hanging residue) iter"+str(iteration)] )  # append the Hydrogen resonance line
#                     print "DEBUG: point 5 appending iter: ", [atom_index, float(Hydrogen_resonance), 0.02, Hydrogen_name, resid, aa_type, "HCNH (C-term hanging residue) iter"+str(iteration)]
#                     atom_index += 1
#                     if N_resonance != None and HN_resonance != None:    # if no N-H was found (e.g. N-term residue) don't write the i-i Sparky lines for this residue, only the xeasy lines
#                         try:
#                             #sparky_fout.write("\t%s-%s-%s-H\t\t%.3f\t%.3f\t%.3f\t%.3f\t\n" % ( aa3to1_dict[aa_type] + residue[3:] + Hydrogen_name, Carbon_name, protein_alignment_list[residue2position_dict[residue]] + str(int(residue[3:])+1) + "N",float(Hydrogen_resonance), float(Carbon_resonance), float(N_resonance), float(HN_resonance)) )
#                             sparky_lines_list.append( [i_residue + Hydrogen_name, Carbon_name, i_residue + "N",float(Hydrogen_resonance), float(Carbon_resonance), float(N_resonance), float(HN_resonance)] )
#                         except IndexError:
#                             #sparky_fout.write("\t%s-%s-%s-H\t\t%.3f\t%.3f\t%.3f\t%.3f\t\n" % ( aa3to1_dict[aa_type] + residue[3:] + Hydrogen_name, Carbon_name, str(int(residue[3:])+1) + "N", float(Hydrogen_resonance), float(Carbon_resonance), float(N_resonance), float(HN_resonance)) )
#                             sparky_lines_list.append( [i_residue + Hydrogen_name, Carbon_name, i_residue + "N", float(Hydrogen_resonance), float(Carbon_resonance), float(N_resonance), float(HN_resonance)] )
#
#
#                 Cnames_list = [cl[1] for cl in Creson_Cname_lines_list]
#                 Cnames_set = set(Cnames_list)
#                 for Carbon_name in Cnames_set:
#                     if Cnames_list.count(Carbon_name) > 1:  # if this carbon occurs in multiple C-H pairs, calculate its average resonance for the XEASY format
#                         Carbon_resonance = round(np.average([cl[0] for cl in Creson_Cname_lines_list if cl[1]==Carbon_name]), 3)
#                     else:
#                         Carbon_resonance = round(np.average([cl[0] for cl in Creson_Cname_lines_list if cl[1]==Carbon_name]), 3)
#                     #xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s\n" % (atom_index, float(Carbon_resonance), 0.2, Carbon_name, int(residue[3:]), aa_type))
#                     xeasy_lines_list.append( [atom_index, float(Carbon_resonance), 0.2, Carbon_name, resid, aa_type, "HCNH (C-term hanging residue) iter"+str(iteration)] )
#                     print "DEBUG: point 6 appending iter: ", [atom_index, float(Carbon_resonance), 0.2, Carbon_name, resid, aa_type, "HCNH (C-term hanging residue) iter"+str(iteration)]
#                     atom_index += 1
#                 # append unmatched TOCSY peaks, if any
#                 if i_residue in residue_unmatched_TOCSY_peaks_dict.keys(): print "DEBUG: point 4, appending unmatched TOCSY peaks for residue", i_residue, residue_unmatched_TOCSY_peaks_dict[i_residue]
#                 xeasy_lines_list, atom_index = append_unmatched_TOCSY_peaks(i_residue, xeasy_lines_list, atom_index, resid, aa_type, residue_unmatched_TOCSY_peaks_dict)
#
#                 # RENAME DEGENERATE PROTONS
#                 new_xeasy_lines_list, new_sparky_lines_list = rename_protons(xeasy_lines_list, sparky_lines_list)
#                 print "DEBUG: writing flanked residues new_xeasy_lines_list=", new_xeasy_lines_list
#                 # WRITE C AND H LINES IN XEASY FORMAT
#                 for xeasy_line in new_xeasy_lines_list:
#                     xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s %s\n" % (xeasy_line[0], xeasy_line[1], xeasy_line[2], xeasy_line[3], xeasy_line[4], xeasy_line[5], xeasy_line[6] ))
#                 # WRITE C AND H LINES IN SPARKY FORMAT
#                 for sparky_line in new_sparky_lines_list:
#                     sparky_fout.write("\t%s-%s-%s-H\t\t%.3f\t%.3f\t%.3f\t%.3f\t\n" % (sparky_line[0], sparky_line[1], sparky_line[2], sparky_line[3], sparky_line[4], sparky_line[5], sparky_line[6] ))
#                 update_HCNH_peaks(i_residue, matched_HCNH_residue_assignments_dict, protein_alignment_list, absolute_matches_alignment_list, args, new_sparky_lines_list)

