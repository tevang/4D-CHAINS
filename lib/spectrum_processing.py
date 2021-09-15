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
from .global_func import *

##################################### FUNCTIONS TO PROCESS THE 3 INPUT SPECTRA (HSQC, 4D TOCSY, 4D NOESY) ############################################



##  FIND NON-OVERLAPPING RESONANCE GROUPS IN THE HSQC SPECTRUM
def find_nonoverlapping_groups_in_root(remaining_root_contents, rtolH, rtolN, keep_closest=False):
    """
        RETURNS:
        nonoverlapping_groups_lines_and_tolerances:    list with the line of each non-overlapping group of the root spectrum along with the rtolH and rtolN
        remaining_groups_lines:                        list with the lines of the Root spectrum that contain overlapping groups with the current torelances
    """
    
    nonoverlapping_groups_lines_and_tolerances = []
    remaining_groups_lines = []
    
    if keep_closest == False:
        for root1_line in remaining_root_contents:
            found_overlap = False
            #print "DEBUG: root1_line=", root1_line
            root1_words_list = root1_line.split()
            try:
                if not root1_words_list[0][-3:] in ['N-H', 'NX-HX']:   # if it's not a backbone amide, skip it
                    continue
                root1_name = root1_words_list[0]
                root1_N_resonance = float(root1_words_list[1])
                root1_HN_resonance = float(root1_words_list[2])
                for root2_line in remaining_root_contents:
                    if root1_line == root2_line:    # in case we are looking at the same line
                        continue
                    #print "DEBUG: root2_line=",root2_line
                    root2_words_list = root2_line.split()
                    try:
                        root2_name = root2_words_list[0]
                        root2_N_resonance = float(root2_words_list[1])
                        root2_HN_resonance = float(root2_words_list[2])
                        if ( (root2_HN_resonance -rtolH) <= root1_HN_resonance <= (root2_HN_resonance +rtolH) ) and ( (root2_N_resonance -rtolN) <= root1_N_resonance <= (root2_N_resonance +rtolN) ):
                            remaining_groups_lines.append(root1_line)   # 
                            found_overlap = True
                            break
                    except (ValueError, IndexError):
                        #print "WARNING: the 4th and 5th elements of the following root file line are not numbers:"
                        #print "Root file line:", root_line
                        continue
                if found_overlap == False:  # if no noverlapping line was found for this line then save it
                    nonoverlapping_groups_lines_and_tolerances.append((root1_line, rtolH, rtolN))
            except (ValueError, IndexError):
                #print "WARNING: the 2nd and 3rd elements of the following root file line are not numbers:"
                #print "Root file line:", root_line
                continue
    
    #if keep_closest == True:
    #    for root1_line in remaining_groups_lines:
    #        matching_lines = []
    #        found_overlap = False
    #        #print "DEBUG: root1_line=", root1_line
    #        root1_words_list = root1_line.split()
    #        try:
    #            if root1_words_list[0][-3:] != "N-H":   # if it's not a backbone amide, skip it
    #                continue
    #            root1_name = root1_words_list[0]
    #            root1_N_resonance = float(root1_words_list[1])
    #            root1_HN_resonance = float(root1_words_list[2])
    #            for root2_line in remaining_root_contents:
    #                if root1_line == root2_line:    # in case we are looking at the same line
    #                    continue
    #                #print "DEBUG: root2_line=",root2_line
    #                root2_words_list = root2_line.split()
    #                try:
    #                    root2_name = root2_words_list[0]
    #                    root2_N_resonance = float(root2_words_list[1])
    #                    root2_HN_resonance = float(root2_words_list[2])
    #                    if ( (root2_HN_resonance -rtolH) <= root1_HN_resonance <= (root2_HN_resonance +rtolH) ) and ( (root2_N_resonance -rtolN) <= root1_N_resonance <= (root2_N_resonance +rtolN) ):
    #                        found_overlap = True
    #                        matching_lines.append(root1_line)
    #                except (ValueError, IndexError):
    #                    #print "WARNING: the 4th and 5th elements of the following root file line are not numbers:"
    #                    #print "Root file line:", root_line
    #                    continue
    #            if found_overlap == False:
    #                delta_previous_root = np.sqrt((previous_root_HN_resonance - query_HN_resonance)**2 + ((previous_root_N_resonance - query_N_resonance)/6)**2)
    #                delta_root = np.sqrt((root_HN_resonance - query_HN_resonance)**2 + ((root_N_resonance - query_N_resonance)/6)**2)
    #                if delta_previous_root > delta_root: 
    #                    nonoverlapping_groups_lines_and_tolerances.append((root1_line, rtolH, rtolN))
    #        except (ValueError, IndexError):
    #            #print "WARNING: the 2nd and 3rd elements of the following root file line are not numbers:"
    #            #print "Root file line:", root_line
    #            continue
    #
    #print "DEBUG: len(nonoverlapping_groups_lines_and_tolerances)=", len(nonoverlapping_groups_lines_and_tolerances)
    print("DEBUG: nonoverlapping_groups_lines_and_tolerances =", nonoverlapping_groups_lines_and_tolerances)
    print("DEBUG: remaining_groups_lines = ", remaining_groups_lines)
    return nonoverlapping_groups_lines_and_tolerances


def load_spectrum_lines(query_fname, spectrum_combo):
    """
    Load the contents of TOCSY or NOESY and remove replicate lines.
    :param query_fname:
    :param spectrum_combo:
    :return:
    """
    with open(query_fname, 'r') as f:
        tmp_query_contents = f.readlines()  # contents of original query_fname (4D TOCSY or 4D NOESY) in 5 column format (name H C N HN)
    query_lines = []
    query_CS_set = set()  # store the numbers here to discard replicate lines
    for line in tmp_query_contents:  # CLEANING THE SPECTRUM FROM IRRELEVANT AND REPLICATE LINES
        word_list = line.split()
        try:
            float(word_list[1]);
            float(word_list[2]);
            float(word_list[3]);
            float(word_list[4]);  # checking if the line contains numbers (resonances)
            if len(word_list) == 6:  # if there is intensity column check if it is a number
                float(word_list[5])
                word_list[5] = str(abs(float(word_list[5])))  # convert all intensities to positive numbers
                CS_values = (float(word_list[1]), float(word_list[2]), float(word_list[3]),
                             float(word_list[4]), float(word_list[5]))
            else:
                CS_values = (float(word_list[1]), float(word_list[2]), float(word_list[3]),
                             float(word_list[4]))
            if not CS_values in query_CS_set:
                query_lines.append(" ".join(word_list) + "\n")
                query_CS_set.add(CS_values)
            else:
                print(bcolors.WARNING + "WARNING: Discarding replicate " + spectrum_combo + " line: " + line + bcolors.ENDC)
        except (IndexError, ValueError):
            print(bcolors.WARNING + "WARNING: Discarding invalid " + spectrum_combo + " line: " + line + bcolors.ENDC)
            print("DEBUG: copy_aaindices_from_root_spectrum_2 point1 word_list=", word_list)

    # remove duplicate lines from query_fname (4D TOCSY or 4D NOESY)
    lines2remove_set = set()
    for qline1 in query_lines:
        try:
            counter = 0
            query1_words_list = qline1.split()
            q1_w1 = float(query1_words_list[1])
            q1_w2 = float(query1_words_list[2])
            q1_w3 = float(query1_words_list[3])
            q1_w4 = float(query1_words_list[4])
            for qline2 in query_lines:
                query2_words_list = qline2.split()
                q2_w1 = float(query2_words_list[1])
                q2_w2 = float(query2_words_list[2])
                q2_w3 = float(query2_words_list[3])
                q2_w4 = float(query2_words_list[4])
                if approx_equal(q1_w1, q2_w1) and approx_equal(q1_w2, q2_w2) and approx_equal(q1_w3,
                                                                                              q2_w3) and approx_equal(
                        q1_w4, q2_w4):
                    counter += 1
                    if counter > 1:
                        lines2remove_set.add(qline2)
                        # print "DEBUG: will remove line",qline2
        except (ValueError, IndexError):
            # print "WARNING: the 2nd and 3rd elements of the following root file line are not numbers:"
            # print "Root file line:", root_line
            continue
    # now remove the duplicate lines
    for qline in lines2remove_set:
        query_lines.remove(qline)

def copy_aaindices_from_root_spectrum_2(root_contents, query_fname, spectrum_combo, TOCSY_sidechain_resonances_list=None):
    
    """
    FUNCTION that finds the closest TOCSY peaks to each HSQC spectrum group (namely, each HSQC line). If the lower
    distance of a TOCSY peak from any HSQC group is greater than (0.04, 0.4) then it will be considered as garbage
    and will be discarded from TOCSYnum.list file.
    """
    tolH = 0.04
    tolN = 0.4

    query_lines = load_spectrum_lines(query_fname, spectrum_combo)

    AAIGsign_Hreson_Creson_tuple_list = []
    for triplet in root_contents:
        root_line = triplet[0]
        #print "DEBUG: root_line=", root_line
        root_words_list = root_line.split()
        try:
            root_name = root_words_list[0]
            root_N_resonance = float(root_words_list[1])
            root_HN_resonance = float(root_words_list[2])
            for q_index in range(0,len(query_lines)):
                query_line = query_lines[q_index]
                #print "DEBUG: query_line=",query_line
                try:
                    query_words_list = query_line.split()
                    query_N_resonance = float(query_words_list[3])
                    query_HN_resonance = float(query_words_list[4])
                    # First remove the garbage that are far away from any peak
                    if ( (query_HN_resonance -tolH) <= root_HN_resonance <= (query_HN_resonance +tolH) ) and ( (query_N_resonance -tolN) <= root_N_resonance <= (query_N_resonance +tolN) ):
                        #print "DEBUG: ?-?-"+root_name+"\t"+query_words_list[1]+"\t"+query_words_list[2]+"\t"+query_words_list[3]+"\t"+query_words_list[4]
                        AAIGsign = root_name
                        if re.search("^\s*\?-\?-\?-\?\s+", query_lines[q_index]):   # if the line does not have an AAIG group already, assign the currect AAIG group to it
                            query_lines[q_index] = re.sub("\?-\?-\?-\?", AAIGsign, query_lines[q_index]).replace("\\", "")
                            Hresonance = query_lines[q_index].split()[1]
                            Cresonance = query_lines[q_index].split()[2]
                            AAIGsign_Hreson_Creson_tuple_list.append((AAIGsign, Hresonance, Cresonance))
                        else: # if it has already been assigned an AAIG group, check if the currect AAIG group has closer N & HN resonances
                            #print "DEBUG: Already assigned AAIG group. Line: ", query_contents[q_index]
                            mo = re.search('^\s*([A-Za-z0-9-]+)\s+', query_lines[q_index])
                            if mo:
                                previous_AAIGsign = mo.group(1)     # the AAIG group that is already assigned to this peak from a previous round
                                for tmp_triplet, tmp_rtolH, tmp_rtolN in root_contents:   # iterate over all AAIG groups until you find previous_AAIGsign, in order to retrieve its N & HN resonances
                                    tmp_root_line = tmp_triplet[0]
                                    tmp_root_words_list = tmp_root_line.split()
                                    try:
                                        tmp_root_name = tmp_root_words_list[0]
                                        tmp_AAIGsign = tmp_root_name
                                        if tmp_AAIGsign == previous_AAIGsign:   # we found the N & HN resonances of the assigned AAIG group
                                            #print "DEBUG: we found the N & HN resonances of the assigned AAIG group:", tmp_root_words_list
                                            previous_root_N_resonance = float(tmp_root_words_list[1])
                                            previous_root_HN_resonance = float(tmp_root_words_list[2])
                                            # if the current AAIG group has N & HN resonances closer to the query line than the previously assigned N & HN resonances, change the AAIG group
                                            #print "DEBUG: query_N_resonance=", query_N_resonance, "previous_root_N_resonance=", previous_root_N_resonance, "root_N_resonance=", root_N_resonance
                                            #print "DEBUG: query_HN_resonance=", query_HN_resonance, "previous_root_HN_resonance=", previous_root_HN_resonance, "root_HN_resonance=", root_HN_resonance
                                            delta_previous_root = np.sqrt((previous_root_HN_resonance - query_HN_resonance)**2 + ((previous_root_N_resonance - query_N_resonance)/6)**2)
                                            delta_root = np.sqrt((root_HN_resonance - query_HN_resonance)**2 + ((root_N_resonance - query_N_resonance)/6)**2)
                                            if delta_previous_root > delta_root:
                                                #print "DEBUG: substituting previous_AAIGsign", previous_AAIGsign, " for ", AAIGsign, " in ", query_contents[q_index]
                                                if len(query_words_list) == 6:
                                                    query_lines[q_index] = " ".join(
                                                        [AAIGsign, query_words_list[1], query_words_list[2],
                                                         query_words_list[3], query_words_list[4],
                                                         query_words_list[5], "\n"])
                                                else:
                                                    query_lines[q_index] = " ".join(
                                                        [AAIGsign, query_words_list[1], query_words_list[2],
                                                         query_words_list[3], query_words_list[4], "\n"])
                                                #query_contents[q_index] = re.sub(previous_AAIGsign, AAIGsign, query_contents[q_index])
                                            # otherwise keep the previous AAIG group
                                            # ATTENTION: THIS IS NOT ENTIRELY CORRECT, BECAUSE THE N RESONANCE MAY BE CLOSER BUT THE H NOT!
                                            break
                                    except (ValueError, IndexError):
                                        #print "WARNING: the 2nd and 3rd elements of the following root file line are not numbers:"
                                        #print "Root file line:", root_line
                                        continue
                            else:
                                print("ERROR: wrong line (modified or not) in ", query_fname)
                                print(query_lines[q_index])
                                sys.exit(1)
                            
                except (ValueError, IndexError):
                    #print "WARNING: the 4th and 5th elements of the following query file line are not numbers:"
                    #print "Query file line:", query_line
                    continue
        except (ValueError, IndexError):
            #print "WARNING: the 2nd and 3rd elements of the following root file line are not numbers:"
            #print "Root file line:", root_line
            continue
    
    print("DEBUG: ", spectrum_combo," contents: ")
    print("".join(query_lines))
    with open(replace_alt(query_fname, [".list",".txt",".dat",".sparky"], "") + "num.list", 'w') as f:
        for line in query_lines:
            try:
                print("DEBUG: saving line:", line)
                word_list = line.split()
                if len(word_list) == 6:
                    appendix = "\t"+word_list[5]
                else:
                    appendix = ""
                if spectrum_combo == "TOCSY":
                    if not "?" in word_list[0]:
                        f.write("?-?-"+word_list[0]+"\t\t"+word_list[1]+"\t"+word_list[2]+"\t"+word_list[3]+"\t"+word_list[4]+appendix+"\n")
                    else:
                        f.write("?-?-?-?\t\t"+word_list[1]+"\t"+word_list[2]+"\t"+word_list[3]+"\t"+word_list[4]+appendix+"\n")
                elif spectrum_combo == "HCNH":
                    if not "?" in word_list[0]:
                        f.write("?-?-"+word_list[0]+"\t\t"+word_list[1]+"\t"+word_list[2]+"\t"+word_list[3]+"\t"+word_list[4]+appendix+"\n")
                    else:
                        f.write("?-?-?-?\t\t"+word_list[1]+"\t"+word_list[2]+"\t"+word_list[3]+"\t"+word_list[4]+appendix+"\n")
            except IndexError:
                print("WARNING: Discarding", spectrum_combo, "line:", line)
                #print "DEBUG: line=", line
                print("DEBUG: copy_aaindices_from_root_spectrum_2 point 2 word_list=", word_list)
                continue
    
    # #print "DEBUG: AAIGsign_Hreson_Creson_tuple_list = ",AAIGsign_Hreson_Creson_tuple_list
    # query_lines = []    # list of the query_fname lines in 3 column format (namely without the root spectrum N & HN resonances) 
    # with open(query_fname.replace(".list","").replace(".txt","").replace(".dat","").replace(".sparky", "")+"_aa-H-C.list", 'w') as f:
    #     for triplet in AAIGsign_Hreson_Creson_tuple_list:
    #         line = "\t?-?-?\t"+triplet[0]+"\t"+triplet[1]+"\t"+triplet[2]+"\n"
    #         f.write(line)
    #         query_lines.append(line)
    
    if spectrum_combo == "TOCSY":
        # save unassigned TOCSY lines that correspond to side-chains. The resonances must be later removed from NOESY
        TOCSY_sidechain_resonances_list = []
        for qline in query_lines:
            try:
                print("DEBUG: qline.split()[0]:",qline.split()[0])
                if qline.split()[0] == "?-?-?-?":
                    TOCSY_sidechain_resonances_list.append(qline)
            except IndexError:
                continue
        for TOCSY_line in TOCSY_sidechain_resonances_list:  # TEMPORARILY REMOVE ALL "?-?-?-?" LINES (NOT ALL OF THEM ARE SIDE CHAIN RESONANCES FROM TOCSY FILE CONTENTS)
            query_lines.remove(TOCSY_line)   # now you won't see TOCSY index groups like "?-?-?-?" in the connectivities and aa type prediction files
        return query_lines, TOCSY_sidechain_resonances_list
    elif spectrum_combo == "HCNH":
        lines2remove_set = set()
        for qline in query_lines:    # FIND ALL LINES WITH "?-?-?-?"
            try:
                #print "DEBUG: qline.split()[0]:",qline.split()[0]
                if qline.split()[0] == "?-?-?-?":
                    lines2remove_set.add(qline)
            except IndexError:
                continue
        # remove the sidechain resonances from NOESY
        #for TOCSY_line in TOCSY_sidechain_resonances_list:
        for qline in lines2remove_set:  # REMOVE ALL LINES WITH "?-?-?-?"
            #print "DEBUG: removing side chain line from NOESY: ",TOCSY_line
            try:
                #query_contents.remove(TOCSY_line)
                query_lines.remove(qline)
            except ValueError:
                print("DEBUG: side chain line not present in NOESY!")
                continue
        return query_lines
    
    pass

def remove_native_peaks_from_NOESY(TOCSY_contents, NOESY_contents, tolH, tolC):
    """
        When you match TOCSY peaks of AAIG(i) to NOESY peaks in order to find the AAIG(i-1), you don't need to use the AAIG(i) in NOESY (these
        are commone between NOESY and TOCSY). This function iterates over all TOCSY AAIG(i), finds the NOESY AAIG(i) and removed the peaks that
        match from NOESY_contents.
There is a problem with the removal of peaks from the same AAIG. Look at this example from nEIt. There are 3 consecutive E in the sequence, E37,E38,E39.

In [62]: [w for w in NOESY_contents if w[0]=='E38']
Out[62]:
[['E38', '4.684', '55.456', '121.017', '8.064', '10353967.0'],
 ['E38', '1.022', '22.189', '120.987', '8.066', '3772312.0'],
 ['E38', '3.054', '42.354', '120.975', '8.069', '3508862.0'],
 ['E38', '4.011', '60.287', '121.009', '8.066', '8719190.0'],
 ['E38', '2.136', '28.789', '120.997', '8.067', '32965242.0'],
 ['E38', '4.207', '59.453', '120.970', '8.066', '41168996.0'],
 ['E38', '2.432', '36.962', '121.004', '8.066', '22151098.0'],
 ['E38', '2.337', '36.957', '121.002', '8.067', '16331055.0'],
 ['E38', '2.251', '29.850', '121.003', '8.067', '60914216.0']]

In [63]: [w for w in TOCSY if w[0] in ['E38']]
Out[63]:
[['E38', '4.010', '60.370', '120.991', '8.063'],
 ['E38', '2.429', '36.937', '120.991', '8.063'],
 ['E38', '2.343', '36.943', '120.986', '8.063'],
 ['E38', '2.133', '28.864', '120.995', '8.061']]

In [64]: [w for w in TOCSY if w[0] in ['E39']]
Out[64]:
[['E39', '4.207', '59.379', '121.895', '8.254'],
 ['E39', '2.430', '36.954', '121.892', '8.253'],
 ['E39', '2.328', '36.986', '121.886', '8.253'],
 ['E39', '2.250', '29.810', '121.908', '8.254']]

As you may notice, there is only one pair of carbons 36 ppm (only one methylene) in NOESY. Normally, it should be two at least because both E38 and E39 have this type of methylene (E37 too). If I remove the TOCSY peaks of E38 from the NOESY list I end up with this:

In [61]: [w for w in clean_NOESY_contents if w[0]=='E38']
Out[61]:
[['E38', '4.684', '55.456', '121.017', '8.064', '10353967.0'],
 ['E38', '1.022', '22.189', '120.987', '8.066', '3772312.0'],
 ['E38', '3.054', '42.354', '120.975', '8.069', '3508862.0'],
 ['E38', '4.207', '59.453', '120.970', '8.066', '41168996.0'],
 ['E38', '2.251', '29.850', '121.003', '8.067', '60914216.0']]

As you can see there are no carbons at 36 ppm to match with E39. This results in lower occupancy ('E38', 2, 4)<--E39.

This was just one scenario (when you have consecutive residues of the same type you lose connectivities). Although the connectivities become fewer, the NH-mapping deteriorates when I remove NOESY peaks from the same AAIG when I do the TOCY->NOESY matching. Therefore I suggest to leave it as it is.
    """
    print("Removing native peaks from NOESY Amino Acid Index Groups.")
    NOESY_contents2remove = []  # the NOESY lines that contains peaks of AAIG(i)
    for TOCSY_wordlist in TOCSY_contents:
        try:
            TOCSY_AAIG=str(TOCSY_wordlist[0])    # residue i
            TOCSY_H_resonance=float(TOCSY_wordlist[1]) # aliphatic H resonance of residue i-1 
            TOCSY_C_resonance=float(TOCSY_wordlist[2]) # aliphatic C (Ca,Cb,Cc,Cg,Ce,etc.) resonance of residue i-1; this C is covalently bonded to the above H
            #print "DEBUG: TOCSY_AAIG=", TOCSY_AAIG
            NOESY_AAIG_wordlists = [words for words in NOESY_contents if words[0]==TOCSY_AAIG]
            for NOESY_wordlist in NOESY_AAIG_wordlists:
                NOESY_w2=float(NOESY_wordlist[1]) # resonance of an aliphatic H close to residue (i-1)
                NOESY_w3=float(NOESY_wordlist[2]) # resonance of an aliphatic C close to residue (i-1);
                                                  # this C is covalently bonded to the above H
                if ( (NOESY_w2-tolH) <= TOCSY_H_resonance <= (NOESY_w2+tolH) ) and \
                        ( (NOESY_w3-tolC) <= TOCSY_C_resonance <= (NOESY_w3+tolC) ):
                    NOESY_contents2remove.append(NOESY_wordlist)
                    # print "DEBUG: removing NOESY native peak:", " ".join(NOESY_wordlist)
                    break
        except (ValueError, IndexError):
            print("WARNING: the 2nd, 3rd or 4th elements of the following TOCSY file line are not numbers:")
            print("TOCSY file line:", TOCSY_wordlist)
            continue
    
    cleaned_NOESY_contents = [wordlist for wordlist in NOESY_contents if not wordlist in NOESY_contents2remove]
    
    return cleaned_NOESY_contents


def get_possible_connectivities(TOCSY_contents, NOESY_lines, tolH, tolC, clean_native_peaks=False):
    """
    FUNCTION that find all the NOESY aa indices that have resonances w2,w3 that match at least one w2,w3 TOCSY resonance pair.
    
    ARGUMENTS:
    TOCSY_contents:     list of lists, sorted by the 2nd element (TOCSY AAIG); these are the contents of TOCSY file in a convenient form.
                        e.g. TOCSY_contents= [['A118', '4.430', '60.596', '127.371', '8.391'], ...]
                        where 'A118' is the i AAIG that has the specified N,HN shifts. The C,H shifts correspond to i-1 AAIG (117).
    NOESY_lines:     list of lists; these are the contents of NOESY file but not sorted. Same format as TOCSY_contents, the label corresponds
                    to the AAIG with these N,HN shifts. E.g.
                    NOESY_contents= [['I88', '0.400', '11.833', '119.817', '7.368', '64604988.0'], ...]
    clean_native_peaks:   from each NOESY AAIG(i), remove the peaks that correspond to residue i. (DANGEROUS)
    
    RETURNS:
    TOCSYAAIG_NOESYAAIG_OccupancyNumOfResonancesList_dict:    a multidimensional dictionary with structure:
    TOCSY AAIG => NOESY AAIG => [occupancy, total resonance number] , namely
    TOCSY AAIG => NOESY AAIG => [number of matched TOCSY w2,w3 resonances in NOESY for that NOESY AAIG, total number of TOCSY w2,w3 resonances for that TOCSY AAIG]
    """
    
    NOESY_contents = [line.split() for line in NOESY_lines]    # split NOESY lines into words to be in same format as TOCSY_contents
    # # print "DEBUG: TOCSY_contents=", TOCSY_contents
    # # print "DEBUG: NOESY_contents=", NOESY_contents
    # # print "DEBUG: tolH=", tolH
    # # print "DEBUG: tolC=", tolC
    # # sys.exit(1)
    if clean_native_peaks:
        # remove NOESY peaks of AAIG(i) to improve the TOCSY->NOESY matching
        NOESY_contents = remove_native_peaks_from_NOESY(TOCSY_contents, NOESY_contents, tolH, tolC)
    
    # dictionary with keys the TOCSY aa indices, and values
    TOCSYAAIG_NOESYAAIG_OccupancyNumOfResonancesList_dict = tree()
    previous_TOCSY_AAIG = None
    matchingNOESYAAIG_occupancy_dict = {}
    Num_of_TOCSY_resonances = 0 # save here the number of TOCSY resonances for a particular AAIG
    for TOCSY_words_list in TOCSY_contents:
        try:
            # discard 1st column
            #TOCSY_AAIG=int(float(TOCSY_words_list[1]))    # residue i
            TOCSY_AAIG=str(TOCSY_words_list[0])    # residue i
            #print "DEBUG: TOCSY_AAIG=", TOCSY_AAIG
            if previous_TOCSY_AAIG != None and TOCSY_AAIG != previous_TOCSY_AAIG: # if we look at a different AAIG in TOCSY, print the matches and occupancies
                                                                                           # of the previous TOCSY AAIG and clear matchingNOESYAAIG_occupancy_dict
                for k,v in list(matchingNOESYAAIG_occupancy_dict.items()):
                    #print "residue i =",previous_TOCSY_AAIG," possible residue i-1 =",k," occupancy",str(v)+"/"+str(Num_of_TOCSY_resonances)
                    TOCSYAAIG_NOESYAAIG_OccupancyNumOfResonancesList_dict[previous_TOCSY_AAIG][k]=[v,Num_of_TOCSY_resonances]
                matchingNOESYAAIG_occupancy_dict = {}
                Num_of_TOCSY_resonances = 0
            TOCSY_H_resonance=float(TOCSY_words_list[1]) # aliphatic H resonance of residue i-1 
            TOCSY_C_resonance=float(TOCSY_words_list[2]) # aliphatic C (Ca,Cb,Cc,Cg,Ce,etc.) resonance of residue i-1; this C is covalently bonded to the above H
            matchingNOESYAAIG_list = []  # list with the NOESY index groups that match the current TOCSY index group
            Num_of_TOCSY_resonances += 1
            for q_index in range(0,len(NOESY_contents)):    # iterate over all NOESY resonances
                try:
                    NOESY_words_list = NOESY_contents[q_index]
                    
                    #NOESY_AAIG =int(float(NOESY_words_list[1]))  # residue (i-1)
                    NOESY_AAIG =str(NOESY_words_list[0])
                    NOESY_w2=float(NOESY_words_list[1]) # resonance of an aliphatic H close to residue (i-1)
                    NOESY_w3=float(NOESY_words_list[2]) # resonance of an aliphatic C close to residue (i-1); this C is covalently bonded to the above H
                    #print "DEBUG: NOESY_AAIG=",NOESY_AAIG,"NOESY_w2=",NOESY_w2,"NOESY_w3=",NOESY_w3
                    if ( (NOESY_w2-tolH) <= TOCSY_H_resonance <= (NOESY_w2+tolH) ) and ( (NOESY_w3-tolC) <= TOCSY_C_resonance <= (NOESY_w3+tolC) ):
                        #print "TOCSY AAIG ",TOCSY_AAIG," w2",TOCSY_H_resonance," w3",TOCSY_C_resonance," matching with NOESY AAIG ",NOESY_AAIG," w2",NOESY_w2," w3",NOESY_w3
                        if NOESY_AAIG in matchingNOESYAAIG_list:  # just a warning that this TOCSY peak mathces multiple peaks from the same NOESY AAIG
                                                                        # in this case matching will be counted only once, it won't make a difference
                            print("WARNING: another peak from NOESY AAIG ",NOESY_AAIG, "was already found matching this peak from TOCSY AAIG ",TOCSY_AAIG,"!!!!!")
                        else:
                            matchingNOESYAAIG_list.append(NOESY_AAIG)
                except (ValueError, IndexError):
                    #print "WARNING: the 4th and 5th elements of the following NOESY file line are not numbers:"
                    #print "NOESY file line:", NOESY_line
                    continue
            
            # populate matchingNOESYAAIG_list using the matches of the current TOCSY line
            for AAIG in matchingNOESYAAIG_list:
                try:
                    matchingNOESYAAIG_occupancy_dict[AAIG] += 1 # if it exist alread increment its occupancy
                except KeyError:
                    matchingNOESYAAIG_occupancy_dict[AAIG] = 1 # if it does not exist initialize occupancy
            previous_TOCSY_AAIG = TOCSY_AAIG
        except (ValueError, IndexError):
            print("WARNING: the 2nd, 3rd or 4th elements of the following TOCSY file line are not numbers:")
            print("TOCSY file line:", TOCSY_words_list)
            continue
    
    # save the possible connectivities of the last TOCSY AAIG
    for k,v in list(matchingNOESYAAIG_occupancy_dict.items()):
        #print "residue i =",previous_TOCSY_AAIG," possible residue i-1 =",k," occupancy",str(v)+"/"+str(Num_of_TOCSY_resonances)
        TOCSYAAIG_NOESYAAIG_OccupancyNumOfResonancesList_dict[previous_TOCSY_AAIG][k]=[v,Num_of_TOCSY_resonances]
    
    return TOCSYAAIG_NOESYAAIG_OccupancyNumOfResonancesList_dict