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

import sys, re, os
from operator import itemgetter
from collections import OrderedDict
from argparse import ArgumentParser
from tabulate import tabulate
from ete3 import Tree
from lib.alignment import *

## Parse command line arguments
def cmdlineparse():
    parser = ArgumentParser(description="command line arguments",
        epilog="EXAMPLE: chain_linker.py -nhmap results_summary.3mers_round1_rst -poolconfile connectivities_all.pool.3mers_round1_rst -multicon -mcincr 0.1 -mcmin 0.5 -patch")
    parser.add_argument("-nhmap", dest="ABSOLUTE_MATCHES_FILE", required=False,
                        help="file with absolute matches from previous run. (default: %(default)s)",
                        metavar="<absolute matches file>")
    parser.add_argument("-poolconfile", dest="POOL_CONNECTIVITIES_FILE", required=False,
                        help="pool connectivities file from the previous round. (default: %(default)s)",
                        metavar="<pool connectivities file>")
    parser.add_argument("-multicon", dest="ALLOW_MULTIPLE_CONNECTIVITIES", required=False, action='store_true', default=False,
                        help="allow linker extension even if multiple alternative connectivities (including the correct one) are present. (default: %(default)s)")
    parser.add_argument("-mcincr", dest="MCUTOFF_INCREMENT", required=False, default=1.0, type=float,
                        help="matching cutoff increment. (default: %(default)s)",
                        metavar="<mcutoff increment>")
    parser.add_argument("-mcmin", dest="MINIMUM_MCUTOFF", required=False, default=0.0, type=float,
                        help="minimum matching cutoff to be used. (default: %(default)s)",
                        metavar="<minimum mcutoff>")
    parser.add_argument("-patch", dest="EXTEND_N_TERMINALS", required=False, action='store_true', default=False,
                        help="place single TAAIG groups at the N-terminals if the respective connectivity exists. (USE WITH CAUTION!). (default: %(default)s)")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s version 1.2')
    args=parser.parse_args()
    return args


args = cmdlineparse()

##########################################################################################################
##                                         LOADING FILES                                                ##
##########################################################################################################

## READ FILE WITH ABSOLUTE MATCHES
absolute_AAIGmatches_alignment_list, protein_alignment_list = \
    Alignment.read_NHmap_file(args.ABSOLUTE_MATCHES_FILE, get_protein_alignment=True)


######################################################### ALL FUNCTION DEFINITIONS ###########################################################

## MODIFYING i_iminus1_dict ACCORDING TO absolute_AAIGmatches_alignment_list
def modify_i_iminus1_dict(i, iminus1):
    global i_iminus1_dict
    
    ## MODIFY i_iminus1_dict
    try:
        for quintuplet in i_iminus1_dict[i]:
            if quintuplet[0] == iminus1:
                print("Setting connectivity of ", i, " to ", quintuplet)
                i_iminus1_dict[i] = [quintuplet]
                break
    except KeyError:    # CLEVER: add a new connection that is not present in the connectivities file
        print("Adding missing connectivity of ", i, " to ", (iminus1, 1, 1, 1, 1))
        i_iminus1_dict[i] = [(iminus1, 1, 1, 1, 1)] # set occupancy and TOCSY resonance number arbitrarily to 1 and 1
    
    # remove iminus1 from the list of values of every other key
    i_list = list(i_iminus1_dict.keys())
    i_list.remove(i)
    for other_i in i_list:
        new_quintuplet_list = []
        for quintuplet in i_iminus1_dict[other_i]:
            if quintuplet[0] != iminus1:
                new_quintuplet_list.append(quintuplet)
        i_iminus1_dict[other_i] = new_quintuplet_list


def insert_linker_to_alignment(linker):
    
    global new_absolute_AAIGmatches_alignment_list
    
    starting_index = new_absolute_AAIGmatches_alignment_list.index(linker[0])
    ending_index = new_absolute_AAIGmatches_alignment_list.index(linker[-1])
    linker_index = 1
    for aln_index in range(starting_index+1, ending_index):
        new_absolute_AAIGmatches_alignment_list[aln_index] = linker[linker_index]
        linker_index += 1

def insert_single_TIG_to_alignment(linker):
    
    global new_absolute_AAIGmatches_alignment_list
    
    ending_index = new_absolute_AAIGmatches_alignment_list.index(linker[-1])
    starting_index = ending_index -1
    new_absolute_AAIGmatches_alignment_list[starting_index] = linker[0]


# NOW DEFINE THE NECESSARY FUNCTIONS TO BUILD THE CHAIN TREES
def build_Tree(starting_TAAIG):
    
    def populate_leaves(Assignment_Tree):
        """
            FUNCTION that adds new branches to the leaves of the Tree.
            ARGUMENTS:
            Assignment_Tree:    The Tree structure with connectivities
            RETURNS:
            (Assignment_Tree, BOOLEAN):    A tuple with elements the input Tree structure with new branches (if applicable), and a BOOLEAN value which is True if the function added
                                           new leaves to the Tree, or False otherwise
        """
        
        global i_iminus1_dict
        number_of_new_leaves = 0
        intersection, multi = None, None    # default values (if None then the connectivities are in the old format)
        for leaf, name in zip(Assignment_Tree.get_leaves(), Assignment_Tree.get_leaf_names()):
            try:
                #for child_tuple in i_iminus1_dict[name]:
                for child_quintuplet in i_iminus1_dict[name]: # WARNING: i_iminus1_normProbabilities_dict contains only matches above the Z-score cutoff
                    NOESYaaindex = child_quintuplet[0]
                    _occupancy = child_quintuplet[1]     # add prefix "_" to discriminate from new child feature "occupancy" that will be added
                    _numOfResonances = child_quintuplet[2]
                    if len(child_quintuplet) == 5:
                        intersection, multi = child_quintuplet[3], child_quintuplet[4]
                    ancestors_list = [ancestor.name for ancestor in leaf.get_ancestors()]
                    if NOESYaaindex in ancestors_list:     # if the current NOESY AAIG is already a node or leaf in the Tree, continue to the next
                        continue
                    #new_child = leaf.add_child(name=NOESYaaindex, dist=float(_occupancy)/_numOfResonances) # add a new brach to the current TOCSY add index (leaf) with length the ratio occupancy/numOfResonances
                    
                    new_child = leaf.add_child(name=NOESYaaindex) # add a new brach to the current TOCSY add index (leaf) with length the respective probability
                    
                    new_child.add_features(occupancy=_occupancy,
                                           numOfResonances=_numOfResonances,
                                           intersection=intersection,
                                           multi=multi)
                    number_of_new_leaves += 1
                    #print "DEBUG: adding connection: ",name,"-->",NOESYaaindex
            except KeyError:
                continue
        
        del leaf, name
        #print Assignment_Tree.get_ascii(show_internal=True, compact=False, attributes=["name", "dist", "occupancy", "numOfResonances"])
        #print Assignment_Tree.get_ascii(show_internal=True, compact=False)
        if number_of_new_leaves > 0:
            return True
        else:
            return False
    
    expand_tree = True
    Assignment_Tree = Tree()
    Root = Assignment_Tree.get_tree_root()
    Root.add_feature("name", starting_TAAIG)
    level = 1
    #sys.stdout.write("Expanding tree from level ")
    while expand_tree:
        #sys.stdout.write(str(level)+" ")
        #sys.stdout.flush()
        expand_tree = populate_leaves(Assignment_Tree)
        level += 1
        if level == linker_length:
            break
    if level < linker_length: # discard chains shorter than the cutoff
        return []   # return an empty chain_list
    # Print the Tree
    #print Assignment_Tree.get_ascii(show_internal=True, compact=False)
    #print Assignment_Tree.get_ascii(show_internal=True, compact=False, attributes=["name", "dist", "occupancy", "numOfResonances"])
    
    #print "\nSaving chains from Tree..."
    
    chain_list = []
    for leaf in Assignment_Tree.get_leaves():
        chain = []
        score = leaf.dist
        chain.append(leaf.name)
        for ancestor in leaf.get_ancestors():
            chain.append(ancestor.name)
        if len(chain) == linker_length:
            chain_list.append(chain)
        del chain
        del ancestor
        del leaf
    del Assignment_Tree
    
    return chain_list

##################################################################### END OF FUNCTION DEFINITIONS #############################################################


print("Loading Cutoff Connectivities from file", args.POOL_CONNECTIVITIES_FILE)
with open(args.POOL_CONNECTIVITIES_FILE, 'r') as f:
    connectivities_file_contents = f.readlines()

iterationNo = 0
for mcutoff in reversed(list(range(int(args.MINIMUM_MCUTOFF * 10), 10, int(args.MCUTOFF_INCREMENT * 10)))):
    mcutoff /= float(10)
    print("################################ SETTING MATCHING CUTOFF TO "+str(mcutoff)+" ################################")
    iterationNo += 1

    i_iminus1_dict = {} # a dictionary containing all possible connectivities of every TOCSY AAIG; has keys the TOCSY aa indices and values lists of quintuplets (tuples)
                        # consisting of (NOESYaaindex, occupancy, numOfResonances)
    
    for line in connectivities_file_contents[1:]:
        word_list = line.split()
        TOCSY_aaindex = word_list[0]
        i_iminus1_dict[TOCSY_aaindex] = []
        values_string = ''.join(word_list[1:])
        elements_string = re.sub('[\(\)\'\"]', '',  values_string ).split("),(")[0]
        elements_list = elements_string.split(",")
        print("DEBUG: elements_list=", elements_list)
        try:
            for aa, occupancy, TOCSY_resonnum, intersection, multi in zip(elements_list[0::5], elements_list[1::5], elements_list[2::5],
                                                                            elements_list[3::5], elements_list[4::5]):
                #print aa, occupancy, TOCSY_resonnum
                if float(occupancy)/float(TOCSY_resonnum) >= mcutoff:
                    i_iminus1_dict[TOCSY_aaindex].append((aa, int(occupancy), int(TOCSY_resonnum), float(intersection), float(multi)))
        except ValueError:  # if this is the old connectivity format
            for aa, occupancy, TOCSY_resonnum in zip(elements_list[0::3], elements_list[1::3], elements_list[2::3]):
                #print aa, occupancy, TOCSY_resonnum
                if float(occupancy)/float(TOCSY_resonnum) >= mcutoff:
                    i_iminus1_dict[TOCSY_aaindex].append((aa, int(occupancy), int(TOCSY_resonnum)))
    
    if iterationNo == 1:    # if this is the first iteration copy absolute_AAIGmatches_alignment_list to new_absolute_AAIGmatches_alignment_list, otherwise let
                            # new_absolute_AAIGmatches_alignment_list retain its contents from the previous iteration
        new_absolute_AAIGmatches_alignment_list = list(absolute_AAIGmatches_alignment_list)
    
    for position in range(1, len(new_absolute_AAIGmatches_alignment_list)):
        if new_absolute_AAIGmatches_alignment_list[position] != '-' and new_absolute_AAIGmatches_alignment_list[position-1] != '-':
            modify_i_iminus1_dict(new_absolute_AAIGmatches_alignment_list[position], new_absolute_AAIGmatches_alignment_list[position-1])
    
    
    ##########################################################################################################
    ##                                         ADDING LINKERS                                               ##
    ##########################################################################################################
    
    ##
    ##  ADDING LINKERS ALLOWING ONLY SINGLE CONNECTIVITIES PER TOCSY INDEX
    ##
    ADDED_LINKER = True
    while ADDED_LINKER:
        print("STARTING NEW CYCLE OF LINKER INSERTION.")
        ADDED_LINKER = False
        
        # FINDING THE INTACT CHAINS IN new_absolute_AAIGmatches_alignment_list
        absolute_match_chains_list = []
        chain = []
        for position in range(len(new_absolute_AAIGmatches_alignment_list)):
            TAAIG = new_absolute_AAIGmatches_alignment_list[position]
            if TAAIG != '-' and TAAIG != 'N/A':
                chain.append(TAAIG)
            else:
                if len(chain) > 0:  # don't save single Tocsy indices because they give no connectivity information
                    absolute_match_chains_list.append(chain)
                chain = []
        if len(chain) > 0:  # save the last chain
            absolute_match_chains_list.append(chain)
        
        # ADDING LINKERS
        for right_chain_index in reversed(list(range(1, len(absolute_match_chains_list)))):
            right_chain = absolute_match_chains_list[right_chain_index]
            left_chain = absolute_match_chains_list[right_chain_index-1]
            starting_TAAIG = right_chain[0]
            ending_TAAIG = left_chain[-1]
            # THESE ARE TO CHECK IF THE TWO CHAINS HAVE A PROLINE BETWEEN THEM
            Cterm_index = new_absolute_AAIGmatches_alignment_list.index(ending_TAAIG)
            Nterm_index = new_absolute_AAIGmatches_alignment_list.index(starting_TAAIG)
            if 'P' in protein_alignment_list[Cterm_index:Nterm_index+1]:    # IF THE IS A PROLINE, MOVE TO THE NEXT 2 CHAINS
                continue
                
            linker = [starting_TAAIG]  # create a linker list that will start (left -> right) from the ending_TAAIG and end at starting_TAAIG
            
            #if not starting_TAAIG in i_iminus1_dict.keys() or len(i_iminus1_dict[starting_TAAIG]) > 1:
            #    continue
            #else:
            #    quintuplet = i_iminus1_dict[starting_TAAIG][0]
            #    linker.append(quintuplet[0]) # insert to the linker list the TOCSY index that is connected to the N-term of the right chain
            
            while 1:
                if linker[0] == ending_TAAIG:
                    break
                if not linker[0] in list(i_iminus1_dict.keys()) or not len(i_iminus1_dict[linker[0]]) == 1:   # takes care also for no connectivities
                    #linker = []
                    break
                if len(linker) != len(set(linker)): # if the linker contains duplicate Tindices exit
                    break
                else:
                    #print "DEBUG: i_iminus1_dict[linker[0]][0]=", i_iminus1_dict[linker[0]]
                    quintuplet = i_iminus1_dict[linker[0]][0]
                    linker.insert(0, quintuplet[0]) # insert to the linker list the TOCSY index that is connected to the N-term of the linker
            
            #print "DEBUG: ending_TAAIG=", ending_TAAIG, "linker", linker[:], "starting_TAAIG=", starting_TAAIG
            if linker[0] == ending_TAAIG and linker[-1] == starting_TAAIG and len(linker[1:-1]) == len(new_absolute_AAIGmatches_alignment_list[Nterm_index+1:Cterm_index]):
                print("Adding the linker ",linker[1:-1]," into the absolute matches alignment")
                ADDED_LINKER = True
                # FILL IN THE GAPS IN absolute_match_chains_list
                insert_linker_to_alignment(linker)
                
                # CLEAN i_iminus1_dict CONNECTIVITIES FROM INSERTED TOCSY INDICES AND KEEP ONLY ONE CONNECTIVITY FOR THE TOCSY INDICES OF THE LINKER
                for lindex, link_TAAIG in enumerate(linker[1:]):   # ommit the N-terminal TAAIG
                        for quintuplet in i_iminus1_dict[link_TAAIG]:
                            if quintuplet[0] == linker[lindex-1][0]:   # keep only the connectivity that agrees with the linker
                                i_iminus1_dict[link_TAAIG] = [quintuplet]
                                break
                        
                for i in list(i_iminus1_dict.keys()):
                    new_quintuplet_set = set()
                    for link_TAAIG in enumerate(linker[1:]):   # ommit the N-term and C-term that 
                        if link_TAAIG == i:    # skip the connectivities of the current link, as they have been modified above
                            continue
                        for quintuplet in i_iminus1_dict[i]:   # clean all other connectivities from the link
                            if link_TAAIG != quintuplet[0]:
                                new_quintuplet_set.add(quintuplet)
                    i_iminus1_dict[i] = list(new_quintuplet_set)
                
                
        
    
    ##
    ##  ADDING LINKERS ALLOWING MULTIPLE CONNECTIVITIES PER TOCSY INDEX
    ##
    if args.ALLOW_MULTIPLE_CONNECTIVITIES:
        ADDED_LINKER = True
        while ADDED_LINKER:
            print("STARTING NEW CYCLE OF LINKER INSERTION ALLOWING MULTIPLE CONNECTIVITIES PER TOCSY GROUP.")
            ADDED_LINKER = False
            
            # FINDING THE INTACT CHAINS IN new_absolute_AAIGmatches_alignment_list
            absolute_match_chains_list = []
            chain = []
            for position in range(len(new_absolute_AAIGmatches_alignment_list)):
                TAAIG = new_absolute_AAIGmatches_alignment_list[position]
                if TAAIG != '-' and TAAIG != 'N/A':
                    chain.append(TAAIG)
                else:
                    if len(chain) > 0:  # don't save single Tocsy indices because they give no connectivity information
                        absolute_match_chains_list.append(chain)
                    chain = []
            if len(chain) > 0:  # save the last chain
                absolute_match_chains_list.append(chain)
            
            # ADDING LINKERS
            #print "DEBUG: absolute_match_chains_list=", absolute_match_chains_list
            added_linker_list = []
            for right_chain_index in reversed(list(range(1, len(absolute_match_chains_list)))):
                right_chain = absolute_match_chains_list[right_chain_index]
                left_chain = absolute_match_chains_list[right_chain_index-1]
                starting_TAAIG = right_chain[0]
                ending_TAAIG = left_chain[-1]
                #print "DEBUG: left_chain=", left_chain
                #print "DEBUG: right_chain=", right_chain
                # THESE ARE TO CHECK IF THE TWO CHAINS HAVE A PROLINE BETWEEN THEM
                Cterm_index = new_absolute_AAIGmatches_alignment_list.index(ending_TAAIG)
                Nterm_index = new_absolute_AAIGmatches_alignment_list.index(starting_TAAIG)
                if 'P' in protein_alignment_list[Cterm_index:Nterm_index+1]:    # IF THERE IS A PROLINE, MOVE TO THE NEXT 2 CHAINS
                    continue
                
                linker_length = new_absolute_AAIGmatches_alignment_list.index(starting_TAAIG) - new_absolute_AAIGmatches_alignment_list.index(ending_TAAIG) + 1
                
                #print "Building Tree starting from amino acid index", starting_TAAIG, ", looking for chains of length ", linker_length
                potential_linker_list = build_Tree(starting_TAAIG)
                
                correct_linker_list = []
                for potential_linker in potential_linker_list:
                    if potential_linker[-1] == starting_TAAIG and potential_linker[0] == ending_TAAIG:
                        if len(potential_linker) != len(set(potential_linker)): # if the linker contains duplicate Tindices, skip it
                            continue
                        correct_linker_list.append(potential_linker)
                #print "DEBUG: correct_linker_list=", correct_linker_list
                if len(correct_linker_list) == 1:   # if only one correct linker was found, then modify the new_absolute_AAIGmatches_alignment_list
                    linker = correct_linker_list[0]
                    # IF ANY OF THE TOCSY INDICES OF THE LINKER IS FOUND MULTIPLE TIMES IN THE ALIGNMENT, DO NOT ADD THE LINKER
                    ADD_THE_LINKER = True
                    for linker_TAAIG in linker:
                        if new_absolute_AAIGmatches_alignment_list.count(linker_TAAIG) > 1:
                            ADD_THE_LINKER = False
                            break
                    # IF ANY OF THE TOCSY INDICES OF THE LINKER IS 
                    if ADD_THE_LINKER and len(linker[1:-1]) > 0:
                        ADDED_LINKER = True
                        added_linker_list.append(linker)
            
            print("DEBUG: added_linker_list=", added_linker_list)
            # REMOVE LINKERS WITH COMMON TOCSY INDEX GROUPS
            all_bridges_TIGs_list=[]
            for linker in added_linker_list:
                all_bridges_TIGs_list.extend(linker[1:-1])
            
            TIGs2remove_set = set() # set with duplicate TOCSY index groups in linkers
            for TIG in all_bridges_TIGs_list:
                if all_bridges_TIGs_list.count(TIG) > 1:
                    TIGs2remove_set.add(TIG)
            
            unique_added_linker_list = []   # linkers with non-common TOCSY index groups
            for linker in added_linker_list:
                KEEP_IT = True
                for TIG in linker[1:-1]:
                    if TIG in TIGs2remove_set:
                        KEEP_IT = False
                        break
                if KEEP_IT == True:
                    unique_added_linker_list.append(linker)
            print("DEBUG: unique_added_linker_list=", unique_added_linker_list)
            
            if len(unique_added_linker_list) == 0:  # IF NO UNIQUE LINKER HAS BEEN ADDED, BREAK THE WHILE LOOP
                break
            
            for linker in unique_added_linker_list:
                # FILL IN THE GAPS IN absolute_match_chains_list
                print("Adding the linker ", linker[1:-1], " into the absolute matches alignment")
                insert_linker_to_alignment(linker)
                # CLEAN i_iminus1_dict CONNECTIVITIES FROM INSERTED TOCSY INDICES AND KEEP ONLY ONE CONNECTIVITY FOR THE TOCSY INDICES OF THE LINKER
                for lindex, link_TAAIG in enumerate(linker[1:]):   # ommit the N-terminal TAAIG
                        for quintuplet in i_iminus1_dict[link_TAAIG]:
                            if quintuplet[0] == linker[lindex-1][0]:   # keep only the connectivity that agrees with the linker
                                i_iminus1_dict[link_TAAIG] = [quintuplet]
                                break
                        
                for i in list(i_iminus1_dict.keys()):
                    new_quintuplet_set = set()
                    for link_TAAIG in enumerate(linker[1:]):   # ommit the N-term and C-term that 
                        if link_TAAIG == i:    # skip the connectivities of the current link, as they have been modified above
                            continue
                        for quintuplet in i_iminus1_dict[i]:   # clean all other connectivities from the link
                            if link_TAAIG != quintuplet[0]:
                                new_quintuplet_set.add(quintuplet)
                    i_iminus1_dict[i] = list(new_quintuplet_set)
                for k,v in list(i_iminus1_dict.items()):
                    print(k, v)
    
    
    #print "DEBUG: new_absolute_AAIGmatches_alignment_list after =", new_absolute_AAIGmatches_alignment_list


##########################################################################################################
# IF REQUSTED, PLACE SINGLE TAAIG GROUPS USING ONLY THEIR CONNECTIVITY WITH THE HANGING N-TERMINAL ENDS #
##########################################################################################################
patched_TIGs_list = []   # residue that have NOESY but not TOCSY and are flanked by two residue with both TOCSY and NOESY
if args.EXTEND_N_TERMINALS:
    
    iterationNo = 0
    for mcutoff in reversed(list(range(int(args.MINIMUM_MCUTOFF * 10), 10, int(args.MCUTOFF_INCREMENT * 10)))):
        mcutoff /= float(10)
        print("################################ SETTING MATCHING CUTOFF TO "+str(mcutoff)+" ################################")
        print("DEBUG: new_absolute_AAIGmatches_alignment_list after =", new_absolute_AAIGmatches_alignment_list)
        iterationNo += 1
    
        i_iminus1_dict = {} # a dictionary containing all possible connectivities of every TOCSY AAIG; has keys the TOCSY aa indices and values lists of quintuplets (tuples)
                            # consisting of (NOESYaaindex, occupancy, numOfResonances)
        
        for line in connectivities_file_contents[1:]:
            word_list = line.split()
            TOCSY_aaindex = word_list[0]
            i_iminus1_dict[TOCSY_aaindex] = []
            values_string = ''.join(word_list[1:])
            elements_string = re.sub('[\(\)\'\"]', '',  values_string ).split("),(")[0]
            elements_list = elements_string.split(",")
            try:
                for aa, occupancy, TOCSY_resonnum, intersection, multi in zip(elements_list[0::5], elements_list[1::5],
                                                                              elements_list[2::5],
                                                                              elements_list[3::5], elements_list[4::5]):
                    # print aa, occupancy, TOCSY_resonnum
                    if float(occupancy)/float(TOCSY_resonnum) >= mcutoff:
                        i_iminus1_dict[TOCSY_aaindex].append(
                            (aa, int(occupancy), int(TOCSY_resonnum), float(intersection), float(multi)))
            except IndexError:  # if this is the old connectivity format
                for aa, occupancy, TOCSY_resonnum in zip(elements_list[0::3], elements_list[1::3], elements_list[2::3]):
                    # print aa, occupancy, TOCSY_resonnum
                    if float(occupancy)/float(TOCSY_resonnum) >= mcutoff:
                        i_iminus1_dict[TOCSY_aaindex].append((aa, int(occupancy), int(TOCSY_resonnum)))
        
        #if iterationNo == 1:    # if this is the first iteration copy absolute_AAIGmatches_alignment_list to new_absolute_AAIGmatches_alignment_list, otherwise let
        #                        # new_absolute_AAIGmatches_alignment_list retain its contents from the previous iteration
        #    new_absolute_AAIGmatches_alignment_list = list(absolute_AAIGmatches_alignment_list)
        
        for position in range(1, len(new_absolute_AAIGmatches_alignment_list)):
            if new_absolute_AAIGmatches_alignment_list[position] != '-' and new_absolute_AAIGmatches_alignment_list[position-1] != '-':
                modify_i_iminus1_dict(new_absolute_AAIGmatches_alignment_list[position], new_absolute_AAIGmatches_alignment_list[position-1])
        
        ##
        ##  ADDING SINGLE LINKERS ALLOWING ONLY SINGLE CONNECTIVITIES PER TOCSY INDEX
        ##
        print("DEBUG: new_absolute_AAIGmatches_alignment_list after =", new_absolute_AAIGmatches_alignment_list)
        print("DEBUG: i_iminus1_dict=", i_iminus1_dict)
        ADDED_LINKER = True
        while ADDED_LINKER:
            print("STARTING NEW CYCLE OF N-TERMINAL ENDS PATCHING.")
            ADDED_LINKER = False
            
            # FINDING THE INTACT CONTINGS IN new_absolute_AAIGmatches_alignment_list
            absolute_match_chains_list = []
            chain = []
            for position in range(len(new_absolute_AAIGmatches_alignment_list)):
                TAAIG = new_absolute_AAIGmatches_alignment_list[position]
                if TAAIG != '-' and TAAIG != 'N/A':
                    chain.append(TAAIG)
                else:
                    if len(chain) > 0:  # don't save single Tocsy indices because they give no connectivity information
                        absolute_match_chains_list.append(chain)
                    chain = []
            if len(chain) > 0:  # save the last chain
                absolute_match_chains_list.append(chain)
            
            # ADDING LINKERS
            for right_chain_index in reversed(list(range(1, len(absolute_match_chains_list)))):
                right_chain = absolute_match_chains_list[right_chain_index]
                left_chain = absolute_match_chains_list[right_chain_index-1]
                starting_TAAIG = right_chain[0]
                ending_TAAIG = left_chain[-1]
                ## THESE ARE TO CHECK IF THE TWO CHAINS HAVE A PROLINE BETWEEN THEM
                Cterm_index = new_absolute_AAIGmatches_alignment_list.index(ending_TAAIG)
                Nterm_index = new_absolute_AAIGmatches_alignment_list.index(starting_TAAIG)
                print("DEBUG: protein_alignment_list[Nterm_index-1:Nterm_index+1] = ", protein_alignment_list[Nterm_index-1:Nterm_index+1])
                if 'P' in protein_alignment_list[Nterm_index-1:Nterm_index+1]:    # IF THERE IS A PROLINE, MOVE TO THE NEXT 2 CHAINS
                    continue
                    
                linker = [starting_TAAIG]  # create a linker list that will start (left -> right) from the ending_TAAIG and end at starting_TAAIG
                
                ##if not starting_TAAIG in i_iminus1_dict.keys() or len(i_iminus1_dict[starting_TAAIG]) > 1:
                ##    continue
                ##else:
                ##    quintuplet = i_iminus1_dict[starting_TAAIG][0]
                ##    linker.append(quintuplet[0]) # insert to the linker list the TOCSY index that is connected to the N-term of the right chain
                
                while 1:
                    if linker[0] == ending_TAAIG:
                        break
                    if not linker[0] in list(i_iminus1_dict.keys()) or not len(i_iminus1_dict[linker[0]]) == 1:   # takes care also for no connectivities
                        #linker = []
                        break
                    if len(linker) != len(set(linker)): # if the linker contains duplicate Tindices exit
                        break
                    else:
                        #print "DEBUG: i_iminus1_dict[linker[0]][0]=", i_iminus1_dict[linker[0]]
                        quintuplet = i_iminus1_dict[linker[0]][0]
                        linker.insert(0, quintuplet[0]) # insert to the linker list the TOCSY index that is connected to the N-term of the linker and break
                        break
                
                print("DEBUG: ending_TAAIG=",  ending_TAAIG, "linker", linker[:], "starting_TAAIG=", starting_TAAIG)
                if starting_TAAIG not in list(i_iminus1_dict.keys()) or linker[0] in new_absolute_AAIGmatches_alignment_list:
                    continue
                patch_quintuplet = [t for t in i_iminus1_dict[linker[-1]] if t[0]==linker[0]]
                print("DEBUG: linker=", linker, "patch_quintuplet=", patch_quintuplet)
                if linker[-1] == starting_TAAIG and len(linker) == 2 and float(patch_quintuplet[0][1])/float(patch_quintuplet[0][2]) >= mcutoff:  # CHANGE ME <== the mcutoff for patching
                    print("Patching ",linker[0]," after ", linker[1], " into the absolute matches alignment")
                    ADDED_LINKER = True
                    # FILL IN THE GAPS IN absolute_match_chains_list
                    insert_single_TIG_to_alignment(linker)
                    patched_TIGs_list.append(linker[0])
                
                    # CLEAN i_iminus1_dict CONNECTIVITIES FROM INSERTED TOCSY INDICES AND KEEP ONLY ONE CONNECTIVITY FOR THE TOCSY INDICES OF THE LINKER
                    print("DEBUG: i_iminus1_dict=", i_iminus1_dict)
                    for lindex, link_TAAIG in enumerate(linker[1:]):   # ommit the N-terminal TAAIG
                            for quintuplet in i_iminus1_dict[link_TAAIG]:
                                if quintuplet[0] == linker[lindex-1][0]:   # keep only the connectivity that agrees with the linker
                                    i_iminus1_dict[link_TAAIG] = [quintuplet]
                                    break
                            
                    for i in list(i_iminus1_dict.keys()): # iterate over all TAAIG groups
                        new_quintuplet_set = set()
                        for link_TAAIG in enumerate(linker):
                            if link_TAAIG == i:    # skip the connectivities of the current link, as they have been modified above
                                continue
                            for quintuplet in i_iminus1_dict[i]:   # clean all other connectivities from the link
                                if link_TAAIG != quintuplet[0]:
                                    new_quintuplet_set.add(quintuplet)
                        i_iminus1_dict[i] = list(new_quintuplet_set)

## WRITING new_absolute_AAIGmatches_alignment_list TO FILE
if args.EXTEND_N_TERMINALS == True:
    outfname = args.ABSOLUTE_MATCHES_FILE+".chainlinkers.patch"
else:
    outfname = args.ABSOLUTE_MATCHES_FILE+".chainlinkers"

Alignment("", 1).write_alignment_to_file(outfname=outfname,
                                         protein_alignment_list=protein_alignment_list,
                                         new_absolute_AAIGmatches_alignment_list=new_absolute_AAIGmatches_alignment_list,
                                         i_iminus1_dict=i_iminus1_dict,
                                         patched_TIGs_list=patched_TIGs_list)
# # CREATE THE CONNECTIVITIES ALIGNMENT
# print "DEBUG: i_iminus1_dict=", i_iminus1_dict
# connectivities_alignment = []
# if i_iminus1_dict != None:
#     for i in range(len(new_absolute_AAIGmatches_alignment_list)-1):
#         TIG_i = new_absolute_matches_alignment_list[i]
#         TIG_iplus1 = new_absolute_matches_alignment_list[i+1]
#         if TIG_i != '-' and TIG_iplus1 != '-':
#             try:
#                 connectivity = [c for c in i_iminus1_dict[TIG_iplus1] if c[0]==TIG_i][0] # connectivity = (TIG_i, occupancy, TOCSY_resonnum)
#                 print "DEBUG: connectivity=", connectivity
#                 if TIG_iplus1 in patched_TIGs_list: # if it has been patched, write "0/0" to distinguish it
#                     connectivities_alignment.append("0/0")
#                 else:
#                     connectivities_alignment.append(str(connectivity[1])+"/"+str(connectivity[2]))
#             except (IndexError, KeyError):
#                 connectivities_alignment.append('-')
#         else:
#             connectivities_alignment.append('-')
#     connectivities_alignment.append('-')
# connectivities_alignment.insert(0, '-')
# connectivities_alignment[0] = 'N/A'
#
# with open(outfname, 'w') as f:
#     f.write("\n\n\nABSOLUTE MATCHES WITH CHAIN LINKERS:\n\n\n" + "\n")
#     print "\n\n\nABSOLUTE MATCHES WITH CHAIN LINKERS:\n\n\n"
#     #print "DEBUG: new_absolute_matches_alignment_list=", new_absolute_matches_alignment_list
#     #print "DEBUG: protein_alignment_list=", protein_alignment_list
#     for start, end in zip(range(0,len(new_absolute_matches_alignment_list),25), range(25,len(new_absolute_matches_alignment_list),25)):
#         headers = protein_alignment_list[start:end]
#         table = [ new_absolute_matches_alignment_list[start:end] ]
#         table.append( connectivities_alignment[start:end] )
#         print tabulate(table, headers, tablefmt="fancy_grid").encode('utf8', 'strict')
#         f.write(tabulate(table, headers, tablefmt="fancy_grid").encode('utf8', 'strict') + "\n")
#         print "\n\n"
#         f.write("\n\n" + "\n")
#     headers = protein_alignment_list[end:]
#     table = [ new_absolute_matches_alignment_list[end:] ]
#     table.append( connectivities_alignment[end:] )
#     print tabulate(table, headers, tablefmt="fancy_grid").encode('utf8', 'strict')
#     f.write(tabulate(table, headers, tablefmt="fancy_grid").encode('utf8', 'strict') + "\n")