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



import sys, re, os
import numpy as np
from scipy.stats.mstats import zscore
from operator import itemgetter
from ordereddict import OrderedDict
from argparse import ArgumentParser
from ete3 import Tree
import gc

code3 = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
code1 = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
aa1to3_dict = dict((c1,c3) for c1,c3 in zip(code1,code3))
aa3to1_dict = dict((c3,c1) for c3,c1 in zip(code3,code1))

def approx_equal(x, y, tolerance=0.001):
    return abs(x-y) <= 0.5 * tolerance * (x + y)

def cmdlineparse():
    parser = ArgumentParser(description="command line arguments",
        epilog="EXAMPLES: \
               If you run the script for the 1st time use the complete connectivities and amino acid type prediction files as pool files:\n \
        ../modify_connectivities_and_aatypes.py -rstfile consensus_alignment.overlapped_chains_common_sequence.iteration7.6mers -poolconfile connectivities_all -allconfile connectivities_all -poolaafile amino_acid_type_prediction_probabilities -allaafile amino_acid_type_prediction_probabilities \
        If you have alread created pool files from previous runs do:\n \
        ../modify_connectivities_and_aatypes.py -rstfile consensus_alignment.overlapped_chains_common_sequence.iteration7.6mers -poolconfile connectivities_all.pool -allconfile connectivities_all -poolaafile amino_acid_type_prediction_probabilities.pool -allaafile amino_acid_type_prediction_probabilities")
    parser.add_argument("-rstfile", dest="ABSOLUTE_MATCHES_FILE", required=True, help="file with absolute matches from previous run to be used as restraints", metavar="<absolute matches restraint file>")
    parser.add_argument("-poolconfile", dest="POOL_CONNECTIVITIES_FILE", required=True, help="pool connectivities file; if specified connectivity calculation from input files will be skipped", metavar="<pool connectivities file>")
    parser.add_argument("-allconfile", dest="COMPLETE_CONNECTIVITIES_FILE", required=True, help="all connectivities file; necessary if -confile specified in order to calculate correct probabilities", metavar="<all connectivities file>")
    parser.add_argument("-poolaafile", dest="POOL_AA_TYPES_FILE", required=True, help="pool amino acid assignment file; if specified amino acid assignment calculation from input files will be skipped", metavar="<pool amino acid assignment file>")
    parser.add_argument("-allaafile", dest="COMPLETE_AA_TYPES_FILE", required=True, help="all amino acid assignment file; necessary if -aafile specified in order to calculate correct probabilities", metavar="<all amino acid assignment file>")
    parser.add_argument("-keepgly", dest="KEEP_ONLY_GLY", required=False, action='store_false', help="if GLY is the top ranked prediction with a probability at least 2 orders of magnitude greater than the 2nd prediction, keep only GLY")
    parser.add_argument("-keepala", dest="KEEP_ONLY_ALA", required=False, action='store_false', help="if ALA is the top ranked prediction with a probability at least 2 orders of magnitude greater than the 2nd prediction, keep only ALA")
    parser.add_argument("-nocorrect", dest="SELF_CORRECTION", required=False, action='store_false', default=True,
                        help="do not do self-correction")
    parser.add_argument("-addconn", dest="ADD_MISSING_CONNECTIVITIES", required=False, action='store_true', default=False,
                        help="add missing connectivities")
    parser.add_argument("-bigchains", dest="BIG_CHAIN_CUTOFF", required=False, type=int, default=None,
                        help="keep only connectivities that form chains of length equal or greater of the specified value. This options is useful when you have a big \
                        protein and you want to reduce the number of possible peptides.", metavar="<big chain cutoff length>")
    parser.add_argument("-mcutoff", dest="RESONANCE_MATCH_CUTOFF", required=False, type=float, default=0.8, help="number 0.0-1.0 saying how much of TOCSY resonances should match in the NOESY in order to consider it a possible match", metavar="<resonance match cutoff>")
    parser.add_argument("-zmcutoff", dest="ZSCORE_RESONANCE_MATCH_CUTOFF", required=False, type=float, default=0.0, help="a real number specifying the lower Z-score for a resonance match to be retained for chain building", metavar="<Z-score resonance match cutoff>")
    parser.add_argument("-maxocc", dest="MAXIMUM_OCCUPANCY_TOLERANCE", required=False, type=int, default=None,
                       help="keep only the connectivities of a TOCSY index group with the highest occupancy or those that differ from the maximum by -maxocc number (recommended: 0).",
                       metavar="<maximum occupancy tolerance>")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s version 0.8')
    args=parser.parse_args()
    return args

args = cmdlineparse()
if not args.KEEP_ONLY_GLY:
    args.KEEP_ONLY_GLY = True
if not args.KEEP_ONLY_ALA:
    args.KEEP_ONLY_ALA = True

protein_alignment_list = []
absolute_matches_alignment_list = []
if args.BIG_CHAIN_CUTOFF == None and args.MAXIMUM_OCCUPANCY_TOLERANCE == None:
    with open(args.ABSOLUTE_MATCHES_FILE, 'r') as f:
        for line in f:
            if line[0:53] == "CONSENSUS CONFIDENCE TOP-SCORED AND ABSOLUTE MATCHES:" or line[0:36] == "ABSOLUTE MATCHES WITH CHAIN LINKERS:":
                try:
                    while True:
                        line = f.next()
                        while line[0] == '\n':
                            line = f.next()
                        line = f.next()
                        protein_alignment_list.extend(line.replace('\xe2\x94\x82', '').split())
                        line = f.next()
                        line = f.next()
                        absolute_matches_alignment_list.extend(line.replace('\xe2\x94\x82', '').replace('*', '-').split())
                        line = f.next()
                        line = f.next()
                        line = f.next()
                        line = f.next()
                except StopIteration:
                    break
    


print "Loading Pool Connectivities from file", args.POOL_CONNECTIVITIES_FILE
i_iminus1_pool_dict = {} # a dictionary containing all possible connectivities of every TOCSY aa index; has keys the TOCSY aa indices and values lists of triplets (tuples) \
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

print "Loading All Connectivities from file", args.COMPLETE_CONNECTIVITIES_FILE
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


if args.MAXIMUM_OCCUPANCY_TOLERANCE == None:
    print "Loading Cutoff Amino Acid Prediction file ..."
    iaaindex_iminus1aaTypesProbPoolTupleList_dict = OrderedDict() # ordereddict with keys the aa index of residue i and values lists of tuples of the form
    with open(args.POOL_AA_TYPES_FILE, 'r') as f:
        aa_type_file_contents = f.readlines()
        for line in aa_type_file_contents[1:]:
            word_list = re.sub('---> \(i-1\) aa type', '', line).split()
            key = word_list[0]
            iaaindex_iminus1aaTypesProbPoolTupleList_dict[key] = []
            values_string = ''.join(word_list[1:])
            elements_string = re.sub('[\(\)\[\]\'\"]', '',  values_string ).split("),(")[0]
            elements_list = elements_string.split(",")
            for aa, Zscore in zip(elements_list[0::2], elements_list[1::2]):
                duplet = (aa, float(Zscore))
                iaaindex_iminus1aaTypesProbPoolTupleList_dict[key].append(duplet)
    
    print "Loading All Amino acid Prediction file ..."
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



if args.BIG_CHAIN_CUTOFF:   # KEEP ONLY CONNECTIVITIES THAT FORM CHAINS OF LENGTH EQUAL OR GREATER OF THE SPECIFIED VALUE.

    
    i2remove_set = set()
    for i in i_iminus1_pool_dict.keys():
        triplet_list = i_iminus1_pool_dict[i]
        new_triplet_list = []
        for triplet in triplet_list:
            if float(triplet[1])/triplet[2] >= args.RESONANCE_MATCH_CUTOFF:
                new_triplet_list.append(triplet)
        if len(new_triplet_list) == 0:
            i2remove_set.add(i)
            continue
        else:
            i_iminus1_pool_dict[i] = new_triplet_list
    for i in i2remove_set:
        del i_iminus1_pool_dict[i]
    
    new_i_iminus1_pool_dict = {}
    for i in i_iminus1_pool_dict.keys():
        prob_list = [float(triplet[1])/triplet[2] for triplet in i_iminus1_pool_dict[i]]
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
        new_i_iminus1_pool_dict[i] = []  # add i to the new dictionary
        for triplet, Zscore in zip(i_iminus1_pool_dict[i], zscore_array):
            if Zscore >= args.ZSCORE_RESONANCE_MATCH_CUTOFF or approx_equal(Zscore, args.ZSCORE_RESONANCE_MATCH_CUTOFF):
                new_i_iminus1_pool_dict[i].append(triplet)
    
    del i_iminus1_pool_dict
    i_iminus1_pool_dict = new_i_iminus1_pool_dict
    del new_i_iminus1_pool_dict
    
    
            
    
    Tindices2keep_set = set()    # a list of lists containing the possible connectivities and the overall score as the last element
    
    def build_Tree():
        global Tindices2keep_set
        
        def populate_leaves(Assignment_Tree):
            """
                FUNCTION that adds new branches to the leaves of the Tree.
                ARGUMENTS:
                Assignment_Tree:    The Tree structure with connectivities
                RETURNS:
                (Assignment_Tree, BOOLEAN):    A tuple with elements the input Tree structure with new branches (if applicable), and a BOOLEAN value which is True if the function added
                                               new leaves to the Tree, or False otherwise
            """
            
            global i_iminus1_pool_dict
            number_of_new_leaves = 0
            for leaf, name in zip(Assignment_Tree.iter_leaves(), Assignment_Tree.iter_leaf_names()):
                try:
                    for child_triplet in i_iminus1_pool_dict[name]: # WARNING: i_iminus1_normProbabilities_dict contains only matches above the Z-score cutoff
                        NOESYaaindex = child_triplet[0]
                        _occupancy = child_triplet[1]     # add prefix "_" to discriminate from new child feature "occupancy" that will be added
                        _numOfResonances = child_triplet[2]
                        ancestors_list = [ancestor.name for ancestor in leaf.get_ancestors()]
                        if NOESYaaindex in ancestors_list:     # if the current NOESY aa index is already a node or leaf in the Tree, continue to the next
                            continue
                        
                        new_child = leaf.add_child(name=NOESYaaindex) # add a new brach to the current TOCSY add index (leaf) with length the respective probability
                        
                        new_child.add_features(occupancy=_occupancy, numOfResonances=_numOfResonances)
                        number_of_new_leaves += 1
                except KeyError:
                    continue
            
            del leaf, name
            if number_of_new_leaves > 0:
                return True
            else:
                return False
        
        expand_tree = True
        Assignment_Tree = Tree()
        Root = Assignment_Tree.get_tree_root()
        Root.add_feature("name", i)
        level = 1
        sys.stdout.write("Expanding tree from level ")
        while expand_tree:
            sys.stdout.write(str(level)+" ")
            sys.stdout.flush()
            expand_tree = populate_leaves(Assignment_Tree)
            level += 1
            if level == args.BIG_CHAIN_CUTOFF:
                break
        if level < args.BIG_CHAIN_CUTOFF: # discard chains shorter than the cutoff
            return
        
        print "\nSaving chains from Tree..."
    
        for leaf in Assignment_Tree.iter_leaves():
            chain = []
            score = leaf.dist
            chain.append(leaf.name)
            for ancestor in leaf.get_ancestors():
                chain.append(ancestor.name)
            if len(chain) >= args.BIG_CHAIN_CUTOFF:
                for Tindex in chain:
                    Tindices2keep_set.add(Tindex)
            del chain
            del ancestor
            del leaf
        del Assignment_Tree
        gc.collect()
    
    
    N = len(i_iminus1_pool_dict.keys())
    index = 1
    for i in i_iminus1_pool_dict.keys():
        print "Building Tree starting from amino acid index",i,"("+str(index)+"/"+str(N)+") ..."
        build_Tree()
        index += 1
    
    
    if len(Tindices2keep_set) == len(i_iminus1_pool_dict.keys()):
        print "NO TOCSY INDEX GROUPS WERE REMOVED FROM CONNECTIVITY FILE !!!"
        sys.exit(0)
    for Tindex in i_iminus1_pool_dict.keys():
        if Tindex not in Tindices2keep_set:
            print "Deleting TOCSY index group", Tindex, " from connectivities pool file."
            del i_iminus1_pool_dict[Tindex]
            del i_iminus1_complete_dict[Tindex]
            del iaaindex_iminus1aaTypesProbTupleList_dict[Tindex]
            del iaaindex_iminus1aaTypesProbPoolTupleList_dict[Tindex]
    
    for Tindex in iaaindex_iminus1aaTypesProbPoolTupleList_dict.keys():
        if args.KEEP_ONLY_GLY:
            duplet_list = iaaindex_iminus1aaTypesProbPoolTupleList_dict[Tindex]
            if len(duplet_list) > 1 and duplet_list[0][0] == 'GLY' and duplet_list[0][1] >= 100 * duplet_list[1][1]:
                print "Setting amino acid type of TOCSY index group ", Tindex, " to ", [duplet_list[0]]
                iaaindex_iminus1aaTypesProbPoolTupleList_dict[Tindex] = [duplet_list[0]]
        if args.KEEP_ONLY_ALA:
            duplet_list = iaaindex_iminus1aaTypesProbPoolTupleList_dict[Tindex]
            if len(duplet_list) > 1 and duplet_list[0][0] == 'ALA' and duplet_list[0][1] >= 100 * duplet_list[1][1]:
                print "Setting amino acid type of TOCSY index group ", Tindex, " to ", [duplet_list[0]]
                iaaindex_iminus1aaTypesProbPoolTupleList_dict[Tindex] = [duplet_list[0]]
    
elif args.MAXIMUM_OCCUPANCY_TOLERANCE != None:  # KEEP ALL CONNECTIVITIES (not only those forming long chains)

    Tindex_maxOccupancy_dict = {}
    for triplet_list in i_iminus1_pool_dict.values():
        for triplet in triplet_list:
            Tindex = triplet[0]
            occupancy = int(triplet[1])
            if Tindex not in Tindex_maxOccupancy_dict.keys():
                Tindex_maxOccupancy_dict[Tindex] = occupancy
            elif Tindex_maxOccupancy_dict[Tindex] < occupancy:
                Tindex_maxOccupancy_dict[Tindex] = occupancy
    
    new_i_iminus1_pool_dict = {}
    for i in i_iminus1_pool_dict.keys():
        new_i_iminus1_pool_dict[i] = []  # add i to the new dictionary
        for triplet in i_iminus1_pool_dict[i]:
            Tindex = triplet[0]
            occupancy = triplet[1]
            if occupancy >= Tindex_maxOccupancy_dict[Tindex] - args.MAXIMUM_OCCUPANCY_TOLERANCE:
                new_i_iminus1_pool_dict[i].append(triplet)
    
    del i_iminus1_pool_dict
    i_iminus1_pool_dict = new_i_iminus1_pool_dict
    del new_i_iminus1_pool_dict
    
    outfile = args.POOL_CONNECTIVITIES_FILE+".mod"
    if args.POOL_CONNECTIVITIES_FILE == args.COMPLETE_CONNECTIVITIES_FILE:
        outfile = args.POOL_CONNECTIVITIES_FILE+".pool.mod"   # avoid ovewriting the complete aa types file
    with open(outfile, 'w') as f:
        f.write("i\tpossible i-1\n")
        for k,v in i_iminus1_pool_dict.items():
            f.write(k+"\t"+', '.join(str(e) for e in v)+"\n")
    
    with open(args.COMPLETE_CONNECTIVITIES_FILE+".mod", 'w') as f:
        f.write("i\tpossible i-1\n")
        for k,v in i_iminus1_complete_dict.items():
            f.write(k+"\t"+', '.join(str(e) for e in v)+"\n")

    sys.exit(0)
    
else:

    
    Nindices_with_connectivities_list = []
    def modify_i_iminus1_pool_dict(i, iminus1):
        global i_iminus1_pool_dict, i_iminus1_complete_dict, Nindices_with_connectivities_list
        
        Nindices_with_connectivities_list.append(iminus1)   # keep record of the added Nindices
        try:
            for triplet in i_iminus1_pool_dict[i]:
                if triplet[0] == iminus1:
                    print "MODIFYING POOLCONFILE: Setting connectivity of ", i, " to ", triplet
                    i_iminus1_pool_dict[i] = [triplet]
                    break
        except KeyError:    # CLEVER: add a new connection that is not present in the connectivities file
            if args.ADD_MISSING_CONNECTIVITIES:
                print "MODIFYING POOLCONFILE: Adding missing connectivity of ", i, " to ", (iminus1, 1, 1)
                i_iminus1_pool_dict[i] = [(iminus1, 1, 1)] # set occupancy and TOCSY resonance number arbitrarily to 1 and 1
        
        i_list = i_iminus1_pool_dict.keys()
        try:
            i_list.remove(i)
        except ValueError:
            print "WARNING: ",i," not in connectivities file! This connectivity file wasn't used to create the results_summary file you supplied!"
        for other_i in i_list:
            new_triplet_list = []
            for triplet in i_iminus1_pool_dict[other_i]:
                if triplet[0] != iminus1:
                    new_triplet_list.append(triplet)
            i_iminus1_pool_dict[other_i] = new_triplet_list
        
        found_aatype = False
        try:
            for triplet in i_iminus1_complete_dict[i]:
                if triplet[0] == iminus1:
                    found_aatype = True
                    break
        except KeyError:    # CLEVER: add a new connection that is not present in the connectivities file
            if args.ADD_MISSING_CONNECTIVITIES:
                print "MODIFYING CONFILE: Adding missing connectivity of ", i, " to ", (iminus1, 1, 1)
                i_iminus1_complete_dict[i] = [(iminus1, 1, 1)] # set occupancy and TOCSY resonance number arbitrarily to 1 and 1
        
        if found_aatype == False and args.ADD_MISSING_CONNECTIVITIES:
            print "MODIFYING CONFILE: Adding missing connectivity of ", i, " to ", (iminus1, 1, 1)
            i_iminus1_complete_dict[i] = [(iminus1, 1, 1)] # set occupancy and TOCSY resonance number arbitrarily to 1 and 1
        
    def add_back_Tindex_to_i_iminus1_pool_dict(i):
        """
            FUNCTION to read back the original connectivities for a particular Tindex(i), but also to place back the wrong Tindex(i-1) to all other positions that was in the
            original connectivities file.
        """
        global i_iminus1_pool_dict, i_iminus1_complete_dict
        
        wrong_iminus1 = i_iminus1_pool_dict[i][0][0]
        for Tindex, connectivities_list in i_iminus1_complete_dict.items():
            for triplet in connectivities_list:
                if wrong_iminus1 in triplet and triplet not in i_iminus1_pool_dict[Tindex]:
                    i_iminus1_pool_dict[Tindex].append(triplet) # add the wrong Tindex(i-1) to every position that existed before filtering
        i_iminus1_pool_dict[i] = i_iminus1_complete_dict[i]
        iaaindex_iminus1aaTypesProbPoolTupleList_dict[i] = iaaindex_iminus1aaTypesProbTupleList_dict[i]
    
    
    def remove_iminus1_from_other_positions(i, iminus1):
        """
            FUNCTION to remove iminus1 from connectivities others that i's
        """
        global i_iminus1_pool_dict
        
        i_list = i_iminus1_pool_dict.keys()
        try:
            i_list.remove(i)
        except ValueError:
            print "WARNING: ",i," not in connectivities file! This connectivity file wasn't used to create the results_summary file you supplied!"
        for other_i in i_list:
            new_triplet_list = []
            for triplet in i_iminus1_pool_dict[other_i]:
                if triplet[0] != iminus1:
                    new_triplet_list.append(triplet)
            i_iminus1_pool_dict[other_i] = new_triplet_list
    
    
    def modify_iaaindex_iminus1aaTypesProbPoolTupleList_dict(Tindex, aatype):
        global iaaindex_iminus1aaTypesProbPoolTupleList_dict, iaaindex_iminus1aaTypesProbTupleList_dict
        
        found_aatype = False
        try:
            for duplet in iaaindex_iminus1aaTypesProbPoolTupleList_dict[Tindex]:
                if duplet[0] == aatype:
                    print "MODIFYING POOLAAFILE: Setting aa type prediction of ", Tindex," to ", [duplet]
                    iaaindex_iminus1aaTypesProbPoolTupleList_dict[Tindex] = [duplet]
                    found_aatype = True
                    break
        except KeyError:    # CLEVER: add a new aa type that is not present in the predictions file
            print "MODIFYING POOLAAFILE: Adding missing aa type prediction of ", Tindex, " to ", (aatype, 1.0)
            iaaindex_iminus1aaTypesProbPoolTupleList_dict[Tindex] = [(aatype, 1.0)]
            
        if found_aatype == False:
            print "MODIFYING POOLAAFILE: Adding missing aa type prediction of ", Tindex, " to ", (aatype, 1.0)
            iaaindex_iminus1aaTypesProbPoolTupleList_dict[Tindex] = [(aatype, 1.0)]    # se the probability arbitrarily to 10
        
        found_aatype = False
        try:
            for duplet in iaaindex_iminus1aaTypesProbTupleList_dict[Tindex]:
                if duplet[0] == aatype:
                    found_aatype = True
                    break
        except KeyError:    # CLEVER: add a new aa type prediction that is not present in the aa type predictions file
            print "MODIFYING AAFILE: Adding missing aa type prediction of ", Tindex, " to ", (aatype, 1.0)
            iaaindex_iminus1aaTypesProbTupleList_dict[Tindex] = [(aatype, 1.0)]
            
        if found_aatype == False:
            print "MODIFYING AAFILE: Adding missing aa type prediction of ", Tindex, " to ", (aatype, 1.0)
            iaaindex_iminus1aaTypesProbTupleList_dict[Tindex] = [(aatype, 1.0)]    # se the Probability arbitrarily to 10
    
    
    absolute_match_chains_list = []
    chain = []
    for position in range(len(absolute_matches_alignment_list)):
        Tindex = absolute_matches_alignment_list[position]
        if Tindex != '-' and Tindex != 'N/A':
            chain.append(Tindex)
        else:
            if len(chain) > 1:  # don't save single Tocsy indices because they give no connectivity information
                absolute_match_chains_list.append(chain)
            chain = []
    
    
    for chain in absolute_match_chains_list:
        N = len(chain)-1    # -1 because the indexing in lists begins from 0
        index = N-1
        while index >=0:
            i = chain[index+1]
            iminus1 = chain[index]
            modify_i_iminus1_pool_dict(i, iminus1)   # keep only iminus1 TOCSY index as value
            index -= 1
    if absolute_matches_alignment_list[1] != '-':
        i_iminus1_pool_dict[absolute_matches_alignment_list[1]] = []
    for chain in absolute_match_chains_list:
        N = len(chain)-1    # -1 because the indexing in lists begins from 0
        index = N-1
        while index >=0:
            i = chain[index+1]
            iminus1 = chain[index]
            remove_iminus1_from_other_positions(i, iminus1)
            index -= 1
    
    
    
    
    for position in range(1, len(absolute_matches_alignment_list[1:])): # ommit first element: 'N/A'
        Tindex = absolute_matches_alignment_list[position]
        if Tindex != '-':
            chain.append(Tindex)
            aatype = protein_alignment_list[position-1]
            modify_iaaindex_iminus1aaTypesProbPoolTupleList_dict(Tindex, aa1to3_dict[aatype])
    
    for i in i_iminus1_pool_dict.keys():
        if len(i_iminus1_pool_dict[i]) == 1 and len(i_iminus1_complete_dict[i]) > 1 and ( not i in absolute_matches_alignment_list or absolute_matches_alignment_list[absolute_matches_alignment_list.index(i)-1] in ['-', 'N/A']):
            self_correct = False
            for Nindex_triplet in i_iminus1_complete_dict[i]:
                if Nindex_triplet != i_iminus1_pool_dict[i][0] and Nindex_triplet[0] not in Nindices_with_connectivities_list:
                    self_correct = True
                    break
            if self_correct == True:
                print "SELF-CORRECTION: found possible mistake in connectivities. Reading to the pools all possible connectivities and aa type predictions of ",i
                add_back_Tindex_to_i_iminus1_pool_dict(i)
        elif len(iaaindex_iminus1aaTypesProbPoolTupleList_dict[i]) == 1 and len(iaaindex_iminus1aaTypesProbTupleList_dict[i]) > 1 and not i in absolute_matches_alignment_list:
            print "SELF-CORRECTION: found possible mistake in aa type predictions. Reading to the pools all possible connectivities and aa type predictions of ",i
            i_iminus1_pool_dict[i] = i_iminus1_complete_dict[i]
            iaaindex_iminus1aaTypesProbPoolTupleList_dict[i] = iaaindex_iminus1aaTypesProbTupleList_dict[i]
    
    for Tindex in iaaindex_iminus1aaTypesProbPoolTupleList_dict.keys():
        if args.KEEP_ONLY_GLY:
            duplet_list = iaaindex_iminus1aaTypesProbPoolTupleList_dict[Tindex]
            if len(duplet_list) > 1 and duplet_list[0][0] == 'GLY' and duplet_list[0][1] >= 100 * duplet_list[1][1]:
                print "Setting amino acid type of TOCSY index group ", Tindex, " to ", [duplet_list[0]]
                iaaindex_iminus1aaTypesProbPoolTupleList_dict[Tindex] = [duplet_list[0]]
        if args.KEEP_ONLY_ALA:
            duplet_list = iaaindex_iminus1aaTypesProbPoolTupleList_dict[Tindex]
            if len(duplet_list) > 1 and duplet_list[0][0] == 'ALA' and duplet_list[0][1] >= 100 * duplet_list[1][1]:
                print "Setting amino acid type of TOCSY index group ", Tindex, " to ", [duplet_list[0]]
                iaaindex_iminus1aaTypesProbPoolTupleList_dict[Tindex] = [duplet_list[0]]
                
    
    

outfile = args.POOL_CONNECTIVITIES_FILE+".mod"
if args.POOL_CONNECTIVITIES_FILE == args.COMPLETE_CONNECTIVITIES_FILE:
    outfile = args.POOL_CONNECTIVITIES_FILE+".pool.mod"   # avoid ovewriting the complete aa types file
with open(outfile, 'w') as f:
    f.write("i\tpossible i-1\n")
    for k,v in i_iminus1_pool_dict.items():
        f.write(k+"\t"+', '.join(str(e) for e in v)+"\n")

with open(args.COMPLETE_CONNECTIVITIES_FILE+".mod", 'w') as f:
    f.write("i\tpossible i-1\n")
    for k,v in i_iminus1_complete_dict.items():
        f.write(k+"\t"+', '.join(str(e) for e in v)+"\n")

with open(args.COMPLETE_AA_TYPES_FILE+".mod", 'w') as f:
    f.write("i aa index\tpossible i-1 aa types\n")
    for i_aaindex in iaaindex_iminus1aaTypesProbTupleList_dict.keys():
        f.write(i_aaindex + " ---> (i-1) aa type " + str(sorted(iaaindex_iminus1aaTypesProbTupleList_dict[i_aaindex], key=itemgetter(1), reverse=True)) + "\n")

outfile = args.POOL_AA_TYPES_FILE+".mod"
if args.POOL_AA_TYPES_FILE == args.COMPLETE_AA_TYPES_FILE:
    outfile = args.POOL_AA_TYPES_FILE+".pool.mod"   # avoid ovewriting the complete aa types file
with open(outfile, 'w') as f:
    f.write("i aa index\tpossible i-1 aa types\n")
    for i_aaindex in iaaindex_iminus1aaTypesProbPoolTupleList_dict.keys():
        f.write(i_aaindex + " ---> (i-1) aa type " + str(sorted(iaaindex_iminus1aaTypesProbPoolTupleList_dict[i_aaindex], key=itemgetter(1), reverse=True)) + "\n")