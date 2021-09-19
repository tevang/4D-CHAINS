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
import gc
import sys
from argparse import ArgumentParser
from operator import itemgetter

import numpy as np

## GLOBAL VARIABLE DEFINITION
from ete3 import Tree
from lib.alignment import Alignment
from lib.bayes.statistics import Probability
from lib.connectivities import Connectivities
from lib.global_func import Debuginfo, save_pickle
from scipy.stats import zscore

code3 = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
code1 = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
aa1to3_dict = dict((c1,c3) for c1,c3 in zip(code1,code3))
aa3to1_dict = dict((c3,c1) for c3,c1 in zip(code3,code1))

def approx_equal(x, y, tolerance=0.001):
    return abs(x-y) <= 0.5 * tolerance * (x + y)

## Parse command line arguments
def cmdlineparse():
    parser = ArgumentParser(description="command line arguments",
        epilog="EXAMPLES: \
               If you run the script for the 1st time use the complete connectivities and amino acid type prediction files as pool files:\n \
        ../modify_connectivities_and_aatypes.py -rstfile consensus_alignment.overlapped_chains_common_sequence.iteration7.6mers -poolconfile connectivities_all -allconfile connectivities_all -poolaafile amino_acid_type_prediction_probabilities -allaafile amino_acid_type_prediction_probabilities \
        If you have alread created pool files from previous runs do:\n \
        ../modify_connectivities_and_aatypes.py -rstfile consensus_alignment.overlapped_chains_common_sequence.iteration7.6mers -poolconfile connectivities_all.pool -allconfile connectivities_all -poolaafile amino_acid_type_prediction_probabilities.pool -allaafile amino_acid_type_prediction_probabilities")
    parser.add_argument("-rstfile", dest="ABSOLUTE_MATCHES_FILE", required=True,
                        help="file with absolute matches from previous run to be used as restraints. (default: %(default)s)",
                        metavar="<absolute matches restraint file>")
    
    parser.add_argument("-poolconfile", dest="POOL_CONNECTIVITIES_FILE", required=True,
                        help="pool connectivities file; if specified connectivity calculation from input files will be skipped. (default: %(default)s)",
                        metavar="<pool connectivities file>")
    
    parser.add_argument("-allconfile", dest="COMPLETE_CONNECTIVITIES_FILE", required=True,
                        help="all connectivities file; necessary if -confile specified in order to calculate correct probabilities. (default: %(default)s)",
                        metavar="<all connectivities file>")
    
    parser.add_argument("-poolaafile", dest="POOL_AA_TYPES_FILE", required=True,
                        help="pool amino acid assignment file; if specified amino acid assignment calculation from input files will be skipped. (default: %(default)s)",
                        metavar="<pool amino acid assignment file>")
    
    parser.add_argument("-allaafile", dest="COMPLETE_AA_TYPES_FILE", required=True,
                        help="all amino acid assignment file; necessary if -aafile specified in order to calculate correct probabilities. (default: %(default)s)",
                        metavar="<all amino acid assignment file>")
    
    parser.add_argument("-keepgly", dest="KEEP_ONLY_GLY", required=False, default=False, action='store_true',
                        help="if GLY is the top ranked prediction with a probability at least 2 orders of magnitude greater than the 2nd \
                             prediction, keep only GLY. (default: %(default)s)")
    parser.add_argument("-keepala", dest="KEEP_ONLY_ALA", required=False, default=False, action='store_true',
                        help="if ALA is the top ranked prediction with a probability at least 2 orders of magnitude greater than the 2nd prediction, \
                             keep only ALA. (default: %(default)s)")
    
    parser.add_argument("-nocorrect", dest="SELF_CORRECTION", required=False, action='store_false', default=True,
                        help="do not do self-correction. (default: %(default)s)")
    
    parser.add_argument("-addconn", dest="ADD_MISSING_CONNECTIVITIES", required=False, action='store_true', default=False,
                        help="add missing connectivities. (default: %(default)s)")
    
    parser.add_argument("-bigchains", dest="BIG_CHAIN_CUTOFF", required=False, type=int, default=None,
                        help="keep only connectivities that form chains of length equal or greater of the specified value. This options is useful \
                             when you have a big protein and you want to reduce the number of possible peptides. (default: %(default)s)",
                        metavar="<big chain cutoff length>")
    
    parser.add_argument("-mcutoff", dest="RESONANCE_MATCH_CUTOFF", required=False, type=float, default=0.8,
                        help="number 0.0-1.0 saying how much of TOCSY resonances should match in the NOESY in order to consider it a possible match. \
                             (default: %(default)s)",
                        metavar="<resonance match cutoff>")
    
    parser.add_argument("-zmcutoff", dest="ZSCORE_RESONANCE_MATCH_CUTOFF", required=False, type=float, default=0.0,
                        help="a real number specifying the lower Z-score for a resonance match to be retained for chain building. (default: %(default)s)",
    metavar="<Z-score resonance match cutoff>")
    
    parser.add_argument("-maxocc", dest="MAXIMUM_OCCUPANCY_TOLERANCE", required=False, type=int, default=None,
                       help="keep only the connectivities of a TOCSY index group with the highest occupancy or those that differ from the maximum \
                            by -maxocc number (recommended: 0). (default: %(default)s)",
                       metavar="<maximum occupancy tolerance>")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s version 0.8')
    args=parser.parse_args()
    return args

args = cmdlineparse()


## READ FILE WITH ABSOLUTE MATCHES
protein_alignment_list = []
absolute_AAIGmatches_alignment_list = []
if args.BIG_CHAIN_CUTOFF == None and args.MAXIMUM_OCCUPANCY_TOLERANCE == None:
    absolute_AAIGmatches_alignment_list, protein_alignment_list = \
        Alignment.read_NHmap_file(args.ABSOLUTE_MATCHES_FILE, get_protein_alignment=True)
    # print("DEBUG: protein_alignment_list=", protein_alignment_list)
    # print("DEBUG: absolute_AAIGmatches_alignment_list=", absolute_AAIGmatches_alignment_list)


print("Loading Pool Connectivities from file", args.POOL_CONNECTIVITIES_FILE)
# a dictionary containing all possible connectivities of every TOCSY AAIG; has keys the TOCSY aa indices and values lists of triplets (tuples) \
# consisting of (NOESYAAIG, occupancy, numOfResonances)
i_iminus1_pool_dict = Connectivities.load_connectivities_from_file(args.POOL_CONNECTIVITIES_FILE)

print("Loading All Connectivities from file", args.COMPLETE_CONNECTIVITIES_FILE)
# Now load the complete connectivities file to calculate correct probabilities
# the same dictionary but with all possible connectivities, including those below args.RESONANCE_MATCH_CUTOFF
i_iminus1_complete_dict = Connectivities.load_connectivities_from_file(args.COMPLETE_CONNECTIVITIES_FILE)

if args.MAXIMUM_OCCUPANCY_TOLERANCE == None:
    print("Loading Cutoff Amino Acid Prediction file ...")
    # OLD CODE
    # iAAIG_iminus1aaTypesProbPoolTupleList_dict = OrderedDict() # ordereddict with keys the AAIG of residue i and values lists of tuples of the form
    #                                                       # (residue i-1 matching aa type, Z-score)
    # with open(args.POOL_AA_TYPES_FILE, 'r') as f:
    #     aa_type_file_contents = f.readlines()
    #     for line in aa_type_file_contents[1:]:
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
    iAAIG_iminus1aaTypesProbPoolTupleList_dict = Probability.load_aatype_probs_from_file(
        fname=args.POOL_AA_TYPES_FILE)  # ordereddict with keys the AAIG of residue i and values lists of tuples of the form
                                                          # (residue i-1 matching aa type, Z-score)
    
    print("Loading All Amino acid Prediction file ...")
    # OLD CODE
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
    iAAIG_iminus1aaTypesProbTupleList_dict = Probability.load_aatype_probs_from_file(
        fname=args.COMPLETE_AA_TYPES_FILE)  # ordereddict with keys the AAIG of residue i and values lists of tuples of the form
                                        # (residue i-1 matching aa type, average probability)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##              FORM CHAINS FROM CONNECTIVITIES               ##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

if args.BIG_CHAIN_CUTOFF:   # KEEP ONLY CONNECTIVITIES THAT FORM CHAINS OF LENGTH EQUAL OR GREATER OF THE SPECIFIED VALUE.

    
    # FIRST REMOVE THE CONNECTIVITIES BELOW THE CUTOFF
    i2remove_set = set()
    for i in list(i_iminus1_pool_dict.keys()):
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
    
    # FILTER i_iminus1_pool_dict ACCORDING TO THE Z-SCORE CUTOFF
    new_i_iminus1_pool_dict = {}
    for i in list(i_iminus1_pool_dict.keys()):
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
            ## ONLY IF THE Z-SCORE OF THE CONNECTIVITY IS GREATER THAN THE CUTOFF
            if Zscore >= args.ZSCORE_RESONANCE_MATCH_CUTOFF or approx_equal(Zscore, args.ZSCORE_RESONANCE_MATCH_CUTOFF):
                new_i_iminus1_pool_dict[i].append(triplet)
    
    del i_iminus1_pool_dict
    i_iminus1_pool_dict = new_i_iminus1_pool_dict
    del new_i_iminus1_pool_dict
    
    
    # NOW DEFINE THE NECESSARY FUNCTIONS TO BUILD THE CHAIN TREES
    #def populate_leaves():
    #    """
    #        FUNCTION that adds new branches to the leaves of the Tree.
    #        ARGUMENTS:
    #        Assignment_Tree:    The Tree structure with connectivities
    #        RETURNS:
    #        (Assignment_Tree, BOOLEAN):    A tuple with elements the input Tree structure with new branches (if applicable), and a BOOLEAN value which is True if the function added
    #                                       new leaves to the Tree, or False otherwise
    #    """
    #    
    #    global i_iminus1_pool_dict, Assignment_Tree
    #    number_of_new_leaves = 0
    #    # ATTENTION: never use Assignment_Tree.iter_leaf_names(), it doesn't return the names in the order
    #    # ATTENTION: corresponding to Assignment_Tree.get_leaves()!!!
    #    for leaf in Assignment_Tree.get_leaves():
    #        try:
    #            #for child_tuple in i_iminus1_pool_dict[name]:
    #            for child_triplet in i_iminus1_pool_dict[name]: # WARNING: i_iminus1_normProbabilities_dict contains only matches above the Z-score cutoff
    #                NOESYAAIG = child_triplet[0]
    #                _occupancy = child_triplet[1]     # add prefix "_" to discriminate from new child feature "occupancy" that will be added
    #                _numOfResonances = child_triplet[2]
    #                ancestors_list = [ancestor.name for ancestor in leaf.get_ancestors()]
    #                if NOESYAAIG in ancestors_list:     # if the current NOESY AAIG is already a node or leaf in the Tree, continue to the next
    #                    continue
    #                #new_child = leaf.add_child(name=NOESYAAIG, dist=float(_occupancy)/_numOfResonances) # add a new brach to the current TOCSY add index (leaf) with length the ratio occupancy/numOfResonances
    #                
    #                new_child = leaf.add_child(name=NOESYAAIG) # add a new brach to the current TOCSY add index (leaf) with length the respective probability
    #                
    #                new_child.add_features(occupancy=_occupancy, numOfResonances=_numOfResonances)
    #                number_of_new_leaves += 1
    #                #print "DEBUG: adding connection: ",name,"-->",NOESYAAIG
    #        except KeyError:
    #            continue
    #    
    #    del leaf, name, new_child, ancestor, ancestors_list, child_triplet, NOESYAAIG, _occupancy, _numOfResonances, number_of_new_leaves
    #    #print Assignment_Tree.get_ascii(show_internal=True, compact=False, attributes=["name", "dist", "occupancy", "numOfResonances"])
    #    #print Assignment_Tree.get_ascii(show_internal=True, compact=False)
    #    if number_of_new_leaves > 0:
    #        return True
    #    else:
    #        return False
            
    
    Tindices2keep_set = set()    # a list of lists containing the possible connectivities and the overall score as the last element
    
    def build_Tree(i, Tindices2keep_set):

        def populate_leaves(Assignment_Tree, i_iminus1_pool_dict):
            """
                FUNCTION that adds new branches to the leaves of the Tree.
                ARGUMENTS:
                Assignment_Tree:    The Tree structure with connectivities
                RETURNS:
                (Assignment_Tree, BOOLEAN):    A tuple with elements the input Tree structure with new branches (if applicable), and a BOOLEAN value which is True if the function added
                                               new leaves to the Tree, or False otherwise
            """
            number_of_new_leaves = 0
            # ATTENTION: never use Assignment_Tree.iter_leaf_names(), it doesn't return the names in the order
            # ATTENTION: corresponding to Assignment_Tree.get_leaves()!!!
            for leaf in Assignment_Tree.get_leaves():
                try:
                    #for child_tuple in i_iminus1_pool_dict[name]:
                    for child_triplet in i_iminus1_pool_dict[leaf.name]: # WARNING: i_iminus1_normProbabilities_dict contains only matches above the Z-score cutoff
                        NOESYAAIG = child_triplet[0]
                        _occupancy = child_triplet[1]     # add prefix "_" to discriminate from new child feature "occupancy" that will be added
                        _numOfResonances = child_triplet[2]
                        ancestors_list = [ancestor.name for ancestor in leaf.get_ancestors()]
                        if NOESYAAIG in ancestors_list:     # if the current NOESY AAIG is already a node or leaf in the Tree, continue to the next
                            continue
                        #new_child = leaf.add_child(name=NOESYAAIG, dist=float(_occupancy)/_numOfResonances) # add a new brach to the current TOCSY add index (leaf) with length the ratio occupancy/numOfResonances
                        
                        new_child = leaf.add_child(name=NOESYAAIG) # add a new brach to the current TOCSY add index (leaf) with length the respective probability
                        
                        new_child.add_features(occupancy=_occupancy, numOfResonances=_numOfResonances)
                        number_of_new_leaves += 1
                        #print "DEBUG: adding connection: ",name,"-->",NOESYAAIG
                except KeyError:
                    continue
            
            del leaf
            #print Assignment_Tree.get_ascii(show_internal=True, compact=False, attributes=["name", "dist", "occupancy", "numOfResonances"])
            #print Assignment_Tree.get_ascii(show_internal=True, compact=False)
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
            expand_tree = populate_leaves(Assignment_Tree, i_iminus1_pool_dict)
            level += 1
            if level == args.BIG_CHAIN_CUTOFF:
                break
        if level < args.BIG_CHAIN_CUTOFF: # discard chains shorter than the cutoff
            return
        # Print the Tree
        #print Assignment_Tree.get_ascii(show_internal=True, compact=False)
        #print Assignment_Tree.get_ascii(show_internal=True, compact=False, attributes=["name", "dist", "occupancy", "numOfResonances"])
        
        print("\nSaving chains from Tree...")
    
        for leaf in Assignment_Tree.get_leaves():
            chain = []
            score = leaf.dist
            chain.append(leaf.name)
            for ancestor in leaf.get_ancestors():
                if ancestor == Root:
                    break
                chain.append(ancestor.name)
            if len(chain) >= args.BIG_CHAIN_CUTOFF:
                for TAAIG in chain:
                    Tindices2keep_set.add(TAAIG)
            del chain
            del ancestor
            del leaf
        del Assignment_Tree
        gc.collect()
    
    
    N = len(list(i_iminus1_pool_dict.keys()))
    index = 1
    for i in list(i_iminus1_pool_dict.keys()):
        print("Building Tree starting from amino acid index",i,"("+str(index)+"/"+str(N)+") ...")
        build_Tree(i, Tindices2keep_set)
        index += 1
    
    
    # CLEAN DICTIONARIES FOR UNUSED TOCSY INDEX GROUPS
    if len(Tindices2keep_set) == len(list(i_iminus1_pool_dict.keys())):
        print("NO TOCSY INDEX GROUPS WERE REMOVED FROM CONNECTIVITY FILE !!!")
        sys.exit(0)
    for TAAIG in list(i_iminus1_pool_dict.keys()):
        if TAAIG not in Tindices2keep_set:
            print("Deleting TOCSY index group", TAAIG, " from connectivities pool file.")
            del i_iminus1_pool_dict[TAAIG]
            del i_iminus1_complete_dict[TAAIG]
            del iAAIG_iminus1aaTypesProbTupleList_dict[TAAIG]
            del iAAIG_iminus1aaTypesProbPoolTupleList_dict[TAAIG]
    
    ## KEEP HIGH FINDELITY GLY AND ALA PREDICTIONS
    for TAAIG in list(iAAIG_iminus1aaTypesProbPoolTupleList_dict.keys()):
        if args.KEEP_ONLY_GLY:
            duplet_list = iAAIG_iminus1aaTypesProbPoolTupleList_dict[TAAIG]
            # if GLY is the top ranked prediction with a probability at least 2 orders of magnitude greater than the 2nd prediction, keep only GLY
            if len(duplet_list) > 1 and duplet_list[0][0] == 'GLY' and duplet_list[0][1] >= 100 * duplet_list[1][1]:
                print("Setting amino acid type of TOCSY index group ", TAAIG, " to ", [duplet_list[0]])
                iAAIG_iminus1aaTypesProbPoolTupleList_dict[TAAIG] = [duplet_list[0]]
        if args.KEEP_ONLY_ALA:
            duplet_list = iAAIG_iminus1aaTypesProbPoolTupleList_dict[TAAIG]
            # if ALA is the top ranked prediction with a probability at least 2 orders of magnitude greater than the 2nd prediction, keep only ALA
            if len(duplet_list) > 1 and duplet_list[0][0] == 'ALA' and duplet_list[0][1] >= 100 * duplet_list[1][1]:
                print("Setting amino acid type of TOCSY index group ", TAAIG, " to ", [duplet_list[0]])
                iAAIG_iminus1aaTypesProbPoolTupleList_dict[TAAIG] = [duplet_list[0]]
    
elif args.MAXIMUM_OCCUPANCY_TOLERANCE != None:  # KEEP ALL CONNECTIVITIES (not only those forming long chains)
    # AND KEEP ONLY THE CONNECTIVITIES WITH THE MAXIMUM OCCUPANCY FOR EACH TOCSY INDEX GROUP

    # Find the maximum occupancy for every TOCSY index group
    TAAIG_maxOccupancy_dict = {}
    for triplet_list in list(i_iminus1_pool_dict.values()):
        for triplet in triplet_list:
            TAAIG = triplet[0]
            occupancy = int(triplet[1])
            if TAAIG not in list(TAAIG_maxOccupancy_dict.keys()):
                TAAIG_maxOccupancy_dict[TAAIG] = occupancy
            elif TAAIG_maxOccupancy_dict[TAAIG] < occupancy:
                TAAIG_maxOccupancy_dict[TAAIG] = occupancy
    
    # FILTER i_iminus1_pool_dict ACCORDING TO THE MAXIMUM OCCUPANCY TOLERANCE
    new_i_iminus1_pool_dict = {}
    for i in list(i_iminus1_pool_dict.keys()):
        new_i_iminus1_pool_dict[i] = []  # add i to the new dictionary
        for triplet in i_iminus1_pool_dict[i]:
            TAAIG = triplet[0]
            occupancy = triplet[1]
            if occupancy >= TAAIG_maxOccupancy_dict[TAAIG] - args.MAXIMUM_OCCUPANCY_TOLERANCE:
                new_i_iminus1_pool_dict[i].append(triplet)
    
    del i_iminus1_pool_dict
    i_iminus1_pool_dict = new_i_iminus1_pool_dict
    del new_i_iminus1_pool_dict
    
    # WRITE CONNECTIVITIES FILES AND EXIT
    # Write modifies Connectivities Files
    outfile = args.POOL_CONNECTIVITIES_FILE+".mod"
    if args.POOL_CONNECTIVITIES_FILE == args.COMPLETE_CONNECTIVITIES_FILE:
        outfile = args.POOL_CONNECTIVITIES_FILE+".pool.mod"   # avoid ovewriting the complete aa types file
    with open(outfile, 'w') as f:
        f.write("i\tpossible i-1\n")
        for k,v in list(i_iminus1_pool_dict.items()):
            #print k,v
            f.write(k+"\t"+', '.join(str(e) for e in v)+"\n")
    
    with open(args.COMPLETE_CONNECTIVITIES_FILE+".mod", 'w') as f:
        f.write("i\tpossible i-1\n")
        for k,v in list(i_iminus1_complete_dict.items()):
            #print k,v
            f.write(k+"\t"+', '.join(str(e) for e in v)+"\n")

    sys.exit(0)
    
else:

    #########################################################################
    #           MODIFY CONNECTIVITIES AND AA TYPE PREDICTIONS               #
    #########################################################################
    
    Nindices_with_connectivities_list = []
    def modify_i_iminus1_pool_dict(i, iminus1, i_iminus1_pool_dict,
                                   i_iminus1_complete_dict, Nindices_with_connectivities_list):
        
        Nindices_with_connectivities_list.append(iminus1)   # keep record of the added Nindices
        ## MODIFY i_iminus1_pool_dict
        try:
            for triplet in i_iminus1_pool_dict[i]:
                if triplet[0] == iminus1:
                    print("MODIFYING POOLCONFILE: Setting connectivity of ", i, " to ", triplet)
                    i_iminus1_pool_dict[i] = [triplet]
                    break
        except KeyError:    # CLEVER: add a new connection that is not present in the connectivities file
            if args.ADD_MISSING_CONNECTIVITIES:
                print("MODIFYING POOLCONFILE: Adding missing connectivity of ", i, " to ", (iminus1, 1, 1))
                i_iminus1_pool_dict[i] = [(iminus1, 1, 1)] # set occupancy and TOCSY resonance number arbitrarily to 1 and 1
        
        # remove iminus1 from the list of values of every other key
        #if i=="I202": print "DEBUG: i_iminus1_pool_dict before =", i_iminus1_pool_dict
        i_list = list(i_iminus1_pool_dict.keys())
        try:
            i_list.remove(i)
        except ValueError:
            print("WARNING: ",i," not in connectivities file! This connectivity file wasn't used to create the results_summary file you supplied!")
        for other_i in i_list:
            new_triplet_list = []
            for triplet in i_iminus1_pool_dict[other_i]:
                if triplet[0] != iminus1:
                    new_triplet_list.append(triplet)
            i_iminus1_pool_dict[other_i] = new_triplet_list
        #if i=="I202": print "DEBUG: i_iminus1_pool_dict after =", i_iminus1_pool_dict
        
        # NOW MODIFY i_iminus1_complete_dict
        found_aatype = False
        try:
            for triplet in i_iminus1_complete_dict[i]:
                if triplet[0] == iminus1:
                    # don't remove all alternative predictions otherwise the probabilities that will be calculated will be wrong!
                    found_aatype = True
                    break
        except KeyError:    # CLEVER: add a new connection that is not present in the connectivities file
            if args.ADD_MISSING_CONNECTIVITIES:
                print("MODIFYING CONFILE: Adding missing connectivity of ", i, " to ", (iminus1, 1, 1))
                i_iminus1_complete_dict[i] = [(iminus1, 1, 1)] # set occupancy and TOCSY resonance number arbitrarily to 1 and 1
        
        # if the connectivity does not exist, add it!
        if found_aatype == False and args.ADD_MISSING_CONNECTIVITIES:
            print("MODIFYING CONFILE: Adding missing connectivity of ", i, " to ", (iminus1, 1, 1))
            i_iminus1_complete_dict[i] = [(iminus1, 1, 1)] # set occupancy and TOCSY resonance number arbitrarily to 1 and 1
        
    def add_back_TAAIG_to_i_iminus1_pool_dict(i, i_iminus1_pool_dict, i_iminus1_complete_dict):
        """
            FUNCTION to read back the original connectivities for a particular TAAIG(i), but also to place back the wrong TAAIG(i-1) to all other positions that was in the
            original connectivities file.
        """
        # PLACE BACK THE WRONG TAAIG(i-1)
        wrong_iminus1 = i_iminus1_pool_dict[i][0][0]
        for TAAIG, connectivities_list in list(i_iminus1_complete_dict.items()):
            for triplet in connectivities_list:
                if wrong_iminus1 in triplet and (TAAIG not in i_iminus1_pool_dict.keys() or
                                                 triplet not in i_iminus1_pool_dict[TAAIG]):
                    i_iminus1_pool_dict[TAAIG].append(triplet) # add the wrong TAAIG(i-1) to every position that existed before filtering
        # AT THE END READ BACK THE ORIGINAL CONNECTIVITIES FOR A PARTICULAR TAAIG(i)
        i_iminus1_pool_dict[i] = i_iminus1_complete_dict[i]
        iAAIG_iminus1aaTypesProbPoolTupleList_dict[i] = iAAIG_iminus1aaTypesProbTupleList_dict[i]
    
    
    def remove_iminus1_from_other_positions(i, iminus1, i_iminus1_pool_dict):
        """
            FUNCTION to remove iminus1 from connectivities others that i's
        """

        i_list = list(i_iminus1_pool_dict.keys())
        try:
            i_list.remove(i)
        except ValueError:
            print("WARNING: ",i," not in connectivities file! This connectivity file wasn't used to create the results_summary file you supplied!")
        for other_i in i_list:
            new_triplet_list = []
            for triplet in i_iminus1_pool_dict[other_i]:
                if triplet[0] != iminus1:
                    new_triplet_list.append(triplet)
            i_iminus1_pool_dict[other_i] = new_triplet_list
    
    
    def modify_iAAIG_iminus1aaTypesProbPoolTupleList_dict(TAAIG, aatype, iAAIG_iminus1aaTypesProbPoolTupleList_dict,
                                                          iAAIG_iminus1aaTypesProbTupleList_dict):
        
        found_aatype = False
        try:
            for duplet in iAAIG_iminus1aaTypesProbPoolTupleList_dict[TAAIG]:
                if duplet[0] == aatype:
                    print("MODIFYING POOLAAFILE: Setting aa type prediction of ", TAAIG," to ", [duplet])
                    iAAIG_iminus1aaTypesProbPoolTupleList_dict[TAAIG] = [duplet]
                    found_aatype = True
                    break
        except KeyError:    # CLEVER: add a new aa type that is not present in the predictions file
            print("MODIFYING POOLAAFILE: Adding missing aa type prediction of ", TAAIG, " to ", (aatype, 1.0))
            iAAIG_iminus1aaTypesProbPoolTupleList_dict[TAAIG] = [(aatype, 1.0)]
            
        # if the aatype prediction does not exist, add it!
        if found_aatype == False:
            print("MODIFYING POOLAAFILE: Adding missing aa type prediction of ", TAAIG, " to ", (aatype, 1.0))
            iAAIG_iminus1aaTypesProbPoolTupleList_dict[TAAIG] = [(aatype, 1.0)]    # se the probability arbitrarily to 10
        
        # Now modify iAAIG_iminus1aaTypesProbTupleList_dict
        found_aatype = False
        try:
            for duplet in iAAIG_iminus1aaTypesProbTupleList_dict[TAAIG]:
                if duplet[0] == aatype:
                    # don't remove all alternative predictions otherwise the probabilities that will be calculated will be wrong!
                    found_aatype = True
                    break
        except KeyError:    # CLEVER: add a new aa type prediction that is not present in the aa type predictions file
            print("MODIFYING AAFILE: Adding missing aa type prediction of ", TAAIG, " to ", (aatype, 1.0))
            iAAIG_iminus1aaTypesProbTupleList_dict[TAAIG] = [(aatype, 1.0)]
            
        # if the aatype prediction does not exist, add it!
        if found_aatype == False:
            print("MODIFYING AAFILE: Adding missing aa type prediction of ", TAAIG, " to ", (aatype, 1.0))
            iAAIG_iminus1aaTypesProbTupleList_dict[TAAIG] = [(aatype, 1.0)]    # se the Probability arbitrarily to 10
    
    
    absolute_match_chains_list = []
    chain = []
    for position in range(len(absolute_AAIGmatches_alignment_list)):
        TAAIG = absolute_AAIGmatches_alignment_list[position]
        if TAAIG != '-' and TAAIG != 'N/A':
            chain.append(TAAIG)
        else:
            if len(chain) > 1:  # don't save single Tocsy indices because they give no connectivity information
                absolute_match_chains_list.append(chain)
            chain = []
    
    #print "DEBUG: i_iminus1_complete_dict before =",
    #for k,v in i_iminus1_complete_dict.items():
    #    print k,v
    #print "DEBUG: before i_iminus1_pool_dict['R435']=", i_iminus1_pool_dict['R435']
    
    for chain in absolute_match_chains_list:
        N = len(chain)-1    # -1 because the indexing in lists begins from 0
        index = N-1
        while index >=0:
            i = chain[index+1]
            iminus1 = chain[index]
            modify_i_iminus1_pool_dict(i, iminus1,
                                       i_iminus1_pool_dict,
                                       i_iminus1_complete_dict,
                                       Nindices_with_connectivities_list)   # keep only iminus1 TOCSY index as value
            index -= 1
    ## FINALY MODIFY THE CONNECTIVITY if resid 2->resid 1, and clean resid 1 from every other position in the pool
    if absolute_AAIGmatches_alignment_list[1] != '-':
        i_iminus1_pool_dict[absolute_AAIGmatches_alignment_list[1]] = []
    # FINALY DO ONE LAST ROUND OF CLEANING the TAAIG(i-1) from every other position (I did one test and didn't see any difference)
    for chain in absolute_match_chains_list:
        N = len(chain)-1    # -1 because the indexing in lists begins from 0
        index = N-1
        while index >=0:
            i = chain[index+1]
            iminus1 = chain[index]
            remove_iminus1_from_other_positions(i, iminus1, i_iminus1_pool_dict)
            index -= 1
    
    #print "DEBUG: after i_iminus1_pool_dict['R435']=", i_iminus1_pool_dict['R435']
    #print "DEBUG: i_iminus1_pool_dict after =",
    #for k,v in i_iminus1_pool_dict.items():
    #    print k,v
    
    
    #print "DEBUG: iAAIG_iminus1aaTypesProbTupleList_dict before =",
    #for k,v in iAAIG_iminus1aaTypesProbTupleList_dict.items():
    #    print k,v
    
    ## MODIFY iAAIG_iminus1aaTypesProbPoolTupleList_dict
    for position in range(1, len(absolute_AAIGmatches_alignment_list[1:])): # ommit first element: 'N/A'
        TAAIG = absolute_AAIGmatches_alignment_list[position]
        if TAAIG != '-':
            chain.append(TAAIG)
            aatype = protein_alignment_list[position-1]
            modify_iAAIG_iminus1aaTypesProbPoolTupleList_dict(TAAIG, aa1to3_dict[aatype],
                                                              iAAIG_iminus1aaTypesProbPoolTupleList_dict,
                                                              iAAIG_iminus1aaTypesProbTupleList_dict)
    
    # SELF-CORRECT PREVIOUS WRONG CONSENSUS ABSOLUTE MATCHES
    #print "DEBUG: i_iminus1_pool_dict=", i_iminus1_pool_dict
    #print "DEBUG: i_iminus1_complete_dict=", i_iminus1_complete_dict
    #print "DEBUG: absolute_AAIGmatches_alignment_list=", absolute_AAIGmatches_alignment_list
    for i in list(i_iminus1_pool_dict.keys()):
        assert i in iAAIG_iminus1aaTypesProbTupleList_dict.keys() and i in iAAIG_iminus1aaTypesProbPoolTupleList_dict.keys(), \
            Debuginfo("FAIL: NOESY AAIG %s has connectivities but not aa-type predictions!" % i, fail=True)
        #print "DEBUG: i=", i, "i_iminus1_pool_dict[i]=", i_iminus1_pool_dict[i], "i_iminus1_complete_dict[i]=", i_iminus1_complete_dict[i]
        ## If the connectivities of i in the pool are different than in the original file, and i is not in the the alignments anymore or it follows a gap
        if len(i_iminus1_pool_dict[i]) == 1 and len(i_iminus1_complete_dict[i]) > 1 and ( not i in absolute_AAIGmatches_alignment_list or absolute_AAIGmatches_alignment_list[absolute_AAIGmatches_alignment_list.index(i)-1] in ['-', 'N/A']):
            self_correct = False
            for Nindex_triplet in i_iminus1_complete_dict[i]:
                if Nindex_triplet != i_iminus1_pool_dict[i][0] and Nindex_triplet[0] not in Nindices_with_connectivities_list:
                    self_correct = True
                    break
            if self_correct == True:
                print("SELF-CORRECTION: found possible mistake in connectivities. Reading to the pools all possible connectivities and aa type predictions of ",i)
                save_pickle('modify_conn.pkl', i, i_iminus1_pool_dict, i_iminus1_complete_dict)
                add_back_TAAIG_to_i_iminus1_pool_dict(i, i_iminus1_pool_dict, i_iminus1_complete_dict)
        elif len(iAAIG_iminus1aaTypesProbPoolTupleList_dict[i]) == 1 and len(iAAIG_iminus1aaTypesProbTupleList_dict[i]) > 1 and not i in absolute_AAIGmatches_alignment_list:
            print("SELF-CORRECTION: found possible mistake in aa type predictions. Reading to the pools all possible connectivities and aa type predictions of ",i)
            i_iminus1_pool_dict[i] = i_iminus1_complete_dict[i]
            iAAIG_iminus1aaTypesProbPoolTupleList_dict[i] = iAAIG_iminus1aaTypesProbTupleList_dict[i]
    
    ## KEEP HIGH FINDELITY GLY AND ALA PREDICTIONS
    for TAAIG in list(iAAIG_iminus1aaTypesProbPoolTupleList_dict.keys()):
        if args.KEEP_ONLY_GLY:
            duplet_list = iAAIG_iminus1aaTypesProbPoolTupleList_dict[TAAIG]
            # if GLY is the top ranked prediction with a probability at least 2 orders of magnitude greater than the 2nd prediction, keep only GLY
            if len(duplet_list) > 1 and duplet_list[0][0] == 'GLY' and duplet_list[0][1] >= 100 * duplet_list[1][1]:
                print("Setting amino acid type of TOCSY index group ", TAAIG, " to ", [duplet_list[0]])
                iAAIG_iminus1aaTypesProbPoolTupleList_dict[TAAIG] = [duplet_list[0]]
        if args.KEEP_ONLY_ALA:
            duplet_list = iAAIG_iminus1aaTypesProbPoolTupleList_dict[TAAIG]
            # if ALA is the top ranked prediction with a probability at least 2 orders of magnitude greater than the 2nd prediction, keep only ALA
            if len(duplet_list) > 1 and duplet_list[0][0] == 'ALA' and duplet_list[0][1] >= 100 * duplet_list[1][1]:
                print("Setting amino acid type of TOCSY index group ", TAAIG, " to ", [duplet_list[0]])
                iAAIG_iminus1aaTypesProbPoolTupleList_dict[TAAIG] = [duplet_list[0]]
                
    
    #print "DEBUG: iAAIG_iminus1aaTypesProbTupleList_dict after =",
    #for k,v in iAAIG_iminus1aaTypesProbTupleList_dict.items():1
    #    print k,v
    
###########################################################################################   
#                               WRITE MODIFIED FILES                                      #
###########################################################################################

# Write modifies Connectivities Files
outfile = args.POOL_CONNECTIVITIES_FILE+".mod"
if args.POOL_CONNECTIVITIES_FILE == args.COMPLETE_CONNECTIVITIES_FILE:
    outfile = args.POOL_CONNECTIVITIES_FILE+".pool.mod"   # avoid ovewriting the complete aa types file
with open(outfile, 'w') as f:
    f.write("i\tpossible i-1\n")
    for k,v in list(i_iminus1_pool_dict.items()):
        #print k,v
        f.write(k+"\t"+', '.join(str(e) for e in v)+"\n")

with open(args.COMPLETE_CONNECTIVITIES_FILE+".mod", 'w') as f:
    f.write("i\tpossible i-1\n")
    for k,v in list(i_iminus1_complete_dict.items()):
        #print k,v
        f.write(k+"\t"+', '.join(str(e) for e in v)+"\n")

# Save ALL amino acid type predictions with the respective probabilities
#with open("amino_acid_type_prediction_probabilities.mcutoff"+str(args.RESONANCE_MATCH_CUTOFF)+".acutoff"+str(args.ASSIGNMENT_CUTOFF)+".zcutoff"+str(args.ZSCORE_ASSIGNMENT_CUTOFF)
          #+".wH"+str(args.H_weight)+".wC"+str(args.C_weight)+".tolH"+str(args.tolH)+".tolC"+str(args.tolC)+".rtolH"+str(args.rtolH)+".rtolN"+str(args.rtolN), 'w') as f:
with open(args.COMPLETE_AA_TYPES_FILE+".mod", 'w') as f:
    f.write("i AAIG\tpossible i-1 aa types\n")
    #print "Probabilities:"
    for i_AAIG in list(iAAIG_iminus1aaTypesProbTupleList_dict.keys()):
        #print i_AAIG,"---> (i-1) aa type",sorted(iAAIG_iminus1aaTypesProbTupleList_dict[i_AAIG], key=itemgetter(1), reverse=True)
        f.write(i_AAIG + " ---> (i-1) aa type " + str(sorted(iAAIG_iminus1aaTypesProbTupleList_dict[i_AAIG], key=itemgetter(1), reverse=True)) + "\n")

# Save amino acid type predictions with the respective Z-scores
#with open("amino_acid_type_prediction_Z-scores.mcutoff"+str(args.RESONANCE_MATCH_CUTOFF)+".acutoff"+str(args.ASSIGNMENT_CUTOFF)+".zcutoff"+str(args.ZSCORE_ASSIGNMENT_CUTOFF)+
          #".wH"+str(args.H_weight)+".wC"+str(args.C_weight)+".tolH"+str(args.tolH)+".tolC"+str(args.tolC)+".rtolH"+str(args.rtolH)+".rtolN"+str(args.rtolN), 'w') as f:
outfile = args.POOL_AA_TYPES_FILE+".mod"
if args.POOL_AA_TYPES_FILE == args.COMPLETE_AA_TYPES_FILE:
    outfile = args.POOL_AA_TYPES_FILE+".pool.mod"   # avoid overwriting the complete aa types file
with open(outfile, 'w') as f:
    f.write("i AAIG\tpossible i-1 aa types\n")
    #print "Z-scores:"
    for i_AAIG in list(iAAIG_iminus1aaTypesProbPoolTupleList_dict.keys()):
        #print i_AAIG,"---> (i-1) aa type",sorted(iAAIG_iminus1aaTypesProbPoolTupleList_dict[i_AAIG], key=itemgetter(1), reverse=True)
        f.write(i_AAIG + " ---> (i-1) aa type " + str(sorted(iAAIG_iminus1aaTypesProbPoolTupleList_dict[i_AAIG], key=itemgetter(1), reverse=True)) + "\n")