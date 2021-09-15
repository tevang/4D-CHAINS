import sys, re, os, pickle, traceback, shutil, bz2, math, gc, collections
import scoop
import numpy as np
from operator import itemgetter
from collections import OrderedDict
from ete3 import Tree
import ftplib
from argparse import ArgumentParser
from scipy.stats.mstats import zscore
from scipy import stats, sqrt
from cluster import HierarchicalClustering
from lib.global_func import *

def form_chains(AAIG_signature_chains, peak_connectivities_mdict, i_iminus1_complete_dict, spectrum):
    # TODO: under development.
    for chain in AAIG_signature_chains:
        AAIG_list = [spectrum.__get_AAIG__(AAIG_signature) for AAIG_signature in AAIG_signature_chains]
        matching_peak_pairs = [[peak_connectivities_mdict[chain[i]][chain[i+1]]] for i in range(len(AAIG_list)-1)]
        # discard the peak distances
        matching_peak_pairs = [[tuple(triplet[:-1]) for triplet in peak_pairs] for peak_pairs in matching_peak_pairs]
        occupancies, intersections, multi_scores = [], [], []
        for i in range(len(AAIG_list) - 1):
            AAIG1_signature, AAIG2_signature = chain[i], chain[i+1]
            for quintuplet in i_iminus1_complete_dict[AAIG1_signature]:
                if quintuplet[0] == AAIG2_signature:
                    occupancies.append(quintuplet[1]/quintuplet[2])
                    intersections.append(quintuplet[3])
                    multi_scores.append(quintuplet[4])
        yield Chain(AAIG_list, matching_peak_pairs, occupancies, intersections, multi_scores)

class Chain():

    def __init__(self,
                 AAIG_list=[],
                 matching_peak_pairs=[],
                 occupancies=[],
                 intersections=[],
                 multi_scores=[],
                 i_iminus1_dict={},
                 i_iminus1_complete_dict={},
                 all_chainScore_set=set()):
        """

        :param AAIG_list:
        :param matching_peak_pairs:
        :param occupancies:
        :param intersections:
        :param multi_scores:
        :param i_iminus1_dict: the dict with selected connectivities that pass the filters of this cycle/iteration
        :param i_iminus1_complete_dict: this must be the dict with ALL POSSIBLE CONNECTIVITIES in order to calculate correct
                                        chain probabilities.
        :param all_chainScore_set:
        """
        self.nodes = AAIG_list
        self.bonds = matching_peak_pairs
        self.occupancies = occupancies
        self.intersections = intersections
        self.multi_scores = multi_scores
        self.i_iminus1_dict = i_iminus1_dict
        self.i_iminus1_complete_dict = i_iminus1_complete_dict
        self.all_chainScore_set = all_chainScore_set     # a set of lists containing the possible connectivities and the overall score as the last element
        self.i_iminus1_normProbabilities_dict = {}  # dict of the form AAIG1->[(AAIG2, probability), (AAIG2, probability), ...];

    def __reverse__(self):
        """
        Reverse the chain object. This is useful to create NOESY chains, where the direction of each connectivity
        is unknown (may be i->i+1 or i-1->i). E.g. if the original chain was:
        nodes = [N45NH, R46NH, F47NH]
        bonds = [[(peakA1, peakB1), (peakA2, peakB2), (peakA3, peakB3)],
                [(peakB1, peakC1), (peakB2, peakC2), (peakB3, peakC3)],
                [(peakC1, peakD1), (peakC2, peakD2)]
                ]
        The reversed chain will be:
        nodes = [F47NH, R46NH, N45NH]
        bonds = [[(peakD1, peakC1), (peakD2, peakC2)],
                [(peakC1, peakB1), (peakC2, peakB2), (peakC3, peakB3)],
                [(peakB1, peakA1), (peakB2, peakA2), (peakB3, peakA3)]
                ]

        :return:
        """
        self.nodes = list(reversed(self.nodes))
        self.bonds = [[tuple(reversed(peakpair)) for peakpair in bond] for bond in reversed(self.bonds)]
        self.occupancies = list(reversed(self.occupancies))
        self.intersections = list(reversed(self.intersections))
        self.multi_scores = list(reversed(self.multi_scores))


################################################# OLD FUNCTIONS TO FORM CHAINS FROM CONNECTIVITIES ########################################################

    def populate_leaves(self, Assignment_Tree):
        """
        Method that adds new branches to the leaves of the Tree.

        :param Assignment_Tree: the Tree structure with connectivities
        :param i_iminus1_dict:
        :return: (Assignment_Tree, BOOLEAN):    A tuple with elements the input Tree structure with new branches (if applicable), and a BOOLEAN value which is True if the function added
                                           new leaves to the Tree, or False otherwise
        """

        number_of_new_leaves = 0
        intersection, multi = None, None  # default values (in case of old connectivities format)
        # ATTENTION: never use Assignment_Tree.iter_leaf_names(), it doesn't return the names in the order
        # ATTENTION: corresponding to Assignment_Tree.get_leaves()!!!
        for leaf in Assignment_Tree.get_leaves():
            try:
                # for child_tuple in i_iminus1_dict[name]:
                for child_triplet, child_duplet in zip(self.i_iminus1_dict[leaf.name], self.i_iminus1_normProbabilities_dict[leaf.name]):  # WARNING: i_iminus1_normProbabilities_dict contains only matches above the Z-score cutoff
                    NOESYaaindex = child_duplet[0]
                    norm_probability = child_duplet[1]  # normalized probability of the connectivity
                    _occupancy = child_triplet[1]  # add prefix "_" to discriminate from new child feature "occupancy" that will be added
                    _numOfResonances = child_triplet[2]
                    if len(child_triplet) == 5:
                        intersection, multi = child_triplet[3], child_triplet[4]
                    ancestors_list = [ancestor.name for ancestor in leaf.get_ancestors()]
                    if NOESYaaindex in ancestors_list:  # if the current NOESY AAIG is already a node or leaf in the Tree, continue to the next
                        continue
                    # new_child = leaf.add_child(name=NOESYaaindex, dist=float(_occupancy)/_numOfResonances) # add a new brach to the current TOCSY add index (leaf) with length the ratio occupancy/numOfResonances

                    new_child = leaf.add_child(name=NOESYaaindex,
                                               dist=norm_probability)  # add a new brach to the current TOCSY add index (leaf) with length the respective probability

                    new_child.add_features(occupancy=_occupancy,
                                           numOfResonances=_numOfResonances,
                                           intersection=intersection,
                                           multi=multi)
                    number_of_new_leaves += 1
                    # print "DEBUG: adding connection: ",name,"-->",NOESYaaindex
            except KeyError:
                continue

        # print Assignment_Tree.get_ascii(show_internal=True, compact=False, attributes=["name", "dist", "occupancy", "numOfResonances"])
        # print Assignment_Tree.get_ascii(show_internal=True, compact=False)
        if number_of_new_leaves > 0:
            return (Assignment_Tree, True)
        else:
            return (Assignment_Tree, False)


    def calc_chain_probs(self):
        """
        Method to convert occupancies to probabilities. It populates self.i_iminus1_normProbabilities_dict.
        """

        # eventually will contain only the matches above the cutoff
        AAIG_weightSum_dict = {}
        for TOCSY_AAIG in list(self.i_iminus1_complete_dict.keys()):
            weight_sum = 0
            for quintuplet in self.i_iminus1_complete_dict[TOCSY_AAIG]:  # use all possible matches to calculate the weight_sum
                weight = quintuplet[4]  # multi score
                weight_sum += weight
            AAIG_weightSum_dict[TOCSY_AAIG] = weight_sum

        for TOCSY_AAIG in list(self.i_iminus1_dict.keys()):  # use only the matches above the cutoff to save their probabilities
            # print "DEBUG: TOCSY_AAIG=", TOCSY_AAIG, "AAIG_weightSum_dict[TOCSY_AAIG] =",AAIG_weightSum_dict[TOCSY_AAIG]
            for quintuplet in self.i_iminus1_dict[TOCSY_AAIG]:
                iminus1_AAIG_name = quintuplet[0]
                weigh = quintuplet[4]  # the multi score
                duplet = (iminus1_AAIG_name, weight / AAIG_weightSum_dict[TOCSY_AAIG])  # recover the probability by dividing by the weight_sum
                try:
                    self.i_iminus1_normProbabilities_dict[TOCSY_AAIG].append(duplet)
                except KeyError:
                    self.i_iminus1_normProbabilities_dict[TOCSY_AAIG] = [duplet]

    def build_Chain_Tree(self, i, MAX_CHAIN_LENGTH, MIN_CHAIN_LENGTH):
        """

        :param i:
        :param MAX_CHAIN_LENGTH:
        :param MIN_CHAIN_LENGTH:
        :return:
        """
        print("Building Tree starting from amino acid index", i, "...")
        expand_tree = True
        Assignment_Tree = Tree()
        Root = Assignment_Tree.get_tree_root()
        Root.add_feature("name", i)
        level = 1
        sys.stdout.write("Expanding tree from level ")
        while expand_tree:
            sys.stdout.write(str(level) + " ")
            sys.stdout.flush()
            Assignment_Tree, expand_tree = self.populate_leaves(Assignment_Tree)
            level += 1
            if level == MAX_CHAIN_LENGTH:
                break
        if level < MIN_CHAIN_LENGTH:  # discard chains with a single AAIG
            return
        # Print the Tree
        # print Assignment_Tree.get_ascii(show_internal=True, compact=False)
        # print Assignment_Tree.get_ascii(show_internal=True, compact=False, attributes=["name", "dist", "occupancy", "numOfResonances"])

        print("\nSaving chains from Tree...")

        for leaf in Assignment_Tree.get_leaves():
            chain = []
            score = leaf.dist
            chain.append(leaf.name)
            for ancestor in leaf.get_ancestors():
                # if ancestor == Root:
                #     break
                # TODO: check why it doesn't add the Root. Is it because it doesn't have property 'dist'?
                chain.append(ancestor.name)
                score *= ancestor.dist  # follow the chain rule for conditional probabilities to calculate the score (probability)
            # ATTENTION: the Tree is pruned in the sense that only connectivities above the Z-score cutoff were used,
            # otherwise the scores of all chains would sum to 1 !
            # Therefore there is no need for normalization of chain scores!
            if MIN_CHAIN_LENGTH <= len(chain) <= MAX_CHAIN_LENGTH:  # save only chains of valid length
                chain.append(score)
                self.all_chainScore_set.add(tuple(chain))
            del chain
            del ancestor
            del leaf
            # Assignment_Tree = None
        del Assignment_Tree
        gc.collect()

    def form_chains(self, MAX_PEPTIDE_LENGTH, MIN_PEPTIDE_LENGTH):
        """
        Main method the create all possible chains from a connectivity dictionary.

        :param i_iminus1_dict:  the connectivity dictionary. It may not contain all possible connectivities.
        :param MAX_PEPTIDE_LENGTH:
        :param MIN_PEPTIDE_LENGTH:
        :return:
        """
        self.calc_chain_probs() # populate self.i_iminus1_normProbabilities_dict
        for i in list(self.i_iminus1_dict.keys()):
            self.build_Chain_Tree(i, MAX_PEPTIDE_LENGTH, MIN_PEPTIDE_LENGTH)

    def find_redundant_chains(self, chainScore_list):
        """
        Method that finds the redundant chains to remove from the pool of chains.

        :param chainScore_list:
        :return:
        """
        progbar = scoop.shared.getConst('PROGBAR')
        all_chainScore_list = scoop.shared.getConst('ALL_CHAINSCORE_LIST')
        chainScores2remove_set = set()
        chain_num = len(chainScore_list)
        for i, chainScore1 in enumerate(chainScore_list):
            chain1 = chainScore1[0:-1]
            if len(set(chain1)) < len(chain1):  # this chains contains multiple times the same AAIG!
                chainScores2remove_set.add(tuple(chainScore1))
                continue    # this will be removed, so no reason to look for overlapping chains
            score1 = chainScore1[-1]
            progbar.set_progress((float(i) / chain_num))
            for chainScore2 in all_chainScore_list:
                chain2 = chainScore2[0:-1]
                score2 = chainScore2[-1]
                if is_sublist(chain1, chain2) and score1 == score2:  # this is the correct one
                    chainScores2remove_set.add(tuple(chainScore1))  # make first the list immutable to be able to insert it into a set

        return chainScores2remove_set

    def remove_redundant_chains(self):
        """

        :return: all_chainScore_list: E.g. [ ('V126NH', 'S99NH', 'A111NH', 'T128NH', 7.187111855209132e-13), ...]
        """
        all_chainScore_list = list(self.all_chainScore_set)  # convert the set to a list to keep track of the sequence that each score belongs to
        self.all_chainScore_set = set()  # delete it to save memory
        progbar = ProgressBar(100)
        all_chainScore_listList_list = chunkIt(all_chainScore_list, scoop.SIZE)
        try:
            scoop.shared.setConst(PROGBAR=progbar)
        except TypeError:
            pass
        scoop.shared.setConst(ALL_CHAINSCORE_LIST=all_chainScore_list)
        results = list(scoop.futures.map(self.find_redundant_chains, all_chainScore_listList_list))
        # results is a list of sets with chains to be removed
        chainScores2remove_set = set()
        for s in results:
            chainScores2remove_set = chainScores2remove_set.union(s)

        # remove subchains to keep only the longest, unique, non-overlapping ones
        for chainScore in chainScores2remove_set:
            all_chainScore_list.remove(tuple(chainScore))

        return all_chainScore_list

    def __add_reversed_chains__(self):
        self.all_chainScore_set = self.all_chainScore_set.union( set(Chain.reverse_chains(self.all_chainScore_set)) )

    def write_chains_file(self, all_chainScore_set, outfname):
        """

        :param all_chainScore_set: can be both a set or a list
        :param outfname:
        :return:
        """
        all_chainScore_list = list(all_chainScore_set)
        all_chainScore_list.sort(key=itemgetter(-1), reverse=True)
        with open(outfname, 'w') as f:
            for chainScore_list in all_chainScore_list:
                chain = chainScore_list[0:-1]
                score = chainScore_list[-1]
                # print "Overall Score = ",score," chain = ", chain
                f.write("Overall Score = " + str(score) + " chain = " + str(chain) + "\n")

    @staticmethod
    def load_chains_file(fname):
        all_chainScore_list = []
        with open(fname, 'r') as f:
            for line in f:
                word_list = line.split()
                chain_score = float(word_list[3])
                values_string = ''.join(word_list[6:])
                Tindices_string = re.sub('[\(\)\'\"]', '', values_string).split("),(")[0]
                chain = Tindices_string.split(",")
                chain.append(chain_score)  # append the chain score to the end of the chain list
                all_chainScore_list.append( tuple(chain) )
        return all_chainScore_list

    @staticmethod
    def reverse_chains(all_chainScore_list):
        all_rev_chainScore_list = []
        for chainScore in all_chainScore_list:
            chain = chainScore[:-1]
            score = chainScore[-1]
            rev_chainScore = list(reversed(chain))
            rev_chainScore.append(score)
            all_rev_chainScore_list.append( tuple(rev_chainScore) )
        return all_rev_chainScore_list

#####################################################################################################################################