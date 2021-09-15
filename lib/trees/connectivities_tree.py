from operator import itemgetter

from ete3 import Tree
from itertools import combinations, permutations
import numpy as np
import gc
from copy import deepcopy
from lib.global_func import *
from lib.proofread.print_functions import *

from scoop import shared, futures


class Connectivities_Tree(object):

    def __init__(self, connections_dict):
        """

        :param connections_dict: TOCSY peak -> [(NOESY peak, distance), (NOESY peak, distance), ...]
        """
        # IMPORTANT: delete AAIG1peaks without matched AAIG2peaks, because they block the expansion of the branches of the tree
        for k,v in list(connections_dict.items()):
            if len(v) == 0:
                del connections_dict[k]
        self.connections_dict = connections_dict

    def get_nonredundant_permutations(self):
        """
        WRONG!!! It misses some non-redundant permutations.
        :return:
        """
        AAIG1peaks_with_common_AAIG2peaks = set()
        for comb in combinations(list(self.connections_dict.keys()), 2):
            AAIG1peak1, AAIG1peak2 = comb
            if AAIG1peak1 in AAIG1peaks_with_common_AAIG2peaks and AAIG1peak2 in AAIG1peaks_with_common_AAIG2peaks:
                continue
            for AAIG2peak1 in self.connections_dict[AAIG1peak1]:
                if AAIG2peak1 in self.connections_dict[AAIG1peak2]:
                    AAIG1peaks_with_common_AAIG2peaks.add(AAIG1peak1)
                    AAIG1peaks_with_common_AAIG2peaks.add(AAIG1peak2)
                    break

        if len(AAIG1peaks_with_common_AAIG2peaks) == 0:
            return [list(self.connections_dict.keys())]

        all_AAIG1peaks = set(self.connections_dict.keys())
        fixed_AAIG1peaks = list(all_AAIG1peaks.difference(AAIG1peaks_with_common_AAIG2peaks))
        valid_permutations = []
        AAIG1peaks_with_common_AAIG2peaks.add('fixed_AAIG1peaks')  # 'fixed_AAIG1peaks' designates the position of the block of fixed AAIG1peaks in the permutation
        for perm in permutations(AAIG1peaks_with_common_AAIG2peaks, len(AAIG1peaks_with_common_AAIG2peaks)+1):
            insert_indx = perm.index('fixed_AAIG1peaks')
            valid_permutations.append(perm[:insert_indx] + fixed_AAIG1peaks + perm[insert_indx+1:])
        return valid_permutations

    def expand_tree(self, perm, length, connections_dict):
        """
        Recall that we do AAIG1-->AAIG2 peak matching to find connectivities.
        AAIG1 can be TOCSY (in case of TOCSY-HCNH) matching, or NOESY (in case of HCNH-HCNH matching).
        AAIG2 can only be NOESY.
        :param perm: a tuple of AAIG1 peaks which will be used
        :param length:  max chain length, namely number of AAIG1 peaks to match with AAIG2 peaks.
        :param connections_dict:
        :return:
        """
        # print("DEBUG: building new tree using AAIG1peak permutation:", [p.__get_CHresons__() for p in perm])
        expand_tree = True
        Connectivity_Tree = Tree()
        Root = Connectivity_Tree.get_tree_root()
        Root.add_features(name="root", distance=0.0)
        level = 0
        for AAIG1peak in perm:
            Connectivity_Tree, expand_tree = self.populate_leaves(Connectivity_Tree, AAIG1peak, connections_dict)
            if expand_tree == False:
                break
            level += 1
            if level == length:
                break

        # print "\nSaving chains from Tree..."
        for leaf in Connectivity_Tree.get_leaves():
            chain = []  # list of matching NOESY Peak objects
            euclidean_dist = leaf.distance ** 2
            chain.append((leaf.AAIG1peak, leaf.name))  # (AAIG1 Peak object, matching AAIG2 Peak object)
            # print("DEBUG: leaf.name=", leaf.name)
            for ancestor in leaf.get_ancestors():
                if ancestor.name == 'root':
                    continue
                # print("DEBUG: ancestor.name=", ancestor.name)
                chain.append((ancestor.AAIG1peak, ancestor.name))
                euclidean_dist += ancestor.distance ** 2
            euclidean_dist = np.sqrt(euclidean_dist) / length
            chain.append(euclidean_dist)
            # print("DEBUG: saving chain=", [(pp[0].__get_CHresons__(), pp[1].__get_CHresons__()) for pp in chain[:-1]],
            #       chain[-1])
            self.all_chainDist_list.append(tuple(chain))
            del chain
            del ancestor
            del leaf

    def expand_tree_parallel(self, perm, length, connections_dict):
        # print "DEBUG: building new tree using AAIG1peak permutation:", perm
        expand_tree = True
        Connectivity_Tree = Tree()
        Root = Connectivity_Tree.get_tree_root()
        Root.add_features(name="root", distance=0.0)
        level = 0
        for AAIG1peak in perm:
            Connectivity_Tree, expand_tree = self.populate_leaves(Connectivity_Tree, AAIG1peak, connections_dict)
            if expand_tree == False:
                break
            level += 1
            if level == length:
                break

        # print "\nSaving chains from Tree..."
        chainDist_list = [] # local list otherwise self.all_chainDist_list won't work in parallel execution!
        for leaf in Connectivity_Tree.get_leaves():
            chain = []  # list of matching NOESY Peak objects
            euclidean_dist = leaf.distance ** 2
            chain.append((leaf.AAIG1peak, leaf.name))  # (TOCSY Peak object, matching NOESY Peak object)
            # print "DEBUG: leaf.name=", leaf.name
            for ancestor in leaf.get_ancestors():
                if ancestor.name == 'root':
                    continue
                # print "DEBUG: ancestor.name=", ancestor.name
                chain.append((ancestor.AAIG1peak, ancestor.name))
                euclidean_dist += ancestor.distance ** 2
            euclidean_dist = np.sqrt(euclidean_dist) / length
            chain.append(euclidean_dist)
            # print "DEBUG: saving chain=", chain
            chainDist_list.append(tuple(chain))
            del chain
            del ancestor
            del leaf
        return chainDist_list

    def form_chains(self, length='max'):
        if length == 'max':
            length = len(list(self.connections_dict.keys()))

        self.all_chainDist_list = []
        # print("DEBUG: list(self.connections_dict.keys())=", [p.__get_CHresons__() for p in self.connections_dict.keys()])
        # print("DEBUG: self.connections_dict=", self.connections_dict)
        # for perm in self.get_nonredundant_permutations():     # WRONG!!!
        for perm in permutations(list(self.connections_dict.keys()), len(list(self.connections_dict.keys()))):
            # print("DEBUG: perm=", perm)
            self.expand_tree(perm, length, self.connections_dict)

    def form_chains_parallel(self, length='max'):
        if length == 'max':
            length = len(list(self.connections_dict.keys()))

        self.all_chainDist_list = []
        all_perms = list(permutations(list(self.connections_dict.keys()), len(list(self.connections_dict.keys()))))

        # Parallel version
        results = list(futures.map(self.expand_tree_parallel,
                                   all_perms,
                                   [length]*len(all_perms),
                                   [self.connections_dict]*len(all_perms)
                                   ))
        self.all_chainDist_list = [r[0] for r in results]
        # print "DEBUG: self.all_chainDist_list=", self.all_chainDist_list

        # # Serial version for Debugging
        # for perm in all_perms:
        #     self.expand_tree(perm, length, self.connections_dict)

    def is_AAIG2peak_in_ancenstors(self, leaf, AAIG2peak):
        ancestors_AAIG2peak_list = [ancestor.name for ancestor in
                                    leaf.get_ancestors()]  # list with all the  currently in the branch starting from this leaf
        ancestors_AAIG2peak_list.append(leaf.name)
        return AAIG2peak in ancestors_AAIG2peak_list

    def populate_leaves(self, Connectivity_Tree, AAIG1peak, connections_dict):
        number_of_new_leaves = 0
        ## ATTENTION: IN PYTHON 3 ALWAYS WHEN YOU ITERATE OVER A CONTAINER AND ADD ELEMENTS IN EACH ITERATION USE A COPY!!!
        ## ATTENTION: NEVER USE Tree.iter_leaves() AGAIN IN PYTHON 3 BECAUSE IT ITERATES OVER EVERY NEW LEAF YOU ADD!!!
        ## ATTENTION: USE INSTEAD Tree.get_leaves() WHICH RETURNS A CONSTANT LIST OF LEAVES.
        leaves = Connectivity_Tree.get_leaves()
        # print("DEBUG: leaves=", leaves)
        for leaf in Connectivity_Tree.get_leaves():
            # if leaf.name== 'root':
            #     print("DEBUG: processing leaf:", leaf.name)
            # else:
            #     print("DEBUG: processing leaf:", leaf.name, leaf.AAIG1peak, leaf.distance)
            try:
                for AAIG2peak, distance in connections_dict[AAIG1peak]:
                    if self.is_AAIG2peak_in_ancenstors(leaf, AAIG2peak):  # if this AAIG2peak is already in the Tree, skip it
                        continue
                    new_child = leaf.add_child(name=AAIG2peak)
                    # print("DEBUG: added new child leaf with features: ", AAIG1peak, AAIG2peak, distance)
                    new_child.add_features(AAIG1peak=AAIG1peak, distance=distance)
                    number_of_new_leaves += 1
            except KeyError:
                continue

        # print(Connectivity_Tree.get_ascii(show_internal=True, compact=False, attributes=["AAIG1peak", "name"]))
        # print Connectivity_Tree.get_ascii(show_internal=True, compact=False)
        if number_of_new_leaves > 0:
            return (Connectivity_Tree, True)
        else:
            return (Connectivity_Tree, False)

    def get_best_chain(self, length='max', maxdist=False):
        """
        Method to find the longest chain with the lowest C,H distance, which represents the maximum possible number of
        matching TOCSY-HCNH peak pairs (in case of TOCSY-HCNH AAIG matching) or matching HCNH-HCNH peak pairs
        (in case of HCNH-HCNH AAIG matching) for a particular AAIG1(TOCSY or NOESY)-AAIG2(NOESY) combination.

        :param length: number of AAIG1 peaks to use
        :param maxdist: if True then is will keep the longest chain with the highest distance
        :return: best_chain:    the optimum list of [(TOCSY Peak object, matching NOESY Peak object), ...] (without the distance)
                                that match with the Peak objects of the current TOCSY AAIG.
        """
        ColorPrint("Searching for best chain", "OKBLUE")
        if length == 'max':
            length = len(list(self.connections_dict.keys()))
            # Debuginfo("DEBUG: list(self.connections_dict.keys()=%s" % list(self.connections_dict.keys()))
        # Debuginfo("DEBUG: length=%i" % length)

        # print("DEBUG: self.all_chainDist_list=", self.all_chainDist_list)
        # NOTE: self.all_chainDist_list = [((AAIG1peakA, matching AAIG2peak), (AAIG1peakB, matching AAIG2peak),..., total C-H distance),
        #                                   ...]
        if len(self.all_chainDist_list) == 0:
            raise Exception("ERROR: empty self.all_chainDist_list!!!")
        valid_chains = []
        while len(valid_chains) == 0:
            assert length > 0, Debuginfo("FAIL: infinite while loop!", fail=True)
            # valid_chains = [chain for chain in self.all_chainDist_list if len(chain)==length+1]  # +1 for the distance
            valid_chains = []
            for chain in self.all_chainDist_list:
                # print("DEBUG: chain=", chain, "len(chain)=", len(chain), "length+1=", length+1)
                if len(chain)==length+1:  # +1 for the distance
                    valid_chains.append(chain)
            valid_chains.sort(key=itemgetter(length), reverse=maxdist)
            length -= 1

        # print("DEBUG: valid_chains=", valid_chains)
        # print("DEBUG: valid_chains=", [([(unfold_peak(sc[0]), unfold_peak(sc[1])) for sc in c[:-1]], c[-1]) for c in valid_chains])
        return valid_chains[0][:-1]     # without the distance