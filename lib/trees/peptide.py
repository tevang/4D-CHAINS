#!/usr/bin/env python

import sys, re, os, pickle, traceback, shutil, bz2, math

from lib.open_func import bcolors
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
from lib.global_func import *
from lib.global_vars import aa3to1_dict
import scoop


class Peptide():

    def __init__(self,
                 total_chain_number,
                 min_peptide_length,
                 AAIG_aaType_Condprob_mdict,
                 AAIG_aaTypesCutoffProbTupleList_dict,
                 AAIG_aaTypesTransformedProbTupleList_dict,
                 aatype_P_dict,
                 connectivity_dict,
                 spectrum_type="TOCSY",
                 workdir="."):    # connectivity_dict is the old "i_iminus1_dict"
        """

        :param total_chain_number:
        :param min_peptide_length:
        :param AAIG_aaType_Condprob_mdict:  all CONDITIONAL probabilities. These are the probabilities that are used for
                                                peptide score calculation. Must be re-calculated in every iteration/cycle
                                                because the number of unassigned residues varies.
        :param AAIG_aaTypesCutoffProbTupleList_dict:    used only to filter the conditional probabilities.
        :param AAIG_aaTypesTransformedProbTupleList_dict:   probabilities after a prespecified transformation (e.g. log). I saw it only used in statistics.py to
                                                            calculate the conditional probabilities.
        :param aatype_P_dict:
        :param connectivity_dict:
        :param spectrum_type:
        :param workdir:
        """
        self.spectrum_type = spectrum_type
        self.peptideScoreChainindex_list = []
        self.peptideProbList_list = []
        self.workdir = workdir

        # SET SHARED VARIABLES BETWEEN THREADS
        shared.setConst(TOTAL_CHAIN_NUMBER=total_chain_number)
        shared.setConst(SPECTRUM_COMBO=spectrum_type)
        shared.setConst(MIN_PEPTIDE_LENGTH=min_peptide_length)
        # print("DEBUG: = AAIG_aaTypesCutoffProbTupleList_dict=", AAIG_aaTypesCutoffProbTupleList_dict)
        # print("DEBUG: AAIG_aaType_Condprob_mdict=", AAIG_aaType_Condprob_mdict)
        shared.setConst(IAAINDEX_IMINUS1AATYPE_CONDPROB_MULTIDICT=AAIG_aaType_Condprob_mdict)
        shared.setConst(IAAINDEX_IMINUS1AATYPESCUTOFFPROBTUPLELIST_DICT=AAIG_aaTypesCutoffProbTupleList_dict)
        shared.setConst(AA3TO1_DICT=aa3to1_dict)
        # print("DEBUG: AAIG_aaTypesTransformedProbTupleList_dict=",
        #       AAIG_aaTypesTransformedProbTupleList_dict)
        shared.setConst(
            IAAINDEX_IMINUS1AATYPESTRANSFORMEDPROBTUPLELIST_DICT=AAIG_aaTypesTransformedProbTupleList_dict)
        # shared.setConst(IAAINDEX_IMINUS1AATYPESPAATYPETUPLELIST_DICT=iAAIG_iminus1aaTypesPaatypeTupleList_dict)
        shared.setConst(AATYPE_P_DICT=aatype_P_dict)
        # print("DEBUG: connectivity_dict=", connectivity_dict)
        shared.setConst(I_IMINUS1_DICT=connectivity_dict)

######################################### OLD FUNCTIONS TO FORM PEPTIDES FROM ONE CHAIN #########################################


    def populate_peptideTree_leaves(self, Peptide_Tree, iAAIG):
        """
            FUNCTION that adds new branches to the leaves of the Tree.

            ARGUMENTS:
            Peptide_Tree:   The Tree structure with connectivities
            iAAIG:       The AAIG of the current element of the chain
            RETURNS:
            (Peptide_Tree, BOOLEAN):    A tuple with elements the input Tree structure with new branches (if applicable), and a BOOLEAN value which is True if the function added
                                           new leaves to the Tree, or False otherwise
        """
        AAIG_aaType_Condprob_mdict = shared.getConst('IAAINDEX_IMINUS1AATYPE_CONDPROB_MULTIDICT')
        AAIG_aaTypesCutoffProbTupleList_dict = shared.getConst('IAAINDEX_IMINUS1AATYPESCUTOFFPROBTUPLELIST_DICT')
        aa3to1_dict = shared.getConst('AA3TO1_DICT')
        #iminus1aa_iAAIGTypesProbTupleList_dict = shared.getConst(iminus1aa_iAAIGTypesProbTupleList_dict=IMINUS1AA_IAAINDEXTYPESPROBTUPLELIST_DICT)
        #iAAIG_iminus1aaTypesPaatypeTupleList_dict = shared.getConst(iAAIG_iminus1aaTypesPaatypeTupleList_dict=IAAINDEX_IMINUS1AATYPESPAATYPETUPLELIST_DICT)
        aatype_P_dict = shared.getConst('AATYPE_P_DICT')

        number_of_new_leaves = 0
        for leaf in Peptide_Tree.get_leaves():
            if leaf.name == "P":    # proline can only be the root of the tree, not a leaf!
                continue
            leaf.add_feature("AAIG", iAAIG)
            try:
                # print("DEBUG: AAIG_aaTypesCutoffProbTupleList_dict[iAAIG]=", AAIG_aaTypesCutoffProbTupleList_dict[iAAIG])
                for child_tuple in AAIG_aaTypesCutoffProbTupleList_dict[iAAIG]:
                    # print("DEBUG: adding child_tuple=", child_tuple)
                    aa_type = child_tuple[0]    # i-1 aa type
                    # Paaiminus1_TAAIGi = calculate_Paaiminus1_TAAIGi(iAAIG, aa_type)    # OBSOLE: calculate the conditional probability on the fly
                    Paaiminus1_TAAIGi = AAIG_aaType_Condprob_mdict[iAAIG][aa_type]
                    try:
                        _ = leaf.add_child(name=aa3to1_dict[aa_type], dist=Paaiminus1_TAAIGi) # add a new brach to the current aa type (leaf) with length the respective P(TAAIGi) weighted by
                                                                                              # the probability of the current transition in the current chain
                    except TypeError:   # when no Paaiminus1_TAAIGi probability could be calculated, don't add this new leaf!
                                        # OCCURS ONLY WHEN RESTRAINTS WERE IMPOSED VIA modify_connectivities_and_aatypes.py SCRIPT
                        # print("DEBUG: Paaiminus1_TAAIGi=", Paaiminus1_TAAIGi)
                        continue
                    #print "DEBUG: adding connection: ",leaf.name,"-->",aa_type
                    number_of_new_leaves += 1
            except KeyError:
                continue

        if number_of_new_leaves > 0:
            return (Peptide_Tree, True)
        else:
            return (Peptide_Tree, False)
        # print Peptide_Tree.get_ascii(show_internal=True, compact=False, attributes=["name", "dist"])
        # print Peptide_Tree.get_ascii(show_internal=True, compact=False)


    def build_Peptide_Tree(self, chainIndex, chainScore):
        """
            Method that takes in one chain and builds and save into files all the possible peptides. Every time you run
            this method, the previous peptide list will be erased but the pickle files will remain.

        :param chainIndex: a unique integer number to be appended at the end of the file which will contain the peptides
        :param chainScore: list containing the AAIGs that consist the chain, and as the last element the frequency (score) of the chain.
                            E.g. ('V126NH', 'S99NH', 'A111NH', 'T128NH', 7.187111855209132e-13)
        :return:
        """
        # geometric progression (see my How_to_run_HRex_with_PLUMED.sh file) ??? The last 2 number were added by intuition cause I
        # couldn't find the function that I used to generate the first 7 numbers!
        resonNumWeight_list = [1.000000, 2.656854, 3.928203, 5.000000, 5.944272, 6.797959, 7.583005, 8.3413055, 8.50813] # max 9 peaks (LYS)

        try:
            total_chain_number = shared.getConst('TOTAL_CHAIN_NUMBER')
            spectrum_combo = shared.getConst('SPECTRUM_COMBO')
            i_iminus1_dict = shared.getConst('I_IMINUS1_DICT')
            min_peptide_length = shared.getConst('MIN_PEPTIDE_LENGTH')
            chain = list(chainScore[0:-1])
            score_of_chain = chainScore[-1]
            print("Predicting possible peptide sequences from chain ("+str(chainIndex+1)+"/"+str(total_chain_number)+"):  ",chain,"...")

            self.peptideScoreChainindex_list = []   # remove all previously formed peptides
            self.peptideProbList_list = []          # remove all previously formed peptides

            Peptide_Tree = Tree()
            Root = Peptide_Tree.get_tree_root()
            iAAIG = chain[-1]    # we start from the C-term of the chains and move twoards the N-term
            Root.add_features(name=iAAIG, AAIG=iAAIG)   # only the root of the Peptide_Tree has the name of the AAIG,
                                                        # the rest of the nodes have the names of possible aa types in the
                                                        # peptide sequence starting from iAAIG

            if len(chain) < min_peptide_length: # if the chain contains fewer aa indices than the specified minimum, do not create the Peptide Tree
                print("Chain ", chain, " has length smaller that the cutoff", min_peptide_length, ". No peptides will be built.")
                return
            sys.stdout.write("Adding branches corresponding to amino acid index ")

            for iAAIG in reversed(chain):
                sys.stdout.write(str(iAAIG)+" ")
                sys.stdout.flush()
                # Print the Tree
                #print Peptide_Tree.get_ascii(show_internal=True, compact=False, attributes=["name", "dist"])
                Peptide_Tree, expand_tree = self.populate_peptideTree_leaves(Peptide_Tree, iAAIG)
                if expand_tree == False:    # if no leaves were added in the last iteration, do not expand the Peptide_Tree anymore
                    break
            print("")    # change line after sys.stdout.write() invokation

            # print("DEBUG: printing peptide tree:")
            # print(Peptide_Tree.get_ascii(show_internal=True, compact=False, attributes=["name", "dist", "AAIG"]))

            for leaf in Peptide_Tree.get_leaves(): # iterate over the bottom nodes (leafs) of the Peptide_Tree, namely the N-teminal ends of the peptides
                peptide = []
                peptide_probability_list = []   # list containing the prediction scores of each aa of the peptide, in the direction N-term -> C-term
                total_peptide_probability = leaf.dist
                peptide.append(leaf.name)
                peptide_probability_list.append(leaf.dist * score_of_chain) # save transformed distance (probability that this aa type is correctly predicted) x the frequency (score) of the source chain
                for ancestor in leaf.get_ancestors():   # multiply the individual probabilities of each prediction to get the total ZscoreRatio that this peptide sequence is correct
                    if ancestor == Root:    # otherwise in Python 3 it will include the Root which is iAAIG, not aa type
                        break
                    peptide.append(ancestor.name)
                    # print("DEBUG: leaf.name=", leaf.name, "--> ancestor.name", ancestor.name,
                    #       "--> ancestor.dist", ancestor.dist)
                    peptide_probability_list.append(ancestor.dist * score_of_chain) # save transformed distance (probability that this aa type is correctly predicted) x the frequency (score) of the source chain
                    total_peptide_probability *= ancestor.dist
                    del ancestor
                # print("DEBUG: peptide=", peptide, "peptide_probability_list=", peptide_probability_list)
                if len(peptide) >= min_peptide_length:   # don't save peptides shorter than the minimum length
                    # print("DEBUG: appending peptide:", peptide)
                    peptide.append(total_peptide_probability * score_of_chain)   # finally multiply the total peptide probability with the score of the source chain
                    peptide.append(chainIndex)  # append the chain index
                    # print("DEBUG: saving peptide=",peptide)
                    # print("DEBUG: saving peptide_probability_list=", peptide_probability_list)
                    self.peptideScoreChainindex_list.append(tuple(peptide))

                    ##~~~~~~~~~~~~~~~~~~~~~~~~ WEIGHTS OF PREDICTION SCORES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    TAAIG_i = chain[-1]   # we move from the leaf to the root, namely from the N-term to the C-term of the chain, so this is the i TOCSY group
                    length = len(peptide_probability_list)  # peptide length
                    # print("DEBUG: i_iminus1_dict[TAAIG_i]=", i_iminus1_dict[TAAIG_i])
                    for j, TAAIG_iminus1 in enumerate(reversed(chain[:-1])):
                        try:
                            # print("DEBUG: TAAIG_iminus1=", TAAIG_iminus1, "i_iminus1_dict[TAAIG_i]=", i_iminus1_dict[TAAIG_i])
                            quintuplet = [t for t in i_iminus1_dict[TAAIG_i] if t[0]==TAAIG_iminus1][0]  # retrieve the connectivity info (occupancy and total num of resonances)
                            if len(quintuplet) == 5:    # score the peptide my the multi value
                                TAAIG_iminus1, occupancy, tot_reson_num, intersection, multi = quintuplet
                                peptide_probability_list[length-(j+2)] *= multi*resonNumWeight_list[tot_reson_num-1]
                            elif len(quintuplet) == 3:    # old connectivities format
                                TAAIG_iminus1, occupancy, tot_reson_num = quintuplet
                                peptide_probability_list[length-(j+2)] *= (occupancy/float(tot_reson_num)) * \
                                                                              resonNumWeight_list[tot_reson_num-1]*occupancy
                                # -1 because counting starts from 0 and another -1 because the C-term aa (root of the tree) does not
                                # have a connectivity
                        except IndexError:
                            if spectrum_combo == "TOCSY" and tot_reson_num-1 > len(resonNumWeight_list):
                                ColorPrint("ERROR: the maximum number of peaks (9) allowed for a TOCSY AAIG has been exceeded!" + \
                                " The TOCSY AAIG " + TAAIG_i + " had " + str(tot_reson_num) + " peaks. Please correct the TOCSY spectrum file " + \
                                "and run 4D-CHAINS again.", "FAIL")
                            pass
                        TAAIG_i = TAAIG_iminus1
                    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                    self.peptideProbList_list.append(peptide_probability_list)
                del leaf
                del peptide
            del Peptide_Tree
            # Save the peptides in a pickle object of the form (peptideScoreChainindex_list, peptideProbList_list)
            if len(self.peptideScoreChainindex_list) == 0 and len(self.peptideProbList_list) == 0:
                ColorPrint("WARNING: saving empty peptide file tmp_peptide_folder/chainIndex_"+str(chainIndex)+'.pickle.bz2',
                           "WARNING")
            save_pickle(self.workdir + '/tmp_peptide_folder/chainIndex_'+str(chainIndex)+'.pickle.bz2',
                        self.peptideScoreChainindex_list,
                        self.peptideProbList_list)
            del self.peptideScoreChainindex_list, self.peptideProbList_list
            gc.collect()

        except:
            type, value, tb = sys.exc_info()
            lines = traceback.format_exception(type, value, tb)
            print((''.join(lines)))
            raise

    def build_all_peptide_trees(self, all_chainScore_list, resume=False, is_parallel=True):

        ColorPrint("Forming peptides in parallel.", "BOLDGREEN")
        # Load previously generated peptides from pickle files
        peptide_file_num = 0
        # while peptide_file_num < total_chain_number:    # CREATE MISSING PEPTIDE TREES
        if os.path.exists(
                self.workdir + '/tmp_peptide_folder/') == True and resume == False:  # if the folder exists and this is not a resumed run, remove it and make a new one
            shutil.rmtree(self.workdir + '/tmp_peptide_folder/', ignore_errors=True)
            os.makedirs(self.workdir + '/tmp_peptide_folder/')
        elif os.path.exists(
                self.workdir + '/tmp_peptide_folder/') == False and resume == False:  # if the folder does not exists, make a new one
            os.makedirs(self.workdir + '/tmp_peptide_folder/')
        elif os.path.exists(self.workdir + '/tmp_peptide_folder/') == False and resume == True:
            print("ERROR: folder tmp_peptide_folder/ with peptide sequence files does not exist! "
                  "The process cannot be resumed!")

        fnames = os.listdir(self.workdir + '/tmp_peptide_folder/')
        fpattern = re.compile('chainIndex_[0-9]+.pickle.bz2')
        peptide_files_list = list(filter(fpattern.search, fnames))
        peptide_file_num = len(peptide_files_list)
        chainIndex_args, chainScore_args = [], []
        for chainIndex, chainScore in enumerate(all_chainScore_list):
            if not 'chainIndex_' + str(chainIndex) + '.pickle.bz2' in peptide_files_list:
                chainIndex_args.append(chainIndex)
                chainScore_args.append(chainScore)

        if is_parallel:
            results = list(futures.map(self.build_Peptide_Tree, chainIndex_args, chainScore_args))  # results is empty
        else:
            for chainIndex, chainScore in zip(chainIndex_args, chainScore_args):
                self.build_Peptide_Tree(chainIndex, chainScore)

    @staticmethod
    def write_peptides_from_pickle_file(peptide_files_list, part, workdir="."):

        progbar = shared.getConst('PROGBAR')
        max_peptide_length = shared.getConst('MAX_PEPTIDE_LENGTH')
        all_chainScore_list = shared.getConst('NEW_ALL_CHAINSCORE_LIST')
        f = open(workdir + "/peptides." + str(max_peptide_length) + "mers.list.chunk" + part, 'w')
        index = 0
        peptide_file_index = 0
        peptide_file_num = len(peptide_files_list)
        for peptide_file in peptide_files_list:
            mo = re.search('chainIndex_([0-9]+).pickle.bz2', peptide_file)
            if mo:
                peptide_file_index += 1
                progbar.set_progress((float(peptide_file_index) / peptide_file_num))
                chainIndex = int(mo.group(1))
                # print("DEBUG: loading tmp_peptide_folder/" + peptide_file)
                try:
                    # peptideScoreChainindexList_list: a list of lists containing the possible peptide sequences and the overall score as the last element
                    # peptideProbList_list: a list of lists containing the probability that each amino acid of the peptide is correctly predicted from the chain
                    peptideScoreChainindexList_list, \
                    peptideProbList_list = load_pickle(workdir + '/tmp_peptide_folder/' + peptide_file)
                except EOFError:
                    Debuginfo("ERROR: file tmp_peptide_folder/" + peptide_file, " seems to be corrupted perhaps due to abrupt termination of the program before the file was properly saved.\
                    Please run again the program without the -resume option.", fail=True)
                    sys.exit(1)
                for peptideScoreChainindexList, peptideProbList in zip(peptideScoreChainindexList_list,
                                                                       peptideProbList_list):
                    chain = all_chainScore_list[chainIndex]
                    peptide = peptideScoreChainindexList[0:-2]
                    Cterm = peptideScoreChainindexList[-3]
                    peptide_score = peptideScoreChainindexList[-2]
                    chainIndex2 = peptideScoreChainindexList[-1]
                    if chainIndex != chainIndex2:
                        Debuginfo("ERROR: the chainIndex of file ", peptide_file, " should be ", chainIndex2, ", not ",
                              chainIndex, fail=True)
                        sys.exit(1)
                    if peptide_score == 0.0:  # do not save peptides containing aa types that are not present in the protein
                        continue
                    # print "DEBUG: chain=", chain, " peptide = ", peptide, "index=", index, "peptideProbList=", peptideProbList
                    f.write("Chain probability = " + str(chain[-1]) + " peptide P_" + str(len(peptide)) + "_" + str(
                        index + 1) + " = " + ''.join(peptide) + " " + '-'.join(chain[0:-1]) + " " + ','.join(
                        map(str, peptideProbList)) + "\n")
                    index += 1
        f.close()

    @staticmethod
    def write_all_peptides_to_files(all_chainScore_list, max_peptide_length, workdir="."):
        ##
        ## Write peptides to files
        ##
        # LOAD ALL WRITEN PEPTIDE FILES
        progbar = ProgressBar(100)
        peptide_file_num = 0
        # while peptide_file_num < total_chain_number:    # CREATE MISSING PEPTIDE TREES
        fnames = os.listdir(workdir + '/tmp_peptide_folder/')
        fpattern = re.compile('chainIndex_[0-9]+.pickle.bz2')
        peptide_files_list = list(filter(fpattern.search, fnames))
        peptide_file_num = len(peptide_files_list)

        print("Saving all the peptide sequences to peptides.list file ...")
        peptide_filesList_list = chunkIt(peptide_files_list, scoop.SIZE)
        parts_list = [str(part) for part in range(1, scoop.SIZE + 1)]
        try:
            shared.setConst(MAX_PEPTIDE_LENGTH=max_peptide_length)
        except TypeError:
            pass
        try:
            shared.setConst(PROGBAR=progbar)
        except TypeError:
            pass
        shared.setConst(NEW_ALL_CHAINSCORE_LIST=all_chainScore_list)
        results = list(futures.map(Peptide.write_peptides_from_pickle_file,
                                   peptide_filesList_list,
                                   parts_list,
                                   [workdir]*len(parts_list)))
        # results will be an empty list
        # Now concatenate all parts of peptides.list and peptides.fasta
        fasta_filehandler = open(workdir + "/peptides." + str(max_peptide_length) + "mers.fasta", 'w')
        list_filehandler = open(workdir + "/peptides." + str(max_peptide_length) + "mers.list", 'w')
        index = 0
        for part in parts_list:
            with open(workdir + "/peptides." + str(max_peptide_length) + "mers.list.chunk" + part, 'r') as fin:
                for line in fin:
                    mo = re.search("^(.* peptide )P_[0-9]+_[0-9]+ = ([A-Z]+)( .*)$", line)
                    if mo:
                        line_part1 = mo.group(1)
                        peptide_seq = mo.group(2)
                        line_part2 = mo.group(3)
                        peptide_name = "P_" + str(len(peptide_seq)) + "_" + str(index)
                        # list_filehandler.write("Chain probability = "+str(chain[-1])+" peptide P_"+str(len(peptide))+"_"+str(index+1) + " = "+''.join(peptide)+" "+'-'.join(chain[0:-1])+" "+','.join(map(str, peptideProbList[:-1]))+"\n")
                        # fasta_filehandler.write(">P_"+str(len(peptide))+"_"+str(index+1)+"\n"+''.join(peptide)+"\n")
                        list_filehandler.write(line_part1 + peptide_name + " = " + peptide_seq + line_part2 + "\n")
                        fasta_filehandler.write(">" + peptide_name + "\n" + peptide_seq + "\n")
                        index += 1
            os.remove("peptides." + str(max_peptide_length) + "mers.list.chunk" + part)
        fasta_filehandler.close()
        list_filehandler.close()
        shutil.rmtree(workdir + '/tmp_peptide_folder/', ignore_errors=True)