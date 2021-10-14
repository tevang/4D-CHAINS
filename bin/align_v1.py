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

#from __future__ import print_function
#from __future__ import unicode_literals
import bz2
import multiprocessing
import pickle
import re
import sys, os
# Import 4D-CHAINS libraries
import traceback
from argparse import ArgumentParser, RawDescriptionHelpFormatter

from lib.alignment import create_consensus_alignemnt, write_consensus_alignment, create_overlapped_chain_alignments, \
    save_overlappingChainsGroups, write_overlapped_chains_alignment, write_consensus_alignment_and_absolute_macthes, \
    Alignment, filter_alignments_using_absolute_matches
from lib.connectivities import Connectivities, print_results_summary, shutil
from lib.global_func import ColorPrint, ProgressBar, shared, split_peptidefile, futures, align_and_read, chunkIt

CHAINS_BIN_DIR = os.path.dirname(os.path.realpath(__file__))
sys.path.append(CHAINS_BIN_DIR)


## Parse command line arguments
def cmdlineparse():
    parser = ArgumentParser(description="command line arguments", formatter_class=RawDescriptionHelpFormatter,
                            epilog="""
EXAMPLE: python -m scoop -n 8 /usr/local/bin/align_v1.py \
-fasta aLP.fasta \
-pseq peptides.5mers.fasta \
-plist peptides.5mers.list \
-rstfile results_summary.6mers_round1""")
    parser.add_argument("-fasta", dest="template_sequence_file", required=True,
                        help="teplate sequence file in fasta format",
                        metavar="<template sequence file>")
    parser.add_argument("-pseq", dest="peptide_fasta_file", required=True, default="protein.fasta",
                        help="peptide sequences file in fasta format",
                        metavar="<peptide sequences file>")
    parser.add_argument("-plist", dest="peptide_list_file", required=True, default="protein.list",
                        help="peptide list file", metavar="<peptide list file>")
    parser.add_argument("-rstfile", dest="ABSOLUTE_MATCHES_FILE", required=False,
                        help="file with absolute matches from previous run to be used as restraints",
                        metavar="<absolute matches restraint file>")
    parser.add_argument("-confile", dest="CONNECTIVITIES_FILE", required=False,
                        help="OPTIONAL: File with connectivities to annotate the resulting table with the NH assignments. It can have the full list of \
                             connectivities or just the pool of connectivities used to generate the peptides.",
                        metavar="<connectivities file>")
    parser.add_argument("-procaln", dest="PROCESSED_ALIGNMENT_FILE", required=False,
                        help="a processed alignment pickle file from a previous run that either failed or was stopped",
                        metavar="<processed alignment pickle file>")
    parser.add_argument("-matrix", dest="MATRIX", required=False, type=str, default="EBLOSUM90",
                        help="matrix file for alignment (default EBLOSUM90); allowed values: EBLOSUM30, EBLOSUM40, EBLOSUM50, EBLOSUM60, \
                             EBLOSUM62-12, EBLOSUM70, EBLOSUM80, EBLOSUM90, EBLOSUM35, EBLOSUM45, EBLOSUM55, EBLOSUM62, EBLOSUM65, EBLOSUM75, \
                             EBLOSUM85, EBLOSUMN. (default %(default)s)",
                        metavar="<matrix file>")
    parser.add_argument("-gapopen", dest="gap_open_penalty", required=False, type=float, default=100.0,
                        help="gap open penalty. (default %(default)s)",
                        metavar="<gap open penalty>")
    parser.add_argument("-gapextend", dest="gap_extension_penalty", required=False, type=float, default=10.0,
                        help="gap entension penalty. (default %(default)s)",
                        metavar="<gap extension penalty>")
    parser.add_argument("-outaln", dest="consensus_alignment_file", required=False, type=str, default="consensus_alignment",
                        help="output consensus alignment file. (default %(default)s)",
                        metavar="<consensus alignment file>")
    parser.add_argument("-cpl", dest="chars_per_line", required=False, type=int, default=162,
                        help="characters per line in the consensus alignment file. (default %(default)s)",
                        metavar="<characters per line>")
    parser.add_argument("-noalign", dest="DO_ALIGN", required=False, action='store_false',
                        help="align the template sequence with the peptides. (default %(default)s)")
    parser.add_argument("-useall", dest="USE_ALL_PEPTIDES_FOR_PREDICTION", required=False, action='store_true',
                        help="use all peptides (not only overlapping which is the default) for prediction score calculation. "
                             "(default %(default)s)")
    parser.add_argument("-inaln", dest="peptide_alignment_file", required=False, type=str,
                        help="input peptide alignment file from previous run. (default %(default)s)",
                        metavar="<peptide alignment file>")
    parser.add_argument("-idcutoff", dest="IDENTITY_CUTOFF", required=False, type=float, default=1.0,
                        help="sequence identity cutoff. (default %(default)s)",
                        metavar="<sequence identity cutoff>")
    parser.add_argument("-minoverlap", dest="OVERLAP_CUTOFF", required=False, type=int,
                        help="minimum number of overlapping TOCSY amino acid indices (default: maximum peptide length -1)",
                        metavar="<overlap cutoff>")
    parser.add_argument("-cpus", dest="CPUs", required=False, type=int, default=multiprocessing.cpu_count(),
                        help="number of processes for multi-threading calculations. (default %(default)s)",
                        metavar="<number of processes>")
    parser.add_argument("-compress", dest="COMPRESS", required=False, action='store_true', default=False,
                        help="compress intermediate files to save disk space with the cost of computation time. (default %(default)s)")
    parser.add_argument("-trim", dest="TRIM_ENDS", required=False, action='store_true', default=False,
                        help="trim the N- or C-term overhanging ends given the condition that there is at least one other contig which \
                             extends at least 2 positions in the respective direction and that it conflicts with the N- or C-term of the \
                             currect contig. (default %(default)s)")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s version 1.2')
    args=parser.parse_args()
    return args




if __name__ == "__main__":
    
    try:
        
        args = cmdlineparse()
        #################################################### SANITY CHECKS ####################################################
        #if args.H_weight <= 0 or args.H_weight > 1:
        #    print "ERROR: -wH must be a number greater than 0 and lower or equal to 1!"
        #elif args.C_weight <= 0 or args.C_weight > 1:
        #    print "ERROR: -wC must be a number greater than 0 and lower or equal to 1!"
        for fname in [args.template_sequence_file, args.peptide_fasta_file, args.peptide_list_file]:
            if len([l for l in open(fname, 'r')]) == 0:
                raise Exception(ColorPrint("ERROR: Input file %s is empty! Decrease the -minlen." % fname, "FAIL"))
        # INITIALIZE args.OVERLAP_CUTOFF
        if not args.OVERLAP_CUTOFF:
            f = open(args.peptide_list_file, 'r')
            for line in f:
                words = line.split()
                if len(words) > 0:
                    max_peptide_length = len(words[7])  # ATTENTION: I assume that all peptides have the same length.
                break
            f.close()
            args.OVERLAP_CUTOFF = max_peptide_length -1
        #######################################################################################################################
        
        if not args.PROCESSED_ALIGNMENT_FILE:
            
            # ALIGN AND LOAD ALIGNMENT IN PARALLEL
            progbar = ProgressBar(100)
            aln_index = 0
            shared.setConst(PROGBAR=progbar)
            shared.setConst(ALN_INDEX=aln_index)
            shared.setConst(ARGS=args)
            shared.setConst(chains_bin_dir=CHAINS_BIN_DIR)
            fasta_file_list = [args.peptide_fasta_file+".chunk"+str(num) for num in range(1, args.CPUs + 1)]
            if args.COMPRESS:
                list_file_list = [args.peptide_list_file+".chunk"+str(num)+".bz2" for num in range(1, args.CPUs + 1)]
            else:
                list_file_list = [args.peptide_list_file+".chunk"+str(num) for num in range(1, args.CPUs + 1)]
            chunkNo_list = list(range(1, args.CPUs + 1))
            
            if args.DO_ALIGN == True:
                split_peptidefile(args)
            else:
                split = False
                for fname in fasta_file_list:
                    if os.path.isfile(fname) == False:
                        split = True
                        break
                split_peptidefile(args)
                
            # alignment_results in a list of tuples of the form (peptide_identity_dict, peptide_score_dict, peptideName2chain_dict, peptide_alignment_dict, peptideName_chainScoreList_dict, max_peptide_length)
            alignment_results = list(futures.map(align_and_read, fasta_file_list, list_file_list, chunkNo_list))
            peptide_identity_dict, peptide_score_dict, peptideName2chain_dict, peptide_alignment_dict, peptideName_chainScoreList_dict, max_peptide_length = {}, {}, {}, {}, {}, 0
            for worker_results in alignment_results:
                peptide_identity_dict.update(worker_results[0])
                peptide_score_dict.update(worker_results[1])
                peptideName2chain_dict.update(worker_results[2])
                peptide_alignment_dict.update(worker_results[3])
                peptideName_chainScoreList_dict.update(worker_results[4])
                if max_peptide_length < worker_results[5]:
                    max_peptide_length = worker_results[5]
            #print "DEBUG: len(peptide_identity_dict)", len(peptide_identity_dict)
            #print "DEBUG: len(peptide_score_dict)", len(peptide_score_dict)
            #print "DEBUG: len(peptideName2chain_dict)", len(peptideName2chain_dict)
            #print "DEBUG: len(peptide_alignment_dict)", len(peptide_alignment_dict)
            #print "DEBUG: len(peptideName_chainScoreList_dict)", len(peptideName_chainScoreList_dict)
            #print "DEBUG max_peptide_length=",max_peptide_length
            
            ## IF ALL WORKERS HAVE FINISHED SUCCESSFULLY, BACKUP THE DATA
            if args.COMPRESS:
                f = bz2.BZ2File('processed_alignment.pickle.bz2', 'wb')
            else:
                f = open('processed_alignment.pickle', 'wb')
            pickle.dump((peptide_identity_dict, peptide_score_dict, peptideName2chain_dict, peptide_alignment_dict, peptideName_chainScoreList_dict, max_peptide_length), f)
            f.close()
        else:   ## IF A PROCESSED ALIGNMENT PICKLE FILE WAS PROVIDED, LOAD IT AND POPULATE ALL RELEVANT DICTIONARIES
            if args.COMPRESS:
                f = bz2.BZ2File(args.PROCESSED_ALIGNMENT_FILE, "rb")
            else:
                f = open(args.PROCESSED_ALIGNMENT_FILE, "rb")
            pickle_object = pickle.load(f)
            f.close()
            peptide_identity_dict = pickle_object[0]
            peptide_score_dict = pickle_object[1]
            peptideName2chain_dict = pickle_object[2]
            peptide_alignment_dict = pickle_object[3]
            peptideName_chainScoreList_dict = pickle_object[4]
            max_peptide_length = pickle_object[5]
            
            
        ordered_consensus_peptide_alignment_list, \
        ordered_consensus_chain_alignment_list, \
        ordered_consensus_TAAIGScore_alignment_list, \
        template_sequence_alignment_string = create_consensus_alignemnt(peptide_identity_dict,
                                                                        peptide_score_dict,
                                                                        peptideName2chain_dict,
                                                                        peptide_alignment_dict,
                                                                        peptideName_chainScoreList_dict,
                                                                        args)
        
        write_consensus_alignment(ordered_consensus_peptide_alignment_list,
                                  ordered_consensus_chain_alignment_list,
                                  ordered_consensus_TAAIGScore_alignment_list,
                                  template_sequence_alignment_string,
                                  args,
                                  1)
        
        refchain_overlappingChainsList_dict, \
        all_overlapped_chain_alignments_list = \
            create_overlapped_chain_alignments(ordered_consensus_peptide_alignment_list,
                                               ordered_consensus_chain_alignment_list,
                                               args)
        # refchain_overlappingChainsList_dict:      ordereddict with keys the reference chain alignment that have at least 1 overlapping chain alignment, and values
        #                                           tuples of the overlapping alignments
        
        ####################################################################################################################################################
        # split all_overlapped_chain_alignments_list into CPUs parts
        overlapped_chain_alignmentsList_list = chunkIt(all_overlapped_chain_alignments_list, args.CPUs)
        progbar = ProgressBar(100)
        try:
            shared.setConst(PROGBAR=progbar)
        except TypeError:
            pass
        shared.setConst(ALL_OVERLAPPED_CHAIN_ALIGNMENTS_LIST_1=all_overlapped_chain_alignments_list)
        overlappingChainsGroup_num = len(all_overlapped_chain_alignments_list)
        tuple_indices_list = list(range(1, overlappingChainsGroup_num +1))
        iterationNo_list = ["1" for i in range(1, overlappingChainsGroup_num +1)]
        tuple_indicesList_list = chunkIt(tuple_indices_list, args.CPUs)
        ## REMOVE CONTIGS OF OVERLAPPING CHAIN ALIGNMENTS THAT ARE INCLUDED IN BIGGER CONTIGS
        print("Removing redundant overlapping chain alignments...")
        # Parallel execution
        results = list(futures.map(save_overlappingChainsGroups,
                                   overlapped_chain_alignmentsList_list,
                                   tuple_indicesList_list,
                                   iterationNo_list))
        # # Serial execution for debugging
        # results = []
        # for o, t, i in zip(overlapped_chain_alignmentsList_list, tuple_indicesList_list, iterationNo_list):
        #     results.append(save_overlappingChainsGroups(o, t, i))
        ####################################################################################################################################################
          
        write_overlapped_chains_alignment(template_sequence_alignment_string, args, 1)
        
        absolute_matches_alignment, \
        position_TAAIG_confidenceScore_mdict = write_consensus_alignment_and_absolute_macthes(template_sequence_alignment_string,
                                                                                                   refchain_overlappingChainsList_dict,
                                                                                                   max_peptide_length,
                                                                                                   ordered_consensus_peptide_alignment_list,
                                                                                                   ordered_consensus_chain_alignment_list,
                                                                                                   ordered_consensus_TAAIGScore_alignment_list,
                                                                                                   args,
                                                                                                   None,
                                                                                                   1)
        
        ### 2nd ITERATION
        if args.ABSOLUTE_MATCHES_FILE:  # USE THE RESTRAINTS, IF PROVIDED, ONLY BEFORE THE 2nd ITERATION
            absolute_matches_alignment = None   # clean absolute_matches_alignment and reinitialize it from input file
            absolute_matches_alignment = Alignment.read_NHmap_file(args.ABSOLUTE_MATCHES_FILE)
        #print "DEBUG: iteration1 absolute_matches_alignment=", absolute_matches_alignment
            
        iter2_ordered_consensus_peptide_alignment_list, \
        iter2_ordered_consensus_chain_alignment_list, \
        iter2_ordered_consensus_TAAIGScore_alignment_list, \
        iter2tmp_position_TAAIG_confidenceScore_mdict, \
        consensus_absolute_matches_alignment = filter_alignments_using_absolute_matches(absolute_matches_alignment,
                                                                                        ordered_consensus_peptide_alignment_list,
                                                                                        ordered_consensus_chain_alignment_list,
                                                                                        ordered_consensus_TAAIGScore_alignment_list,
                                                                                        position_TAAIG_confidenceScore_mdict)
        
        del ordered_consensus_peptide_alignment_list, ordered_consensus_chain_alignment_list, ordered_consensus_TAAIGScore_alignment_list
        #print "DEBUG: iter2_ordered_consensus_chain_alignment_list=", iter2_ordered_consensus_chain_alignment_list
        write_consensus_alignment(iter2_ordered_consensus_peptide_alignment_list, iter2_ordered_consensus_chain_alignment_list, iter2_ordered_consensus_TAAIGScore_alignment_list,
                                  template_sequence_alignment_string, args, 2)
            
        iter2_refchain_overlappingChainsList_dict, iter2_all_overlapped_chain_alignments_list = create_overlapped_chain_alignments(iter2_ordered_consensus_peptide_alignment_list, iter2_ordered_consensus_chain_alignment_list, args)
        del refchain_overlappingChainsList_dict, all_overlapped_chain_alignments_list
        
        ####################################################################################################################################################
        # split iter2_all_overlapped_chain_alignments_list into CPUs parts
        overlapped_chain_alignmentsList_list = chunkIt(iter2_all_overlapped_chain_alignments_list, args.CPUs)
        shared.setConst(ALL_OVERLAPPED_CHAIN_ALIGNMENTS_LIST_2=iter2_all_overlapped_chain_alignments_list)
        overlappingChainsGroup_num = len(iter2_all_overlapped_chain_alignments_list)
        tuple_indices_list = list(range(1, overlappingChainsGroup_num +1))
        iterationNo_list = ["2" for i in range(1, overlappingChainsGroup_num +1)]
        tuple_indicesList_list = chunkIt(tuple_indices_list, args.CPUs)
        results = list(futures.map(save_overlappingChainsGroups, overlapped_chain_alignmentsList_list, tuple_indicesList_list, iterationNo_list))
        ####################################################################################################################################################
        
        write_overlapped_chains_alignment(template_sequence_alignment_string, args, 2)
        
        iter2_absolute_matches_alignment, iter2_position_TAAIG_confidenceScore_mdict = write_consensus_alignment_and_absolute_macthes(template_sequence_alignment_string,
                iter2_refchain_overlappingChainsList_dict, max_peptide_length, iter2_ordered_consensus_peptide_alignment_list, iter2_ordered_consensus_chain_alignment_list,
                iter2_ordered_consensus_TAAIGScore_alignment_list, args, iter2tmp_position_TAAIG_confidenceScore_mdict, 2)
        #print "DEBUG: iter2_absolute_matches_alignment=", iter2_absolute_matches_alignment
        
        
        ### 3rd ITERATION, no update of confidence scores from the previous iterations
        #iter3_ordered_consensus_peptide_alignment_list, iter3_ordered_consensus_chain_alignment_list, iter3_ordered_consensus_TAAIGScore_alignment_list, iter3tmp_position_TAAIG_confidenceScore_mdict, consensus_absolute_matches_alignment = filter_alignments_using_absolute_matches(iter2_absolute_matches_alignment,
        #                    iter2_ordered_consensus_peptide_alignment_list, iter2_ordered_consensus_chain_alignment_list, iter2_ordered_consensus_TAAIGScore_alignment_list, iter2_position_TAAIG_confidenceScore_mdict)
        #
        #del iter2_ordered_consensus_peptide_alignment_list, iter2_ordered_consensus_chain_alignment_list, iter2_ordered_consensus_TAAIGScore_alignment_list
        ##print "DEBUG: iter3_ordered_consensus_chain_alignment_list=", iter3_ordered_consensus_chain_alignment_list
        #write_consensus_alignment(iter3_ordered_consensus_peptide_alignment_list, iter3_ordered_consensus_chain_alignment_list, iter3_ordered_consensus_TAAIGScore_alignment_list,
        #                          template_sequence_alignment_string, 3)
        #    
        #iter3_unique_overlapped_chain_alignments_set, iter3_refchain_overlappingChainsList_dict = create_overlapped_chain_alignments(iter3_ordered_consensus_chain_alignment_list)
        ##del unique_overlapped_chain_alignments_set, refchain_overlappingChainsList_dict 
        #
        #write_overlapped_chains_alignment(template_sequence_alignment_string, iter3_unique_overlapped_chain_alignments_set, 3)
        #
        #iter3_absolute_matches_alignment, iter3_position_TAAIG_confidenceScore_mdict = write_consensus_alignment_and_absolute_macthes(template_sequence_alignment_string, iter3_unique_overlapped_chain_alignments_set,
        #        iter3_refchain_overlappingChainsList_dict, max_peptide_length, iter3_ordered_consensus_peptide_alignment_list, iter3_ordered_consensus_chain_alignment_list,
        #        iter3_ordered_consensus_TAAIGScore_alignment_list, iter3tmp_position_TAAIG_confidenceScore_mdict, 3)
        ##print "DEBUG: iter3_absolute_matches_alignment=", iter3_absolute_matches_alignment
        
        
        ## 3rd ITERATION, update of confidence scores from the previous iterations
        iter3_ordered_consensus_peptide_alignment_list, iter3_ordered_consensus_chain_alignment_list, iter3_ordered_consensus_TAAIGScore_alignment_list, iter3tmp_position_TAAIG_confidenceScore_mdict, consensus_absolute_matches_alignment = filter_alignments_using_absolute_matches(iter2_absolute_matches_alignment,
                    iter2_ordered_consensus_peptide_alignment_list, iter2_ordered_consensus_chain_alignment_list, iter2_ordered_consensus_TAAIGScore_alignment_list, iter2_position_TAAIG_confidenceScore_mdict)
        
        del iter2_ordered_consensus_peptide_alignment_list, iter2_ordered_consensus_chain_alignment_list, iter2_ordered_consensus_TAAIGScore_alignment_list
        #print "DEBUG: iter2_ordered_consensus_chain_alignment_list=", iter2_ordered_consensus_chain_alignment_list
        write_consensus_alignment(iter3_ordered_consensus_peptide_alignment_list, iter3_ordered_consensus_chain_alignment_list, iter3_ordered_consensus_TAAIGScore_alignment_list,
                                  template_sequence_alignment_string, args, 3)
        
        iter3_refchain_overlappingChainsList_dict, iter3_all_overlapped_chain_alignments_list = create_overlapped_chain_alignments(iter3_ordered_consensus_peptide_alignment_list, iter3_ordered_consensus_chain_alignment_list, args)
        del iter2_refchain_overlappingChainsList_dict, iter2_all_overlapped_chain_alignments_list 
        
        ####################################################################################################################################################
        # split iter3_all_overlapped_chain_alignments_list into CPUs parts
        overlapped_chain_alignmentsList_list = chunkIt(iter3_all_overlapped_chain_alignments_list, args.CPUs)
        shared.setConst(ALL_OVERLAPPED_CHAIN_ALIGNMENTS_LIST_3=iter3_all_overlapped_chain_alignments_list)
        overlappingChainsGroup_num = len(iter3_all_overlapped_chain_alignments_list)
        tuple_indices_list = list(range(1, overlappingChainsGroup_num +1))
        iterationNo_list = ["3" for i in range(1, overlappingChainsGroup_num +1)]
        tuple_indicesList_list = chunkIt(tuple_indices_list, args.CPUs)
        results = list(futures.map(save_overlappingChainsGroups, overlapped_chain_alignmentsList_list, tuple_indicesList_list, iterationNo_list))
        ####################################################################################################################################################
        
        write_overlapped_chains_alignment(template_sequence_alignment_string, args, 3)
        
        iter3_absolute_matches_alignment, iter3_position_TAAIG_confidenceScore_mdict = write_consensus_alignment_and_absolute_macthes(template_sequence_alignment_string,
                iter3_refchain_overlappingChainsList_dict, max_peptide_length, iter3_ordered_consensus_peptide_alignment_list, iter3_ordered_consensus_chain_alignment_list,
                iter3_ordered_consensus_TAAIGScore_alignment_list, args, iter3tmp_position_TAAIG_confidenceScore_mdict, 3)
        #print "DEBUG: iter3_absolute_matches_alignment=", iter3_absolute_matches_alignment
    
    
        #### 4th ITERATION
        iter4_ordered_consensus_peptide_alignment_list, iter4_ordered_consensus_chain_alignment_list, iter4_ordered_consensus_TAAIGScore_alignment_list, iter4tmp_position_TAAIG_confidenceScore_mdict, consensus_absolute_matches_alignment = filter_alignments_using_absolute_matches(iter3_absolute_matches_alignment,
                    iter3_ordered_consensus_peptide_alignment_list, iter3_ordered_consensus_chain_alignment_list, iter3_ordered_consensus_TAAIGScore_alignment_list, iter3_position_TAAIG_confidenceScore_mdict)
        
        del iter3_ordered_consensus_peptide_alignment_list, iter3_ordered_consensus_chain_alignment_list, iter3_ordered_consensus_TAAIGScore_alignment_list
        #print "DEBUG: iter3_ordered_consensus_chain_alignment_list=", iter3_ordered_consensus_chain_alignment_list
        write_consensus_alignment(iter4_ordered_consensus_peptide_alignment_list, iter4_ordered_consensus_chain_alignment_list, iter4_ordered_consensus_TAAIGScore_alignment_list,
                                  template_sequence_alignment_string, args, 4)
        
        iter4_refchain_overlappingChainsList_dict, iter4_all_overlapped_chain_alignments_list = create_overlapped_chain_alignments(iter4_ordered_consensus_peptide_alignment_list, iter4_ordered_consensus_chain_alignment_list, args)
        del  iter3_refchain_overlappingChainsList_dict, iter3_all_overlapped_chain_alignments_list 
        
        ####################################################################################################################################################
        # split all_overlapped_chain_alignments_list into CPUs parts
        overlapped_chain_alignmentsList_list = chunkIt(iter4_all_overlapped_chain_alignments_list, args.CPUs)
        shared.setConst(ALL_OVERLAPPED_CHAIN_ALIGNMENTS_LIST_4=iter4_all_overlapped_chain_alignments_list)
        overlappingChainsGroup_num = len(iter4_all_overlapped_chain_alignments_list)
        tuple_indices_list = list(range(1, overlappingChainsGroup_num +1))
        iterationNo_list = ["4" for i in range(1, overlappingChainsGroup_num +1)]
        tuple_indicesList_list = chunkIt(tuple_indices_list, args.CPUs)
        results = list(futures.map(save_overlappingChainsGroups, overlapped_chain_alignmentsList_list, tuple_indicesList_list, iterationNo_list))
        ####################################################################################################################################################
        
        write_overlapped_chains_alignment(template_sequence_alignment_string, args, 4)
        
        iter4_absolute_matches_alignment, iter4_position_TAAIG_confidenceScore_mdict = write_consensus_alignment_and_absolute_macthes(template_sequence_alignment_string,
                iter4_refchain_overlappingChainsList_dict, max_peptide_length, iter4_ordered_consensus_peptide_alignment_list, iter4_ordered_consensus_chain_alignment_list,
                iter4_ordered_consensus_TAAIGScore_alignment_list, args, iter4tmp_position_TAAIG_confidenceScore_mdict, 4)
        #print "DEBUG: iter4_absolute_matches_alignment=", iter4_absolute_matches_alignment
        
        
        ## 5th ITERATION
        iter5_ordered_consensus_peptide_alignment_list, iter5_ordered_consensus_chain_alignment_list, iter5_ordered_consensus_TAAIGScore_alignment_list, iter5tmp_position_TAAIG_confidenceScore_mdict, consensus_absolute_matches_alignment = filter_alignments_using_absolute_matches(iter4_absolute_matches_alignment,
                    iter4_ordered_consensus_peptide_alignment_list, iter4_ordered_consensus_chain_alignment_list, iter4_ordered_consensus_TAAIGScore_alignment_list, iter4_position_TAAIG_confidenceScore_mdict)
        
        del iter4_ordered_consensus_peptide_alignment_list, iter4_ordered_consensus_chain_alignment_list, iter4_ordered_consensus_TAAIGScore_alignment_list
        #print "DEBUG: iter4_ordered_consensus_chain_alignment_list=", iter4_ordered_consensus_chain_alignment_list
        write_consensus_alignment(iter5_ordered_consensus_peptide_alignment_list, iter5_ordered_consensus_chain_alignment_list, iter5_ordered_consensus_TAAIGScore_alignment_list,
                                  template_sequence_alignment_string, args, 5)
        
        iter5_refchain_overlappingChainsList_dict, iter5_all_overlapped_chain_alignments_list = create_overlapped_chain_alignments(iter5_ordered_consensus_peptide_alignment_list, iter5_ordered_consensus_chain_alignment_list, args)
        del iter4_refchain_overlappingChainsList_dict, iter4_all_overlapped_chain_alignments_list 
        
        ####################################################################################################################################################
        # split iter5_all_overlapped_chain_alignments_list into CPUs parts
        overlapped_chain_alignmentsList_list = chunkIt(iter5_all_overlapped_chain_alignments_list, args.CPUs)
        shared.setConst(ALL_OVERLAPPED_CHAIN_ALIGNMENTS_LIST_5=iter5_all_overlapped_chain_alignments_list)
        overlappingChainsGroup_num = len(iter5_all_overlapped_chain_alignments_list)
        tuple_indices_list = list(range(1, overlappingChainsGroup_num +1))
        iterationNo_list = ["5" for i in range(1, overlappingChainsGroup_num +1)]
        tuple_indicesList_list = chunkIt(tuple_indices_list, args.CPUs)
        results = list(futures.map(save_overlappingChainsGroups, overlapped_chain_alignmentsList_list, tuple_indicesList_list, iterationNo_list))
        ####################################################################################################################################################
        
        write_overlapped_chains_alignment(template_sequence_alignment_string, args, 5)
        
        iter5_absolute_matches_alignment, iter5_position_TAAIG_confidenceScore_mdict = write_consensus_alignment_and_absolute_macthes(template_sequence_alignment_string,
                iter5_refchain_overlappingChainsList_dict, max_peptide_length, iter5_ordered_consensus_peptide_alignment_list, iter5_ordered_consensus_chain_alignment_list,
                iter5_ordered_consensus_TAAIGScore_alignment_list, args, iter5tmp_position_TAAIG_confidenceScore_mdict, 5)
        #print "DEBUG: iter5_absolute_matches_alignment=", iter5_absolute_matches_alignment
        
        ## 6th ITERATION
        iter6_ordered_consensus_peptide_alignment_list, iter6_ordered_consensus_chain_alignment_list, iter6_ordered_consensus_TAAIGScore_alignment_list, iter6tmp_position_TAAIG_confidenceScore_mdict, consensus_absolute_matches_alignment = filter_alignments_using_absolute_matches(iter5_absolute_matches_alignment,
                    iter5_ordered_consensus_peptide_alignment_list, iter5_ordered_consensus_chain_alignment_list, iter5_ordered_consensus_TAAIGScore_alignment_list, iter5_position_TAAIG_confidenceScore_mdict)
        
        del iter5_ordered_consensus_peptide_alignment_list, iter5_ordered_consensus_chain_alignment_list, iter5_ordered_consensus_TAAIGScore_alignment_list
        #print "DEBUG: iter5_ordered_consensus_chain_alignment_list=", iter5_ordered_consensus_chain_alignment_list
        write_consensus_alignment(iter6_ordered_consensus_peptide_alignment_list, iter6_ordered_consensus_chain_alignment_list, iter6_ordered_consensus_TAAIGScore_alignment_list,
                                  template_sequence_alignment_string, args, 6)
        
        iter6_refchain_overlappingChainsList_dict, iter6_all_overlapped_chain_alignments_list = create_overlapped_chain_alignments(iter6_ordered_consensus_peptide_alignment_list, iter6_ordered_consensus_chain_alignment_list, args)
        del iter5_refchain_overlappingChainsList_dict, iter5_all_overlapped_chain_alignments_list 
        
        ####################################################################################################################################################
        # split iter6_all_overlapped_chain_alignments_list into CPUs parts
        overlapped_chain_alignmentsList_list = chunkIt(iter6_all_overlapped_chain_alignments_list, args.CPUs)
        shared.setConst(ALL_OVERLAPPED_CHAIN_ALIGNMENTS_LIST_6=iter6_all_overlapped_chain_alignments_list)
        overlappingChainsGroup_num = len(iter6_all_overlapped_chain_alignments_list)
        tuple_indices_list = list(range(1, overlappingChainsGroup_num +1))
        iterationNo_list = ["6" for i in range(1, overlappingChainsGroup_num +1)]
        tuple_indicesList_list = chunkIt(tuple_indices_list, args.CPUs)
        results = list(futures.map(save_overlappingChainsGroups, overlapped_chain_alignmentsList_list, tuple_indicesList_list, iterationNo_list))
        ####################################################################################################################################################
        
        write_overlapped_chains_alignment(template_sequence_alignment_string, args, 6)
        
        iter6_absolute_matches_alignment, iter6_position_TAAIG_confidenceScore_mdict = write_consensus_alignment_and_absolute_macthes(template_sequence_alignment_string,
                iter6_refchain_overlappingChainsList_dict, max_peptide_length, iter6_ordered_consensus_peptide_alignment_list, iter6_ordered_consensus_chain_alignment_list,
                iter6_ordered_consensus_TAAIGScore_alignment_list, args, iter6tmp_position_TAAIG_confidenceScore_mdict, 6)
        #print "DEBUG: iter6_absolute_matches_alignment=", iter6_absolute_matches_alignment
        
        ## 7th ITERATION
        iter7_ordered_consensus_peptide_alignment_list, iter7_ordered_consensus_chain_alignment_list, iter7_ordered_consensus_TAAIGScore_alignment_list, iter7tmp_position_TAAIG_confidenceScore_mdict, consensus_absolute_matches_alignment = filter_alignments_using_absolute_matches(iter6_absolute_matches_alignment,
                    iter6_ordered_consensus_peptide_alignment_list, iter6_ordered_consensus_chain_alignment_list, iter6_ordered_consensus_TAAIGScore_alignment_list, iter6_position_TAAIG_confidenceScore_mdict)
        
        del iter6_ordered_consensus_peptide_alignment_list, iter6_ordered_consensus_chain_alignment_list, iter6_ordered_consensus_TAAIGScore_alignment_list
        #print "DEBUG: iter6_ordered_consensus_chain_alignment_list=", iter6_ordered_consensus_chain_alignment_list
        write_consensus_alignment(iter7_ordered_consensus_peptide_alignment_list, iter7_ordered_consensus_chain_alignment_list, iter7_ordered_consensus_TAAIGScore_alignment_list,
                                  template_sequence_alignment_string, args, 7)
        
        iter7_refchain_overlappingChainsList_dict, iter7_all_overlapped_chain_alignments_list = create_overlapped_chain_alignments(iter7_ordered_consensus_peptide_alignment_list, iter7_ordered_consensus_chain_alignment_list, args)
        del iter6_refchain_overlappingChainsList_dict, iter6_all_overlapped_chain_alignments_list 
        
        ####################################################################################################################################################
        # split iter7_all_overlapped_chain_alignments_list into CPUs parts
        overlapped_chain_alignmentsList_list = chunkIt(iter7_all_overlapped_chain_alignments_list, args.CPUs)
        shared.setConst(ALL_OVERLAPPED_CHAIN_ALIGNMENTS_LIST_7=iter7_all_overlapped_chain_alignments_list)
        overlappingChainsGroup_num = len(iter7_all_overlapped_chain_alignments_list)
        tuple_indices_list = list(range(1, overlappingChainsGroup_num +1))
        iterationNo_list = ["7" for i in range(1, overlappingChainsGroup_num +1)]
        tuple_indicesList_list = chunkIt(tuple_indices_list, args.CPUs)
        results = list(futures.map(save_overlappingChainsGroups, overlapped_chain_alignmentsList_list, tuple_indicesList_list, iterationNo_list))
        ####################################################################################################################################################
        
        write_overlapped_chains_alignment(template_sequence_alignment_string, args, 7)
        
        iter7_absolute_matches_alignment, iter7_position_TAAIG_confidenceScore_mdict = write_consensus_alignment_and_absolute_macthes(template_sequence_alignment_string,
                iter7_refchain_overlappingChainsList_dict, max_peptide_length, iter7_ordered_consensus_peptide_alignment_list, iter7_ordered_consensus_chain_alignment_list,
                iter7_ordered_consensus_TAAIGScore_alignment_list, args, iter7tmp_position_TAAIG_confidenceScore_mdict, iteration=7)
        #print "DEBUG: iter7_absolute_matches_alignment=", iter7_absolute_matches_alignment
        
        ## WRITE RESULTS SUMMARY (CONSENSUS CONFIDENCE SCORE AND ABSOLUTE MATCHES)
        final_ordered_consensus_peptide_alignment_list, final_ordered_consensus_chain_alignment_list, final_ordered_consensus_TAAIGScore_alignment_list, final_position_TAAIG_confidenceScore_mdict, consensus_absolute_matches_alignment = filter_alignments_using_absolute_matches(iter7_absolute_matches_alignment,
                    iter7_ordered_consensus_peptide_alignment_list, iter7_ordered_consensus_chain_alignment_list, iter7_ordered_consensus_TAAIGScore_alignment_list, iter7_position_TAAIG_confidenceScore_mdict)
        
        # LOAD CONNECTIVITES IF PROVIDED
        i_iminus1_dict = None   # default value
        if args.CONNECTIVITIES_FILE:
            i_iminus1_dict = Connectivities.load_connectivities_from_file(args.CONNECTIVITIES_FILE)

        print_results_summary(template_sequence_alignment_string,
                              consensus_absolute_matches_alignment,
                              final_position_TAAIG_confidenceScore_mdict,
                              final_ordered_consensus_chain_alignment_list,
                              i_iminus1_dict,
                              SPECTRUM_COMBO='TOCSY-HCNH')
        
        # remove intermediate files
        shutil.rmtree('tmp_overlapped_chain_alignments_folder/', ignore_errors=True)
        fnames = os.listdir("./")
        fpattern = re.compile("peptide_alignment.needle.chunk.*")
        files_list = list(filter(fpattern.search, fnames))
        for fname in files_list:
            os.remove(fname)
    
    except:
        type, value, tb = sys.exc_info()
        lines = traceback.format_exception(type, value, tb)
        print((''.join(lines)))
        raise
    