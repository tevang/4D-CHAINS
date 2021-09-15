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
"""
This is the new version of "form_peptides.py" script which work only with the 4D-NOESY spectrum.
TODO: to be updated.
At the end you should have two files at each protein folder:

(i) "%s_HCNH-HCNH_connectivities.list" % protein: a text file with the connectivities.

(ii) "%s_HCNH-HCNH_connectivities.pkl" % protein: a pickle file which contains the 'i_iminus1_complete_dict' and
'NOESY_NOESY_peak_connectivities_mdict'. 'i_iminus1_complete_dict' dictionary contains all possible connectivities of
every NOESY AAIG1; has keys the AAIG signatures and values lists of quintuplets (tuples) consisting of
(NOESY AAIG2 signature, # matching peaks, total # peaks, intersection, multi score).
'NOESY_NOESY_peak_connectivities_mdict' is a 2D dictionary of the form:
# AAIG1 -> matching AAIG2 -> [(matching NOESY peak1, matching NOESY peak 2, distance), ...], where peaks are Peak objects
from lib/peak.py.

# Always launch the script with 'scoop' for parallelization. E.g.
python -m scoop -n 8 ./find_NOESY_connectivities.py
"""

from lib.hcnh_connectivities.__init__ import *
from lib.hcnh_connectivities._write_methods import write_spectrum_info, write_NOESY_NOESY_connectivities
from lib.trees.chain import *
from lib.trees.peptide import *
from lib.bayes.statistics import *


##~~~~~~~~~~~~~~~~~~~~~~~~ GLOBAL PARAMETERS. Set the right paths on your computer ~~~~~~~~~~~~~~~~~~~~##
proofread=True          # Will print in red which correct connectivities were removed by MULTI_RATIO or
                        # REVERSE CONNECTIVITY criterion.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##


## Parse command line arguments
def cmdlineparse():
    parser = ArgumentParser(description="command line arguments", formatter_class=RawDescriptionHelpFormatter,
            epilog="""
    EXAMPLE: 
       python -m scoop -n 8 ./find_NOESY_connectivities.py -hsqc TDhsqc.list -tocsy TDtocsy.list -noesy 
       TDnoesy.list -tolH 0.01 -tolC 0.1 -rtolH 0.01 -rtolN 0.1 -mcutoff 1.0 -acutoff 1.0 
       -wH 0.5 -wC 1
                   """)
    parser.add_argument("-workdir", dest="WORKDIR", required=False, default=None,
                        help="the working directory",
                        metavar="<working directory>")
    parser.add_argument("-fasta", dest="FASTA_FILE", required=True,
                        help="teplate sequence file in fasta format", metavar="<template sequence file>")
    parser.add_argument("-hsqc", dest="HSQC_FILE", required=True,
                        help="2D N-H HSQC root spectrum", metavar="<2D N-H HSQC input file>")
    parser.add_argument("-hcnh", dest="HCNH_FILE", required=True,
                        help="4D NOESY (HCNOENH) file", metavar="<4D HCNH NOESY input file>")
    parser.add_argument("-hnnh", dest="HNNH_FILE", required=False,
                        help="4D NOESY (HNNOENH) file. It is not necessary but recommended.",
                        metavar="<4D HNNH NOESY input file>")
    parser.add_argument("-tolH", dest="tolH", required=False, type=float, default=0.04,
                        help="tolerance for the proton resonance when matching NOESY peaks to TOCSY peaks. (default %(default)s)",
                        metavar="<proton tolerance>")
    parser.add_argument("-tolC", dest="tolC", required=False, type=float, default=0.4,
                        help="tolerance for the carbon resonance when matching NOESY peaks to TOCSY peaks. (default %(default)s)",
                        metavar="<carbon tolerance>")
    # parser.add_argument("-rtolH", dest="rtolH", required=False, type=float, default=0.02,
    # help="tolerance for the proton resonance when matching peaks between the root spectum and TOCSY/NOESY", metavar="<proton tolerance>")
    # parser.add_argument("-rtolN", dest="rtolN", required=False, type=float, default=0.2,
    # help="tolerance for the nitrogen resonance when matching peaks between the root spectum and TOCSY/NOESY", metavar="<carbon tolerance>")
    parser.add_argument("-update", dest="DOWNLOAD_CS_HISTOGRAMS", required=False, action='store_true',
                        help="download the latest chemical shift histgrams from BMRB")
    parser.add_argument("-wH", dest="H_weight", required=False, type=float, default=0.1,
                        help="weight in (0.0-1.0] to apply on aa type prediction from aliphatic H resonances. (default %(default)s)",
                        metavar="<H weight>")
    parser.add_argument("-wC", dest="C_weight", required=False, type=float, default=1.0,
                        help="weight in (0.0, 1.0] to apply on aa type prediction from aliphatic C resonances. (default %(default)s)",
                        metavar="<C weight>")
    # parser.add_argument("-examples", dest='EXAMPLES', action='store_true', required=False, help="Usage and Examples")
    # parser.add_argument("-mcutoff", dest="RESONANCE_MATCH_CUTOFF", required=False, type=float, default=0.8,
    #                     help="number 0.0-1.0 saying how much of TOCSY resonances should match in the NOESY in order to consider "
    #                          "it a possible match. (default %(default)s)",
    #                     metavar="<resonance match cutoff>")   # NOT SURE IF APPLICABLE IN NOESY!
    parser.add_argument("-zmcutoff", dest="ZSCORE_RESONANCE_MATCH_CUTOFF", required=False, type=float, default=0.0,
                        help="a real number specifying the lower Z-score for a resonance match to be retained for chain building."
                             " (default %(default)s)",
                        metavar="<Z-score resonance match cutoff>")
    parser.add_argument("-mratio", dest="MULTI_RATIO", required=False, type=float, default=10.0,
                        help="keeps all connected AAIG2s that have 'multi' score up to MULTI_RATIO times lower than the maximum"
                            "'multi' score of the connectivities of AAIG1. The higher the MULTI_RATIO, the more connectivities"
                            "are retained.")
    parser.add_argument("-revconn", "--reverse-conn", dest="REVERSE_CONN", required=False, action='store_true', default=False,
                        help="EXPERIMENTAL FEATURE. If not both AAIG1-AAIG2 and AAIG2-AAIG1 connectivites exist, then remove any of them from the "
                             "connectivities list.When MULTI_RATIO is low it removes correct connectivities. Leave it as False "
                             "for the time being. Default %(default)s")
    # parser.add_argument("-acutoff", dest="ASSIGNMENT_CUTOFF", required=False, type=float, default=None,
    #                     help="number 0.0-1.0 saying how much of TOCSY resonances should match with a particular aa type in order to be considered, e.g. \
    #                     if TOCSY resonances are 5, 0.8 will mean that at least 4 resonances must match. (default %(default)s)",
    #                     metavar="<aa type assignment cutoff>")    # NOT SURE IF APPLICABLE IN NOESY!
    parser.add_argument("-zacutoff", dest="ZSCORE_ASSIGNMENT_CUTOFF", required=False, type=float, default=-1.0,
                        help="a real number specifying the lower Z-score for an aa type prediction to be considered as valid. (default %(default)s)",
                        metavar="<Z-score aa type assignment cutoff>")
    # parser.add_argument("-moreaa", dest="ALLOW_SINGLE_CH_PAIRS", required=False, action='store_true', help="allow aa type predictions from single C-H resonance pairs. Because it increases memory consumption and lowers findelity, it should be used in the last rounds of assignment.")
    parser.add_argument("-maxlen", dest="MAX_PEPTIDE_LENGTH", required=False, type=int, default=6,
                        help="The maximum peptide length (high values require more memory). Default: %(default)s.",
                        metavar="<max peptide length>")
    parser.add_argument("-minlen", dest="MIN_PEPTIDE_LENGTH", required=False, type=int, default=3,
                        help="the minimum peptide length (low values increase noise). Default: %(default)s",
                        metavar="<min peptide length>")
    parser.add_argument("-resoncut", dest="JUST_CARBON_MATCH_CUTOFF", required=False, type=int, default=3,
                        help="the minimum number of C-H resonance pairs for a TOCSY index group to start predicting aa types from both C-H or C \
                             only resonances. (default %(default)s)",
                        metavar="<Carbon prediction cutoff>")

    # parser.add_argument("-confile", dest="CONNECTIVITIES_FILE", required=False, help="connectivities file; if specified connectivity calculation from input files will be skipped", metavar="<connectivities file>")
    parser.add_argument("-poolconfile", dest="POOL_CONNECTIVITIES_FILE", required=False, default=None,
                        help="connectivities pool file; necessary if -confile specified in order to calculate correct probabilities",
                        metavar="<connectivities pool file>")
    parser.add_argument("-allconfile", dest="COMPLETE_CONNECTIVITIES_FILE", required=False, default=None,
                        help="all connectivities file; necessary if -confile specified in order to calculate correct probabilities",
                        metavar="<all connectivities file>")
    # parser.add_argument("-aafile", dest="AA_TYPES_FILE", required=False, help="amino acid assignment file; if specified amino acid assignment calculation from input files will be skipped", metavar="<amino acid assignment file>")
    parser.add_argument("-poolaafile", dest="POOL_AA_TYPES_FILE", required=False, default=None,
                        help="pool amino acid assignment file; if specified amino acid assignment calculation from input files will be skipped",
                        metavar="<pool amino acid assignment file>")
    parser.add_argument("-allaafile", dest="COMPLETE_AA_TYPES_FILE", required=False, default=None,
                        help="all amino acid assignment file; necessary if -aafile specified in order to calculate correct probabilities",
                        metavar="<all amino acid assignment file>")
    parser.add_argument("-chainfile", dest="NON_REDUNDANT_CHAINS_FILE", required=False, default=None,
                        help="non-redundant chains file; if specified, calculation of chains from connectivities will be omitted",
                        metavar="<non-redundant chains file>")
    parser.add_argument("-zmin", dest="MIN_NUM_OF_PREDICTIONS", required=False, type=int, default=4,
                        help="minimum number of aa type predictions required to apply the Z-score cutoff. The fewer the predictions \
                        the more inaccurate is Z-score. (default %(default)s)",
                        metavar="<minimum number of aa types prediction for Z-score filtering>")
    parser.add_argument("-transform", dest="TRANSFORM_TYPE", required=False, type=str, default='None',
                        help="type of mathematical transform to apply on the probabilities P[AAIG(i)|aatype(i-1)]. Allowed values \
                             are: \"None\", \"log\", \"log10\", \"boxcox_pearson\", \"boxcox_mle\". (default %(default)s)",
                        metavar="<type of mathematical transform>")
    parser.add_argument("-resume", dest="RESUME", required=False, action='store_true', default=False,
                        help="resume previous run by loading all peptide sequences saved in tmp_peptide_folder/. By default tmp_peptide_folder/ \
                             will be cleaned and new peptide sequences will be writen inside it. (default %(default)s)")
    parser.add_argument("-log", dest="LOG_TRANSFORM", required=False, action='store_true', default=False,
                        help="convert aa type prediction probabilities to logarithmic scale and then calculate Z-scores, if the \
                             min(probability)/max(probability) > 1000. (default %(default)s)")
    parser.add_argument("-skipchains", dest="SKIP_CHAINS", required=False, action='store_true', default=False,
                        help="Generate connectivies and amino acid type predictions, but skip chain and peptide formation. "
                             "(default %(default)s)")
    parser.add_argument("-proofread", dest="PROOFREAD", required=False, action='store_true', default=False,
                        help="Will print in red which correct connectivities were removed by MULTI_RATIO or"
                        " REVERSE CONNECTIVITY criterion. Default: %(default)s.")
    parser.add_argument("-delpred", dest="DELETE_AA_TYPE_PREDICTIONS", required=False, action='store_true',
                        default=False,
                        help="delete aa type predictions with probabilities which are 1000, 10000 or 100000 times lower than the highest, if \
                             the highest is >10e-10, >10e-20, <=10e-20, respectively. (default %(default)s)")
    parser.add_argument("-classifier", dest="CLASSIFIER_FILE", required=True, default=None, type=str,
                        help="path to the AA-type Classifier.",
                        metavar="<path to classifier>")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s version 1.0')
    args = parser.parse_args()
    return args


############################################### FUNCTION DEFINITIONS ##########################################
def create_NOESY_NOESY_connectivity_dictionaries(con,
                                                 MULTI_RATIO=1.0,
                                                 calc_intersection=True,
                                                 uniquify=True,
                                                 selected_AAIG_connectivities_dict=None,
                                                 reverse_conn=True,
                                                 proofread=False,
                                                 remove_matched_AAIGs_without_CA=True,
                                                 min_peak_num=1):
    """
    This method calculates HCNH-HCNH connectivities, returns them as dictionaries but also writes them into a text file
    named "%s_HCNH-HCNH_connectivities.list" % protein and a pickle file named "%s_HCNH-HCNH_connectivities.pkl" % protein
    which contains the 'i_iminus1_complete_dict' and 'NOESY_NOESY_peak_connectivities_mdict'.
    'i_iminus1_complete_dict' dictionary contains all possible connectivities of every NOESY AAIG1; has keys the
    AAIG signatures and values lists of quintuplets (tuples) consisting of (NOESY AAIG2 signature, # matching peaks,
    total # peaks, intersection, multi score).
    'NOESY_NOESY_peak_connectivities_mdict' is a 2D dictionary of the form: # AAIG1 -> matching AAIG2 ->
    [(matching NOESY peak1, matching NOESY peak 2, distance), ...], where peaks are Peak objects from lib/peak.py.


    :param con: HCNH_Connectivities() object
    :param tolH: proton tolerance in ppm to match AAIG1 peaks with AAIG2 peaks
    :param tolC: carbon tolerance in ppm to match AAIG1 peaks with AAIG2 peaks
    :param MULTI_RATIO: keep all connected AAIG2s that have 'multi' score up to MULTI_RATIO times lower than the maximum
                        'multi' score of the connectivities of AAIG1. The higher the MULTI_RATIO, the more connectivities
                        are retained.
    :param calc_intersection: keep it True for the moment
    :param uniquify: keep it True for the moment
    :param selected_AAIG_connectivities_dict: keep it None for the moment.
    :return:
    """

    con.match_NOESY_to_NOESY_Peaks(calc_intersection=calc_intersection,
                                  uniquify=uniquify,
                                  selected_AAIG_connectivities_dict=selected_AAIG_connectivities_dict,
                                  remove_matched_AAIGs_without_CA=remove_matched_AAIGs_without_CA,
                                  min_peak_num=min_peak_num)
    i_iminus1_complete_dict = con.__get_NOESY_NOESY_connectivities_dict__()  # all possible connectivities in dictionary form

    # Filter the connectivities by applying the MULTI_RATIO and the reverse_conn criteria.
    i_iminus1_dict = con.__get_filtered_NOESY_NOESY_connectivities_dict__(MULTI_RATIO=MULTI_RATIO,
                                                                          reverse_conn=reverse_conn,
                                                                          proofread=proofread)
    return i_iminus1_dict, i_iminus1_complete_dict, con


def check_input_files_format():
    """
    Method to check input files' format for correctness.

    :return:
    """
    global args

    # Check NOESY file
    query_contents = []
    HAS_INTENSITIES = False
    with open(args.HCNH_FILE, 'r') as f:
        contents = f.readlines()
        for line in contents:  # CLEANING THE SPECTRUM FROM IRRELEVANT LINES
            word_list = line.split()
            try:
                if not word_list[0] == '?-?-?-?':
                    continue
                float(word_list[1]);
                float(word_list[2]);
                float(word_list[3]);
                float(word_list[4]);  # checking if the line contains numbers (resonances)
                if len(word_list) == 6:
                    float(word_list[5])
                    HAS_INTENSITIES = True
                query_contents.append(line)
            except (IndexError, ValueError):
                # print "WARNING: Discarding", spectrum_type, "line:", line
                print(
                    bcolors.FAIL + "DEBUG check_input_files_format: this NOESY line will be discarded: " + line + bcolors.ENDC)
                continue
    # Now check if all the lines contain intensity
    if HAS_INTENSITIES:
        for word_list in query_contents:
            if len(word_list) < 6:
                print(
                    bcolors.FAIL + "ERROR: Please check your input NOESY file. Not all lines have intensity!!!" + bcolors.ENDC)
                sys.exit(1)

################################################################################################################

if __name__ == "__main__":

    try:

        ############################################ SANITY CHECKS ############################################
        args = cmdlineparse()
        print("Input argument values:")
        for arg in vars(args):
            print(arg, "=", getattr(args, arg))
        if args.H_weight <= 0 or args.H_weight > 1:
            print("ERROR: -wH must be a number greater than 0 and lower or equal to 1!")
            sys.exit(1)
        if args.C_weight <= 0 or args.C_weight > 1:
            print("ERROR: -wC must be a number greater than 0 and lower or equal to 1!")
            sys.exit(1)
        # if (args.ASSIGNMENT_CUTOFF != None and args.POOL_AA_TYPES_FILE != None) or (
        #         args.ASSIGNMENT_CUTOFF != None and args.COMPLETE_AA_TYPES_FILE != None):
        #     print("ERROR: Usage of -poolaafile or -allaafile is incompatible with -acutoff argument! Either remove -poolaafile and -allaafile to make new amino acid \
        #     type predictions from scratch, or remove -acutoff to use the amino acid prediction files you provided as input.")
        #     sys.exit(1)   # NOT SURE IF APPLICABLE IN NOESY!
        # if args.ASSIGNMENT_CUTOFF == None and args.POOL_AA_TYPES_FILE == None and args.COMPLETE_AA_TYPES_FILE == None:
        #     print("ERROR: you must either provide an aa type assignment cutoff with -acutoff argument (recommented value: 1.0), or provide a pool amino acid assignment file with \
        #     argument -poolaafile and file with all the aa type assignment with -allaafile argument.")
        #     sys.exit(1)   # NOT SURE IF APPLICABLE IN NOESY!

        # Launch input file check
        check_input_files_format()

        # if args.EXAMPLES:
        #    print "USAGE: python HCTOCSYNH_to_HCNOENH_tree-based.py -tocsy <HCTOCSYNH file> -noesy <HCNOENH file> [-stdH <stdev1>] [-stdC <stdev2>] [-update]"
        #    print "Example: ./4D_assignment_tree-based.py -hsqc TDhsqc.list -tocsy TDtocsy.list -noesy TDnoesy.list -stdH 0.03 -stdC 0.3 -rstdH 0.02 -rstdC 0.2 -cutoff 1.0"
        #    sys.exit(0)
        if args.DOWNLOAD_CS_HISTOGRAMS:
            download_CS_histograms()

        if not args.WORKDIR:
            args.WORKDIR = os.getcwd()
        elif not os.path.exists(args.WORKDIR):
            os.mkdir(args.WORKDIR)
        ############################################## END OF SANITY CHECKS ###################################

        # ENTER THE WORKDIR
        os.chdir(args.WORKDIR)

        # LOAD THE HSQC SPECTRUM AND GRADUALLY REDUCE TOLERANCES TO FIND MORE NON-OVERLAPPING GROUPS IN ROOT SPECTRUM
        HSQC_spec = HSQC_spectrum(args.HSQC_FILE, fasta=args.FASTA_FILE)

        ## COPY AAIGs FROM ROOT SPECTRUM (HSQC->NOESY matching)
        # In this scenario we don't have a TOCSY spectrum to find the side-chain peaks, so we pass an empty list
        # to 'TOCSY_sidechain_resonances_list='.
        NOESY_lines, _ = \
            HSQC_spec.copy_AAIGnames_from_HSQC_spectrum(spectrum4D_fname=args.HCNH_FILE,
                                                        spectrum4D_type="HCNH",
                                                        TOCSY_sidechain_resonances_list=[])
        ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
        ##        FIND OR LOAD HCNH-HCNH CONNECTIVITIES             ##
        ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
        noesy_con = HCNH_Connectivities()
        noesy_con.load_HCNH_spectrum(HSQC_spec.numlist_files['HCNH'])
        write_spectrum_info(noesy_con.HCNH_spec)  # write statistics about TOCSY and NOESY spectra
        noesy_con.__scale_intensities__()  # scale the intensities wrt the whole spectrum maximum
        if not args.POOL_CONNECTIVITIES_FILE:
            ## FIND HCNH-HCNH CONNECTIVITIES
            ColorPrint("Calculating HCNH-HCNH connectivities including native peaks...", "BOLDGREEN")
            # The function below will also write
            i_iminus1_dict, \
            i_iminus1_complete_dict,\
            noesy_con = \
                create_NOESY_NOESY_connectivity_dictionaries(noesy_con,
                                                             MULTI_RATIO=args.MULTI_RATIO,
                                                             reverse_conn=args.REVERSE_CONN,
                                                             proofread=proofread)
            write_NOESY_NOESY_connectivities(noesy_con.NOESY_NOESY_peak_connectivities_mdict,
                                             noesy_con.AAIG_connectivities_complete_dict, "connectivities_all",
                                             selected_AAIG_connectivities_dict=i_iminus1_complete_dict)
            write_NOESY_NOESY_connectivities(noesy_con.NOESY_NOESY_peak_connectivities_mdict,
                                             noesy_con.AAIG_connectivities_complete_dict,
                                             "connectivities_mratio_" + str(args.MULTI_RATIO),
                                             selected_AAIG_connectivities_dict=i_iminus1_dict)
        else:
            print("Loading Connectivities from file", args.POOL_CONNECTIVITIES_FILE)
            # Now load the pool connectivities file to calculate correct probabilities
            # same dictionary with the pool of possible connectivities remained after filtering using absolute consensus matches of a previous run
            i_iminus1_pool_dict = noesy_con.load_connectivities_from_file(args.POOL_CONNECTIVITIES_FILE)
            # Now load the complete connectivities file to calculate correct probabilities
            # the same dictionary but with all possible connectivities, including those below args.RESONANCE_MATCH_CUTOFF
            i_iminus1_complete_dict = noesy_con.load_connectivities_from_file(args.COMPLETE_CONNECTIVITIES_FILE)

            # FIND NEW CONNECTIVITIES FOR THOSE AAIGs THAT WERE NOT CORRECTED
            # FIRST QUICK RECONNAISSANCE ROUND: find missing connectivities
            ColorPrint("FIRST QUICK RECONNAISSANCE ROUND: find missing connectivities.", "BOLDGREEN")
            new_i_iminus1_dict, new_i_iminus1_complete_dict, noesy_con \
                                                = create_NOESY_NOESY_connectivity_dictionaries(noesy_con,
                                                                                               calc_intersection=False,
                                                                                               uniquify=False)  # new_i_iminus1_dict is useless here!
            selected_i_iminus1_dict = {}  # dict with the missing connectivities to be calculated later
            for AAIG_signature in list(new_i_iminus1_complete_dict.keys()):
                if AAIG_signature not in list(
                        i_iminus1_complete_dict.keys()):  # if there was not connectivity for this AAIG add it to the dictionaries
                    selected_i_iminus1_dict[AAIG_signature] = new_i_iminus1_complete_dict[AAIG_signature]
                elif len(new_i_iminus1_complete_dict[AAIG_signature]) > len(i_iminus1_complete_dict[
                                                                           AAIG_signature]):  # if there were found extra connectivities using the new tolerances
                    if not AAIG_signature in list(i_iminus1_pool_dict.keys()):  # probably all NOESY AAIGs which had connectivities have been mapped
                        # to the sequence
                        # TODO: Maybe it would make more sense to add only the new connenctivites that did not exist in
                        # TODO: i_iminus1_complete_dict to selected_i_iminus1_dict.
                        selected_i_iminus1_dict[AAIG_signature] = new_i_iminus1_complete_dict[
                            AAIG_signature]  # update the connectivities
                    elif len(i_iminus1_pool_dict[AAIG_signature]) == 1 and len(i_iminus1_complete_dict[AAIG_signature]) > 1:
                        are_all_unique = True  # do all AAIGs in the complete dict have a single connectivity in the pool
                        for tmp_quintuplet in i_iminus1_complete_dict[AAIG_signature]:
                            tmp_NOESYAAIG = tmp_quintuplet[0]
                            if tmp_NOESYAAIG in list(i_iminus1_pool_dict.keys()) and len(i_iminus1_pool_dict[tmp_NOESYAAIG]) > 1:
                                are_all_unique = False
                                break
                        if are_all_unique == False:  # this means that i_iminus1_pool_dict[AAIG_signature] has been modified and not cleaned from AAIGs used elsewhere
                            continue  # leave the connectivity of this AAIG_signature unchanged
                        selected_i_iminus1_dict[AAIG_signature] = new_i_iminus1_complete_dict[AAIG_signature]
                    elif len(i_iminus1_pool_dict[AAIG_signature]) == 1 and len(i_iminus1_complete_dict[AAIG_signature]) == 1:
                        selected_i_iminus1_dict[AAIG_signature] = new_i_iminus1_complete_dict[
                            AAIG_signature]  # update the connectivities

            ColorPrint("SECOND ROUND: calculate intersections of missing connectivities and add them.", "BOLDGREEN")
            noesy_con.clean_connectivities()  # VERY IMPORTANT!!!
            # 2D-hist intersections will be calculated only for the missing connectivities in selected_i_iminus1_dict
            new_i_iminus1_dict, \
            new_i_iminus1_complete_dict,\
            noesy_con = \
                        create_NOESY_NOESY_connectivity_dictionaries(noesy_con,
                                                                       MULTI_RATIO=args.MULTI_RATIO,
                                                                       reverse_conn=args.REVERSE_CONN,
                                                                       proofread=proofread,
                                                                       selected_AAIG_connectivities_dict=selected_i_iminus1_dict)  # new_i_iminus1_dict is useless here!
            for AAIG_signature in list(new_i_iminus1_complete_dict.keys()):
                if AAIG_signature not in list(i_iminus1_complete_dict.keys()):  # if there was not connectivity for this AAIG add it to the dictionaries
                    i_iminus1_complete_dict[AAIG_signature] = new_i_iminus1_complete_dict[AAIG_signature]  # update the connectivities
                    i_iminus1_pool_dict[AAIG_signature] = new_i_iminus1_complete_dict[AAIG_signature]
                elif len(new_i_iminus1_complete_dict[AAIG_signature]) > len(i_iminus1_complete_dict[AAIG_signature]):  # if there were found extra connectivities using the new tolerances
                    if not AAIG_signature in list(i_iminus1_pool_dict.keys()):  # probably all NOESY AAIGs which had connectivities have been mapped
                        # to the sequence
                        # TODO: Maybe it would make more sense to add only the new connenctivites that did not exist in
                        # TODO: i_iminus1_complete_dict to i_iminus1_pool_dict and i_iminus1_complete_dict.
                        i_iminus1_complete_dict[AAIG_signature] = new_i_iminus1_complete_dict[AAIG_signature]  # update the connectivities
                        i_iminus1_pool_dict[AAIG_signature] = new_i_iminus1_complete_dict[AAIG_signature]
                    elif len(i_iminus1_pool_dict[AAIG_signature]) == 1 and len(i_iminus1_complete_dict[AAIG_signature]) > 1:
                        are_all_unique = True  # do all AAIGs in the complete dict have a single connectivity in the pool
                        for tmp_quintuplet in i_iminus1_complete_dict[AAIG_signature]:
                            tmp_NOESYAAIG = tmp_quintuplet[0]
                            if tmp_NOESYAAIG in list(i_iminus1_pool_dict.keys()) and len(i_iminus1_pool_dict[tmp_NOESYAAIG]) > 1:
                                are_all_unique = False
                                break
                        if are_all_unique == False:  # this means that i_iminus1_pool_dict[AAIG_signature] has been modified and not cleaned from AAIGs used elsewhere
                            continue  # leave the connectivity of this AAIG_signature unchanged
                        i_iminus1_complete_dict[AAIG_signature] = new_i_iminus1_complete_dict[AAIG_signature]  # update the connectivities
                        i_iminus1_pool_dict[AAIG_signature] = new_i_iminus1_complete_dict[AAIG_signature]
                    elif len(i_iminus1_pool_dict[AAIG_signature]) == 1 and len(i_iminus1_complete_dict[AAIG_signature]) == 1:
                        i_iminus1_complete_dict[AAIG_signature] = new_i_iminus1_complete_dict[AAIG_signature]  # update the connectivities
                        i_iminus1_pool_dict[AAIG_signature] = new_i_iminus1_complete_dict[AAIG_signature]

            # CALCULATING CONNECTIVITIES FROM THE POOL OF CONNECTIVITIES USING THE args.MULTI_RATIO
            # print("Finding Connectivities within the -mratio from the Pool of Connectivities...")
            # # connectivities_mdict is a multidimensional dictionary with structure:
            # # NOESY AAIG => NOESY AAIG => [number of matched NOESY w2,w3 resonances in NOESY for that NOESY AAIG,
            # # total number of NOESY w2,w3 resonances for that NOESY AAIG]
            # i_iminus1_dict = {}  # a dictionary containing all possible connectivities of every NOESY AAIG above 'multi' threshold;
            # # has keys the NOESY AAIG signatures and values lists of quintuplets (tuples) \
            # # consisting of (NOESY AAIG_signature, occupancy, numOfResonances)
            # for AAIG_signature in list(i_iminus1_pool_dict.keys()):
            #     # Iterate over all NOESY aa indices and, if applicable, keep those that match in all resonance pairs w2,w3
            #     max_multi = np.max([float(q[4]) for q in i_iminus1_pool_dict[AAIG_signature]])
            #     for quintuplet in i_iminus1_pool_dict[AAIG_signature]:
            #         multi = float(quintuplet[4])
            #         if max_multi / multi <= args.MULTI_RATIO and multi > 10e-6:  # adjustable threshold
            #             try:
            #                 i_iminus1_dict[AAIG_signature].append(quintuplet)
            #             except KeyError:
            #                 i_iminus1_dict[AAIG_signature] = [(quintuplet)]
            # TODO: modify this function to take i_iminus1_pool_dict as input and not use the self.AAIG_connectivities_complete_dict
            # Filter the connectivities by applying the MULTI_RATIO and the reverse_conn criteria.

            i_iminus1_dict = noesy_con.__get_filtered_NOESY_NOESY_connectivities_dict__(MULTI_RATIO=args.MULTI_RATIO,
                                                                                  reverse_conn=args.REVERSE_CONN,
                                                                                  proofread=proofread,
                                                                                  pool_AAIG_connectivities_dict=i_iminus1_pool_dict)
            
            noesy_con.write_NOESY_NOESY_connectivities_from_dict(i_iminus1_dict,
                                                                 "connectivities_mratio_" + str(args.MULTI_RATIO))
            noesy_con.write_NOESY_NOESY_connectivities_from_dict(i_iminus1_complete_dict,
                                                       "connectivities_all")
            save_pickle("all_connectivities_dicts.pkl", i_iminus1_complete_dict, i_iminus1_pool_dict,
                        new_i_iminus1_complete_dict, new_i_iminus1_dict, i_iminus1_dict)
        ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
        ##              FORM CHAINS FROM CONNECTIVITIES               ##
        ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
        ColorPrint("Forming chains.", "BOLDGREEN")
        chain = Chain(i_iminus1_dict=i_iminus1_dict,
                      i_iminus1_complete_dict=i_iminus1_complete_dict)
        chain.form_chains(MAX_PEPTIDE_LENGTH=args.MAX_PEPTIDE_LENGTH,
                          MIN_PEPTIDE_LENGTH=args.MIN_PEPTIDE_LENGTH)
        chain.__add_reversed_chains__()   # Comment this out if you was only i->i-1 direction
        ColorPrint("Saving all the chains to chains.list file ...", "BOLDBLUE")
        chain.write_chains_file(chain.all_chainScore_set, "%s/chains.list" % args.WORKDIR)
        ColorPrint("Removing redundant chains and saving the rest to file chains.non-redundant.list...", "BOLDBLUE")
        all_chainScore_list = chain.remove_redundant_chains()
        # EXAMPLE: all_chainScore_list = [ ('V126NH', 'S99NH', 'A111NH', 'T128NH', 7.187111855209132e-13), ...]
        chain.write_chains_file(all_chainScore_list, "%s/chains.non-redundant.list" % args.WORKDIR)

        ## CALCULATE AA-TYPE PROBABILITIES
        prob = Probability(fasta_file=args.FASTA_FILE, classifier_file=args.CLASSIFIER_FILE)
        if not args.POOL_AA_TYPES_FILE and not args.COMPLETE_AA_TYPES_FILE:
            ColorPrint("\nPredicting the possible amino acid types of each AAIG.", "BOLDGREEN")
            prob.calc_NOESY_aa_type_conditional_probs(noesy_spec=noesy_con.HCNH_spec,
                                                      transform_type="None")
            # The Pool in the very beginning contains ALL predictions!
            AAIG_aaTypesProbPoolTupleList_dict = prob.AAIG_aaTypesProbTupleList_dict
            # Save ALL amino acid type predictions with the respective probabilities
            Probability.save_aa_type_predictions_to_file(prob.AAIG_aaTypesProbTupleList_dict,
                                                         "amino_acid_type_prediction_probabilities",
                                                         spectrum_type="NOESY")
        else:
        ################################################################################################################
        ##  LOAD AMINO ACID TYPE PREDICTIONS FOR FASTER DEBUGGING
        ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            print("\nLoading amino acid prediction files ...\n")
            AAIG_aaTypesProbPoolTupleList_dict = Probability.load_aatype_probs_from_file(
                                                                    fname=args.POOL_AA_TYPES_FILE)
            prob.calc_NOESY_aa_type_conditional_probs(noesy_spec=noesy_con.HCNH_spec,
                                                      COMPLETE_AA_TYPES_FILE=args.COMPLETE_AA_TYPES_FILE,
                                                      transform_type="None")

        ## FORM PEPTIDES
        # AAIG_aaTypesCutoffProbTupleList_dict = prob.get_pruned_condprob_dict(prob_cutoff=0.1)    # <== CHANGE ME
        AAIG_aaTypesCutoffProbTupleList_dict,\
            AAIG_aaTypesZscoreTupleList_dict = \
            Probability.get_pruned_prob_dicts_by_zscore(AAIG_aaTypesProbPoolTupleList_dict,
                                                    ZSCORE_ASSIGNMENT_CUTOFF=args.ZSCORE_ASSIGNMENT_CUTOFF,
                                                    DELETE_AA_TYPE_PREDICTIONS=args.DELETE_AA_TYPE_PREDICTIONS,
                                                    LOG_TRANSFORM=args.LOG_TRANSFORM,
                                                    MIN_NUM_OF_PREDICTIONS=args.MIN_NUM_OF_PREDICTIONS)

        ## SAVE THE AA-TYPE PROBABILITIES TO FILES FOR EASY DEBUGGING
        # Save amino acid type predictions that are above the Z-Score cutoff, along with the respective probabilities
        Probability.save_aa_type_predictions_to_file(AAIG_aaTypesCutoffProbTupleList_dict,
                                                     "amino_acid_type_prediction_probabilities_above_cutoff.zacutoff" + str(
                                                         args.ZSCORE_ASSIGNMENT_CUTOFF),
                                                     spectrum_type="NOESY")
        # Save amino acid type predictions with the respective Z-scores
        Probability.save_aa_type_predictions_to_file(AAIG_aaTypesZscoreTupleList_dict,
                                                     "amino_acid_type_prediction_Z-scores.zacutoff" + str(
                                                         args.ZSCORE_ASSIGNMENT_CUTOFF),
                                                     spectrum_type="NOESY")
        # Save ALL amino acid type predictions with the respective conditional probabilities
        Probability.save_aa_type_predictions_to_file(prob.AAIG_aaType_Condprob_mdict,
                                                     fname="amino_acid_type_prediction_conditional_probabilities",
                                                     spectrum_type="NOESY")

        # For DEBUGGING
        save_pickle("%s/aaType_prob_dicts.pkl" % args.WORKDIR, AAIG_aaTypesZscoreTupleList_dict, AAIG_aaTypesCutoffProbTupleList_dict)
        """
        # For proof-reading:
        for k, v in AAIG_aaTypesZscoreTupleList_dict.items(): 
            aa_type = get_aatype_from_AAIG_signature(k) 
            print(np.any([aa_type==d[0] for d in v])) 
        """
        pept = Peptide(total_chain_number=len(all_chainScore_list),
                            min_peptide_length=args.MIN_PEPTIDE_LENGTH,
                            AAIG_aaType_Condprob_mdict=prob.AAIG_aaType_Condprob_mdict,
                            AAIG_aaTypesCutoffProbTupleList_dict=AAIG_aaTypesCutoffProbTupleList_dict,
                            AAIG_aaTypesTransformedProbTupleList_dict=prob.AAIG_aaTypesTransformedProbTupleList_dict,
                            aatype_P_dict=prob.aatype_P_dict,
                            connectivity_dict=i_iminus1_dict,
                            spectrum_type="NOESY",
                            workdir=args.WORKDIR)
        pept.build_all_peptide_trees(all_chainScore_list, resume=args.RESUME, is_parallel=True) # default RESUME=False
        # Repeat peptide building for a second time to create whatever was left out accindentally.
        pept.build_all_peptide_trees(all_chainScore_list, resume=True, is_parallel=True)
        # Finally extract the peptide from the pickle files and write them to text files
        Peptide.write_all_peptides_to_files(all_chainScore_list, args.MAX_PEPTIDE_LENGTH, workdir=args.WORKDIR)

    except:
        type, value, tb = sys.exc_info()
        lines = traceback.format_exception(type, value, tb)
        print((''.join(lines)))
        raise
