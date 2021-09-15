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

import os
from argparse import ArgumentParser, RawDescriptionHelpFormatter, SUPPRESS

# Import 4D-CHAINS libraries
# CHAINS_BIN_DIR = os.path.dirname(os.path.realpath(__file__))
from lib.alignment import *
from lib.onlynoesy_csa import *
from lib.probhist import *
from lib.global_func import *
from lib.csa_for_tocsy import *

## Parse command line arguments
def cmdlineparse():
    parser = ArgumentParser(description="""The script to assign chemical shifts from only the NOESY spectrum.""",
        epilog="""
EXAMPLE: 
cs_assignment_onlyNOESY.py 
-nhmap 4DCHAINS_NHmap 
-rstart 89 
-noesy Tudor_NOESYnum.list 
-probprod 
-2dhist
        """)
    parser.add_argument("-nhmap", dest="NHmap_FILE", required=False, 
                        help="file with NH-mappings from previous run",
                        metavar="<NH-mapping file>")
    parser.add_argument("-fasta", dest="FASTA", required=False, type=str,
                        help="Fasta sequence file of the protein. Necessary only if the NHmap file comes from the "
                             "GrapH algorithm",
                        metavar="<fasta>")
    parser.add_argument("-rstart", dest="FIRST_RESIDUE_NUMBER", required=False, default=1, type=int, 
                        help="The number of the first residue in the protein sequence inside the NHmap_FILE (default: 1)",
                        metavar="<first residue number>")
    parser.add_argument("-noesy", dest="HCNH_FILE", required=True, help="the 4D HCNH file (*num.list) produced by 4D_assignment_parallel.py", metavar="<4D HCNH input file>")
    parser.add_argument("-o", dest="OUT_fname", required=False, default=None,
                        help="output file name (assigned 4D TOCSY in sparky format)", metavar="<output file>")
    parser.add_argument("-wH", dest="H_weight", required=False, type=float, default=0.1,
                        help="weight in (0.0-1.0] to apply on aa type prediction from aliphatic H resonances. (default: %(default)s)",
                        metavar="<H weight>")
    parser.add_argument("-wC", dest="C_weight", required=False, type=float, default=1.0,
                        help="weight in (0.0, 1.0] to apply on aa type prediction from aliphatic C resonances. (default: %(default)s)",
                        metavar="<C weight>")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s version 0.9')
    parser.add_argument("-probprod", dest="PROB_PRODUCT", required=False, action='store_true', default=False,
                        help="select the best C-H type assignment combination based on the product of probabilities of the individual C-H assignments")
    parser.add_argument("-2dhist", dest="USE_2D_HISTOGRAMS", required=False, action='store_true', default=False,
                        help="use 2D BMRB histograms for aa type prediction")
    parser.add_argument("-wcthres1", dest="WITHIN_CGROUP_THRESHOLD_ITER1", required=False, default=0.0, type=float,
                        help="The first Carbon type must have probability greater this ratio from the second one of the SAME CARBON GROUP in order to be kept in the final HCNH assignments")
    parser.add_argument("-bcthres1", dest="BETWEEN_CGROUP_THRESHOLD_ITER1", required=False, default=10.0, type=float,
                        help="The same Carbon type of a Carbon group must have probability greater this ratio from the probability of the SAME CARBON TYPE in any other Carbon group in order to be kept in the final HCNH assignments")
    parser.add_argument("-wcthres2", dest="WITHIN_CGROUP_THRESHOLD_ITER2", required=False, default=0.0, type=float,
                        help="same as -wcthres1, but for the 2nd iteration")
    parser.add_argument("-bcthres2", dest="BETWEEN_CGROUP_THRESHOLD_ITER2", required=False, default=10.0, type=float,
                        help="same as -bcthres1, but for the 2nd iteration")
    parser.add_argument("-wcthres3", dest="WITHIN_CGROUP_THRESHOLD_ITER3", required=False, default=0.0, type=float,
                        help="same as -wcthres1, but for the 3rd iteration")
    parser.add_argument("-bcthres3", dest="BETWEEN_CGROUP_THRESHOLD_ITER3", required=False, default=0.0, type=float,
                        help="same as -bcthres1, but for the 3rd iteration")
    parser.add_argument("-wcthres4", dest="WITHIN_CGROUP_THRESHOLD_ITER4", required=False, default=0.0, type=float,
                        help="same as -wcthres1, but for the 4th iteration")
    parser.add_argument("-bcthres4", dest="BETWEEN_CGROUP_THRESHOLD_ITER4", required=False, default=0.0, type=float,
                        help="same as -bcthres1, but for the 4th iteration")
    parser.add_argument("-wcthres5", dest="WITHIN_CGROUP_THRESHOLD_ITER5", required=False, default=0.0, type=float,
                        help="same as -wcthres1, but for the 5th iteration")
    parser.add_argument("-bcthres5", dest="BETWEEN_CGROUP_THRESHOLD_ITER5", required=False, default=0.0, type=float,
                        help="same as -bcthres1, but for the 5th iteration")
    parser.add_argument("-percentile1", dest="PERCENTILE_ITER1", required=False, default=0.85, type=float,
                        help="This arguments is used as a theshold to discard low probability C-H type predictions in iteration 1. Can be 0.9, 0.85, 0.8, 0.75, or 0.0 (no filtering). A 0.9 percentile of value X means that the top 10 percent of the probability density values were above X.")
    parser.add_argument("-percentile2", dest="PERCENTILE_ITER2", required=False, default=0.9, type=float,
                        help="This arguments is used as a theshold to discard low probability C-H type predictions in iteration 2. Can be 0.9, 0.85, 0.8, 0.75, or 0.0 (no filtering). A 0.9 percentile of value X means that the top 10 percent of the probability density values were above X.")
    parser.add_argument("-percentile3", dest="PERCENTILE_ITER3", required=False, default=0.8, type=float,
                        help="This arguments is used as a theshold to discard low probability C-H type predictions in iteration 3. Can be 0.9, 0.85, 0.8, 0.75, or 0.0 (no filtering). A 0.9 percentile of value X means that the top 10 percent of the probability density values were above X.")
    parser.add_argument("-percentile4", dest="PERCENTILE_ITER4", required=False, default=0.8, type=float,
                        help="This arguments is used as a theshold to discard low probability C-H type predictions in iteration 4. Can be 0.9, 0.85, 0.8, 0.75, or 0.0 (no filtering). A 0.9 percentile of value X means that the top 10 percent of the probability density values were above X.")
    parser.add_argument("-percentile5", dest="PERCENTILE_ITER5", required=False, default=0.8, type=float,
                        help="This arguments is used as a theshold to discard low probability C-H type predictions in iteration 5. Can be 0.9, 0.85, 0.8, 0.75, or 0.0 (no filtering). A 0.9 percentile of value X means that the top 10 percent of the probability density values were above X.")
    parser.add_argument("-usertocsy", dest="user_TOCSY_FILE", required=False,
                        help="4D TOCSY (HCTOCSYNH) file in Sparky format with atom assignments made by the user", metavar="<Sparky 4D TOCSY input file with user-made atom assignments>")
    parser.add_argument("-usernoesy", dest="user_HCNH_FILE", required=False,
                        help="4D HCNH (NOESY) file in Sparky format with atom assignments made by the user", metavar="<Sparky 4D HCNH input file with user-made atom assignments>")
    parser.add_argument("-int1", dest="USE_INTENSITIES_ITERATION1", required=False, action='store_true', default=False,
                        help="use peak intensities in iteration 1")
    parser.add_argument("-int2", dest="USE_INTENSITIES_ITERATION2", required=False, action='store_true', default=False,
                        help="use peak intensities in iteration 2")
    parser.add_argument("-int3", dest="USE_INTENSITIES_ITERATION3", required=False, action='store_true', default=False,
                        help="use peak intensities in iteration 3")
    parser.add_argument("-int4", dest="USE_INTENSITIES_ITERATION4", required=False, action='store_true', default=False,
                        help="use peak intensities in iteration 4")
    parser.add_argument("-int5", dest="USE_INTENSITIES_ITERATION5", required=False, action='store_true', default=False,
                        help="use peak intensities in iteration 5")
    parser.add_argument("-ithres1", dest="INTENSITY_THRESHOLD_ITERATION1", required=False, type=float, default=0.1,
                        help="use relative peak intensity threshold (normalized wrt residue not the whole spectrum)")
    parser.add_argument("-ithres2", dest="INTENSITY_THRESHOLD_ITERATION2", required=False, type=float, default=0.0,
                        help="use relative peak intensity threshold (normalized wrt residue not the whole spectrum)")
    parser.add_argument("-ithres3", dest="INTENSITY_THRESHOLD_ITERATION3", required=False, type=float, default=0.1,
                        help="use relative peak intensity threshold (normalized wrt residue not the whole spectrum)")
    parser.add_argument("-ithres4", dest="INTENSITY_THRESHOLD_ITERATION4", required=False, type=float, default=0.0,
                        help="use relative peak intensity threshold (normalized wrt residue not the whole spectrum)")
    parser.add_argument("-ithres5", dest="INTENSITY_THRESHOLD_ITERATION5", required=False, type=float, default=0.0,
                        help="use relative peak intensity threshold (normalized wrt residue not the whole spectrum)")
    parser.add_argument("-itrans1", dest="INTENSITY_TRANSFORM_TYPE_ITER1", required=False, type=int, default=4,
                        help="transform the intensities in iteration1. Can be 1: do nothing; 2: (default: 1)")
    parser.add_argument("-itrans2", dest="INTENSITY_TRANSFORM_TYPE_ITER2", required=False, type=int, default=2,
                        help="transform the intensities in iteration2. Can be 1: do nothing; 2: (default: 1)")
    parser.add_argument("-itrans3", dest="INTENSITY_TRANSFORM_TYPE_ITER3", required=False, type=int, default=3,
                        help="transform the intensities in iteration3. Can be 1: do nothing; 2: (default: 1)")
    parser.add_argument("-itrans4", dest="INTENSITY_TRANSFORM_TYPE_ITER4", required=False, type=int, default=2,
                        help="transform the intensities in iteration4. Can be 1: do nothing; 2: (default: 1)")
    parser.add_argument("-itrans5", dest="INTENSITY_TRANSFORM_TYPE_ITER5", required=False, type=int, default=1,
                        help="transform the intensities in iteration5. Can be 1: do nothing; 2: (default: 1)")
    parser.add_argument("-probmodel", dest="PROBABILITY_MODEL", required=False, type=int, default=2,
                        help="If '1' the probability of each peak will be given by [wH*1Dhist(H)+wC*1Dhist(C)]/(wH+wC). If '2' then by 1Dhist(H)*1Dhist(C)")
    parser.add_argument("-probmode", dest="CONSENSUS_CGROUP_PROB_MODE", required=False, default=5, type=int,
                        help="""The way to calculate the total score of a set of chemical shift assignment (default: 3).
                    The following values control how to calculate the consensus probability of each C-group. The total score will be the
                    product of this consensus C-group probabilities.
                    1: average;
                    2: sqrt(prob1*prob2)    ; geometric mean
                    3: (prob2-prob1)/(log(prob2)-log(prob1))    ; logarithmic mean
                    4: (prob1*prob2)/((prob1+prob2)/2.) ;composite average 
                    5: prob1 + prob2    ; simple sum
                        """, metavar="<way to calculate the cs assignment score>")
    parser.add_argument("-debugresidues", dest="DEBUG_RESIDUES", required=False, type=str, default=None,
                        help="Only for debugging! You can give specific residues separated by ',' to keep, the rest will"
                             " be discarded from the input files. E.g. -debugresidues 'V124,V125' .",
                        metavar="<debug residues>")
    # OPTIONAL ARGUMENTS
    parser.add_argument("-ithresA", dest="ithresA", type=float, default=0.1, help=SUPPRESS) # intensity threshold for ALA CB-QB methyl for both iteration 1 & 3
    parser.add_argument("-ithresT", dest="ithresT", type=float, default=0.1, help=SUPPRESS) # intensity threshold for THR CG2-QG2 methyl for both iteration 1 & 3
    parser.add_argument("-ithresV", dest="ithresV", type=float, default=0.1, help=SUPPRESS) # intensity threshold for VAL CG1-QG1 & CG2-QG2 methyls for both iteration 1 & 3
    parser.add_argument("-ithresI", dest="ithresI", type=float, default=0.1, help=SUPPRESS) # intensity threshold for ILE CD1-QD1 & CG2-QG2 methyl for both iteration 1 & 3
    parser.add_argument("-ithresL", dest="ithresL", type=float, default=None, help=SUPPRESS) # intensity threshold for LEU CD1-QD1 & CD2-QD2 methyls for both iteration 1 & 3
    parser.add_argument("-pthresA", dest="pthresA", type=float, default=0.8, help=SUPPRESS) # intensity threshold for ALA CB-QB methyl for both iteration 1 & 3
    parser.add_argument("-pthresT", dest="pthresT", type=float, default=0.8, help=SUPPRESS) # intensity threshold for THR CG2-QG2 methyl for both iteration 1 & 3
    parser.add_argument("-pthresV", dest="pthresV", type=float, default=0.8, help=SUPPRESS) # intensity threshold for VAL CG1-QG1 & CG2-QG2 methyls for both iteration 1 & 3
    parser.add_argument("-pthresI", dest="pthresI", type=float, default=0.8, help=SUPPRESS) # intensity threshold for ILE CD1-QD1 & CG2-QG2 methyl for both iteration 1 & 3
    parser.add_argument("-pthresL", dest="pthresL", type=float, default=None, help=SUPPRESS) # intensity threshold for LEU CD1-QD1 & CD2-QD2 methyls for both iteration 1 & 3
    
    # Alternatively you can define one ithres and pthres for all Methyls. If args.ithresMethyl and args.pthresMethyl != None then the above 10
    # values will be overriden!
    parser.add_argument("-ithresMethyl", dest="ithresMethyl", type=float, default=0.1, help=SUPPRESS) # intensity threshold for all Methyls
    parser.add_argument("-pthresMethyl", dest="pthresMethyl", type=float, default=0.8, help=SUPPRESS) # intensity threshold for all Methyls
    
    args = parser.parse_args()
    return args

if __name__ == '__main__':
############################################### SANITY CHECKS #########################################

    args = cmdlineparse()
    print("Input argument values:")
    for arg in vars(args):
        print(arg, "=", getattr(args, arg))

    try:

        if args.ithresMethyl:
            args.ithresA, args.ithresT, args.ithresV, args.ithresI, args.ithresL = args.ithresMethyl, args.ithresMethyl, args.ithresMethyl, args.ithresMethyl, args.ithresMethyl
        if args.pthresMethyl:
            args.pthresA, args.pthresT, args.pthresV, args.pthresI, args.pthresL = args.pthresMethyl, args.pthresMethyl, args.pthresMethyl, args.pthresMethyl, args.pthresMethyl

        # Make file paths absolute (not relative, since you work with symlinks)
        args.NHmap_FILE = os.path.abspath(args.NHmap_FILE)
        if is_graph_NHmap(args.NHmap_FILE):
            assert args.FASTA, ColorPrint("ERROR: the input NHmap file comes from the GrapH algorithm, therefore "
                                          "you must provide the FASTA sequence of the protein (-fasta).", "FAIL")
            print_graph_results_summary(graph_output=args.NHmap_FILE,
                                        fasta=args.FASTA,
                                        out_fname=args.NHmap_FILE + "_4DCHAINS_format",
                                        silent=True)
            args.NHmap_FILE = args.NHmap_FILE + "_4DCHAINS_format"
        args.HCNH_FILE = os.path.abspath(args.HCNH_FILE)

############################################ END OF SANITY CHECKS #########################################


        ##########################################################################################################
        ##                                         LOADING FILES                                                ##
        ##########################################################################################################

        ncsa = HCNH_CSA(args, spectrum_combo="HCNH-HCNH")

        ########################################################## END OF FILE LOADING #########################################################


        ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ITERATION 1 only for methyl groups ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
        ncsa.onlynoesy_csa_iter1()

        ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ITERATION 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        ncsa.onlynoesy_csa_iter2()

        ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ITERATION 3 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        ncsa.onlynoesy_csa_iter3()

        ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ITERATION 4 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        ncsa.onlynoesy_csa_iter4()

        ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ITERATION 5 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        ## In iteration 5 only PHE-->CD-QD and TYR-->CD-QD,CE-QE are assigned
        ncsa.onlynoesy_csa_iter5()

        ##
        ## Write in HCNH Sparky format the i-1 and i+1 peaks, too. Write also unassigned peaks with just
        # the HSQC AAIG signature.
        ##
        matched_HCNH_residue_assignments_dict = load_pickle("iter5.pkl")[0]
        annotate_HCNH_file(args.HCNH_FILE,
                            ncsa.ali.absolute_AAIGmatches_alignment_list,
                            ncsa.ali.absolute_matches_alignment_list,
                            matched_HCNH_residue_assignments_dict,
                            ncsa.HCNH_residue_peak_intensity_mdict,
                            out_fname="only4DNOESY_assignedall.sparky")

    except:
        type, value, tb = sys.exc_info()
        lines = traceback.format_exception(type, value, tb)
        print((''.join(lines)))
        raise