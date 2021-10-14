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
from argparse import ArgumentParser
from collections import defaultdict, OrderedDict
from operator import itemgetter

import numpy as np
from lib.alignment import Alignment
from lib.csa_for_noesy import group_carbons, rename_protons

 
## Parse command line arguments
from lib.csa_for_tocsy import get_all_assignments_from_H_C_resonpair_2Dhist, get_aatypes_from_H_C_resonpair, \
    get_presenceprob_from_all_H_C_resonpairs
from lib.global_func import Debuginfo, approx_equal_proportional
from lib.global_vars import aa_carbonBondedHydrogensDict_dict, aatype_carbon_methylHydrogens_mdict, \
    aa_CarbonListwithDegenerateH_dict, aa3to1_dict, aa1to3_dict
from lib.hcnh_spectrum import HCNH_Spectrum
from lib.probhist import ProbHist_Loader
from lib.tocsy_spectrum import TOCSY_Spectrum


def cmdlineparse():
    parser = ArgumentParser(description="command line arguments",
        epilog="EXAMPLE: cs_assignment.py -nhmap results_summary.3mers_round1_rst.chainlinkers.patch -rstart 200 -tocsy tocsyHCNH.21.4.2016num.list -probprod -2dhist -noesy noesyHCNH.21.4.2016num.list")
    parser.add_argument("-nhmap", dest="NHmap_FILE", required=True,
                        help="file with absolute matches from previous run",
                        metavar="<absolute matches file>")
    parser.add_argument("-rstart", dest="FIRST_RESIDUE_NUMBER", required=True, default=1, type=int, 
                        help="The number of the first residue in the protein sequence inside the NHmap_FILE (default: 1)",
                        metavar="<first residue number>")
    parser.add_argument("-tocsy", dest="TOCSY_fname", required=True,
                        help="4D TOCSY (HCTOCSYNH) file", metavar="<4D TOCSY input file>")
    parser.add_argument("-noesy", dest="NOESY_FILE", required=False,
                        help="optionally the 4D NOESY file in case the -patch argument was used when running chain_linker.py",
                        metavar="<4D NOESY input file>")
    parser.add_argument("-o", dest="OUT_fname", required=True, default=None,
                        help="output file name (assigned 4D TOCSY in sparky format). (default: like -nhmap but with a different file extension)",
                        metavar="<output file>")
    parser.add_argument("-wH", dest="H_weight", required=False, type=float, default=0.1,
                        help="weight in (0.0-1.0] to apply on aa type prediction from aliphatic H resonances. (default: %(default)s)",
                        metavar="<H weight>")
    parser.add_argument("-wC", dest="C_weight", required=False, type=float, default=1.0,
                        help="weight in (0.0, 1.0] to apply on aa type prediction from aliphatic C resonances. (default: %(default)s)",
                        metavar="<C weight>")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s version 1.2')
    parser.add_argument("-probprod", dest="PROB_PRODUCT", required=False, action='store_true', default=False,
                        help="select the best C-H type assignment combination based on the product of probabilities of the individual C-H assignments. (default: %(default)s)")
    parser.add_argument("-probmodel", dest="PROBABILITY_MODEL", required=False, type=int, default=2,
                        help="If '1' the probability of each peak will be given by [wH*1Dhist(H)+wC*1Dhist(C)]/(wH+wC). If '2' then by \
                             1Dhist(H)*1Dhist(C). (default: %(default)s)")
    parser.add_argument("-probmode", dest="CONSENSUS_CGROUP_PROB_MODE", required=False, default=2, type=int,
                        help="""The way to calculate the total score of a set of chemical shift assignment. Can be:
                    0: just multiply all the probabilities of the individual peaks
                    The following values control how to calculate the consensus probability of each C-group. The total score will be the
                    product of this consensus C-group probabilities.
                    1: average;
                    2: sqrt(prob1*prob2)    ; geometric mean
                    3: (prob2-prob1)/(log(prob2)-log(prob1))    ; logarithmic mean
                    4: (prob1*prob2)/((prob1+prob2)/2.) ;composite average 
                    5: (log(prob1)+log(prob2))/2.
                     (default: %(default)s)
                        """, metavar="<way to calculate the cs assignment score>")
    parser.add_argument("-2dhist", dest="USE_2D_HISTOGRAMS", required=False, action='store_true', default=False,
                        help="use 2D BMRB histograms for aa type prediction. (default: %(default)s)")
    args=parser.parse_args()
    return args


args = cmdlineparse()
print("Input argument values:")
for arg in vars(args):
    print(arg, "=", getattr(args, arg))
            
# Make file paths absolute (not real, since you work with symlinks)
args.NHmap_FILE = os.path.abspath(args.NHmap_FILE)
args.TOCSY_fname = os.path.abspath(args.TOCSY_fname)
args.NOESY_FILE = os.path.abspath(args.NOESY_FILE)


##########################################################################################################
##                                         LOADING FILES                                                ##
##########################################################################################################

# LOAD PROBABILITY HISTOGRAMS
histload = ProbHist_Loader()
histload.load_1Dhistograms()
histload.load_2Dhistograms()


## READ FILE WITH ABSOLUTE MATCHES
# absolute_AAIGmatches_alignment_list contains the absolutely matched HSQC AAIGs, namely the names given to each
# N,HN peak in the HSQC. It does not necessarily contain real residue names!

absolute_AAIGmatches_alignment_list, protein_alignment_list = \
    Alignment.read_NHmap_file(args.NHmap_FILE, get_protein_alignment=True)

# Load TOCSY contents
# OLD WAY: TOCSY_AAIGsign_CSlist_dict = read_spectrum_file(args.TOCSY_fname)
tocsy_spec = TOCSY_Spectrum(TOCSY_fname=args.TOCSY_fname)   # don't read into a dataframe because it ignores non 'N-H' peaks
tocsy_spec.find_AAIG2residue_correspondaces(args.NHmap_FILE, starting_resid=args.FIRST_RESIDUE_NUMBER)  # DO NOT RENAME AAIG signatures!
# create local variables from tocsy_spec variable to save space
residue2TAAIGsign_dict = tocsy_spec.residue2AAIGsign_dict    # OrderedDict with keys: residue (resname{3 letters}+resid) -> assigned TOCSY AAIG signature
                                        # the TAAIG corresponds to the i+1 residue, as that residue "sees" the TOCSY peaks of i residue
residue2position_dict = tocsy_spec.residue2position_dict   # OrderedDict with keys: residue (resname+resid) -> position in protein sequence (starting from 0)
print("DEBUG: residue2TAAIGsign_dict=", residue2TAAIGsign_dict)

TOCSY_AAIGsign_CSlist_dict = defaultdict(list)  # j_AAIG_signature -> [[Hreson, Creson, Nreson, HNreson], [...], ...]
                                            # E.g.  {'X10NH': [['4.775', '56.946', '120.151', '8.920'],
                                            # ['3.006', '32.245', '120.153', '8.920'], ...]
for p in tocsy_spec.__get_all_peaks__():
    # ATTENTION: the aa type must be in 3-letter code. Some keys in TOCSY_AAIGsign_CSlist_dict will be real residue codes,
    # ATTENTION: while others will be virtual (e.g. 'X45' without 'NXHX').
    # try:
    #     residue = aa1to3_dict[p.j_AAIG_signature[0]] + remove_NH_suffix(p.j_AAIG_signature[1:])
    # except KeyError:
    #     residue = remove_NH_suffix(p.j_AAIG_signature)
    TOCSY_AAIGsign_CSlist_dict[p.j_AAIG_signature].append([p.Hreson, p.Creson, p.j_Nreson, p.j_HNreson])
print("DEBUG: TOCSY_AAIGsign_CSlist_dict=", TOCSY_AAIGsign_CSlist_dict)
residue_CSlist_dict = OrderedDict() # ordereddict with keys the resname+resid -> the list of the associated H-C-N-HN resonances. E.g.
# I202 --> [['4.998', '52.252', '125.616', '9.165'], ['2.638', '38.088', '125.620', '9.167'], ['2.116', '38.070', '125.609', '9.165']]
patched_residues_list = []  # list with the residues (e.g. E104, not TOCSY index group)that were added by -patch option in the chain_linker.py (they have no TOCSY peaks but have NOESY)
for residue in list(residue2TAAIGsign_dict.keys()):
    print("DEBUG: residue=", residue) 
    if residue2TAAIGsign_dict[residue] == None:
        print("DEBUG: residue2TAAIGsign_dict[residue]=", residue2TAAIGsign_dict[residue])
        continue
    try:
        residue_CSlist_dict[residue] = TOCSY_AAIGsign_CSlist_dict[residue2TAAIGsign_dict[residue]]
    except KeyError:    # in case this residue from the alignment has no TOCSY peaks (was added using -patch in the chain_linker.py), save it
        patched_residues_list.append(residue)
        continue

print("DEBUG: residue_CSlist_dict =", residue_CSlist_dict)
print("DEBUG: patched_residues_list=", patched_residues_list)
NOESY_residue_CSlist_dict = OrderedDict() # ordereddict with keys the resname+reside -> the list of the associated H-C-N-HN resonances. E.g.
# I202 --> [['4.998', '52.252', '125.616', '9.165'], ['2.638', '38.088', '125.620', '9.167'], ['2.116', '38.070', '125.609', '9.165']]
if args.NOESY_FILE:
    # OLD WAY: NOESY_AAIG_CSlist_dict = read_spectrum_file(args.NOESY_FILE)  # same as TOCSY_AAIGsign_CSlist_dict, but for the NOESY;
    #                                                                 # it includes the right peaks as well as peaks of nearby residues
    noesy_spec = HCNH_Spectrum(HCNH_FILE=args.NOESY_FILE)   # don't read into a dataframe because it ignores non 'N-H' peaks
    NOESY_AAIG_CSlist_dict = defaultdict(list)    # same as TOCSY_AAIGsign_CSlist_dict, but for the NOESY;
                                                        # it includes peaks of AAIG(i) as well as peaks of nearby residues
                                                        # E.g.  {'X10': [['4.775', '56.946', '120.151', '8.920'],
                                                        # ['3.006', '32.245', '120.153', '8.920'], ...]
    for p in noesy_spec.__get_all_peaks__():
        NOESY_AAIG_CSlist_dict[p.j_AAIG_signature].append([p.Hreson, p.Creson, p.j_Nreson, p.j_HNreson])

    for residue in list(residue2TAAIGsign_dict.keys()):
        if residue2TAAIGsign_dict[residue] == None:
            continue
        elif residue in patched_residues_list:
            try:
                NOESY_residue_CSlist_dict[residue] = NOESY_AAIG_CSlist_dict[residue2TAAIGsign_dict[residue]]
            except KeyError:    # in case this residue from the alignment has no TOCSY and no NOESY peaks (was added using -patch in the chain_linker.py), print error!
                raise KeyError(Debuginfo("ERROR: residue %s was patched but has no NOESY peaks!!!" %
                                         residue2TAAIGsign_dict[residue], fail=True))

    # Also rewrite the HCNH num.list files with the real residue names
    noesy_spec.rename_AAIGs_to_residues(args.NHmap_FILE, starting_resid=args.FIRST_RESIDUE_NUMBER)
    noesy_spec.write_sparky_list(os.path.realpath(args.NOESY_FILE).replace("num.list", "realnum.list"))

print("DEBUG: NOESY_residue_CSlist_dict=", NOESY_residue_CSlist_dict)
##print "DEBUG: TOCSY_AAIGsign_CSlist_dict =", TOCSY_AAIGsign_CSlist_dict
#for k,v in TOCSY_AAIGsign_CSlist_dict.items():
#    print k, "-->", v

##print "DEBUG: residue_CSlist_dict =", residue_CSlist_dict
#for k,v in residue_CSlist_dict.items():
#    print k, "-->", v


residue_residueTypesProbTupleList_dict = OrderedDict() # ordereddict with keys the AAIG of residue i and values lists of tuples of the form (residue i-1 matching aa type, average probability)
Num_of_TOCSY_resonances = 0 # save here the number of TOCSY resonances for a particular AAIG
residue = None
possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list = [] # list of lists of the form (aa type, Carbon name, Hydrogen name, weighted average probability, TOCSY_reson_index, H resonance, C resonance, carbon group)
previous_aaindex = ""
atom_index = 1
if args.OUT_fname:
    xeasy_fout_name = args.OUT_fname+".xeasy"
    sparky_fout = open(args.OUT_fname+".sparky", 'w+')
else:
    xeasy_fout_name = args.NHmap_FILE+".xeasy"
    sparky_fout = open(args.NHmap_FILE+".sparky", 'w')
xeasy_fout = open(xeasy_fout_name, 'w')

residue_list = list(residue_CSlist_dict.keys())
carbon_groups_list = []  # list of lists of the form [H resonance, C resonance, N resonance, HN resonance, carbon group]
for res_index, residue in enumerate(residue_list):
    #if residue == "GLN136":
    #try:
    print("Assigning chemical shifts to residue ", residue)
    aa_type = residue[0:3]  # aa type in 3-letter code
    #print "DEBUG: residue_CSlist_dict[residue] =", residue_CSlist_dict[residue]
    TOCSYindex_ResonancesTuple_dict = {}
    for TOCSY_Hreson,TOCSY_Creson,TOCSY_Nreson,TOCSY_HNreson in residue_CSlist_dict[residue]:
        Num_of_TOCSY_resonances += 1
        # print "DEBUG: TOCSY_Creson=",TOCSY_Creson,"TOCSY_Hreson=",TOCSY_Hreson
        # valid_matches_list: list of lists of the form: [aa type,
                                                        # Carbon name,
                                                        # Hydrogen name,
                                                        # weighted average probability,
                                                        # TOCSY_reson_index,
                                                        # H resonance,
                                                        # C resonance,
                                                        # carbon group]
        if args.USE_2D_HISTOGRAMS == True:
            valid_matches_list = get_all_assignments_from_H_C_resonpair_2Dhist(aa_type,
                                                                               float(TOCSY_Hreson),
                                                                               float(TOCSY_Creson),
                                                                               Num_of_TOCSY_resonances,
                                                                               histload,
                                                                               useonlyCarbon=False) # currently useonlyCarbon=False not supported by select_correct_H_C_resonpairs4
        else:
            valid_matches_list = get_aatypes_from_H_C_resonpair(aa_type,
                                                                float(TOCSY_Hreson),
                                                                float(TOCSY_Creson),
                                                                Num_of_TOCSY_resonances)
        TOCSYindex_ResonancesTuple_dict[Num_of_TOCSY_resonances] = (TOCSY_Hreson, TOCSY_Creson, TOCSY_Nreson, TOCSY_HNreson)
        #if TOCSY_aaindex == "Y130":
        print("DEBUG: valid_matches_list=",valid_matches_list)
        possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list.extend(valid_matches_list)
        print("point 0")
        carbon_groups_list.append([TOCSY_Hreson, TOCSY_Creson, TOCSY_Nreson, TOCSY_HNreson, None])
        print("point 1")
        
    # GROUP THE CARBON RESONANCE OF EACH TAAIG
    if len(possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list) == 0: # if no atom assignment could be made, skip this aa!
        print("WARNING: No CS assignment could be made for residue", residue, ", therefore I skip it!")
        # reset these variables for the next residue 
        possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list = [] # list of tuples of the form (aa type, probability, Hreson, Creson, TOCSY_reson_index) containing all matching sets of H,C resonances
        Num_of_TOCSY_resonances = 0
        continue
    print("point 2")
    new_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list = group_carbons(possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list, iteration=None)
    
    if len(new_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list) == 0: # if carbon grouping excluded all peaks, skip this aa!
        print("WARNING: No CS assignment could be made for residue", residue, ", therefore I skip it!")
        # reset these variables for the next residue 
        possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list = [] # list of tuples of the form (aa type, probability, Hreson, Creson, TOCSY_reson_index) containing all matching sets of H,C resonances
        Num_of_TOCSY_resonances = 0
        continue
    
    print("point 3")
    aatype_CHnucleiType_presenceProbSum_mdict, \
    aatype_CHnucleiType_TresonIndex_mdict = get_presenceprob_from_all_H_C_resonpairs(new_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list,
                                                                                        Num_of_TOCSY_resonances,
                                                                                        args)
    #print "DEBUG: previous_TOCSY_aaindex=",previous_TOCSY_aaindex,"residue_residueTypesProbTupleList_dict[previous_TOCSY_aaindex]=",residue_residueTypesProbTupleList_dict[previous_TOCSY_aaindex]
    print("point 4")
    possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list = [] # list of tuples of the form (aa type, probability, Hreson, Creson, TOCSY_reson_index) containing all matching sets of H,C resonances
    
    print("point 5")
    Num_of_TOCSY_resonances = 0
    
    if aatype_CHnucleiType_presenceProbSum_mdict == {} and aatype_CHnucleiType_TresonIndex_mdict == {}: # if no atom assignment could be made, skip this aa! (TEMPORARY FIX)
        print("WARNING: No CS assignment could be made for residue", residue, ", therefore I skip it!")
        continue

#except (ValueError, IndexError):
#    print "WARNING: the 3rd and 4th elements of the following TOCSY file line are not numbers:"
#    print "TOCSY_Hreson, TOCSY_Creson:", TOCSY_Hreson, TOCSY_Creson
#    sys.exit(1)
    
    # Write the Chemical Shifts Assignments in XEASY format for Rosetta
    xeasy_lines_list, sparky_lines_list = [], []    # lines to save the xeasy and sparky lines for modification
    HNresons_list, Nresons_list = [], []
    Creson_Cname_lines_list = [] # list with the names of the printed carbons to avoid double printing (e.g. CB, CG, CD, CE)
    print("DEBUG: aatype_CHnucleiType_TresonIndex_mdict[aa_type].keys() =", list(aatype_CHnucleiType_TresonIndex_mdict[aa_type].keys()))
    for CHnucleiType, TresonIndex  in list(aatype_CHnucleiType_TresonIndex_mdict[aa_type].items()):
        ResonancesTuple = TOCSYindex_ResonancesTuple_dict[TresonIndex]
        #print "DEBUG: aa_type=", aa_type, "CHnucleiType=", CHnucleiType, "TresonIndex=", TresonIndex, "ResonancesTuple=", ResonancesTuple
        Carbon_name = CHnucleiType.split('_')[0]
        Hydrogen_name = CHnucleiType.split('_')[1]
        # Rename un-resolved Geminal protons to QA, QB, QG, QD, QE.
        if len(aa_carbonBondedHydrogensDict_dict[aa_type][Carbon_name]) > 1:    # if this carbon has geminal protons
            print("DEBUG: aa_type=", aa_type, "Carbon_name=", Carbon_name)
            print("DEBUG: aatype_CHnucleiType_TresonIndex_mdict[aa_type].keys()=", list(aatype_CHnucleiType_TresonIndex_mdict[aa_type].keys()))
            print("DEBUG: aa_carbonBondedHydrogensDict_dict[aa_type][Carbon_name]=", aa_carbonBondedHydrogensDict_dict[aa_type][Carbon_name])
            if [Carbon_name+"_" in x for x in list(aatype_CHnucleiType_TresonIndex_mdict[aa_type].keys())].count(True) < len(aa_carbonBondedHydrogensDict_dict[aa_type][Carbon_name]):
                print("DEBUG entered condition 1!")
                if aa_type in list(aatype_carbon_methylHydrogens_mdict.keys()) and Carbon_name in list(aatype_carbon_methylHydrogens_mdict[aa_type].keys()):
                    print("DEBUG entered condition 2!")
                    new_Hydrogen_name = aatype_carbon_methylHydrogens_mdict[aa_type][Carbon_name]
                    Hydrogen_name = new_Hydrogen_name
                #else:
                #    print "DEBUG entered condition 3!"
                #    new_Hydrogen_name = 'Q'+Hydrogen_name[1]
                #    Hydrogen_name = new_Hydrogen_name
        # CONDITIONS:
        # 1) 
        if aa_type in list(aatype_carbon_methylHydrogens_mdict.keys()) and Carbon_name in list(aatype_carbon_methylHydrogens_mdict[aa_type].keys()) and (not aa_type in list(aa_CarbonListwithDegenerateH_dict.keys()) or Carbon_name in aa_CarbonListwithDegenerateH_dict[aa_type]):
            print("DEBUG entered condition 4! aa_type = ", aa_type, "Carbon_name=", Carbon_name)
            new_Hydrogen_name = aatype_carbon_methylHydrogens_mdict[aa_type][Carbon_name]
            Hydrogen_name = new_Hydrogen_name
        Carbon_resonance = ResonancesTuple[1]
        Hydrogen_resonance = ResonancesTuple[0]
        Nreson = ResonancesTuple[2]
        HNreson = ResonancesTuple[3]
        HNresons_list.append(float(ResonancesTuple[3]))
        Nresons_list.append(float(ResonancesTuple[2]))
        
        #print "\t"+str(atom_index)+"\t"+str(Carbon_resonance)+"\t0.2\t"+str(Carbon_name)+"\t"+str(residue[3:])
        Creson_Cname_lines_list.append([float(Carbon_resonance), Carbon_name])  
            
        #print "\t"+str(atom_index)+"\t"+str(Hydrogen_resonance)+"\t0.02\t"+str(Hydrogen_name)+"\t"+str(residue[3:])
        #xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s\n" % (atom_index, float(Hydrogen_resonance), 0.02, Hydrogen_name, int(residue[3:]), aa_type))
        xeasy_lines_list.append( [atom_index, float(Hydrogen_resonance), 0.02, Hydrogen_name, int(residue[3:]), aa_type] )  # append the Hydrogen resonance line
        atom_index += 1
        try:
            #sparky_fout.write("\t%s-%s-%s-H\t\t%.3f\t%.3f\t%.3f\t%.3f\t\n" % ( aa3to1_dict[aa_type] + residue[3:] + Hydrogen_name, Carbon_name, protein_alignment_list[residue2position_dict[residue]] + str(int(residue[3:])+1) + "N",float(Hydrogen_resonance), float(Carbon_resonance), float(Nreson), float(HNreson)) )
            sparky_lines_list.append( [aa3to1_dict[aa_type] + residue[3:] + Hydrogen_name, Carbon_name, protein_alignment_list[residue2position_dict[residue]] + str(int(residue[3:])+1) + "N",float(Hydrogen_resonance), float(Carbon_resonance), float(Nreson), float(HNreson)] )
        except IndexError:
            #sparky_fout.write("\t%s-%s-%s-H\t\t%.3f\t%.3f\t%.3f\t%.3f\t\n" % ( aa3to1_dict[aa_type] + residue[3:] + Hydrogen_name, Carbon_name, str(int(residue[3:])+1) + "N", float(Hydrogen_resonance), float(Carbon_resonance), float(Nreson), float(HNreson)) )
            sparky_lines_list.append( [aa3to1_dict[aa_type] + residue[3:] + Hydrogen_name, Carbon_name, str(int(residue[3:])+1) + "N", float(Hydrogen_resonance), float(Carbon_resonance), float(Nreson), float(HNreson)] )
        
    
    Cnames_list = [cl[1] for cl in Creson_Cname_lines_list]
    Cnames_set = set(Cnames_list)
    for Carbon_name in Cnames_set:
        if Cnames_list.count(Carbon_name) > 1:  # if this carbon occurs in multiple C-H pairs, calculate its average resonance for the XEASY format
            Carbon_resonance = np.average([cl[0] for cl in Creson_Cname_lines_list if cl[1]==Carbon_name])
        else:
            Carbon_resonance = np.average([cl[0] for cl in Creson_Cname_lines_list if cl[1]==Carbon_name])
        #xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s\n" % (atom_index, float(Carbon_resonance), 0.2, Carbon_name, int(residue[3:]), aa_type))
        xeasy_lines_list.append( [atom_index, float(Carbon_resonance), 0.2, Carbon_name, int(residue[3:]), aa_type] )
        atom_index += 1
    
    average_HNreson = np.mean(HNresons_list)
    stdev_HNreson = np.std(HNresons_list)
    average_Nreson = np.mean(Nresons_list)
    stdev_Nreson = np.std(Nresons_list)
    
    # RENAME DEGENERATE PROTONS
    new_xeasy_lines_list, new_sparky_lines_list = rename_protons(xeasy_lines_list, sparky_lines_list)
    # WRITE C AND H LINES IN XEASY FORMAT
    for xeasy_line in new_xeasy_lines_list:
        xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s\n" % (xeasy_line[0], xeasy_line[1], xeasy_line[2], xeasy_line[3], xeasy_line[4], xeasy_line[5]) )
    # WRITE C AND H LINES IN SPARKY FORMAT
    for sparky_line in new_sparky_lines_list:
        sparky_fout.write("\t%s-%s-%s-H\t\t%.3f\t%.3f\t%.3f\t%.3f\t\n" % (sparky_line[0], sparky_line[1], sparky_line[2], sparky_line[3], sparky_line[4], sparky_line[5], sparky_line[6]))
    
    #print "\t"+str(atom_index)+"\t"+str(float(average_HNreson))+"\t0.02\tH\t"+str(int(residue[3:])+1)
    
    # FINALLY WRITE N AND H LINES IN XEASY FORMAT
    try:
        xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s\n" % (atom_index, float(average_HNreson), 0.02, "H", int(residue[3:])+1, aa1to3_dict[protein_alignment_list[residue2position_dict[residue]]])) # The N-HN of i sees the C-H of i-1
        atom_index += 1
        #print "\t"+str(atom_index)+"\t"+str(float(average_Nreson))+"\t0.2\tN\t"+str(int(residue[3:])+1)
        xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s\n" % (atom_index, average_Nreson, stdev_Nreson, "N", int(residue[3:])+1, aa1to3_dict[protein_alignment_list[residue2position_dict[residue]]]))   # The N-HN of i sees the C-H of i-1
        atom_index += 1
    except IndexError:
        xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s\n" % (atom_index, float(average_HNreson), 0.02, "H", int(residue[3:])+1, aa1to3_dict[protein_alignment_list[residue2position_dict[residue]]])) # The N-HN of i sees the C-H of i-1
        atom_index += 1
        #print "\t"+str(atom_index)+"\t"+str(float(average_Nreson))+"\t0.2\tN\t"+str(int(residue[3:])+1)
        xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s\n" % (atom_index, average_Nreson, stdev_Nreson, "N", int(residue[3:])+1, aa1to3_dict[protein_alignment_list[residue2position_dict[residue]]]))   # The N-HN of i sees the C-H of i-1
        atom_index += 1
    
    # IN ADDITION WRITE N AND H LINES IN HSQC FILE
    


## NOW PRINT THE N AND HN RESONANCES OF THE PATCHED RESIDUES, IF APPLICABLE
residue_residueTypesProbTupleList_dict = OrderedDict() # ordereddict with keys the AAIG of residue i and values lists of tuples of the form (residue i-1 matching aa type, average probability)
Num_of_NOESY_resonances = 0 # save here the number of TOCSY resonances for a particular AAIG
residue = None
possible_aatype_prob_C_H_resonpair_NOESYindex_list_list = [] # list of lists of the form (aa type, Carbon name, Hydrogen name, weighted average probability, NOESY_reson_index, H resonance, C resonance, carbon group)
previous_aaindex = ""
residue_list = list(NOESY_residue_CSlist_dict.keys())
carbon_groups_list = []  # list of lists of the form [H resonance, C resonance, N resonance, HN resonance, carbon group]
for res_index, residue in enumerate(residue_list):
    #if residue == "GLN136":
    #try:
        print("Assigning chemical shifts to residue ", residue)
        aa_type = residue[0:3]  # aa type in 3-letter code
        #print "DEBUG: NOESY_residue_CSlist_dict[residue] =", NOESY_residue_CSlist_dict[residue]
        NOESYindex_ResonancesTuple_dict = {}
        for NOESY_Hreson,NOESY_Creson,NOESY_Nreson,NOESY_HNreson in NOESY_residue_CSlist_dict[residue]:
            Num_of_NOESY_resonances += 1
            #print "DEBUG: NOESY_Creson=",NOESY_Creson,"NOESY_Hreson=",NOESY_Hreson
            NOESYindex_ResonancesTuple_dict[Num_of_NOESY_resonances] = (NOESY_Hreson, NOESY_Creson, NOESY_Nreson, NOESY_HNreson)
        Num_of_NOESY_resonances = 0
    
    #except (ValueError, IndexError):
    #    print "WARNING: the 3rd and 4th elements of the following NOESY file line are not numbers:"
    #    print "NOESY_Hreson, NOESY_Creson:", NOESY_Hreson, NOESY_Creson
    #    sys.exit(1)
        
        # Write the Chemical Shifts Assignments in XEASY format for Rosetta
        HNresons_list, Nresons_list = [], []
        for TresonIndex  in list(NOESYindex_ResonancesTuple_dict.keys()):
            ResonancesTuple = NOESYindex_ResonancesTuple_dict[TresonIndex]
            Nreson = ResonancesTuple[2]
            HNreson = ResonancesTuple[3]
            HNresons_list.append(float(ResonancesTuple[3]))
            Nresons_list.append(float(ResonancesTuple[2]))
        
        average_HNreson = np.mean(HNresons_list)
        stdev_HNreson = np.std(HNresons_list)
        average_Nreson = np.mean(Nresons_list)
        stdev_Nreson = np.std(Nresons_list)
        
        # FINALLY WRITE N AND H LINES IN XEASY FORMAT
        try:
            xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s\n" % (atom_index, float(average_HNreson), 0.02, "H", int(residue[3:])+1, aa1to3_dict[protein_alignment_list[residue2position_dict[residue]]])) # The N-HN of i sees the C-H of i-1
            atom_index += 1
            print("\t"+str(atom_index)+"\t"+str(float(average_HNreson))+"\t0.02\tH\t"+str(int(residue[3:])+1))
            print("\t"+str(atom_index)+"\t"+str(float(average_Nreson))+"\t0.2\tN\t"+str(int(residue[3:])+1))
            xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s\n" % (atom_index, average_Nreson, stdev_Nreson, "N", int(residue[3:])+1, aa1to3_dict[protein_alignment_list[residue2position_dict[residue]]]))   # The N-HN of i sees the C-H of i-1
            atom_index += 1
        except IndexError:
            xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s\n" % (atom_index, float(average_HNreson), 0.02, "H", int(residue[3:])+1, aa1to3_dict[protein_alignment_list[residue2position_dict[residue]]])) # The N-HN of i sees the C-H of i-1
            atom_index += 1
            print("\t"+str(atom_index)+"\t"+str(float(average_HNreson))+"\t0.02\tH\t"+str(int(residue[3:])+1))
            print("\t"+str(atom_index)+"\t"+str(float(average_Nreson))+"\t0.2\tN\t"+str(int(residue[3:])+1))
            xeasy_fout.write("\t%g\t%.3f\t%.3f\t%s\t%g\t# %s\n" % (atom_index, average_Nreson, stdev_Nreson, "N", int(residue[3:])+1, aa1to3_dict[protein_alignment_list[residue2position_dict[residue]]]))   # The N-HN of i sees the C-H of i-1
            atom_index += 1

xeasy_fout.close()
sparky_fout.close()

# Sort xeasy to make sure all the lines for each residue are together
with open(xeasy_fout_name, 'r') as fin:
    contents = []
    for line in fin:
        word_list = line.split()
        word_list[4] = int(word_list[4])
        contents.append(word_list)
    contents.sort(key=itemgetter(4,3))
    with open('tmp.xeasy', 'w') as fout:
        aindex = 1
        for line in contents:
            line[0] = str(aindex)
            fout.write("\t%s\t%s\t%s\t%s\t%g\t%s %s\n" % tuple(line))
            aindex += 1
os.remove(xeasy_fout_name)
os.rename('tmp.xeasy', xeasy_fout_name)

## Append the unlabeled lines to Sparky file
# READ TOCSY FILE
def read_whole_spectrum_file(query_fname):
    """ Function to load all TOCSY file contents, including lines ?-?-?-?. """
    with open(query_fname, 'r') as f:
        tmp_query_contents=f.readlines()    # contents of original query_fname (4D TOCSY or 4D NOESY) in 5 column format (name H C N HN)
    query_contents=[]   # list with the line_lists of TOCSY file
    for line in tmp_query_contents: # CLEANING THE SPECTRUM FROM IRRELEVANT LINES
        word_list = line.split()
        try:
            float(word_list[1]); float(word_list[2]); float(word_list[3]); float(word_list[4]); # checking if the line contains numbers (resonances)
            #print "DEBUG: appending line:", line
            query_contents.append(line)
        except (IndexError, ValueError):
            print("WARNING: Discarding TOCSY line:", line)
    
    # remove duplicate lines from query_fname (4D TOCSY or 4D NOESY)
    lines2remove_set = set()
    for qline1 in query_contents:
        try:
            counter = 0
            query1_words_list = qline1.split()
            q1_w1=float(query1_words_list[1])
            q1_w2=float(query1_words_list[2])
            q1_w3=float(query1_words_list[3])
            q1_w4=float(query1_words_list[4])
            for qline2 in query_contents:
                query2_words_list = qline2.split()
                q2_w1=float(query2_words_list[1])
                q2_w2=float(query2_words_list[2])
                q2_w3=float(query2_words_list[3])
                q2_w4=float(query2_words_list[4])
                if approx_equal_proportional(q1_w1, q2_w1) and approx_equal_proportional(q1_w2, q2_w2) and approx_equal_proportional(q1_w3, q2_w3) and approx_equal_proportional(q1_w4, q2_w4):
                    counter += 1
                    if counter > 1:
                        lines2remove_set.add(qline2)
                        #print "DEBUG: will remove line",qline2
        except (ValueError, IndexError):
            #print "WARNING: the 2nd and 3rd elements of the following root file line are not numbers:"
            #print "Root file line:", root_line
            continue
    # now remove the duplicate lines
    for qline in lines2remove_set:
        query_contents.remove(qline)
    
    return query_contents


TOCSY_line_list = read_whole_spectrum_file(args.TOCSY_fname)    # read the original TOCSY, including the "?-?-?-?" lines
if args.OUT_fname:
    sparky_fout = open(args.OUT_fname+".sparky", 'a+')
else:
    sparky_fout = open(args.NHmap_FILE+".sparky", 'a+')   # append to the annotated Sparky the unlabeled lines
sparky_lines = sparky_fout.readlines()
for TOCSY_line in TOCSY_line_list:
    TOCSY_lineList = TOCSY_line.split()
    # print "DEBUG: TOCSY_lineList=", TOCSY_lineList
    FOUND = False
    TOCSY_label = TOCSY_lineList[0]
    TOCSY_Hreson = float(TOCSY_lineList[1])
    TOCSY_Creson = float(TOCSY_lineList[2])
    TOCSY_Nreson = float(TOCSY_lineList[3])
    TOCSY_HNreson = float(TOCSY_lineList[4])
    for line in sparky_lines:
        # try:
        # print "DEBUG: line=", line
        word_list = line.split()
        Hreson = float(word_list[1])
        Creson = float(word_list[2])
        Nreson = float(word_list[3])
        HNreson = float(word_list[4])
        if Hreson == TOCSY_Hreson and Creson == TOCSY_Creson and Nreson == TOCSY_Nreson and HNreson == TOCSY_HNreson:
            # print "DEBUG: matching lines:"
            print(line, end=' ')
            print(TOCSY_line)
            FOUND = True
            break
        # except ValueError:  # if this sparky line does not comply  with the format (e.g. "# PATCHED RESIDUES: A213 R304 S345 G360"), skip it
        #     continue
    
    if FOUND == False:  # if this line hasn't been write into the annotated TOCSY Sparky file, append it
        sparky_fout.write("\t%s\t\t%.3f\t%.3f\t%.3f\t%.3f\t\n" % (TOCSY_label, TOCSY_Hreson, TOCSY_Creson, TOCSY_Nreson, TOCSY_HNreson))


# BEFORE YOU CLOSE SPARKY FILE APPEND AT THE END AS A COMMENT THE PATCHED N-TERMINAL RESIDUES
actual_patched_residues_list = []   # patched_residues_list has the i-1 residues, actual_patched_residues_list has the i. I.e. N212 and A213,
                                    # only the A213 is the actually patched residue, N212 is not in the alignment
for residue in patched_residues_list:
    resid = int(residue[3:])    # recall that we made residue variable by concatenating the aa type (3-letter) and the resid
    actual_patched_residues_list.append(protein_alignment_list[resid+1 -args.FIRST_RESIDUE_NUMBER] + str(resid+1))
sparky_fout.write("# PATCHED RESIDUES: " + " ".join(actual_patched_residues_list) + "\n")
sparky_fout.close()

