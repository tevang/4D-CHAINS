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


#!/usr/bin/env python2.7

import sys, re, os, csv
from operator import itemgetter
from ordereddict import OrderedDict
from argparse import ArgumentParser
from tabulate import tabulate
from ete3 import Tree
from scipy import stats
import gc
import collections
import numpy as np
from cluster import HierarchicalClustering
from scipy.stats.mstats import zscore
from scipy.stats import norm

def tree(): # function to create multidimensional dictionaries
    return collections.defaultdict(tree)

HOME_DIR = os.path.dirname(os.path.realpath(__file__))
code3 = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
code1 = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
aa1to3_dict = dict((c1,c3) for c1,c3 in zip(code1,code3))
aa3to1_dict = dict((c3,c1) for c3,c1 in zip(code3,code1))

allowed_aa_atoms_dict = {
"ALA" : ["HA", "HB", "CA", "CB", "N", "H"],
"ARG" : ["HA", "HB2", "HB3", "HG2", "HG3", "HD2", "HD3", "CA", "CB", "CG", "CD", "N", "H"],
"ASP" : ["HA", "HB2", "HB3", "CA", "CB", "N", "H"],
"ASN" : ["HA", "HB2", "HB3", "CA", "CB", "N", "H"],
"CYS" : ["HA", "HB2", "HB3", "CA", "CB", "N", "H"],
"GLU" : ["HA", "HB2", "HB3", "HG2", "HG3", "CA", "CB", "CG", "N", "H"],
"GLN" : ["HA", "HB2", "HB3", "HG2", "HG3", "CA", "CB", "CG", "N", "H"],
"GLY" : ["HA2", "HA3", "CA", "N", "H"],
"HIS" : ["HA", "HB2", "HB3", "CA", "CB", "N", "H"],
"ILE" : ["HA", "HB", "HG12", "HG13", "HG2", "HD1", "CA", "CB", "CG1", "CG2", "CD1", "N", "H"],
"LEU" : ["HA", "HB2", "HB3", "HG", "HD1", "HD2", "CA", "CB", "CG", "CD1", "CD2", "N", "H"],
"LYS" : ["HA", "HB2", "HB3", "HG2", "HG3", "HD2", "HD3", "HE2", "HE3", "CA", "CB", "CG", "CD", "CE", "N", "H"],
"MET" : ["HA", "HB2", "HB3", "HG2", "HG3", "HE2", "HE3", "CA", "CB", "CG", "N", "H", "CE"],
"PHE" : ["HA", "HB2", "HB3", "CA", "CB", "N", "H"],
"PRO" : ["HA", "HB2", "HB3", "HG2", "HG3", "HD2", "HD3", "CA", "CB", "CG", "CD", "N"], # Prolines are not detected by the method at position "i" due to lack of HN hydrogen
"SER" : ["HA", "HB2", "HB3", "CA", "CB", "N", "H"],
"THR" : ["HA", "HB", "HG2", "CA", "CB", "CG2", "N", "H"],
"TRP" : ["HA", "HB2", "HB3", "CA", "CB", "N", "H"],
"TYR" : ["HA", "HB2", "HB3", "HD1", "HD2", "HE1", "HE2", "CA", "CB", "CD1", "CD2", "CE1", "CE2", "N", "H"],
"VAL" : ["HA", "HB", "HG1", "HG2", "CA", "CB", "CG1", "CG2", "N", "H"]
}

aa_alternativeAtomNames_dict = tree()   # the first element in the list (dict value) is always the conventional atom Name
aa_alternativeAtomNames_dict["ILE"]["CG2"] = ["HG2", "HG21", "HG22", "HG23"]
aa_alternativeAtomNames_dict["ILE"]["CD1"] = ["HD1", "HD11", "HD12", "HD13"]
aa_alternativeAtomNames_dict["LEU"]["CD2"] = ["HD2", "HD21", "HD22", "HD23"]
aa_alternativeAtomNames_dict["LEU"]["CD1"] = ["HD1", "HD11", "HD12", "HD13"]
aa_alternativeAtomNames_dict["VAL"]["CG2"] = ["HG2", "HG21", "HG22", "HG23"]
aa_alternativeAtomNames_dict["VAL"]["CG1"] = ["HG1", "HG11", "HG12", "HG13"]
aa_alternativeAtomNames_dict["THR"]["CG2"] = ["HG2", "HG21", "HG22", "HG23"]
aa_alternativeAtomNames_dict["ALA"]["CB"] = ["HB", "HB1", "HB2", "HB3"]
aa_alternativeAtomNames_dict["MET"]["CE"] = ["HE", "HE1", "HE2", "HE3"]


aa_CHpair_2Dhist_multidict = tree()
aa_CHpair_2Dhist_multidict["ALA"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_multidict["ALA"]["CB-HB"] = [[], []]

aa_CHpair_2Dhist_multidict["ARG"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_multidict["ARG"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_multidict["ARG"]["CB-HB3"] = [[], []]
aa_CHpair_2Dhist_multidict["ARG"]["CG-HG2"] = [[], []]
aa_CHpair_2Dhist_multidict["ARG"]["CG-HG3"] = [[], []]
aa_CHpair_2Dhist_multidict["ARG"]["CD-HD2"] = [[], []]
aa_CHpair_2Dhist_multidict["ARG"]["CD-HD3"] = [[], []]

aa_CHpair_2Dhist_multidict["ASP"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_multidict["ASP"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_multidict["ASP"]["CB-HB3"] = [[], []]

aa_CHpair_2Dhist_multidict["ASN"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_multidict["ASN"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_multidict["ASN"]["CB-HB3"] = [[], []]

aa_CHpair_2Dhist_multidict["CYS"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_multidict["CYS"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_multidict["CYS"]["CB-HB3"] = [[], []]

aa_CHpair_2Dhist_multidict["GLU"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_multidict["GLU"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_multidict["GLU"]["CB-HB3"] = [[], []]
aa_CHpair_2Dhist_multidict["GLU"]["CG-HG2"] = [[], []]
aa_CHpair_2Dhist_multidict["GLU"]["CG-HG3"] = [[], []]

aa_CHpair_2Dhist_multidict["GLN"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_multidict["GLN"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_multidict["GLN"]["CB-HB3"] = [[], []]
aa_CHpair_2Dhist_multidict["GLN"]["CG-HG2"] = [[], []]
aa_CHpair_2Dhist_multidict["GLN"]["CG-HG3"] = [[], []]

aa_CHpair_2Dhist_multidict["GLY"]["CA-HA2"] = [[], []]
aa_CHpair_2Dhist_multidict["GLY"]["CA-HA3"] = [[], []]

aa_CHpair_2Dhist_multidict["HIS"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_multidict["HIS"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_multidict["HIS"]["CB-HB3"] = [[], []]

aa_CHpair_2Dhist_multidict["ILE"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_multidict["ILE"]["CB-HB"] = [[], []]
aa_CHpair_2Dhist_multidict["ILE"]["CG1-HG12"] = [[], []]
aa_CHpair_2Dhist_multidict["ILE"]["CG1-HG13"] = [[], []]
aa_CHpair_2Dhist_multidict["ILE"]["CG2-HG2"] = [[], []]
aa_CHpair_2Dhist_multidict["ILE"]["CD1-HD1"] = [[], []]

aa_CHpair_2Dhist_multidict["LEU"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_multidict["LEU"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_multidict["LEU"]["CB-HB3"] = [[], []]
aa_CHpair_2Dhist_multidict["LEU"]["CG-HG"] = [[], []]
aa_CHpair_2Dhist_multidict["LEU"]["CD1-HD1"] = [[], []]
aa_CHpair_2Dhist_multidict["LEU"]["CD2-HD2"] = [[], []]

aa_CHpair_2Dhist_multidict["LYS"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_multidict["LYS"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_multidict["LYS"]["CB-HB3"] = [[], []]
aa_CHpair_2Dhist_multidict["LYS"]["CG-HG2"] = [[], []]
aa_CHpair_2Dhist_multidict["LYS"]["CG-HG3"] = [[], []]
aa_CHpair_2Dhist_multidict["LYS"]["CD-HD2"] = [[], []]
aa_CHpair_2Dhist_multidict["LYS"]["CD-HD3"] = [[], []]
aa_CHpair_2Dhist_multidict["LYS"]["CE-HE2"] = [[], []]
aa_CHpair_2Dhist_multidict["LYS"]["CE-HE3"] = [[], []]

aa_CHpair_2Dhist_multidict["MET"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_multidict["MET"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_multidict["MET"]["CB-HB3"] = [[], []]
aa_CHpair_2Dhist_multidict["MET"]["CG-HG2"] = [[], []]
aa_CHpair_2Dhist_multidict["MET"]["CG-HG3"] = [[], []]
aa_CHpair_2Dhist_multidict["MET"]["CE-HE"] = [[], []]

aa_CHpair_2Dhist_multidict["PHE"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_multidict["PHE"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_multidict["PHE"]["CB-HB3"] = [[], []]
aa_CHpair_2Dhist_multidict["PHE"]["CD1-HD1"] = [[], []]
aa_CHpair_2Dhist_multidict["PHE"]["CD2-HD2"] = [[], []]
aa_CHpair_2Dhist_multidict["PHE"]["CE1-HE1"] = [[], []]
aa_CHpair_2Dhist_multidict["PHE"]["CE2-HE2"] = [[], []]
aa_CHpair_2Dhist_multidict["PHE"]["CZ-HZ"] = [[], []]

aa_CHpair_2Dhist_multidict["PRO"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_multidict["PRO"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_multidict["PRO"]["CB-HB3"] = [[], []]
aa_CHpair_2Dhist_multidict["PRO"]["CG-HG2"] = [[], []]
aa_CHpair_2Dhist_multidict["PRO"]["CG-HG3"] = [[], []]
aa_CHpair_2Dhist_multidict["PRO"]["CD-HD2"] = [[], []]
aa_CHpair_2Dhist_multidict["PRO"]["CD-HD3"] = [[], []]

aa_CHpair_2Dhist_multidict["SER"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_multidict["SER"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_multidict["SER"]["CB-HB3"] = [[], []]

aa_CHpair_2Dhist_multidict["THR"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_multidict["THR"]["CB-HB"] = [[], []]
aa_CHpair_2Dhist_multidict["THR"]["CG2-HG2"] = [[], []]

aa_CHpair_2Dhist_multidict["TRP"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_multidict["TRP"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_multidict["TRP"]["CB-HB3"] = [[], []]

aa_CHpair_2Dhist_multidict["TYR"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_multidict["TYR"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_multidict["TYR"]["CB-HB3"] = [[], []]
aa_CHpair_2Dhist_multidict["TYR"]["CD1-HD1"] = [[], []]
aa_CHpair_2Dhist_multidict["TYR"]["CD2-HD2"] = [[], []]
aa_CHpair_2Dhist_multidict["TYR"]["CE1-HE1"] = [[], []]
aa_CHpair_2Dhist_multidict["TYR"]["CE2-HE2"] = [[], []]

aa_CHpair_2Dhist_multidict["VAL"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_multidict["VAL"]["CB-HB"] = [[], []]
aa_CHpair_2Dhist_multidict["VAL"]["CG1-HG1"] = [[], []]
aa_CHpair_2Dhist_multidict["VAL"]["CG2-HG2"] = [[], []]


aa_CHpairs2beMerged_dict = {}
aa_CHpairs2beMerged_dict["ALA"] = []
aa_CHpairs2beMerged_dict["ARG"] = [
    ["CB-HB2", "CB-HB3"],
    ["CG-HG2", "CG-HG3"],
    ["CD-HD2", "CD-HD3"]
]
aa_CHpairs2beMerged_dict["ASP"] = [
    ["CB-HB2", "CB-HB3"]
]
aa_CHpairs2beMerged_dict["ASN"] = [
    ["CB-HB2", "CB-HB3"]
]
aa_CHpairs2beMerged_dict["CYS"] = [
    ["CB-HB2", "CB-HB3"]
]
aa_CHpairs2beMerged_dict["GLU"] = [
    ["CB-HB2", "CB-HB3"],
    ["CG-HG2", "CG-HG3"]
]
aa_CHpairs2beMerged_dict["GLN"] = [
    ["CB-HB2", "CB-HB3"],
    ["CG-HG2", "CG-HG3"]
]
aa_CHpairs2beMerged_dict["GLY"] = [
    ["CA-HA2", "CA-HA3"]
]
aa_CHpairs2beMerged_dict["HIS"] = [
    ["CB-HB2", "CB-HB3"]
]
aa_CHpairs2beMerged_dict["ILE"] = [
    ["CG1-HG12", "CG1-HG13"],
]
aa_CHpairs2beMerged_dict["LEU"] = [
    ["CB-HB2", "CB-HB3"],
    ["CD1-HD1", "CD2-HD2"]
]
aa_CHpairs2beMerged_dict["LYS"] = [
    ["CB-HB2", "CB-HB3"],
    ["CG-HG2", "CG-HG3"],
    ["CD-HD2", "CD-HD3"],
    ["CE-HE2", "CE-HE3"]
]
aa_CHpairs2beMerged_dict["MET"] = [
    ["CB-HB2", "CB-HB3"],
    ["CG-HG2", "CG-HG3"]
]
aa_CHpairs2beMerged_dict["PHE"] = [
    ["CB-HB2", "CB-HB3"]
]
aa_CHpairs2beMerged_dict["PRO"] = [
    ["CB-HB2", "CB-HB3"],
    ["CG-HG2", "CG-HG3"],
    ["CD-HD2", "CD-HD3"]
]
aa_CHpairs2beMerged_dict["SER"] = [
    ["CB-HB2", "CB-HB3"]
]
aa_CHpairs2beMerged_dict["THR"] = []
aa_CHpairs2beMerged_dict["TRP"] = [
    ["CB-HB2", "CB-HB3"]
]
aa_CHpairs2beMerged_dict["TYR"] = [
    ["CB-HB2", "CB-HB3"]
]
aa_CHpairs2beMerged_dict["VAL"] = [
    ["CG1-HG1", "CG2-HG2"]
]


CO_CA_CB_aa_SecStruct_list = []
CB_CG_aa_SecStruct_list = []
CG_CD_aa_SecStruct_list = []
CO_CA_CB_SS_prevAA_curAA_nextAA_list = []

VASCO_naming_exceptions_multidict = tree()
VASCO_naming_exceptions_multidict["ALA"]["HB*"] = "HB"
VASCO_naming_exceptions_multidict["ILE"]["HD1*"] = "HD1"
VASCO_naming_exceptions_multidict["ILE"]["HG2*"] = "HG2"
VASCO_naming_exceptions_multidict["LEU"]["HD1*"] = "HD1"
VASCO_naming_exceptions_multidict["LEU"]["HD2*"] = "HD2"
VASCO_naming_exceptions_multidict["THR"]["HG2*"] = "HG2"
VASCO_naming_exceptions_multidict["VAL"]["HG1*"] = "HG1"
VASCO_naming_exceptions_multidict["VAL"]["HG2*"] = "HG2"
VASCO_naming_exceptions_multidict["MET"]["HE*"] = "HE2"

def cmdlineparse():
    parser = ArgumentParser(description="command line arguments",
        epilog="EXAMPLE: ./create_2D_histograms.py -bmrbfile BMRB_data/Atom_chem_shift.proteins_with_sidechains.csv -selfile BMRB_data/Rescue2_BMRB_selected_sequences.txt")
    parser.add_argument("-bmrbfile", dest="BMRB_ENTRIES_FILE", required=False, type=str, help="csv file with all the chemical shifts deposited in BMRB",
                        metavar="<BMRB entries file>")
    parser.add_argument("-vascodir", dest="VASCO_DIR", required=False, type=str, help="folder with all the VASCO database files",
                        metavar="<VASCO database directory>")
    parser.add_argument("-selfile", dest="SELECTED_ENTRIES_FILE", required=False, type=str, 
                        help="file with a list of BMRB entries to be used to make 2D histograms", metavar="<selected BMRB entries files>")
    parser.add_argument("-wH", dest="H_weight", required=False, type=float, default=0.1, help="weight in (0.0-1.0] to apply on aa type prediction from aliphatic H resonances (default: 0.1)", metavar="<H weight>")
    parser.add_argument("-wC", dest="C_weight", required=False, type=float, default=1.0, help="weight in (0.0, 1.0] to apply on aa type prediction from aliphatic C resonances (default: 1.0)", metavar="<C weight>")
    parser.add_argument("-resH", dest="H_bin_length", required=False, type=float, default=0.04,
                        help="histogram bin length for the carbon in ppm", metavar="<H bin length>")
    parser.add_argument("-resC", dest="C_bin_length", required=False, type=float, default=0.2,
                        help="histogram bin length for the hydrogen in ppm", metavar="<C bin length>")
    parser.add_argument("-nomerge", dest="MERGE_EQUIVALENT_CH_PAIRS", required=False, default=True, action='store_false', help="merge equivalent C-H pairs")
    args=parser.parse_args()
    return args


args = cmdlineparse()


def drange(start, stop, step):
    r = start
    while r < stop:
        yield r
        r += step

i0 = drange(-9.98, 34, args.H_bin_length)
Hedges_list = [round(x,2) for x in i0]
i0 = drange(0.1, 200, args.C_bin_length)
Cedges_list = [round(x,1) for x in i0]


def p_adjust_bh(p):
    """Benjamini-Hochberg p-value correction for multiple hypothesis testing."""
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = len(p) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]


def discard_outliers(Hshifts, Cshifts, H_thresh=15.0, C_thresh=10.0):
    """
    Returns a boolean array with True if points are outliers and False 
    otherwise.

    Parameters:
    -----------
        points : An numobservations by numdimensions array of observations
        thresh : The modified z-score to use as a threshold. Observations with
            a modified z-score (based on the median absolute deviation) greater
            than this value will be classified as outliers.

    Returns:
    --------
        mask : A numobservations-length boolean array.

    References:
    ----------
        Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
        Handle Outliers", The ASQC Basic References in Quality Control:
        Statistical Techniques, Edward F. Mykytka, Ph.D., Editor. 
    """
    
    
    Hzscores = zscore(Hshifts)
    Hp_values = norm.sf(abs(Hzscores))*2
    adj_Hp_values = p_adjust_bh(Hp_values)
    final_discard_H = adj_Hp_values < 10e-100
    
    Czscores = zscore(Cshifts)
    Cp_values = norm.sf(abs(Czscores))*2
    adj_Cp_values = p_adjust_bh(Cp_values)
    final_discard_C = adj_Cp_values < 10e-10
    
    if len(Cshifts.shape) == 1:
        Cshifts = Cshifts[:,None]
    if len(Hshifts.shape) == 1:
        Hshifts = Hshifts[:,None]
    
    outliers = final_discard_C | final_discard_H    # boolean array, True if this value either in Cshifts or Hshifts is an outlier
    i=0
    new_Cshifts_list, new_Hshifts_list = [], []
    for C, H in zip(Cshifts, Hshifts):
        if outliers[i] == False:
            new_Cshifts_list.append(C[0])
            new_Hshifts_list.append(H[0])
        i += 1

    return np.array(new_Hshifts_list), np.array(new_Cshifts_list)


def read_BMRB_file():
    global args
    
    print "Loading BMRB file with chemical shifts ..."
    reader=csv.reader(open(args.BMRB_ENTRIES_FILE,"rb"),delimiter=',')
    line_list=list(reader)

    if args.SELECTED_ENTRIES_FILE:
        selected_entries_list = []
        with open(args.SELECTED_ENTRIES_FILE) as f:
            for entry in f:
                selected_entries_list.append(entry.strip())
        selected_lines_list = [l for l in line_list if l[-2] in selected_entries_list]
        return selected_lines_list
    
    return line_list


def parse_aromatic_resonances(resname, atomname_shift_dict):
    """
    FUNCTION to modify the atomname_shift_dict by removing TYR CD-QD, CE-QE and PHE CD-QD that are not the same in both equivalent
    carbons and protons.
    """
    
    new_atomname_shift_dict = {}
    if not resname in ["TYR", "PHE"]:   # if it's not aromatic, no parsing is needed
        return atomname_shift_dict
    elif resname == "TYR":
        CD_shifts = [shift for atomname, shift in atomname_shift_dict.items() if atomname in ["CD1", "CD2"] ]
        HD_shifts = [shift for atomname, shift in atomname_shift_dict.items() if atomname in ["HD1", "HD2"] ]
        CE_shifts = [shift for atomname, shift in atomname_shift_dict.items() if atomname in ["CE1", "CE2"] ]
        HE_shifts = [shift for atomname, shift in atomname_shift_dict.items() if atomname in ["HE1", "HE2"] ]
        CD_reson, QD_reson, CE_reson, QE_reson = None, None, None, None
        if len(set(CD_shifts))==1 and len(set(HD_shifts))==1:     # if CD1 and CD2 or HD1 and HD2 are not different, save them
            CD_reson = CD_shifts[0]
            QD_reson = HD_shifts[0]
        if len(set(CE_shifts))==1 and len(set(HE_shifts))==1:     # if CE1 and CE2 or HE1 and HE2 are not different, save them
            CE_reson = CE_shifts[0]
            QE_reson = HE_shifts[0]
        for atomname, shift in atomname_shift_dict.items():
            if atomname in ["CD1", "CD2"] and CD_reson != None:   # if CD_reson != None then QD_reson != None too
                new_atomname_shift_dict["CD1"] = CD_reson
                new_atomname_shift_dict["CD2"] = CD_reson
            elif atomname in ["HD1", "HD2"] and CD_reson != None:
                new_atomname_shift_dict["HD1"] = QD_reson
                new_atomname_shift_dict["HD2"] = QD_reson
            elif atomname in ["CE1", "CE2"] and CE_reson != None:   # if CE_reson != None then QE_reson != None too
                new_atomname_shift_dict["CE1"] = CE_reson
                new_atomname_shift_dict["CE2"] = CE_reson
            elif atomname in ["HE1", "HE2"] and CE_reson != None:
                new_atomname_shift_dict["HE1"] = QE_reson
                new_atomname_shift_dict["HE2"] = QE_reson
            elif not atomname in ["CD1", "CD2", "HD1", "HD2", "CE1", "CE2", "HE1", "HE2"]:
                new_atomname_shift_dict[atomname] = shift
    
    elif resname == "PHE":
        CD_shifts = [shift for atomname, shift in atomname_shift_dict.items() if atomname in ["CD1", "CD2"] ]
        HD_shifts = [shift for atomname, shift in atomname_shift_dict.items() if atomname in ["HD1", "HD2"] ]
        CE_shifts = [shift for atomname, shift in atomname_shift_dict.items() if atomname in ["CE1", "CE2"] ]
        HE_shifts = [shift for atomname, shift in atomname_shift_dict.items() if atomname in ["HE1", "HE2"] ]
        CD_reson, QD_reson, CE_reson, QE_reson = None, None, None, None
        if len(set(CD_shifts))==1 and len(set(HD_shifts))==1:     # if CD1 and CD2 or HD1 and HD2 are not different, save them
            CD_reson = CD_shifts[0]
            QD_reson = HD_shifts[0]
        if len(set(CE_shifts))==1 and len(set(HE_shifts))==1:     # if CE1 and CE2 or HE1 and HE2 are not different, save them
            CE_reson = CE_shifts[0]
            QE_reson = HE_shifts[0]
        for atomname, shift in atomname_shift_dict.items():
            if atomname in ["CD1", "CD2"] and CD_reson != None:   # if CD_reson != None then QD_reson != None too
                new_atomname_shift_dict["CD1"] = CD_reson
                new_atomname_shift_dict["CD2"] = CD_reson
            elif atomname in ["HD1", "HD2"] and CD_reson != None:
                new_atomname_shift_dict["HD1"] = QD_reson
                new_atomname_shift_dict["HD2"] = QD_reson
            elif atomname in ["CE1", "CE2"] and CE_reson != None:   # if CE_reson != None then QE_reson != None too
                new_atomname_shift_dict["CE1"] = CE_reson
                new_atomname_shift_dict["CE2"] = CE_reson
            elif atomname in ["HE1", "HE2"] and CE_reson != None:
                new_atomname_shift_dict["HE1"] = QE_reson
                new_atomname_shift_dict["HE2"] = QE_reson
            elif not atomname in ["CD1", "CD2", "HD1", "HD2", "CE1", "CE2", "HE1", "HE2"]:
                new_atomname_shift_dict[atomname] = shift
        
    return new_atomname_shift_dict


def duplicate_degenerate_methylenes(resname, atomname_shift_dict):
    """
        FUNCTION to duplicate the resonances of the CH2 protons that exist only once. These correspond to overlapping peaks.
    """
    global aa_CHpairs2beMerged_dict
    
    try:
        for CHpair_list in aa_CHpairs2beMerged_dict[resname]:
            Cname = CHpair_list[0].split('-')[0]
            Hname1 = CHpair_list[0].split('-')[1]
            Hname2 = CHpair_list[1].split('-')[1]
            if Hname1 in atomname_shift_dict.keys() and not Hname2 in atomname_shift_dict.keys():
                atomname_shift_dict[Hname2] = atomname_shift_dict[Hname1]
            elif Hname2 in atomname_shift_dict.keys() and not Hname1 in atomname_shift_dict.keys():
                atomname_shift_dict[Hname1] = atomname_shift_dict[Hname2]
    except KeyError:
        return
    
    
def save_all_shift_pairs(resname, raw_atomname_shift_dict, previous_resname, next_resname):
    """
        FUNCTION to save all C-H shift pairs of a particular residue. Basically updates aa_CHpair_2Dhist_multidict.
    """
    global aa_CHpair_2Dhist_multidict, CO_CA_CB_list
    saved_alternative_Hnames_list = []  # list with the alternative Hnames that were saved; with the purpose to avoid saving them multiple times
    atomname_shift_dict = parse_aromatic_resonances(resname, raw_atomname_shift_dict)   # take care of aromatic C-H pairs
    if args.MERGE_EQUIVALENT_CH_PAIRS == True:
        duplicate_degenerate_methylenes(resname, atomname_shift_dict)
    
    for CH_pair in aa_CHpair_2Dhist_multidict[resname].keys():
        [Cname, Hname] = CH_pair.split("-")
        if Cname in atomname_shift_dict.keys() and Hname in atomname_shift_dict.keys():
            aa_CHpair_2Dhist_multidict[resname][Cname+"-"+Hname][0].append(round(float(atomname_shift_dict[Cname]), 4))
            aa_CHpair_2Dhist_multidict[resname][Cname+"-"+Hname][1].append(round(float(atomname_shift_dict[Hname]), 4))
                #         equivalent_CH_pair = [CH for CH in CHpair_list if CH_pair!=CH][0]   # we know that len(CHpair_list)=2 and we need only the one

        elif resname in aa_alternativeAtomNames_dict.keys() and Cname in aa_alternativeAtomNames_dict[resname].keys() and Hname in aa_alternativeAtomNames_dict[resname][Cname]:
            for alt_Hname in aa_alternativeAtomNames_dict[resname][Cname][1:]:  # skip the first which is the conventional Hname
                if Cname in atomname_shift_dict.keys() and alt_Hname in atomname_shift_dict.keys():
                    conv_Hname = aa_alternativeAtomNames_dict[resname][Cname][0] # replace the alternative Hname with the conventional Hname
                    if conv_Hname in saved_alternative_Hnames_list: # if we have saved already this conv_Hname with an alternative name, skip this
                        continue
                    saved_alternative_Hnames_list.append(conv_Hname)
                    aa_CHpair_2Dhist_multidict[resname][Cname+"-"+conv_Hname][0].append(round(float(atomname_shift_dict[Cname]), 4))
                    aa_CHpair_2Dhist_multidict[resname][Cname+"-"+conv_Hname][1].append(round(float(atomname_shift_dict[alt_Hname]), 4))
                    #             equivalent_CH_pair = [CH for CH in CHpair_list if Cname+"-"+conv_Hname!=CH][0]   # we know that len(CHpair_list)=2 and we need only the one
                    break # it is uncessessary to check the rest of the alternative Hnames since we have already save one
    
    try:
        CO_CA_CB_aa_SecStruct_list.append((atomname_shift_dict["C"], atomname_shift_dict["CA"], atomname_shift_dict["CB"], resname, atomname_shift_dict["SS"]))
    except KeyError:
        pass
    try:
        CB_CG_aa_SecStruct_list.append((atomname_shift_dict["CB"], atomname_shift_dict["CG"], resname, atomname_shift_dict["SS"]))
    except KeyError:
        pass
    try:
        CB_CG_aa_SecStruct_list.append((atomname_shift_dict["CB"], atomname_shift_dict["CG1"], resname, atomname_shift_dict["SS"]))
    except KeyError:
        pass
    try:
        CB_CG_aa_SecStruct_list.append((atomname_shift_dict["CB"], atomname_shift_dict["CG2"], resname, atomname_shift_dict["SS"]))
    except KeyError:
        pass
    
    try:
        if resname == "ILE":
            CG_CD_aa_SecStruct_list.append((atomname_shift_dict["CG1"], atomname_shift_dict["CD1"], resname, atomname_shift_dict["SS"]))
        elif resname == "LEU":
            if "CD1" in atomname_shift_dict.keys() and "CG" in atomname_shift_dict.keys():
                CG_CD_aa_SecStruct_list.append((atomname_shift_dict["CG"], atomname_shift_dict["CD1"], resname, atomname_shift_dict["SS"]))
            if "CD2" in atomname_shift_dict.keys() and "CG" in atomname_shift_dict.keys():
                CG_CD_aa_SecStruct_list.append((atomname_shift_dict["CG"], atomname_shift_dict["CD2"], resname, atomname_shift_dict["SS"]))
        else:
            CG_CD_aa_SecStruct_list.append((atomname_shift_dict["CG"], atomname_shift_dict["CD"], resname, atomname_shift_dict["SS"]))
    except KeyError:
        pass
    
    if previous_resname != None:
        try:
            CO_CA_CB_SS_prevAA_curAA_nextAA_list.append((atomname_shift_dict["C"], atomname_shift_dict["CA"], atomname_shift_dict["CB"], atomname_shift_dict["SS"],
                                                     previous_resname, resname, next_resname))
        except KeyError:
            pass


def smoothen_distribution(proton, carbon, Hedges, Cedges):
    """
    proton: hydrogen chemical shift array
    carbon: carbon chemical shift array. Both proton & carbon must have the same length.
    RETURN:
    Z0:     array with the smoothed probability distribution.
    """
    
    xmin = proton.min()
    xmax = proton.max()
    ymin = carbon.min()
    ymax = carbon.max()
    binNox = len(Hedges)-1 + 0j  # convert the bin number in x-axis to complex number to emulate e.g. "100j"
    binNoy = len(Cedges)-1 + 0j    # convert the bin number in y-axis to complex number to emulate e.g. "100j"
    X, Y = np.mgrid[xmin:xmax:binNox, ymin:ymax:binNoy]
    positions = np.vstack([X.ravel(), Y.ravel()])   # stack the two arrays vertically (one in each row); creates a 2x10000 array
    values = np.vstack([proton, carbon])    # place the proton observations on the 1st row and the carbon observations of the 2nd; creates a 2xN array (N is the number of observations)
    kernel = stats.gaussian_kde(values) # apply
    Z = np.reshape(kernel(positions).T, X.shape)    # convert the 1D array with 10,000 values to a 2D array with 100x100 values
    Z0= Z/Z.sum()   # normalize the distribution (sum of all values must be 1, like in a 2D probability distribution)
    
    return Z0
    
    import matplotlib.pyplot as plt
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.imshow(np.rot90(Z0), cmap=plt.cm.gist_earth_r, extent=[xmin, xmax, ymin, ymax])
    ax.plot(proton, carbon, 'k.', markersize=2)
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])
    plt.show()


def write_coordinates(raw_x, raw_y, fname):
    
    global HOME_DIR
    
    with open(HOME_DIR + "/BMRB_data/" + fname, 'w') as f:
        for i in xrange(raw_x.shape[0]):
            f.write(str(raw_x[i])+"\t"+str(raw_y[i])+"\n")
    


carbon = None
hydrogen = None
previous_resid = -1000
previous_resname = None
previous2_resname = None
previous_Entry_ID = None
Assigned_chem_shift_list_ID = None
chain_ID = None
previous_chain_ID = None
atomname_shift_dict = {}    # dict of the form: atom type -> chemical shift values (this dict applies to one residue only every time)
protein_sequence = []
counter = 0
res_counter = 0
if args.BMRB_ENTRIES_FILE:
    cs_lines = read_BMRB_file()
    for line in cs_lines:
        resid = line[5]
        resname = line[6]
        atomname = line[7]
        atomtype = line[8]
        cs = line[10]
        cs_stdev = line[11]
        Entry_ID = line[-2]
        Assigned_chem_shift_list_ID = line[-1]
        if (previous_resid != -1000 and resid != previous_resid) or (Entry_ID != previous_Entry_ID and previous_Entry_ID != None and Assigned_chem_shift_list_ID != previous_Assigned_chem_shift_list_ID):
            save_all_shift_pairs(previous_resname, atomname_shift_dict, previous2_resname, resname) # ATTENTION: THIS IS WRONG previous2_resname !!!
            atomname_shift_dict = {}
            res_counter += 1
        atomname_shift_dict[atomname] = cs
        previous_resid = resid
        previous_resname = resname
        previous_Entry_ID = Entry_ID
        previous_Assigned_chem_shift_list_ID = Assigned_chem_shift_list_ID
        if previous2_resname != previous_resname:
            previous2_resname = previous_resname
        counter += 1
    
    print "Read ", counter, "lines and ", res_counter, " residues from file", args.BMRB_ENTRIES_FILE

elif args.VASCO_DIR:
    VASCO_full_path = os.path.realpath(args.VASCO_DIR)
    print "Reading VASCO Database from ", VASCO_full_path
    fnames=os.listdir(VASCO_full_path)
    fpattern = re.compile('bmr[0-9]+.[a-z0-9]{4}.vasco')
    VASCO_files_list = filter(fpattern.search, fnames)
    for VASCO_file in VASCO_files_list:
        print "DEBUG: reading file ", VASCO_file
        mo = re.search('bmr([0-9]+).([a-z0-9]{4}).vasco', VASCO_file)
        if mo:
            Entry_ID = mo.group(1)
            with open(VASCO_full_path + "/" + VASCO_file, 'r') as f:
                for line in f:
                    word_list = line.split()
                    try:
                        if word_list[0] == "#":
                            continue
                    except IndexError:
                        continue
                    resid = word_list[2]
                    resname = word_list[3]
                    atomname = word_list[5]
                    if resname in VASCO_naming_exceptions_multidict.keys() and atomname in VASCO_naming_exceptions_multidict[resname].keys():
                        atomname = VASCO_naming_exceptions_multidict[resname][atomname]
                    atomtype = word_list[7]
                    cs = word_list[8]
                    cs_stdev = word_list[9]
                    chain_ID = word_list[1]
                    ss = word_list[4]
                    if not int(previous_resid) in [int(resid)-1, int(resid)] or previous_Entry_ID != Entry_ID or previous_chain_ID != chain_ID:    # if this is a new protein, or there was a gap
                        protein_sequence = []           # reset the protein_sequence list
                    if (previous_resid != -1000 and resid != previous_resid) or (chain_ID != previous_chain_ID):
                        protein_sequence.append(resname)
                        try:
                            previous2_resname = protein_sequence[-3]
                        except IndexError:
                            previous2_resname = None
                        save_all_shift_pairs(previous_resname, atomname_shift_dict, previous2_resname, resname)
                        atomname_shift_dict = {}
                        res_counter += 1
                    atomname_shift_dict[atomname] = cs
                    atomname_shift_dict["SS"] = ss
                    previous_resid = resid
                    previous_resname = resname
                    previous_Entry_ID = Entry_ID
                    previous_chain_ID = chain_ID
                    counter += 1
                



with open(HOME_DIR + "/BMRB_data/" + "CO_CA_CB_aa_SS_resonances.dat", 'w') as f:
    for triplet in CO_CA_CB_aa_SecStruct_list:
        f.write(triplet[0] + "\t" + triplet[1] + "\t" + triplet[2] + "\t" + triplet[3] + "\t" + triplet[4] + "\n")

with open(HOME_DIR + "/BMRB_data/" + "CB_CG_aa_SS_resonances.dat", 'w') as f:
    for triplet in CB_CG_aa_SecStruct_list:
        f.write(triplet[0] + "\t" + triplet[1] + "\t" + triplet[2] + "\t" + triplet[3] + "\n")

with open(HOME_DIR + "/BMRB_data/" + "CG_CD_aa_SS_resonances.dat", 'w') as f:
    for triplet in CG_CD_aa_SecStruct_list:
        f.write(triplet[0] + "\t" + triplet[1] + "\t" + triplet[2] + "\t" + triplet[3] + "\n")

with open(HOME_DIR + "/BMRB_data/" + "CO_CA_CB_SS_prevAA_curAA_nextAA.dat", 'w') as f:
    for group7 in CO_CA_CB_SS_prevAA_curAA_nextAA_list:
        f.write(group7[0] + "\t" + group7[1] + "\t" + group7[2] + "\t" + group7[3] + "\t" + group7[4] + "\t" + group7[5] + "\t" + group7[6] + "\n")



if args.MERGE_EQUIVALENT_CH_PAIRS == True:
    for resname in aa_CHpairs2beMerged_dict.keys():
        for CH_pairs in aa_CHpairs2beMerged_dict[resname]:
            CH_pair1 = CH_pairs[0]
            CH_pair2 = CH_pairs[1]
            H_list1 = aa_CHpair_2Dhist_multidict[resname][CH_pair1][1]
            C_list1 = aa_CHpair_2Dhist_multidict[resname][CH_pair1][0]
            H_list2 = aa_CHpair_2Dhist_multidict[resname][CH_pair2][1]
            C_list2 = aa_CHpair_2Dhist_multidict[resname][CH_pair2][0]
            if not (resname == "LEU" and CH_pair1 in ['CD1-HD1', 'CD2-HD2'] or resname == "VAL" and CH_pair1 in ["CG1-HG1", "CG2-HG2"]) and (
                len(H_list1) != len(H_list2) or len(C_list1) != len(C_list2) ):
                print "ERROR: the C and H histograms of ", resname, CH_pairs, "do not have equal sizes!!!"
                print "DEBUG: len(H_list1)=", len(H_list1), "len(H_list2)=", len(H_list2), "len(C_list1)=", len(C_list1), "len(C_list2)=", len(C_list2)
                sys.exit(1)
            elif resname == "LEU" and CH_pair1 in ['CD1-HD1', 'CD2-HD2'] or resname == "VAL" and CH_pair1 in ["CG1-HG1", "CG2-HG2"]:
                aa_CHpair_2Dhist_multidict[resname][CH_pair1][1] = H_list1 + H_list2
                aa_CHpair_2Dhist_multidict[resname][CH_pair2][1] = H_list1 + H_list2
                aa_CHpair_2Dhist_multidict[resname][CH_pair1][0] = C_list1 + C_list2
                aa_CHpair_2Dhist_multidict[resname][CH_pair2][0] = C_list1 + C_list2
            else:
                for index in xrange(len(H_list1)):
                    if H_list1[index] != H_list2[index]:
                        aa_CHpair_2Dhist_multidict[resname][CH_pair1][1].append(H_list2[index])
                        aa_CHpair_2Dhist_multidict[resname][CH_pair2][1].append(H_list1[index])
                        aa_CHpair_2Dhist_multidict[resname][CH_pair1][0].append(C_list2[index])
                        aa_CHpair_2Dhist_multidict[resname][CH_pair2][0].append(C_list1[index])

for aa in aa_CHpair_2Dhist_multidict.keys():
        for CH_pair in aa_CHpair_2Dhist_multidict[aa].keys():
                print "DEBUG: aa=", aa, "CH_pair=", CH_pair, "len(aa_CHpair_2Dhist_multidict[aa][CH_pair])=", len(aa_CHpair_2Dhist_multidict[aa][CH_pair][1])
                print "DEBUG: aa=", aa, "CH_pair=", CH_pair, "np.mean(aa_CHpair_2Dhist_multidict[aa][CH_pair])=", np.mean(aa_CHpair_2Dhist_multidict[aa][CH_pair][1])
                raw_x = np.array(aa_CHpair_2Dhist_multidict[aa][CH_pair][1])  # hydrogen CS (x-axis coordinates)
                raw_y = np.array(aa_CHpair_2Dhist_multidict[aa][CH_pair][0])  # carbon CS (y-axis coordinates)
                write_coordinates(raw_x, raw_y, aa+"_"+CH_pair+"_correlated_2D_coords.txt")
                x, y = discard_outliers(raw_x, raw_y)  # discard outliers from C and H chemical shift arrays
                print "Writing file "+ aa+"_"+CH_pair+"_correlated_2D_coords.no_outliers.txt"
                write_coordinates(x, y, aa+"_"+CH_pair+"_correlated_2D_coords.no_outliers.txt")
                xmin = min(x)
                xmax = max(x)
                ymin = min(y)
                ymax = max(y)
                adjusted_Cedges_list = [e for e in Cedges_list if e > ymin-0.4 and e < ymax+0.4]
                adjusted_Hedges_list = [e for e in Hedges_list if e > xmin-0.08 and e < xmax+0.08]
                H, xedges, yedges = np.histogram2d(x, y, bins=(adjusted_Hedges_list, adjusted_Cedges_list))
                Hnorm = H/float(H.sum())   # normalize the histogram to get probabilities
                (xdim,ydim) = Hnorm.shape
                print "Writing file "+ aa+"_"+CH_pair+"_correlated_2Dhist.txt"
                with open(HOME_DIR + "/BMRB_data/" + aa+"_"+CH_pair+"_correlated_2Dhist.txt", 'w') as f:
                   for i in xrange(xdim):
                       for j in xrange(ydim):
                           f.write("%.2f" % xedges[i] + "\t" + "%.2f" % yedges[j] + "\t" + str(Hnorm[i][j]) + "\n")
                Hsmooth = smoothen_distribution(x, y, xedges, yedges)
                (xdim,ydim) = Hsmooth.shape
                print "Writing file "+ aa+"_"+CH_pair+"_correlated_2Dhist.smoothed.txt"
                with open(HOME_DIR + "/BMRB_data/" + aa+"_"+CH_pair+"_correlated_2Dhist.smoothed.txt", 'w') as f:
                    for i in xrange(xdim):
                        for j in xrange(ydim):
                            f.write("%.2f" %  xedges[i] + "\t" + "%.2f" % yedges[j] + "\t" + str(Hsmooth[i][j]) + "\n")

for aa in aa_CHpair_2Dhist_multidict.keys():
    for CH_pair in aa_CHpair_2Dhist_multidict[aa].keys():
        print "DEBUG: writing 1D and 2D uncorrelated histograms for ", aa, CH_pair
        raw_x = np.array(aa_CHpair_2Dhist_multidict[aa][CH_pair][1])  # hydrogen CS (x-axis coordinates)
        raw_y = np.array(aa_CHpair_2Dhist_multidict[aa][CH_pair][0])  # carbon CS (y-axis coordinates)
        x, y = discard_outliers(raw_x, raw_y)  # discard outliers from C and H chemical shift arrays
        print "DEBUG: x=", x
        xmin = min(x)
        xmax = max(x)
        print "DEBUG: y=", y
        ymin = min(y)
        ymax = max(y)
        adjusted_Cedges_list = [e for e in Cedges_list if e > ymin-0.4 and e < ymax+0.4]
        adjusted_Hedges_list = [e for e in Hedges_list if e > xmin-0.08 and e < xmax+0.08]
        Hh, xedges = np.histogram(x, bins=adjusted_Hedges_list, density=False)  # make 1D histogram for Hydrogen
        Hh_norm = Hh/float(Hh.sum())   # normalize the histogram to get probabilities
        print "DEBUG: Hh=", Hh
        print "DEBUG: xedges=", xedges.tolist()
        xdim = Hh_norm.shape[0]
        print "DEBUG: xdim=", xdim
        print "DEBUG: Hh_norm=", Hh_norm.tolist()
        Hhist_1D = []
        for i in xrange(xdim):
            Hhist_1D.append((round(xedges[i], 2), Hh_norm[i]) )
        print "DEBUG: Hhist_1D=", Hhist_1D
        xedge = Hhist_1D[0][0]-args.H_bin_length
        while round(xedge, 2) >= -9.98:
            Hhist_1D.insert(0, (round(xedge, 2), 0.00))
            xedge -= args.H_bin_length
        xedge = Hhist_1D[-1][0]+args.H_bin_length
        while round(xedge, 2) <= 33.98:
            Hhist_1D.append((round(xedge, 2), 0.00))
            xedge += args.H_bin_length
        print "Writing file "+ aa+"_"+CH_pair.split('-')[1]+"_hist.txt"
        with open(HOME_DIR + "/BMRB_data/" + aa+"_"+CH_pair.split('-')[1]+"_hist.txt", 'w') as f:
            for i in xrange(len(Hhist_1D)):
                f.write(str(Hhist_1D[i][0]) + "\t" + str(Hhist_1D[i][1]) + "\n")
        Hc, yedges = np.histogram(y, bins=adjusted_Cedges_list, density=False)  # make 1D histogram for Carbon
        Hc_norm = Hc/float(Hc.sum())   # normalize the histogram to get probabilities
        print "DEBUG: Hc=", Hc
        ydim = Hc_norm.shape[0]
        print "DEBUG: ydim=", ydim
        print "DEBUG: yedges=", yedges.tolist()
        print "DEBUG: Hc_norm=", Hc_norm.tolist()
        Hcist_1D = []
        for i in xrange(ydim):
            Hcist_1D.append((round(yedges[i], 2), Hc_norm[i]) )
        print "DEBUG: Hcist_1D=", Hcist_1D
        yedge = Hcist_1D[0][0]-args.C_bin_length
        while round(yedge, 2) >= 0.1:
            Hcist_1D.insert(0, (round(yedge, 2), 0.00))
            yedge -= args.C_bin_length
        yedge = Hcist_1D[-1][0]+args.C_bin_length
        while round(yedge, 2) <= 199.9:
            Hcist_1D.append((round(yedge, 2), 0.00))
            yedge += args.C_bin_length
        print "Writing file "+ aa+"_"+CH_pair.split('-')[0]+"_hist.txt"
        with open(HOME_DIR + "/BMRB_data/" + aa+"_"+CH_pair.split('-')[0]+"_hist.txt", 'w') as f:
            for i in xrange(len(Hcist_1D)):
                f.write(str(Hcist_1D[i][0]) + "\t" + str(Hcist_1D[i][1]) + "\n")
        weighted_H = np.zeros([xdim, ydim])
        H = np.zeros([xdim, ydim])
        for i in xrange(xdim):
                for j in xrange(ydim):
                    weighted_uncorr_prob = (args.H_weight * Hh_norm[i] + args.C_weight * Hc_norm[j])/(args.H_weight + args.C_weight)
                    weighted_H[i][j] = weighted_uncorr_prob
                    uncorr_prob = Hh_norm[i] * Hc_norm[j]
                    H[i][j] = uncorr_prob
        weighted_Hnorm = weighted_H/float(weighted_H.sum())   # normalize the histogram to get probabilities
        Hnorm = H/float(H.sum())   # normalize the histogram to get probabilities
        print "Writing file "+ aa+"_"+CH_pair+"_weighted_uncorrelated_2Dhist.txt"
        with open(HOME_DIR + "/BMRB_data/" + aa+"_"+CH_pair+"_weighted_uncorrelated_2Dhist.txt", 'w') as f:
            for i in xrange(xdim):
                for j in xrange(ydim):
                    f.write(str(xedges[i]) + "\t" + str(yedges[j]) + "\t" + str(weighted_Hnorm[i][j]) + "\n")
        print "Writing file "+ aa+"_"+CH_pair+"_uncorrelated_2Dhist.txt"
        with open(HOME_DIR + "/BMRB_data/" + aa+"_"+CH_pair+"_uncorrelated_2Dhist.txt", 'w') as f:
            for i in xrange(xdim):
                for j in xrange(ydim):
                    f.write(str(xedges[i]) + "\t" + str(yedges[j]) + "\t" + str(Hnorm[i][j]) + "\n")



# awk 'BEGIN{x='$(head -1 $fname | awk '{print $1}')'}{if($1 != x){print ""; x=$1};print $0}' $fname > ${fname}_for_GNUPLOT
#done