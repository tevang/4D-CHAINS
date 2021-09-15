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

import sys, re, os, csv
from operator import itemgetter
from collections import OrderedDict
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
from lib.probhist import *
from lib.global_func import *

## Set global variables
HOME_DIR = os.getcwd()
code3 = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
code1 = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
aa1to3_dict = dict((c1,c3) for c1,c3 in zip(code1,code3))
aa3to1_dict = dict((c3,c1) for c3,c1 in zip(code3,code1))

allowed_aa_atoms_dict = {
"ALA" : ["HA", "HB", "CA", "CB", "N", "H"],
"ARG" : ["HA", "HB2", "HB3", "HG2", "HG3", "HD2", "HD3", "CA", "CB", "CG", "CD", "N", "H"],
"ASP" : ["HA", "HB2", "HB3", "CA", "CB", "N", "H"],
"ASN" : ["HA", "HB2", "HB3", "CA", "CB", "N", "H", "ND2", "HD21", "HD22"],  # including side chain N-H
"CYS" : ["HA", "HB2", "HB3", "CA", "CB", "N", "H"],
"GLU" : ["HA", "HB2", "HB3", "HG2", "HG3", "CA", "CB", "CG", "N", "H"],
"GLN" : ["HA", "HB2", "HB3", "HG2", "HG3", "CA", "CB", "CG", "N", "H", "NE2", "HE21", "HE22"],  # including side chain N-H
"GLY" : ["HA2", "HA3", "CA", "N", "H"],
"HIS" : ["HA", "HB2", "HB3", "CA", "CB", "N", "H"],
"ILE" : ["HA", "HB", "HG12", "HG13", "HG2", "HD1", "CA", "CB", "CG1", "CG2", "CD1", "N", "H"],
"LEU" : ["HA", "HB2", "HB3", "HG", "HD1", "HD2", "CA", "CB", "CG", "CD1", "CD2", "N", "H"],
"LYS" : ["HA", "HB2", "HB3", "HG2", "HG3", "HD2", "HD3", "HE2", "HE3", "CA", "CB", "CG", "CD", "CE", "N", "H"],
"MET" : ["HA", "HB2", "HB3", "HG2", "HG3", "HE", "CA", "CB", "CG", "CE", "N", "H"],
"PHE" : ["HA", "HB2", "HB3", "CA", "CB", "N", "H"],
"PRO" : ["HA", "HB2", "HB3", "HG2", "HG3", "HD2", "HD3", "CA", "CB", "CG", "CD", "N"], # Prolines are not detected by the method at position "i" due to lack of HN hydrogen
"SER" : ["HA", "HB2", "HB3", "CA", "CB", "N", "H"],
"THR" : ["HA", "HB", "HG2", "CA", "CB", "CG2", "N", "H"],
"TRP" : ["HA", "HB2", "HB3", "CA", "CB", "N", "H", "NE1", "HE1"],   # including side chain N-H
"TYR" : ["HA", "HB2", "HB3", "CA", "CB", "N", "H"],
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


# a multidict with keys the amino acid --> carbon-hydrogen pair and values two arrays, the carbon chemical shifts and the respective hydrogen chemical shifts
# this multidict holds all the information in the input BMRB database and will be used to calculate the histograms
aa_CHpair_2Dhist_mdict = tree()
aa_CHpair_2Dhist_mdict["ALA"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_mdict["ALA"]["CB-HB"] = [[], []]

#"ARG" : ["HA", "HB2", "HB3", "HG2", "HG3", "HD2", "HD3", "CA", "CB", "CG", "CD"],
aa_CHpair_2Dhist_mdict["ARG"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_mdict["ARG"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_mdict["ARG"]["CB-HB3"] = [[], []]
aa_CHpair_2Dhist_mdict["ARG"]["CG-HG2"] = [[], []]
aa_CHpair_2Dhist_mdict["ARG"]["CG-HG3"] = [[], []]
aa_CHpair_2Dhist_mdict["ARG"]["CD-HD2"] = [[], []]
aa_CHpair_2Dhist_mdict["ARG"]["CD-HD3"] = [[], []]

#"ASP" : ["HA", "HB2", "HB3", "CA", "CB"],
aa_CHpair_2Dhist_mdict["ASP"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_mdict["ASP"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_mdict["ASP"]["CB-HB3"] = [[], []]

#"ASN" : ["HA", "HB2", "HB3", "CA", "CB"],
aa_CHpair_2Dhist_mdict["ASN"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_mdict["ASN"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_mdict["ASN"]["CB-HB3"] = [[], []]

#"CYS" : ["HA", "HB2", "HB3", "CA", "CB"],
aa_CHpair_2Dhist_mdict["CYS"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_mdict["CYS"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_mdict["CYS"]["CB-HB3"] = [[], []]

#"GLU" : ["HA", "HB2", "HB3", "HG2", "HG3", "CA", "CB", "CG"],
aa_CHpair_2Dhist_mdict["GLU"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_mdict["GLU"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_mdict["GLU"]["CB-HB3"] = [[], []]
aa_CHpair_2Dhist_mdict["GLU"]["CG-HG2"] = [[], []]
aa_CHpair_2Dhist_mdict["GLU"]["CG-HG3"] = [[], []]

#"GLN" : ["HA", "HB2", "HB3", "HG2", "HG3", "CA", "CB", "CG"],
aa_CHpair_2Dhist_mdict["GLN"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_mdict["GLN"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_mdict["GLN"]["CB-HB3"] = [[], []]
aa_CHpair_2Dhist_mdict["GLN"]["CG-HG2"] = [[], []]
aa_CHpair_2Dhist_mdict["GLN"]["CG-HG3"] = [[], []]

#
aa_CHpair_2Dhist_mdict["GLY"]["CA-HA2"] = [[], []]
aa_CHpair_2Dhist_mdict["GLY"]["CA-HA3"] = [[], []]

#"HIS" : ["HA", "HB2", "HB3", "CA", "CB"],
aa_CHpair_2Dhist_mdict["HIS"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_mdict["HIS"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_mdict["HIS"]["CB-HB3"] = [[], []]

#"ILE" : ["HA", "HB", "HG12", "HG13", "HG2", "HD1", "CA", "CB", "CG1", "CG2", "CD1"],
aa_CHpair_2Dhist_mdict["ILE"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_mdict["ILE"]["CB-HB"] = [[], []]
aa_CHpair_2Dhist_mdict["ILE"]["CG1-HG12"] = [[], []]
aa_CHpair_2Dhist_mdict["ILE"]["CG1-HG13"] = [[], []]
aa_CHpair_2Dhist_mdict["ILE"]["CG2-HG2"] = [[], []]
aa_CHpair_2Dhist_mdict["ILE"]["CD1-HD1"] = [[], []]

#"LEU" : ["HA", "HB2", "HB3", "HG", "HD1", "HD2", "CA", "CB", "CG", "CD1", "CD2"],
aa_CHpair_2Dhist_mdict["LEU"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_mdict["LEU"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_mdict["LEU"]["CB-HB3"] = [[], []]
aa_CHpair_2Dhist_mdict["LEU"]["CG-HG"] = [[], []]
aa_CHpair_2Dhist_mdict["LEU"]["CD1-HD1"] = [[], []]
aa_CHpair_2Dhist_mdict["LEU"]["CD2-HD2"] = [[], []]

#"LYS" : ["HA", "HB2", "HB3", "HG2", "HG3", "HD2", "HD3", "HE2", "HE3", "CA", "CB", "CG", "CD", "CE"],
aa_CHpair_2Dhist_mdict["LYS"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_mdict["LYS"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_mdict["LYS"]["CB-HB3"] = [[], []]
aa_CHpair_2Dhist_mdict["LYS"]["CG-HG2"] = [[], []]
aa_CHpair_2Dhist_mdict["LYS"]["CG-HG3"] = [[], []]
aa_CHpair_2Dhist_mdict["LYS"]["CD-HD2"] = [[], []]
aa_CHpair_2Dhist_mdict["LYS"]["CD-HD3"] = [[], []]
aa_CHpair_2Dhist_mdict["LYS"]["CE-HE2"] = [[], []]
aa_CHpair_2Dhist_mdict["LYS"]["CE-HE3"] = [[], []]

#"MET" : ["HA", "HB2", "HB3", "HG2", "HG3", "CA", "CB", "CG"],
aa_CHpair_2Dhist_mdict["MET"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_mdict["MET"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_mdict["MET"]["CB-HB3"] = [[], []]
aa_CHpair_2Dhist_mdict["MET"]["CG-HG2"] = [[], []]
aa_CHpair_2Dhist_mdict["MET"]["CG-HG3"] = [[], []]
aa_CHpair_2Dhist_mdict["MET"]["CE-HE"] = [[], []]
# aa_CHpair_2Dhist_mdict["MET"]["CE-HE1"] = [[], []]
# aa_CHpair_2Dhist_mdict["MET"]["CE-HE2"] = [[], []]
# aa_CHpair_2Dhist_mdict["MET"]["CE-HE3"] = [[], []]

#"PHE" : ["HA", "HB2", "HB3", "CA", "CB"],
aa_CHpair_2Dhist_mdict["PHE"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_mdict["PHE"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_mdict["PHE"]["CB-HB3"] = [[], []]
aa_CHpair_2Dhist_mdict["PHE"]["CD1-HD1"] = [[], []]
aa_CHpair_2Dhist_mdict["PHE"]["CD2-HD2"] = [[], []]
aa_CHpair_2Dhist_mdict["PHE"]["CE1-HE1"] = [[], []]
aa_CHpair_2Dhist_mdict["PHE"]["CE2-HE2"] = [[], []]
aa_CHpair_2Dhist_mdict["PHE"]["CZ-HZ"] = [[], []]

#"PRO" : ["HA", "HB2", "HB3", "HG2", "HG3", "HD2", "HD3", "CA", "CB", "CG", "CD"],
aa_CHpair_2Dhist_mdict["PRO"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_mdict["PRO"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_mdict["PRO"]["CB-HB3"] = [[], []]
aa_CHpair_2Dhist_mdict["PRO"]["CG-HG2"] = [[], []]
aa_CHpair_2Dhist_mdict["PRO"]["CG-HG3"] = [[], []]
aa_CHpair_2Dhist_mdict["PRO"]["CD-HD2"] = [[], []]
aa_CHpair_2Dhist_mdict["PRO"]["CD-HD3"] = [[], []]

#"SER" : ["HA", "HB2", "HB3", "CA", "CB"],
aa_CHpair_2Dhist_mdict["SER"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_mdict["SER"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_mdict["SER"]["CB-HB3"] = [[], []]

#"THR" : ["HA", "HB", "HG2", "CA", "CB", "CG2"],
aa_CHpair_2Dhist_mdict["THR"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_mdict["THR"]["CB-HB"] = [[], []]
aa_CHpair_2Dhist_mdict["THR"]["CG2-HG2"] = [[], []]

#"TRP" : ["HA", "HB2", "HB3", "CA", "CB"],
aa_CHpair_2Dhist_mdict["TRP"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_mdict["TRP"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_mdict["TRP"]["CB-HB3"] = [[], []]

#"TYR" : ["HA", "HB2", "HB3", "CA", "CB"],
aa_CHpair_2Dhist_mdict["TYR"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_mdict["TYR"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_mdict["TYR"]["CB-HB3"] = [[], []]
aa_CHpair_2Dhist_mdict["TYR"]["CD1-HD1"] = [[], []]
aa_CHpair_2Dhist_mdict["TYR"]["CD2-HD2"] = [[], []]
aa_CHpair_2Dhist_mdict["TYR"]["CE1-HE1"] = [[], []]
aa_CHpair_2Dhist_mdict["TYR"]["CE2-HE2"] = [[], []]

#"VAL" : ["HA", "HB", "HG1", "HG2", "CA", "CB", "CG1", "CG2"]
aa_CHpair_2Dhist_mdict["VAL"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_mdict["VAL"]["CB-HB"] = [[], []]
aa_CHpair_2Dhist_mdict["VAL"]["CG1-HG1"] = [[], []]
aa_CHpair_2Dhist_mdict["VAL"]["CG2-HG2"] = [[], []]

# Backbone & sidechain N-H
["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
aa_CHpair_2Dhist_mdict["ALA"]["N-H"] = [[], []]
aa_CHpair_2Dhist_mdict["ARG"]["N-H"] = [[], []]
aa_CHpair_2Dhist_mdict["ASN"]["N-H"] = [[], []]
aa_CHpair_2Dhist_mdict["ASN"]["ND2-HD21"] = [[], []]
aa_CHpair_2Dhist_mdict["ASN"]["ND2-HD22"] = [[], []]
aa_CHpair_2Dhist_mdict["ASP"]["N-H"] = [[], []]
aa_CHpair_2Dhist_mdict["CYS"]["N-H"] = [[], []]
aa_CHpair_2Dhist_mdict["GLU"]["N-H"] = [[], []]
aa_CHpair_2Dhist_mdict["GLN"]["N-H"] = [[], []]
aa_CHpair_2Dhist_mdict["GLN"]["NE2-HE21"] = [[], []]
aa_CHpair_2Dhist_mdict["GLN"]["NE2-HE22"] = [[], []]
aa_CHpair_2Dhist_mdict["GLY"]["N-H"] = [[], []]
aa_CHpair_2Dhist_mdict["HIS"]["N-H"] = [[], []]
aa_CHpair_2Dhist_mdict["ILE"]["N-H"] = [[], []]
aa_CHpair_2Dhist_mdict["LEU"]["N-H"] = [[], []]
aa_CHpair_2Dhist_mdict["LYS"]["N-H"] = [[], []]
aa_CHpair_2Dhist_mdict["MET"]["N-H"] = [[], []]
aa_CHpair_2Dhist_mdict["PHE"]["N-H"] = [[], []]
aa_CHpair_2Dhist_mdict["SER"]["N-H"] = [[], []]
aa_CHpair_2Dhist_mdict["THR"]["N-H"] = [[], []]
aa_CHpair_2Dhist_mdict["TRP"]["N-H"] = [[], []]
aa_CHpair_2Dhist_mdict["TRP"]["NE1-HE1"] = [[], []]
aa_CHpair_2Dhist_mdict["TYR"]["N-H"] = [[], []]
aa_CHpair_2Dhist_mdict["VAL"]["N-H"] = [[], []]



# a dict with keys the amino acid --> carbons and values list of the covalently bonded hydrogens
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
    ["CB-HB2", "CB-HB3"],
    ["ND2-HD21", "ND2-HD22"]    # side chain N-H
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
    ["CG-HG2", "CG-HG3"],
    ['NE2-HE21', "NE2-HE22"]    # side chain N-H
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
    ["CB-HB2", "CB-HB3"],
]
aa_CHpairs2beMerged_dict["TYR"] = [
    ["CB-HB2", "CB-HB3"]
]
aa_CHpairs2beMerged_dict["VAL"] = [
    ["CG1-HG1", "CG2-HG2"]
]


# List of tuples of the form (CO, CA, CB, aa_type, secondary structure) chemical shifts; to be used to measure the correlation between these resonances
CO_CA_CB_aa_SecStruct_list = []
CB_CG_aa_SecStruct_list = []
CG_CD_aa_SecStruct_list = []
CO_CA_CB_SS_prevAA_curAA_nextAA_list = []

# naming exceptions for VASCO: aa type -> C-H pair -> VASCO name
VASCO_naming_exceptions_mdict = tree()
VASCO_naming_exceptions_mdict["ALA"]["HB*"] = "HB"
VASCO_naming_exceptions_mdict["ILE"]["HD1*"] = "HD1"
VASCO_naming_exceptions_mdict["ILE"]["HG2*"] = "HG2"
VASCO_naming_exceptions_mdict["LEU"]["HD1*"] = "HD1"
VASCO_naming_exceptions_mdict["LEU"]["HD2*"] = "HD2"
VASCO_naming_exceptions_mdict["THR"]["HG2*"] = "HG2"
VASCO_naming_exceptions_mdict["VAL"]["HG1*"] = "HG1"
VASCO_naming_exceptions_mdict["VAL"]["HG2*"] = "HG2"
VASCO_naming_exceptions_mdict["MET"]["HE*"] = "HE2"

## Parse command line arguments
def cmdlineparse():
    parser = ArgumentParser(description="command line arguments",
        epilog="""
        EXAMPLE: 
        create_2D_histograms.py \\
        -bmrbfile BMRB_data/Atom_chem_shift.proteins_with_sidechains.csv \\
        -selfile BMRB_data/Rescue2_BMRB_selected_sequences.txt
        """)
    parser.add_argument("-bmrbfile", dest="BMRB_ENTRIES_FILE", required=False, type=str, help="csv file with all the chemical shifts deposited in BMRB",
                        metavar="<BMRB entries file>")
    parser.add_argument("-vascodir", dest="VASCO_DIR", required=False, type=str, help="folder with all the VASCO database files",
                        metavar="<VASCO database directory>")
    parser.add_argument("-newvascodir", dest="NEW_VASCO_DIR", required=False, type=str,
                        help="folder with all the new VASCO database files",
                        metavar="<new VASCO database directory>")
    parser.add_argument("-newvascofile", dest="NEW_VASCO_FILE", required=False, type=str,
                        help="file from the new VASCO database files. Useful for debugging.",
                        metavar="<new VASCO database file>")    # NOT IN USE YET
    parser.add_argument("-selfile", dest="SELECTED_ENTRIES_FILE", required=False, type=str,
                        help="file with a list of BMRB entries to be used to make 2D histograms", metavar="<selected BMRB entries files>")
    parser.add_argument("-wH", dest="H_weight", required=False, type=float, default=0.1, help="weight in (0.0-1.0] to apply on aa type prediction from aliphatic H resonances (default: 0.1)", metavar="<H weight>")
    parser.add_argument("-wC", dest="C_weight", required=False, type=float, default=1.0, help="weight in (0.0, 1.0] to apply on aa type prediction from aliphatic C resonances (default: 1.0)", metavar="<C weight>")
    parser.add_argument("-resH", dest="H_bin_length", required=False, type=float, default=0.04,
                        help="histogram bin length for the carbon in ppm. (default: %(default)s)", metavar="<H bin length>")
    parser.add_argument("-resC", dest="C_bin_length", required=False, type=float, default=0.2,
                        help="histogram bin length for the hydrogen in ppm. (default: %(default)s)", metavar="<C bin length>")
    parser.add_argument("-nomerge", dest="MERGE_EQUIVALENT_CH_PAIRS", required=False, default=True, action='store_false',
                        help="Merge equivalent C-H pairs, e.g. for all CB methylenes, CB-HB2 and CB-HB3 will have identical 2D-histograms.")
    parser.add_argument("-nofilter", dest="NO_FILTER", required=False, default=False, action='store_true',
                        help="keep carbon resonances that have 0 probability in the 1D carbon histogram (NOT RECOMMENDED).")
    parser.add_argument("-probthres", dest="PROBABILITY_THRESHOLD", required=False, type=float, default=1e-20,
                        help="whichever carbon resonance has 1D-histogram probability below this threshold will be discarded. (default %(default)s). "
                             "This arguments is pointless if you use -nofilter.",
                        metavar="<>")
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
    
    ##print "DEBUG: Cshifts_list=", Cshifts.tolist()
    ##print "DEBUG: Hshifts_list=", Hshifts.tolist()
    ## Find the indices of C shift outliers
    if len(Cshifts.shape) == 1:
        Cshifts = Cshifts[:,None]
    #median = np.median(Cshifts, axis=0)
    #diff = np.sum((Cshifts - median)**2, axis=-1)
    #diff = np.sqrt(diff)
    #med_abs_deviation = np.median(diff)
    ##print "DEBUG: diff=", diff
    ##print "DEBUG: med_abs_deviation=", med_abs_deviation
    #modified_z_score = 0.6745 * diff / med_abs_deviation
    #upper_discard_C = modified_z_score > abs(C_thresh)    # upper outliers
    #lower_discard_C = modified_z_score < -1*abs(C_thresh)    # lower outliers
    #final_discard_C = upper_discard_C | lower_discard_C # all the outliers
    #
    ## Find the indices of H shift outliers
    if len(Hshifts.shape) == 1:
        Hshifts = Hshifts[:,None]
    #median = np.median(Hshifts, axis=0)
    #diff = np.sum((Hshifts - median)**2, axis=-1)
    #diff = np.sqrt(diff)
    #med_abs_deviation = np.median(diff)
    ##print "DEBUG: diff=", diff
    ##print "DEBUG: med_abs_deviation=", med_abs_deviation
    #modified_z_score = 0.6745 * diff / med_abs_deviation
    #upper_discard_H = modified_z_score > abs(H_thresh)    # upper outliers
    #lower_discard_H = modified_z_score < -1*abs(H_thresh)    # lower outliers
    #final_discard_H = upper_discard_H | lower_discard_H
    
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
    
    print("Loading BMRB file with chemical shifts ...")
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
    
    #print "DEBUG: atomname_shift_dict=", atomname_shift_dict
    new_atomname_shift_dict = {}
    if not resname in ["TYR", "PHE"]:   # if it's not aromatic, no parsing is needed
        return atomname_shift_dict
    elif resname == "TYR":
        CD_shifts = [shift for atomname, shift in list(atomname_shift_dict.items()) if atomname in ["CD1", "CD2"] ]
        HD_shifts = [shift for atomname, shift in list(atomname_shift_dict.items()) if atomname in ["HD1", "HD2"] ]
        CE_shifts = [shift for atomname, shift in list(atomname_shift_dict.items()) if atomname in ["CE1", "CE2"] ]
        HE_shifts = [shift for atomname, shift in list(atomname_shift_dict.items()) if atomname in ["HE1", "HE2"] ]
        CD_reson, QD_reson, CE_reson, QE_reson = None, None, None, None
        if len(set(CD_shifts))==1 and len(set(HD_shifts))==1:     # if CD1 and CD2 or HD1 and HD2 are not different, save them
            CD_reson = CD_shifts[0]
            QD_reson = HD_shifts[0]
        if len(set(CE_shifts))==1 and len(set(HE_shifts))==1:     # if CE1 and CE2 or HE1 and HE2 are not different, save them
            CE_reson = CE_shifts[0]
            QE_reson = HE_shifts[0]
        # NOW POPULATE new_atomname_shift_dict
        for atomname, shift in list(atomname_shift_dict.items()):
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
        CD_shifts = [shift for atomname, shift in list(atomname_shift_dict.items()) if atomname in ["CD1", "CD2"] ]
        HD_shifts = [shift for atomname, shift in list(atomname_shift_dict.items()) if atomname in ["HD1", "HD2"] ]
        CE_shifts = [shift for atomname, shift in list(atomname_shift_dict.items()) if atomname in ["CE1", "CE2"] ]
        HE_shifts = [shift for atomname, shift in list(atomname_shift_dict.items()) if atomname in ["HE1", "HE2"] ]
        CD_reson, QD_reson, CE_reson, QE_reson = None, None, None, None
        if len(set(CD_shifts))==1 and len(set(HD_shifts))==1:     # if CD1 and CD2 or HD1 and HD2 are not different, save them
            CD_reson = CD_shifts[0]
            QD_reson = HD_shifts[0]
        if len(set(CE_shifts))==1 and len(set(HE_shifts))==1:     # if CE1 and CE2 or HE1 and HE2 are not different, save them
            CE_reson = CE_shifts[0]
            QE_reson = HE_shifts[0]
        # NOW POPULATE new_atomname_shift_dict
        for atomname, shift in list(atomname_shift_dict.items()):
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
        FUNCTION to duplicate the resonances of the CH2 protons and the C-H resonances of LEU CD2 or VAL CG2 that exist only once.
        These correspond to overlapping peaks.
    """
    global aa_CHpairs2beMerged_dict

    # print "DEBUG: before duplication resname", resname, "atomname_shift_dict=", atomname_shift_dict
    try:
        for CHpair_list in aa_CHpairs2beMerged_dict[resname]:
            Cname = CHpair_list[0].split('-')[0]
            Hname1 = CHpair_list[0].split('-')[1]
            Hname2 = CHpair_list[1].split('-')[1]
            if Hname1 in list(atomname_shift_dict.keys()) and not Hname2 in list(atomname_shift_dict.keys()):
                # print "DEBUG: duplicating", resname, Hname1, atomname_shift_dict[Hname1]
                atomname_shift_dict[Hname2] = atomname_shift_dict[Hname1]
            elif Hname2 in list(atomname_shift_dict.keys()) and not Hname1 in list(atomname_shift_dict.keys()):
                # print "DEBUG: duplicating", resname, Hname2, atomname_shift_dict[Hname2]
                atomname_shift_dict[Hname1] = atomname_shift_dict[Hname2]
    except KeyError:
        return
    # print "DEBUG: after duplication resname", resname, "atomname_shift_dict=", atomname_shift_dict


def even_duplicate_degenerate_methylenes(resname, atomname_shift_dict):
    """
        FUNCTION to remove one of the duplicate resonances of the CH2 protons and the C-H resonances of LEU CD2 or VAL CG2. These
        correspond to case where the user did not find the second methylene proton and just copied the resonance of the first in
        the assignment file.
    """
    global aa_CHpairs2beMerged_dict

    # print "DEBUG: before even resname", resname, "atomname_shift_dict=", atomname_shift_dict
    try:
        for CHpair_list in aa_CHpairs2beMerged_dict[resname]:
            Cname = CHpair_list[0].split('-')[0]
            Hname1 = CHpair_list[0].split('-')[1]
            Hname2 = CHpair_list[1].split('-')[1]
            if Hname1 in list(atomname_shift_dict.keys()) and Hname2 in list(atomname_shift_dict.keys()) and atomname_shift_dict[
                Hname1] == atomname_shift_dict[Hname2]:
                # print "DEBUG: duplicating", resname, Hname1, atomname_shift_dict[Hname1]
                del atomname_shift_dict[Hname2]
            if (
                    resname == "LEU" and Cname == "CD1" and Hname1 in list(atomname_shift_dict.keys()) and Hname2 in list(atomname_shift_dict.keys()) and
                    atomname_shift_dict[Hname1] == atomname_shift_dict[Hname2] and
                    "CD1" in list(atomname_shift_dict.keys()) and "CD2" in list(atomname_shift_dict.keys()) and atomname_shift_dict[
                "CD1"] == atomname_shift_dict["CD2"]):
                del atomname_shift_dict["CD2"], atomname_shift_dict[Hname2]
            if (
                    resname == "VAL" and Cname == "CG1" and Hname1 in list(atomname_shift_dict.keys()) and Hname2 in list(atomname_shift_dict.keys()) and
                    atomname_shift_dict[Hname1] == atomname_shift_dict[Hname2] and
                    "CG1" in list(atomname_shift_dict.keys()) and "CG2" in list(atomname_shift_dict.keys()) and atomname_shift_dict[
                "CG1"] == atomname_shift_dict["CG2"]):
                del atomname_shift_dict["CG2"], atomname_shift_dict[Hname2]
    except KeyError:
        return

    # print "DEBUG: after even resname", resname, "atomname_shift_dict=", atomname_shift_dict

# Load 1D histograms for the following function
histload = ProbHist_Loader()
histload.load_1Dhistograms()
histcalc = ProbHist_Calculator()
def remove_C_outliers(resname, atomname_shift_dict):
    """
        using 1D Carbon histograms.
    """
    global histload, histcalc, args

    for atomname in list(atomname_shift_dict.keys()):
        if atomname[0] == 'C' and atomname in histload.aa_carbon_binDensityList_mdict[resname]:
            Creson = atomname_shift_dict[atomname]
            [bin_array, density_array] = histload.aa_carbon_binDensityList_mdict[resname][atomname]
            probability = histcalc.get_probability_from_histogram(Creson, bin_array, density_array)
            if probability < args.PROBABILITY_THRESHOLD:
                # print "DEBUG: removing outlier:", resname, atomname, atomname_shift_dict[atomname]
                del atomname_shift_dict[atomname]

    return atomname_shift_dict


def remove_replicate_Cresons(atomname_shift_dict, exclude_C=[]):
    """
        Remove all carbons of the same residue that have the same resonance value.
    """

    Creson_list = [s for a, s in list(atomname_shift_dict.items()) if a[0] == 'C']
    for atomname, shift in list(atomname_shift_dict.items()):
        if atomname[0] == 'C' and Creson_list.count(shift) > 1 and not atomname in exclude_C:
            del atomname_shift_dict[atomname]

    return atomname_shift_dict

# def save_all_shift_pairs(resname, raw_atomname_shift_dict, previous_resname, next_resname):
#     """
#         OLD FUNCTION to save all C-H shift pairs of a particular residue. Basically updates aa_CHpair_2Dhist_mdict.
#     """
#     global aa_CHpair_2Dhist_mdict, CO_CA_CB_aa_SecStruct_list, CB_CG_aa_SecStruct_list
#     global CG_CD_aa_SecStruct_list, CO_CA_CB_SS_prevAA_curAA_nextAA_list
#
#     # print "DEBUG: resname=", resname, "raw_atomname_shift_dict=", raw_atomname_shift_dict
#     saved_alternative_Hnames_list = []  # list with the alternative Hnames that were saved; with the purpose to avoid saving them multiple times
#     #if resname == "VAL" and "CA" in atomname_shift_dict.keys():
#     #    print "DEBUG: VAL atomname_shift_dict[CA]", atomname_shift_dict["CA"]
#     #    print "DEBUG: ",resname,"atomname_shift_dict=", atomname_shift_dict
#     #print "DEBUG: previous resname=", previous_resname, "resname=", resname, "next_resname=", next_resname
#     atomname_shift_dict = parse_aromatic_resonances(resname, raw_atomname_shift_dict)   # take care of aromatic C-H pairs
#     if args.MERGE_EQUIVALENT_CH_PAIRS == True:
#         duplicate_degenerate_methylenes(resname, atomname_shift_dict)
#     #atomname_shift_dict = raw_atomname_shift_dict
#
#     for CH_pair in aa_CHpair_2Dhist_mdict[resname].keys():
#         [Cname, Hname] = CH_pair.split("-")
#         if Cname in atomname_shift_dict.keys() and Hname in atomname_shift_dict.keys():
#             # if resname == "GLY" and Cname == "CA" and Hname in ["HA2", "HA3"]:
#             #    print "\n\nDEBUG: saving ", float(atomname_shift_dict[Cname]), "and", float(atomname_shift_dict[Hname]), " to ", resname, Cname+"-"+Hname
#             # print "DEBUG: saving ", float(atomname_shift_dict[Cname]), "and", float(atomname_shift_dict[Hname]), " to ", resname, Cname+"-"+Hname
#             aa_CHpair_2Dhist_mdict[resname][Cname+"-"+Hname][0].append(round(float(atomname_shift_dict[Cname]), 4))
#             aa_CHpair_2Dhist_mdict[resname][Cname+"-"+Hname][1].append(round(float(atomname_shift_dict[Hname]), 4))
#             # if args.MERGE_EQUIVALENT_CH_PAIRS == True:
#                 # for CHpair_list in aa_CHpairs2beMerged_dict[resname]:
#                 #     if CH_pair in CHpair_list:
#                 #         equivalent_CH_pair = [CH for CH in CHpair_list if CH_pair!=CH][0]   # we know that len(CHpair_list)=2 and we need only the one
#                 #         # print "DEBUG: saving ", CH_pair, " and its equivalent ", equivalent_CH_pair
#                 #         [eq_Cname, eq_Hname] = equivalent_CH_pair.split("-")
#                 #         # print "DEBUG: saving ", float(atomname_shift_dict[Cname]), "and", float(atomname_shift_dict[Hname]), " to ", resname, eq_Cname+"-"+eq_Hname
#                 #         aa_CHpair_2Dhist_mdict[resname][eq_Cname+"-"+eq_Hname][0].append(round(float(atomname_shift_dict[Cname]), 4))
#                 #         aa_CHpair_2Dhist_mdict[resname][eq_Cname+"-"+eq_Hname][1].append(round(float(atomname_shift_dict[Hname]), 4))
#
#         elif resname in aa_alternativeAtomNames_dict.keys() and Cname in aa_alternativeAtomNames_dict[resname].keys() and Hname in aa_alternativeAtomNames_dict[resname][Cname]:
#             for alt_Hname in aa_alternativeAtomNames_dict[resname][Cname][1:]:  # skip the first which is the conventional Hname
#                 if Cname in atomname_shift_dict.keys() and alt_Hname in atomname_shift_dict.keys():
#                     conv_Hname = aa_alternativeAtomNames_dict[resname][Cname][0] # replace the alternative Hname with the conventional Hname
#                     if conv_Hname in saved_alternative_Hnames_list: # if we have saved already this conv_Hname with an alternative name, skip this
#                         continue
#                     saved_alternative_Hnames_list.append(conv_Hname)
#                     # print "DEBUG: found existing alternative Hname ", Hname, " and saved it as ", conv_Hname
#                     # print "\n\nDEBUG: saving ", float(atomname_shift_dict[Cname]), "and", float(atomname_shift_dict[alt_Hname]), " to ", resname, Cname+"-"+conv_Hname
#                     aa_CHpair_2Dhist_mdict[resname][Cname+"-"+conv_Hname][0].append(round(float(atomname_shift_dict[Cname]), 4))
#                     aa_CHpair_2Dhist_mdict[resname][Cname+"-"+conv_Hname][1].append(round(float(atomname_shift_dict[alt_Hname]), 4))
#                     # if args.MERGE_EQUIVALENT_CH_PAIRS == True:
#                     #     for CHpair_list in aa_CHpairs2beMerged_dict[resname]:
#                     #         if Cname+"-"+conv_Hname in CHpair_list:
#                     #             equivalent_CH_pair = [CH for CH in CHpair_list if Cname+"-"+conv_Hname!=CH][0]   # we know that len(CHpair_list)=2 and we need only the one
#                     #             #print "DEBUG: saving ", Cname+"-"+conv_Hname, " and its equivalent ", equivalent_CH_pair
#                     #             [eq_Cname, eq_Hname] = equivalent_CH_pair.split("-")
#                     #             #if resname == "VAL" and Cname == "CG2":
#                     #             #    print "DEBUG: saving ", float(atomname_shift_dict[Cname]), "and", float(atomname_shift_dict[alt_Hname]), " to ", resname, eq_Cname+"-"+eq_Hname
#                     #             aa_CHpair_2Dhist_mdict[resname][eq_Cname+"-"+eq_Hname][0].append(round(float(atomname_shift_dict[Cname]), 4))
#                     #             aa_CHpair_2Dhist_mdict[resname][eq_Cname+"-"+eq_Hname][1].append(round(float(atomname_shift_dict[alt_Hname]), 4))
#                     break # it is uncessessary to check the rest of the alternative Hnames since we have already save one
#
#     # Save the CO, CA, CB shift to measure their correlation at the end
#     try:
#         #print "DEBUG: atomname_shift_dict.keys()=", atomname_shift_dict.keys()
#         #print "DEBUG: appending:", atomname_shift_dict["C"], atomname_shift_dict["CA"], atomname_shift_dict["CB"], resname, atomname_shift_dict["SS"]
#         CO_CA_CB_aa_SecStruct_list.append((atomname_shift_dict["C"], atomname_shift_dict["CA"], atomname_shift_dict["CB"], resname, atomname_shift_dict["SS"]))
#     except KeyError:
#         pass
#     # Save the CB, CG shift to measure their correlation at the end
#     try:
#         CB_CG_aa_SecStruct_list.append((atomname_shift_dict["CB"], atomname_shift_dict["CG"], resname, atomname_shift_dict["SS"]))
#     except KeyError:
#         pass
#     try:
#         CB_CG_aa_SecStruct_list.append((atomname_shift_dict["CB"], atomname_shift_dict["CG1"], resname, atomname_shift_dict["SS"]))
#     except KeyError:
#         pass
#     try:
#         CB_CG_aa_SecStruct_list.append((atomname_shift_dict["CB"], atomname_shift_dict["CG2"], resname, atomname_shift_dict["SS"]))
#     except KeyError:
#         pass
#
#     # Save the CG, CD shifts to measure their correlation
#     try:
#         if resname == "ILE":
#             CG_CD_aa_SecStruct_list.append((atomname_shift_dict["CG1"], atomname_shift_dict["CD1"], resname, atomname_shift_dict["SS"]))
#         elif resname == "LEU":
#             if "CD1" in atomname_shift_dict.keys() and "CG" in atomname_shift_dict.keys():
#                 CG_CD_aa_SecStruct_list.append((atomname_shift_dict["CG"], atomname_shift_dict["CD1"], resname, atomname_shift_dict["SS"]))
#             if "CD2" in atomname_shift_dict.keys() and "CG" in atomname_shift_dict.keys():
#                 CG_CD_aa_SecStruct_list.append((atomname_shift_dict["CG"], atomname_shift_dict["CD2"], resname, atomname_shift_dict["SS"]))
#         else:
#             CG_CD_aa_SecStruct_list.append((atomname_shift_dict["CG"], atomname_shift_dict["CD"], resname, atomname_shift_dict["SS"]))
#     except KeyError:
#         pass
#
#     # Save CO, CA, CB, SS, previous aa, current aa, next aa
#     if previous_resname != None:
#         try:
#             CO_CA_CB_SS_prevAA_curAA_nextAA_list.append((atomname_shift_dict["C"], atomname_shift_dict["CA"], atomname_shift_dict["CB"], atomname_shift_dict["SS"],
#                                                      previous_resname, resname, next_resname))
#         except KeyError:
#             pass


def save_all_shift_pairs(resname, raw_atomname_shift_dict, previous_resname, next_resname):
    """
        NEW FUNCTION to save all C-H shift pairs of a particular residue. Basically updates aa_CHpair_2Dhist_mdict.
    """
    global aa_CHpair_2Dhist_mdict, CO_CA_CB_aa_SecStruct_list, CB_CG_aa_SecStruct_list
    global CG_CD_aa_SecStruct_list, CO_CA_CB_SS_prevAA_curAA_nextAA_list

    # print "DEBUG: resname=", resname, "raw_atomname_shift_dict=", raw_atomname_shift_dict
    saved_alternative_Hnames_list = []  # list with the alternative Hnames that were saved; with the purpose to avoid saving them multiple times
    # if resname == "VAL" and "CA" in atomname_shift_dict.keys():
    #    print "DEBUG: VAL atomname_shift_dict[CA]", atomname_shift_dict["CA"]
    #    print "DEBUG: ",resname,"atomname_shift_dict=", atomname_shift_dict
    # print "DEBUG: previous resname=", previous_resname, "resname=", resname
    atomname_shift_dict = parse_aromatic_resonances(resname, raw_atomname_shift_dict)  # take care of aromatic C-H pairs

    if args.NEW_VASCO_DIR or args.NEW_VASCO_FILE:
        assign_ambiguous_degenerate_methylenes(resname, atomname_shift_dict)

    if args.MERGE_EQUIVALENT_CH_PAIRS == True:
        duplicate_degenerate_methylenes(resname, atomname_shift_dict)
    else:
        even_duplicate_degenerate_methylenes(resname, atomname_shift_dict)
    # atomname_shift_dict = raw_atomname_shift_dict

    # DELETE REPLICATED CARBONS
    if resname in ["PHE", "TYR"]:
        atomname_shift_dict = remove_replicate_Cresons(atomname_shift_dict, exclude_C=["CD1", "CD2", "CE1", "CE2"])
    else:
        atomname_shift_dict = remove_replicate_Cresons(atomname_shift_dict)


    new_atomname_shift_dict = {}  # tmp storage of the assigned C-H pairs of this record
    for CH_pair in list(aa_CHpair_2Dhist_mdict[resname].keys()):
        [Cname, Hname] = CH_pair.split("-")
        # print "DEBUG: Cname=", Cname, "Hname=", Hname
        if Cname in list(atomname_shift_dict.keys()) and Hname in list(atomname_shift_dict.keys()):
            # print "DEBUG: resname", resname, "Cname", Cname, "atomname_shift_dict[Cname]=", atomname_shift_dict[Cname], "Hname", Hname, "atomname_shift_dict[Hname]=", atomname_shift_dict[Hname]
            if Cname in list(new_atomname_shift_dict.keys()):  # if aleady a carbon CS, average it
                new_atomname_shift_dict[Cname] = (new_atomname_shift_dict[Cname] + round(
                    float(atomname_shift_dict[Cname]), 4)) / 2.0
            else:
                new_atomname_shift_dict[Cname] = round(float(atomname_shift_dict[Cname]), 4)
            new_atomname_shift_dict[Hname] = round(float(atomname_shift_dict[Hname]), 4)
            # print "DEBUG: new_atomname_shift_dict=", new_atomname_shift_dict

        elif resname in list(aa_alternativeAtomNames_dict.keys()) and Cname in list(aa_alternativeAtomNames_dict[
            resname].keys()) and Hname in aa_alternativeAtomNames_dict[resname][Cname]:
            for alt_Hname in aa_alternativeAtomNames_dict[resname][Cname][
                             1:]:  # skip the first which is the conventional Hname
                if Cname in list(atomname_shift_dict.keys()) and alt_Hname in list(atomname_shift_dict.keys()):
                    conv_Hname = aa_alternativeAtomNames_dict[resname][Cname][
                        0]  # replace the alternative Hname with the conventional Hname
                    if conv_Hname in saved_alternative_Hnames_list:  # if we have saved already this conv_Hname with an alternative name, skip this
                        continue
                    saved_alternative_Hnames_list.append(conv_Hname)
                    if Cname in list(new_atomname_shift_dict.keys()):  # if aleady a carbon CS, average it
                        new_atomname_shift_dict[Cname] = (new_atomname_shift_dict[Cname] + round(
                            float(atomname_shift_dict[Cname]), 4)) / 2.0
                    else:
                        new_atomname_shift_dict[Cname] = round(float(atomname_shift_dict[Cname]), 4)
                    # print "DEBUG: alt_Hname", alt_Hname, "conv_Hname=", conv_Hname, "atomname_shift_dict[alt_Hname]=", atomname_shift_dict[alt_Hname]
                    new_atomname_shift_dict[conv_Hname] = round(float(atomname_shift_dict[alt_Hname]), 4)
                    break  # it is uncessessary to check the rest of the alternative Hnames since we have already save one

    # Finally save the N, H resonances
    if 'N' in list(atomname_shift_dict.keys()):
        new_atomname_shift_dict['N'] = atomname_shift_dict['N']
    if 'H' in list(atomname_shift_dict.keys()):
        new_atomname_shift_dict['H'] = atomname_shift_dict['H']

    if not args.NO_FILTER:
        # Remove Carbon outliers from atomname_shift_dict
        new_atomname_shift_dict = remove_C_outliers(resname, new_atomname_shift_dict)

    # Finally append the contents of new_atomname_shift_dict to aa_CHpair_2Dhist_mdict
    for CH_pair in list(aa_CHpair_2Dhist_mdict[resname].keys()):
        [Cname, Hname] = CH_pair.split("-")
        if not Cname in list(new_atomname_shift_dict.keys()) or not Hname in list(new_atomname_shift_dict.keys()):
            continue
        aa_CHpair_2Dhist_mdict[resname][Cname + "-" + Hname][0].append(round(float(new_atomname_shift_dict[Cname]), 4))
        aa_CHpair_2Dhist_mdict[resname][Cname + "-" + Hname][1].append(round(float(new_atomname_shift_dict[Hname]), 4))

    # # Save the CO, CA, CB shift to measure their correlation at the end
    # try:
    #     # print "DEBUG: atomname_shift_dict.keys()=", atomname_shift_dict.keys()
    #     #print "DEBUG: appending:", atomname_shift_dict["C"], atomname_shift_dict["CA"], atomname_shift_dict["CB"], resname, atomname_shift_dict["SS"]
    #     CO_CA_CB_aa_SecStruct_list.append((atomname_shift_dict["C"], atomname_shift_dict["CA"], atomname_shift_dict["CB"], resname, atomname_shift_dict["SS"]))
    # except KeyError:
    #     pass
    # # Save the CB, CG shift to measure their correlation at the end
    # try:
    #     CB_CG_aa_SecStruct_list.append((atomname_shift_dict["CB"], atomname_shift_dict["CG"], resname, atomname_shift_dict["SS"]))
    # except KeyError:
    #     pass
    # try:
    #     CB_CG_aa_SecStruct_list.append((atomname_shift_dict["CB"], atomname_shift_dict["CG1"], resname, atomname_shift_dict["SS"]))
    # except KeyError:
    #     pass
    # try:
    #     CB_CG_aa_SecStruct_list.append((atomname_shift_dict["CB"], atomname_shift_dict["CG2"], resname, atomname_shift_dict["SS"]))
    # except KeyError:
    #     pass
    #
    # # Save the CG, CD shifts to measure their correlation
    # try:
    #     if resname == "ILE":
    #         CG_CD_aa_SecStruct_list.append((atomname_shift_dict["CG1"], atomname_shift_dict["CD1"], resname, atomname_shift_dict["SS"]))
    #     elif resname == "LEU":
    #         if "CD1" in atomname_shift_dict.keys() and "CG" in atomname_shift_dict.keys():
    #             CG_CD_aa_SecStruct_list.append((atomname_shift_dict["CG"], atomname_shift_dict["CD1"], resname, atomname_shift_dict["SS"]))
    #         if "CD2" in atomname_shift_dict.keys() and "CG" in atomname_shift_dict.keys():
    #             CG_CD_aa_SecStruct_list.append((atomname_shift_dict["CG"], atomname_shift_dict["CD2"], resname, atomname_shift_dict["SS"]))
    #     else:
    #         CG_CD_aa_SecStruct_list.append((atomname_shift_dict["CG"], atomname_shift_dict["CD"], resname, atomname_shift_dict["SS"]))
    # except KeyError:
    #     pass
    #
    # # Save CO, CA, CB, SS, previous aa, current aa, next aa
    # if previous_resname != None:
    #     try:
    #         CO_CA_CB_SS_prevAA_curAA_nextAA_list.append((atomname_shift_dict["C"], atomname_shift_dict["CA"], atomname_shift_dict["CB"], atomname_shift_dict["SS"],
    #                                                  previous_resname, resname, next_resname))
    #     except KeyError:
    #         pass


def write_coordinates(raw_x, raw_y, fname):
    
    global HOME_DIR

    if not os.path.exists(HOME_DIR + "/Probability_Histograms"):
        os.mkdir(HOME_DIR + "/Probability_Histograms")

    with open(HOME_DIR + "/Probability_Histograms/" + fname, 'w') as f:
        for i in range(raw_x.shape[0]):
            f.write(str(raw_x[i])+"\t"+str(raw_y[i])+"\n")


def assign_ambiguous_degenerate_methylenes(resname, atomname_shift_dict):
    """
        FUNCTION for newVASCO to assign one resonance to the 1st and second resonance to the second equivalent proton (whenever applicable).
        These correspond to cases where the user did not find the second methylene proton and just copied the resonance of the
        first in the assignment file.
    """
    global aa_CHpairs2beMerged_dict

    # print "DEBUG assign_ambiguous_degenerate_methylenes: atomname_shift_dict=", atomname_shift_dict
    if resname == "LEU":
        if 'CD1' in list(atomname_shift_dict.keys()) and type(atomname_shift_dict['CD1']) == list:
            Cresons = atomname_shift_dict['CD1']
            if len(set(Cresons)) == 2:
                atomname_shift_dict['CD1'] = Cresons[0]
                atomname_shift_dict['CD2'] = Cresons[1]
            elif len(set(Cresons)) == 1:
                atomname_shift_dict['CD1'] = Cresons[0]
                if 'CD2' in list(atomname_shift_dict.keys()):
                    del atomname_shift_dict['CD2']
        elif 'CD2' in list(atomname_shift_dict.keys()) and type(atomname_shift_dict['CD2']) == list:
            Cresons = atomname_shift_dict['CD2']
            if len(set(Cresons)) == 2:
                atomname_shift_dict['CD1'] = Cresons[0]
                atomname_shift_dict['CD2'] = Cresons[1]
            elif len(set(Cresons)) == 1:
                atomname_shift_dict['CD1'] = Cresons[0]
                if 'CD2' in list(atomname_shift_dict.keys()):
                    del atomname_shift_dict['CD2']
        # TAKE CARE OF METHYL PROTONS
        methyl1_proton_reson_list = []
        for alt_Hname in ['HD11', 'HD12', 'HD13']:
            if alt_Hname in list(atomname_shift_dict.keys()):
                if type(atomname_shift_dict[alt_Hname]) == list:
                    methyl1_proton_reson_list.extend(atomname_shift_dict[alt_Hname])
                else:
                    methyl1_proton_reson_list.append(atomname_shift_dict[alt_Hname])
        if len(set(methyl1_proton_reson_list)) == 2:
            atomname_shift_dict['HD1'] = methyl1_proton_reson_list[0]
            atomname_shift_dict['HD2'] = methyl1_proton_reson_list[1]
        elif len(set(methyl1_proton_reson_list)) == 1:
            atomname_shift_dict['HD1'] = methyl1_proton_reson_list[0]
        # DELETE ALTERNATIVE PROTON NAMES FROM THE DICT
        for alt_Hname in ['HD11', 'HD12', 'HD13', 'HD21', 'HD22', 'HD23']:  # delete the alternative names from the dict
            if alt_Hname in list(atomname_shift_dict.keys()):
                del atomname_shift_dict[alt_Hname]
    elif resname == "VAL":
        if 'CG1' in list(atomname_shift_dict.keys()) and type(atomname_shift_dict['CG1']) == list:
            Cresons = atomname_shift_dict['CG1']
            if len(set(Cresons)) == 2:
                atomname_shift_dict['CG1'] = Cresons[0]
                atomname_shift_dict['CG2'] = Cresons[1]
            elif len(set(Cresons)) == 1:
                atomname_shift_dict['CG1'] = Cresons[0]
                if 'CG2' in list(atomname_shift_dict.keys()):
                    del atomname_shift_dict['CG2']
        elif 'CG2' in list(atomname_shift_dict.keys()) and type(atomname_shift_dict['CG2']) == list:
            Cresons = atomname_shift_dict['CG2']
            if len(set(Cresons)) == 2:
                atomname_shift_dict['CG1'] = Cresons[0]
                atomname_shift_dict['CG2'] = Cresons[1]
            elif len(set(Cresons)) == 1:
                atomname_shift_dict['CG1'] = Cresons[0]
                if 'CG2' in list(atomname_shift_dict.keys()):
                    del atomname_shift_dict['CG2']
        # TAKE CARE OF METHYL PROTONS
        methyl1_proton_reson_list = []
        for alt_Hname in ['HG11', 'HG12', 'HG13']:
            if alt_Hname in list(atomname_shift_dict.keys()):
                if type(atomname_shift_dict[alt_Hname]) == list:
                    methyl1_proton_reson_list.extend(atomname_shift_dict[alt_Hname])
                else:
                    methyl1_proton_reson_list.append(atomname_shift_dict[alt_Hname])
        if len(set(methyl1_proton_reson_list)) == 2:
            atomname_shift_dict['HG1'] = methyl1_proton_reson_list[0]
            atomname_shift_dict['HG2'] = methyl1_proton_reson_list[1]
        elif len(set(methyl1_proton_reson_list)) == 1:
            atomname_shift_dict['HG1'] = methyl1_proton_reson_list[0]
        # DELETE ALTERNATIVE PROTON NAMES FROM THE DICT
        for alt_Hname in ['HG11', 'HG12', 'HG13', 'HG21', 'HG22', 'HG23']:  # delete the alternative names from the dict
            if alt_Hname in list(atomname_shift_dict.keys()):
                del atomname_shift_dict[alt_Hname]

    # print "DEBUG: before ambiguous resname", resname, "atomname_shift_dict=", atomname_shift_dict
    try:
        for CHpair_list in aa_CHpairs2beMerged_dict[resname]:
            Cname = CHpair_list[0].split('-')[0]
            Hname1 = CHpair_list[0].split('-')[1]
            Hname2 = CHpair_list[1].split('-')[1]
            if Hname1 in list(atomname_shift_dict.keys()) and type(atomname_shift_dict[Hname1]) == list:
                # print "DEBUG: duplicating", atomname_shift_dict[Hname1], atomname_shift_dict[Hname2]
                cs = atomname_shift_dict[Hname1]
                if len(set(cs)) == 2:
                    atomname_shift_dict[Hname1] = cs[0]
                    atomname_shift_dict[Hname2] = cs[1]
                elif len(set(cs)) == 1:
                    atomname_shift_dict[Hname1] = cs[0]
                    if Hname2 in list(atomname_shift_dict.keys()):
                        del atomname_shift_dict[Hname2]
            elif Hname2 in list(atomname_shift_dict.keys()) and type(atomname_shift_dict[Hname2]) == list:
                cs = atomname_shift_dict[Hname2]
                if len(set(cs)) == 2:
                    atomname_shift_dict[Hname1] = cs[0]
                    atomname_shift_dict[Hname2] = cs[1]
                elif len(set(cs)) == 1:
                    atomname_shift_dict[Hname1] = cs[0]
                    if Hname2 in list(atomname_shift_dict.keys()):
                        del atomname_shift_dict[Hname2]
    except KeyError:
        return
    # print "DEBUG: after ambiguous resname", resname, "atomname_shift_dict=", atomname_shift_dict


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
# LOAD AND PARSE THE BMRB DATA
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
        #print "DEBUG: reading line:", line
        if (previous_resid != -1000 and resid != previous_resid) or (Entry_ID != previous_Entry_ID and previous_Entry_ID != None and Assigned_chem_shift_list_ID != previous_Assigned_chem_shift_list_ID):
            #print "DEBUG: saving resid=", previous_resid, "Entry_ID", previous_Entry_ID, "Assigned_chem_shift_list_ID=", previous_Assigned_chem_shift_list_ID
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
    
    print("Read ", counter, "lines and ", res_counter, " residues from file", args.BMRB_ENTRIES_FILE)

# TODO: this needs probably update. Look at process_VASCO.py
elif args.VASCO_DIR:
    VASCO_full_path = os.path.realpath(args.VASCO_DIR)
    print("Reading VASCO Database from ", VASCO_full_path)
    fnames=os.listdir(VASCO_full_path)
    fpattern = re.compile('bmr[0-9]+.[a-z0-9]{4}.vasco')
    VASCO_files_list = list(filter(fpattern.search, fnames))
    for VASCO_file in VASCO_files_list:
        print("DEBUG: reading file ", VASCO_file)
        mo = re.search('bmr([0-9]+).([a-z0-9]{4}).vasco', VASCO_file)
        if mo:
            Entry_ID = mo.group(1)
            with open(VASCO_full_path + "/" + VASCO_file, 'r') as f:
                for line in f:
                    #print "DEBUG: line=", line
                    word_list = line.split()
                    try:
                        if word_list[0] == "#":
                            continue
                    except IndexError:
                        continue
                    resid = word_list[2]
                    resname = word_list[3]
                    atomname = word_list[5]
                    if resname in list(VASCO_naming_exceptions_mdict.keys()) and atomname in list(VASCO_naming_exceptions_mdict[resname].keys()):
                        atomname = VASCO_naming_exceptions_mdict[resname][atomname]
                    if not resname in list(allowed_aa_atoms_dict.keys()):
                        print("WARNING: skipping line with unknown amino acid:", line)
                        continue
                    atomtype = word_list[7]
                    cs = word_list[8]
                    cs_stdev = word_list[9]
                    chain_ID = word_list[1]
                    ss = word_list[4]
                    #print "DEBUG: previous_resid", previous_resid, "resid", resid
                    if not int(previous_resid) in [int(resid)-1, int(resid)] or previous_Entry_ID != Entry_ID or previous_chain_ID != chain_ID:    # if this is a new protein, or there was a gap
                        protein_sequence = []           # reset the protein_sequence list
                    #print "DEBUG: reading line:", line
                    if (previous_resid != -1000 and resid != previous_resid) or (chain_ID != previous_chain_ID):
                        protein_sequence.append(resname)
                        try:
                            #print "DEBUG: protein_sequence = ", protein_sequence
                            previous2_resname = protein_sequence[-3]
                        except IndexError:
                            previous2_resname = None
                        #print "DEBUG: saving resid=", previous_resid, "chain_ID", previous_chain_ID
                        save_all_shift_pairs(previous_resname, atomname_shift_dict, previous2_resname, resname)
                        atomname_shift_dict = {}
                        res_counter += 1
                    atomname_shift_dict[atomname] = cs
                    #print "Secondary structure = ", ss
                    atomname_shift_dict["SS"] = ss
                    previous_resid = resid
                    previous_resname = resname
                    previous_Entry_ID = Entry_ID
                    previous_chain_ID = chain_ID
                    counter += 1
                
elif args.NEW_VASCO_DIR or args.NEW_VASCO_FILE:
    # LOAD AND PARSE NEW VASCO FILES
    if args.NEW_VASCO_FILE:
        Vfiles = [args.NEW_VASCO_FILE]
    else:
        Vfiles = list_files(args.NEW_VASCO_DIR, "bmr[0-9]+\.[0-9a-zA-Z]{4}.*\.cosh", full_path=True)
    for i,VASCO_FILE in enumerate(Vfiles):
        print("Reading new VASCO file #", i, ":", VASCO_FILE)
        residue_list = []  # list of residues in the order they appear in the VASCO file
        mo = re.search('bmr([0-9]+).([a-z0-9]{4})_[1-9].cosh', VASCO_FILE)
        if mo:
            Entry_ID = mo.group(1)
            with open(VASCO_FILE, 'r') as f:
                for line in f:
                    #print "DEBUG: line=", line
                    word_list = line.split()
                    try:
                        if line[0:4] != "ATOM":
                            continue
                    except IndexError:
                        continue
                    resid = line[22:26].rstrip()
                    resname = line[17:20].rstrip()
                    atomname = line[12:16].rstrip()
                    if resname in list(VASCO_naming_exceptions_mdict.keys()) and atomname in list(VASCO_naming_exceptions_mdict[resname].keys()):
                        atomname = VASCO_naming_exceptions_mdict[resname][atomname]
                    if not resname in list(allowed_aa_atoms_dict.keys()):
                        print("WARNING: skipping line with unknown amino acid:", line)
                        continue
                    atomtype = line[12:16].rstrip()
                    try:
                        cs = float(word_list[-1])
                    except ValueError:
                        continue
                    try:
                        cs2 = float(word_list[-2])
                        cs = [cs, cs2]
                    except ValueError:
                        pass
                    chain_ID = line[21]
                    ss = word_list[9]       # DANGEROUS: think of another way!
                    # print "DEBUG: previous_resid", previous_resid, "resid", resid
                    if not int(previous_resid) in [int(resid)-1, int(resid)] or previous_Entry_ID != Entry_ID or previous_chain_ID != chain_ID:    # if this is a new protein, or there was a gap
                        protein_sequence = []           # reset the protein_sequence list
                    # print "DEBUG: reading line:", line
                    if ((previous_resid != -1000 and resid != previous_resid) or
                        previous_resname != resname or (chain_ID != previous_chain_ID)):
                        protein_sequence.append(resname)
                        try:
                            # print "DEBUG: protein_sequence = ", protein_sequence
                            previous2_resname = protein_sequence[-3]
                        except IndexError:
                            previous2_resname = None
                        # print "DEBUG: saving resid=", previous_resid, "resname=", previous_resname, "chain_ID=", previous_chain_ID
                        if previous_resname != None:
                            previous_residue = aa3to1_dict[previous_resname]+str(previous_resid)
                            residue_list.append(previous_residue)
                            atomname_shift_dict = save_all_shift_pairs(previous_resname, atomname_shift_dict, previous2_resname, resname)
                            # print "DEBUG: writing atomname_shift_dict=", atomname_shift_dict
                        atomname_shift_dict = {}
                        res_counter += 1
                    atomname_shift_dict[atomname] = cs
                    #print "Secondary structure = ", ss
                    atomname_shift_dict["SS"] = ss
                    previous_resid = resid
                    previous_resname = resname
                    previous_Entry_ID = Entry_ID
                    previous_chain_ID = chain_ID
                    counter += 1

### FOR DEBUGGING
##for k1 in aa_CHpair_2Dhist_mdict.keys():
##    for k2 in aa_CHpair_2Dhist_mdict[k1].keys():
##        print k1, k2, aa_CHpair_2Dhist_mdict[k1][k2]


# # WRITE CO, CA, CB resonances to a file
# with open(HOME_DIR + "/BMRB_data/" + "CO_CA_CB_aa_SS_resonances.dat", 'w') as f:
#     for triplet in CO_CA_CB_aa_SecStruct_list:
#         f.write(triplet[0] + "\t" + triplet[1] + "\t" + triplet[2] + "\t" + triplet[3] + "\t" + triplet[4] + "\n")
#
# # WRITE CB, CG resonances to a file
# with open(HOME_DIR + "/BMRB_data/" + "CB_CG_aa_SS_resonances.dat", 'w') as f:
#     for triplet in CB_CG_aa_SecStruct_list:
#         f.write(triplet[0] + "\t" + triplet[1] + "\t" + triplet[2] + "\t" + triplet[3] + "\n")
#
# with open(HOME_DIR + "/BMRB_data/" + "CG_CD_aa_SS_resonances.dat", 'w') as f:
#     for triplet in CG_CD_aa_SecStruct_list:
#         f.write(triplet[0] + "\t" + triplet[1] + "\t" + triplet[2] + "\t" + triplet[3] + "\n")
#
# with open(HOME_DIR + "/BMRB_data/" + "CO_CA_CB_SS_prevAA_curAA_nextAA.dat", 'w') as f:
#     for group7 in CO_CA_CB_SS_prevAA_curAA_nextAA_list:
#         f.write(group7[0] + "\t" + group7[1] + "\t" + group7[2] + "\t" + group7[3] + "\t" + group7[4] + "\t" + group7[5] + "\t" + group7[6] + "\n")

#sys.exit(1)

# Copy the statistics of HE2 of MET to HE2, because in VASCO both protons are named "HE*" and we replaced "HE*" with "HE2" already
# aa_CHpair_2Dhist_mdict["MET"]["CE-HE3"] = aa_CHpair_2Dhist_mdict["MET"]["CE-HE2"]
# aa_CHpair_2Dhist_mdict["MET"]["CE-HE1"] = aa_CHpair_2Dhist_mdict["MET"]["CE-HE2"]

# BEFORE YOU CREATE ANY HISTOGRAMS, MERGE THE CHEMICAL SHIFT COORDINATES OF THE CH2 PAIRS IN ORDER TO HAVE
# IDENTICAL HISTOGRAMS
if args.MERGE_EQUIVALENT_CH_PAIRS == True:
    for resname in list(aa_CHpairs2beMerged_dict.keys()):
        for CH_pairs in aa_CHpairs2beMerged_dict[resname]:
            CH_pair1 = CH_pairs[0]
            CH_pair2 = CH_pairs[1]
            # save the resonances into separate lists because when you extend the 1st then you will add the the values of the 2nd along with those of the 1st
            H_list1 = aa_CHpair_2Dhist_mdict[resname][CH_pair1][1]
            C_list1 = aa_CHpair_2Dhist_mdict[resname][CH_pair1][0]
            H_list2 = aa_CHpair_2Dhist_mdict[resname][CH_pair2][1]
            C_list2 = aa_CHpair_2Dhist_mdict[resname][CH_pair2][0]
            if not (resname == "LEU" and CH_pair1 in ['CD1-HD1', 'CD2-HD2'] or resname == "VAL" and CH_pair1 in ["CG1-HG1", "CG2-HG2"]) and (
                len(H_list1) != len(H_list2) or len(C_list1) != len(C_list2) ):
                print("ERROR: the C and H histograms of ", resname, CH_pairs, "do not have equal sizes!!!")
                print("DEBUG: len(H_list1)=", len(H_list1), "len(H_list2)=", len(H_list2), "len(C_list1)=", len(C_list1), "len(C_list2)=", len(C_list2))
                sys.exit(1)
            elif resname == "LEU" and CH_pair1 in ['CD1-HD1', 'CD2-HD2'] or resname == "VAL" and CH_pair1 in ["CG1-HG1", "CG2-HG2"]:
                aa_CHpair_2Dhist_mdict[resname][CH_pair1][1] = H_list1 + H_list2
                aa_CHpair_2Dhist_mdict[resname][CH_pair2][1] = H_list1 + H_list2
                # append also the Carbon resonances, otherwise C and H will not have the same length to make a 2Dhist
                aa_CHpair_2Dhist_mdict[resname][CH_pair1][0] = C_list1 + C_list2
                aa_CHpair_2Dhist_mdict[resname][CH_pair2][0] = C_list1 + C_list2
            else:
                for index in range(len(H_list1)):
                    if H_list1[index] != H_list2[index]:
                        aa_CHpair_2Dhist_mdict[resname][CH_pair1][1].append(H_list2[index])
                        aa_CHpair_2Dhist_mdict[resname][CH_pair2][1].append(H_list1[index])
                        # append also the Carbon resonances, otherwise C and H will not have the same length to make a 2Dhist
                        aa_CHpair_2Dhist_mdict[resname][CH_pair1][0].append(C_list2[index])
                        aa_CHpair_2Dhist_mdict[resname][CH_pair2][0].append(C_list1[index])

# CREATE THE CORRELATED 2D C-H HISTOGRAMS
for aa in list(aa_CHpair_2Dhist_mdict.keys()):
    # if aa in ["PHE"]:
    for CH_pair in list(aa_CHpair_2Dhist_mdict[aa].keys()):
        # if aa +"_"+ CH_pair in ["PHE_CD1-HD1"]:
        # make 2D histogram for this C-H pair
        print("DEBUG: aa=", aa, "CH_pair=", CH_pair, "len(aa_CHpair_2Dhist_mdict[aa][CH_pair])=", len(aa_CHpair_2Dhist_mdict[aa][CH_pair][1]))
        print("DEBUG: aa=", aa, "CH_pair=", CH_pair, "np.mean(aa_CHpair_2Dhist_mdict[aa][CH_pair])=", np.mean(aa_CHpair_2Dhist_mdict[aa][CH_pair][1]))
        raw_x = np.array(aa_CHpair_2Dhist_mdict[aa][CH_pair][1])  # hydrogen CS (x-axis coordinates)
        raw_y = np.array(aa_CHpair_2Dhist_mdict[aa][CH_pair][0])  # carbon CS (y-axis coordinates)
        write_coordinates(raw_x, raw_y, aa+"_"+CH_pair+"_correlated_2D_coords.txt")
        x, y = discard_outliers(raw_x, raw_y)  # discard outliers from C and H chemical shift arrays
        ColorPrint("Writing file "+ aa+"_"+CH_pair+"_correlated_2D_coords.no_outliers.txt", "BOLDGREEN")
        write_coordinates(x, y, aa+"_"+CH_pair+"_correlated_2D_coords.no_outliers.txt")
        #x,y = raw_x, raw_y
        #continue
        #print "DEBUG: raw_x=", raw_x.tolist()
        # print "DEBUG: x=", x.tolist()
        xmin = min(x)
        xmax = max(x)
        #print "DEBUG: raw_y=", raw_y.tolist()
        # print "DEBUG: y=", y.tolist()
        ymin = min(y)
        ymax = max(y)
        adjusted_Cedges_list = [e for e in Cedges_list if e > ymin-0.4 and e < ymax+0.4]    # truncated the full-length edge list to fit this C value distribution but leave 0.4 ppm extra space at the ends
        adjusted_Hedges_list = [e for e in Hedges_list if e > xmin-0.08 and e < xmax+0.08]
        #print "DEBUG: Cedges_list=", Cedges_list
        #print "DEBUG: Hedges_list=", Hedges_list
        H, xedges, yedges = np.histogram2d(x, y, bins=(adjusted_Hedges_list, adjusted_Cedges_list))
        # print "DEBUG: xedges=", xedges.tolist()
        # print "DEBUG: yedges=", yedges.tolist()
        Hnorm = H/float(H.sum())   # normalize the histogram to get probabilities
        (xdim,ydim) = Hnorm.shape
        #print "DEBUG: Hnorm.shape", Hnorm.shape
        #print "DEBUG: len(xedges)=", len(xedges), "len(yedges)=", len(yedges)
        #np.set_printoptions(threshold='nan')
        #print "DEBUG: Hnorm="
        #print(Hnorm)
        # save the normalized 2D histogram into a file (by convention x-axis must be H and y-axis C)
        ColorPrint("Writing file "+ aa+"_"+CH_pair+"_correlated_2Dhist.txt", "BOLDGREEN")
        with open(HOME_DIR + "/Probability_Histograms/" + aa+"_"+CH_pair+"_correlated_2Dhist.txt", 'w') as f:
           for i in range(xdim):
               for j in range(ydim):
                   f.write("%.2f" % xedges[i] + "\t" + "%.2f" % yedges[j] + "\t" + str(Hnorm[i][j]) + "\n")
        # SMOOTHEN THE PROBABILITY DISTRIBUTION AND SAVE THE RESPECTIVE HISTOGRAM IN A DIFFERENT FILE
        Hsmooth = histcalc.smoothen_distribution(x, y, xedges, yedges)
        #print "DEBUG: Hsmooth.shape", Hsmooth.shape
        #print "DEBUG: Hsmooth="
        #print(Hsmooth)
        (xdim,ydim) = Hsmooth.shape
        ColorPrint("Writing file "+ aa+"_"+CH_pair+"_correlated_2Dhist.smoothed.txt", "BOLDGREEN")
        with open(HOME_DIR + "/Probability_Histograms/" + aa+"_"+CH_pair+"_correlated_2Dhist.smoothed.txt", 'w') as f:
            for i in range(xdim):
                for j in range(ydim):
                    f.write("%.2f" %  xedges[i] + "\t" + "%.2f" % yedges[j] + "\t" + str(Hsmooth[i][j]) + "\n")

sys.exit(1)
# CREATE THE UNCORRELATED 2D C-H HISTOGRAMS
for aa in list(aa_CHpair_2Dhist_mdict.keys()):
    # if not aa in ['LEU', 'VAL']:
    #     continue
    for CH_pair in list(aa_CHpair_2Dhist_mdict[aa].keys()):
        # make 2D histogram for this C-H pair
        print("DEBUG: writing 1D and 2D uncorrelated histograms for ", aa, CH_pair)
        raw_x = np.array(aa_CHpair_2Dhist_mdict[aa][CH_pair][1])  # hydrogen CS (x-axis coordinates)
        raw_y = np.array(aa_CHpair_2Dhist_mdict[aa][CH_pair][0])  # carbon CS (y-axis coordinates)
        x, y = discard_outliers(raw_x, raw_y)  # discard outliers from C and H chemical shift arrays
        #x, y = raw_x, raw_y
        print("DEBUG: x=", x)
        xmin = min(x)
        xmax = max(x)
        print("DEBUG: y=", y)
        ymin = min(y)
        ymax = max(y)
        adjusted_Cedges_list = [e for e in Cedges_list if e > ymin-0.4 and e < ymax+0.4]
        adjusted_Hedges_list = [e for e in Hedges_list if e > xmin-0.08 and e < xmax+0.08]
        #print "DEBUG: Cedges_list=", Cedges_list
        #print "DEBUG: Hedges_list=", Hedges_list
        Hh, xedges = np.histogram(x, bins=adjusted_Hedges_list, density=False)  # make 1D histogram for Hydrogen
        Hh_norm = Hh/float(Hh.sum())   # normalize the histogram to get probabilities
        print("DEBUG: Hh=", Hh)
        print("DEBUG: xedges=", xedges.tolist())
        xdim = Hh_norm.shape[0]
        print("DEBUG: xdim=", xdim)
        print("DEBUG: Hh_norm=", Hh_norm.tolist())
        Hhist_1D = []
        for i in range(xdim):
            Hhist_1D.append((round(xedges[i], 2), Hh_norm[i]) )
        print("DEBUG: Hhist_1D=", Hhist_1D)
        xedge = Hhist_1D[0][0]-args.H_bin_length
        while round(xedge, 2) >= -9.98:
            Hhist_1D.insert(0, (round(xedge, 2), 0.00))
            xedge -= args.H_bin_length
        xedge = Hhist_1D[-1][0]+args.H_bin_length
        while round(xedge, 2) <= 33.98:
            Hhist_1D.append((round(xedge, 2), 0.00))
            xedge += args.H_bin_length
        print("Writing file "+ aa+"_"+CH_pair.split('-')[1]+"_hist.txt")
        with open(HOME_DIR + "/Probability_Histograms/" + aa+"_"+CH_pair.split('-')[1]+"_hist.txt", 'w') as f:
            for i in range(len(Hhist_1D)):
                f.write(str(Hhist_1D[i][0]) + "\t" + str(Hhist_1D[i][1]) + "\n")
        Hc, yedges = np.histogram(y, bins=adjusted_Cedges_list, density=False)  # make 1D histogram for Carbon
        Hc_norm = Hc/float(Hc.sum())   # normalize the histogram to get probabilities
        print("DEBUG: Hc=", Hc)
        ydim = Hc_norm.shape[0]
        print("DEBUG: ydim=", ydim)
        print("DEBUG: yedges=", yedges.tolist())
        print("DEBUG: Hc_norm=", Hc_norm.tolist())
        Hcist_1D = []
        for i in range(ydim):
            Hcist_1D.append((round(yedges[i], 2), Hc_norm[i]) )
        print("DEBUG: Hcist_1D=", Hcist_1D)
        yedge = Hcist_1D[0][0]-args.C_bin_length
        while round(yedge, 2) >= 0.1:
            Hcist_1D.insert(0, (round(yedge, 2), 0.00))
            yedge -= args.C_bin_length
        yedge = Hcist_1D[-1][0]+args.C_bin_length
        while round(yedge, 2) <= 199.9:
            Hcist_1D.append((round(yedge, 2), 0.00))
            yedge += args.C_bin_length
        print("Writing file "+ aa+"_"+CH_pair.split('-')[0]+"_hist.txt")
        with open(HOME_DIR + "/Probability_Histograms/" + aa+"_"+CH_pair.split('-')[0]+"_hist.txt", 'w') as f:
            for i in range(len(Hcist_1D)):
                f.write(str(Hcist_1D[i][0]) + "\t" + str(Hcist_1D[i][1]) + "\n")
        # create the 2D histogram from two 1D histograms
        weighted_H = np.zeros([xdim, ydim])
        H = np.zeros([xdim, ydim])
        for i in range(xdim):
                for j in range(ydim):
                    weighted_uncorr_prob = (args.H_weight * Hh_norm[i] + args.C_weight * Hc_norm[j])/(args.H_weight + args.C_weight)
                    weighted_H[i][j] = weighted_uncorr_prob
                    # uncorr_prob = (Hh_norm[i] + Hc_norm[j])/2.0
                    uncorr_prob = Hh_norm[i] * Hc_norm[j]
                    H[i][j] = uncorr_prob
        weighted_Hnorm = weighted_H/float(weighted_H.sum())   # normalize the histogram to get probabilities
        Hnorm = H/float(H.sum())   # normalize the histogram to get probabilities
        # save the normalized 2D histogram into a file (by convention x-axis must be H and y-axis C)
        print("Writing file "+ aa+"_"+CH_pair+"_weighted_uncorrelated_2Dhist.txt")
        with open(HOME_DIR + "/Probability_Histograms/" + aa+"_"+CH_pair+"_weighted_uncorrelated_2Dhist.txt", 'w') as f:
            for i in range(xdim):
                for j in range(ydim):
                    f.write(str(xedges[i]) + "\t" + str(yedges[j]) + "\t" + str(weighted_Hnorm[i][j]) + "\n")
        print("Writing file "+ aa+"_"+CH_pair+"_uncorrelated_2Dhist.txt")
        with open(HOME_DIR + "/Probability_Histograms/" + aa+"_"+CH_pair+"_uncorrelated_2Dhist.txt", 'w') as f:
            for i in range(xdim):
                for j in range(ydim):
                    f.write(str(xedges[i]) + "\t" + str(yedges[j]) + "\t" + str(Hnorm[i][j]) + "\n")



## CONVERT THE DATA INTO GNUPLOT FORMAT
#for fname in $(ls *correlated_2Dhist.smoothed.txt)
#do
# awk 'BEGIN{x='$(head -1 $fname | awk '{print $1}')'}{if($1 != x){print ""; x=$1};print $0}' $fname > ${fname}_for_GNUPLOT
#done
## TO SAVE ALL THE 2D HISTOGRAMS AUTOMATICALLY
#for fname in $(ls *correlated_2Dhist.smoothed.txt);
#do
#name=$(perl -pi -e "s/\.txt//" <<< $fname)
##perl -pi -e "s/[A-Z_0-9-]+correlated_2Dhist/${name}/" plot_2D_histogram.gp
#perl -pi -e "s/basename=\".+correlated_2Dhist\"/basename=\"${name}\"/" plot_2D_histogram.gp
#gnuplot -e "basename='$name'" plot_2D_histogram.gp
#done