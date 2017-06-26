#!/usr/bin/env python2.7
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

HOME_DIR = os.path.dirname(os.path.realpath(__file__))
code3 = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
code1 = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
aa1to3_dict = dict((c1,c3) for c1,c3 in zip(code1,code3))
aa3to1_dict = dict((c3,c1) for c3,c1 in zip(code3,code1))

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

def cmdlineparse():
    parser = ArgumentParser(description="command line arguments",
        epilog="EXAMPLE: ./VASCO2sparky_format.py -vascofile BMRB_data/VASCO/vasco/bmr7114.2klf.vasco")
    parser.add_argument("-vascofile", dest="VASCO_FILE", required=False, type=str, help="vasco file",
                        metavar="<BMRB entries file>")
    parser.add_argument("-outfile", dest="OUT_FILE", required=False, type=str, default=None, help="output file name (default input file name with .sparky extension)",
                        metavar="<output file name>")
    args=parser.parse_args()
    return args


args = cmdlineparse()
if args.OUT_FILE == None:
    args.OUT_FILE = args.VASCO_FILE.replace(".vasco", ".sparky")

def save_all_shift_pairs(resname, raw_atomname_shift_dict, previous_resname, next_resname):
    """
        FUNCTION to save all C-H shift pairs of a particular residue. Basically updates aa_CHpair_2Dhist_multidict.
    """
    global aa_CHpair_2Dhist_multidict, CO_CA_CB_list
    saved_alternative_Hnames_list = []  # list with the alternative Hnames that were saved; with the purpose to avoid saving them multiple times
    
    for CH_pair in aa_CHpair_2Dhist_multidict[resname].keys():
        [Cname, Hname] = CH_pair.split("-")
        if Cname in atomname_shift_dict.keys() and Hname in atomname_shift_dict.keys():
            aa_CHpair_2Dhist_multidict[resname][Cname+"-"+Hname][0].append(round(float(atomname_shift_dict[Cname]), 4))
            aa_CHpair_2Dhist_multidict[resname][Cname+"-"+Hname][1].append(round(float(atomname_shift_dict[Hname]), 4))

        elif resname in aa_alternativeAtomNames_dict.keys() and Cname in aa_alternativeAtomNames_dict[resname].keys() and Hname in aa_alternativeAtomNames_dict[resname][Cname]:
            for alt_Hname in aa_alternativeAtomNames_dict[resname][Cname][1:]:  # skip the first which is the conventional Hname
                if Cname in atomname_shift_dict.keys() and alt_Hname in atomname_shift_dict.keys():
                    conv_Hname = aa_alternativeAtomNames_dict[resname][Cname][0] # replace the alternative Hname with the conventional Hname
                    if conv_Hname in saved_alternative_Hnames_list: # if we have saved already this conv_Hname with an alternative name, skip this
                        continue
                    saved_alternative_Hnames_list.append(conv_Hname)
                    aa_CHpair_2Dhist_multidict[resname][Cname+"-"+conv_Hname][0].append(round(float(atomname_shift_dict[Cname]), 4))
                    aa_CHpair_2Dhist_multidict[resname][Cname+"-"+conv_Hname][1].append(round(float(atomname_shift_dict[alt_Hname]), 4))
                    break # it is uncessessary to check the rest of the alternative Hnames since we have already save one

def save_residue_shifts(resname, resid, atomname_shift_dict):
    """
    write each residue shifts to the following format:
    X500HA-CA-N-H                number                number                number                number
    """
    print "DEBUG: resname=", resname, "resid=", resid, "atomname_shift_dict=", atomname_shift_dict
    global aa_CHpair_2Dhist_multidict, aa_alternativeAtomNames_dict, residue_HCpair_N_H_resonance_multidict
    saved_alternative_Hnames_list = []  # list with the alternative Hnames that were saved; with the purpose to avoid saving them multiple times
    residue = aa3to1_dict[resname]+str(resid)
    print "DEBUG: atomname_shift_dict=", atomname_shift_dict
    
    for CH_pair in aa_CHpair_2Dhist_multidict[resname].keys():
        [Cname, Hname] = CH_pair.split("-")
        if Cname in atomname_shift_dict.keys() and Hname in atomname_shift_dict.keys():
            residue_HCpair_N_H_resonance_multidict[residue][Hname+"-"+Cname] = [atomname_shift_dict[Hname], atomname_shift_dict[Cname]]
            
        elif resname in aa_alternativeAtomNames_dict.keys() and Cname in aa_alternativeAtomNames_dict[resname].keys() and Hname in aa_alternativeAtomNames_dict[resname][Cname]:
            for alt_Hname in aa_alternativeAtomNames_dict[resname][Cname][1:]:  # skip the first which is the conventional Hname
                if Cname in atomname_shift_dict.keys() and alt_Hname in atomname_shift_dict.keys():
                    conv_Hname = aa_alternativeAtomNames_dict[resname][Cname][0] # replace the alternative Hname with the conventional Hname
                    if conv_Hname in saved_alternative_Hnames_list: # if we have saved already this conv_Hname with an alternative name, skip this
                        continue
                    saved_alternative_Hnames_list.append(conv_Hname)
                    residue_HCpair_N_H_resonance_multidict[residue][conv_Hname+"-"+Cname] = [atomname_shift_dict[alt_Hname], atomname_shift_dict[Cname]]
                    break # it is uncessessary to check the rest of the alternative Hnames since we have already save one
    if "N" in atomname_shift_dict.keys():
        residue_HCpair_N_H_resonance_multidict[residue]["N"] = atomname_shift_dict["N"]
    if "H" in atomname_shift_dict.keys():
        residue_HCpair_N_H_resonance_multidict[residue]["H"] = atomname_shift_dict["H"]
    

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
residue_HCpair_N_H_resonance_multidict = tree() # residue -> atom name -> resonance
residue_list = [] # list of residues in the order they appear in the VASCO file
print "Reading VASCO file ", args.VASCO_FILE
mo = re.search('bmr([0-9]+).([a-z0-9]{4}).vasco', args.VASCO_FILE)
if mo:
    Entry_ID = mo.group(1)
    with open(args.VASCO_FILE, 'r') as f:
        for line in f:
            word_list = line.split()
            try:
                if word_list[0] in ["", "#"]:
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
            print "DEBUG: previous_resid", previous_resid, "resid", resid
            if not int(previous_resid) in [int(resid)-1, int(resid)] or previous_Entry_ID != Entry_ID or previous_chain_ID != chain_ID:    # if this is a new protein, or there was a gap
                protein_sequence = []           # reset the protein_sequence list
            print "DEBUG: reading line:", line
            if (previous_resid != -1000 and resid != previous_resid) or (chain_ID != previous_chain_ID):
                protein_sequence.append(resname)
                try:
                    print "DEBUG: protein_sequence = ", protein_sequence
                    previous2_resname = protein_sequence[-3]
                except IndexError:
                    previous2_resname = None
                print "DEBUG: saving resid=", previous_resid, "chain_ID", previous_chain_ID
                if previous_resname != None:
                    previous_residue = aa3to1_dict[previous_resname]+str(previous_resid)
                    residue_list.append(previous_residue)
                    save_residue_shifts(previous_resname, previous_resid, atomname_shift_dict)
                atomname_shift_dict = {}
                res_counter += 1
            atomname_shift_dict[atomname] = cs
            atomname_shift_dict["SS"] = ss
            previous_resid = resid
            previous_resname = resname
            previous_Entry_ID = Entry_ID
            previous_chain_ID = chain_ID
            counter += 1

outfile = open(args.OUT_FILE, 'w')
for residue in residue_list:
    if not "N" in residue_HCpair_N_H_resonance_multidict[residue].keys() or not "H" in residue_HCpair_N_H_resonance_multidict[residue].keys():
        continue
    HC_pairs_list = [k for k in residue_HCpair_N_H_resonance_multidict[residue].keys() if not k in ["N", "H"]]
    for HC_pair in HC_pairs_list:
        outfile.write(residue+HC_pair+"-N-H\t"+str(residue_HCpair_N_H_resonance_multidict[residue][HC_pair][0])+"\t"+\
            str(residue_HCpair_N_H_resonance_multidict[residue][HC_pair][1])+"\t"+\
            str(residue_HCpair_N_H_resonance_multidict[residue]["N"])+"\t"+\
            str(residue_HCpair_N_H_resonance_multidict[residue]["H"])+"\n")
outfile.close()