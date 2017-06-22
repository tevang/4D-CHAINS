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

import sys, re, os
import numpy as np
from ordereddict import OrderedDict
from argparse import ArgumentParser


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
"MET" : ["HA", "HB2", "HB3", "HG2", "HG3", "HE", "CA", "CB", "CG", "CE", "N", "H"],
"PHE" : ["HA", "HB2", "HB3", "CA", "CB", "N", "H"],
"PRO" : ["HA", "HB2", "HB3", "HG2", "HG3", "HD2", "HD3", "CA", "CB", "CG", "CD", "N"], # Prolines are not detected by the method at position "i" due to lack of HN hydrogen
"SER" : ["HA", "HB2", "HB3", "CA", "CB", "N", "H"],
"THR" : ["HA", "HB", "HG2", "CA", "CB", "CG2", "N", "H"],
"TRP" : ["HA", "HB2", "HB3", "CA", "CB", "N", "H"],
"TYR" : ["HA", "HB2", "HB3", "CA", "CB", "N", "H"],
"VAL" : ["HA", "HB", "HG1", "HG2", "CA", "CB", "CG1", "CG2", "N", "H"]
}


NMR_detectable_heavy_atoms = {
"ALA" : ["N", "CA", "CB"],
"ARG" : ["N", "CA", "CB", "CG", "CD"],
"ASN" : ["N", "CA", "CB", "ND2"],
"ASP" : ["N", "CA", "CB"],
"CYS" : ["N", "CA", "CB"],
"GLN" : ["N", "CA", "CB", "CG", "NE2"],
"GLU" : ["N", "CA", "CB", "CG"],
"GLY" : ["N", "CA"],
"HIS" : ["N", "CA", "CB", "CD2", "CE1"],
"ILE" : ["N", "CA", "CB", "CG1", "CG2", "CD1"],
"LEU" : ["N", "CA", "CB", "CG", "CD1", "CD2"],
"LYS" : ["N", "CA", "CB", "CG", "CD", "CE"],
"MET" : ["N", "CA", "CB", "CG", "CE"],
"PHE" : ["N", "CA", "CB", "CD", "CE", "CZ"],
"PRO" : ["CA", "CB", "CG", "CD"],
"SER" : ["N", "CA", "CB"],
"THR" : ["N", "CA", "CB", "CG2"],
"TRP" : ["N", "CA", "CB", "CD1", "NE1", "CZ2", "CH2", "CZ3", "CE3"],
"TYR" : ["N", "CA", "CB", "CD", "CE"],
"VAL" : ["N", "CA", "CB", "CG1", "CG2"]
}


methyl_carbons_dict = {}   # resname->list of allowed carbons in iteration1
methyl_carbons_dict["LEU"] = ["CD1", "CD2"]
methyl_carbons_dict["VAL"] = ["CG1", "CG2"]
methyl_carbons_dict["ILE"] = ["CD1", "CG2"]
methyl_carbons_dict["THR"] = ["CG2"]
methyl_carbons_dict["ALA"] = ["CB"]
methyl_carbons_dict["MET"] = ["CE"]


def cmdlineparse():
    parser = ArgumentParser(description="command line arguments",
        epilog="EXAMPLE: calculate_assignment_statistics.py -asgnfile results_summary.3mers_round1_rst.chainlinkers.xeasy -tseq ../Nab3.fasta \\\
               or to calculate statistics from multiple assignment files: ")
    parser.add_argument("-asgnfile", dest="ASSIGNMENT_FILE", required=False, action='append',
                        help="file with chemical shift assignments in xeasy format (may be TOCSY or NOESY)", metavar="<xeasy assignment file>")
    parser.add_argument("-tseq", dest="template_sequence_file", required=True, action='append',
                        help="protein sequences file in fasta format", metavar="<protein sequences file>")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s version 0.9')
    args=parser.parse_args()
    return args


args = cmdlineparse()

def get_resname_from_resid(resid, start):
    global template_sequence_list
    
    resname = aa1to3_dict[template_sequence_list[int(resid)-int(start)]]
    return resname
    

TOTAL_expected_number_of_carbons = 0
TOTAL_num_of_assigned_carbons_TOCSY = 0
TOTAL_num_of_correct_carbons_TOCSY = 0
TOTAL_num_of_wrong_carbons_TOCSY = 0

TOTAL_num_of_assigned_carbons_NOESY = 0
TOTAL_num_of_correct_carbons_NOESY = 0
TOTAL_num_of_wrong_carbons_NOESY = 0

for ASSIGNMENT_FILE,template_sequence_file in zip(args.ASSIGNMENT_FILE, args.template_sequence_file):
    print "\n\n############################# Analyzing files ", ASSIGNMENT_FILE," and ", template_sequence_file, " ################################\n\n"
    
    template_sequence_string = ""
    with open(template_sequence_file, 'r') as f:
        for line in f:
            if not re.match("^>", line):
                template_sequence_string += re.sub(r"[^A-Z]", "", line)
    template_sequence_list = list(template_sequence_string)
    
    num_of_allowed_heavy_atoms = 0  # number of heavy atoms that can be detected (C + N, including aromatic C+N and ASN, GLN sidechain N)
    expected_number_of_methyl_carbons = 0
    expected_number_of_nonmethyl_carbons = 0
    num_of_allowed_aliphatic_carbons = 0
    num_of_allowed_all_atoms = 0    # number of all atoms that can be detected (C+H+N)
    num_of_backbone_heavy_atoms = 0
    total_num_of_NMR_detectable = 0 # total number of NMR detectable heavy atoms in the protein (C+N, including aromatic C and N)
    aatype_occurence_dict = {}
    for resname in allowed_aa_atoms_dict.keys():
        resnum = template_sequence_list.count(aa3to1_dict[resname])  # number of occurence of this aa type in the protein
        aatype_occurence_dict[resname] = resnum
        heavy_atoms_list = [atom for atom in allowed_aa_atoms_dict[resname] if atom[0]!='H']
        aliphatic_carbons_list = [atom for atom in allowed_aa_atoms_dict[resname] if atom[0]!='H' and atom[0]!='N']
        num_of_allowed_heavy_atoms += resnum * len(heavy_atoms_list)
        num_of_allowed_aliphatic_carbons += resnum * len(aliphatic_carbons_list)
        num_of_allowed_all_atoms += resnum * len(allowed_aa_atoms_dict[resname])
        total_num_of_NMR_detectable += resnum * len(NMR_detectable_heavy_atoms[resname])
        if resname in methyl_carbons_dict.keys():
            expected_number_of_methyl_carbons += resnum * len(methyl_carbons_dict[resname])
            expected_number_of_nonmethyl_carbons += resnum * len(aliphatic_carbons_list) -1*resnum * len(methyl_carbons_dict[resname])
        else:
            expected_number_of_nonmethyl_carbons += resnum * len(aliphatic_carbons_list)
    
    
    num_of_assigned_allowed_heavy_atoms = 0  # number of heavy atoms detected by this 4D methods that were assigned
    num_of_assigned_methyl_carbons_TOCSY, num_of_assigned_methyl_carbons_NOESY = 0, 0  # number of methyl carbons detected by this 4D methods that were assigned
    num_of_correct_methyl_carbons_TOCSY, num_of_wrong_methyl_carbons_TOCSY, num_of_correct_methyl_carbons_NOESY, num_of_wrong_methyl_carbons_NOESY = 0, 0, 0, 0
    num_of_assigned_nonmethyl_carbons_TOCSY, num_of_assigned_nonmethyl_carbons_NOESY = 0, 0  # number of non-methyl carbons detected by this 4D methods that were assigned
    num_of_correct_nonmethyl_carbons_TOCSY, num_of_wrong_nonmethyl_carbons_TOCSY, num_of_correct_nonmethyl_carbons_NOESY, num_of_wrong_nonmethyl_carbons_NOESY = 0, 0, 0, 0
    num_of_assigned_NMR_detectable_atoms = 0    # number of all NMR detectable heavy atoms that were assigned
    num_of_assigned_N = 0 # number of backbone N that were assigned
    residue_assignedAtoms_dict = {} # dict of the form residue (resname + resid) -> list of assigned heavy atoms
    with open(ASSIGNMENT_FILE, 'r') as f:
        contents = f.readlines()
        start = int(contents[0].split()[4])
        for line in contents:
            word_list = line.split()
            if len(word_list) == 5: # assume this is a manual assignment file
                line = line.strip() + " # XXX TOCSY-NOESY <CORRECT>"
            resid = word_list[4]
            try:
                resname = word_list[6]
            except IndexError:
                resname = get_resname_from_resid(resid, start)
            atom = word_list[3]
            if not atom in allowed_aa_atoms_dict[resname]:  # to exclude atomatic atoms
                continue
            residue = aa3to1_dict[resname] + str(resid)
            try:
                residue_assignedAtoms_dict[residue].append(atom)
            except (ValueError, KeyError):
                residue_assignedAtoms_dict[residue] = [atom]
            if atom[0] != 'H' and atom in allowed_aa_atoms_dict[resname]:
                num_of_assigned_allowed_heavy_atoms += 1
            if atom[0] != 'H' and atom in NMR_detectable_heavy_atoms[resname]:
                num_of_assigned_NMR_detectable_atoms += 1
            if atom == 'N':
                num_of_assigned_N += 1
            if resname in methyl_carbons_dict.keys() and atom in methyl_carbons_dict[resname]:
                if "TOCSY" in line:
                    num_of_assigned_methyl_carbons_TOCSY += 1
                    if "<CORRECT>" in line:
                        num_of_correct_methyl_carbons_TOCSY += 1
                    elif "<WRONG>" in line:
                        num_of_wrong_methyl_carbons_TOCSY += 1
                    elif "<CHECK MANUALLY>" in line or "<NOT ASSIGNED>" in line:
                        print "Please check manually if the following methyl Carbon assignment is correct:"
                        print line
                elif "NOESY" in line:
                    num_of_assigned_methyl_carbons_NOESY += 1
                    if "<CORRECT>" in line:
                        num_of_correct_methyl_carbons_NOESY += 1
                    elif "<WRONG>" in line:
                        num_of_wrong_methyl_carbons_NOESY += 1
                    elif "<CHECK MANUALLY>" in line or "<NOT ASSIGNED>" in line:
                        print "Please check manually if the following methyl Carbon assignment is correct:"
                        print line
            if atom[0] == 'C' and not (resname in methyl_carbons_dict.keys() and atom in methyl_carbons_dict[resname]):
                if "TOCSY" in line:
                    num_of_assigned_nonmethyl_carbons_TOCSY += 1
                    if "<CORRECT>" in line:
                        num_of_correct_nonmethyl_carbons_TOCSY += 1
                    elif "<WRONG>" in line:
                        num_of_wrong_nonmethyl_carbons_TOCSY += 1
                    elif "<CHECK MANUALLY>" in line or "<NOT ASSIGNED>" in line:
                        print "Please check manually if the following non-methyl Carbon assignment is correct:"
                        print line
                elif "NOESY" in line:
                    num_of_assigned_nonmethyl_carbons_NOESY += 1
                    if "<CORRECT>" in line:
                        num_of_correct_nonmethyl_carbons_NOESY += 1
                    elif "<WRONG>" in line:
                        num_of_wrong_nonmethyl_carbons_NOESY += 1
                    elif "<CHECK MANUALLY>" in line or "<NOT ASSIGNED>" in line:
                        print "Please check manually if the following non-methyl Carbon assignment is correct:"
                        print line
    
    total_num_of_N = len(template_sequence_list)-template_sequence_list.count('P')-1 # minus prolines, minus the N-term residue
    
    
    print "STATISTICS:"
    print "ALL NMR DETECTABLE HEAVY ATOMS ASSIGNED:", str(num_of_assigned_NMR_detectable_atoms)+"/"+str(total_num_of_NMR_detectable), "("+ str(100*float(num_of_assigned_NMR_detectable_atoms)/total_num_of_NMR_detectable), "%)"
    print "ALLOWED HEAVY ATOMS ASSIGNED: "+str(num_of_assigned_allowed_heavy_atoms)+"/"+str(num_of_allowed_heavy_atoms)+" (" + \
            str(100*float(num_of_assigned_allowed_heavy_atoms)/num_of_allowed_heavy_atoms), "%)"
    print ""
    print "ALLOWED METHYL CARBONS ASSIGNED in TOCSY: "+str(num_of_assigned_methyl_carbons_TOCSY)+"/"+str(expected_number_of_methyl_carbons)+" (" + \
            str(100*float(num_of_assigned_methyl_carbons_TOCSY)/expected_number_of_methyl_carbons), " %)"
    print "CORRECT/WRONG METHYL CARBONS in TOCSY: "+str(num_of_correct_methyl_carbons_TOCSY)+"/"+str(num_of_wrong_methyl_carbons_TOCSY)+" ("+ \
            str(100*float(num_of_correct_methyl_carbons_TOCSY)/(expected_number_of_methyl_carbons))+"%/"+str(100*float(num_of_wrong_methyl_carbons_TOCSY)/(expected_number_of_methyl_carbons))+"%)"
    print ""
    print "ALLOWED METHYL CARBONS ASSIGNED in NOESY: "+str(num_of_assigned_methyl_carbons_NOESY)+"/"+str(expected_number_of_methyl_carbons)+" (" + \
            str(100*float(num_of_assigned_methyl_carbons_NOESY)/expected_number_of_methyl_carbons), "%)"
    print "CORRECT/WRONG METHYL CARBONS in NOESY: "+str(num_of_correct_methyl_carbons_NOESY)+"/"+str(num_of_wrong_methyl_carbons_NOESY)+" ("+ \
            str(100*float(num_of_correct_methyl_carbons_NOESY)/(expected_number_of_methyl_carbons))+"%/"+str(100*float(num_of_wrong_methyl_carbons_NOESY)/(expected_number_of_methyl_carbons))+"%)"
    print ""
    print "ALLOWED NON-METHYL CARBONS ASSIGNED in TOCSY: "+str(num_of_assigned_nonmethyl_carbons_TOCSY)+"/"+str(expected_number_of_nonmethyl_carbons)+" (" + \
            str(100*float(num_of_assigned_nonmethyl_carbons_TOCSY)/expected_number_of_nonmethyl_carbons), "%)"
    print "CORRECT/WRONG NON-METHYL CARBONS in TOCSY: "+str(num_of_correct_nonmethyl_carbons_TOCSY)+"/"+str(num_of_wrong_nonmethyl_carbons_TOCSY)+" ("+ \
            str(100*float(num_of_correct_nonmethyl_carbons_TOCSY)/(expected_number_of_nonmethyl_carbons))+"%/"+str(100*float(num_of_wrong_nonmethyl_carbons_TOCSY)/(expected_number_of_nonmethyl_carbons))+"%)"
    print ""
    print "ALLOWED NON-METHYL CARBONS ASSIGNED in NOESY: "+str(num_of_assigned_nonmethyl_carbons_NOESY)+"/"+str(expected_number_of_nonmethyl_carbons)+" (" + \
            str(100*float(num_of_assigned_nonmethyl_carbons_NOESY)/expected_number_of_nonmethyl_carbons), "%)"
    print "CORRECT/WRONG NON-METHYL CARBONS in NOESY: "+str(num_of_correct_nonmethyl_carbons_NOESY)+"/"+str(num_of_wrong_nonmethyl_carbons_NOESY)+" ("+ \
            str(100*float(num_of_correct_nonmethyl_carbons_NOESY)/(expected_number_of_nonmethyl_carbons))+"%/"+str(100*float(num_of_wrong_nonmethyl_carbons_NOESY)/(expected_number_of_nonmethyl_carbons))+"%)"
    print ""
    
    
    total_expected_number_of_carbons = expected_number_of_methyl_carbons + expected_number_of_nonmethyl_carbons
    total_num_of_assigned_carbons_TOCSY = num_of_assigned_methyl_carbons_TOCSY + num_of_assigned_nonmethyl_carbons_TOCSY
    total_num_of_correct_carbons_TOCSY = num_of_correct_methyl_carbons_TOCSY + num_of_correct_nonmethyl_carbons_TOCSY
    total_num_of_wrong_carbons_TOCSY = num_of_wrong_methyl_carbons_TOCSY + num_of_wrong_nonmethyl_carbons_TOCSY
    
    total_num_of_assigned_carbons_NOESY = num_of_assigned_methyl_carbons_NOESY + num_of_assigned_nonmethyl_carbons_NOESY
    total_num_of_correct_carbons_NOESY = num_of_correct_methyl_carbons_NOESY + num_of_correct_nonmethyl_carbons_NOESY
    total_num_of_wrong_carbons_NOESY = num_of_wrong_methyl_carbons_NOESY + num_of_wrong_nonmethyl_carbons_NOESY
    
    print "TOTAL ALLOWED CARBONS ASSIGNED in TOCSY: "+str(total_num_of_assigned_carbons_TOCSY)+"/"+str(total_expected_number_of_carbons)+" (" + \
            str(100*float(total_num_of_assigned_carbons_TOCSY)/total_expected_number_of_carbons), " %)"
    print "TOTAL CORRECT/WRONG CARBONS in TOCSY: "+str(total_num_of_correct_carbons_TOCSY)+"/"+str(total_num_of_wrong_carbons_TOCSY)+" ("+ \
            str(100*float(total_num_of_correct_carbons_TOCSY)/(total_expected_number_of_carbons))+"%/"+str(100*float(total_num_of_wrong_carbons_TOCSY)/(total_expected_number_of_carbons))+"%)"
    print ""
    print "TOTAL ALLOWED CARBONS ASSIGNED in NOESY: "+str(total_num_of_assigned_carbons_NOESY)+"/"+str(total_expected_number_of_carbons)+" (" + \
            str(100*float(total_num_of_assigned_carbons_NOESY)/total_expected_number_of_carbons), "%)"
    print "TOTAL CORRECT/WRONG CARBONS in NOESY: "+str(total_num_of_correct_carbons_NOESY)+"/"+str(total_num_of_wrong_carbons_NOESY)+" ("+ \
            str(100*float(total_num_of_correct_carbons_NOESY)/(total_expected_number_of_carbons))+"%/"+str(100*float(total_num_of_wrong_carbons_NOESY)/(total_expected_number_of_carbons))+"%)"
    print ""
    
    
    TOTAL_expected_number_of_carbons += total_expected_number_of_carbons
    TOTAL_num_of_assigned_carbons_TOCSY += total_num_of_assigned_carbons_TOCSY
    TOTAL_num_of_correct_carbons_TOCSY += total_num_of_correct_carbons_TOCSY
    TOTAL_num_of_wrong_carbons_TOCSY += total_num_of_wrong_carbons_TOCSY
    
    TOTAL_num_of_assigned_carbons_NOESY += total_num_of_assigned_carbons_NOESY
    TOTAL_num_of_correct_carbons_NOESY += total_num_of_correct_carbons_NOESY
    TOTAL_num_of_wrong_carbons_NOESY += total_num_of_wrong_carbons_NOESY
    
    print "BACKBONE NITROGENS ASSIGNED:", str(num_of_assigned_N)+"/"+str(total_num_of_N), "("+str(100*float(num_of_assigned_N)/total_num_of_N)+" %)"
    
    
    # # REDEFINITION OF of allowed_aa_atoms_dict, excluding CD2 from LEU and CG2 from VAL to avoid double counting !
    #         print "There is no ", aa_type, " in the protein! Moving to next aa type."
    #                 if resname != aa_type:
    #         y_pos = np.arange(len([a for a in allowed_aa_atoms_dict[aa_type] if a=="H" or a[0] != "H"]))
    #         plt.yticks(y_pos, [a for a in allowed_aa_atoms_dict[aa_type] if a=="H" or a[0] != "H"])
        

print "\n\n############################# TOTAL STATISTICS ################################\n\n"

print ""

print "TOTAL ALLOWED CARBONS ASSIGNED in TOCSY: "+str(TOTAL_num_of_assigned_carbons_TOCSY)+"/"+str(TOTAL_expected_number_of_carbons)+" (" + \
        str(100*float(TOTAL_num_of_assigned_carbons_TOCSY)/TOTAL_expected_number_of_carbons), " %)"
print "TOTAL CORRECT/WRONG CARBONS in TOCSY: "+str(TOTAL_num_of_correct_carbons_TOCSY)+"/"+str(TOTAL_num_of_wrong_carbons_TOCSY)+" ("+ \
        str(100*float(TOTAL_num_of_correct_carbons_TOCSY)/(TOTAL_expected_number_of_carbons))+"%/"+str(100*float(TOTAL_num_of_wrong_carbons_TOCSY)/(TOTAL_expected_number_of_carbons))+"%)"
print ""
print "TOTAL ALLOWED CARBONS ASSIGNED in NOESY: "+str(TOTAL_num_of_assigned_carbons_NOESY)+"/"+str(TOTAL_expected_number_of_carbons)+" (" + \
        str(100*float(TOTAL_num_of_assigned_carbons_NOESY)/TOTAL_expected_number_of_carbons), "%)"
print "TOTAL CORRECT/WRONG CARBONS in NOESY: "+str(TOTAL_num_of_correct_carbons_NOESY)+"/"+str(TOTAL_num_of_wrong_carbons_NOESY)+" ("+ \
        str(100*float(TOTAL_num_of_correct_carbons_NOESY)/(TOTAL_expected_number_of_carbons))+"%/"+str(100*float(TOTAL_num_of_wrong_carbons_NOESY)/(TOTAL_expected_number_of_carbons))+"%)"
print ""
    