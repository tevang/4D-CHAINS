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



import sys, re, os
import numpy as np
from ordereddict import OrderedDict
from argparse import ArgumentParser
from operator import itemgetter
import collections, copy
def tree(): # function to create multidimensional dictionaries
    return collections.defaultdict(tree)

HOME_DIR = os.path.dirname(os.path.realpath(__file__))
code3 = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
code1 = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
aa1to3_dict = dict((c1,c3) for c1,c3 in zip(code1,code3))
aa3to1_dict = dict((c3,c1) for c3,c1 in zip(code3,code1))


aatype_carbon_nondegenerateHlist_multidict = tree() # names of non-degenerate methylene protons (in that case they are named HG2, HG3, etc.)
aatype_carbon_nondegenerateHlist_multidict['G']['CA'] = ['HA2', 'HA3']

aatype_carbon_nondegenerateHlist_multidict['D']['CB'] = ['HB2', 'HB3']
aatype_carbon_nondegenerateHlist_multidict['E']['CB'] = ['HB2', 'HB3']
aatype_carbon_nondegenerateHlist_multidict['H']['CB'] = ['HB2', 'HB3']
aatype_carbon_nondegenerateHlist_multidict['K']['CB'] = ['HB2', 'HB3']
aatype_carbon_nondegenerateHlist_multidict['L']['CB'] = ['HB2', 'HB3']
aatype_carbon_nondegenerateHlist_multidict['N']['CB'] = ['HB2', 'HB3']
aatype_carbon_nondegenerateHlist_multidict['F']['CB'] = ['HB2', 'HB3']
aatype_carbon_nondegenerateHlist_multidict['P']['CB'] = ['HB2', 'HB3']
aatype_carbon_nondegenerateHlist_multidict['Q']['CB'] = ['HB2', 'HB3']
aatype_carbon_nondegenerateHlist_multidict['R']['CB'] = ['HB2', 'HB3']
aatype_carbon_nondegenerateHlist_multidict['S']['CB'] = ['HB2', 'HB3']
aatype_carbon_nondegenerateHlist_multidict['C']['CB'] = ['HB2', 'HB3']
aatype_carbon_nondegenerateHlist_multidict['M']['CB'] = ['HB2', 'HB3']
aatype_carbon_nondegenerateHlist_multidict['W']['CB'] = ['HB2', 'HB3']
aatype_carbon_nondegenerateHlist_multidict['Y']['CB'] = ['HB2', 'HB3']

aatype_carbon_nondegenerateHlist_multidict['R']['CG'] = ['HG2', 'HG3']
aatype_carbon_nondegenerateHlist_multidict['Q']['CG'] = ['HG2', 'HG3']
aatype_carbon_nondegenerateHlist_multidict['E']['CG'] = ['HG2', 'HG3']
aatype_carbon_nondegenerateHlist_multidict['I']['CG1'] = ['HG12', 'HG13']
aatype_carbon_nondegenerateHlist_multidict['K']['CG'] = ['HG2', 'HG3']
aatype_carbon_nondegenerateHlist_multidict['M']['CG'] = ['HG2', 'HG3']
aatype_carbon_nondegenerateHlist_multidict['P']['CG'] = ['HG2', 'HG3']

aatype_carbon_nondegenerateHlist_multidict['R']['CD'] = ['HD2', 'HD3']
aatype_carbon_nondegenerateHlist_multidict['K']['CD'] = ['HD2', 'HD3']
aatype_carbon_nondegenerateHlist_multidict['P']['CD'] = ['HD2', 'HD3']

aatype_carbon_nondegenerateHlist_multidict['K']['CE'] = ['HE2', 'HE3']



def report_statistics(fpaths, norm_intensities_multidict):
    
    methyls_multidict = tree()   # resname->list of allowed carbons in iteratio1
    methyls_multidict["L"]["CD1"] = []
    methyls_multidict["L"]["CD2"] = []
    methyls_multidict["V"]["CG1"] = []
    methyls_multidict["V"]["CG2"] = []
    methyls_multidict["I"]["CD1"] = []
    methyls_multidict["I"]["CG2"] = []
    methyls_multidict["T"]["CG2"] = []
    methyls_multidict["A"]["CB"] = []
    methyls_multidict["M"]["CE"] = []
    
    methyl_ratios_dict = {}
    methyl_ratios_dict["L"] = []
    methyl_ratios_dict["V"] = []

    ratios_multidict = tree()
    KR_intensities = tree()
    KR_intensities["K"]["CD"] = []
    KR_intensities["K"]["CE"] = []
    KR_intensities["R"]["CD"] = []
    for aa in aatype_carbon_nondegenerateHlist_multidict.keys():
        for carbon in aatype_carbon_nondegenerateHlist_multidict[aa].keys():
            H2 = aatype_carbon_nondegenerateHlist_multidict[aa][carbon][0]
            H3 = aatype_carbon_nondegenerateHlist_multidict[aa][carbon][1]
            for fpath in fpaths:
                for AA in norm_intensities_multidict[fpath].keys():
                    if AA == aa:
                        for RESID in norm_intensities_multidict[fpath][AA].keys():
                            for C in norm_intensities_multidict[fpath][AA][RESID].keys():
                                if C == carbon:
                                    intensities = norm_intensities_multidict[fpath][AA][RESID][C].values()
                                    if len(intensities) == 2:
                                        intensities.sort(reverse=True)
                                        ratio = intensities[0]/float(intensities[1])
                                        try:
                                            ratios_multidict[aa][carbon].append(ratio)
                                        except (AttributeError, KeyError):
                                            ratios_multidict[aa][carbon] = [ratio]
                                        if aa in ["K", "R"] and C in ["CD", "CE"]:
                                            KR_intensities[aa][C].extend(intensities)
                                    elif len(intensities) == 1:
                                        if aa in ["K", "R"] and C in ["CD", "CE"]:
                                            KR_intensities[aa][C].extend(intensities)
    
    for aa in ratios_multidict.keys():
        for carbon in ratios_multidict[aa].keys():
            total = str(len(ratios_multidict[aa][carbon]))
            mu = np.mean(ratios_multidict[aa][carbon])
            stdev = np.std(ratios_multidict[aa][carbon])
            print "Average intensity ratio of aa "+aa+" and methylene carbon "+carbon+" is :", mu, "+-", stdev, " (#"+total+")"
    
    for aa in KR_intensities.keys():
        for C in KR_intensities[aa].keys():
            total = str(len(KR_intensities[aa][C]))
            mu = np.mean(KR_intensities[aa][C])
            stdev = np.std(KR_intensities[aa][C])
            print "Average intensity of aa "+aa+" and methylene carbon "+C+" is :", mu, "+-", stdev, " (#"+total+")"
    
    
    for aa in ["A", "I", "T", "M"]:
        for carbon in methyls_multidict[aa].keys():
            for fpath in fpaths:
                for AA in norm_intensities_multidict[fpath].keys():
                    if AA == aa:
                        for RESID in norm_intensities_multidict[fpath][AA].keys():
                            for C in norm_intensities_multidict[fpath][AA][RESID].keys():
                                if C == carbon:
                                    intensity = norm_intensities_multidict[fpath][AA][RESID][C].values()[0]   # there is only 1 proton and thus only 1 value
                                    methyls_multidict[aa][carbon].append(intensity)
                                    
    for aa in ["A", "I", "T", "M"]:
        for carbon in methyls_multidict[aa].keys():
            total = str(len(methyls_multidict[aa][carbon]))
            mu = np.mean(methyls_multidict[aa][carbon])
            stdev = np.std(methyls_multidict[aa][carbon])
            print "Average intensity of aa "+aa+" and methyl carbon "+carbon+" is :", mu, "+-", stdev, " (#"+total+")"
            
    equivalent_carbons_dict = { "CD1":"CD2", "CD2":"CD1", "CG1":"CG2", "CG2":"CG1"}
    for aa in ["L", "V"]:
        for carbon in methyls_multidict[aa].keys():
            for fpath in fpaths:
                for AA in norm_intensities_multidict[fpath].keys():
                    if AA == aa:
                        for RESID in norm_intensities_multidict[fpath][AA].keys():
                            for C in norm_intensities_multidict[fpath][AA][RESID].keys():
                                if C == carbon:
                                    intensity = norm_intensities_multidict[fpath][AA][RESID][C].values()[0]   # there is only 1 proton and thus only 1 value
                                    methyls_multidict[aa][carbon].append(intensity)
                                    if not equivalent_carbons_dict[C] in norm_intensities_multidict[fpath][AA][RESID].keys():
                                        methyls_multidict[aa][equivalent_carbons_dict[C]].append(intensity)
                                    else:   # if both equivalent carbons exist, save their intensity ratio
                                        intensity2 = norm_intensities_multidict[fpath][AA][RESID][equivalent_carbons_dict[C]].values()[0]
                                        if intensity >= intensity2:
                                            methyl_ratios_dict[aa].append(intensity/float(intensity2))
                                        else:
                                            methyl_ratios_dict[aa].append(intensity2/float(intensity))
                                    
    for aa, C_pair in zip(["L", "V"], ["CD1_CD2", "CG1_CG2"]):
        C1 = C_pair.split("_")[0]
        C2 = C_pair.split("_")[1]
        total = str(len(methyls_multidict[aa][C1] + methyls_multidict[aa][C2]))
        C1_array = np.array(methyls_multidict[aa][C1], dtype=float)
        C2_array = np.array(methyls_multidict[aa][C2], dtype=float)
        C1C2_array = np.vstack((C1_array, C2_array))
        C1C2_array.sort(axis=0)
        mean_ratio = np.mean(methyl_ratios_dict[aa])
        ratio_stdev = np.std(methyl_ratios_dict[aa])
        mu = np.mean(methyls_multidict[aa][C1] + methyls_multidict[aa][C2])
        stdev = np.std(methyls_multidict[aa][C1] + methyls_multidict[aa][C2])
        print "Average intensity of aa "+aa+" and equivalent methyl carbons "+ C1 +" and "+ C2 +" is :", mu, "+-", stdev, " (#"+total+")"
        print "Average intensity ratio of aa "+aa+" and equivalent methyl carbons "+ C1 +" and "+ C2 +" is :", mean_ratio, "+-", ratio_stdev



def find_missing_assignments(norm_intensities_multidict1, norm_intensities_multidict2):
    
    fpaths = set(norm_intensities_multidict1.keys() + norm_intensities_multidict2.keys())
    for fpath in fpaths:
        AAs = set(norm_intensities_multidict1[fpath].keys() + norm_intensities_multidict2[fpath].keys())
        for AA in AAs:
            RESIDs = set(norm_intensities_multidict1[fpath][AA].keys() + norm_intensities_multidict2[fpath][AA].keys())
            for RESID in RESIDs:
                CNAMEs = set(norm_intensities_multidict1[fpath][AA][RESID].keys() + norm_intensities_multidict2[fpath][AA][RESID].keys())
                for CNAME in CNAMEs:
                    HNAMEs = set(norm_intensities_multidict1[fpath][AA][RESID][CNAME].keys() + norm_intensities_multidict2[fpath][AA][RESID][CNAME].keys())
                    for HNAME in HNAMEs:
                        if not HNAME in norm_intensities_multidict1[fpath][AA][RESID][CNAME].keys() and HNAME in norm_intensities_multidict2[fpath][AA][RESID][CNAME].keys():
                            print "Missing intra-residue NOE peak assignments from", fpath.split("_")[0], AA+RESID + HNAME + "-" + CNAME 
                        elif not HNAME in norm_intensities_multidict2[fpath][AA][RESID][CNAME].keys() and HNAME in norm_intensities_multidict1[fpath][AA][RESID][CNAME].keys():
                            print "Missing inter-residue i+1<-i NOE peak assignments from", fpath.split("_")[0], AA+RESID + HNAME + "-" + CNAME 



os.chdir(".")
fpaths=os.listdir(".")
fpattern = re.compile("^.*.list.curated.normintensities$")
fpaths=filter(fpattern.search, fpaths)
norm_intensities_multidict1 = tree()
for fpath in fpaths:
    with open(fpath, 'r') as f:
        for line in f:
            components = line.split()
            mo = re.search('^([A-Z])([0-9]+)([HQM][ABGDE][1-3]{0,2})-(C[ABGDE][12]{0,1})-N-H$', components[0])
            if mo:
                aa = mo.group(1)
                resid = mo.group(2)
                Hname = mo.group(3)
                Cname = mo.group(4)
                norm_intensity = float(components[6])
                norm_intensities_multidict1[fpath][aa][resid][Cname][Hname] = norm_intensity


    
print "\n\t\t\tREPORT INTRA-RESIDUE NOE INTENSITY STATISTICS\n"
report_statistics(fpaths, norm_intensities_multidict1)


os.chdir(".")
fpaths=os.listdir(".")
fpattern = re.compile("^.*.list.curated.normintensities$")
fpaths=filter(fpattern.search, fpaths)
norm_intensities_multidict2 = tree()
for fpath in fpaths:
    with open(fpath, 'r') as f:
        for line in f:
            components = line.split()
            mo = re.search('^([A-Z])([0-9]+)([HQM][ABGDE][1-3]{0,2})-(C[ABGDE][12]{0,1})-([A-Z])([0-9]+)N-H$', components[0])
            if mo:
                i_aa = mo.group(1)
                i_resid = mo.group(2)
                i_Hname = mo.group(3)
                i_Cname = mo.group(4)
                iplus1_aa = mo.group(5)
                iplus1_resid = mo.group(6)
                if not int(iplus1_resid) == float(i_resid) + 1:
                    continue
                norm_intensity = float(components[6])
                norm_intensities_multidict2[fpath][i_aa][i_resid][i_Cname][i_Hname] = norm_intensity


print "\n\t\t\tREPORT INTER-RESIDUE i<-i+1 NOE INTENSITY STATISTICS\n"
report_statistics(fpaths, norm_intensities_multidict2)

os.chdir(".")
fpaths=os.listdir(".")
fpattern = re.compile("^.*.list.curated.normintensities$")
fpaths=filter(fpattern.search, fpaths)
norm_intensities_multidict3 = tree()
for fpath in fpaths:
    with open(fpath, 'r') as f:
        for line in f:
            components = line.split()
            mo = re.search('^([A-Z])([0-9]+)([HQM][ABGDE][1-3]{0,2})-(C[ABGDE][12]{0,1})-([A-Z])([0-9]+)N-H$', components[0])
            if mo:
                i_aa = mo.group(1)
                i_resid = mo.group(2)
                i_Hname = mo.group(3)
                i_Cname = mo.group(4)
                iplus1_aa = mo.group(5)
                iplus1_resid = mo.group(6)
                if not int(iplus1_resid) == float(i_resid) - 1:
                    continue
                norm_intensity = float(components[6])
                norm_intensities_multidict3[fpath][i_aa][i_resid][i_Cname][i_Hname] = norm_intensity


print "\n\t\t\tREPORT INTER-RESIDUE i<-i-1 NOE INTENSITY STATISTICS\n"
report_statistics(fpaths, norm_intensities_multidict3)

# find_missing_assignments(norm_intensities_multidict1, norm_intensities_multidict2)