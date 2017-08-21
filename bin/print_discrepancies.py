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




import sys
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from collections import OrderedDict
from operator import itemgetter
class tree(OrderedDict):
    def __missing__(self, key):
        self[key] = type(self)()
        return self[key]


code3 = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
code1 = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
aa1to3_dict = dict((c1,c3) for c1,c3 in zip(code1,code3))
aa3to1_dict = dict((c3,c1) for c3,c1 in zip(code3,code1))

def cmdlineparse():
    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter, description="""
DESCRIPTION:
   Prints discrepancies between two assignment files and write a file for ploting the correlation and overlays of assignment resonances. 
                            """,
    epilog="""
EXAMPLES:
    
    collect_MMGBSA_results.py -cf compound_IC50.filtered.txt
    """)
    
    parser.add_argument("-n1", dest="Nassignments1", required=False, type=str, default=None,
                        help="Nitrogen assignments file 1.",
                        metavar="<>")
    parser.add_argument("-n2", dest="Nassignments2", required=False, type=str, default=None,
                        help="Nitrogen assignments file 2.",
                        metavar="<>")
    parser.add_argument("-c1", dest="Cassignments1", required=False, type=str, default=None,
                        help="Carbon assignments file 1.",
                        metavar="<>")
    parser.add_argument("-c2", dest="Cassignments2", required=False, type=str, default=None,
                        help="Carbon assignments file 2.",
                        metavar="<>")
    parser.add_argument("-h1", dest="Hassignments1", required=False, type=str, default=None,
                        help="Proton assignments file 1.",
                        metavar="<>")
    parser.add_argument("-h2", dest="Hassignments2", required=False, type=str, default=None,
                        help="Proton assignments file 2.",
                        metavar="<>")
    parser.add_argument("-ntol", dest="Ntol", required=False, type=float, default=0.4,
                        help="Nitrogen tolerance to find discrepances between the two files.",
                        metavar="<>")
    parser.add_argument("-ctol", dest="Ctol", required=False, type=float, default=0.4,
                        help="Carbon tolerance to find discrepances between the two files.",
                        metavar="<>")
    parser.add_argument("-htol", dest="Htol", required=False, type=float, default=0.04,
                        help="Proton tolerance to find discrepances between the two files.",
                        metavar="<>")
    
    args=parser.parse_args()
    return args


args = cmdlineparse()


def select_matches(duplets1, duplets2):
    
    print "DEBUG: duplets1=", duplets1
    print "DEBUG: duplets2=", duplets2
    
    matches_dict = {}
    if len(duplets1) == 1:
        reson1 = duplets1[0][1]
        diff_list = []
        for d2 in duplets2:
            reson2 = d2[1]
            diff_list.append(abs(reson1-reson2))
        minindex = diff_list.index(min(diff_list))
        matching_peak = duplets2[minindex]
        matches_dict[duplets1[0]] = matching_peak
    elif len(duplets2) == 1:
        reson2 = duplets2[0][1]
        diff_list = []
        for d1 in duplets1:
            reson1 = d1[1]
            diff_list.append(abs(reson1-reson2))
        minindex = diff_list.index(min(diff_list))
        matching_peak = duplets1[minindex]
        matches_dict[matching_peak] = duplets2[0]
    if len(duplets1) == 2 and len(duplets2) == 2:
        reson11 = duplets1[0][1]
        reson12 = duplets1[1][1]
        reson21 = duplets2[0][1]
        reson22 = duplets2[1][1]
        diff1121 = abs(reson11- reson21)
        diff1122 = abs(reson11- reson22)    
        if diff1121 < diff1122:
            matches_dict[duplets1[0]] = duplets2[0]
            matches_dict[duplets1[1]] = duplets2[1]
        else:
            matches_dict[duplets1[1]] = duplets2[0]
            matches_dict[duplets1[0]] = duplets2[1]
        
    return matches_dict


if args.Nassignments1 and args.Nassignments2: 

    Ndata1 = tree()
    Ndata2 = tree()
    resid2aa_dict = {}
    with open(args.Nassignments1, 'r') as f:
        for line1 in f:
            words1 = line1.split()
            if len(words1) == 0: continue
            if len(words1[0]) == 1: words1[0] = aa1to3_dict[words1[0]]
            resid2aa_dict[int(words1[1])] = words1[0]
            Ndata1[int(words1[1])] = float(words1[3])
    with open(args.Nassignments2, 'r') as f:
        for line2 in f:
            words2 = line2.split()
            if len(words2) == 0: continue
            if len(words2[0]) == 1: words2[0] = aa1to3_dict[words2[0]]
            resid2aa_dict[int(words2[1])] = words2[0]
            Ndata2[int(words2[1])] = float(words2[3])                        
    all_resids = list(set(Ndata1.keys() + Ndata2.keys()))
    all_resids = sorted(all_resids)
    
    fout = open("N_for_plotting.dat", 'w')
    for resid in all_resids:
        if not resid in Ndata1.keys() or not resid in Ndata2.keys():
            print "Skipping ", resid2aa_dict[resid], "resid", resid
            continue
        diff = abs(Ndata1[resid]-Ndata2[resid])
        if diff > args.Ntol:
            print "DISCREPANCY: Matching ", resid2aa_dict[resid], "resid", resid, "N :", Ndata1[resid], Ndata2[resid], "diff=", diff
        else:
            print "Matching ", resid2aa_dict[resid], "resid", resid, "N :", Ndata1[resid], Ndata2[resid], "diff=", diff
        fout.write(str(resid) + " N " + str(Ndata1[resid]) + " N " + str(Ndata2[resid]) + "\n")
    fout.close()


if args.Cassignments1 and args.Cassignments2: 

    Cdata1 = tree()
    Cdata2 = tree()
    resid2aa_dict = {}
    with open(args.Cassignments1, 'r') as f:
        for line1 in f:
            words1 = line1.split()
            if len(words1) == 0: continue
            if len(words1[0]) == 1: words1[0] = aa1to3_dict[words1[0]]
            resid2aa_dict[int(words1[1])] = words1[0]
            Cdata1[int(words1[1])][words1[2]] = float(words1[3])
    with open(args.Cassignments2, 'r') as f:
        for line2 in f:
            words2 = line2.split()
            if len(words2) == 0: continue
            if len(words2[0]) == 1: words2[0] = aa1to3_dict[words2[0]]
            resid2aa_dict[int(words2[1])] = words2[0]
            Cdata2[int(words2[1])][words2[2]] = float(words2[3])                        
    all_resids = list(set(Cdata1.keys() + Cdata2.keys()))
    all_resids = sorted(all_resids)
    
    fout = open("C_for_plotting.dat", 'w')
    for resid in all_resids:
        if not resid in Cdata1.keys() or not resid in Cdata2.keys():
            print "Skipping ", resid2aa_dict[resid], "resid", resid
            continue
        for C in set(Cdata1[resid].keys() + Cdata2[resid].keys()):
            if not C in Cdata1[resid].keys() or not C in Cdata2[resid].keys():
                print "Skipping ", resid2aa_dict[resid], "resid", resid, "carbon", C
                continue
            if resid2aa_dict[resid] in ["TYR", "PHE"] and C in ['CE1', 'CE2', 'CD1', 'CD2']:
                print "Skipping aromatic carbon:", resid2aa_dict[resid], "resid", resid, "carbon", C
                continue
            if resid2aa_dict[resid] in ["LEU", "VAL"] and C in ['CG1', 'CG2', 'CD1', 'CD2']:    # do the assignment of equivalent carbons later
                continue
            diff = abs(Cdata1[resid][C]-Cdata2[resid][C])
            if diff > args.Ctol:
                print "DISCREPANCY: Matching ", resid2aa_dict[resid], "resid", resid, C, ":", Cdata1[resid][C], Cdata2[resid][C], "diff=", diff
            else:
                print "Matching ", resid2aa_dict[resid], "resid", resid, C, ":", Cdata1[resid][C], Cdata2[resid][C], "diff=", diff
            fout.write(str(resid) + " " + C + " " + str(Cdata1[resid][C]) + " " + C + " " + str(Cdata2[resid][C]) + "\n")
        if resid2aa_dict[resid] in ["LEU", "VAL"]:
            duplets1 = [(C, Cdata1[resid][C]) for C in ['CG1', 'CG2', 'CD1', 'CD2'] if C in Cdata1[resid].keys()]
            duplets2 = [(C, Cdata2[resid][C]) for C in ['CG1', 'CG2', 'CD1', 'CD2'] if C in Cdata2[resid].keys()]
            if len(duplets1) == 0 or len(duplets2) == 0:
                continue
            matches_dict = select_matches(duplets1, duplets2)
            for duplet1 in matches_dict.keys():
                duplet2 = matches_dict[duplet1]
                diff = abs(duplet1[1]-duplet2[1])
                if diff > args.Ctol:
                    print "DISCREPANCY: Matching ", resid2aa_dict[resid], "resid", resid, ":", duplet1, duplet2, "diff=", diff
                else:
                    print "Matching ", resid2aa_dict[resid], "resid", resid, ":", duplet1, duplet2, "diff=", diff
                fout.write(str(resid) + " " + duplet1[0] + " " + str(duplet1[1]) + " " + duplet2[0] + " " + str(duplet2[1]) + "\n")
    fout.close()
        

if args.Hassignments1 and args.Hassignments2: 
    Hdata1 = tree()
    Hdata2 = tree()
    resid2aa_dict = {}
    with open(args.Hassignments1, 'r') as f:
        for line1 in f:
            words1 = line1.split()
            if len(words1) == 0: continue
            if len(words1[0]) == 1: words1[0] = aa1to3_dict[words1[0]]
            resid2aa_dict[int(words1[1])] = words1[0]
            if words1[2] == 'H':
                Hdata1[int(words1[1])]['H'] = [ (words1[2], float(words1[3])) ]
                continue
            print words1
            if not words1[1] in Hdata1.keys() and not words1[2][1] in Hdata1[int(words1[1])].keys():
                Hdata1[int(words1[1])][words1[2][1]] = [ (words1[2], float(words1[3])) ] 
            else:
                Hdata1[int(words1[1])][words1[2][1]].append( (words1[2], float(words1[3])) )
    with open(args.Hassignments2, 'r') as f:
        for line2 in f:
            words2 = line2.split()
            if len(words2) == 0: continue
            if len(words2[0]) == 1: words2[0] = aa1to3_dict[words2[0]]
            resid2aa_dict[int(words2[1])] = words2[0]
            if words2[2] == 'H':
                Hdata2[int(words2[1])]['H'] = [ (words2[2], float(words2[3])) ]
                continue
            if not words2[1] in Hdata2.keys() and not words2[2][1] in Hdata2[int(words2[1])].keys():
                Hdata2[int(words2[1])][words2[2][1]] = [ (words2[2], float(words2[3])) ]
            else:
                Hdata2[int(words2[1])][words2[2][1]].append( (words2[2], float(words2[3])) )
                        
    all_resids = list(set(Hdata1.keys() + Hdata2.keys()))
    all_resids = sorted(all_resids)
    for k1 in Hdata1:
        for k2 in Hdata1[k1].keys():
            Hdata1[k1][k2].sort(key=itemgetter(1))
    for k1 in Hdata2:
        for k2 in Hdata2[k1].keys():
            Hdata2[k1][k2].sort(key=itemgetter(1))
    
    fout = open("H_for_plotting.dat", 'w')
    for resid in all_resids:
        if not resid in Hdata1.keys() or not resid in Hdata2.keys():
            print "Skipping ", resid2aa_dict[resid], "resid", resid
            continue
        for H in set(Hdata1[resid].keys() + Hdata2[resid].keys()):
            if not H in Hdata1[resid].keys() or not H in Hdata2[resid].keys():
                print "Skipping ", resid2aa_dict[resid], "resid", resid, "proton", H
                continue
            if resid2aa_dict[resid] in ["TYR", "PHE"] and H in ['E', 'D']:
                print "Skipping aromatic proton:", resid2aa_dict[resid], "resid", resid, "proton", H
                continue
            matches_dict = select_matches(Hdata1[resid][H], Hdata2[resid][H])
            for duplet1 in matches_dict.keys():
                duplet2 = matches_dict[duplet1]
                diff = abs(duplet1[1]-duplet2[1])
                if diff > args.Htol:
                    print "DISCREPANCY: Matching ", resid2aa_dict[resid], "resid", resid, ":", duplet1, duplet2, "diff=", diff
                else:
                    print "Matching ", resid2aa_dict[resid], "resid", resid, ":", duplet1, duplet2, "diff=", diff
                fout.write(str(resid) + " " + duplet1[0] + " " + str(duplet1[1]) + " " + duplet2[0] + " " + str(duplet2[1]) + "\n")
    fout.close()
            