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



import sys, re, os, cPickle, traceback, shutil, bz2, math
from scoop import futures, shared
import numpy as np
from operator import itemgetter
from ordereddict import OrderedDict
from ete3 import Tree
import ftplib
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from scipy.stats.mstats import zscore
from scipy import stats, sqrt
import collections
import gc
from cluster import HierarchicalClustering
def tree(): # function to create multidimensional dictionaries
    return collections.defaultdict(tree)
CHAINS_BIN_DIR = os.path.dirname(os.path.realpath(__file__))
CHAINS_LIB_DIR = CHAINS_BIN_DIR[:-3] + "lib"
sys.path.append(CHAINS_BIN_DIR)
sys.path.append(CHAINS_LIB_DIR)
from global_func import *
from cs import guess_starting_resid, get_protein_sequence, read_NHmap


def cmdlineparse():
    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter, description="""
This script will read in the {N-H}-HSQC file in sparky list format. If there are labels inside that comply with the pattern '[A-Z][0-9]+N-H'
(e.g. A13N-H), the script will check if these labels are valid given the protein sequence. If yes, it will keep them and rename all other
peaks as 'X[0-9]+'. If notm then it will ignore all existing labels and relabel all peaks.
                            """,
                            epilog="EXAMPLE: ")
    parser.add_argument("-root", dest="ROOT_fname", required=True, help="2D N-H HSQC root spectrum", metavar="<2D N-H HSQC input file>")
    parser.add_argument("-absfile", dest="ABSOLUTE_MATCHES_FILE", required=False,
                        help="OPTIONAL: table with NH assignments",
                        metavar="<absolute matches file>")
    parser.add_argument("-rstart", dest="STARTING_RESID", required=False, type=int, default=None,
                        help="OPTIONAL: the resid of the first residue in the sequence. (default: %(default)s)",
                        metavar="<absolute matches file>")
    parser.add_argument("-tseq", dest="FASTA_fname", required=False, type=str,
                        help="teplate sequence file in fasta format", metavar="<template sequence file>")
    parser.add_argument("-o", dest="OUT_fname", required=False, type=str,
                        help="output file name (annotated {N-H}-HSQC).", metavar="<output file>")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s version 1.0')
    args=parser.parse_args()
    return args

args = cmdlineparse()
  

if args.ABSOLUTE_MATCHES_FILE:
    if not args.STARTING_RESID:
        print bcolors.WARNING + "WARNING: you have not specified the starting residue number with the option -rstart. I will try \
to guess it from the labels." + bcolors.ENDC  
    protein_sequence_list, absolute_RIGmatches_alignment_list  = read_NHmap(args.ABSOLUTE_MATCHES_FILE)
    
    TIG2residue_dict = {}
    for i, aa,TIG in zip(range(len(protein_sequence_list)), protein_sequence_list, absolute_RIGmatches_alignment_list):
        if not TIG in ['N/A', '-']:
            TIG2residue_dict[TIG] = aa + str(i + args.STARTING_RESID)
    
    with open(args.ROOT_fname, 'r') as f:
        contents = f.readlines()
        index = 1
        for i in range(len(contents)):
            line = contents[i]
            mo = re.match(".*\s+([0-9.]+)\s+([0-9.]+)\s*", line)    # the last 2 columns must be numbers
            if mo:
                Nreson = mo.group(1)
                Hreson = mo.group(2)
                mp = re.match("^\s*([A-Za-z-]+[0-9]+)+N-H\s+[0-9.]+\s+[0-9.]+\s*", line)
                if mp:
                    TIG = mp.group(1)
                    if TIG in TIG2residue_dict.keys():
                        residue = TIG2residue_dict[TIG]
                    else:
                        residue = TIG
                    contents[i] = "\t%s\t%s\t%s\n" % (str(residue)+"N-H", Nreson, Hreson)
                    index += 1
    
else:
    if not args.FASTA_fname:
        print bcolors.FAIL + "You must also provide a fasta sequence file using the option -tseq !" + bcolors.ENDC
        sys.exit(1)
    protein_sequence_list = get_protein_sequence(args.FASTA_fname)
    args.STARTING_RESID = guess_starting_resid(args.ROOT_fname, args.FASTA_fname, NHmap=None)
    with open(args.ROOT_fname, 'r') as f:
        contents = f.readlines()
        index = 1
        for i in range(len(contents)):
            line = contents[i]
            mo = re.match(".*\s+([0-9.]+)\s+([0-9.]+)\s*", line) # the last 2 columns must be numbers
            if mo:
                Nreson = mo.group(1)
                Hreson = mo.group(2)
                mp = re.match("^\s*([A-Z])([0-9]+)N-H\s+[0-9.]+\s+[0-9.]+\s*$", line)
                if mp: # if the labels matches the conditions, keep it
                    aa_type = mp.group(1)
                    resid = int(mp.group(2))
                    try:
                        if not args.STARTING_RESID: raise ValueError
                        if protein_sequence_list[resid - args.STARTING_RESID] == aa_type:
                            continue
                    except ValueError:
                        contents[i] = "\t%s\t%s\t%s\n" % ("X"+str(index)+"N-H", Nreson, Hreson)
                        index += 1
                else:
                    contents[i] = "\t%s\t%s\t%s\n" % ("X"+str(index)+"N-H", Nreson, Hreson)
                    index += 1

if args.OUT_fname:
    out_fname = args.OUT_fname
else:
    out_fname = args.ROOT_fname+"_annotated"
with open(out_fname, 'w') as f:
    for line in contents:
        f.write(line)