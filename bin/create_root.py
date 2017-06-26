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
from argparse import ArgumentParser
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


def cmdlineparse():
    parser = ArgumentParser(description="command line arguments",
                            epilog="EXAMPLE: ")
    parser.add_argument("-tocsy", dest="TOCSY_fname", required=True, help="4D TOCSY (HCTOCSYNH) file", metavar="<4D TOCSY input file>")
    parser.add_argument("-noesy", dest="NOESY_fname", required=True, help="4D NOESY (HCNOENH) file", metavar="<4D NOESY input file>")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s version 1.0')
    args=parser.parse_args()
    return args

args = cmdlineparse()

NHreson_set = set()
for fname in [args.TOCSY_fname, args.NOESY_fname]:
    with open(fname, 'r') as f:
        tmp_query_contents=f.readlines()    # contents of original query_fname (4D TOCSY or 4D NOESY) in 5 column format (name H C N HN)
    for line in tmp_query_contents: # CLEANING THE SPECTRUM FROM IRRELEVANT LINES
        word_list = line.split()
        try:
            float(word_list[1]); float(word_list[2]); float(word_list[3]); float(word_list[4]); # checking if the line contains numbers (resonances)
            if len(word_list) == 6:    # if there is intensity column check if it is a number
                float(word_list[5])
                word_list[5] = str(abs(float(word_list[5])))    # convert all intensities to positive numbers
            NHreson_set.add( (float(word_list[3]), float(word_list[4])) )
        except (IndexError, ValueError):
            print "WARNING: Discarding from file", fname, " the line:", line

NHreson_list = list(NHreson_set)
NHreson_list.sort(key=itemgetter(0,1))
index = 1
froot = open("4DCHAINS_root.list", 'w')
for NH in NHreson_list:
    froot.write("\t" + "X"+str(index) + "N-H\t" + str(NH[0]) + "\t" + str(NH[1]) + "\n")
    index += 1
froot.close()