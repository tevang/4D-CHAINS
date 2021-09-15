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

from argparse import ArgumentParser, RawDescriptionHelpFormatter
# Import 4D-CHAINS libraries
from lib.global_func import *


## Parse command line arguments
def cmdlineparse():
    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter, description="""
DESCRIPTION:
This is a script to create artificiall HSQC files from 4D-TOCSY and 4D-NOESY. The prerequisite is 
both spectra to have a single N,HN value pair for all the peaks of each AAIG.

TODO: include side chain N-H.

                            """,
                            epilog="""
EXAMPLE1:

    """)
    parser = ArgumentParser(description="command line arguments",
                            epilog="EXAMPLE: ")
    parser.add_argument("-tocsy", dest="TOCSY_fname", required=True,
                        help="4D TOCSY (HCTOCSYNH) file",
                        metavar="<4D TOCSY input file>")
    parser.add_argument("-noesy", dest="NOESY_FILE", required=True,
                        help="4D NOESY (HCNOENH) file",
                        metavar="<4D NOESY input file>")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s version 1.0')
    args=parser.parse_args()
    return args

args = cmdlineparse()

NHreson_set = set()
# Read 4D TOCSY and 4D NOESY lists
for fname in [args.TOCSY_fname, args.NOESY_FILE]:
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
            print("WARNING: Discarding from file", fname, " the line:", line)
            # print "DEBUG: copy_aaindices_from_root_spectrum_2 point1 word_list=", word_list

NHreson_list = list(NHreson_set)
NHreson_list.sort(key=itemgetter(0,1))
index = 1
froot = open("4DCHAINS_root.list", 'w')
for NH in NHreson_list:
    # TODO: include side chain N-H
    froot.write("\t" + "X"+str(index) + "NX-HX\t" + str(NH[0]) + "\t" + str(NH[1]) + "\n")
    index += 1
froot.close()