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
import os
from thirdparty.b2bTools.singleSeq.DynaMine.Predictor import *
from lib.fasta import *
from lib.global_func import *

## Parse command line arguments
def cmdlineparse():
    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter,
                            description="""
This script runs DynaMine on a protein sequence.
""",
        epilog="""
EXAMPLE:


""")
    parser.add_argument("-fasta", dest="FASTA", required=False, type=str, default=None,
                        help="fasta file",
                        metavar="<fasta file>")
    parser.add_argument("-fastapattern", dest="FASTA_PATTERN", required=False, default=".*\.fasta",
                        help="Instead of -fasta, use this argument to give the file name pattern of the .fasta files.")
    # parser.add_argument("-outfile", dest="OUT_FILE", required=False, type=str, default=None,
    #                     help="output file name. It have the format:"
    #                          "protein_ID    cluster_ID",
    #                     metavar="<output file>")

    args=parser.parse_args()
    return args

args = cmdlineparse()

if args.FASTA:
    fname_list = [args.FASTA]
else:
    fname_list = list_files(".", pattern=args.FASTA_PATTERN)

seqs = []
for fname in fname_list:
    fasta = FASTA(fname, full_seq=True)
    seqs.append( (fasta.protname, fasta.protseq_str) )

dm = DynaMine()
start = 0
for end in range(0, len(seqs), 10):
    dm.predictSeqs(seqs[start:end])
    print("\nPreds=", dm.allPredictions)
    start = end
