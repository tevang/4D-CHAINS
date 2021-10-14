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

## Parse command line arguments
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from lib.hsqc_spectrum import HSQC_spectrum


def cmdlineparse():
    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter, description=
"""
This script will read in the {N-H}-HSQC file in sparky list format. If there are labels inside that comply with the pattern '[A-Z][0-9]+N-H'
(e.g. A13N-H), the script will check if these labels are valid given the protein sequence. If yes, it will keep them and rename all other
peaks as 'X[0-9]+'. If not, then it will ignore all existing labels and relabel all peaks.
""",
                            epilog="""
EXAMPLE:

annotate_root.py -hsqc nEIt_HSQC.list -rstart 1 -fasta nEIt.fasta -o nEIt_HSQCnum.list

                            """
                                   )
    parser.add_argument("-hsqc", dest="HSQC_FILE", required=True,
                        help="2D N-H HSQC root spectrum",
                        metavar="<2D N-H HSQC input file>")
    parser.add_argument("-nhmap", dest="ABSOLUTE_MATCHES_FILE", required=False,
                        help="OPTIONAL: table with NH assignments to protein residues.",
                        metavar="<absolute matches file>")
    parser.add_argument("-rstart", dest="STARTING_RESID", required=False, type=int, default=None,
                        help="OPTIONAL: the resid of the first residue in the sequence. (default: %(default)s)",
                        metavar="<absolute matches file>")
    parser.add_argument("-fasta", dest="FASTA_fname", required=False, type=str,
                        help="teplate sequence file in fasta format", metavar="<template sequence file>")
    parser.add_argument("-o", dest="OUT_fname", required=False, type=str,
                        help="output file name (annotated {N-H}-HSQC).", metavar="<output file>")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s version 1.2')
    args=parser.parse_args()
    return args

args = cmdlineparse()
  
HSQC_spectrum.annotate_HSQC(args.HSQC_FILE,
                            args.FASTA_fname,
                            args.STARTING_RESID,
                            args.ABSOLUTE_MATCHES_FILE,
                            args.OUT_fname)

# if args.ABSOLUTE_MATCHES_FILE:  # TRY TO GUESS STARTING RESID FROM AAIG LABELS IN HSQC
#     if args.STARTING_RESID == None:
#         print(bcolors.WARNING + "WARNING: you have not specified the starting residue number with the option -rstart. I will try \
# to guess it from the labels." + bcolors.ENDC)
#     ## READ FILE WITH ABSOLUTE MATCHES
#     absolute_AAIGmatches_alignment_list, protein_sequence_list = \
#         Alignment.read_NHmap_file(args.ABSOLUTE_MATCHES_FILE, get_protein_alignment=True)
#
#     TIG2residue_dict = {}
#     for i, aa,TIG in zip(list(range(len(protein_sequence_list))), protein_sequence_list, absolute_AAIGmatches_alignment_list):
#         if not TIG in ['N/A', '-']:
#             TIG2residue_dict[TIG] = aa + str(i + args.STARTING_RESID)
#     # print "DEBUG: TIG2residue_dict=", TIG2residue_dict
#
#     with open(args.HSQC_FILE, 'r') as f:
#         contents = f.readlines()
#         index = 1
#         for i in range(len(contents)):
#             line = contents[i]
#             mo = re.match(".*\s+([0-9.]+)\s+([0-9.]+)\s*", line)    # the last 2 columns must be numbers
#             if mo:
#                 Nreson = mo.group(1)
#                 Hreson = mo.group(2)
#                 mp = re.match("^\s*([A-Za-z-]+[0-9]+)+N-H\s+[0-9.]+\s+[0-9.]+\s*", line)    # consider only N-H mapped AAIGs, ignore NX-HX, etc.
#                 if mp:
#                     TIG = mp.group(1)
#                     if TIG in list(TIG2residue_dict.keys()):
#                         residue = TIG2residue_dict[TIG]
#                     else:
#                         residue = TIG
#                     contents[i] = "\t%s\t%s\t%s\n" % (str(residue)+"N-H", Nreson, Hreson)
#                     index += 1
#
# else:
#     used_labels = []    # AAIG labels that have been assigned to HSQC peaks already
#     if not args.FASTA_fname:
#         print(bcolors.FAIL + "You must also provide a fasta sequence file using the option -fasta !" + bcolors.ENDC)
#         sys.exit(1)
#     protein_sequence_list = get_protein_sequence(args.FASTA_fname)
#     # print "DEBUG: protein_sequence_list=", protein_sequence_list
#     args.STARTING_RESID = HSQC_spectrum.guess_starting_resid(HSQC_FILE=args.HSQC_FILE, fasta=args.FASTA_fname, NHmap=None)
#     # If no NH-mapping table was given, consider labels meeting the conditions as real labels from the user and relabel the rest of the peaks
#     with open(args.HSQC_FILE, 'r') as f:
#         contents = f.readlines()
#         index = 1
#         for i in range(len(contents)):
#             line = contents[i]
#             mo = re.match(".*\s+([0-9.]+)\s+([0-9.]+)\s*", line) # the last 2 columns must be numbers
#             if mo:
#                 Nreson = mo.group(1)
#                 Hreson = mo.group(2)
#                 mp = re.match("^\s*([A-Z])([0-9]+)NX?-HX?\s+[0-9.]+\s+[0-9.]+\s*$", line)
#                 mq = re.match("^\s*([0-9]+)NX?-HX?\s+[0-9.]+\s+[0-9.]+\s*$", line)
#                 if mp: # if the label matches the conditions, keep it
#                     aa_type = mp.group(1)
#                     resid = int(mp.group(2))
#                     try:
#                         if not args.STARTING_RESID: raise ValueError
#                         # print "DEBUG: resid=", resid, "args.STARTING_RESID=", args.STARTING_RESID
#                         if protein_sequence_list[resid - args.STARTING_RESID] == aa_type:
#                             continue
#                     except ValueError:
#                         contents[i] = "\t%s\t%s\t%s\n" % ("X"+str(index)+"NX-HX", Nreson, Hreson)
#                         index += 1
#                     except IndexError:  # resid too high or too small, probably a sidechain, leave it as it is
#                         continue
#                 # TODO: add conditions about the side chains
#                 elif mq:    # if the label is of the form "[0-9]+NX?-HX?", then keep it to ease the user in the re-assignment
#                     user_index = mq.group(1)
#                     if "X" + user_index in used_labels:
#                         raise KeyError("Label "+user_index+" has been used already for another HSQC peak! Go back to the HSQC file"
#                                         "and change the label in the following line:\n"+line)
#                     contents[i] = "\t%s\t%s\t%s\n" % ("X" + user_index + "NX-HX", Nreson, Hreson)
#                 else:
#                     contents[i] = "\t%s\t%s\t%s\n" % ("X"+str(index)+"NX-HX", Nreson, Hreson)
#                     index += 1
#
# # Now write the modified contents to a new file
# if args.OUT_fname:
#     out_fname = args.OUT_fname
# else:
#     out_fname = args.HSQC_FILE+"_annotated"
# with open(out_fname, 'w') as f:
#     for line in contents:
#         f.write(line)