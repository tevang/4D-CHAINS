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

import sys, re, os, traceback, shutil, bz2, math, multiprocessing
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
from cs import guess_starting_resid



def cmdlineparse():
    parser = ArgumentParser(description="""
The script will create the following files in the work directory:



tocsyNEWnum.list
noesyNEWnum.list
4DCHAINS_HSQC.list
4DCHAINS_NHmap
4DCHAINS_assigned_NH_HSQC.sparky
4DTOCSY_assignedall.xeasy
4DTOCSY_assignedall.sparky
4DNOESY_assignedall.sparky
4DNOESY_assignedall.proofread.xeasy


only4DNOESY_assignedall.sparky
only4DNOESY_assignedall.xeasy


                            """,
                            formatter_class=RawDescriptionHelpFormatter,
                            epilog="""
EXAMPLES:

4Dchains.py -protocol 4DCHAINS_protocol.txt

""" )
    parser.add_argument("-protocol", dest="protocol_file", required=False, default=None,
                        help="protocol file", metavar="<protocol file>")
    parser.add_argument("-workdir", dest="WORKDIR", required=False, type=str, default=os.getcwd(),
                        help="the work directory (default: the current directory: %(default)s )",
                        metavar="<work directory>")
    parser.add_argument("-cpus", dest="CPUs", required=False, type=int, default=multiprocessing.cpu_count(),
                        help="number of processes for multi-threading calculations",
                        metavar="<number of processes>")
    parser.add_argument("-exhaustiveness", dest="EXHAUSTIVENESS", required=False, type=int, default=1,
                        help="""the exhaustuveness level of the NH resonance mapping tp the protein sequence. 1: 4 cycles are conducted, each one
with lower cutoffs than the previous; 2: 10 cycles are conducted; 3: 16 cycles are conducted. (default: %(default)s)""",
                        metavar="<exhaustiveness level>")
    parser.add_argument("-writeprotocol", dest="WRITE_PROTOCOL", required=False, action='store_true', default=False,
                        help="Write the default protocol file given the selected exhaustiveness level and exits. Use this to generate a protocol\
                             file that you will edit and use it later for assignment calculations.")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s version 1.0')
    args=parser.parse_args()
    return args


args = cmdlineparse()


if args.EXHAUSTIVENESS == 1:
    directives = {
        "mcutoff": ['1.0','1.0','0.8','0.6','0.6'],
        "zmcutoff": ['0.0','0.0','-1.0','0.0','-1.0'],
        "zacutoff": ['-0.5','-1.0','-1.0','-0.5','-1.0'],
        "first_length": [6]*5,
        "last_length":  [4]*5,
        'poolconfile': [None]*5,
        'allconfile': [None]*5,
        'poolaafile': [None]*5,
        'allaafile': [None]*5,
        "rst_from_prev_cycle": [False] + [True]*4
    }
elif args.EXHAUSTIVENESS == 2:
    directives = {
        "mcutoff": ['1.0','1.0','0.8','0.6','0.6','0.5','0.5','0.5','1.0','1.0'],
        "zmcutoff": ['0.0','0.0','-1.0','0.0','-1.0','0.0','0.0','-1.0','0.0','0.0'],
        "zacutoff": ['-0.5','-1.0','-1.0','-0.5','-1.0','0.0','-0.5','-1.0','-0.5','-1.0'],
        "first_length": [6]*8 + [3]*2,
        "last_length":  [4]*8 + [3]*2,
        'poolconfile': [None]*10,
        'allconfile': [None]*10,
        'poolaafile': [None]*10,
        'allaafile': [None]*10,
        "rst_from_prev_cycle": [False] + [True]*9
    }
elif args.EXHAUSTIVENESS == 3:
    directives = {
        "mcutoff": ['1.0','1.0','1.0','0.8','0.8','0.6','0.6','0.6','0.6','0.5','0.5','0.5','0.5','1.0','1.0','1.0'],
        "zmcutoff": ['0']*4 + ['-1.0','0.0','0.0','-1.0','-1.0','0.0','0.0','-1.0','-1.0','0.0','0.0','0.0'],
        "zacutoff": ['-0.5','-0.5','-1.0','-1.0','-1.0','0.0','-1.0','0.0','-1.0','0.0','-1.0','-1.0','-100.0','0.0','-0.5','-1.0'],
        "first_length": [6]*13 + [3]*3,
        "last_length":  [4]*13 + [3]*3,
        'poolconfile': [None]*16,
        'allconfile': [None]*16,
        'poolaafile': [None]*16,
        'allaafile': [None]*16,
        "rst_from_prev_cycle": [False] + [True]*15
    }

directives['tolH'] = '0.04'
directives['tolC'] = '0.4'
directives['doNHmapping'] = True
directives['doassign4DTOCSY'] = True
directives['doassign4DNOESY'] = True
directives['doassignonly4DNOESY'] = False
directives['onlyCACB'] = False
directives['patch'] = False

directives['rstart'] = None
directives['4DTOCSY'] = ""
directives['4DNOESY'] = ""
directives['4DTOCSY_assignedNH'] = ""
directives['4DTOCSY_assignedall'] = ""
directives['4DNOESY_assignedNH'] = ""
directives['user_4DTOCSY_assignedall'] = ""
directives['user_4DNOESY_assignedall'] = ""
directives['fasta'] = ""
directives['HSQC'] = ""
directives['rst_file'] = ""
directives['con_file'] = ""
directives['aa_file'] = ""
directives['NHmap'] = ""


external_strlist_directives = [
                       'mcutoff',
                       'zacutoff',
                       'zmcutoff'
                       ]
external_intlist_directives = [
                       'first_length',
                       'last_length'
                       ]
external_booleanlist_directives = [
                        'rst_from_prev_cycle'
                        ]
external_value_directives = [
                        'rstart',
                        'tolH',
                        'tolC',
                        ]

external_booleanvalue_directives = [
                        'doNHmapping',
                        'doassign4DTOCSY',
                        'doassign4DNOESY',
                        'doassignonly4DNOESY',
                        ]

external_file_directives = [
    'fasta',
    'HSQC',
    '4DTOCSY',
    '4DNOESY',
    'NHmap',
    '4DTOCSY_assignedNH',
    '4DTOCSY_assignedall',
    '4DNOESY_assignedNH',
    'rst_file',
    'con_file',
    'aa_file',
    'user_4DTOCSY_assignedall',
    'user_4DNOESY_assignedall',
                        ]

compulsory_directives = [
                        'mcutoff',
                        'zacutoff',
                        'zmcutoff',
                        'first_length',
                        'last_length',
                        'doNHmapping',
                        'doassign4DTOCSY',
                        'doassign4DNOESY',
                        '4DTOCSY',
                        '4DNOESY',
                        'fasta',
                        'HSQC'
                        ]
optional_directives = set(directives.keys()).difference(compulsory_directives)


all_directives = external_intlist_directives + external_strlist_directives + external_booleanlist_directives + \
external_value_directives + external_booleanvalue_directives + external_file_directives


def str_to_bool(s):
    if s == 'True':
         return True
    elif s == 'False':
         return False
    else:
         raise ValueError(s) # evil ValueError that doesn't tell you what the wrong value was

if args.protocol_file and not args.WRITE_PROTOCOL:  # LOAD PROTOCOL FILE
    with open(args.protocol_file, 'r') as f:
        for line in f:
            words = line.split('=')
            if len(words) == 0 or not words[0] in all_directives:
                continue
            directive = words[0]
            if directive in external_strlist_directives:
                values = [v.strip() for v in words[1].split(',')]
                directives[directive] = values
            elif directive in external_intlist_directives:
                values = [int(v.strip()) for v in words[1].split(',')]
                directives[directive] = values
            elif directive in external_booleanvalue_directives:
                value = words[1].strip()
                directives[directive] = str_to_bool(value)
            elif directive in external_value_directives:
                value = words[1].strip()
                directives[directive] = value
            elif directive in external_booleanlist_directives:
                values = [str_to_bool(v.strip()) for v in words[1].split(',')]
                directives[directive] = values
            elif directive in external_file_directives:
                value = words[1].strip()
                if len(value) > 0:
                   value = os.path.abspath(value) 
                directives[directive] = value
    
    if len(set([len(directives[k]) for k in external_strlist_directives + external_intlist_directives + external_booleanlist_directives])) != 1:
        print bcolors.FAIL + "ERROR: all directives in file "+ args.protocol_file+" must have the same number of values separated by ',' !!!" + bcolors.ENDC
        sys.exit(1)

else:   # otherwise write an example protocol file with the current parameter values
    if len(set([len(directives[k]) for k in external_strlist_directives + external_intlist_directives + external_booleanlist_directives])) != 1:
        print bcolors.FAIL + "ERROR: all directives in file "+ args.protocol_file+" must have the same number of values separated by ',' !!!" + bcolors.ENDC
        sys.exit(1)
    with open("4DCHAINS_protocol.txt", 'w') as f:
        f.write("# COMPULSORY DIRECTIVES\n\n")
        f.write("# input files\n")
        for k in external_file_directives:
            if k in compulsory_directives:
                v = directives[k]
                if v == None: v = ""
                f.write(str(k) + "=" + str(v) + "\n")
        f.write("\n# parameters\n")
        for k in external_booleanvalue_directives + external_value_directives:
            if k in compulsory_directives:
                v = directives[k]
                if v == None: v = ""
                f.write(str(k) + "=" + str(v) + "\n")
        for k in external_intlist_directives + external_strlist_directives + external_booleanlist_directives:
            if k in compulsory_directives:
                v = directives[k]
                f.write(str(k) + "=" + ",".join([str(w) for w in v]) + "\n")
        
        f.write("\n# OPTIONAL DIRECTIVES\n\n")
        f.write("# input files\n")
        for k in external_file_directives:
            if k in optional_directives:
                v = directives[k]
                if v == None: v = ""
                f.write(str(k) + "=" + str(v) + "\n")
        f.write("\n# parameters\n")
        for k in external_value_directives + external_booleanvalue_directives:
            if k in optional_directives:
                v = directives[k]
                if v == None: v = ""
                f.write(str(k) + "=" + str(v) + "\n")
        for k in external_intlist_directives + external_strlist_directives + external_booleanlist_directives:
            if k in optional_directives:
                v = directives[k]
                f.write(str(k) + "=" + ",".join([str(w) for w in v]) + "\n")

if args.WRITE_PROTOCOL:
    print bcolors.BOLDGREEN + "Just wrote a sample protocol file ('4DCHAINS_protocol.txt') and exiting." + bcolors.ENDBOLD
    sys.exit(0)


def generate_input_data():
    
    global directives
    
    run_commandline("python2.7 -m scoop -n 4 "+CHAINS_BIN_DIR+"/4D_assignment_tree-based_parallel.py \
    -root "+directives['HSQC']+" \
    -tocsy "+directives['4DTOCSY']+" \
    -noesy "+directives['4DNOESY']+" \
    -tolH "+directives['tolH']+" -tolC "+directives['tolC']+" \
    -mcutoff 1.0 \
    -zmcutoff 1.0 \
    -acutoff 1.0 \
    -zacutoff 1.0 \
    -maxlen 3 -minlen 3 \
    -tseq "+directives['fasta']+" \
    -probprod -2dhist -log -delpred -cgrpprob 2 -skipchains") #>& round0.log
    
    
    shutil.copyfile("connectivities_all", "connectivities_all.pool")
    shutil.copyfile("amino_acid_type_prediction_probabilities", "amino_acid_type_prediction_probabilities.pool")
    directives['4DTOCSY_assignedNH'] = directives['4DTOCSY'].replace(".list", "num.list")
    directives['4DNOESY_assignedNH'] = directives['4DNOESY'].replace(".list", "num.list")

def multiround_doNHmapping(cycle=1, tolH='0.04', tolC='0.4', mcutoff='0.8', zmcutoff='0.0', zacutoff='-1.0',
                         poolconfile="", allconfile="", poolaafile="", allaafile="",
                              first_length=6, last_length=4, rst_file=None):
    
    global directives,args
    
    args.CPUs = str(args.CPUs)
    
    cur_round = '1'
    if rst_file:
        restraints="-rstfile " + rst_file
    else:
        restraints=""    
    pept_length = str(first_length)
    
    run_commandline("python2.7 -m scoop -n "+args.CPUs+" "+CHAINS_BIN_DIR+"/4D_assignment_tree-based_parallel.py \
    -root "+directives['HSQC']+" \
    -tocsy "+directives['4DTOCSY']+" \
    -noesy "+directives['4DNOESY']+" \
    -tolH "+tolH+" -tolC "+tolC+" \
    -mcutoff "+mcutoff+" \
    -zmcutoff "+zmcutoff+" \
    -zacutoff "+zacutoff+" \
    -maxlen "+pept_length+" -minlen "+pept_length+" \
    -tseq "+directives['fasta']+" \
    -poolconfile "+poolconfile+" \
    -allconfile "+allconfile+" \
    -poolaafile "+poolaafile+" \
    -allaafile "+allaafile+" \
    -probprod -log -delpred")
    
    run_commandline("python2.7 -m scoop -n "+args.CPUs+" "+CHAINS_BIN_DIR+"/align_parallel.py \
    -tseq "+directives['fasta']+" \
    -pseq peptides."+pept_length+"mers.fasta \
    -plist peptides."+pept_length+"mers.list "+restraints)
    
    for fname in ["consensus_alignment.iteration7", "consensus_alignment.overlapped_chains.iteration7",
                  "consensus_alignment.overlapped_chains_common_sequence.iteration7", "consensus_alignment.prediction_scores.iteration7",
                  "results_summary"]:
        os.rename(fname, fname+"."+pept_length+"mers_round"+cur_round)
    
    
    run_commandline(CHAINS_BIN_DIR+"/modify_connectivities_and_aatypes.py \
        -rstfile results_summary."+pept_length+"mers_round"+cur_round+" \
        -poolconfile "+poolconfile+" \
        -allconfile "+allconfile+" \
        -poolaafile "+poolaafile+" \
        -allaafile "+allaafile+" \
        -keepgly -keepala")
        
    os.rename(poolconfile+".mod", "connectivities_all.pool."+pept_length+"mers_round"+cur_round)
    os.rename(allconfile+".mod", "connectivities_all."+pept_length+"mers_round"+cur_round)
    os.rename(poolaafile+".mod", "amino_acid_type_prediction_probabilities.pool."+pept_length+"mers_round"+cur_round)
    os.rename(allaafile+".mod", "amino_acid_type_prediction_probabilities."+pept_length+"mers_round"+cur_round)
        
    
    for pept_length in reversed(range(int(last_length), int(first_length))):
        
        prev_round = str(cur_round)
        cur_round = str(int(cur_round)+1)
        prev_pept_length = str(pept_length+1)
        pept_length = str(pept_length)
        
        run_commandline("python2.7 -m scoop -n "+args.CPUs+" "+CHAINS_BIN_DIR+"/4D_assignment_tree-based_parallel.py \
        -root "+directives['HSQC']+" \
        -tocsy "+directives['4DTOCSY']+" \
        -noesy "+directives['4DNOESY']+" \
        -mcutoff "+mcutoff+" \
        -zmcutoff "+zmcutoff+" \
        -zacutoff "+zacutoff+" \
        -maxlen "+pept_length+" -minlen "+pept_length+" \
        -tseq "+directives['fasta']+" \
        -poolconfile connectivities_all.pool."+prev_pept_length+"mers_round"+prev_round+" \
        -poolaafile amino_acid_type_prediction_probabilities.pool."+prev_pept_length+"mers_round"+prev_round+" \
        -allconfile connectivities_all."+prev_pept_length+"mers_round"+prev_round+" \
        -allaafile amino_acid_type_prediction_probabilities."+prev_pept_length+"mers_round"+prev_round+" \
        -probprod -log -delpred")
        
        run_commandline("python2.7 -m scoop -n "+args.CPUs+" "+CHAINS_BIN_DIR+"/align_parallel.py \
        -tseq "+directives['fasta']+" \
        -pseq peptides."+pept_length+"mers.fasta \
        -plist peptides."+pept_length+"mers.list \
        -rstfile results_summary."+prev_pept_length+"mers_round"+prev_round)
        
        for fname in ["consensus_alignment.iteration7", "consensus_alignment.overlapped_chains.iteration7",
                      "consensus_alignment.overlapped_chains_common_sequence.iteration7",
                      "consensus_alignment.prediction_scores.iteration7", "results_summary"]:
            os.rename(fname, fname+"."+pept_length+"mers_round"+cur_round)
    
    
        run_commandline(CHAINS_BIN_DIR+"/modify_connectivities_and_aatypes.py \
        -rstfile results_summary."+pept_length+"mers_round"+cur_round+" \
        -poolconfile connectivities_all.pool."+prev_pept_length+"mers_round"+prev_round+" \
        -poolaafile amino_acid_type_prediction_probabilities.pool."+prev_pept_length+"mers_round"+prev_round+" \
        -allconfile connectivities_all."+prev_pept_length+"mers_round"+prev_round+" \
        -allaafile amino_acid_type_prediction_probabilities."+prev_pept_length+"mers_round"+prev_round+" \
        -keepgly -keepala") #-addconn
        
        os.rename("connectivities_all.pool."+prev_pept_length+"mers_round"+prev_round+".mod", "connectivities_all.pool."+pept_length+"mers_round"+cur_round)
        os.rename("amino_acid_type_prediction_probabilities.pool."+prev_pept_length+"mers_round"+prev_round+".mod", "amino_acid_type_prediction_probabilities.pool."+pept_length+"mers_round"+cur_round)
        os.rename("connectivities_all."+prev_pept_length+"mers_round"+prev_round+".mod", "connectivities_all."+pept_length+"mers_round"+cur_round)
        os.rename("amino_acid_type_prediction_probabilities."+prev_pept_length+"mers_round"+prev_round+".mod", "amino_acid_type_prediction_probabilities."+pept_length+"mers_round"+cur_round)



def chain_linker(absfile, poolconfile, patch=False):

    extra_args = ""
    if patch:
        extra_args = " -patch"
    run_commandline(CHAINS_BIN_DIR+"/chain_linker.py \
    -absfile "+absfile+" \
    -poolconfile "+poolconfile+" \
    -multicon \
    -mcincr 0.1 \
    -mcmin 0.5 " + extra_args)



ALLOW_PROOFREADING = True ; # If the correct rstart has not been specified, or the program could not guess it, then proofreading will be deactivated
os.chdir(args.WORKDIR)
if not os.path.exists('4DCHAINS_workdir'):
    os.mkdir('4DCHAINS_workdir')
os.chdir('4DCHAINS_workdir')

if directives['doNHmapping'] == True:
    print bcolors.BOLDGREEN + "Annotating {N-H}-HSQC root file." + bcolors.ENDBOLD
    if not directives['rstart']:    # if not gived -rstart, try to guess it from the labels in the HSQC file
        directives['rstart'] = str(guess_starting_resid(directives['HSQC'], directives['fasta']))
        if directives['rstart'] == 'None':
            directives['rstart'] = '1';
            ALLOW_PROOFREADING = False
    run_commandline(CHAINS_BIN_DIR+"/annotate_root.py -root " + directives['HSQC']+ " -rstart " + directives['rstart'] + " -tseq " \
                    + directives['fasta'] + " -o " + directives['HSQC'] + "_annotated", logname="annotate_root.log")
    directives['HSQC'] = directives['HSQC'] + "_annotated"
    
    total_cycles = len(directives['mcutoff'])
    directives['rst_file'] = [None]*total_cycles
    directives['poolconfile'] =  [None]*total_cycles
    directives['allconfile'] =  [None]*total_cycles
    directives['poolaafile'] =  [None]*total_cycles
    directives['allaafile'] =  [None]*total_cycles
    
    
    if directives['con_file'] and directives['aa_file']:
        print bcolors.BOLDGREEN + "Using connectivities file", directives['con_file'], "and aa-type predictions file", directives['aa_file'] + bcolors.ENDBOLD
        directives['con_file'] = os.path.realpath(directives['con_file'])
        directives['aa_file'] = os.path.realpath(directives['aa_file'])
        shutil.copy(directives['con_file'], directives['con_file']+".pool")
        shutil.copy(directives['aa_file'], directives['aa_file']+".pool")
        directives['poolconfile'][0] = directives['con_file'] + ".pool"
        directives['allconfile'][0] = directives['con_file']
        directives['poolaafile'][0] = directives['aa_file'] + ".pool"
        directives['allaafile'][0] = directives['aa_file']
    else:
        print bcolors.BOLDGREEN + "Generating connectivities and amino-acid type predictions for each TOCSY index group." + bcolors.ENDBOLD
        generate_input_data() # temporarily deactivated for debugging
        directives['poolconfile'][0] = args.WORKDIR+"/4DCHAINS_workdir/connectivities_all.pool"
        directives['allconfile'][0] = args.WORKDIR+"/4DCHAINS_workdir/connectivities_all"
        directives['poolaafile'][0] = args.WORKDIR+"/4DCHAINS_workdir/amino_acid_type_prediction_probabilities.pool"
        directives['allaafile'][0] = args.WORKDIR+"/4DCHAINS_workdir/amino_acid_type_prediction_probabilities"
    
    
    for cycle in range(1, total_cycles+1):
        
        print bcolors.BOLDGREEN + "Entering NH-mapping cycle", cycle, bcolors.ENDBOLD
        
        if cycle > 1 and directives['rst_from_prev_cycle'][cycle-1]:
            suffix = str(directives['last_length'][cycle-2]) + "mers_round" + str(directives['first_length'][cycle-2]-directives['last_length'][cycle-2]+1)
            directives['rst_file'][cycle-1] = args.WORKDIR+"/4DCHAINS_workdir/cycle"+str(cycle-1)+"/results_summary." + suffix
            directives['poolconfile'][cycle-1] = args.WORKDIR+"/4DCHAINS_workdir/cycle"+str(cycle-1)+"/connectivities_all.pool." + suffix
            directives['allconfile'][cycle-1] = args.WORKDIR+"/4DCHAINS_workdir/cycle"+str(cycle-1)+"/connectivities_all." + suffix
            directives['poolaafile'][cycle-1] = args.WORKDIR+"/4DCHAINS_workdir/cycle"+str(cycle-1)+"/amino_acid_type_prediction_probabilities.pool." + suffix
            directives['allaafile'][cycle-1] = args.WORKDIR+"/4DCHAINS_workdir/cycle"+str(cycle-1)+"/amino_acid_type_prediction_probabilities." + suffix
        
        if os.path.exists("cycle"+str(cycle)):
            shutil.rmtree("cycle"+str(cycle))
        os.mkdir("cycle"+str(cycle))
        os.chdir("cycle"+str(cycle))
        multiround_doNHmapping(cycle=cycle,
                            tolH=directives['tolH'],
                            tolC=directives['tolC'],
                            mcutoff=directives['mcutoff'][cycle-1],
                            zmcutoff=directives['zmcutoff'][cycle-1],
                            zacutoff=directives['zacutoff'][cycle-1],
                            poolconfile=directives['poolconfile'][cycle-1],
                            allconfile=directives['allconfile'][cycle-1],
                            poolaafile=directives['poolaafile'][cycle-1],
                            allaafile=directives['allaafile'][cycle-1],
                            first_length=directives['first_length'][cycle-1],
                            last_length=directives['last_length'][cycle-1],
                            rst_file=directives['rst_file'][cycle-1])
        
        os.chdir("../")
    
    
    print bcolors.BOLDGREEN + "Linking individual contigs in the alignment using only the available connectivities (not the amino-acid type predictions)." + bcolors.ENDBOLD
    os.chdir("cycle"+str(cycle))
    suffix = str(directives['last_length'][cycle-1]) + "mers_round" + str(directives['first_length'][cycle-1]-directives['last_length'][cycle-1]+1)
    chain_linker("results_summary." + suffix, "connectivities_all.pool."+ suffix )
    
    suffix = str(directives['last_length'][cycle-1]) + "mers_round" + str(directives['first_length'][cycle-1]-directives['last_length'][cycle-1]+1)
    NHmap_fname = "results_summary." + suffix+".chainlinkers"
    if directives['patch']:
        NHmap_fname += ".patch"
    directives['NHmap'] = os.path.abspath(NHmap_fname)
    print bcolors.BOLDGREEN + "Transfering assignments to {N-H}-HSQC root file." + bcolors.ENDBOLD
    run_commandline(CHAINS_BIN_DIR+"/annotate_root.py -root "+directives['HSQC']+" -absfile "+directives['NHmap']+" -rstart " +\
                    directives['rstart'] + " -o " + args.WORKDIR + "/4DCHAINS_assigned_NH_HSQC.sparky",
                    logname="annotate_root.log")
    
    os.chdir(args.WORKDIR)
    if os.path.exists('4DCHAINS_NHmap'):
        os.unlink('4DCHAINS_NHmap')
    os.symlink(directives['NHmap'], "4DCHAINS_NHmap")
    
if directives['doassign4DTOCSY'] == True:
    print bcolors.BOLDGREEN + "Assigning 4D-TOCSY peaks." + bcolors.ENDBOLD
    os.chdir(args.WORKDIR + "/4DCHAINS_workdir")
    if not os.path.exists('4DTOCSY_cs_assignment'):
        os.mkdir('4DTOCSY_cs_assignment')
    else:
        shutil.rmtree('4DTOCSY_cs_assignment')
        os.mkdir('4DTOCSY_cs_assignment')
    os.chdir('4DTOCSY_cs_assignment')
    
    if not directives['NHmap']:
        print bcolors.WARNING + "WARNING: you did not specify NHmap file ('NHmap' directive is missing)! Trying to find the '4DCHAINS_NHmap' file\
 under your workdir." + bcolors.ENDC
        if not os.path.exists(args.WORKDIR + "/4DCHAINS_NHmap"):
            print bcolors.FAIL + "ERROR: '4DCHAINS_NHmap' doesn't exist! Please specify the 'NHmap' directive in your protocol!" + bcolors.ENDC
            sys.exit(1)
        else:
             directives['NHmap'] = args.WORKDIR + "/4DCHAINS_NHmap"
    os.symlink(directives['NHmap'], 'NHmap')
    if not directives['4DTOCSY_assignedNH']:
        directives['4DTOCSY_assignedNH'] = directives['4DTOCSY'].replace(".list", "num.list")
    if not directives['4DNOESY_assignedNH']:
        directives['4DNOESY_assignedNH'] = directives['4DNOESY'].replace(".list", "num.list")
    os.symlink(directives['4DTOCSY_assignedNH'], '4DTOCSY_assignedNH')
    os.symlink(directives['4DNOESY_assignedNH'], '4DNOESY_assignedNH')
    
    if not directives['rstart'] and directives['HSQC'] :    # if not gived -rstart, try to guess it from the labels in the HSQC file
        directives['rstart'] = str(guess_starting_resid(directives['HSQC'], "", directives['NHmap']))
        if directives['rstart'] == 'None':
            directives['rstart'] = '1'
            ALLOW_PROOFREADING = False
    elif not directives['rstart']:
        directives['rstart'] = '1'
    
    run_commandline(CHAINS_BIN_DIR+"/cs_assignment.py \
                    -absfile NHmap \
                    -rstart "+directives['rstart']+" \
                    -tocsy 4DTOCSY_assignedNH \
                    -probprod \
                    -2dhist \
                    -noesy 4DNOESY_assignedNH \
                    -probmode 2 -o 4DTOCSY_assignedall", "cs_assignment.log")
    
    directives['4DTOCSY_assignedall'] = os.path.abspath('4DTOCSY_assignedall.sparky')
    if os.path.exists(args.WORKDIR + '/4DTOCSY_assignedall.sparky'):    # copy the sparky file
        os.unlink(args.WORKDIR + '/4DTOCSY_assignedall.sparky')
    shutil.copyfile('4DTOCSY_assignedall.sparky', args.WORKDIR + '/4DTOCSY_assignedall.sparky')
    if os.path.exists(args.WORKDIR + '/4DTOCSY_assignedall.xeasy'):     # copy the xeasy file
        os.unlink(args.WORKDIR + '/4DTOCSY_assignedall.xeasy')
    shutil.copyfile('4DTOCSY_assignedall.xeasy', args.WORKDIR + '/4DTOCSY_assignedall.xeasy')

if directives['doassign4DNOESY'] == True:
    print bcolors.BOLDGREEN + "Assigning 4D-NOESY peaks." + bcolors.ENDBOLD
    os.chdir(args.WORKDIR + "/4DCHAINS_workdir")
    if not os.path.exists('4DNOESY_cs_assignment'):
        os.mkdir('4DNOESY_cs_assignment')
    else:
        shutil.rmtree('4DNOESY_cs_assignment')
        os.mkdir('4DNOESY_cs_assignment')
    os.chdir('4DNOESY_cs_assignment')

    if not directives['NHmap']:
        print bcolors.WARNING + "WARNING: you did not specify NHmap file ('NHmap' directive is missing)! Trying to find the '4DCHAINS_NHmap' file\
 under your workdir." + bcolors.ENDC
        if not os.path.exists(args.WORKDIR + "/4DCHAINS_NHmap"):
            print bcolors.FAIL + "ERROR: '4DCHAINS_NHmap' doesn't exist! Please specify the 'NHmap' directive in your protocol!" + bcolors.ENDC
            sys.exit(1)
        else:
             directives['NHmap'] = args.WORKDIR + "/4DCHAINS_NHmap"
    os.symlink(directives['NHmap'], 'NHmap')
    
    if not directives['4DTOCSY_assignedall']:
        print bcolors.WARNING + "WARNING: you did not specify an assigned 4D-TOCSY file file ('4DTOCSY_assignedall' directive is missing)! Trying to \
find the '4DTOCSY_assignedall.sparky' file under your workdir." + bcolors.ENDC
        if not os.path.exists(args.WORKDIR + "/4DTOCSY_assignedall.sparky"):
            print bcolors.FAIL + "ERROR: '4DTOCSY_assignedall.sparky' doesn't exist! Please specify the '4DTOCSY_assignedall' directive in your protocol!" + bcolors.ENDC
            sys.exit(1)
        else:
             directives['4DTOCSY_assignedall'] = args.WORKDIR + "/4DTOCSY_assignedall.sparky"
    os.symlink(directives['4DTOCSY_assignedall'], '4DTOCSY_assignedall.sparky')
    
    if not directives['4DNOESY_assignedNH']:
        directives['4DNOESY_assignedNH'] = directives['4DNOESY'].replace(".list", "num.list")
    os.symlink(directives['4DNOESY_assignedNH'], '4DNOESY_assignedNH')
    
    proofread = ""
    if directives['user_4DTOCSY_assignedall']:
        os.symlink(directives['user_4DTOCSY_assignedall'], 'user_4DTOCSY_assignedall.sparky')
        proofread += " -usertocsy user_4DTOCSY_assignedall.sparky"
    if directives['user_4DNOESY_assignedall']:
        os.symlink(directives['user_4DNOESY_assignedall'], 'user_4DNOESY_assignedall')
        proofread += " -usernoesy user_4DNOESY_assignedall"
   
    if not directives['rstart'] and directives['HSQC'] :    # if not gived -rstart, try to guess it from the labels in the HSQC file
        directives['rstart'] = str(guess_starting_resid(directives['HSQC'], "", directives['NHmap']))
        if directives['rstart'] == 'None':
            directives['rstart'] = '1'
            ALLOW_PROOFREADING = False
    elif not directives['rstart']:
        directives['rstart'] = '1'
    
    if ALLOW_PROOFREADING == False:
        proofread = ""
    
    run_commandline(CHAINS_BIN_DIR+"/cs_assignment_NOESY.py \
                    -absfile NHmap \
                    -rstart "+directives['rstart']+" \
                    -atocsy 4DTOCSY_assignedall.sparky \
                    -noesy 4DNOESY_assignedNH \
                    -probprod \
                    -2dhist \
                    -int1 -int2 -int3 -int4 \
                    -ithres1 0.1 \
                    -ithres2 0.1 \
                    -ithres3 0.1 \
                    -ithres4 0.1 \
                    -wcthres1 0 \
                    -bcthres1 0 \
                    -wcthres2 0 \
                    -bcthres2 0 \
                    -wcthres1 0 \
                    -bcthres1 0 \
                    -wcthres2 0 \
                    -bcthres2 0 \
                    -percentile1 0 \
                    -percentile2 0.8 \
                    -percentile3 0.9 \
                    -percentile4 0.9 \
                    -percentile5 0.8 \
                    -itrans1 4 -itrans2 1 -itrans3 4 -itrans4 1 -itrans5 1 -o 4DNOESY_assignedall " + proofread, "cs_assignment_NOESY.log")
    
    if os.path.exists(args.WORKDIR + '/4DNOESY_assignedall.sparky'):    # copy the sparky file
        os.unlink(args.WORKDIR + '/4DNOESY_assignedall.sparky')
    shutil.copyfile('4DNOESY_assignedall.sparky', args.WORKDIR + '/4DNOESY_assignedall.sparky')
    if os.path.exists('4DNOESY_assignedall.proofread.xeasy'):
        if os.path.exists(args.WORKDIR + '/4DNOESY_assignedall.proofread.xeasy'):     # copy the xeasy file
            os.unlink(args.WORKDIR + '/4DNOESY_assignedall.proofread.xeasy')
        shutil.copyfile('4DNOESY_assignedall.proofread.xeasy', args.WORKDIR + '/4DNOESY_assignedall.proofread.xeasy')
    else:
        if os.path.exists(args.WORKDIR + '/4DNOESY_assignedall.xeasy'):     # copy the xeasy file
            os.unlink(args.WORKDIR + '/4DNOESY_assignedall.xeasy')
        shutil.copyfile('4DNOESY_assignedall.xeasy', args.WORKDIR + '/4DNOESY_assignedall.xeasy')
    

elif directives['doassignonly4DNOESY'] == True:
    print bcolors.BOLDGREEN + "Assigning 4D-NOESY peaks without using 4D-TOCSY assignments." + bcolors.ENDBOLD
    os.chdir(args.WORKDIR + "/4DCHAINS_workdir")
    if not os.path.exists('only4DNOESY_cs_assignment'):
        os.mkdir('only4DNOESY_cs_assignment')
    else:
        shutil.rmtree('only4DNOESY_cs_assignment')
        os.mkdir('only4DNOESY_cs_assignment')
    os.chdir('only4DNOESY_cs_assignment')

    if not directives['NHmap']:
        print bcolors.WARNING + "WARNING: you did not specify NHmap file ('NHmap' directive is missing)! Trying to find the '4DCHAINS_NHmap' file\
 under your workdir." + bcolors.ENDC
        if not os.path.exists(args.WORKDIR + "/4DCHAINS_NHmap"):
            print bcolors.FAIL + "ERROR: '4DCHAINS_NHmap' doesn't exist! Please specify the 'NHmap' directive in your protocol!" + bcolors.ENDC
            sys.exit(1)
        else:
             directives['NHmap'] = args.WORKDIR + "/4DCHAINS_NHmap"
    os.symlink(directives['NHmap'], 'NHmap')
    if not directives['4DNOESY_assignedNH']:
        directives['4DNOESY_assignedNH'] = directives['4DNOESY'].replace(".list", "num.list")
    os.symlink(directives['4DNOESY_assignedNH'], '4DNOESY_assignedNH')
    
    proofread = ""
    if directives['user_4DNOESY_assignedall']:
        os.symlink(directives['user_4DNOESY_assignedall'], 'user_4DNOESY_assignedall')
        proofread += " -usernoesy user_4DNOESY_assignedall"
    if directives['user_4DTOCSY_assignedall']:
        os.symlink(directives['user_4DTOCSY_assignedall'], 'user_4DTOCSY_assignedall.sparky')
        proofread += " -usertocsy user_4DTOCSY_assignedall.sparky"
   
    if not directives['rstart'] and directives['HSQC'] :    # if not gived -rstart, try to guess it from the labels in the HSQC file
        directives['rstart'] = str(guess_starting_resid(directives['HSQC'], "", directives['NHmap']))
        if directives['rstart'] == 'None':
            directives['rstart'] = '1'
            ALLOW_PROOFREADING = False
    elif not directives['rstart']:
        directives['rstart'] = '1'
    
    if ALLOW_PROOFREADING == False:
        proofread = ""
    
    run_commandline(CHAINS_BIN_DIR+"/cs_assignment_onlyNOESY.py \
                    -absfile NHmap \
                    -rstart "+directives['rstart']+" \
                    -noesy 4DNOESY_assignedNH \
                    -probprod \
                    -2dhist \
                    -int1 -int2 -int3 -int4 -int5 \
                    -ithres1 0.1 -ithres2 0.1 -ithres3 0.1 -ithres4 0 -ithres5 0 \
                    -percentile1 0.85 -percentile2 0.9 -percentile3 0.8 -percentile4 0.8 -percentile5 0.8 \
                    -itrans1 4 -itrans2 2 -itrans3 3 -itrans4 2 -itrans5 1 -o only4DNOESY_assignedall \
                    " + proofread, "cs_assignment_onlyNOESY.log")
    
    if os.path.exists(args.WORKDIR + '/only4DNOESY_assignedall.sparky'):    # copy the sparky file
        os.unlink(args.WORKDIR + '/only4DNOESY_assignedall.sparky')
    shutil.copyfile('only4DNOESY_assignedall.sparky', args.WORKDIR + '/only4DNOESY_assignedall.sparky')
    
    if os.path.exists('only4DNOESY_assignedall.proofread.xeasy'):
        if os.path.exists(args.WORKDIR + '/only4DNOESY_assignedall.proofread.xeasy'):     # copy the xeasy file
            os.unlink(args.WORKDIR + '/only4DNOESY_assignedall.proofread.xeasy')
        shutil.copyfile('only4DNOESY_assignedall.proofread.xeasy', args.WORKDIR + '/only4DNOESY_assignedall.proofread.xeasy')
    else:
        if os.path.exists(args.WORKDIR + '/only4DNOESY_assignedall.xeasy'):     # copy the xeasy file
            os.unlink(args.WORKDIR + '/only4DNOESY_assignedall.xeasy')
        shutil.copyfile('only4DNOESY_assignedall.xeasy', args.WORKDIR + '/only4DNOESY_assignedall.xeasy')
    