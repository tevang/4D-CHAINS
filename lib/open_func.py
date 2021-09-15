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



import sys, os, subprocess, time
from scoop import futures, shared

# Import 4D-CHAINS libraries
CHAINS_LIB_DIR = os.path.dirname(os.path.realpath(__file__))
CHAINS_BIN_DIR = CHAINS_LIB_DIR[:-3] + "bin"
sys.path.append(CHAINS_BIN_DIR)
sys.path.append(CHAINS_LIB_DIR)
from .global_func import *


class bcolors:
    ## See also https://pypi.python.org/pypi/blessings/   https://pypi.python.org/pypi/colorama
    HEADER = '\033[1m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[32m'
    BOLDGREEN = '\x1b[32;1m'
    BOLDBLUE = '\x1b[34;1m'
    BOLDRED = '\x1b[91;1m'
    WARNING = '\033[43m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    ENDBOLD = '\x1b[0m'

    def disable(self):
        self.HEADER = ''
        self.OKBLUE = ''
        self.OKGREEN = ''
        self.BOLDGREEN = ''
        self.BOLDRED = ''
        self.WARNING = ''
        self.FAIL = ''
        self.ENDC = ''
        self.ENDBOLD = ''
        

class ProgressBar:   
    
    BARLENGTH = 10
    
    def __init__(self,barLength):
        self.BARLENGTH = barLength # initialize the length of the progress bar
        
    def set_progress(self, progress):
        # update_progress() : Displays or updates a console progress bar
        ## Accepts a float between 0 and 1. Any int will be converted to a float.
        ## A value under 0 represents a 'halt'.
        ## A value at 1 or bigger represents 100%
        status = ""
        if isinstance(progress, int):
            progress = float(progress)
        if not isinstance(progress, float):
            progress = 0
            status = "error: progress var must be float\r\n"
        if progress < 0:
            progress = 0
            status = "Halt...\r\n"
        if progress >= 1:
            progress = 1
            status = "Done...\r\n"
        block = int(round(self.BARLENGTH*progress))
        text = "\rPercent: [{0}] {1}% {2}".format( "#"*block + "-"*(self.BARLENGTH-block), progress*100, status)
        sys.stdout.write(text)
        sys.stdout.flush()
        

# Example of using pebble instead of scoop for multi-threading calculations
"""
        from pebble import ProcessPool
        from concurrent.futures import TimeoutError
        import copy_reg
        import types
        import multiprocessing

        def _pickle_method(m):
            if m.im_self is None:
                return getattr, (m.im_class, m.im_func.func_name)
            else:
                return getattr, (m.im_self, m.im_func.func_name)

        copy_reg.pickle(types.MethodType, _pickle_method)

        results = []

        with ProcessPool(max_workers=8) as pool:
            future = pool.map(self.expand_tree,
                              all_perms,
                              [length]*len(all_perms),
                              [self.connections_dict]*len(all_perms),
                              timeout=20)

            iterator = future.result()

            # iterate over all results, if a computation timed out
            # print it and continue to the next result
            while True:
                try:
                    result = next(iterator)
                    results.append(result)
                except StopIteration:
                    break
                except TimeoutError as error:
                    print "function took longer than %d seconds" % error.args[1]

        print results
"""

def align_and_read(peptide_fasta_file, peptide_list_file, chunk_No):
    """
        MAIN FUNCTION to align the peptides to the protein sequence and load the alignment file ("*.needle").
    """
    
    try:
        
        progbar = shared.getConst('PROGBAR')
        aln_index = shared.getConst('ALN_INDEX')
        args = shared.getConst('ARGS')
        CHAINS_BIN_DIR = shared.getConst('chains_bin_dir')
        
        ## measure lines of peptides list for Progressbar
        result = subprocess.check_output(['wc', '-l', peptide_fasta_file])
        line_num = int(result.split()[0])
        peptide_num = line_num/2
    
        if args.DO_ALIGN:
            os.environ['EMBOSS_ACDROOT'] = CHAINS_BIN_DIR + "/../databases/acd"
            ## ALIGN THE PEPTIDES WITH THE TEMPLATE
            print("DEBUG: " , CHAINS_BIN_DIR+"/needle", args.template_sequence_file, peptide_fasta_file, "-datafile", CHAINS_BIN_DIR+"/../databases/acd/"+args.MATRIX, \
                                    "-gapopen", str(args.gap_open_penalty), "-gapextend", str(args.gap_extension_penalty), "-outfile", \
                                    "peptide_alignment.needle.chunk"+str(chunk_No))
            result = subprocess.call([CHAINS_BIN_DIR+"/needle", args.template_sequence_file, peptide_fasta_file, "-datafile", CHAINS_BIN_DIR+"/../databases/acd/"+args.MATRIX,
                                    "-gapopen", str(args.gap_open_penalty), "-gapextend", str(args.gap_extension_penalty), "-outfile",
                                    "peptide_alignment.needle.chunk"+str(chunk_No)])
        
        ## ANALYZE ALIGNMENT FILE
        #progbar = ProgressBar(100)
        peptide_lengthsum_dict = {}
        peptide_identity_dict = {}
        peptide_similarity_dict = {}
        peptide_gapnum_dict = {}
        peptide_score_dict = {}
        peptide_alignment_dict = {}
        #aln_index = 0
        max_peptide_length = 0
        peptides2remove_set = set()
        peptide_name = None
        if args.peptide_alignment_file:
            print("Loading alignment file "+args.peptide_alignment_file+".chunk"+str(chunk_No)+"...")
            f = open(args.peptide_alignment_file+".chunk"+str(chunk_No), 'r')
        else:
            print("Loading alignment file "+"peptide_alignment.needle.chunk"+str(chunk_No)+"...")
            f = open("peptide_alignment.needle.chunk"+str(chunk_No), 'r')
        time.sleep(5)   # wait for 5 sec for the workers to open all alignment files and then initialize the progress bar
        for line in f:
            #print line
            if line[0:18] == "# Aligned_sequence":
                line = next(f)
                try:
                    template_name = line.split()[2]
                except IndexError:
                    template_name = ""
                line = next(f)
                peptide_name = line.split()[2]
                #print "DEBUG: peptide_name = ",peptide_name
                peptide_length = int(peptide_name.split('_')[1])
                if peptide_length > max_peptide_length:
                    max_peptide_length = peptide_length
                line = next(f); line = next(f); line = next(f); line = next(f); # discard next 4 lines
                line = next(f) # read Length line
            #elif line[0:9] == "# Length:":
                mo = re.search('# Length:\s+([0-9]+)', line)
                if mo:
                    LENGTH_SUM = int(mo.group(1))
                    peptide_lengthsum_dict[peptide_name] = LENGTH_SUM
            #elif line[0:10] == "# Identity":
                line = next(f) # read Identity line
                mo = re.search('# Identity:\s+([0-9]+)/.*', line)
                if mo:
                    IDENTITY_NUM = float(mo.group(1))
                    peptide_identity_dict[peptide_name] = IDENTITY_NUM/peptide_length
            #elif line[0:12] == "# Similarity":
                line = next(f) # read Similarity line
                mo = re.search('# Similarity:\s+([0-9]+)/.*', line)
                if mo:
                    SIMILARITY_NUM = float(mo.group(1))
                    peptide_similarity_dict[peptide_name] = SIMILARITY_NUM/peptide_length
            #elif line[0:6] == "# Gaps":
                line = next(f) # read Gaps line
                mo = re.search('# Gaps:\s+([0-9]+)/.*', line)
                if mo:
                    GAP_NUM = int(mo.group(1))
                    peptide_gapnum_dict[peptide_name] = GAP_NUM
            #elif line[0:7] == "# Score":
                line = next(f) # read Score line
                mo = re.search('# Score: ([0-9.]+)', line)
                if mo:
                    SCORE = float(mo.group(1))
                    peptide_score_dict[peptide_name] = SCORE
            elif line[0] != "#" and peptide_name != None:
                try:
                    try:
                        template_alignment = ""
                        peptide_alignment = ""
                        line = next(f)
                        while line[0] != "#":
                            if re.match(template_name+"\s+[0-9]+\s+[A-Z-]+\s+[0-9]+", line):
                                template_alignment += line.split()[2]
                            elif re.match(peptide_name+"\s+[0-9]+\s+[A-Z-]+\s+[0-9]+", line):
                                peptide_alignment += line.split()[2]
                            line = next(f)
                    except StopIteration:
                        pass
                    if IDENTITY_NUM/peptide_length >= args.IDENTITY_CUTOFF: # save the alignments only of the identity is above the cutoff
                        if peptide_alignment[-1] == '-':    # DISCARD ALIGNMENTS WHERE AN EXTRA C-TERMINAL TINDEX IS NEEDED
                            peptide_alignment_dict[peptide_name] = (template_alignment, peptide_alignment)
                    aln_index += 1
                    progbar.set_progress((float(aln_index)/peptide_num))
                except NameError:
                    pass
            try:
                if peptide_identity_dict[peptide_name] < args.IDENTITY_CUTOFF:  # IF THE SEQUENCE IDENTITY WAS BELOW THE SPECIFIED CUTOFF DELETE ALL ENTRIES TO SAVE MEMEORY
                    peptides2remove_set.add(peptide_name)
                    try:
                        del peptide_identity_dict[peptide_name]
                    except KeyError:
                        pass
                    try:
                        del peptide_lengthsum_dict[peptide_name]
                    except KeyError:
                        pass
                    try:
                        del peptide_similarity_dict[peptide_name]
                    except KeyError:
                        pass
                    try:
                        del peptide_gapnum_dict[peptide_name]
                    except KeyError:
                        pass
                    try:
                        del peptide_score_dict[peptide_name]
                    except KeyError:
                        pass
            except KeyError:
                        pass
            
        ### FINALLY TO MAKE SURE ALL DICTIONARIES CONTAIN THE SAME KEYS CLEAN THE REDUNDANT ONES
        for peptide_name in list(peptide_score_dict.keys()):
            if not peptide_name in list(peptide_alignment_dict.keys()):
                del peptide_identity_dict[peptide_name]
                del peptide_lengthsum_dict[peptide_name]
                del peptide_similarity_dict[peptide_name]
                del peptide_gapnum_dict[peptide_name]
                del peptide_score_dict[peptide_name]
        #keys1 = set(peptide_identity_dict.keys())
        #keys2 = set(peptide_lengthsum_dict.keys())
        #keys3 = set(peptide_similarity_dict.keys())
        #keys4 = set(peptide_gapnum_dict.keys())
        #keys5 = set(peptide_score_dict.keys())
        #keys6 = set(peptide_alignment_dict.keys())
        #peptides2keep_set = set.intersection(keys1, keys2, keys3, keys4, keys5, keys6)
        #peptides_not_removed_set = peptides2keep_set.intersection(peptides2remove_set)
        #for peptide_name in peptides2remove_set:
        #    print "DEBUG: removing peptide", peptide_name,"from peptide_alignment_dict"
        #    del peptide_alignment_dict[peptide_name]
        
        print("")    # change line after progress bar
        
        # INITIALIZE args.OVERLAP_CUTOFF
        if not args.OVERLAP_CUTOFF:
            args.OVERLAP_CUTOFF = max_peptide_length -1
        
        peptideName2chain_dict = {}
        peptideName_chainScoreList_dict = {}
        if args.COMPRESS: 
            f = bz2.BZ2File(peptide_list_file, 'r')
        else:
            f = open(peptide_list_file, 'r')
        for line in f:
            word_list = line.split()
            peptide_name = word_list[5]
            if peptide_name in peptides2remove_set: # ommit loading peptides with sequence identity below the cutoff
                continue
            peptide_sequence = word_list[7]
            chain = word_list[8]
            peptideName_chainScoreList_dict[peptide_name] = list(map(float, word_list[9].split(",")))
            peptideName2chain_dict[peptide_name] = chain
        f.close()
            
        #print "DEBUG: identity",len(peptide_identity_dict)
        #print "DEBUG: similarity",len(peptide_similarity_dict)
        #print "DEBUG: peptide_score_dict=",len(peptide_score_dict)
        #print "DEBUG: peptide_alignment_dict=", len(peptide_alignment_dict)
        
        # remove intermediate peptide fast and list files
        os.remove(peptide_fasta_file)
        os.remove(peptide_list_file)
        return (peptide_identity_dict, peptide_score_dict, peptideName2chain_dict, peptide_alignment_dict, peptideName_chainScoreList_dict, max_peptide_length)
    
    except:
        type, value, tb = sys.exc_info()
        lines = traceback.format_exception(type, value, tb)
        print((''.join(lines)))
        raise
    

def split_peptidefile(args):
    """
        FUNCTION to divide the input files (peptides.*.list and peptides.*.fasta) into args.CPUs parts to parallelize the calculations.
    """
    
    print("Splitting input files into "+str(args.CPUs)+" parts to speed up calculations...")
    result = subprocess.check_output(['wc', '-l', args.peptide_fasta_file])
    line_num = int(result.split()[0])
    peptide_num = line_num/2
    offset = 2*int(peptide_num/args.CPUs)
    stop_points_list = list(range(offset, line_num+1, offset))
    stop_points_list[-1] = line_num
    stop_points_iterator = iter(stop_points_list)
    lineNo2stop = next(stop_points_iterator)
    chunk_index = 1
    fasta_fout = open(args.peptide_fasta_file+".chunk"+str(chunk_index), 'w')
    with open(args.peptide_fasta_file, 'r') as f:
        for lineNo, line  in enumerate(f):
            lineNo += 1
            fasta_fout.write(line)
            if lineNo == lineNo2stop:
                fasta_fout.close()
                try:
                    lineNo2stop = next(stop_points_iterator)
                    chunk_index += 1
                    fasta_fout = open(args.peptide_fasta_file+".chunk"+str(chunk_index), 'w')
                except StopIteration:   # now write the lines that were left into the last chunk
                    break
        fasta_fout.close()
    
    chunk_index = 1
    offset = int(peptide_num/args.CPUs)
    result = subprocess.check_output(['wc', '-l', args.peptide_list_file])
    line_num = int(result.split()[0])
    stop_points_list = list(range(offset, line_num+1, offset))
    stop_points_list[-1] = line_num
    stop_points_iterator = iter(stop_points_list)
    lineNo2stop = next(stop_points_iterator)
    if args.COMPRESS:
        list_fout = bz2.BZ2File(args.peptide_list_file+".chunk"+str(chunk_index)+".bz2", 'w')
    else:
        list_fout = open(args.peptide_list_file+".chunk"+str(chunk_index), 'w')
    with open(args.peptide_list_file, 'r') as f:
        for lineNo, line  in enumerate(f):
            lineNo += 1
            list_fout.write(line)
            if lineNo == lineNo2stop:
                list_fout.close()
                try:
                    lineNo2stop = next(stop_points_iterator)
                    chunk_index += 1
                    if args.COMPRESS:
                        list_fout = bz2.BZ2File(args.peptide_list_file+".chunk"+str(chunk_index)+".bz2", 'w')
                    else:
                        list_fout = open(args.peptide_list_file+".chunk"+str(chunk_index), 'w')
                except StopIteration:
                    break
        list_fout.close()


