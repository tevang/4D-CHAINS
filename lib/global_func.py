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


import bz2
import collections
import ftplib
import pickle
import re
from collections import OrderedDict
from decimal import getcontext, Decimal, InvalidOperation
from subprocess import call
import numpy as np

from .open_func import *
import regex
from inspect import getframeinfo, stack

import regex

from .open_func import *

getcontext().prec = 7


class tree(OrderedDict):
    def __missing__(self, key):
        self[key] = type(self)()
        return self[key]

class bcolors:
    ## See also https://pypi.python.org/pypi/blessings/   https://pypi.python.org/pypi/colorama
    HEADER = '\033[1m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[32m'
    BOLDGREEN = '\x1b[32;1m'
    BOLDBLUE = '\x1b[34;1m'
    FAIL = '\x1b[91;1m'
    WARNING = '\033[43m'
    OKRED = '\033[91m'
    ENDC = '\033[0m'
    ENDBOLD = '\x1b[0m'

    def disable(self):
        self.HEADER = ''
        self.OKBLUE = ''
        self.OKGREEN = ''
        self.BOLDGREEN = ''
        self.OKRED = ''
        self.WARNING = ''
        self.FAIL = ''
        self.ENDC = ''
        self.ENDBOLD = ''

class ColorPrint():

    def __init__(self, message, type):
        prefix = {
            "HEADER": bcolors.HEADER,
            "OKBLUE": bcolors.OKBLUE,
            "OKGREEN": bcolors.OKGREEN,
            "BOLDGREEN": bcolors.BOLDGREEN,
            "BOLDBLUE": bcolors.BOLDBLUE,
            "OKRED": bcolors.OKRED,
            "WARNING": bcolors.WARNING,
            "FAIL": bcolors.FAIL,
        }
        suffix = {
            "HEADER" : bcolors.ENDBOLD,
            "OKBLUE" : bcolors.ENDC,
            "OKGREEN" : bcolors.ENDC,
            "BOLDGREEN" : bcolors.ENDBOLD,
            "BOLDBLUE" : bcolors.ENDBOLD,
            "OKRED" : bcolors.ENDC,
            "WARNING" : bcolors.ENDC,
            "FAIL" : bcolors.ENDBOLD,
        }
        print(prefix[type] + message + suffix[type])

def Debuginfo(message, fail=False):
    caller = getframeinfo(stack()[1][0])
    if fail:
        ColorPrint("%s:%d - %s" % (caller.filename, caller.lineno, message), "FAIL")
    else:
        print("%s:%d - %s" % (caller.filename, caller.lineno, message))

## OLD FUNCTION DEFINITION, tolerance is proportional to the difference between the two values
def approx_equal_proportional(x, y, tolerance=0.001):
    return abs(x-y) <= 0.5 * tolerance * (x + y)

def approx_equal(x, y, tolerance=0.001):
    return abs(float(x)-float(y)) <= tolerance


def list_files(folder, pattern, full_path=False):
    """
        FUNCTION to list the files in 'folder' that match the 'pattern'.
    """

    fpaths = os.listdir(folder)
    fpattern = re.compile(pattern)
    file_list = list(filter(fpattern.search, fpaths))
    if full_path:
        file_list = [folder + "/" + f for f in file_list]
    return file_list

def download_CS_histograms():
    """
    FUNCTION to download the most recent chemical shift histograms from BMRB 
    """
    
    global aa3to1_dict, aa1to3_dict
    
    if not os.path.exists(CHAINS_BIN_DIR+"/../BMRB_data"):
            os.makedirs(CHAINS_BIN_DIR+"/../BMRB_data")
            
    ftp = ftplib.FTP("ftp.bmrb.wisc.edu")
    ftp.login("anonymous", "")
    ftp.cwd("/pub/bmrb/statistics/chem_shifts/selected/aasel/")
    fnames = ftp.nlst()
    for aa in list(aa3to1_dict.keys()):
        fpattern = re.compile(aa+"_[A-Z0-9]+_hist.txt$")
        file_list=list(filter(fpattern.search, fnames))
        print("Downloading latest "+aa+" chemical shift histograms from BMRB...")
        for filename in file_list:
            #print "Downloading file ",filename
            ftp.retrbinary("RETR " + filename ,open(CHAINS_BIN_DIR+"/../BMRB_data/"+filename, 'wb').write)

def is_H_in_low_occupancy_region(H):
    """
        FUNCTION to assess whether an H resonance lies in the low occupancy regions of BMRB histograms. If true then only Carbon-based aa type predictions will be pursued.
        
        ARGUMENTS:
        H:          H resonance
        
        threshold 0.0001:
        [-2, '-1.23', '-1.17', '-1.15', '-1.11', '-1.11', '-1.07', '-1.07', '6.61', '6.69', '6.75', '6.97', '7.01', '8.01', '8.03', '8.69', '8.73', '8.89', '8.93', '17.01']
        threshold 0.002:
        [-2, '-1.01', '-0.99', '-0.27', '-0.25', '-0.11', '5.85', '5.85', '5.89', '8.01', '8.03', 17.01]
        threshold 0.001:
        [-2, '-1.01', '-0.99', '-0.37', '-0.33', '-0.27', '-0.23', '-0.23', '6.07', '6.09', '6.13', '8.01', '8.03', 17.01]
    """
    #if H <= -1.17 or 6.61 <= H <= 6.69 or 6.75 <= H <= 6.97 or 7.01 <= H:
    if H <= 0.5 or H >= 6.00:
        return True
    else:
        return False


def run_commandline(commandline, logname="log", append=False, return_out=False, error_keywords=[], skip_fail=False,
                    verbose=True):
    """
        FUNCTION to run a single command on the UNIX shell. The worker will only receive an index from network.
    """
    if append:
        fout = open(logname, 'a')
    else:
        fout = open(logname, 'w')
    if verbose:
        print("Running commandline:", commandline)
    return_code = call(commandline, stdout=fout, stderr=fout, shell=True, executable='/bin/bash')

    if (return_code != 0):
        print(ColorPrint("ERROR, THE FOLLOWING COMMAND FAILED TO RUN:", "FAIL"))
        print(commandline)
        print("return_code=", return_code)
        fout.close()
        print("Output:")
        with open(logname, 'r') as f:
            contents = f.readlines()
            for line in contents:
                print(line)
        if not skip_fail:
            raise Exception()
    fout.close()

    if len(error_keywords) > 0:
        with open(logname, 'r') as f:
            contents = f.readlines()
            for line in contents:
                for word in error_keywords:
                    if word in line:
                        ColorPrint("ERROR, THE FOLLOWING COMMAND FAILED TO RUN:", "FAIL")
                        print(commandline)
                        ColorPrint("COMMAND OUTPUT:", "WARNING")
                        for line in contents:
                            print(line)
                        raise Exception()

    if return_out:
        with open(logname, 'r') as f:
            contents = f.readlines()
            return contents

def minmax_scale(X, minmax=[], final_range=[0., 1.]):
    """
    This is the original implementation in Scikit-Learn:
    X_std = (X - X.min(axis=0)) / (X.max(axis=0) - X.min(axis=0))
    X_scaled = X_std * (max - min) + min
    :param X:   an iterable of values to be scaled
    :param minmax:    the min and max values that the variable can take (not necessarily min(X) and max(X)).
    :param final_range: the final value range that you wish X to have
    :return:    X scaled in the "final_range"
    """
    if isinstance(X, (Decimal, float, int)):
        return Decimal('1.0')
    else:
        X = np.array(X)

    if minmax:
        # min, max = float(minmax[0]), float(minmax[1])   # to cope with Decimal
        min, max = minmax[0], minmax[1]   # to cope with Decimal
    else:
        min, max = X.min(), X.max()

    X_std = (X - min) / (max - min)
    if type(min) == Decimal:
        final_range = [Dec(0.), Dec(1.)]  # convert them to Decimals for consistency in calculations
    min_range, max_range = final_range
    X_scaled = X_std * (max_range - min_range) + min_range
    return X_scaled


def chunkIt(seq, num, weights=[]):
    """
    Method to split a list into a specified number of approximately equal sublists.
    :param seq: input iterable
    :param num: number of chunks to split seq
    :param weights: must be a list of integers of the length num
    :return:
    """
    if num == 1:
        return [seq]
    elif num > len(seq):
        return [[s] for s in seq]
    if not weights:
        weights = [1]*num

    assert len(weights) == num, Debuginfo("ERROR: weights must be a list of integers of length num.", fail=True)

    quantum = len(seq) / np.sum(weights) ;
    out = []
    last = 0.0

    for i in range(num):
        out.append(seq[int(last):int(last + quantum*weights[i])])
        last += quantum*weights[i]
    out[i].extend(seq[int(last):])    # append the remaining (if any) to the last chunk

    return out


def multiplier_for_int(val):
    """
    Finds the power of 10 that is necessary to convert all the values of val to integers.
    :param val: float or Decimal object.
    :return:
    """
    decimals = 1  # with 8 it runs out of memory
    alpha = 10 ** decimals
    while alpha * val < 1.0:  # we don't want 0 intensity because that peak will have no density in the 2D-hist!
        decimals += 1
        alpha = 10 ** decimals
    return alpha


def find_nearest_index(array, value):
    idx = (np.abs(array-value)).argmin()
    return idx

def flatten(l):
    """
    FUNCTION to flatten any Iterable object.
    """
    for el in l:
        if isinstance(el, collections.Iterable) and not isinstance(el, str):
            for sub in flatten(el):
                yield sub
        else:
            yield el

def writelist2file(List, fname, header=None):
    """
        FUNCTION to write the contents of a list into a file, each element to a different line.
    """
    with open(fname, 'w') as f:
        if header:
            f.write(header)
        for l in List:
            if type(l) == str:
                if l[-1] != '\n':
                    l += '\n'
                f.write(l)
            elif type(l) in [int, float]:
                f.write(str(l) + "\n")
            else:
                l = [str(i) for i in l] # convert all elements to string
                f.write(" ".join(l) + "\n")

def writelists2file(fname, *lists):
    """
        FUNCTION to write the contents of multiple lists of strings and numbers in filename fname.
        E.g. if l1 = ["a", "b", "c"], l2 = [25.6, 50.5, 100.3], then fname will contain:
        a 25.6
        b 50.5
        c 100.3
    """
    N = len(lists[0])  # asume that all lists have the same size
    with open(fname, 'w') as f:
        for i in range(N):
            line = lists[0][i] + " "
            for l in lists[1:]:
                line += " " + str(l[i])
            f.write(line + "\n")
            f.flush

def replace_multi(text, dict):
    "Method to do multiple replacements in a string."
    for i, j in dict.items():
        text = text.replace(i, j)
    return text

def replace_alt(text, alt_txtlist, new_txt):
    "Method to do replace multiple alternative texts with one specific text."
    for txt in alt_txtlist:
        text = text.replace(txt, new_txt)
    return text

def remove_NH_suffix(text):
    return replace_alt(text, ['N-H', "ND2-HD21", "ND2-HD22", "NE2-HE21", "NE2-HE22", "NE1-HE1", 'NE-HE', 'NX-HX',
                              'NH', "ND2HD21", "ND2HD22", "NE2HE21", "NE2HE22", "NE1HE1", 'NEHE', 'NXHX'], '')

def dash_to_NH(label):
    if '-' in label:    # if a dash already exists in the label
        return label
    elif '?' in label:
        return "?-?"

    label = list(label)
    for i in reversed(list(range(len(label)))):
        if 'H' in label[i]:
            label.insert(i, '-')
            break
    return "".join(label)

def split_AAIG_signature(AAIG_signature):
    """
    Method to split the AAIG signature into AAIG and amide name. If the signature contains '-' then it will
    return the AAIG name and the amide name with '-'. If not then the returned amide will not contain '-' either.
    :param AAIG_signature:
    :return:
    """
    Nnames = ['N', 'ND2', 'ND2', 'NE2', 'NE2', 'NE', 'NE1', '?', 'NX']  # all alternative valid N names
    HNnames = ['H', 'HD21', 'HD22', 'HE21', 'HE22', 'HE', 'HE1', '?', 'HX']  # and the respective HN names
    if '-' in AAIG_signature:
        p = regex.compile(r"^([^\s]+)(\L<Nname>-\L<HNname>)$", regex.X, Nname=Nnames, HNname=HNnames)
    else:
        p = regex.compile(r"^([^\s]+)(\L<Nname>\L<HNname>)$", regex.X, Nname=Nnames, HNname=HNnames)
    m = p.search(AAIG_signature)
    if m:
        AAIG_name = m.group(1)
        NH_name = m.group(2)
    else:
        AAIG_name = '?'
        NH_name = '?'
    return AAIG_name, NH_name

def split_reversed_AAIG_signature(AAIG_signature):
    """
    Method to split the AAIG signature into AAIG and amide name. If the signature contains '-' then it will
    return the AAIG name and the amide name with '-'. If not then the returned amide will not contain '-' either.
    :param AAIG_signature:
    :return:
    """
    Nnames = ['N', 'ND2', 'ND2', 'NE2', 'NE2', 'NE', 'NE1', '?', 'NX']  # all alternative valid N names
    HNnames = ['H', 'HD21', 'HD22', 'HE21', 'HE22', 'HE', 'HE1', '?', 'HX']  # and the respective HN names
    if '-' in AAIG_signature:
        p = regex.compile(r"^([^\s]+)(\L<HNname>-\L<Nname>)$", regex.X, Nname=Nnames, HNname=HNnames)
    else:
        p = regex.compile(r"^([^\s]+)(\L<HNname>\L<Nname>)$", regex.X, Nname=Nnames, HNname=HNnames)
    m = p.search(AAIG_signature)
    if m:
        AAIG_name = m.group(1)
        NH_name = m.group(2)
    else:
        AAIG_name = '?'
        NH_name = '?'
    return AAIG_name, NH_name

def get_NH_name(AAIG_signature):
    """
    Get the amide name from the AAIG signature.
    :param AAIG_signature:
    :return:
    """
    return split_AAIG_signature(AAIG_signature)[1]

def get_aatype_from_AAIG_signature(AAIG_signature):
    AAIG_name = remove_NH_suffix(AAIG_signature)
    return aa1to3_dict[AAIG_name[0]]

def is_sidechain(AAIG_signature):
    """
    Does the given AAIG signature belong to a side chain?
    :param AAIG_signature:
    :return True/False:
    """
    return get_NH_name(AAIG_signature) in ['ND2HD21', 'ND2-HD21', 'ND2HD22', 'ND2-HD22', 'NE2HE21', 'NE2-HE21',
                                           'NE2HE22', 'NE2-HE22', 'NEHE', 'NE-HE', 'NE1HE1', 'NE1-HE1']
def is_valid_signature(AAIG_signature):
    # TODO: maybe is imperfect.
    if re.match("^[A-Z][0-9]+[XNDE12]+[XHDE12]+$", AAIG_signature) or \
            re.match("^[A-Z][0-9]+[XNDE12]+\-[XHDE12]+$", AAIG_signature):
        return True
    else:
        return False

def get_resid(AAIG_signature):
    """
    Get the residue ID from the AAIG signature.
    :param AAIG_signature:
    :return:
    """
    return int(split_AAIG_signature(AAIG_signature)[0][1:])

def reverse_NH_suffix(AAIG_signature):
    """
    For i_AAIG, my naming is "N-HN" while the correct is "HN-N".
    :return:
    """
    if AAIG_signature.split('-')[1][0] == 'N':
        AAIG_name, NH_name = split_reversed_AAIG_signature(AAIG_signature)
    else:
        AAIG_name, NH_name = split_AAIG_signature(AAIG_signature)
    if NH_name == '?':
        return AAIG_signature
    new_NH_name = NH_name.split('-')[1] + "-" + NH_name.split('-')[0]
    return AAIG_name + new_NH_name

def Dec(number, decpoints=4):
    """
    Method to create a read fixed decimal point float number. Beware that every time you do operations between 2 such numbers
    then the decimal points will change (e.g. v1*v2 will have 8 decimal points if v1 & v2 had 4 decimal points).
    :param number: float, int, Decimal or string
    :param decpoints:
    :return:
    """
    try:
        # NOTE: ...if the length of the coefficient after the quantize operation would be greater than precision,
        # NOTE: then an InvalidOperation is signaled. SOLUTION: increase getcontext().prec
        return Decimal(number).quantize(Decimal('1.' + '0'*decpoints))
    except InvalidOperation:
        Debuginfo("FAIL: number " + str(number) + " %s return decimal.InvalidOperation. "
                                                  "Try a smaller 'decpoints' input parameters." % type(number),
                              fail=True)

def Decarange(start, end, step):
    """
    Method like numpy.arange but for Decimals.
    :param start:
    :param end:
    :param step:
    :return:
    """
    start = Dec(start)
    end = Dec(end)
    step = Dec(step)
    current = start
    series = []
    while current < end:
        series.append(current)
        current += step
    return np.array(series, dtype=Decimal)

def save_pickle(fname, *kargs):
    """
        Method to save a variable number of objects into a pickled file.
    """
    with bz2.BZ2File(fname, 'wb') as f:
        pickler = pickle.Pickler(f)
        for item in kargs:
            pickler.dump(item)
        # newData = cPickle.dumps(kargs, 1)
        # f.write(newData)

def load_pickle(fname, pred_num=1000000):
    """
        Method to load a variable number of objects from a pickled file.
        It returns a list with the objects loaded from the pickle file.
    """
    object_list = []
    with bz2.BZ2File(fname, "rb") as pickled_file:
        unpickler = pickle.Unpickler(pickled_file)
        for i in range(pred_num):
            try:
                object_list.append(unpickler.load())
            except EOFError:
                break
        # pickle_object = cPickle.load(pickled_file)
    return object_list


def is_sublist(lst1, lst2):
    try:
        ii = lst2.index(lst1[0])
    except ValueError:
        return False

    if (lst2[ii:ii + len(lst1)] == lst1) and (len(lst2) > len(lst1)):
        return True
    else:
        return False

def print_dict_contents(dictionary):
    """
    Print the contents of a dictionary. Can work with multi-dimensional dictionaries, too.
    :param dictionary: a dict or tree() object
    :return:
    """
    for k,v in list(dictionary.items()):
        sys.stdout.write(str(k) + " --> ")
        try:
            print_dict_contents(v)
        except AttributeError:
            print(v)


def get_residue_from_AAIGsignature(AAIGsignature, absolute_AAIGmatches_alignment_list, absolute_matches_alignment_list):
    """
    Applies only to N-H AAIGs (N-H mapped)
    :param AAIGsignature:
    :param absolute_AAIGmatches_alignment_list:
    :param absolute_matches_alignment_list:
    :return:
    """
    try:
        index = absolute_AAIGmatches_alignment_list.index(AAIGsignature)
        return absolute_matches_alignment_list[index]
    except ValueError:  # if not in the NH table, return the original RIG
        return AAIGsignature


def get_resid_from_residue(residue):

    resid = int(residue[1:])
    return resid