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

"""
ATTENTION: This scripts is OUTDATED and need major update!!!
"""

import sys, re, os, pickle, traceback, shutil, bz2, math
from scoop import futures, shared
import numpy as np
from operator import itemgetter
from collections import OrderedDict
from ete3 import Tree
import ftplib
from argparse import ArgumentParser
from scipy.stats.mstats import zscore
from scipy import stats, sqrt
import collections
import gc
from cluster import HierarchicalClustering
from lib.open_func import ProgressBar
from lib.global_func import *
from lib.trees.chain import *

## Set global variables
HOME_DIR = os.path.dirname(os.path.realpath(__file__))
code3 = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
code1 = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
aa1to3_dict = dict((c1,c3) for c1,c3 in zip(code1,code3))
aa3to1_dict = dict((c3,c1) for c3,c1 in zip(code3,code1))
allowed_aa_atoms_dict = {
"ALA" : ["HA", "HB", "CA", "CB", "N", "H"],
"ARG" : ["HA", "HB2", "HB3", "CA", "CB" "N", "H"],
"ASP" : ["HA", "HB2", "HB3", "CA", "CB", "N", "H"],
"ASN" : ["HA", "HB2", "HB3", "CA", "CB", "N", "H"],
"CYS" : ["HA", "HB2", "HB3", "CA", "CB", "N", "H"],
"GLU" : ["HA", "HB2", "HB3", "CA", "CB", "N", "H"],
"GLN" : ["HA", "HB2", "HB3", "CA", "CB", "N", "H"],
"GLY" : ["HA2", "HA3", "CA", "N", "H"],
"HIS" : ["HA", "HB2", "HB3", "CA", "CB", "N", "H"],
"ILE" : ["HA", "HB", "CA", "CB", "N", "H"],
"LEU" : ["HA", "HB2", "HB3", "CA", "CB", "N", "H"],
"LYS" : ["HA", "HB2", "HB3", "CA", "CB", "N", "H"],
"MET" : ["HA", "HB2", "HB3", "CA", "CB", "N", "H"],
"PHE" : ["HA", "HB2", "HB3", "CA", "CB", "N", "H"],
"PRO" : ["HA", "HB2", "HB3", "CA", "CB", "N"], # Prolines are not detected by the method at position "i" due to lack of HN hydrogen
"SER" : ["HA", "HB2", "HB3", "CA", "CB", "N", "H"],
"THR" : ["HA", "HB", "CA", "CB", "N", "H"],
"TRP" : ["HA", "HB2", "HB3", "CA", "CB", "N", "H"],
"TYR" : ["HA", "HB2", "HB3", "CA", "CB", "N", "H"],
"VAL" : ["HA", "HB", "CA", "CB", "N", "H"]
}
aatype_maxH_C_pairs_dict = {
"ALA" : 2,
"ARG" : 3,
"ASP" : 3,
"ASN" : 3,
"CYS" : 3,
"GLU" : 3,
"GLN" : 3,
"GLY" : 2,
"HIS" : 3,
"ILE" : 2,
"LEU" : 3,
"LYS" : 3,
"MET" : 3,
"PHE" : 3,
"PRO" : 3, # Prolines are not detected by the method at position "i" due to lack of HN hydrogen
"SER" : 3,
"THR" : 2,
"TRP" : 3,
"TYR" : 3,
"VAL" : 2
}


# a dict with keys the amino acid --> carbons and values list of the covalently bonded hydrogens
aa_carbonBondedHydrogensDict_dict = {}
aa_carbonBondedHydrogensDict_dict["ALA"] = {
    "CA" : ["HA"],
    "CB" : ["HB"]
}
#"ARG" : ["HA", "HB2", "HB3", "HG2", "HG3", "HD2", "HD3", "CA", "CB", "CG", "CD"],
aa_carbonBondedHydrogensDict_dict["ARG"] = {
    "CA" : ["HA"],
    "CB" : ["HB2", "HB3"]
}
#"ASP" : ["HA", "HB2", "HB3", "CA", "CB"],
aa_carbonBondedHydrogensDict_dict["ASP"] = {
    "CA" : ["HA"],
    "CB" : ["HB2", "HB3"]
}
#"ASN" : ["HA", "HB2", "HB3", "CA", "CB"],
aa_carbonBondedHydrogensDict_dict["ASN"] = {
    "CA" : ["HA"],
    "CB" : ["HB2", "HB3"]
}
#"CYS" : ["HA", "HB2", "HB3", "CA", "CB"],
aa_carbonBondedHydrogensDict_dict["CYS"] = {
    "CA" : ["HA"],
    "CB" : ["HB2", "HB3"]
}
#"GLU" : ["HA", "HB2", "HB3", "HG2", "HG3", "CA", "CB", "CG"],
aa_carbonBondedHydrogensDict_dict["GLU"] = {
    "CA" : ["HA"],
    "CB" : ["HB2", "HB3"]
}
#"GLN" : ["HA", "HB2", "HB3", "HG2", "HG3", "CA", "CB", "CG"],
aa_carbonBondedHydrogensDict_dict["GLN"] = {
    "CA" : ["HA"],
    "CB" : ["HB2", "HB3"]
}
aa_carbonBondedHydrogensDict_dict["GLY"] = {
    "CA" : ["HA2", "HA3"]
}
#"HIS" : ["HA", "HB2", "HB3", "CA", "CB"],
aa_carbonBondedHydrogensDict_dict["HIS"] = {
    "CA" : ["HA"],
    "CB" : ["HB2", "HB3"]
}
#"ILE" : ["HA", "HB", "HG12", "HG13", "HG2", "HD1", "CA", "CB", "CG1", "CG2", "CD1"],
aa_carbonBondedHydrogensDict_dict["ILE"] = {
    "CA" : ["HA"],
    "CB" : ["HB"]
}
#"LEU" : ["HA", "HB2", "HB3", "HG", "HD1", "HD2", "CA", "CB", "CG", "CD1", "CD2"],
aa_carbonBondedHydrogensDict_dict["LEU"] = {
    "CA" : ["HA"],
    "CB" : ["HB2", "HB3"]
}
#"LYS" : ["HA", "HB2", "HB3", "HG2", "HG3", "HD2", "HD3", "HE2", "HE3", "CA", "CB", "CG", "CD", "CE"],
aa_carbonBondedHydrogensDict_dict["LYS"] = {
    "CA" : ["HA"],
    "CB" : ["HB2", "HB3"]
}
#"MET" : ["HA", "HB2", "HB3", "HG2", "HG3", "CA", "CB", "CG"],
aa_carbonBondedHydrogensDict_dict["MET"] = {
    "CA" : ["HA"],
    "CB" : ["HB2", "HB3"]
}
#"PHE" : ["HA", "HB2", "HB3", "CA", "CB"],
aa_carbonBondedHydrogensDict_dict["PHE"] = {
    "CA" : ["HA"],
    "CB" : ["HB2", "HB3"]
}
#"PRO" : ["HA", "HB2", "HB3", "HG2", "HG3", "HD2", "HD3", "CA", "CB", "CG", "CD"],
aa_carbonBondedHydrogensDict_dict["PRO"] = {
    "CA" : ["HA"],
    "CB" : ["HB2", "HB3"]
}
#"SER" : ["HA", "HB2", "HB3", "CA", "CB"],
aa_carbonBondedHydrogensDict_dict["SER"] = {
    "CA" : ["HA"],
    "CB" : ["HB2", "HB3"]
}
#"THR" : ["HA", "HB", "HG2", "CA", "CB", "CG2"],
aa_carbonBondedHydrogensDict_dict["THR"] = {
    "CA" : ["HA"],
    "CB" : ["HB"]
}
#"TRP" : ["HA", "HB2", "HB3", "CA", "CB"],
aa_carbonBondedHydrogensDict_dict["TRP"] = {
    "CA" : ["HA"],
    "CB" : ["HB2", "HB3"]
}
#"TYR" : ["HA", "HB2", "HB3", "CA", "CB"],
aa_carbonBondedHydrogensDict_dict["TYR"] = {
    "CA" : ["HA"],
    "CB" : ["HB2", "HB3"]
}
#"VAL" : ["HA", "HB", "HG1", "HG2", "CA", "CB", "CG1", "CG2"]
aa_carbonBondedHydrogensDict_dict["VAL"] = {
    "CA" : ["HA"],
    "CB" : ["HB"]
}

# same as aa_carbonBondedHydrogensDict_dict but only with the carbons that have geminal protons
aa_carbonBondedGeminalHydrogensDict_dict = {}
for aa_type in list(aa_carbonBondedHydrogensDict_dict.keys()):
    for carbon, proton_list in list(aa_carbonBondedHydrogensDict_dict[aa_type].items()):
        if len(proton_list) == 2:   # only for carbons with 2 geminal protons !
            try:
                aa_carbonBondedGeminalHydrogensDict_dict[aa_type].append(proton_list)
            except KeyError:
                aa_carbonBondedGeminalHydrogensDict_dict[aa_type] = [proton_list]



## PROBABILITIES OF THE PRESENCE OF A CHEMICAL SHIFT OF A PARTICULAR NUCLEUS TYPE 
Prob_CS = tree()
## ATENTION: "H" probabilies are wrong but are not used up to this version!
Prob_CS["ALA"]["HA"] = 0.627588
Prob_CS["ALA"]["HB"] = 0.485091
Prob_CS["ALA"]["CA"] = 0.720575
Prob_CS["ALA"]["CB"] = 0.681273
Prob_CS["ALA"]["N"] = 0.775434
Prob_CS["ALA"]["H"] = 0.828109
Prob_CS["ARG"]["HA"] = 0.651418
Prob_CS["ARG"]["HB2"] = 0.61372
Prob_CS["ARG"]["HB3"] = 0.560987
Prob_CS["ARG"]["HG2"] = 0.550487
Prob_CS["ARG"]["HG3"] = 0.493159
Prob_CS["ARG"]["HD2"] = 0.540319
Prob_CS["ARG"]["HD3"] = 0.476462
Prob_CS["ARG"]["CA"] = 0.71567
Prob_CS["ARG"]["CB"] = 0.668594
Prob_CS["ARG"]["CG"] = 0.436393
Prob_CS["ARG"]["CD"] = 0.442028
Prob_CS["ARG"]["N"] = 0.766385
Prob_CS["ARG"]["H"] = 0.841117
Prob_CS["ASP"]["HA"] = 0.631589
Prob_CS["ASP"]["HB2"] = 0.617238
Prob_CS["ASP"]["HB3"] = 0.567902
Prob_CS["ASP"]["CA"] = 0.728253
Prob_CS["ASP"]["CB"] = 0.691
Prob_CS["ASP"]["N"] = 0.78801
Prob_CS["ASP"]["H"] = 0.834681
Prob_CS["ASN"]["HA"] = 0.630765
Prob_CS["ASN"]["HB2"] = 0.616623
Prob_CS["ASN"]["HB3"] = 0.573352
Prob_CS["ASN"]["CA"] = 0.700944
Prob_CS["ASN"]["CB"] = 0.665637
Prob_CS["ASN"]["N"] = 0.743466
Prob_CS["ASN"]["H"] = 0.808263
Prob_CS["CYS"]["HA"] = 0.71314
Prob_CS["CYS"]["HB2"] = 0.713565
Prob_CS["CYS"]["HB3"] = 0.672861
Prob_CS["CYS"]["CA"] = 0.565749
Prob_CS["CYS"]["CB"] = 0.537025
Prob_CS["CYS"]["N"] = 0.623102
Prob_CS["CYS"]["H"] = 0.820489
Prob_CS["GLU"]["HA"] = 0.643911
Prob_CS["GLU"]["HB2"] = 0.607481
Prob_CS["GLU"]["HB3"] = 0.550605
Prob_CS["GLU"]["HG2"] = 0.566361
Prob_CS["GLU"]["HG3"] = 0.506864
Prob_CS["GLU"]["CA"] = 0.74446
Prob_CS["GLU"]["CB"] = 0.696098
Prob_CS["GLU"]["CG"] = 0.47692
Prob_CS["GLU"]["N"] = 0.80193
Prob_CS["GLU"]["H"] = 0.844035
Prob_CS["GLN"]["HA"] = 0.646441
Prob_CS["GLN"]["HB2"] = 0.610609
Prob_CS["GLN"]["HB3"] = 0.560065
Prob_CS["GLN"]["HG2"] = 0.573652
Prob_CS["GLN"]["HG3"] = 0.512974
Prob_CS["GLN"]["CA"] = 0.741299
Prob_CS["GLN"]["CB"] = 0.695634
Prob_CS["GLN"]["CG"] = 0.480245
Prob_CS["GLN"]["N"] = 0.78849
Prob_CS["GLN"]["H"] = 0.836707
Prob_CS["GLY"]["HA2"] = 0.621721
Prob_CS["GLY"]["HA3"] = 0.562958
Prob_CS["GLY"]["CA"] = 0.681443
Prob_CS["GLY"]["N"] = 0.719233
Prob_CS["GLY"]["H"] = 0.779058
Prob_CS["HIS"]["HA"] = 0.470238
Prob_CS["HIS"]["HB2"] = 0.454692
Prob_CS["HIS"]["HB3"] = 0.42737
Prob_CS["HIS"]["CA"] = 0.525354
Prob_CS["HIS"]["CB"] = 0.495373
Prob_CS["HIS"]["N"] = 0.545207
Prob_CS["HIS"]["H"] = 0.594266
Prob_CS["ILE"]["HA"] = 0.647827
Prob_CS["ILE"]["HB"] = 0.611565
Prob_CS["ILE"]["HG12"] = 0.582129
Prob_CS["ILE"]["HG13"] = 0.53749
Prob_CS["ILE"]["HG2"] = 0.496295
Prob_CS["ILE"]["HD1"] = 0.506604
Prob_CS["ILE"]["CA"] = 0.749527
Prob_CS["ILE"]["CB"] = 0.704023
Prob_CS["ILE"]["CG1"] = 0.481033
Prob_CS["ILE"]["CG2"] = 0.511658
Prob_CS["ILE"]["CD1"] = 0.520396
Prob_CS["ILE"]["N"] = 0.800226
Prob_CS["ILE"]["H"] = 0.848307
Prob_CS["LEU"]["HA"] = 0.64082
Prob_CS["LEU"]["HB2"] = 0.614473
Prob_CS["LEU"]["HB3"] = 0.565649
Prob_CS["LEU"]["HG"] = 0.527574
Prob_CS["LEU"]["HD1"] = 0.513923
Prob_CS["LEU"]["HD2"] = 0.502025
Prob_CS["LEU"]["CA"] = 0.739859
Prob_CS["LEU"]["CB"] = 0.695436
Prob_CS["LEU"]["CG"] = 0.447941
Prob_CS["LEU"]["CD1"] = 0.497406
Prob_CS["LEU"]["CD2"] = 0.474578
Prob_CS["LEU"]["N"] = 0.790436
Prob_CS["LEU"]["H"] = 0.841981
Prob_CS["LYS"]["HA"] = 0.632944
Prob_CS["LYS"]["HB2"] = 0.592301
Prob_CS["LYS"]["HB3"] = 0.536489
Prob_CS["LYS"]["HG2"] = 0.535574
Prob_CS["LYS"]["HG3"] = 0.475269
Prob_CS["LYS"]["HD2"] = 0.475505
Prob_CS["LYS"]["HD3"] = 0.410943
Prob_CS["LYS"]["HE2"] = 0.467366
Prob_CS["LYS"]["HE3"] = 0.397022
Prob_CS["LYS"]["CA"] = 0.697505
Prob_CS["LYS"]["CB"] = 0.652287
Prob_CS["LYS"]["CG"] = 0.432478
Prob_CS["LYS"]["CD"] = 0.40896
Prob_CS["LYS"]["CE"] = 0.396342
Prob_CS["LYS"]["N"] = 0.753123
Prob_CS["LYS"]["H"] = 0.820389
Prob_CS["MET"]["HA"] = 0.569117
Prob_CS["MET"]["HB2"] = 0.530263
Prob_CS["MET"]["HB3"] = 0.482514
Prob_CS["MET"]["HG2"] = 0.491627
Prob_CS["MET"]["HG3"] = 0.446494
Prob_CS["MET"]["CA"] = 0.66296
Prob_CS["MET"]["CB"] = 0.614469
Prob_CS["MET"]["CG"] = 0.395561
Prob_CS["MET"]["N"] = 0.682714
Prob_CS["MET"]["H"] = 0.718516
Prob_CS["PHE"]["HA"] = 0.628538
Prob_CS["PHE"]["HB2"] = 0.612522
Prob_CS["PHE"]["HB3"] = 0.574238
Prob_CS["PHE"]["CA"] = 0.718078
Prob_CS["PHE"]["CB"] = 0.677572
Prob_CS["PHE"]["N"] = 0.772927
Prob_CS["PHE"]["H"] = 0.827693
Prob_CS["PRO"]["HA"] = 0.587267
Prob_CS["PRO"]["HB2"] = 0.568503
Prob_CS["PRO"]["HB3"] = 0.529848
Prob_CS["PRO"]["HG2"] = 0.51334
Prob_CS["PRO"]["HG3"] = 0.460657
Prob_CS["PRO"]["HD2"] = 0.530976
Prob_CS["PRO"]["HD3"] = 0.493155
Prob_CS["PRO"]["CA"] = 0.651203
Prob_CS["PRO"]["CB"] = 0.612526
Prob_CS["PRO"]["CG"] = 0.432489
Prob_CS["PRO"]["CD"] = 0.435353
Prob_CS["PRO"]["N"] = 0.0214249
Prob_CS["PRO"]["H"] = 0.00036084
Prob_CS["SER"]["HA"] = 0.595354
Prob_CS["SER"]["HB2"] = 0.566237
Prob_CS["SER"]["HB3"] = 0.50678
Prob_CS["SER"]["CA"] = 0.678518
Prob_CS["SER"]["CB"] = 0.63219
Prob_CS["SER"]["N"] = 0.702867
Prob_CS["SER"]["H"] = 0.756586
Prob_CS["THR"]["HA"] = 0.633859
Prob_CS["THR"]["HB"] = 0.584473
Prob_CS["THR"]["HG2"] = 0.484363
Prob_CS["THR"]["CA"] = 0.715949
Prob_CS["THR"]["CB"] = 0.666764
Prob_CS["THR"]["CG2"] = 0.473236
Prob_CS["THR"]["N"] = 0.771457
Prob_CS["THR"]["H"] = 0.826343
Prob_CS["TRP"]["HA"] = 0.629356
Prob_CS["TRP"]["HB2"] = 0.619159
Prob_CS["TRP"]["HB3"] = 0.580376
Prob_CS["TRP"]["CA"] = 0.675044
Prob_CS["TRP"]["CB"] = 0.638108
Prob_CS["TRP"]["N"] = 0.735426
Prob_CS["TRP"]["H"] = 0.82062
Prob_CS["TYR"]["HA"] = 0.64056
Prob_CS["TYR"]["HB2"] = 0.622185
Prob_CS["TYR"]["HB3"] = 0.586077
Prob_CS["TYR"]["CA"] = 0.70362
Prob_CS["TYR"]["CB"] = 0.656173
Prob_CS["TYR"]["N"] = 0.761541
Prob_CS["TYR"]["H"] = 0.835459
Prob_CS["VAL"]["HA"] = 0.65401
Prob_CS["VAL"]["HB"] = 0.612094
Prob_CS["VAL"]["HG1"] = 0.527849
Prob_CS["VAL"]["HG2"] = 0.521265
Prob_CS["VAL"]["CA"] = 0.752713
Prob_CS["VAL"]["CB"] = 0.70053
Prob_CS["VAL"]["CG1"] = 0.516292
Prob_CS["VAL"]["CG2"] = 0.497813
Prob_CS["VAL"]["N"] = 0.809193
Prob_CS["VAL"]["H"] = 0.854516

## AVERAGE HISTOGRAM PROBABILITIES OF A PARTICULAR NUCLEUS TYPE 
aa_nucleus_averageHistProb_mdict = tree()
aa_nucleus_averageHistProb_mdict["CYS"]["HA"] =  0.005
aa_nucleus_averageHistProb_mdict["CYS"]["HB3"] =  0.0048309178744
aa_nucleus_averageHistProb_mdict["CYS"]["HB2"] =  0.00507614213198
aa_nucleus_averageHistProb_mdict["GLN"]["HA"] =  0.00490196078431
aa_nucleus_averageHistProb_mdict["GLN"]["HG2"] =  0.00636942675159
aa_nucleus_averageHistProb_mdict["GLN"]["HG3"] =  0.00641025641026
aa_nucleus_averageHistProb_mdict["GLN"]["HB3"] =  0.00662251655629
aa_nucleus_averageHistProb_mdict["GLN"]["HB2"] =  0.00675675675676
aa_nucleus_averageHistProb_mdict["HIS"]["HA"] =  0.00531914893617
aa_nucleus_averageHistProb_mdict["HIS"]["HB3"] =  0.00555555555556
aa_nucleus_averageHistProb_mdict["HIS"]["HB2"] =  0.0059880239521
aa_nucleus_averageHistProb_mdict["SER"]["HA"] =  0.00471698113208
aa_nucleus_averageHistProb_mdict["SER"]["HB3"] =  0.00584795321637
aa_nucleus_averageHistProb_mdict["SER"]["HB2"] =  0.00613496932515
aa_nucleus_averageHistProb_mdict["VAL"]["HB"] =  0.00552486187845
aa_nucleus_averageHistProb_mdict["VAL"]["HA"] =  0.00471698113208
aa_nucleus_averageHistProb_mdict["VAL"]["HG2"] =  0.0066222783277
aa_nucleus_averageHistProb_mdict["VAL"]["HG1"] =  0.00662251655629
aa_nucleus_averageHistProb_mdict["MET"]["HA"] =  0.00515463917526
aa_nucleus_averageHistProb_mdict["MET"]["HG2"] =  0.00680272108844
aa_nucleus_averageHistProb_mdict["MET"]["HG3"] =  0.00625
aa_nucleus_averageHistProb_mdict["MET"]["HB3"] =  0.00625
aa_nucleus_averageHistProb_mdict["MET"]["HB2"] =  0.00645161290323
aa_nucleus_averageHistProb_mdict["ASN"]["HA"] =  0.00564971751412
aa_nucleus_averageHistProb_mdict["ASN"]["HB3"] =  0.00537634408602
aa_nucleus_averageHistProb_mdict["ASN"]["HB2"] =  0.00546448087432
aa_nucleus_averageHistProb_mdict["PRO"]["HD3"] =  0.00523560209424
aa_nucleus_averageHistProb_mdict["PRO"]["HD2"] =  0.00549450549451
aa_nucleus_averageHistProb_mdict["PRO"]["HG2"] =  0.00537634408602
aa_nucleus_averageHistProb_mdict["PRO"]["HG3"] =  0.00561797752809
aa_nucleus_averageHistProb_mdict["PRO"]["HA"] =  0.00518134715026
aa_nucleus_averageHistProb_mdict["PRO"]["HB3"] =  0.00549450549451
aa_nucleus_averageHistProb_mdict["PRO"]["HB2"] =  0.00561797752809
aa_nucleus_averageHistProb_mdict["LYS"]["HD3"] =  0.00641025641026
aa_nucleus_averageHistProb_mdict["LYS"]["HD2"] =  0.00649350649351
aa_nucleus_averageHistProb_mdict["LYS"]["HE2"] =  0.00724637681159
aa_nucleus_averageHistProb_mdict["LYS"]["HE3"] =  0.00714285714286
aa_nucleus_averageHistProb_mdict["LYS"]["HG2"] =  0.00649350649351
aa_nucleus_averageHistProb_mdict["LYS"]["HG3"] =  0.00632911392405
aa_nucleus_averageHistProb_mdict["LYS"]["HA"] =  0.00471698113208
aa_nucleus_averageHistProb_mdict["LYS"]["HB3"] =  0.00588235294118
aa_nucleus_averageHistProb_mdict["LYS"]["HB2"] =  0.00609756097561
aa_nucleus_averageHistProb_mdict["THR"]["HB"] =  0.00558659217877
aa_nucleus_averageHistProb_mdict["THR"]["HA"] =  0.00529100529101
aa_nucleus_averageHistProb_mdict["THR"]["HG2"] =  0.00724637681159
aa_nucleus_averageHistProb_mdict["PHE"]["HA"] =  0.00485436893204
aa_nucleus_averageHistProb_mdict["PHE"]["HB3"] =  0.00549450549451
aa_nucleus_averageHistProb_mdict["PHE"]["HB2"] =  0.00581395348837
aa_nucleus_averageHistProb_mdict["ALA"]["HB"] =  0.00645161290323
aa_nucleus_averageHistProb_mdict["ALA"]["HA"] =  0.00423728813559
aa_nucleus_averageHistProb_mdict["GLY"]["HA2"] =  0.00425531914894
aa_nucleus_averageHistProb_mdict["GLY"]["HA3"] =  0.00444444444444
aa_nucleus_averageHistProb_mdict["ASP"]["HA"] =  0.00561797752809
aa_nucleus_averageHistProb_mdict["ASP"]["HB3"] =  0.00571428571429
aa_nucleus_averageHistProb_mdict["ASP"]["HB2"] =  0.00602409638554
aa_nucleus_averageHistProb_mdict["ILE"]["HD1"] =  0.00675675675676
aa_nucleus_averageHistProb_mdict["ILE"]["HG2"] =  0.00699300699301
aa_nucleus_averageHistProb_mdict["ILE"]["HG12"] =  0.00564773446129
aa_nucleus_averageHistProb_mdict["ILE"]["HG13"] =  0.00534703677389
aa_nucleus_averageHistProb_mdict["ILE"]["HB"] =  0.00632911392405
aa_nucleus_averageHistProb_mdict["ILE"]["HA"] =  0.00480769230769
aa_nucleus_averageHistProb_mdict["LEU"]["HD2"] =  0.00591715976331
aa_nucleus_averageHistProb_mdict["LEU"]["HD1"] =  0.00591715976331
aa_nucleus_averageHistProb_mdict["LEU"]["HA"] =  0.0049504950495
aa_nucleus_averageHistProb_mdict["LEU"]["HG"] =  0.00581395348837
aa_nucleus_averageHistProb_mdict["LEU"]["HB3"] =  0.00473933649289
aa_nucleus_averageHistProb_mdict["LEU"]["HB2"] =  0.00492610837438
aa_nucleus_averageHistProb_mdict["ARG"]["HD3"] =  0.00666666666667
aa_nucleus_averageHistProb_mdict["ARG"]["HD2"] =  0.00694444444444
aa_nucleus_averageHistProb_mdict["ARG"]["HG2"] =  0.00606060606061
aa_nucleus_averageHistProb_mdict["ARG"]["HG3"] =  0.00613496932515
aa_nucleus_averageHistProb_mdict["ARG"]["HA"] =  0.00469483568075
aa_nucleus_averageHistProb_mdict["ARG"]["HB3"] =  0.00591715976331
aa_nucleus_averageHistProb_mdict["ARG"]["HB2"] =  0.00595238095238
aa_nucleus_averageHistProb_mdict["TRP"]["HA"] =  0.00549450549451
aa_nucleus_averageHistProb_mdict["TRP"]["HB3"] =  0.00684931506849
aa_nucleus_averageHistProb_mdict["TRP"]["HB2"] =  0.00714285714286
aa_nucleus_averageHistProb_mdict["GLU"]["HA"] =  0.00478468899522
aa_nucleus_averageHistProb_mdict["GLU"]["HG2"] =  0.00724637681159
aa_nucleus_averageHistProb_mdict["GLU"]["HG3"] =  0.00751879699248
aa_nucleus_averageHistProb_mdict["GLU"]["HB3"] =  0.00709219858156
aa_nucleus_averageHistProb_mdict["GLU"]["HB2"] =  0.00704225352113
aa_nucleus_averageHistProb_mdict["TYR"]["HA"] =  0.00471698113208
aa_nucleus_averageHistProb_mdict["TYR"]["HB3"] =  0.00546448087432
aa_nucleus_averageHistProb_mdict["TYR"]["HB2"] =  0.00555555555556
aa_nucleus_averageHistProb_mdict["CYS"]["CB"] =  0.00337837837838
aa_nucleus_averageHistProb_mdict["CYS"]["CA"] =  0.00540540540541
aa_nucleus_averageHistProb_mdict["GLN"]["CB"] =  0.00581395348837
aa_nucleus_averageHistProb_mdict["GLN"]["CA"] =  0.00684931506849
aa_nucleus_averageHistProb_mdict["GLN"]["CG"] =  0.0077519379845
aa_nucleus_averageHistProb_mdict["ILE"]["CG1"] =  0.00564971751412
aa_nucleus_averageHistProb_mdict["ILE"]["CB"] =  0.00578034682081
aa_nucleus_averageHistProb_mdict["ILE"]["CA"] =  0.00571428571429
aa_nucleus_averageHistProb_mdict["ILE"]["CG2"] =  0.0062893081761
aa_nucleus_averageHistProb_mdict["ILE"]["CD1"] =  0.00653594771242
aa_nucleus_averageHistProb_mdict["SER"]["CB"] =  0.00602409638554
aa_nucleus_averageHistProb_mdict["SER"]["CA"] =  0.00588235294118
aa_nucleus_averageHistProb_mdict["VAL"]["CG1"] =  0.00763358778626
aa_nucleus_averageHistProb_mdict["VAL"]["CB"] =  0.00595238095238
aa_nucleus_averageHistProb_mdict["VAL"]["CA"] =  0.00574712643678
aa_nucleus_averageHistProb_mdict["VAL"]["CG2"] =  0.00740740740741
aa_nucleus_averageHistProb_mdict["LYS"]["CB"] =  0.0054347826087
aa_nucleus_averageHistProb_mdict["LYS"]["CA"] =  0.00636942675159
aa_nucleus_averageHistProb_mdict["LYS"]["CG"] =  0.00793650793651
aa_nucleus_averageHistProb_mdict["LYS"]["CE"] =  0.00892857142857
aa_nucleus_averageHistProb_mdict["LYS"]["CD"] =  0.00757575757576
aa_nucleus_averageHistProb_mdict["PRO"]["CB"] =  0.00689655172414
aa_nucleus_averageHistProb_mdict["PRO"]["CA"] =  0.00694444444444
aa_nucleus_averageHistProb_mdict["PRO"]["CG"] =  0.00847457627119
aa_nucleus_averageHistProb_mdict["PRO"]["CD"] =  0.00826446280992
aa_nucleus_averageHistProb_mdict["GLY"]["CA"] =  0.00735294117647
aa_nucleus_averageHistProb_mdict["THR"]["CB"] =  0.00505050505051
aa_nucleus_averageHistProb_mdict["THR"]["CA"] =  0.00546448087432
aa_nucleus_averageHistProb_mdict["THR"]["CG2"] =  0.00819672131148
aa_nucleus_averageHistProb_mdict["PHE"]["CB"] =  0.00571428571429
aa_nucleus_averageHistProb_mdict["PHE"]["CA"] =  0.00578034682081
aa_nucleus_averageHistProb_mdict["ALA"]["CB"] =  0.00595238095238
aa_nucleus_averageHistProb_mdict["ALA"]["CA"] =  0.00609756097561
aa_nucleus_averageHistProb_mdict["HIS"]["CB"] =  0.00584795321637
aa_nucleus_averageHistProb_mdict["HIS"]["CA"] =  0.00595238095238
aa_nucleus_averageHistProb_mdict["MET"]["CB"] =  0.00595238095238
aa_nucleus_averageHistProb_mdict["MET"]["CA"] =  0.00671140939597
aa_nucleus_averageHistProb_mdict["MET"]["CG"] =  0.00847457627119
aa_nucleus_averageHistProb_mdict["ASP"]["CB"] =  0.00578034682081
aa_nucleus_averageHistProb_mdict["ASP"]["CA"] =  0.00617283950617
aa_nucleus_averageHistProb_mdict["GLU"]["CB"] =  0.00564971751412
aa_nucleus_averageHistProb_mdict["GLU"]["CA"] =  0.00649350649351
aa_nucleus_averageHistProb_mdict["GLU"]["CG"] =  0.00675675675676
aa_nucleus_averageHistProb_mdict["LEU"]["CB"] =  0.00540540540541
aa_nucleus_averageHistProb_mdict["LEU"]["CA"] =  0.00621118012422
aa_nucleus_averageHistProb_mdict["LEU"]["CG"] =  0.00671140939597
aa_nucleus_averageHistProb_mdict["LEU"]["CD1"] =  0.00684931506849
aa_nucleus_averageHistProb_mdict["LEU"]["CD2"] =  0.00704225352113
aa_nucleus_averageHistProb_mdict["ARG"]["CB"] =  0.00578034682081
aa_nucleus_averageHistProb_mdict["ARG"]["CA"] =  0.00613496932515
aa_nucleus_averageHistProb_mdict["ARG"]["CG"] =  0.00793650793651
aa_nucleus_averageHistProb_mdict["ARG"]["CD"] =  0.00934579439252
aa_nucleus_averageHistProb_mdict["TRP"]["CB"] =  0.00699300699301
aa_nucleus_averageHistProb_mdict["TRP"]["CA"] =  0.00641025641026
aa_nucleus_averageHistProb_mdict["ASN"]["CB"] =  0.00591715976331
aa_nucleus_averageHistProb_mdict["ASN"]["CA"] =  0.00689655172414
aa_nucleus_averageHistProb_mdict["TYR"]["CB"] =  0.0062893081761
aa_nucleus_averageHistProb_mdict["TYR"]["CA"] =  0.00621118012422


# a multidict with keys the amino acid --> carbon-hydrogen pair and values two arrays, the carbon chemical shifts and the respective hydrogen chemical shifts
# this multidict holds all the information in the input BMRB database and will be used to calculate the histograms
aa_CHpair_2Dhist_mdict = tree()
aa_CHpair_2Dhist_mdict["ALA"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_mdict["ALA"]["CB-HB"] = [[], []]

#"ARG" : ["HA", "HB2", "HB3", "HG2", "HG3", "HD2", "HD3", "CA", "CB", "CG", "CD"],
aa_CHpair_2Dhist_mdict["ARG"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_mdict["ARG"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_mdict["ARG"]["CB-HB3"] = [[], []]

#"ASP" : ["HA", "HB2", "HB3", "CA", "CB"],
aa_CHpair_2Dhist_mdict["ASP"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_mdict["ASP"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_mdict["ASP"]["CB-HB3"] = [[], []]

#"ASN" : ["HA", "HB2", "HB3", "CA", "CB"],
aa_CHpair_2Dhist_mdict["ASN"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_mdict["ASN"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_mdict["ASN"]["CB-HB3"] = [[], []]

#"CYS" : ["HA", "HB2", "HB3", "CA", "CB"],
aa_CHpair_2Dhist_mdict["CYS"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_mdict["CYS"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_mdict["CYS"]["CB-HB3"] = [[], []]

#"GLU" : ["HA", "HB2", "HB3", "HG2", "HG3", "CA", "CB", "CG"],
aa_CHpair_2Dhist_mdict["GLU"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_mdict["GLU"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_mdict["GLU"]["CB-HB3"] = [[], []]

#"GLN" : ["HA", "HB2", "HB3", "HG2", "HG3", "CA", "CB", "CG"],
aa_CHpair_2Dhist_mdict["GLN"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_mdict["GLN"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_mdict["GLN"]["CB-HB3"] = [[], []]

#
aa_CHpair_2Dhist_mdict["GLY"]["CA-HA2"] = [[], []]
aa_CHpair_2Dhist_mdict["GLY"]["CA-HA3"] = [[], []]

#"HIS" : ["HA", "HB2", "HB3", "CA", "CB"],
aa_CHpair_2Dhist_mdict["HIS"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_mdict["HIS"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_mdict["HIS"]["CB-HB3"] = [[], []]

#"ILE" : ["HA", "HB", "HG12", "HG13", "HG2", "HD1", "CA", "CB", "CG1", "CG2", "CD1"],
aa_CHpair_2Dhist_mdict["ILE"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_mdict["ILE"]["CB-HB"] = [[], []]

#"LEU" : ["HA", "HB2", "HB3", "HG", "HD1", "HD2", "CA", "CB", "CG", "CD1", "CD2"],
aa_CHpair_2Dhist_mdict["LEU"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_mdict["LEU"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_mdict["LEU"]["CB-HB3"] = [[], []]

#"LYS" : ["HA", "HB2", "HB3", "HG2", "HG3", "HD2", "HD3", "HE2", "HE3", "CA", "CB", "CG", "CD", "CE"],
aa_CHpair_2Dhist_mdict["LYS"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_mdict["LYS"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_mdict["LYS"]["CB-HB3"] = [[], []]

#"MET" : ["HA", "HB2", "HB3", "HG2", "HG3", "CA", "CB", "CG"],
aa_CHpair_2Dhist_mdict["MET"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_mdict["MET"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_mdict["MET"]["CB-HB3"] = [[], []]

#"PHE" : ["HA", "HB2", "HB3", "CA", "CB"],
aa_CHpair_2Dhist_mdict["PHE"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_mdict["PHE"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_mdict["PHE"]["CB-HB3"] = [[], []]

#"PRO" : ["HA", "HB2", "HB3", "HG2", "HG3", "HD2", "HD3", "CA", "CB", "CG", "CD"],
aa_CHpair_2Dhist_mdict["PRO"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_mdict["PRO"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_mdict["PRO"]["CB-HB3"] = [[], []]

#"SER" : ["HA", "HB2", "HB3", "CA", "CB"],
aa_CHpair_2Dhist_mdict["SER"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_mdict["SER"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_mdict["SER"]["CB-HB3"] = [[], []]

#"THR" : ["HA", "HB", "HG2", "CA", "CB", "CG2"],
aa_CHpair_2Dhist_mdict["THR"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_mdict["THR"]["CB-HB"] = [[], []]

#"TRP" : ["HA", "HB2", "HB3", "CA", "CB"],
aa_CHpair_2Dhist_mdict["TRP"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_mdict["TRP"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_mdict["TRP"]["CB-HB3"] = [[], []]

#"TYR" : ["HA", "HB2", "HB3", "CA", "CB"],
aa_CHpair_2Dhist_mdict["TYR"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_mdict["TYR"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_mdict["TYR"]["CB-HB3"] = [[], []]

#"VAL" : ["HA", "HB", "HG1", "HG2", "CA", "CB", "CG1", "CG2"]
aa_CHpair_2Dhist_mdict["VAL"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_mdict["VAL"]["CB-HB"] = [[], []]


def approx_equal(x, y, tolerance=0.001):
    return abs(x-y) <= 0.5 * tolerance * (x + y)

def download_CS_histograms():
    """
    FUNCTION to download the most recent chemical shift histograms from BMRB 
    """
    
    global aa3to1_dict, aa1to3_dict
    
    if not os.path.exists(HOME_DIR+"/BMRB_data"):
            os.makedirs(HOME_DIR+"/BMRB_data")
            
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
            ftp.retrbinary("RETR " + filename ,open(HOME_DIR+"/BMRB_data/"+filename, 'wb').write)

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


## Parse command line arguments
def cmdlineparse():
    parser = ArgumentParser(description="command line arguments", epilog="EXAMPLE: 4D_assignment_tree-based.py -hsqc TDhsqc.list -tocsy TDtocsy.list -noesy TDnoesy.list -tolH 0.01 -tolC 0.1 -rtolH 0.01 -rtolN 0.1 -mcutoff 1.0 -acutoff 1.0 -wH 0.5 -wC 1")
    parser.add_argument("-fasta", dest="template_sequence_file", required=True, help="teplate sequence file in fasta format", metavar="<template sequence file>")
    parser.add_argument("-hsqc", dest="HSQC_FILE", required=True, help="2D N-H HSQC root spectrum", metavar="<2D N-H HSQC input file>")
    parser.add_argument("-tocsy", dest="TOCSY_fname", required=True, help="4D TOCSY (HCTOCSYNH) file", metavar="<4D TOCSY input file>")
    parser.add_argument("-noesy", dest="NOESY_FILE", required=True, help="4D NOESY (HCNOENH) file", metavar="<4D NOESY input file>")
    parser.add_argument("-tolH", dest="tolH", required=False, type=float, default=0.03, help="tolerance for the proton resonance when matching NOESY peaks to TOCSY peaks", metavar="<proton tolerance>")
    parser.add_argument("-tolC", dest="tolC", required=False, type=float, default=0.3, help="tolerance for the carbon resonance when matching NOESY peaks to TOCSY peaks", metavar="<carbon tolerance>")
    #parser.add_argument("-rtolH", dest="rtolH", required=False, type=float, default=0.02, help="tolerance for the proton resonance when matching peaks between the root spectum and TOCSY/NOESY", metavar="<proton tolerance>")
    #parser.add_argument("-rtolN", dest="rtolN", required=False, type=float, default=0.2, help="tolerance for the nitrogen resonance when matching peaks between the root spectum and TOCSY/NOESY", metavar="<carbon tolerance>")
    parser.add_argument("-update", dest="DOWNLOAD_CS_HISTOGRAMS", required=False, action='store_true', help="download the latest chemical shift histgrams from BMRB")
    parser.add_argument("-wH", dest="H_weight", required=False, type=float, default=1.0, help="weight in (0.0-1.0] to apply on aa type prediction from aliphatic H resonances", metavar="<H weight>")
    parser.add_argument("-wC", dest="C_weight", required=False, type=float, default=1.0, help="weight in (0.0, 1.0] to apply on aa type prediction from aliphatic C resonances", metavar="<C weight>")
    #parser.add_argument("-examples", dest='EXAMPLES', action='store_true', required=False, help="Usage and Examples")
    parser.add_argument("-mcutoff", dest="RESONANCE_MATCH_CUTOFF", required=False, type=float, default=0.8, help="number 0.0-1.0 saying how much of TOCSY resonances should match in the NOESY in order to consider it a possible match", metavar="<resonance match cutoff>")
    parser.add_argument("-zmcutoff", dest="ZSCORE_RESONANCE_MATCH_CUTOFF", required=False, type=float, default=0.0, help="a real number specifying the lower Z-score for a resonance match to be retained for chain building", metavar="<Z-score resonance match cutoff>")
    parser.add_argument("-acutoff", dest="ASSIGNMENT_CUTOFF", required=False, type=float, default=None, help="number 0.0-1.0 saying how much of TOCSY resonances should match with a particular aa type in order to be considered, e.g. \
                        if TOCSY resonances are 5, 0.8 will mean that at least 4 resonances must match", metavar="<aa type assignment cutoff>")
    parser.add_argument("-zacutoff", dest="ZSCORE_ASSIGNMENT_CUTOFF", required=False, type=float, default=-1.0, help="a real number specifying the lower Z-score for an aa type prediction to be considered as valid", metavar="<Z-score aa type assignment cutoff>")
    #parser.add_argument("-moreaa", dest="ALLOW_SINGLE_CH_PAIRS", required=False, action='store_true', help="allow aa type predictions from single C-H resonance pairs. Because it increases memory consumption and lowers findelity, it should be used in the last rounds of assignment.")
    parser.add_argument("-maxlen", dest="MAX_PEPTIDE_LENGTH", required=False, type=int, default=8, help="the maximum peptide length (high values require more memery; default: 8)", metavar="<max peptide length>")
    parser.add_argument("-minlen", dest="MIN_PEPTIDE_LENGTH", required=False, type=int, default=3, help="the minimum peptide length (low values increase noise; default: 2)", metavar="<min peptide length>")
    parser.add_argument("-resoncut", dest="JUST_CARBON_MATCH_CUTOFF", required=False, type=int, default=3,
                        help="the minimum number of C-H resonance pairs for a TOCSY index group to start predicting aa types from both C-H or C only resonances",
                        metavar="<Carbon prediction cutoff>")
    
    #parser.add_argument("-confile", dest="CONNECTIVITIES_FILE", required=False, help="connectivities file; if specified connectivity calculation from input files will be skipped", metavar="<connectivities file>")
    parser.add_argument("-poolconfile", dest="POOL_CONNECTIVITIES_FILE", required=False, default=None, help="connectivities pool file; necessary if -confile specified in order to calculate correct probabilities", metavar="<connectivities pool file>")
    parser.add_argument("-allconfile", dest="COMPLETE_CONNECTIVITIES_FILE", required=False, default=None, help="all connectivities file; necessary if -confile specified in order to calculate correct probabilities", metavar="<all connectivities file>")
    #parser.add_argument("-aafile", dest="AA_TYPES_FILE", required=False, help="amino acid assignment file; if specified amino acid assignment calculation from input files will be skipped", metavar="<amino acid assignment file>")
    parser.add_argument("-poolaafile", dest="POOL_AA_TYPES_FILE", required=False, default=None, help="pool amino acid assignment file; if specified amino acid assignment calculation from input files will be skipped", metavar="<pool amino acid assignment file>")
    parser.add_argument("-allaafile", dest="COMPLETE_AA_TYPES_FILE", required=False, default=None, help="all amino acid assignment file; necessary if -aafile specified in order to calculate correct probabilities", metavar="<all amino acid assignment file>")
    parser.add_argument("-chainfile", dest="NON_REDUNDANT_CHAINS_FILE", required=False, default=None, help="non-redundant chains file; if specified, calculation of chains from connectivities will be omitted", metavar="<non-redundant chains file>")
    parser.add_argument("-zmin", dest="MIN_NUM_OF_PREDICTIONS", required=False, type=int, default=4, help="minimum number of aa type predictions required to apply the Z-score cutoff. The fewer the predictions \
                        the more inaccurate is Z-score (default: 4)", metavar="<minimum number of aa types prediction for Z-score filtering>")
    parser.add_argument("-transform", dest="TRANSFORM_TYPE", required=False, type=str, default='None', help="type of mathematical transform to apply on the probabilities P[TAAIG(i)|aatype(i-1)]. Allowed values are: \"None\", \"log\", \"log10\", \"boxcox_pearson\", \"boxcox_mle\"",
                        metavar="<type of mathematical transform>")
    parser.add_argument("-resume", dest="RESUME", required=False, action='store_true', default=False,
                        help="resume previous run by loading all peptide sequences saved in tmp_peptide_folder/. By default tmp_peptide_folder/ will be cleaned and new peptide sequences will be writen inside it.")
    parser.add_argument("-probprod", dest="PROB_PRODUCT", required=False, action='store_true', default=False,
                        help="select the best C-H type assignment combination based on the product of probabilities of the individual C-H assignments")
    parser.add_argument("-probmodel", dest="PROBABILITY_MODEL", required=False, type=int, default=2,
                        help="If '1' the probability of each peak will be given by [wH*1Dhist(H)+wC*1Dhist(C)]/(wH+wC). If '2' then by 1Dhist(H)*1Dhist(C)")
    parser.add_argument("-cgrpprob", dest="CONSENSUS_CGROUP_PROB_MODE", required=False, default=2, type=int,
                        help="""The way to calculate the total score of a set of chemical shift assignment (default: 2). Can be:
                        0: just multiply all the probabilities of the individual peaks
                        The following values control how to calculate the consensus probability of each C-group. The total score will be the
                        product of this consensus C-group probabilities.
                        1: average;
                        2: sqrt(prob1*prob2)    ; geometric mean
                        3: (prob2-prob1)/(log(prob2)-log(prob1))    ; logarithmic mean
                        4: (prob1*prob2)/((prob1+prob2)/2.) ;composite average 
                        5: prob1*prob2    ; product of probabilities
                        When the C-group contains only one peak, the respective probability will be prob1 in all cases.
                            """, metavar="<way to calculate the cs assignment score>")
    parser.add_argument("-2dhist", dest="USE_2D_HISTOGRAMS", required=False, action='store_true', default=False,
                        help="use 2D BMRB histograms for aa type prediction")
    parser.add_argument("-log", dest="LOG", required=False, action='store_true', default=False,
                        help="convert aa type prediction probabilities to logarithmic scale and then calculate Z-scores, if the min(probability)/max(probability) > 1000")
    parser.add_argument("-delpred", dest="DELETE_AA_TYPE_PREDICTIONS", required=False, action='store_true', default=False,
                        help="delete aa type predictions with probabilities which are 1000, 10000 or 100000 times lower than the highest, if the highest is >10e-10, >10e-20, <=10e-20, respectively.")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s version 1.0')
    args=parser.parse_args()
    return args


##
##  CHECK INPUT FILES FORMAT FOR CORRECTNESS
##

def check_input_files_format():
    global args
    
    # Check NOESY file
    query_contents = []
    HAS_INTENSITIES = False
    with open(args.NOESY_FILE, 'r') as f:
        contents = f.readlines()
        for line in contents: # CLEANING THE SPECTRUM FROM IRRELEVANT LINES
            word_list = line.split()
            try:
                if not word_list[0] == '?-?-?-?':
                    continue
                float(word_list[1]); float(word_list[2]); float(word_list[3]); float(word_list[4]); # checking if the line contains numbers (resonances)
                if len(word_list) == 6:
                    float(word_list[5])
                    HAS_INTENSITIES = True
                query_contents.append(line)
            except (IndexError, ValueError):
                # print "WARNING: Discarding", source_spectrum, "line:", line
                print("DEBUG check_input_files_format: this NOESY line will be discarded:", line)
                continue
    # Now check if all the lines contain intensity
    if HAS_INTENSITIES:
        for word_list in query_contents:
            if len(word_list) < 6:
                print("ERROR: Please check your input NOESY file. Not all lines have intensity!!!")
                sys.exit(1)


##
##  FIND NON-OVERLAPPING RESONANCE GROUPS IN THE HSQC SPECTRUM
##

def find_nonoverlapping_groups_in_root(remaining_root_contents, rtolH, rtolN, keep_closest=False):
    """
        RETURNS:
        nonoverlapping_groups_lines_and_tolerances:    list with the line of each non-overlapping group of the root spectrum along with the rtolH and rtolN
        remaining_groups_lines:                        list with the lines of the Root spectrum that contain overlapping groups with the current torelances
    """
    
    nonoverlapping_groups_lines_and_tolerances = []
    remaining_groups_lines = []
    
    if keep_closest == False:
        for root1_line in remaining_root_contents:
            found_overlap = False
            #print "DEBUG: root1_line=", root1_line
            root1_words_list = root1_line.split()
            try:
                if root1_words_list[0][-3:] != "N-H":   # if it's not a backbone amide, skip it
                    continue
                root1_name = root1_words_list[0]
                root1_N_resonance = float(root1_words_list[1])
                root1_HN_resonance = float(root1_words_list[2])
                for root2_line in remaining_root_contents:
                    if root1_line == root2_line:    # in case we are looking at the same line
                        continue
                    #print "DEBUG: root2_line=",root2_line
                    root2_words_list = root2_line.split()
                    try:
                        root2_name = root2_words_list[0]
                        root2_N_resonance = float(root2_words_list[1])
                        root2_HN_resonance = float(root2_words_list[2])
                        if ( (root2_HN_resonance -rtolH) <= root1_HN_resonance <= (root2_HN_resonance +rtolH) ) and ( (root2_N_resonance -rtolN) <= root1_N_resonance <= (root2_N_resonance +rtolN) ):
                            remaining_groups_lines.append(root1_line)   # 
                            found_overlap = True
                            break
                    except (ValueError, IndexError):
                        #print "WARNING: the 4th and 5th elements of the following root file line are not numbers:"
                        #print "Root file line:", root_line
                        continue
                if found_overlap == False:  # if no noverlapping line was found for this line then save it
                    nonoverlapping_groups_lines_and_tolerances.append((root1_line, rtolH, rtolN))
            except (ValueError, IndexError):
                #print "WARNING: the 2nd and 3rd elements of the following root file line are not numbers:"
                #print "Root file line:", root_line
                continue
    
    #if keep_closest == True:
    #    for root1_line in remaining_groups_lines:
    #        matching_lines = []
    #        found_overlap = False
    #        #print "DEBUG: root1_line=", root1_line
    #        root1_words_list = root1_line.split()
    #        try:
    #            if root1_words_list[0][-3:] != "N-H":   # if it's not a backbone amide, skip it
    #                continue
    #            root1_name = root1_words_list[0]
    #            root1_N_resonance = float(root1_words_list[1])
    #            root1_HN_resonance = float(root1_words_list[2])
    #            for root2_line in remaining_root_contents:
    #                if root1_line == root2_line:    # in case we are looking at the same line
    #                    continue
    #                #print "DEBUG: root2_line=",root2_line
    #                root2_words_list = root2_line.split()
    #                try:
    #                    root2_name = root2_words_list[0]
    #                    root2_N_resonance = float(root2_words_list[1])
    #                    root2_HN_resonance = float(root2_words_list[2])
    #                    if ( (root2_HN_resonance -rtolH) <= root1_HN_resonance <= (root2_HN_resonance +rtolH) ) and ( (root2_N_resonance -rtolN) <= root1_N_resonance <= (root2_N_resonance +rtolN) ):
    #                        found_overlap = True
    #                        matching_lines.append(root1_line)
    #                except (ValueError, IndexError):
    #                    #print "WARNING: the 4th and 5th elements of the following root file line are not numbers:"
    #                    #print "Root file line:", root_line
    #                    continue
    #            if found_overlap == False:
    #                delta_previous_root = np.sqrt((previous_root_HN_resonance - query_HN_resonance)**2 + ((previous_root_N_resonance - query_N_resonance)/6)**2)
    #                delta_root = np.sqrt((root_HN_resonance - query_HN_resonance)**2 + ((root_N_resonance - query_N_resonance)/6)**2)
    #                if delta_previous_root > delta_root: 
    #                    nonoverlapping_groups_lines_and_tolerances.append((root1_line, rtolH, rtolN))
    #        except (ValueError, IndexError):
    #            #print "WARNING: the 2nd and 3rd elements of the following root file line are not numbers:"
    #            #print "Root file line:", root_line
    #            continue
    #
    #print "DEBUG: len(nonoverlapping_groups_lines_and_tolerances)=", len(nonoverlapping_groups_lines_and_tolerances)
    print("DEBUG: nonoverlapping_groups_lines_and_tolerances =", nonoverlapping_groups_lines_and_tolerances)
    print("DEBUG: remaining_groups_lines = ", remaining_groups_lines)
    return nonoverlapping_groups_lines_and_tolerances


def copy_aaindices_from_root_spectrum_2(root_contents, query_fname, spectrum_type, TOCSY_sidechain_resonances_list=None):
    
    """
    FUNCTION that finds the closest TOCSY peaks to each HSQC spectrum group (namely, each HSQC line). If the lower distance of a TOCSY peak from any HSQC group is
    greater than (0.04, 0.4) then it will be considered as garbage and will be discarded from TOCSYnum.list file.
    """
    tolH = 0.04
    tolN = 0.4
    with open(query_fname, 'r') as f:
        tmp_query_contents=f.readlines()    # contents of original query_fname (4D TOCSY or 4D NOESY) in 5 column format (name H C N HN)
    query_contents=[]
    for line in tmp_query_contents: # CLEANING THE SPECTRUM FROM IRRELEVANT LINES
        word_list = line.split()
        try:
            float(word_list[1]); float(word_list[2]); float(word_list[3]); float(word_list[4]); # checking if the line contains numbers (resonances)
            if len(word_list) == 6:    # if there is intensity column check if it is a number
                float(word_list[5])
                word_list[5] = str(abs(float(word_list[5])))    # convert all intensities to positive numbers
            query_contents.append(" ".join(word_list)+"\n")
        except (IndexError, ValueError):
            print("WARNING: Discarding", spectrum_type, "line:", line)
            print("DEBUG: copy_aaindices_from_root_spectrum_2 point1 word_list=", word_list)
    
    # remove duplicate lines from query_fname (4D TOCSY or 4D NOESY)
    lines2remove_set = set()
    for qline1 in query_contents:
        try:
            counter = 0
            query1_words_list = qline1.split()
            q1_w1=float(query1_words_list[1])
            q1_w2=float(query1_words_list[2])
            q1_w3=float(query1_words_list[3])
            q1_w4=float(query1_words_list[4])
            for qline2 in query_contents:
                query2_words_list = qline2.split()
                q2_w1=float(query2_words_list[1])
                q2_w2=float(query2_words_list[2])
                q2_w3=float(query2_words_list[3])
                q2_w4=float(query2_words_list[4])
                if approx_equal(q1_w1, q2_w1) and approx_equal(q1_w2, q2_w2) and approx_equal(q1_w3, q2_w3) and approx_equal(q1_w4, q2_w4):
                    counter += 1
                    if counter > 1:
                        lines2remove_set.add(qline2)
                        #print "DEBUG: will remove line",qline2
        except (ValueError, IndexError):
            #print "WARNING: the 2nd and 3rd elements of the following root file line are not numbers:"
            #print "Root file line:", root_line
            continue
    # now remove the duplicate lines
    for qline in lines2remove_set:
        query_contents.remove(qline)
    
    aaindex_Hresonance_Cresonance_tuple_list = []
    for triplet in root_contents:
        root_line = triplet[0]
        #print "DEBUG: root_line=", root_line
        root_words_list = root_line.split()
        try:
            root_name = root_words_list[0]
            root_N_resonance = float(root_words_list[1])
            root_HN_resonance = float(root_words_list[2])
            for q_index in range(0,len(query_contents)):
                query_line = query_contents[q_index]
                #print "DEBUG: query_line=",query_line
                try:
                    query_words_list = query_line.split()
                    query_N_resonance = float(query_words_list[3])
                    query_HN_resonance = float(query_words_list[4])
                    # First remove the garbage that are far away from any peak
                    if ( (query_HN_resonance -tolH) <= root_HN_resonance <= (query_HN_resonance +tolH) ) and ( (query_N_resonance -tolN) <= root_N_resonance <= (query_N_resonance +tolN) ):
                        #print "DEBUG: ?-?-"+root_name+"\t"+query_words_list[1]+"\t"+query_words_list[2]+"\t"+query_words_list[3]+"\t"+query_words_list[4]
                        aa_index = re.sub("N-H", "", root_name)
                        if re.search("^\s*\?-\?-\?-\?\s+", query_contents[q_index]):   # if the line does not have an AAIG group already, assign the currect AAIG group to it
                            query_contents[q_index] = re.sub("\?-\?-\?-\?", aa_index, query_contents[q_index]).replace("\\", "")
                            Hresonance = query_contents[q_index].split()[1]
                            Cresonance = query_contents[q_index].split()[2]
                            aaindex_Hresonance_Cresonance_tuple_list.append((aa_index, Hresonance, Cresonance))
                        else: # if it has already been assigned an AAIG group, check if the currect AAIG group has closer N & HN resonances
                            #print "DEBUG: Already assigned AAIG group. Line: ", query_contents[q_index]
                            mo = re.search('^\s*([A-Za-z0-9]+)\s+', query_contents[q_index])
                            if mo:
                                previous_aa_index = mo.group(1)     # the AAIG group that is already assigned to this peak from a previous round
                                for tmp_triplet in root_contents:   # iterate over all AAIG groups until you find previous_aa_index, in order to retrieve its N & HN resonances
                                    tmp_root_line = tmp_triplet[0]
                                    tmp_rtolH = tmp_triplet[1]
                                    tmp_rtolN = tmp_triplet[2]
                                    tmp_root_words_list = tmp_root_line.split()
                                    try:
                                        tmp_root_name = tmp_root_words_list[0]
                                        tmp_aa_index = re.sub("N-H", "", tmp_root_name)
                                        if tmp_aa_index == previous_aa_index:   # we found the N & HN resonances of the assigned AAIG group
                                            #print "DEBUG: we found the N & HN resonances of the assigned AAIG group:", tmp_root_words_list
                                            previous_root_N_resonance = float(tmp_root_words_list[1])
                                            previous_root_HN_resonance = float(tmp_root_words_list[2])
                                            # if the current AAIG group has N & HN resonances closer to the query line than the previously assigned N & HN resonances, change the AAIG group
                                            #print "DEBUG: query_N_resonance=", query_N_resonance, "previous_root_N_resonance=", previous_root_N_resonance, "root_N_resonance=", root_N_resonance
                                            #print "DEBUG: query_HN_resonance=", query_HN_resonance, "previous_root_HN_resonance=", previous_root_HN_resonance, "root_HN_resonance=", root_HN_resonance
                                            delta_previous_root = np.sqrt((previous_root_HN_resonance - query_HN_resonance)**2 + ((previous_root_N_resonance - query_N_resonance)/6)**2)
                                            delta_root = np.sqrt((root_HN_resonance - query_HN_resonance)**2 + ((root_N_resonance - query_N_resonance)/6)**2)
                                            if delta_previous_root > delta_root:
                                                #print "DEBUG: substituting previous_aa_index", previous_aa_index, " for ", aa_index, " in ", query_contents[q_index]
                                                if len(query_words_list) == 6:
                                                    query_contents[q_index] = " ".join([aa_index, query_words_list[1], query_words_list[2], query_words_list[3], query_words_list[4], query_words_list[5], "\n"])
                                                else:
                                                    query_contents[q_index] = " ".join([aa_index, query_words_list[1], query_words_list[2], query_words_list[3], query_words_list[4], "\n"])
                                                #query_contents[q_index] = re.sub(previous_aa_index, aa_index, query_contents[q_index])
                                            # otherwise keep the previous AAIG group
                                            # ATTENTION: THIS IS NOT ENTIRELY CORRECT, BECAUSE THE N RESONANCE MAY BE CLOSER BUT THE H NOT!
                                            break
                                    except (ValueError, IndexError):
                                        #print "WARNING: the 2nd and 3rd elements of the following root file line are not numbers:"
                                        #print "Root file line:", root_line
                                        continue
                            else:
                                print("ERROR: wrong line (modified or not) in ", query_fname)
                                print(query_contents[q_index])
                                sys.exit(1)
                            
                except (ValueError, IndexError):
                    #print "WARNING: the 4th and 5th elements of the following query file line are not numbers:"
                    #print "Query file line:", query_line
                    continue
        except (ValueError, IndexError):
            #print "WARNING: the 2nd and 3rd elements of the following root file line are not numbers:"
            #print "Root file line:", root_line
            continue
    
    print("DEBUG: ", spectrum_type," contents: ")
    print("".join(query_contents))
    with open(query_fname.replace(".list", "num.list"), 'w') as f:
        for line in query_contents:
            try:
                print("DEBUG: saving line:", line)
                word_list = line.split()
                if len(word_list) == 6:
                    appendix = "\t"+word_list[5]
                else:
                    appendix = ""
                if spectrum_type == "TOCSY":
                    if not "?" in word_list[0]:
                        f.write("?-?-"+word_list[0]+"N-H\t\t"+word_list[1]+"\t"+word_list[2]+"\t"+word_list[3]+"\t"+word_list[4]+appendix+"\n")
                    else:
                        f.write("?-?-?-?\t\t"+word_list[1]+"\t"+word_list[2]+"\t"+word_list[3]+"\t"+word_list[4]+appendix+"\n")
                elif spectrum_type == "NOESY":
                    if not "?" in word_list[0]:
                        f.write("?-?-"+word_list[0]+"N-H\t\t"+word_list[1]+"\t"+word_list[2]+"\t"+word_list[3]+"\t"+word_list[4]+appendix+"\n")
                    else:
                        f.write("?-?-?-?\t\t"+word_list[1]+"\t"+word_list[2]+"\t"+word_list[3]+"\t"+word_list[4]+appendix+"\n")
            except IndexError:
                print("WARNING: Discarding", spectrum_type, "line:", line)
                #print "DEBUG: line=", line
                print("DEBUG: copy_aaindices_from_root_spectrum_2 point 2 word_list=", word_list)
                continue
    
    #print "DEBUG: aaindex_Hresonance_Cresonance_tuple_list = ",aaindex_Hresonance_Cresonance_tuple_list
    query_lines = []    # list of the query_fname lines in 3 column format (namely without the root spectrum N & HN resonances) 
    with open(query_fname.replace(".list", "_aa-H-C.list"), 'w') as f:
        for triplet in aaindex_Hresonance_Cresonance_tuple_list:
            line = "\t?-?-?\t"+triplet[0]+"\t"+triplet[1]+"\t"+triplet[2]+"\n"
            f.write(line)
            query_lines.append(line)
    
    if spectrum_type == "TOCSY":
        # save unassigned TOCSY lines that correspond to sidechains. The resonances must be later removed from NOESY
        TOCSY_sidechain_resonances_list = []
        for qline in query_contents:
            try:
                print("DEBUG: qline.split()[0]:",qline.split()[0])
                if qline.split()[0] == "?-?-?-?":
                    TOCSY_sidechain_resonances_list.append(qline)
            except IndexError:
                continue
        for TOCSY_line in TOCSY_sidechain_resonances_list:  # TEMPORARILY REMOVE ALL "?-?-?-?" LINES (NOT ALL OF THEM ARE SIDE CHAIN RESONANCES FROM TOCSY FILE CONTENTS)
            query_contents.remove(TOCSY_line)   # now you won't see TOCSY index groups like "?-?-?-?" in the connectivities and aa type prediction files
        return query_contents, TOCSY_sidechain_resonances_list
    elif spectrum_type == "NOESY":
        lines2remove_set = set()
        for qline in query_contents:    # FIND ALL LINES WITH "?-?-?-?"
            try:
                #print "DEBUG: qline.split()[0]:",qline.split()[0]
                if qline.split()[0] == "?-?-?-?":
                    lines2remove_set.add(qline)
            except IndexError:
                continue
        # remove the sidechain resonances from NOESY
        #for TOCSY_line in TOCSY_sidechain_resonances_list:
        for qline in lines2remove_set:  # REMOVE ALL LINES WITH "?-?-?-?"
            #print "DEBUG: removing side chain line from NOESY: ",TOCSY_line
            try:
                #query_contents.remove(TOCSY_line)
                query_contents.remove(qline)
            except ValueError:
                print("DEBUG: side chain line not present in NOESY!")
                continue
        return query_contents
    
    pass


def get_possible_connectivities(TOCSY_contents, NOESY_lines):
    """
    FUNCTION that find all the NOESY aa indices that have resonances w2,w3 that match at least one w2,w3 TOCSY resonance pair.
    
    ARGUMENTS:
    TOCSY_contents:     list of lists, sorted by the 2nd element (TOCSY AAIG); these are the contents of TOCSY file in a convenient form
    NOESY_lines:     list of lists; these are the contents of NOESY file but not sorted
    
    RETURNS:
    TOCSYaaindex_NOESYaaindex_OccupancyNumOfResonancesList_dict:    a multidimensional dictionary with structure:
    TOCSY AAIG => NOESY AAIG => [occupancy, total resonance number] , namely
    TOCSY AAIG => NOESY AAIG => [number of matched TOCSY w2,w3 resonances in NOESY for that NOESY AAIG, total number of TOCSY w2,w3 resonances for that TOCSY AAIG]
    """
    global args
    # dictionary with keys the TOCSY aa indices, and values
    TOCSYaaindex_NOESYaaindex_OccupancyNumOfResonancesList_dict = tree()
    previous_TOCSY_aaindex = None
    matchingNOESYaaindex_occupancy_dict = {}
    Num_of_TOCSY_resonances = 0 # save here the number of TOCSY resonances for a particular AAIG
    for TOCSY_words_list in TOCSY_contents:
        try:
            # discard 1st column
            #TOCSY_aaindex=int(float(TOCSY_words_list[1]))    # residue i
            TOCSY_aaindex=str(TOCSY_words_list[0])    # residue i
            #print "DEBUG: TOCSY_aaindex=", TOCSY_aaindex
            if previous_TOCSY_aaindex != None and TOCSY_aaindex != previous_TOCSY_aaindex: # if we look at a different AAIG in TOCSY, print the matches and occupancies
                                                                                           # of the previous TOCSY AAIG and clear matchingNOESYaaindex_occupancy_dict
                for k,v in list(matchingNOESYaaindex_occupancy_dict.items()):
                    #print "residue i =",previous_TOCSY_aaindex," possible residue i-1 =",k," occupancy",str(v)+"/"+str(Num_of_TOCSY_resonances)
                    TOCSYaaindex_NOESYaaindex_OccupancyNumOfResonancesList_dict[previous_TOCSY_aaindex][k]=[v,Num_of_TOCSY_resonances]
                matchingNOESYaaindex_occupancy_dict = {}
                Num_of_TOCSY_resonances = 0
            TOCSY_H_resonance=float(TOCSY_words_list[1]) # aliphatic H resonance of residue i-1 
            TOCSY_C_resonance=float(TOCSY_words_list[2]) # aliphatic C (Ca,Cb,Cc,Cg,Ce,etc.) resonance of residue i-1; this C is covalently bonded to the above H
            matchingNOESYaaindex_list = []  # list with the NOESY index groups that match the current TOCSY index group
            Num_of_TOCSY_resonances += 1
            for q_index in range(0,len(NOESY_lines)):    # iterate over all NOESY resonances
                NOESY_line = NOESY_lines[q_index]
                try:
                    NOESY_words_list = NOESY_line.split()
                    
                    #NOESY_aaindex =int(float(NOESY_words_list[1]))  # residue (i-1)
                    NOESY_aaindex =str(NOESY_words_list[0])
                    NOESY_w2=float(NOESY_words_list[1]) # resonance of an aliphatic H close to residue (i-1)
                    NOESY_w3=float(NOESY_words_list[2]) # resonance of an aliphatic C close to residue (i-1); this C is covalently bonded to the above H
                    #print "DEBUG: NOESY_aaindex=",NOESY_aaindex,"NOESY_w2=",NOESY_w2,"NOESY_w3=",NOESY_w3
                    if ( (NOESY_w2-args.tolH) <= TOCSY_H_resonance <= (NOESY_w2+args.tolH) ) and ( (NOESY_w3-args.tolC) <= TOCSY_C_resonance <= (NOESY_w3+args.tolC) ):
                        #print "TOCSY AAIG ",TOCSY_aaindex," w2",TOCSY_H_resonance," w3",TOCSY_C_resonance," matching with NOESY AAIG ",NOESY_aaindex," w2",NOESY_w2," w3",NOESY_w3
                        if NOESY_aaindex in matchingNOESYaaindex_list:
                            print("WARNING: NOESY AAIG ",NOESY_aaindex, "was already found matching TOCSY AAIG ",TOCSY_aaindex,"!!!!!")
                        else:
                            matchingNOESYaaindex_list.append(NOESY_aaindex)
                except (ValueError, IndexError):
                    #print "WARNING: the 4th and 5th elements of the following NOESY file line are not numbers:"
                    #print "NOESY file line:", NOESY_line
                    continue
            
            # populate matchingNOESYaaindex_list using the matches of the current TOCSY line
            for aaindex in matchingNOESYaaindex_list:
                try:
                    matchingNOESYaaindex_occupancy_dict[aaindex] += 1 # if it exist alread increment its occupancy
                except KeyError:
                    matchingNOESYaaindex_occupancy_dict[aaindex] = 1 # if it does not exist initialize occupancy
            previous_TOCSY_aaindex = TOCSY_aaindex
        except (ValueError, IndexError):
            print("WARNING: the 2nd, 3rd or 4th elements of the following TOCSY file line are not numbers:")
            print("TOCSY file line:", TOCSY_words_list)
            continue
    
    # print the possible connectivities of the last TOCSY AAIG
    for k,v in list(matchingNOESYaaindex_occupancy_dict.items()):
        #print "residue i =",previous_TOCSY_aaindex," possible residue i-1 =",k," occupancy",str(v)+"/"+str(Num_of_TOCSY_resonances)
        TOCSYaaindex_NOESYaaindex_OccupancyNumOfResonancesList_dict[previous_TOCSY_aaindex][k]=[v,Num_of_TOCSY_resonances]
    
    return TOCSYaaindex_NOESYaaindex_OccupancyNumOfResonancesList_dict


################################################# FUNCTIONS TO FORM CHAINS FROM CONNECTIVITIES ########################################################
def populate_leaves(Assignment_Tree):
    """
        FUNCTION that adds new branches to the leaves of the Tree.
        ARGUMENTS:
        Assignment_Tree:    The Tree structure with connectivities
        RETURNS:
        (Assignment_Tree, BOOLEAN):    A tuple with elements the input Tree structure with new branches (if applicable), and a BOOLEAN value which is True if the function added
                                       new leaves to the Tree, or False otherwise
    """
    
    global i_iminus1_dict
    number_of_new_leaves = 0
    # ATTENTION: never use Assignment_Tree.iter_leaf_names(), it doesn't return the names in the order
    # ATTENTION: corresponding to Assignment_Tree.get_leaves()!!!
    for leaf in Assignment_Tree.get_leaves():
        try:
            #for child_tuple in i_iminus1_dict[name]:
            for child_triplet, child_duplet in zip(i_iminus1_dict[name], i_iminus1_normProbabilities_dict[name]): # WARNING: i_iminus1_normProbabilities_dict contains only matches above the Z-score cutoff
                NOESYaaindex = child_duplet[0]
                norm_probability = child_duplet[1]   # normalized probability of the connectivity 
                _occupancy = child_triplet[1]     # add prefix "_" to discriminate from new child feature "occupancy" that will be added
                _numOfResonances = child_triplet[2]
                ancestors_list = [ancestor.name for ancestor in leaf.get_ancestors()]
                if NOESYaaindex in ancestors_list:     # if the current NOESY AAIG is already a node or leaf in the Tree, continue to the next
                    continue
                #new_child = leaf.add_child(name=NOESYaaindex, dist=float(_occupancy)/_numOfResonances) # add a new brach to the current TOCSY add index (leaf) with length the ratio occupancy/numOfResonances
                
                new_child = leaf.add_child(name=NOESYaaindex, dist=norm_probability) # add a new brach to the current TOCSY add index (leaf) with length the respective probability
                
                new_child.add_features(occupancy=_occupancy, numOfResonances=_numOfResonances)
                number_of_new_leaves += 1
                #print "DEBUG: adding connection: ",name,"-->",NOESYaaindex
        except KeyError:
            continue
    
    # print Assignment_Tree.get_ascii(show_internal=True, compact=False, attributes=["name", "dist", "occupancy", "numOfResonances"])
    # print Assignment_Tree.get_ascii(show_internal=True, compact=False)
    if number_of_new_leaves > 0:
        return (Assignment_Tree, True)
    else:
        return (Assignment_Tree, False)
        

def build_Chain_Tree(i):
    
    global all_chainScore_set
    
    print("Building Tree starting from amino acid index",i,"...")
    expand_tree = True
    Assignment_Tree = Tree()
    Root = Assignment_Tree.get_tree_root()
    Root.add_feature("name", i)
    level = 1
    sys.stdout.write("Expanding tree from level ")
    while expand_tree:
        sys.stdout.write(str(level)+" ")
        sys.stdout.flush()
        Assignment_Tree, expand_tree = populate_leaves(Assignment_Tree)
        level += 1
        if level == args.MAX_PEPTIDE_LENGTH:
            break
    if level < args.MIN_PEPTIDE_LENGTH: # discard chains with a single AAIG
        return
    # Print the Tree
    # print Assignment_Tree.get_ascii(show_internal=True, compact=False)
    # print Assignment_Tree.get_ascii(show_internal=True, compact=False, attributes=["name", "dist", "occupancy", "numOfResonances"])
    
    print("\nSaving chains from Tree...")

    for leaf in Assignment_Tree.get_leaves():
        chain = []
        score = leaf.dist
        chain.append(leaf.name)
        for ancestor in leaf.get_ancestors():
            chain.append(ancestor.name)
            score *= ancestor.dist  # follow the chain rule for conditional probabilities to calculate the score (probability)
        # ATTENTION: the Tree is pruned in the sense that only connectivities above the Z-score cutoff were used, otherwise the scores of all chains would sum to 1 !
        # Therefore there is no need for normalization of chain scores!
        chain.append(score)
        all_chainScore_set.add(tuple(chain))
        del chain
        del ancestor
        del leaf
        #Assignment_Tree = None
    del Assignment_Tree
    gc.collect()


def is_sublist(lst1,lst2):
    try:
        ii = lst2.index(lst1[0])
    except ValueError:
        return False

    if (lst2[ii:ii+len(lst1)] == lst1) and (len(lst2)>len(lst1)):
        return True
    else:
        return False
#####################################################################################################################################



################################################# AUTOMATIC ASSIGNMENT OF AMINO ACID TYPES TO TOCSY RESONANCES FUNCTION DEFINITIONS ####################################################
def get_probability_from_histogram(resonance, bin_array, density_array):
    previous_hist_bin = bin_array[0]
    # print "DEBUG get_probability_from_histogram: resonance=", resonance
    # print "DEBUG get_probability_from_histogram: bin_array=", bin_array.tolist()
    # print "DEBUG get_probability_from_histogram: density_array=", density_array.tolist()
    for index, hist_bin in enumerate(np.nditer(bin_array[1:], order='K')):  # start counting from the 2nd element because we have appended "0" at the beginning of bin_list
        #print "DEBUG: previous_hist_bin=", previous_hist_bin, "resonance=", resonance,"hist_bin=",hist_bin
        if previous_hist_bin <= resonance and resonance < hist_bin:
            #print "DEBUG: index=",index, "previous_hist_bin=",previous_hist_bin,"resonance=",resonance,"hist_bin=",hist_bin
            probability = float(density_array[index])   # checked it and works correctly !
            break
        previous_hist_bin = hist_bin
    return probability


def get_probability_from_2Dhistogram(H_resonance, C_resonance, x_bin_array, y_bin_array, ):
    """
        FUNCTION to find the probability from a 2D histogram.
    """
    x_bin_list = list(set(x_bin_array))
    x_bin_list.sort()
    x_bin_length = round(x_bin_list[1]- x_bin_list[0], 2)
    y_bin_list = list(set(y_bin_array))
    y_bin_list.sort()
    y_bin_length = round(y_bin_list[1]- y_bin_list[0], 2)
    #print "DEBUG: x_bin_array=", x_bin_array.tolist()
    #print "DEBUG: y_bin_array=", y_bin_array.tolist()
    #print "DEBUG: probability_array=", probability_array
    index = 0
    probability = 0.0
    for x_hist_bin, y_hist_bin in zip(np.nditer(x_bin_array, order='K'), np.nditer(y_bin_array, order='K')):
        #print "DEBUG: x_hist_bin=", x_hist_bin, "H_resonance=", H_resonance,"x_hist_bin+x_bin_length=",x_hist_bin+x_bin_length
        #print "DEBUG: y_hist_bin=", y_hist_bin, "C_resonance=", C_resonance,"y_hist_binyx_bin_length=",y_hist_bin+y_bin_length
        if x_hist_bin <= H_resonance and H_resonance < x_hist_bin+x_bin_length and y_hist_bin <= C_resonance and C_resonance < y_hist_bin+y_bin_length:
            #print "DEBUG: index=",index, "x_hist_bin=", x_hist_bin, "H_resonance=", H_resonance,"x_hist_bin+x_bin_length=",x_hist_bin+x_bin_length
            #print "DEBUG: index=",index, "y_hist_bin=", y_hist_bin, "C_resonance=", C_resonance,"y_hist_binyx_bin_length=",y_hist_bin+y_bin_length
            probability = float(probability_array[index])   # I checked it and works correctly !
            break
        index += 1
    #print "DEBUG: probability=", probability
    return probability


def get_aatypes_from_H_C_resonpair_2Dhist(H_resonance, C_resonance, TOCSY_reson_index):
    """
        FUNCTION to find all possible amino acid types and the respective probabilities from a pair of aliphatic H,C resonances.
        ARGUMENTS:
        TOCSY_reson_index:  an index indicating the position of this H-C resonance pair in the TOCSY group it belongs to
        RETURNS:
        matches_list: list of lists of the form [aa type, Carbon name, Hydrogen name, weighted average probability, TOCSY_reson_index, H resonance, C resonance, carbon group]
    """
    global aa_CHpair_binProbabilityList_mdict, aa_carbon_binDensityList_mdict
    
    #print "DEBUG: H_resonance=",H_resonance,"C_resonance=",C_resonance
    matches_list = [] # list of tuples of the form (aa type, Carbon name, Hydrogen name, overall probability, )
    for aa in list(aa_CHpair_binProbabilityList_mdict.keys()):  # iterate over all amino acids
        carbon_matches_list = []    # list of tuples of the form (aa type,carbon,probability); here aa should be the same in all the tuples of the list!
        
        for carbon in list(aa_carbon_binDensityList_mdict[aa].keys()):    # iterate over all carbons of the current amino acid
            #print "DEBUG: carbon=", carbon
            bin_array, density_array = aa_carbon_binDensityList_mdict[aa][carbon]
            probability = get_probability_from_histogram(C_resonance, bin_array, density_array)
            #print "DEBUG: aa=",aa,"carbon=",carbon,"probability=",probability
            if probability > 0.0:
                carbon_matches_list.append((aa,carbon,probability))
        
        for CH_pair in list(aa_CHpair_binProbabilityList_mdict[aa].keys()):
            x_bin_array, y_bin_array, probability_array = aa_CHpair_binProbabilityList_mdict[aa][CH_pair]
            probability = get_probability_from_2Dhistogram(H_resonance, C_resonance, x_bin_array, y_bin_array, probability_array)
            # Keep only the amino acids for which were found both the carbon and the respective covalently bonded hydrogen
            # The probability of each prediction is given by the 2D histogram
            if probability > 0.0:
                C_name, H_name = CH_pair.split("-")
                matches_list.append([aa, C_name, H_name, probability, TOCSY_reson_index, H_resonance, C_resonance, None])  # add also the index of this C-H pair in the TOCSY group it
            elif probability == 0.0:    # if the probability of the 2D hist at this point is 0, use only the Carbon 1D histogram to get the probability
                C_name, H_name = CH_pair.split("-")
                #print "DEBUG: C_name=", C_name, "carbon_matches_list=", carbon_matches_list
                carbon_match = [match for match in carbon_matches_list if match[1]==C_name]
                if len(carbon_match) == 0:  # if the hist(C) of this C_name was zero, don't save it 
                    continue
                #print "DEBUG: len(carbon_match) is ", len(carbon_match)," and it should be 1 !"
                #print "DEBUG: carbon_match=", carbon_match
                weighted_average_probability = -1 * carbon_match[0][2] #  but first make it negative to distiguish it from weighted average probabilities
                #print "DEBUG: only the Carbon 1D histogram used:", [aa, C_name, H_name, weighted_average_probability, TOCSY_reson_index, H_resonance, C_resonance, None]
                matches_list.append([aa, C_name, H_name, weighted_average_probability, TOCSY_reson_index, H_resonance, C_resonance, None]) 
            
    #print "DEBUG: returning matches_list=", matches_list
    return matches_list


def get_aatypes_from_H_C_resonpair(H_resonance, C_resonance, TOCSY_reson_index):
    """
        FUNCTION to find all possible amino acid types and the respective probabilities from a pair of aliphatic H,C resonances.
        ARGUMENTS:
        TOCSY_reson_index:  an index indicating the position of this H-C resonance pair in the TOCSY group it belongs to
        RETURNS:
        valid_matches_list: list of lists of the form [aa type, Carbon name, Hydrogen name, weighted average probability, TOCSY_reson_index, H resonance, C resonance, carbon group]
    """
    global aa_carbon_binDensityList_mdict, aa_hydrogen_binDensityList_mdict, aa_carbonBondedHydrogensDict_dict
    
    def is_valid_hydrogen(aa_type, C_name, H_name, C_bondedH_dict):
        """
            FUNCTION to check whether a given hydrogen is covalently bonded to the given carbon of a given aa type
            ARGUMENTS:
            C_bondedH_dict:     dictonary with keys the carbon names of a given amino acid and values the respective covalently hydrogen names
        """
        
        if H_name in C_bondedH_dict[C_name]:
            return True
        
        return False
    
    # print "DEBUG: H_resonance=",H_resonance,"C_resonance=",C_resonance
    hydrogen_matches_list = []   # list of tuples of the form (aa type,hydrogen,probability)
    carbon_matches_list = []    # list of tuples of the form (aa type,carbon,probability)
    for aa in list(aa_hydrogen_binDensityList_mdict.keys()):  # iterate over all amino acids
        
        for hydrogen in list(aa_hydrogen_binDensityList_mdict[aa].keys()):    # iterate over all hydrogens of the current amino acid
            bin_array, density_array = aa_hydrogen_binDensityList_mdict[aa][hydrogen]
            probability = get_probability_from_histogram(H_resonance, bin_array, density_array)
            if probability > 0.0:
                hydrogen_matches_list.append((aa,hydrogen,probability))
            elif probability == 0.0:    # consider all cases where hist(H)=0
                hydrogen_matches_list.append((aa,hydrogen, 0.0))
        
        for carbon in list(aa_carbon_binDensityList_mdict[aa].keys()):    # iterate over all carbons of the current amino acid
            bin_array, density_array = aa_carbon_binDensityList_mdict[aa][carbon]
            probability = get_probability_from_histogram(C_resonance, bin_array, density_array)
            #print "DEBUG: aa=",aa,"carbon=",carbon,"probability=",probability
            if probability > 0.0:
                carbon_matches_list.append((aa,carbon,probability))
    
    # Keep only the amino acids for which were found both the carbon and the respective covalently bonded hydrogen
    # The probability of each prediction will be (Carbon_probability + Hydrogen_probability) /2
    matches_list = [] # list of tuples of the form (aa type, Carbon name, Hydrogen name, overall probability, )
    for hydrogen_match in hydrogen_matches_list:    # iterate over all hydrogen matches found in the previous step
        for carbon_match in carbon_matches_list:    # iterate over all carbons matches found in the previous step
            if hydrogen_match[0] == carbon_match[0]:    # if they belong to the same amino acid
                aa_type = hydrogen_match[0]
                H_name = hydrogen_match[1]
                C_name = carbon_match[1]
                if hydrogen_match[2] > 0.0:
                    if args.PROBABILITY_MODEL == 1: # treat the H and C probabilities as dependent events
                        # print "DEBUG: treating the H and C probabilities as dependent events."
                        if aa == "LEU":
                            weighted_average_probability = (1.0 * hydrogen_match[2] + 0.1 * carbon_match[2])/float((1.0 + 0.1))
                        else:
                            weighted_average_probability = (args.H_weight * hydrogen_match[2] + args.C_weight * carbon_match[2])/float((args.H_weight + args.C_weight))
                    elif args.PROBABILITY_MODEL == 2:
                        # print  "DEBUG: treating the H and C probabilities as independent events."
                        weighted_average_probability = hydrogen_match[2] * carbon_match[2]
                elif hydrogen_match[2] == 0.0:  # if hist(H)=0, consider only hist(C) but make it negative to distiguish it from weighted average probabilities
                    weighted_average_probability = -1 * carbon_match[2]
                matches_list.append((aa_type, C_name, H_name, weighted_average_probability))
    
    # TODO: calculate Z-scores from the histogram probabilities of the remaining unmatched carbons and ignore the hydrogens for those C with Z-score > 0.
    valid_matches_list = []    # the same list but contains only the most probable C-H matching pair from each aminoacid. 
    previous_aatype = None
    for quartet in matches_list:
        #print "DEBUG: quartet=",quartet
        aa_type = quartet[0]
        C_name = quartet[1]
        H_name = quartet[2]
        if is_valid_hydrogen(aa_type, C_name, H_name, aa_carbonBondedHydrogensDict_dict[aa_type]) == False: # if this carbon is not covalently bonded to this hydrogen, skip it
            continue
        #OBSOLETE # if aa_type == previous_aatype:  # elegant way to keep only the (aatype, C_name, H_name) with the highest probability, e.g. we may have:
        #OBSOLETE # ('VAL', 'CB', 'HB', 1.4736221632773356e-05), ('VAL', 'CG2', 'HG2', 0.01135282415016002), ('VAL', 'CG1', 'HG1', 0.01433564978569348) but we keep only the 3rd
        #OBSOLETE # othewise we would keep 'CG2_HG2', 'CG1_HG1', 'CB_HB' for just one C-H resonance pair!
        #OBSOLETE #    continue    
        valid_matches_list.append([aa_type, C_name, H_name, quartet[3], TOCSY_reson_index, H_resonance, C_resonance, None])  # add also the index of this H-C pair in the TOCSY group it
        # belongs to, but also the H and C resonances for future usage
        previous_aatype = aa_type
    
    #print "DEBUG: matches_list=", matches_list
    #print "DEBUG: valid_matches_list=", valid_matches_list
    return valid_matches_list


def select_correct_H_C_resonpairs(aatypeResonpairMatchesTuple_list):
    """
        FUNCTION  to match the TOCSY resonance pairs to as much as possible C-H resonance pairs of a PARTICULAR aatype. I may need to used a Genetic Algorith to do that
        if the current implementation proves to be error prone.
        
        ARGUMETS:
        aatypeResonpairMatchesTuple_list:   list of lists of the form [aa_type, C, H, weighted average probability, selected TOCSY reson index, H resonance, C resonance, clusterID] where aa_type is the same in the whole list (e.g. "PRO")
        
        RETURNS:
        correct_C_H_resonpair_matches_list: list of list of the form [aa_type, C, H, weighted average probability, selected TOCSY reson index], where aa_type is the same in the whole list (e.g. "PRO"). This list contains
                                            the final, uniq nucleus type assignment for the current amino acid. E.g.
        [['MET', 'CA', 'HA', 0.01993300618212495, 2, 3.843, 55.13, 1], ['MET', 'CB', 'HB2', 0.012634347058708313, 4, 1.266, 31.504, 2],
        ['MET', 'CG', 'HG2', 0.02281033340649995, 1, 1.508, 31.481, 2], ['MET', 'CB', 'HB3', 0.009737867955406978, 3, 1.911, 31.223, 3],
        ['MET', 'CG', 'HG3', 0.01607381733870664, 5, 0.403, 31.186, 3]]
    """
    global aa_carbonBondedHydrogensDict_dict
    
    def do_carbons_match(match, correct_C_H_resonpair_matches_list):
        """
            FUNCTION to check if the carbon resonances of geminal protons match. E.g. If he have added 
        """
        C_name = match[1]
        C_resonance = match[6]
        for correct_match in correct_C_H_resonpair_matches_list:
            # if the carbon name is already in the correct matches but its resonance is different, return false
            #print "DEBUG: C_name=",C_name, "correct_match[1]=", correct_match[1], "C_resonance=", C_resonance, "correct_match[6]=", correct_match[6]
            if C_name == correct_match[1] and approx_equal(C_resonance, correct_match[6], 0.3) == False:
                #print "DEBUG: conflict found, do_carbons_match() returning False!"
                return False
        
        return True # otherwise return true
        
        
    #print "DEBUG: entered function select_correct_H_C_resonpairs() with aatypeResonpairMatchesTuple_list=",aatypeResonpairMatchesTuple_list
    aa_type = aatypeResonpairMatchesTuple_list[0][0]    # recall that we have only one aa type in this script
    correct_C_H_resonpair_matches_list = [] # list of lists, see the Function definition
    used_TOCSY_reson_indices_list = []  # list of assigned peaks (TOCSY reson indices)
    used_C_H_pairs_list = []    # list of assigned C and H nucleus types
    
    # FIRST TAKE CARE OF CARBON GROUPING AND THEIR PROTONS (e.g. CB, CB2, CB3)
    clustIDs_list = [x[7] for x in aatypeResonpairMatchesTuple_list]    # list of cluster ID of each possible assignment
    clustIDs_set = set(clustIDs_list)
    # FIRST OF ALL SORT THE CLUSTERS ACCORDING TO THE TOTAL PROBABILITY OF THEIR CARBONS
    clustID_highestTotProb_dict = {}
    for clustID in clustIDs_set:
        TresonIndex_set = set([x[4] for x in aatypeResonpairMatchesTuple_list if x[7]==clustID])    # set with all different peaks (TOCSY reson indices) of the current cluster ID
        # NOW USE THE Treson indices WITH MULTIPLE TYPES OF CARBONS
        #print "DEBUG: clustIDs_list.count(clustID) =", clustIDs_list.count(clustID), " len(TresonIndex_set)=", len(TresonIndex_set)
        #if clustIDs_list.count(clustID) > 1 and len(TresonIndex_set) > 1:    # if there are 2 or more possible assignments and at least 2 different peaks within the same cluster ID
        nucleusAssignmentsforCurrentClustID_list = [x for x in aatypeResonpairMatchesTuple_list if x[7]==clustID]    # save all the members of this cluster ID
        # Recall that nucleusAssignmentsforCurrentClustID_list is already sorted by lines 1,2,3, BUT WE MUST ...
        nucleusAssignmentsforCurrentClustID_list.sort(key=itemgetter(3), reverse=True)  # ... SORT BY THE PREDICTION PROBABILITY
        
        # Make a list of the carbon types in each Treson index of this cluster and its maximum probability
        C_totalProb_dict = {}   # the total probability of each C type in this cluster. Use this to decide the assignment
        carbon_set = set([x[1] for x in nucleusAssignmentsforCurrentClustID_list])
        TresonIndex_set = set([x[4] for x in nucleusAssignmentsforCurrentClustID_list])
        for TRI in TresonIndex_set:
            nucleusAssignmentsforTRI_list = [x for x in nucleusAssignmentsforCurrentClustID_list if x[4]]
            nucleusAssignmentsforTRI_list.sort(key=itemgetter(3), reverse=True)  # sort by the prediction probability
            C_prob_dict = {}
            for prediction_list in nucleusAssignmentsforTRI_list:
                C = prediction_list[1]
                if C not in list(C_prob_dict.keys()):
                    C_prob_dict[C] = prediction_list[3]
            #print "DEBUG: TRI = ", TRI
            for C in carbon_set:
                if C not in list(C_prob_dict.keys()): # if this C type was not in possible predictions if this TRI, set its total probability to 0 (it will be completely excluded!)
                    C_totalProb_dict[C] = 0
                    continue
                try:
                    C_totalProb_dict[C] *= C_prob_dict[C]
                except KeyError:
                    C_totalProb_dict[C] = C_prob_dict[C]
        
        #print "DEBUG: C_totalProb_dict=", C_totalProb_dict
        sorted_C_totalProb_list = sorted(list(C_totalProb_dict.items()), key=itemgetter(1), reverse=True)   # sorted C_totalProb_list, e.g. [(CB, 0.9), (CG, 0.7)]
        
        clustID_highestTotProb_dict[clustID] = sorted_C_totalProb_list[0][1]    # save the total probability of the first element only, since the list was sorted
    # AND THIS IS THE LIST OF CLUSTER IDS SORTED BY THEIR HIGHEST TOTAL PROBABILITY
    #print "DEBUG: clustID_highestTotProb_dict = ", clustID_highestTotProb_dict
    sorted_clustIDs_list = [x[0] for x in sorted(list(clustID_highestTotProb_dict.items()), key=itemgetter(1), reverse=True)]
    #print "DEBUG: sorted_clustIDs_list=", sorted_clustIDs_list
    for clustID in clustIDs_set:
        if clustID not in sorted_clustIDs_list:
            sorted_clustIDs_list.append(clustID)
    #print "DEBUG: after addition of missing clusters sorted_clustIDs_list=", sorted_clustIDs_list
    
    
    CONTINUE_1ST_GROUPING = True # if no new assignments are saved, set to False to stop the iterations
    while CONTINUE_1ST_GROUPING:
        CONTINUE_1ST_GROUPING = False
        for clustID in sorted_clustIDs_list:
            TresonIndex_set = set([x[4] for x in aatypeResonpairMatchesTuple_list if x[7]==clustID])    # set with all different peaks (TOCSY reson indices) of the current cluster ID
            if clustIDs_list.count(clustID) > 1 and len(TresonIndex_set) > 1:    # if there are 2 or more possible assignments and at least 2 different peaks within the same cluster ID
                nucleusAssignmentsforCurrentClustID_list = [x for x in aatypeResonpairMatchesTuple_list if x[7]==clustID]    # save all the members of this cluster ID
                # Recall that nucleusAssignmentsforCurrentClustID_list is already sorted by lines 1,2,3
                # BUT FIRST CHECK IF ANY OF THE Treson INDICES CONTAINS A SINGLE TYPE OF CARBON. IN THAT CASE, THAT CARBON HAS PRIORITY IN THE ASSIGNMENT.
                for TresonIndex in TresonIndex_set:
                    carbon_list = [assignment[1] for assignment in nucleusAssignmentsforCurrentClustID_list if assignment[4]==TresonIndex]  # save all types of C in this Treson index
                    carbon_set = set(carbon_list)
                    if len(carbon_set) == 1:    # CONDITION: if there is only a single type of carbon within the Treson index 
                        assignment_list = [assignment for assignment in nucleusAssignmentsforCurrentClustID_list if assignment[4]==TresonIndex]  # save all types of C in this Treson index
                        assignment_list.sort(key=itemgetter(3), reverse=True)  # sort by the prediction probability
                        assignment = assignment_list[0] # USE ONLY THE ASSIGNMENT WITH THE HIGHEST PREDICTION PROBABILITY
                        C = assignment[1]
                        CARBON_OF_THIS_CLUSTER = C  # all assignments of this clustID must contain the this type of Carbon!!!
                        H = assignment[2]
                        # CONDITIONS:
                        # 1) TOCSY reson index has not been used again
                        # 2) this H must be covalently bonded to this C
                        # 3) this C,H pair has not been used again
                        if not TresonIndex in used_TOCSY_reson_indices_list and H in aa_carbonBondedHydrogensDict_dict[aa_type][C] and not C+"_"+H in used_C_H_pairs_list:
                            #print "1st GROUPING EXAMPLE: nucleusAssignmentsforCurrentClustID_list=", nucleusAssignmentsforCurrentClustID_list
                            correct_C_H_resonpair_matches_list.append(assignment)
                            used_TOCSY_reson_indices_list.append(TresonIndex)
                            used_C_H_pairs_list.append(C+"_"+H)
                            CONTINUE_1ST_GROUPING = True # since a new assignment has been saved, iterate the cycle once more
                            # CLEAN aatypeResonpairMatchesTuple_list from other possible assignments of the same clustID that do not contain this C type
                            for x in aatypeResonpairMatchesTuple_list:
                                # CONDITIONS:
                                # 1) the possible assignment must belong to the current clusterID
                                # 2) if the Carbon type in the assignment is not the CARBON_OF_THIS_CLUSTER, remove the assignment
                                if x[7] == clustID and x[1] != CARBON_OF_THIS_CLUSTER:
                                    aatypeResonpairMatchesTuple_list.remove(x)
        #print "DEBUG: after 1st grouping correct_C_H_resonpair_matches_list=", correct_C_H_resonpair_matches_list
    
    for clustID in sorted_clustIDs_list:
        TresonIndex_set = set([x[4] for x in aatypeResonpairMatchesTuple_list if x[7]==clustID])    # set with all different peaks (TOCSY reson indices) of the current cluster ID
        # NOW USE THE Treson indices WITH MULTIPLE TYPES OF CARBONS
        if clustIDs_list.count(clustID) > 1 and len(TresonIndex_set) > 1:    # if there are 2 or more possible assignments and at least 2 different peaks within the same cluster ID
            nucleusAssignmentsforCurrentClustID_list = [x for x in aatypeResonpairMatchesTuple_list if x[7]==clustID]    # save all the members of this cluster ID
            # Recall that nucleusAssignmentsforCurrentClustID_list is already sorted by lines 1,2,3, BUT WE MUST ...
            nucleusAssignmentsforCurrentClustID_list.sort(key=itemgetter(3), reverse=True)  # ... SORT BY THE PREDICTION PROBABILITY
            
            # Make a list of the carbon types in each Treson index of this cluster and its maximum probability
            C_totalProb_dict = {}   # the total probability of each C type in this cluster. Use this to decide the assignment
            carbon_set = set([x[1] for x in nucleusAssignmentsforCurrentClustID_list])
            TresonIndex_set = set([x[4] for x in nucleusAssignmentsforCurrentClustID_list])
            for TRI in TresonIndex_set:
                nucleusAssignmentsforTRI_list = [x for x in nucleusAssignmentsforCurrentClustID_list if x[4]]
                nucleusAssignmentsforTRI_list.sort(key=itemgetter(3), reverse=True)  # sort by the prediction probability
                C_prob_dict = {}
                for prediction_list in nucleusAssignmentsforTRI_list:
                    C = prediction_list[1]
                    if C not in list(C_prob_dict.keys()):
                        C_prob_dict[C] = prediction_list[3]
                
                for C in carbon_set:
                    if C not in list(C_prob_dict.keys()): # if this C type was not in possible predictions if this TRI, set its total probability to 0 (it will be completely excluded!)
                        C_totalProb_dict[C] = 0
                        continue
                    try:
                        C_totalProb_dict[C] *= C_prob_dict[C]
                    except KeyError:
                        C_totalProb_dict[C] = C_prob_dict[C]
            
            sorted_C_totalProb_list = sorted(list(C_totalProb_dict.items()), key=itemgetter(1), reverse=True)   # sorted C_totalProb_list, e.g. [(CB, 0.9), (CG, 0.7)]
            #print "DEBUG: clustID=", clustID, "sorted_C_totalProb_list=", sorted_C_totalProb_list
            for C_prob in sorted_C_totalProb_list:
                CARBON_OF_THIS_CLUSTER = C_prob[0] # all assignments of this clustID must contain this type of Carbon!!!
                
                for assignment in nucleusAssignmentsforCurrentClustID_list: # assignmend are sorted by probability; assignment example: ['MET', 'CG', 'HG3', 0.02342287175549055, 4, 1.266, 31.504, 2]
                    C = assignment[1]
                    if C != CARBON_OF_THIS_CLUSTER: # if this possible assignment contains a different C type, skip it
                        continue
                    H = assignment[2]
                    TresonIndex = assignment[4]
                    # CONDITIONS:
                    # 1) TOCSY reson index has not been used again
                    # 2) this H must be covalently bonded to this C
                    # 3) this C,H pair has not been used again
                    if not TresonIndex in used_TOCSY_reson_indices_list and H in aa_carbonBondedHydrogensDict_dict[aa_type][C] and not C+"_"+H in used_C_H_pairs_list:
                        #print "2nd GROUPING EXAMPLE: nucleusAssignmentsforCurrentClustID_list=", nucleusAssignmentsforCurrentClustID_list
                        correct_C_H_resonpair_matches_list.append(assignment)
                        used_TOCSY_reson_indices_list.append(TresonIndex)
                        used_C_H_pairs_list.append(C+"_"+H)
                        # CLEAN aatypeResonpairMatchesTuple_list from other possible assignments of the same clustID that do not contain this C type
                        for x in aatypeResonpairMatchesTuple_list:
                            # CONDITIONS:
                            # 1) the possible assignment must belong to the current clusterID
                            # 2) if the Carbon type in the assignment is not the CARBON_OF_THIS_CLUSTER, remove the assignment
                            if x[7] == clustID and x[1] != CARBON_OF_THIS_CLUSTER:
                                aatypeResonpairMatchesTuple_list.remove(x)
        #print "DEBUG: after 2nd grouping correct_C_H_resonpair_matches_list=", correct_C_H_resonpair_matches_list
    
    #print "DEBUG: after grouping all clustIDs correct_C_H_resonpair_matches_list=", correct_C_H_resonpair_matches_list
    
    
    # NOW TAKE CARE OF UNGROUPED CARBONS, SPECIFICALLY CLUSTERS WITH ONLY ONE Treson index AND ONLY ONE POSSIBLE ASSIGNMENT
    for clustID in sorted_clustIDs_list:
        TresonIndex_set = set([x[4] for x in aatypeResonpairMatchesTuple_list if x[7]==clustID])    # set with all different peaks (TOCSY reson indices) of the current cluster ID
        # NOW USE THE Treson indices WITH MULTIPLE TYPES OF CARBONS
        if clustIDs_list.count(clustID) == 1 and len(TresonIndex_set) == 1:    # if there is only 1 possible assignment and only 1 peak within the same cluster ID
            nucleusAssignmentsforCurrentClustID_list = [x for x in aatypeResonpairMatchesTuple_list if x[7]==clustID]    # save the single member of this cluster ID
            # SINCE IT IS ONLY ONE POSSIBLE ASSIGNMENT, WE DON'T NEED TO SORT nucleusAssignmentsforCurrentClustID_list
            assignment = nucleusAssignmentsforCurrentClustID_list[0]
            C = assignment[1]
            CARBON_OF_THIS_CLUSTER = C  # all assignments of this clustID must contain the this type of Carbon!!!
            H = assignment[2]
            TresonIndex = assignment[4]
            # CONDITIONS:
            # 1) TOCSY reson index has not been used again
            # 2) this H must be covalently bonded to this C
            # 3) this C,H pair has not been used again
            if not TresonIndex in used_TOCSY_reson_indices_list and H in aa_carbonBondedHydrogensDict_dict[aa_type][C] and not C+"_"+H in used_C_H_pairs_list:
                #print "3rd GROUPING EXAMPLE: nucleusAssignmentsforCurrentClustID_list=", nucleusAssignmentsforCurrentClustID_list
                correct_C_H_resonpair_matches_list.append(assignment)
                used_TOCSY_reson_indices_list.append(TresonIndex)
                used_C_H_pairs_list.append(C+"_"+H)
                # CLEAN aatypeResonpairMatchesTuple_list from other possible assignments of the same clustID that do not contain this C type
                for x in aatypeResonpairMatchesTuple_list:
                    # CONDITIONS:
                    # 1) the possible assignment must belong to the current clusterID
                    # 2) if the Carbon type in the assignment is not the CARBON_OF_THIS_CLUSTER, remove the assignment
                    if x[7] == clustID and x[1] != CARBON_OF_THIS_CLUSTER:
                        aatypeResonpairMatchesTuple_list.remove(x)
    #print "DEBUG: after clusters with a single TRI and single assignment correct_C_H_resonpair_matches_list=", correct_C_H_resonpair_matches_list
    
    # THEN, TAKE CARE OF UNGROUPED CARBONS, SPECIFICALLY CLUSTERS WITH ONLY ONE Treson index BUT MULTIPLE POSSIBLE ASSIGNMENTS
    for clustID in sorted_clustIDs_list:
        TresonIndex_set = set([x[4] for x in aatypeResonpairMatchesTuple_list if x[7]==clustID])    # set with all different peaks (TOCSY reson indices) of the current cluster ID
        # NOW USE THE Treson indices WITH MULTIPLE TYPES OF CARBONS
        if clustIDs_list.count(clustID) > 1 and len(TresonIndex_set) == 1:    # if there is only 1 possible assignment and only 1 peak within the same cluster ID
            nucleusAssignmentsforCurrentClustID_list = [x for x in aatypeResonpairMatchesTuple_list if x[7]==clustID]    # save the single member of this cluster ID
            TresonIndex = TresonIndex_set.pop() # substract the single TresonIndex of this clustID
            assignment_list = [assignment for assignment in nucleusAssignmentsforCurrentClustID_list if assignment[4]==TresonIndex]  # save all types of C in this Treson index
            assignment_list.sort(key=itemgetter(3), reverse=True)  # SORT BY THE PREDICTION PROBABILITY
            assignment = assignment_list[0]                        # USE ONLY THE ASSIGNMENT WITH THE HIGHEST PREDICTION PROBABILITY
            C = assignment[1]
            CARBON_OF_THIS_CLUSTER = C  # all assignments of this clustID must contain the this type of Carbon!!!
            H = assignment[2]
            TresonIndex = assignment[4]
            # CONDITIONS:
            # 1) TOCSY reson index has not been used again
            # 2) this H must be covalently bonded to this C
            # 3) this C,H pair has not been used again
            if not TresonIndex in used_TOCSY_reson_indices_list and H in aa_carbonBondedHydrogensDict_dict[aa_type][C] and not C+"_"+H in used_C_H_pairs_list:
                #print "4th GROUPING EXAMPLE: nucleusAssignmentsforCurrentClustID_list=", nucleusAssignmentsforCurrentClustID_list
                correct_C_H_resonpair_matches_list.append(assignment)
                used_TOCSY_reson_indices_list.append(TresonIndex)
                used_C_H_pairs_list.append(C+"_"+H)
                # CLEAN aatypeResonpairMatchesTuple_list from other possible assignments of the same clustID that do not contain this C type
                for x in aatypeResonpairMatchesTuple_list:
                    # CONDITIONS:
                    # 1) the possible assignment must belong to the current clusterID
                    # 2) if the Carbon type in the assignment is not the CARBON_OF_THIS_CLUSTER, remove the assignment
                    if x[7] == clustID and x[1] != CARBON_OF_THIS_CLUSTER:
                        aatypeResonpairMatchesTuple_list.remove(x)
    #print "DEBUG: after clusters with a single TRI but multiple assignments correct_C_H_resonpair_matches_list=", correct_C_H_resonpair_matches_list
    
    # FINALLY TAKE CARE OF WHATEVER IS LEFT IN aatypeResonpairMatchesTuple_list. DON'T KNOW IF THIS PART IS NECESSARY ANY MORE!
    C_H_occupancy_dict = {} # dict with keys C_H pairs and values the number of TOCSY peaks of the currect group this C_H pair matched
    C_H_matches_dict = {}
    for C, H_list in list(aa_carbonBondedHydrogensDict_dict[aa_type].items()):
        for H in H_list:
            #print "DEBUG:"+C+"_"+H
            if C+"_"+H in used_C_H_pairs_list:  # if this C,H pair has been used again, skip it
                continue
            matchesTuples_list = [t for t in aatypeResonpairMatchesTuple_list if t[1]==C and t[2]==H]    # keep all the tuples matching the current C-H pair
            C_H_occupancy_dict[C+"_"+H] = len(matchesTuples_list)
            C_H_matches_dict[C+"_"+H] = matchesTuples_list
    
    #print "DEBUG: C_H_occupancy_dict = ", C_H_occupancy_dict
    #print "DEBUG: C_H_matches_dict = ", C_H_matches_dict
    sorted_C_H_occupancyTuple_list = sorted(list(C_H_occupancy_dict.items()), key=itemgetter(1))  # the C_H_occupancyTuple tuples sorted by occupancy, so ambiguous assignments (higher occupancy) are left for the end
    #print "DEBUG: sorted_C_H_occupancyTuple_list=", sorted_C_H_occupancyTuple_list
    for C_H_occupancy_tuple in sorted_C_H_occupancyTuple_list:
        C_H_pair = C_H_occupancy_tuple[0]
        C = C_H_pair.split("_")[0]
        H = C_H_pair.split("_")[1]
        for match in C_H_matches_dict[C_H_pair]:   # we don't need to sort here the alternative C-H resonances according to the weighted average prob,
                                # because it has been done already at "possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list.sort(key=itemgetter(0,1,2,3), reverse=True)"
            #print "DEBUG: match = ", match
            #print "DEBUG: correct_C_H_resonpair_matches_list=", correct_C_H_resonpair_matches_list
            TOCSY_reson_index = match[4]
            # CONDITIONS:
            # 1) This Treson index has not been assigned yet
            # 2) This condition was useful before the incorporation of carbon grouping. I.e. it checked whether HB2-CB and HB3-CB had the same carbon resonance (+- 0.3)
            if not TOCSY_reson_index in used_TOCSY_reson_indices_list and do_carbons_match(match, correct_C_H_resonpair_matches_list):
                #print "5th GROUPING EXAMPLE: nucleusAssignmentsforCurrentClustID_list=", nucleusAssignmentsforCurrentClustID_list
                correct_C_H_resonpair_matches_list.append(match)
                used_TOCSY_reson_indices_list.append(TOCSY_reson_index)
                break   # move to the next C-H resonance pair
    
    #print "DEBUG: final correct_C_H_resonpair_matches_list=",correct_C_H_resonpair_matches_list
    return correct_C_H_resonpair_matches_list


def populate_leaves_for_correct_H_C_resonpairs2(Assignment_Tree, possible_aatype_prob_C_H_resonpair_TOCSYindex_assignID_list_list, TRI):
    """
        FUNCTION that adds new branches to the leaves of the Tree.
        ARGUMENTS:
        Assignment_Tree:    The Tree structure with connectivities
        RETURNS:
        (Assignment_Tree, BOOLEAN):    A tuple with elements the input Tree structure with new branches (if applicable), and a BOOLEAN value which is True if the function added
                                       new leaves to the Tree, or False otherwise
    """
    nucleusAssignmentsforCurrentTRI_list = [x for x in possible_aatype_prob_C_H_resonpair_TOCSYindex_assignID_list_list if x[4]==TRI] # save all the members of this cluster ID
    #TresonIndex_set = set([x[4] for x in possible_aatype_prob_C_H_resonpair_TOCSYindex_assignID_list_list if x[7]==clustID])    # set with all different peaks (TOCSY reson indices) of the current cluster ID
    assignments_list = []
    #if len(TresonIndex_set) == 1:
    assignments_list = nucleusAssignmentsforCurrentTRI_list
    number_of_new_leaves = 0
    # ATTENTION: never use Assignment_Tree.iter_leaf_names(), it doesn't return the names in the order
    # ATTENTION: corresponding to Assignment_Tree.get_leaves()!!!
    for leaf in Assignment_Tree.get_leaves():
        try:
            for assignment in assignments_list:
                _Ctype = assignment[1]
                _Htype = assignment[2]
                _Prob = assignment[3]
                _TRI = assignment[4]
                _H = assignment[5]
                _C = assignment[6]
                _clustID = assignment[7]
                assignID = assignment[8]
                ancestors_list = [ancestor.Ctype + "_" + ancestor.Htype for ancestor in leaf.get_ancestors()]
                if _Ctype + "_" + _Htype in ancestors_list:   # if this C-H type is already in the Tree, skip it
                    continue
                new_child = leaf.add_child(name=assignID) # add a new brach to the current TOCSY add index (leaf) with length the respective probability and return it
                new_child.add_features(Ctype=_Ctype, Htype=_Htype, Prob=_Prob, TRI=_TRI, H=_H, C=_C, clustID=_clustID)
                number_of_new_leaves += 1
        except KeyError:
            continue
    
    # print Assignment_Tree.get_ascii(show_internal=True, compact=False, attributes=["name", "Ctype", "Htype", "TRI"])
    # print Assignment_Tree.get_ascii(show_internal=True, compact=False)
    if number_of_new_leaves > 0:
        return (Assignment_Tree, True)
    else:
        return (Assignment_Tree, False)


def get_prediction_score(prob_clustID_tuple_list, mode=1):
    """
    For each C-group (aka clustID), let the individual probabilities are prob1 and prob2, this function will calculate a consensus probability
    in one of the following ways (determined by 'mode' variable), and at the end return the product if these consensus probabilities.
    mode:   1: average;
            2: sqrt(prob1*prob2)    ; geometric mean
            3: (prob2-prob1)/(log(prob2)-log(prob1))    ; logarithmic mean
            4: (prob1*prob2)/((prob1+prob2)/2.) ;composite average
            5: prob1*prob2    ; product of probabilities
    """
    
    def get_Cgroup_consensus_probability(prob_list, mode=1):
        """
        ARGUMENTS:
        prob_list:  list of the peak probabilities of this C-group
        mode:       controls how to calculate the consensus probability of this C-group.
                    1: average;
                    2: sqrt(prob1*prob2)    ; geometric mean
                    3: (prob2-prob1)/(log(prob2)-log(prob1))    ; logarithmic mean
                    4: (prob1*prob2)/((prob1+prob2)/2.) ;composite average 
                    5: prob1*prob2    ; product of probabilities
        """
        consensus_probability = None
        
        if len(prob_list) == 1:
            prob1 = prob_list[0]
            print("DEBUG get_Cgroup_consensus_probability: only one probability is available prob1=", prob1, "mode=", mode)
            consensus_probability = prob1
        elif len(prob_list) == 2:
            prob1 = prob_list[0]
            prob2 = prob_list[1]
            print("DEBUG get_Cgroup_consensus_probability: prob1=", prob1, "prob2=", prob2, "mode=", mode)
            if mode == 1:   # arithmetic average
                consensus_probability = (prob1 + prob2)/2.
            elif mode == 2: # geometric average
                consensus_probability = math.sqrt(prob1*prob2)
            elif mode == 3: # logarithmic average
                if prob1 == 0 and prob2 == 0:
                    consensus_probability = 0
                elif prob1 == prob2:
                    consensus_probability = prob1
                else:
                    consensus_probability = (prob2-prob1)/(math.log10(prob2)-math.log10(prob1))
            elif mode == 4: # composite average
                consensus_probability = (prob1*prob2)/((prob1+prob2)/2.)
            elif mode == 5: # product of probabilities
                consensus_probability = prob1*prob2
        elif len(prob_list) > 2:
            print("ERROR: this Carbon group contains more than two peaks!!!!")
            print("DEBUG: prob_list=", prob_list)
            sys.exit(1)
        
        print("DEBUG: consensus_probability=", consensus_probability)
        return consensus_probability
    
    prediction_score = 1 # product of the consensus probabilites of each C-group
    Cgroup_set = set([t[1] for t in prob_clustID_tuple_list])
    for Cgroup in Cgroup_set:
        prob_set = set([t[0] for t in prob_clustID_tuple_list if t[1]==Cgroup]) # get the peak probabilities of this C-group
        prob_list = list(prob_set)
        prediction_score *= get_Cgroup_consensus_probability(prob_list, mode)
    
    print("DEBUG: prediction_score=", prediction_score)
    return prediction_score


def select_correct_H_C_resonpairs2(aatypeResonpairMatchesTuple_list):
    """
        NEW FUNCTION  to match the TOCSY resonance pairs to as much as possible C-H resonance pairs of a PARTICULAR aatype. This functions selects the best C-H type assignment
        combination based on the product of probabilities of the individual C-H assignments.
        
        ARGUMETS:
        aatypeResonpairMatchesTuple_list:   list of lists of the form [aa_type, C, H, weighted average probability, selected TOCSY reson index, H resonance, C resonance, clusterID] where aa_type is the same in the whole list (e.g. "PRO")
        
        RETURNS:
        correct_C_H_resonpair_matches_list: list of lists of the form [aa_type, C, H, weighted average probability, selected TOCSY reson index], where aa_type is the same in the whole list (e.g. "PRO"). This list contains
                                            the final, uniq nucleus type assignment for the current amino acid. E.g.
        [['MET', 'CA', 'HA', 0.01993300618212495, 2, 3.843, 55.13, 1], ['MET', 'CB', 'HB2', 0.012634347058708313, 4, 1.266, 31.504, 2],
        ['MET', 'CG', 'HG2', 0.02281033340649995, 1, 1.508, 31.481, 2], ['MET', 'CB', 'HB3', 0.009737867955406978, 3, 1.911, 31.223, 3],
        ['MET', 'CG', 'HG3', 0.01607381733870664, 5, 0.403, 31.186, 3]]
        aa_type:                            recall that we run this function for only one aa type!
        prediction_score:                   the total score of this amino acid type
    """
    global aa_carbonBondedHydrogensDict_dict
    
    possible_aatype_prob_C_H_resonpair_TOCSYindex_assignID_list_list = []   # same as possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list but with an assignment extra index at the end
    for i,x in enumerate(aatypeResonpairMatchesTuple_list):
        possible_aatype_prob_C_H_resonpair_TOCSYindex_assignID_list_list.append([x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], i])
    print("DEBUG select_correct_H_C_resonpairs2: possible_aatype_prob_C_H_resonpair_TOCSYindex_assignID_list_list=", possible_aatype_prob_C_H_resonpair_TOCSYindex_assignID_list_list)
    
    clustIDs_list = [x[7] for x in possible_aatype_prob_C_H_resonpair_TOCSYindex_assignID_list_list]    # list of cluster ID of each possible assignment
    clustIDs_set = set(clustIDs_list)
    TRI_set = set([x[4] for x in possible_aatype_prob_C_H_resonpair_TOCSYindex_assignID_list_list])
    TRI_num = len(TRI_set)
    print("DEBUG: TRI_set =", TRI_set)
    
    print("Building Tree of possible C-H type assignment combinations")
    print("DEBUG: possible_aatype_prob_C_H_resonpair_TOCSYindex_assignID_list_list=", possible_aatype_prob_C_H_resonpair_TOCSYindex_assignID_list_list)
    aa_type = possible_aatype_prob_C_H_resonpair_TOCSYindex_assignID_list_list[0][0]
    expand_tree = True
    Assignment_Tree = Tree()
    Root = Assignment_Tree.get_tree_root()
    Root.add_features(name=aa_type, Ctype="None", Htype="None")
    level = 0
    sys.stdout.write("Expanding tree from level ")
    for TRI in TRI_set:
        sys.stdout.write("level="+str(level)+" ")
        sys.stdout.flush()
        Assignment_Tree, expand_tree = populate_leaves_for_correct_H_C_resonpairs2(Assignment_Tree, possible_aatype_prob_C_H_resonpair_TOCSYindex_assignID_list_list, TRI)
        if expand_tree == False:
            break
        level += 1
        if level == TRI_num:
            break
    # if level < TRI_num: # discard combinations with less than the number of peaks
    #    print "ERROR: NOT ALL C-GROUPS WERE USED IN THE ASSIGNMENT OF ", aa_type
    #    sys.exit(1)
    
    
    #print "\nSaving chains from Tree..."
    all_chainScore_list = []
    for leaf in Assignment_Tree.get_leaves():
        chain = []
        prob_product = leaf.Prob
        prob_clustID_tuple_list = [(leaf.Prob, leaf.clustID)]
        assignID = int(leaf.name)
        chain.append(possible_aatype_prob_C_H_resonpair_TOCSYindex_assignID_list_list[assignID])
        print("DEBUG: leaf.name=", leaf.name)
        for ancestor in leaf.get_ancestors():
            print("DEBUG: ancestor.name=", ancestor.name)
            try:
                assignID = int(ancestor.name)
            except ValueError:
                continue
            chain.append(possible_aatype_prob_C_H_resonpair_TOCSYindex_assignID_list_list[assignID])
            prob_product *= ancestor.Prob
            prob_clustID_tuple_list.append((ancestor.Prob, ancestor.clustID))
        if args.CONSENSUS_CGROUP_PROB_MODE == 0:
            chain.append(prob_product)
        else:
            chain.append(get_prediction_score(prob_clustID_tuple_list, mode=args.CONSENSUS_CGROUP_PROB_MODE))
        print("DEBUG: saving chain=", chain)

        all_chainScore_list.append(tuple(chain))
        del chain
        del ancestor
        del leaf
    del Assignment_Tree
    gc.collect()
    
    
    correct_chainScore_list = []
    print("DEBUG: all_chainScore_list=", all_chainScore_list)
    carbon_exceptions_list = ["CD1", "CD2"] # carbons that may belong to the same group but have different names
    aatype_Ctype_exceptions_dict = {}   # not in use yet
    aatype_Ctype_exceptions_dict['GLY'] = ["CA"]
    for x in all_chainScore_list:
        DISCARD = False
        # discard the prediction if it does not contain all the TRIs
        if len(x) < len(TRI_set)+1 :
            DISCARD = True
        # check the C-groups
        Ctype2ClustID_dict = {}
        for clustID in clustIDs_set:
            Ctype_list = list(set([a[1] for a in x[:-1] if a[7]==clustID]))    # C-type(s) of this cluster
            if len(Ctype_list) > 1 and not Ctype_list[0] in carbon_exceptions_list and not Ctype_list[1] in carbon_exceptions_list:
                DISCARD = True
                break
        # check that each C-group contains a different C-type
        Ctype_set = set([a[1] for a in x[:-1]])
        for Ctype in Ctype_set:
            if len(set([a[7] for a in x[:-1] if a[1]==Ctype])) > 1 and (a[0] != 'GLY' and a[1] != 'CA'):
                DISCARD = True
                break
        # check that each C-type belongs to only one C-group
        Ctype_set = set([a[1] for a in x[:-1]])
        for Ctype in Ctype_set:
            Cgroups_list = list(set([a[-2] for a in x[:-1] if a[1]==Ctype]))    # C-groups of this cluster
            if len(Cgroups_list) > 1:
                DISCARD = True
                break
        # check for H-type duplication
        Htype_set = set([a[2] for a in x[:-1]])
        Htype_list = [a[2] for a in x[:-1]]
        if len(Htype_set) !=  len(Htype_list):
            DISCARD = True
        
        if DISCARD:
            continue
        else:
            correct_chainScore_list.append(x)
    
    
    correct_chainScore_list.sort(key=itemgetter(len(TRI_set)), reverse=True) # sort probability by descending order; the highest probability combination must be 1st
    print("DEBUG: all possible predictions correct_chainScore_list=", correct_chainScore_list)
    
    print("\nDEBUG: Best Prediction:")
    try:
        for x in correct_chainScore_list[0][:-1]:
            if x[3] <0:
                print("DEBUG: Best Prediction contains C-based probability!!!")
            print(x)
    except IndexError:
        print("DEBUG: correct_chainScore_list=", correct_chainScore_list)
        #sys.exit(1)    # temporarily inactive
        return [], aa_type, 0.0
    return correct_chainScore_list[0][:-1], aa_type, correct_chainScore_list[0][-1]



def get_combined_weighted_and_presence_probability(AA_type, aatype_probsum_dict, aatype_CHnucleiType_presenceProbSum_mdict):
    """
    FUNCTION to calculate the probability of an amino acid type AA_type given an assignment, as in eq. 6 of RESCUE2 paper. But instead of multiplying all C & H probabilities,
    we multiply the weighted averages of C-H pair probabilities.
    
    ARGUMENTS:
    aatype_probsum_dict:    in this case it is a dictionary with the product of the weighted average probabilities of all C-H pairs assigned.
    aatype_CHnucleiType_presenceProbSum_mdict:  multidimensiona dict of the form aa type --> C_H pair --> weighted average presence probabilities as in eq. 6 of RESCUE2 paper.
                                                    It misses the absence probabilities, which must be computed afterwards.
    """
    global aa_carbonBondedHydrogensDict_dict, Prob_CS
    
    carbon_bondexHydrogensList_dict = aa_carbonBondedHydrogensDict_dict[AA_type]
    
    for C_name in list(carbon_bondexHydrogensList_dict.keys()):
        for H_name in carbon_bondexHydrogensList_dict[C_name]:
            #print "DEBUG: checking ", C_name+"_"+H_name, "of aa type", AA_type
            if not C_name+"_"+H_name in list(aatype_CHnucleiType_presenceProbSum_mdict[AA_type].keys()):
                # if there is no weighted presence probability for the current C-H pair, save the respective absence probability
                aatype_CHnucleiType_presenceProbSum_mdict[AA_type][C_name+"_"+H_name] = 1 - (args.H_weight * Prob_CS[AA_type][H_name] + args.C_weight * Prob_CS[AA_type][C_name])/float((args.H_weight + args.C_weight))
                #print "DEBUG: Assigning absence to ", C_name+"_"+H_name, ":1 - *(",args.H_weight,"*", Prob_CS[AA_type][H_name], "+", args.C_weight, "*", Prob_CS[AA_type][C_name], ")/float((args.H_weight + args.C_weight)) =", aatype_CHnucleiType_presenceProbSum_mdict[AA_type][C_name+"_"+H_name]
    
    #print "DEBUG: aatype_CHnucleiType_presenceProbSum_mdict[AA_type]=", aatype_CHnucleiType_presenceProbSum_mdict[AA_type].items()
    for weighted_presence_probabilty in list(aatype_CHnucleiType_presenceProbSum_mdict[AA_type].values()):
        try:
            weighted_presence_probabilities_product *= weighted_presence_probabilty
        except NameError:
            weighted_presence_probabilities_product = weighted_presence_probabilty
    
    #print "DEBUG: final probability=", weighted_presence_probabilities_product * aatype_probsum_dict[AA_type]
    return weighted_presence_probabilities_product * aatype_probsum_dict[AA_type]


def consider_only_Carbon_matches(possible_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list, only_missing=True, remove_Cbased_predictions=False):
    """
        FUNCTION to check if there are C-H resonance pairs in the currect TAAIG group for every aa type, and to allow C-H type prediction only
        from C resonances when one of the following conditions is met:
                1) The H resonance is negative.
                2) No aa type prediction can be made for a particular TAAIG group.
        
        ARGUMENTS:
        possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list:   list of lists of the form (aa type, Carbon name, Hydrogen name, weighted average probability, TOCSY_reson_index, H resonance, C resonance, carbon group)
                                                                        This list includes also negative probabilities wherever only the Carbon resonance match was considered.
        possible_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list:   list of tuples of the form (aa type, Carbon name, Hydrogen name, weighted average probability, TOCSY_reson_index, H resonance, C resonance)
        only_missing:   if True only the C hist probabilities will be considered in C-H resonance cases where no C-H type prediction could be made for a particular aa. If False,
                        then only the C hist probabilities will be considered in every C-H resonance case. E.g. If we have the following predictions:
                        ('THR', 'CB', 'HB', -0.007423257864062593, 1, 4.239, 61.662), ('THR', 'CA', 'HA', 0.022600159959100728, 1, 4.239, 61.662)
                        ('THR', 'CB', 'HB', -0.020221544866552674, 2, 1.634, 68.57), ('THR', 'CA', 'HA', -0.0014439440299790283, 2, 1.634, 68.57)
                        ('THR', 'CG2', 'HG2', 0.001264928291792422, 3, 0.836, 26.318)
                        If only_missing=True, the allowed predictions will become:
                        ('THR', 'CA', 'HA', 0.022600159959100728, 1, 4.239, 61.662)
                        ('THR', 'CB', 'HB', -0.020221544866552674, 2, 1.634, 68.57), ('THR', 'CA', 'HA', -0.0014439440299790283, 2, 1.634, 68.57)
                        ('THR', 'CG2', 'HG2', 0.001264928291792422, 3, 0.836, 26.318)
                        But if only_missing=False, they will remain the same:
                        ('THR', 'CB', 'HB', -0.007423257864062593, 1, 4.239, 61.662), ('THR', 'CA', 'HA', 0.022600159959100728, 1, 4.239, 61.662)
                        ('THR', 'CB', 'HB', -0.020221544866552674, 2, 1.634, 68.57), ('THR', 'CA', 'HA', -0.0014439440299790283, 2, 1.634, 68.57)
                        ('THR', 'CG2', 'HG2', 0.001264928291792422, 3, 0.836, 26.318)
    """
    #print "DEBUG: possible_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list=", possible_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list
    possible_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list.sort(key=itemgetter(0,4), reverse=True)  # sort by aa type and Treson index
    previous_aatype = None
    aatype_rest7_dict = {}  # dict of the form: aa type -> (Carbon name, Hydrogen name, weighted average probability, TOCSY_reson_index, H resonance, C resonance)
    for group8 in possible_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list:
        aa_type = group8[0]
        if aa_type != previous_aatype and previous_aatype != None:
            aatype_rest7_dict[aa_type] = [(group8[1:])]
        elif aa_type == previous_aatype or previous_aatype == None:
            try:
                aatype_rest7_dict[aa_type].append((group8[1:]))
            except KeyError:
                aatype_rest7_dict[aa_type] = [(group8[1:])]
            
        previous_aatype = aa_type    
    
    new_possible_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list = [] # possible_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list after filtering
    for aa_type in list(aatype_rest7_dict.keys()):
        #print "DEBUG: aa_type=", aa_type, "len(TresonIndex_set)=", len(TresonIndex_set)
        previous_TresonIndex = None
        C_H_type_prediction_list = []   # list with the valid 
        for group7 in aatype_rest7_dict[aa_type]:
            #print "DEBUG: aa_type=", aa_type, "group7", group7
            if remove_Cbased_predictions == True and group7[2] < 0:
                continue
            if only_missing == False:   # if instructed add both C- and C-H-based aa type predictions
                new_possible_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list.append( (aa_type, group7[0], group7[1], abs(group7[2]), group7[3], group7[4], group7[5], group7[6]) )
            
            elif only_missing == True:
                # Recall the TresonIndex numbers are sorted in aatype_rest7_dict as a result of sorting possible_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list
                #print "DEBUG: previous_TresonIndex=", previous_TresonIndex, "group7[3]=", group7[3]
                if previous_TresonIndex != group7[3] and previous_TresonIndex != None:
                    # do the appropriate changes and save the C-H type predictions
                    positive_C_H_type_prediction_list = [g7 for g7 in C_H_type_prediction_list if g7[2] > 0]    # list of the C-H type predictions with positive weighted average probabilities
                    #print "DEBUG: positive_C_H_type_prediction_list=", positive_C_H_type_prediction_list
                    for C_H_type_prediction in C_H_type_prediction_list:
                        # ADD ALL POSITIVE PREDICTIONS AND THOSE NEGATIVE PREDICTIONS WITH PROBABILITY ABOVE THE AVERAGE PROBABILITY OF ALL POSITIVE PREDICTIONS, ONLY IF THE H RESONACNE
                        # IS WITHIN THE LOW OCCUPANCY REGIONS
                        if C_H_type_prediction[3] > 0 or ( (C_H_type_prediction[3] < 0 and is_H_in_low_occupancy_region(C_H_type_prediction[5]) ) and ( len(positive_C_H_type_prediction_list) == 0 or abs(C_H_type_prediction[3]) >= np.average([x[3] for x in positive_C_H_type_prediction_list]) ) ):
                            lst = list(C_H_type_prediction)
                            if (C_H_type_prediction[3] < 0 and is_H_in_low_occupancy_region(C_H_type_prediction[5])):
                                lst[5] = str(lst[5]) + "C"   # MARK THE H RESONANCE TO KNOW THAT IT WAS NOT USED FOR AA TYPE PREDICTION (marking C resonance will raise error)
                            lst[3] = abs(C_H_type_prediction[3])
                            C_H_type_prediction = tuple(lst)
                            new_possible_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list.append(C_H_type_prediction)
                    
                    # now reinitialize C_H_type_prediction_list
                    C_H_type_prediction_list = [ (aa_type, group7[0], group7[1], group7[2], group7[3], group7[4], group7[5], group7[6]) ]
                
                elif previous_TresonIndex == group7[3] or previous_TresonIndex == None:
                    C_H_type_prediction_list.append( (aa_type, group7[0], group7[1], group7[2], group7[3], group7[4], group7[5], group7[6]) )
        
            previous_TresonIndex = group7[3]
        
        if only_missing == True:  # NOW PROCESS THE LAST TresonIndex
            # do the appropriate changes and save the C-H type predictions
            positive_C_H_type_prediction_list = [g7 for g7 in C_H_type_prediction_list if g7[2] > 0]    # list of the C-H type predictions with positive weighted average probabilities
            #print "DEBUG: positive_C_H_type_prediction_list=", positive_C_H_type_prediction_list
            for C_H_type_prediction in C_H_type_prediction_list:
                # ADD ALL POSITIVE PREDICTIONS AND THOSE NEGATIVE PREDICTIONS WITH PROBABILITY ABOVE THE AVERAGE PROBABILITY OF ALL POSITIVE PREDICTIONS, ONLY IF THE H RESONACNE
                # IS NEGATIVE
                if C_H_type_prediction[3] > 0 or ( (C_H_type_prediction[3] < 0 and is_H_in_low_occupancy_region(C_H_type_prediction[5])) and ( len(positive_C_H_type_prediction_list) == 0 or abs(C_H_type_prediction[3]) >= np.average([x[3] for x in positive_C_H_type_prediction_list]) ) ):
                    lst = list(C_H_type_prediction)
                    if (C_H_type_prediction[3] < 0 and is_H_in_low_occupancy_region(C_H_type_prediction[5])):
                        lst[5] = str(lst[5]) + "C"   # MARK THE H RESONANCE TO KNOW THAT IT WAS NOT USED FOR AA TYPE PREDICTION (marking C resonance will raise error)
                    lst[3] = abs(C_H_type_prediction[3])
                    C_H_type_prediction = tuple(lst)
                    new_possible_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list.append(C_H_type_prediction)
        
    #print "DEBUG: new_possible_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list="
    #for t in new_possible_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list:
    #    print t
    
    return new_possible_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list


def get_aatypes_from_all_H_C_resonpairs(raw_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list, Num_of_TOCSY_resonances):
    """
        FUNCTION to return the matched amino acid types for a particular TOCSY AAIG by using all its H,C resonance pairs.
        
        ARGUMENTS:
        raw_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list:   list of lists of the form (aa type, Carbon name, Hydrogen name, weighted average probability, TOCSY_reson_index, H resonance, C resonance, carbon group)
                                                                        This list includes also negative probabilities wherever only the Carbon resonance match was considered.
        
        RETURN:
        matching_aaTypesProbTuple_list:     list of tuples of the form (aa type, average probability)
    """
    global aatype_maxH_C_pairs_dict, Prob_CS
    
    #print "DEBUG: raw_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list=", raw_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list
    # If the last arguments in True, it will remove all C-based predictions
    possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list = consider_only_Carbon_matches(raw_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list, True, False)
    #print "DEBUG get_aatypes_from_all_H_C_resonpairs: possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list=",possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list
    # The following sorting ensures that the alternative C-H pairs are sorted according to their weighted average probability 
    possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list.sort(key=itemgetter(0,1,2,3), reverse=True)  # sort by every column except the last 4
    
    ## The following block keeps only the matches for a particular aatype
    TRI_set = set([x[4] for x in possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list])
    TRI_num = len(TRI_set)  # the number of different peaks available for this spin system
    previous_aatype = None
    aatypeResonpairMatchesTuple_list = []
    selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list = []
    aatype_predictionScore_dict = {}     # dict aa type -> total prediction score of this aa type. Used only if args.PROB_PRODUCT == True !
    for group8 in possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list:
        #print "DEBUG: group8=", group8
        aa_type = group8[0]
        weighted_average_probability = group8[3]    # = [wH * hist(H) + wC * hist(C)]/2
        if aa_type != previous_aatype and previous_aatype != None:
            #print "DEBUG: aa_type=",aa_type,"previous_aatype=",previous_aatype,"aatypeResonpairMatchesTuple_list=",aatypeResonpairMatchesTuple_list
            if args.PROB_PRODUCT == True:
                new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list, aatype2save, aatype_predictionScore = select_correct_H_C_resonpairs2(aatypeResonpairMatchesTuple_list)
                aatype_predictionScore_dict[aatype2save] = aatype_predictionScore
            else:
                new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list = select_correct_H_C_resonpairs(aatypeResonpairMatchesTuple_list)
            # CHECK IF ALL THE PEAKS HAVE BEEN ASSIGNED TO C-H TYPES FOR THE CURRENT AA TYPE. IF NOT DO NOT SAVE THIS AA TYPE AS POSSIBLE PREDICTION FOR THE CURRENT TOCSY index group
            if len(new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list) >= (args.ASSIGNMENT_CUTOFF * TRI_num):
                selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list.extend(new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list)
            else:
                pass
                #print "DEBUG: the assigned C-H pairs for aa type", previous_aatype, "were below the requested threshold (",args.ASSIGNMENT_CUTOFF * TRI_num,")"
                #print "DEBUG: new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list=", new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list
            #print "DEBUG: selected correct H-C resonpairs!"
            aatypeResonpairMatchesTuple_list = []
        
        aatypeResonpairMatchesTuple_list.append(group8)
        previous_aatype = aa_type
    #print "DEBUG get_aatypes_from_all_H_C_resonpairs: aa_type=",aa_type,"previous_aatype=",previous_aatype,"aatypeResonpairMatchesTuple_list=",aatypeResonpairMatchesTuple_list
    if len(selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list) > 0:
    # SELECT THE OPTIMUM C-H COMBINATIONS FOR THE LAST aa_type
        if args.PROB_PRODUCT == True:
            new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list, aatype2save, aatype_predictionScore = select_correct_H_C_resonpairs2(aatypeResonpairMatchesTuple_list)  # do the last aa_type
            aatype_predictionScore_dict[aatype2save] = aatype_predictionScore
        else:
            new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list = select_correct_H_C_resonpairs(aatypeResonpairMatchesTuple_list)  # do the last aa_type
        # CHECK IF ALL THE PEAKS HAVE BEEN ASSIGNED TO C-H TYPES FOR THE CURRENT AA TYPE. IF NOT DO NOT SAVE THIS AA TYPE AS POSSIBLE PREDICTION FOR THE CURRENT TOCSY index group
        if len(new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list) >= (args.ASSIGNMENT_CUTOFF * TRI_num):
            selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list.extend(new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list)
        else:
            pass
            #print "DEBUG: the assigned C-H pairs for aa type", previous_aatype, "were below the requested threshold (",args.ASSIGNMENT_CUTOFF * TRI_num,")"
            #print "DEBUG: new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list=", new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list
        #print "DEBUG: selected correct H-C resonpairs!"
    # IF NO C-H TYPE ASSIGNMENTS WERE MADE, REPEAT THE PROCEDURE INCLUDING ONLY CARBON RESONANCES
    elif len(selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list) == 0:
        possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list = consider_only_Carbon_matches(raw_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list, False, False)
        # The following sorting ensures that the alternative C-H pairs are sorted according to their weighted average probability 
        possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list.sort(key=itemgetter(0,1,2,3), reverse=True)  # sort by every column except the last 3
        ## The following block keeps only the matches for a particular aatype
        previous_aatype = None
        aatypeResonpairMatchesTuple_list = []
        selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list = []
        #print "DEBUG: 2032 possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list=", possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list
        for group8 in possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list:
            #print "DEBUG: group8=", group8
            aa_type = group8[0]
            weighted_average_probability = group8[3]    # = [wH * hist(H) + wC * hist(C)]/2
            if aa_type != previous_aatype and previous_aatype != None:
                #print "DEBUG: aa_type=",aa_type,"previous_aatype=",previous_aatype,"aatypeResonpairMatchesTuple_list=",aatypeResonpairMatchesTuple_list
                if args.PROB_PRODUCT == True:
                    new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list, aatype2save, aatype_predictionScore = select_correct_H_C_resonpairs2(aatypeResonpairMatchesTuple_list)
                    aatype_predictionScore_dict[aatype2save] = aatype_predictionScore
                else:
                    new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list = select_correct_H_C_resonpairs(aatypeResonpairMatchesTuple_list)
                # CHECK IF ALL THE PEAKS HAVE BEEN ASSIGNED TO C-H TYPES FOR THE CURRENT AA TYPE. IF NOT DO NOT SAVE THIS AA TYPE AS POSSIBLE PREDICTION FOR THE CURRENT TOCSY index group
                if len(new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list) >= (args.ASSIGNMENT_CUTOFF * TRI_num):
                    selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list.extend(new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list)
                else:
                    pass
                    #print "DEBUG: the assigned C-H pairs for aa type", previous_aatype, "were below the requested threshold (",args.ASSIGNMENT_CUTOFF * TRI_num,")"
                    #print "DEBUG: new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list=", new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list
                #print "DEBUG: selected correct H-C resonpairs!"
                aatypeResonpairMatchesTuple_list = []
            
            aatypeResonpairMatchesTuple_list.append(group8)
            previous_aatype = aa_type
        #print "DEBUG: aa_type=",aa_type,"previous_aatype=",previous_aatype,"aatypeResonpairMatchesTuple_list=",aatypeResonpairMatchesTuple_list
        if args.PROB_PRODUCT == True:
            new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list, aatype2save, aatype_predictionScore = select_correct_H_C_resonpairs2(aatypeResonpairMatchesTuple_list)  # do the last aa_type
            aatype_predictionScore_dict[aatype2save] = aatype_predictionScore
        else:
            new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list = select_correct_H_C_resonpairs(aatypeResonpairMatchesTuple_list)  # do the last aa_type
        # CHECK IF ALL THE PEAKS HAVE BEEN ASSIGNED TO C-H TYPES FOR THE CURRENT AA TYPE. IF NOT DO NOT SAVE THIS AA TYPE AS POSSIBLE PREDICTION FOR THE CURRENT TOCSY index group
        if len(new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list) >= (args.ASSIGNMENT_CUTOFF * TRI_num):
            selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list.extend(new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list)
        else:
            pass
            #print "DEBUG: the assigned C-H pairs for aa type", previous_aatype, "were below the requested threshold (",args.ASSIGNMENT_CUTOFF * TRI_num,")"
            #print "DEBUG: new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list=", new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list
        #print "DEBUG: selected correct H-C resonpairs!"
    
    #print "DEBUG: selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list=", selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list
    aatype_occupancy_dict = {}  # dict with key the aa types and values how many times they were found during the H,C resonance analysis
    aatype_probsum_dict = {}    # dictionary with the aa type --> sum of probabilities of individual matched H,C resonance pairs
    aatype_CHnucleiType_presenceProbSum_mdict = tree()    # multidimensional dict: aa type --> C,H nuclei pair type --> probability of presence or absence (eq 6 in RESCUE2 paper)
    aatype_CHnucleiType_TresonIndex_mdict = tree()    # multidimensional dict: aa type --> C,H nuclei pair type --> TOCSY_reson_index
    aa_type_set = set()
    aatype_validC_H_pairsList_dict = {}
    #print "DEBUG: selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list=", selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list
    for group8 in selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list:
        #print "DEBUG: group8=",group8
        aa_type = group8[0]
        C_name = group8[1]
        H_name = group8[2]
        weighted_average_probability = group8[3]
        TOCSY_reson_index = group8[4]
        #TOCSY_reson_index = group5[4]  # we don't need this
        #if (aa_type, C_name, H_name) == previous_aatype_C_H_tuple:  # elegant way to keep only the (aatype, C_name, H_name) with the highest weighted_average_probability, e.g. we may have:
        #    # ('CYS', 'CA', 'HA', 0.012993371769822381), ('CYS', 'CA', 'HA', 0.002004997884084288), ('CYS', 'CA', 'HA', 0.0018243946693470656) but we keep only the 1st 
        #    continue
        try:
            if (C_name+"_"+H_name in aatype_validC_H_pairsList_dict[aa_type]):  # if we have saved twice this aatype,C,H combination print ERROR & exit
                print("ERROR: ",C_name+"_"+H_name,"already found !!!")
                sys.exit(1)
            aatype_validC_H_pairsList_dict[aa_type].append(C_name+"_"+H_name)
            aa_type_set.add(aa_type)
            aatype_occupancy_dict[aa_type] += 1
            ## THERE 2 DIFFENT WAYS TO CALCULATE THE OVERALL WEIGHTED AVERAGE PROBABILITY
            #aatype_probsum_dict[aa_type] += weighted_average_probability    # 1st WAY
            aatype_probsum_dict[aa_type] *= weighted_average_probability    # 2nd WAY
            aatype_CHnucleiType_presenceProbSum_mdict[aa_type][C_name+"_"+H_name] = (args.H_weight * Prob_CS[aa_type][H_name] + args.C_weight * Prob_CS[aa_type][C_name])/float((args.H_weight + args.C_weight))  # 2nd WAY
            aatype_CHnucleiType_TresonIndex_mdict[aa_type][C_name+"_"+H_name] = TOCSY_reson_index
        except KeyError:
            aatype_validC_H_pairsList_dict[aa_type] = [C_name+"_"+H_name]
            #if args.ALLOW_SINGLE_CH_PAIRS:
            aa_type_set.add(aa_type)
            aatype_occupancy_dict[aa_type] = 1
            aatype_probsum_dict[aa_type] = weighted_average_probability     # COMMON IN 1st & 2nd WAY
            aatype_CHnucleiType_presenceProbSum_mdict[aa_type][C_name+"_"+H_name] = (args.H_weight * Prob_CS[aa_type][H_name] + args.C_weight * Prob_CS[aa_type][C_name])/float((args.H_weight + args.C_weight))  # 2nd WAY
            aatype_CHnucleiType_TresonIndex_mdict[aa_type][C_name+"_"+H_name] = TOCSY_reson_index
            
        previous_aatype_C_H_tuple = (aa_type, C_name, H_name)
    
    matching_aaTypesProbTuple_list = []
    #print "DEBUG: Num_of_TOCSY_resonances=", Num_of_TOCSY_resonances
    #print "DEBUG: args.ASSIGNMENT_CUTOFF=", args.ASSIGNMENT_CUTOFF
    #print "DEBUG: aatype_occupancy_dict=",aatype_occupancy_dict
    #print "DEBUG: aatype_validC_H_pairsList_dict=",aatype_validC_H_pairsList_dict
    #print "DEBUG: aa_type_set=", aa_type_set
    for AA_type in aa_type_set:
        if aatype_occupancy_dict[AA_type] > aatype_maxH_C_pairs_dict[AA_type]:    # if we found more C-H pairs than the maximum number for that aa type, ommit it
            #print "DEBUG: ommiting AA_type=",AA_type,"due to higher occupancy (",aatype_occupancy_dict[AA_type],") than the maximum number of allowed C-H pairs (",aatype_maxH_C_pairs_dict[AA_type],")!"
            continue
        # consider this aa type as a match ONLY IF A PORTION OF THE TOCSY RESONANCES MATCH AND if ...
        # The latter excludes cases were we had 5 TOCSY resonances, we found 4 mathcing Val, but Val has only 4 C-H pairs in maximum
        # But when we have 1 TOCSY resonance, we found matching Cys, but Cys has 3 aatype_maxH_C_pairs_dict
        #print "DEBUG: AA_type=",AA_type,"aatype_occupancy_dict[AA_type]=", aatype_occupancy_dict[AA_type], ">= args.ASSIGNMENT_CUTOFF * Num_of_TOCSY_resonances=",args.ASSIGNMENT_CUTOFF * Num_of_TOCSY_resonances
        #print "DEBUG: aatype_maxH_C_pairs_dict[AA_type]=", aatype_maxH_C_pairs_dict[AA_type], ">= Num_of_TOCSY_resonances=", Num_of_TOCSY_resonances
        if (aatype_occupancy_dict[AA_type] >= (args.ASSIGNMENT_CUTOFF * Num_of_TOCSY_resonances)) and (aatype_maxH_C_pairs_dict[AA_type] >= Num_of_TOCSY_resonances):   
            try:
                #print "DEBUG: Saving AA_type",AA_type," as valid prediction."
                ## THERE ARE A FEW ALTERNATIVE SCORING SCHEMES HERE:
                if args.PROBABILITY_MODEL == 1:
                    if args.PROB_PRODUCT == True:
                        score = aatype_predictionScore_dict[AA_type]
                    else:
                        score = get_combined_weighted_and_presence_probability(AA_type, aatype_probsum_dict, aatype_CHnucleiType_presenceProbSum_mdict)
                elif args.PROBABILITY_MODEL == 2:
                    if args.PROB_PRODUCT == True:
                        score = aatype_predictionScore_dict[AA_type]
                    else:
                        score = aatype_probsum_dict[AA_type]
                #score = aatype_probsum_dict[AA_type]/float(Num_of_TOCSY_resonances) # divide by the total num of TOCSY resonance pairs; 1st WAY
                # CYS GET WORSE & Z-SCORES MORE NEGATIVE # score = aatype_probsum_dict[AA_type]/float(aatype_occupancy_dict[AA_type]) # divide by the prob sum by the occupancy of that aatype
                # BAD! # score = float(aatype_maxH_C_pairs_dict[AA_type]) * aatype_probsum_dict[AA_type]/(aatype_occupancy_dict[AA_type]**2)
                matching_aaTypesProbTuple_list.append( (AA_type, score) )
            except TypeError:
                print(aatype_probsum_dict[aa_type], Num_of_TOCSY_resonances)
                sys.exit(1)
    
    
    # FINALLY IF NO AA TYPE PREDICTIONS WERE MADE, REPEAT THE PROCEDURE INCLUDING ONLY CARBON RESONANCES
    if len(matching_aaTypesProbTuple_list) == 0:
        print("NO AA TYPE PREDICTIONS WERE MADE, THUS ACTIVATING CARBON-BASED PREDICTIONS!")
        possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list = consider_only_Carbon_matches(raw_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list, False, False)
        # The following sorting ensures that the alternative C-H pairs are sorted according to their weighted average probability 
        possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list.sort(key=itemgetter(0,1,2,3), reverse=True)  # sort by every column except the last 3
        ## The following block keeps only the matches for a particular aatype
        previous_aatype = None
        aatypeResonpairMatchesTuple_list = []
        selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list = []
        for group8 in possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list:
            #print "DEBUG: group8=", group8
            aa_type = group8[0]
            weighted_average_probability = group8[3]    # = [wH * hist(H) + wC * hist(C)]/2
            if aa_type != previous_aatype and previous_aatype != None:
                #print "DEBUG: aa_type=",aa_type,"previous_aatype=",previous_aatype,"aatypeResonpairMatchesTuple_list=",aatypeResonpairMatchesTuple_list
                if args.PROB_PRODUCT == True:
                    # selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list.extend( select_correct_H_C_resonpairs2(aatypeResonpairMatchesTuple_list) )
                    new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list, aatype2save, aatype_predictionScore = select_correct_H_C_resonpairs2(aatypeResonpairMatchesTuple_list)  # do the last aa_type
                    selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list.extend(new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list)
                    aatype_predictionScore_dict[aatype2save] = aatype_predictionScore
                else:
                    selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list.extend( select_correct_H_C_resonpairs(aatypeResonpairMatchesTuple_list) )
                aatypeResonpairMatchesTuple_list = []
            
            aatypeResonpairMatchesTuple_list.append(group8)
            previous_aatype = aa_type
        #print "DEBUG: aa_type=",aa_type,"previous_aatype=",previous_aatype,"aatypeResonpairMatchesTuple_list=",aatypeResonpairMatchesTuple_list
        if args.PROB_PRODUCT == True:
            # selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list.extend( select_correct_H_C_resonpairs2(aatypeResonpairMatchesTuple_list) )  # do the last aa_type
            new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list, aatype2save, aatype_predictionScore = select_correct_H_C_resonpairs2(aatypeResonpairMatchesTuple_list)  # do the last aa_type
            selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list.extend(new_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list)
            aatype_predictionScore_dict[aatype2save] = aatype_predictionScore
        else:
            selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list.extend( select_correct_H_C_resonpairs(aatypeResonpairMatchesTuple_list) )  # do the last aa_type
    
        #print "DEBUG: selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list=", selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list
        aatype_occupancy_dict = {}  # dict with key the aa types and values how many times they were found during the H,C resonance analysis
        aatype_probsum_dict = {}    # dictionary with the aa type --> sum of probabilities of individual matched H,C resonance pairs
        aatype_CHnucleiType_presenceProbSum_mdict = tree()    # multidimensional dict: aa type --> C,H nuclei pair type --> probability of presence or absence (eq 6 in RESCUE2 paper)
        aatype_CHnucleiType_TresonIndex_mdict = tree()    # multidimensional dict: aa type --> C,H nuclei pair type --> TOCSY_reson_index
        aa_type_set = set()
        aatype_validC_H_pairsList_dict = {}
        print("DEBUG: selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list=", selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list)
        for group8 in selected_aatype_prob_C_H_resonpair_TOCSYindex_tuple_list:
            #print "DEBUG: group8=",group8
            aa_type = group8[0]
            C_name = group8[1]
            H_name = group8[2]
            weighted_average_probability = group8[3]
            TOCSY_reson_index = group8[4]
            #if (aa_type, C_name, H_name) == previous_aatype_C_H_tuple:  # elegant way to keep only the (aatype, C_name, H_name) with the highest weighted_average_probability, e.g. we may have:
            #    # ('CYS', 'CA', 'HA', 0.012993371769822381), ('CYS', 'CA', 'HA', 0.002004997884084288), ('CYS', 'CA', 'HA', 0.0018243946693470656) but we keep only the 1st 
            #    continue
            try:
                if (C_name+"_"+H_name in aatype_validC_H_pairsList_dict[aa_type]):  # if we have saved twice this aatype,C,H combination print ERROR & exit
                    print("ERROR: ",C_name+"_"+H_name,"already found !!!")
                    sys.exit(1)
                aatype_validC_H_pairsList_dict[aa_type].append(C_name+"_"+H_name)
                aa_type_set.add(aa_type)
                aatype_occupancy_dict[aa_type] += 1
                ## THERE 2 DIFFENT WAYS TO CALCULATE THE OVERALL WEIGHTED AVERAGE PROBABILITY
                #aatype_probsum_dict[aa_type] += weighted_average_probability    # 1st WAY
                aatype_probsum_dict[aa_type] *= weighted_average_probability    # 2nd WAY
                aatype_CHnucleiType_presenceProbSum_mdict[aa_type][C_name+"_"+H_name] = (args.H_weight * Prob_CS[aa_type][H_name] + args.C_weight * Prob_CS[aa_type][C_name])/float((args.H_weight + args.C_weight))  # 2nd WAY
                aatype_CHnucleiType_TresonIndex_mdict[aa_type][C_name+"_"+H_name] = TOCSY_reson_index
            except KeyError:
                aatype_validC_H_pairsList_dict[aa_type] = [C_name+"_"+H_name]
                #if args.ALLOW_SINGLE_CH_PAIRS:
                aa_type_set.add(aa_type)
                aatype_occupancy_dict[aa_type] = 1
                aatype_probsum_dict[aa_type] = weighted_average_probability     # COMMON IN 1st & 2nd WAY
                aatype_CHnucleiType_presenceProbSum_mdict[aa_type][C_name+"_"+H_name] = (args.H_weight * Prob_CS[aa_type][H_name] + args.C_weight * Prob_CS[aa_type][C_name])/float((args.H_weight + args.C_weight))  # 2nd WAY
                aatype_CHnucleiType_TresonIndex_mdict[aa_type][C_name+"_"+H_name] = TOCSY_reson_index
                
            previous_aatype_C_H_tuple = (aa_type, C_name, H_name)
        
        matching_aaTypesProbTuple_list = []
        #print "DEBUG: Num_of_TOCSY_resonances=", Num_of_TOCSY_resonances
        #print "DEBUG: args.ASSIGNMENT_CUTOFF=", args.ASSIGNMENT_CUTOFF
        #print "DEBUG: aatype_occupancy_dict=",aatype_occupancy_dict
        #print "DEBUG: aatype_validC_H_pairsList_dict=",aatype_validC_H_pairsList_dict
        #print "DEBUG: aa_type_set=", aa_type_set
        for AA_type in aa_type_set:
            if aatype_occupancy_dict[AA_type] > aatype_maxH_C_pairs_dict[AA_type]:    # if we found more C-H pairs than the maximum number for that aa type, ommit it
                #print "DEBUG: ommiting AA_type=",AA_type,"due to higher occupancy (",aatype_occupancy_dict[AA_type],") than the maximum number of allowed C-H pairs (",aatype_maxH_C_pairs_dict[AA_type],")!"
                continue
            # consider this aa type as a match ONLY IF A PORTION OF THE TOCSY RESONANCES MATCH AND if ...
            # The latter excludes cases were we had 5 TOCSY resonances, we found 4 mathcing Val, but Val has only 4 C-H pairs in maximum
            # But when we have 1 TOCSY resonance, we found matching Cys, but Cys has 3 aatype_maxH_C_pairs_dict
            #print "DEBUG: AA_type=",AA_type,"aatype_occupancy_dict[AA_type]=", aatype_occupancy_dict[AA_type], ">= args.ASSIGNMENT_CUTOFF * Num_of_TOCSY_resonances=",args.ASSIGNMENT_CUTOFF * Num_of_TOCSY_resonances
            #print "DEBUG: aatype_maxH_C_pairs_dict[AA_type]=", aatype_maxH_C_pairs_dict[AA_type], ">= Num_of_TOCSY_resonances=", Num_of_TOCSY_resonances
            if (aatype_occupancy_dict[AA_type] >= (args.ASSIGNMENT_CUTOFF * Num_of_TOCSY_resonances)) and (aatype_maxH_C_pairs_dict[AA_type] >= Num_of_TOCSY_resonances):   
                try:
                    #print "DEBUG: Saving AA_type",AA_type," as valid prediction."
                    ## THERE ARE A FEW ALTERNATIVE SCORING SCHEMES HERE:
                    if args.PROBABILITY_MODEL == 1:
                        if args.PROB_PRODUCT == True:
                            score = aatype_predictionScore_dict[AA_type]
                        else:
                            score = get_combined_weighted_and_presence_probability(AA_type, aatype_probsum_dict, aatype_CHnucleiType_presenceProbSum_mdict)
                    elif args.PROBABILITY_MODEL == 2:
                        if args.PROB_PRODUCT == True:
                            score = aatype_predictionScore_dict[AA_type]
                        else:
                            score = aatype_probsum_dict[AA_type]
                    #score = aatype_probsum_dict[AA_type]/float(Num_of_TOCSY_resonances) # divide by the total num of TOCSY resonance pairs; 1st WAY
                    # CYS GET WORSE & Z-SCORES MORE NEGATIVE # score = aatype_probsum_dict[AA_type]/float(aatype_occupancy_dict[AA_type]) # divide by the prob sum by the occupancy of that aatype
                    # BAD! # score = float(aatype_maxH_C_pairs_dict[AA_type]) * aatype_probsum_dict[AA_type]/(aatype_occupancy_dict[AA_type]**2)
                    matching_aaTypesProbTuple_list.append( (AA_type, score) )
                except TypeError:
                    print(aatype_probsum_dict[aa_type], Num_of_TOCSY_resonances)
                    sys.exit(1)
    
    
    return matching_aaTypesProbTuple_list
########################################################## END OF FUNCTION DEFINITIONS #########################################################


################################################ PEPTIDE TREES FUNCTION DEFITION ##################################################

def populate_peptideTree_leaves(Peptide_Tree, iaaindex):
    """
        FUNCTION that adds new branches to the leaves of the Tree.
        
        ARGUMENTS:
        Peptide_Tree:   The Tree structure with connectivities
        iaaindex:       The AAIG of the current element of the chain
        RETURNS:
        (Peptide_Tree, BOOLEAN):    A tuple with elements the input Tree structure with new branches (if applicable), and a BOOLEAN value which is True if the function added
                                       new leaves to the Tree, or False otherwise
    """
    
    def calculate_Paaiminus1_TAAIGi(TAAIGi, aaiminus1):
        # aatype_P_dict:                                  probability P[aatype(i-1)], dict with the P[aatype(i-1)] for every aatype(i-1)

        aatype_P_dict = shared.getConst('AATYPE_P_DICT')
        iaaindex_iminus1aaTypesTransformedProbTupleList_dict = shared.getConst('IAAINDEX_IMINUS1AATYPESTRANSFORMEDPROBTUPLELIST_DICT')
        
        PTAAIGi = 0    # P[TAAIG(i)] = Sum{P[TAAIG(i)|aatype(i-1)]}
        for duplet in iaaindex_iminus1aaTypesTransformedProbTupleList_dict[TAAIGi]:
            #print "DEBUG: duplet=", duplet
            aatype = duplet[0]
            PTAAIGi += duplet[1]
            if aatype == aaiminus1:
                PTAAIGi_aaiminus1 = duplet[1]  # this is P[TAAIG(i)|aatype(i-1)]
        
        #print "DEBUG: aaiminus1=",aaiminus1
        #print "DEBUG: aatype_P_dict=", aatype_P_dict
        Paaiminus1 = aatype_P_dict[aaiminus1]
        #print "DEBUG: PTAAIGi_aaiminus1=",PTAAIGi_aaiminus1
        #print "DEBUG: Paaiminus1=",Paaiminus1
        #print "DEBUG: PTAAIGi=", PTAAIGi
        Paaiminus1_TAAIGi = PTAAIGi_aaiminus1 * Paaiminus1 / PTAAIGi
        #print "DEBUG: Paaiminus1_TAAIGi=", Paaiminus1_TAAIGi
        return Paaiminus1_TAAIGi
    
    iaaindex_iminus1aaTypesCutoffProbTupleList_dict = shared.getConst('IAAINDEX_IMINUS1AATYPESCUTOFFPROBTUPLELIST_DICT')
    aa3to1_dict = shared.getConst('AA3TO1_DICT')
    #iminus1aa_iaaindexTypesProbTupleList_dict = shared.getConst(iminus1aa_iaaindexTypesProbTupleList_dict=IMINUS1AA_IAAINDEXTYPESPROBTUPLELIST_DICT)
    #iaaindex_iminus1aaTypesPaatypeTupleList_dict = shared.getConst(iaaindex_iminus1aaTypesPaatypeTupleList_dict=IAAINDEX_IMINUS1AATYPESPAATYPETUPLELIST_DICT)
    aatype_P_dict = shared.getConst('AATYPE_P_DICT')
        
    number_of_new_leaves = 0
    for leaf in Peptide_Tree.get_leaves():
        if leaf.name == "P":    # proline can only be the root of the tree, not a leaf!
            continue
        leaf.add_feature("aaindex", iaaindex)
        try:
        #print "DEBUG: iaaindex_iminus1aaTypesCutoffProbTupleList_dict[iaaindex]=", iaaindex_iminus1aaTypesCutoffProbTupleList_dict[iaaindex]
            for child_tuple in iaaindex_iminus1aaTypesCutoffProbTupleList_dict[iaaindex]:
                #print "DEBUG: adding child_tuple=", child_tuple
                aa_type = child_tuple[0]    # i-1 aa type
                Paaiminus1_TAAIGi = calculate_Paaiminus1_TAAIGi(iaaindex, aa_type)
                new_child = leaf.add_child(name=aa3to1_dict[aa_type], dist=Paaiminus1_TAAIGi) # add a new brach to the current aa type (leaf) with length the respective P(TAAIGi) weighted by
                                                                                     # the probability of of the current transition in the current chain
                #print "DEBUG: adding connection: ",leaf.name,"-->",aa_type
                number_of_new_leaves += 1
        except KeyError:
            continue
    
    if number_of_new_leaves > 0:
        return (Peptide_Tree, True)
    else:
        return (Peptide_Tree, False)
    # print Peptide_Tree.get_ascii(show_internal=True, compact=False, attributes=["name", "dist"])
    # print Peptide_Tree.get_ascii(show_internal=True, compact=False)


def build_Peptide_Tree(chainIndex, chainScore):
    """
        FUNCTION to build and save into files all the possible peptides derived for a specified chain.
        
        ARGUMENTS:
        chainScore:     list containing to TOCSY index groups of the chain, and as the last element the frequency (score) of the chain
        chainIndex:     a unique integer number to be appended at the end of the file which will contain the peptides
    """
    # geometric progression (see my H-REX with PLUMED instructions)
    resonNumWeight_list = [1.000000, 2.656854, 3.928203, 5.000000, 5.944272, 6.797959, 7.583005]
    
    
    try:
        total_chain_number = shared.getConst('TOTAL_CHAIN_NUMBER')
        i_iminus1_dict = shared.getConst('I_IMINUS1_DICT')
        args = shared.getConst('ARGS')
        chain = list(chainScore[0:-1])
        score_of_chain = chainScore[-1]
        print("Predicting possible peptide sequences from chain ("+str(chainIndex+1)+"/"+str(total_chain_number)+"):  ",chain,"...")
        
        peptideScoreChainindex_list = []
        peptideProbList_list = []
        
        Peptide_Tree = Tree()
        Root = Peptide_Tree.get_tree_root()
        iaaindex = chain[-1]    # we start from the C-term of the chains and move twoards the N-term
        Root.add_features(name=iaaindex, aaindex=iaaindex)  # only the root of the Peptide_Tree has the name of the AAIG, the rest of the nodes have the names of possible aa types in the
                                            # peptide sequence starting from iaaindex
        
        if len(chain) < args.MIN_PEPTIDE_LENGTH: # if the chain contains fewer aa indices than the specified minimum, do not create the Peptide Tree
            print("Chain ", chain, " has length smaller that the cutoff", args.MIN_PEPTIDE_LENGTH, ". No peptides will be built.")
            return
        sys.stdout.write("Adding branches corresponding to amino acid index ")
        
        for iaaindex in reversed(chain):
            sys.stdout.write(str(iaaindex)+" ")
            sys.stdout.flush()
            # Print the Tree
            #print Peptide_Tree.get_ascii(show_internal=True, compact=False, attributes=["name", "dist"])
            Peptide_Tree, expand_tree = populate_peptideTree_leaves(Peptide_Tree, iaaindex)
            if expand_tree == False:    # if no leaves were added in the last iteration, do not expand the Peptide_Tree anymore
                break
        print("")    # change line after sys.stdout.write() invokation
        
        
        # print "DEBUG: printing peptide tree:"
        # print Peptide_Tree.get_ascii(show_internal=True, compact=False, attributes=["name", "dist", "aaindex"])
        for leaf in Peptide_Tree.get_leaves(): # iterate over the bottom nodes (leafs) of the Peptide_Tree, namely the N-teminal ends of the peptides
            peptide = []
            peptide_probability_list = []   # list containing the prediction scores of each aa of the peptide, in the direction N-term -> C-term
            total_peptide_probability = leaf.dist
            peptide.append(leaf.name)
            peptide_probability_list.append(leaf.dist * score_of_chain) # save transformed distance (probability that this aa type is correctly predicted) x the frequency (score) of the source chain
            seq_length = 0
            for ancestor in leaf.get_ancestors():   # multiply the individual probabilities of each prediction to get the total ZscoreRatio that this peptide sequence is correct
                seq_length += 1
                peptide.append(ancestor.name)
                peptide_probability_list.append(ancestor.dist * score_of_chain) # save transformed distance (probability that this aa type is correctly predicted) x the frequency (score) of the source chain
                total_peptide_probability *= ancestor.dist
                del ancestor
            if seq_length >= args.MIN_PEPTIDE_LENGTH:   # don't save peptides shorter than the minimum length
                #print "DEBUG: appending peptide:",peptide
                peptide.append(total_peptide_probability * score_of_chain)   # finally multiply the total peptide probability with the score of the source chain
                peptide.append(chainIndex)  # append the chain index
                #print "DEBUG: saving peptide=",peptide
                #print "DEBUG: saving peptide_probability_list=", peptide_probability_list
                peptideScoreChainindex_list.append(tuple(peptide))
                
                ##~~~~~~~~~~~~~~~~~~~~~~~~ WEIGHTS OF PREDICTION SCORES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                TAAIG_i = chain[-1]   # we move from the leaf to the root, namely from the N-term to the C-term of the chain, so this is the i TOCSY group
                length = len(peptide_probability_list)
                #print "DEBUG: i_iminus1_dict[TAAIG_i]=", i_iminus1_dict[TAAIG_i]
                for j, TAAIG_iminus1 in enumerate(reversed(chain[:-1])):
                    #print "DEBUG: TAAIG_iminus1=", TAAIG_iminus1
                    (TAAIG_iminus1, occupancy, tot_reson_num) = [t for t in i_iminus1_dict[TAAIG_i] if t[0]==TAAIG_iminus1][0]  # retrieve the connectivity info (occupancy and total num of resonances)
                    peptide_probability_list[length-(j+2)] *= (occupancy/float(tot_reson_num)) * resonNumWeight_list[tot_reson_num-1] * occupancy  # -1 because counting starts from 0 and another -1 because the C-term aa (root of the tree) does not have a connectivity
                    TAAIG_i = TAAIG_iminus1
                ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                
                peptideProbList_list.append(peptide_probability_list)
            del leaf
            del peptide
        del Peptide_Tree
        # save a pickle object of the form (peptideScoreChainindex_list, peptideProbList_list)
        if len(peptideScoreChainindex_list) == 0 and len(peptideProbList_list) == 0:
            print("WARNING: saving empty peptide file tmp_peptide_folder/chainIndex_"+str(chainIndex)+'.pickle.bz2')
        with bz2.BZ2File('tmp_peptide_folder/chainIndex_'+str(chainIndex)+'.pickle.bz2', 'wb') as f:
            pickle.dump((peptideScoreChainindex_list, peptideProbList_list), f)
        del peptideScoreChainindex_list, peptideProbList_list
        gc.collect()
    
    except:
        type, value, tb = sys.exc_info()
        lines = traceback.format_exception(type, value, tb)
        print((''.join(lines)))
        raise

################################################ END OF FUNCTION DEFITION ##################################################

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##              CALCULATE BAYESIAN STATISTICS                        #
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

def transform(probability_array):
    """
        FUNCTION to do a mathematical transform on a value
    """
    if args.TRANSFORM_TYPE == "None":
        return probability_array
    elif args.TRANSFORM_TYPE == "log":
        return np.log(probability_array)
    elif args.TRANSFORM_TYPE == "log10":
        return np.log10(probability_array)
    elif args.TRANSFORM_TYPE == "boxcox_pearson":
        lmax_pearsonr = stats.boxcox_normmax(probability_array)
        prob_pearson = stats.boxcox(probability_array, lmbda=lmax_pearsonr)
        return prob_pearson
    elif args.TRANSFORM_TYPE == "boxcox_mle":
        prob_mle, lmax_mle = stats.boxcox(probability_array)
        return prob_mle


def chunkIt(seq, num):
    """
        split a list into a specified number of approximately equal sublists.
    """
    avg = len(seq) / float(num)
    out = []
    last = 0.0
  
    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg
  
    return out


def remove_redundant_chains(chainScore_list):
    
    progbar = shared.getConst('PROGBAR')
    all_chainScore_list = shared.getConst('ALL_CHAINSCORE_LIST')
    chainScores2remove_set = set()
    chain_num = len(chainScore_list)
    for i, chainScore1 in enumerate(chainScore_list):
        chain1 = chainScore1[0:-1]
        score1 = chainScore1[-1]
        progbar.set_progress((float(i)/chain_num))
        for chainScore2 in all_chainScore_list:
            chain2 = chainScore2[0:-1]
            score2 = chainScore2[-1]
            if Chain.is_sublist(chain1, chain2) and score1 == score2: # this is the correct one
                chainScores2remove_set.add(tuple(chainScore1))   # make first the list immutable to be able to insert it into a set
    
    return chainScores2remove_set


def save_peptides_to_file(peptide_files_list, part):
    
    progbar = shared.getConst('PROGBAR')
    args = shared.getConst('ARGS')
    all_chainScore_list = shared.getConst('NEW_ALL_CHAINSCORE_LIST')
    f = open("peptides."+str(args.MAX_PEPTIDE_LENGTH)+"mers.list.chunk"+part, 'w')
    index = 0
    peptide_file_index = 0
    peptide_file_num = len(peptide_files_list)
    for peptide_file in peptide_files_list:
        mo = re.search('chainIndex_([0-9]+).pickle.bz2', peptide_file)
        if mo:
            peptide_file_index += 1
            progbar.set_progress((float(peptide_file_index)/peptide_file_num))
            chainIndex = int(mo.group(1))
            #print "DEBUG: loading tmp_peptide_folder/" + peptide_file
            try:
                with bz2.BZ2File('tmp_peptide_folder/' + peptide_file, "rb") as pickled_file:
                    pickle_object = pickle.load(pickled_file)
            except EOFError:
                print("ERROR: file tmp_peptide_folder/" + peptide_file, " seems to be corrupted perhaps due to abrupt termination of the program before the file was properly saved.\
                Please run again the program without the -resume option.")
                sys.exit(1)
            peptideScoreChainindexList_list = pickle_object[0]  # a list of lists containing the possible peptide sequences and the overall score as the last element
            peptideProbList_list = pickle_object[1] # a list of lists containing the probability that each amino acid of the peptide is correctly predicted from the chain
            for peptideScoreChainindexList, peptideProbList in zip(peptideScoreChainindexList_list, peptideProbList_list):
                chain = all_chainScore_list[chainIndex]
                peptide = peptideScoreChainindexList[0:-3]
                Cterm = peptideScoreChainindexList[-3]
                peptide_score = peptideScoreChainindexList[-2]
                chainIndex2 = peptideScoreChainindexList[-1]
                if chainIndex != chainIndex2:
                    print("ERROR: the chainIndex of file ", peptide_file, " should be ", chainIndex2, ", not ", chainIndex)
                    sys.exit(1)
                if peptide_score == 0.0:    # do not save peptides containing aa types that are not present in the protein
                    continue
                #print "DEBUG: chain=", chain, " peptide = ", peptide, "index=", index, "peptideProbList=", peptideProbList
                f.write("Chain probability = "+str(chain[-1])+" peptide P_"+str(len(peptide))+"_"+str(index+1) + " = "+''.join(peptide)+" "+'-'.join(chain[0:-1])+" "+','.join(map(str, peptideProbList[:-1]))+"\n")
                index += 1
    f.close()


def group_carbons(possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list):
    """
        ARGUMENTS:
        possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list:    list of lists of the form (aa type, Carbon name, Hydrogen name, weighted average probability, TOCSY_reson_index, H resonance, C resonance, carbon group)
        RETURNS:
        The same as the input but with the correct carbon groups (last element)
    """
    global aa_carbonBondedGeminalHydrogensDict_dict
    
    #print "DEBUG group_carbons: possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list=", possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list
    aa_types_set = set([x[0] for x in possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list])
    updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list = []    # same as possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list but with the cluster IDs for each aa type
    for aa_type in aa_types_set:
        singleAAtype_prob_C_H_resonpair_TOCSYindex_list_list = [x for x in possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list if x[0]==aa_type]
        if not aa_type in list(aa_carbonBondedGeminalHydrogensDict_dict.keys()):  # this aa does not have geminal protons, skip clustering
            clustID = 0
            for index in range(len(singleAAtype_prob_C_H_resonpair_TOCSYindex_list_list)):
                clustID += 1
                singleAAtype_prob_C_H_resonpair_TOCSYindex_list_list[index][7] = clustID
            # save the updated list entries with the correct cluster IDs
            updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list.extend(singleAAtype_prob_C_H_resonpair_TOCSYindex_list_list)
            continue    # continue with the next possible aa type
        
        new_aaindex_carbonGroupsTuple_dict = {}
        #print "DEBUG: before singleAAtype_prob_C_H_resonpair_TOCSYindex_list_list=", singleAAtype_prob_C_H_resonpair_TOCSYindex_list_list
        carbonCS_list = [l[6] for l in singleAAtype_prob_C_H_resonpair_TOCSYindex_list_list]
        #print "DEBUG: carbonCS_list=", carbonCS_list
        cl = HierarchicalClustering(carbonCS_list , lambda x,y: float(abs(x-y)))
        
        CONTINUE_CLUSTERING = True
        cutoff = 0.300001       # <== CHANGE ME (arbitrary starting value)
        while CONTINUE_CLUSTERING:  # continue clustering by lowering the cutoff until no more than 2 carbons resonances are within each cluster
            CONTINUE_CLUSTERING = False
            cluster_list = cl.getlevel(cutoff)    # <== CHANGE ME (the tolerance for carbon)
            #print "DEBUG: cluster_list=", cluster_list
            for cluster in cluster_list:
                try:
                    if len(set(cluster)) > 2:
                        CONTINUE_CLUSTERING = True
                        cutoff -= 0.02
                        break
                except TypeError:   # in case only one nucleus prediction is available (and only one peak)
                    CONTINUE_CLUSTERING = False
            
        for index in range(len(singleAAtype_prob_C_H_resonpair_TOCSYindex_list_list)):
            carbon_resonance = singleAAtype_prob_C_H_resonpair_TOCSYindex_list_list[index][6]
            for clustID, clust in enumerate(cluster_list):
                #print "DEBUG: clustID=", clustID, "clust=", clust
                if len(cluster_list) == 1 and carbon_resonance == clust:    # to avoid "TypeError: argument of type 'float' is not iterable"
                    singleAAtype_prob_C_H_resonpair_TOCSYindex_list_list[index][7] = clustID + 1  # because enumeration starts from 0
                elif carbon_resonance in clust:
                    singleAAtype_prob_C_H_resonpair_TOCSYindex_list_list[index][7] = clustID + 1  # because enumeration starts from 0
        # save the updated list entries with the correct cluster IDs
        updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list.extend(singleAAtype_prob_C_H_resonpair_TOCSYindex_list_list)
        #print "DEBUG group_carbons: appended carbon clusters with ", singleAAtype_prob_C_H_resonpair_TOCSYindex_list_list
    
    #print "DEBUG group_carbons: updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list=", updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list
    return updated_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list


if __name__ == "__main__":
    
    try:
        
        args = cmdlineparse()
        if args.H_weight <= 0 or args.H_weight > 1:
            print("ERROR: -wH must be a number greater than 0 and lower or equal to 1!")
            sys.exit(1)
        if args.C_weight <= 0 or args.C_weight > 1:
            print("ERROR: -wC must be a number greater than 0 and lower or equal to 1!")
            sys.exit(1)
        if (args.ASSIGNMENT_CUTOFF != None and args.POOL_AA_TYPES_FILE != None) or (args.ASSIGNMENT_CUTOFF != None and args.COMPLETE_AA_TYPES_FILE != None):
            print("ERROR: Usage of -poolaafile or -allaafile is incompatible with -acutoff argument! Either remove -poolaafile and -allaafile to make new amino acid \
            type predictions from scratch, or remove -acutoff to use the amino acid prediction files you provided as input.")
            sys.exit(1)
        if args.ASSIGNMENT_CUTOFF == None and args.POOL_AA_TYPES_FILE == None and args.COMPLETE_AA_TYPES_FILE == None:
            print("ERROR: you must either provide an aa type assignment cutoff with -acutoff argument (recommented value: 1.0), or provide a pool amino acid assignment file with \
            argument -poolaafile and file with all the aa type assignment with -allaafile argument.")
            sys.exit(1)
        
        # Launch input file check    
        check_input_files_format()
        
        #if args.EXAMPLES:
        #    print "USAGE: python HCTOCSYNH_to_HCNOENH_tree-based.py -tocsy <HCTOCSYNH file> -noesy <HCNOENH file> [-stdH <stdev1>] [-stdC <stdev2>] [-update]"
        #    print "Example: ./4D_assignment_tree-based.py -hsqc TDhsqc.list -tocsy TDtocsy.list -noesy TDnoesy.list -stdH 0.03 -stdC 0.3 -rstdH 0.02 -rstdC 0.2 -cutoff 1.0"
        #    sys.exit(0)
        if args.DOWNLOAD_CS_HISTOGRAMS:
            download_CS_histograms()
        
        
        # GRADUALLY REDUCE TOLERANCES TO FIND MORE NON-OVERLAPPING GROUPS IN HSQC SPECTRUM
        with open(args.HSQC_FILE, 'r') as f:
            root_contents=f.readlines()
        remaining_root_contents = []    # list of Root spectrum lines with overlapping groups
        nonoverlapping_root_contents_and_tolerances = []    # list wiht the lines of the Root spectrum that contain non-overlapping groups and the associated rtoH, rtolN
        for rtolH, rtolN in zip([0.02, 0.01, 0.01, 0.005], [0.2, 0.1, 0.05, 0.015]):
        #for rtolH, rtolN in zip([0.02, 0.01, 0.01], [0.2, 0.1, 0.05]):
            #print "DEBUG: len(root_contents)=", len(root_contents)
            # FIND NON-OVERLAPPING RESONANCE GROUPS IN THE HSQC SPECTRUM
            nonoverlapping_root_contents_and_tolerances.extend(find_nonoverlapping_groups_in_root(root_contents, rtolH, rtolN))
            for triplet in nonoverlapping_root_contents_and_tolerances:
                line = triplet[0]
                try:
                    root_contents.remove(line)
                except ValueError:
                    continue
        
        #sys.exit(1)
        ## AT THE END GROUP WHICHEVER PEAKS IS LEFT TO THE GROUP IN THE HSQC SPECTRUCT THAT HAS THE NEAREST H,N RESONANCE DISTANCE
        #nonoverlapping_root_contents_and_tolerances.extend(find_nonoverlapping_groups_in_root(root_contents, rtolH, rtolN, True))
        #for triplet in nonoverlapping_root_contents_and_tolerances:
        #    line = triplet[0]
        #    try:
        #        root_contents.remove(line)
        #    except ValueError:
        #        continue
        
        #print "DEBUG: len(nonoverlapping_root_contents_and_tolerances)=", len(nonoverlapping_root_contents_and_tolerances)
        #for triplet in nonoverlapping_root_contents_and_tolerances:
        #    print triplet
        
        
        ##
        ## COPY AAIGs FROM HSQC SPECTRUM
        ##
        ## TOCSY_lines contains the AAIG of i, the H,C resonances of residue i-1 and the N,HN resonances of residue i
        print("DEBUG: final nonoverlapping_root_contents_and_tolerances=", nonoverlapping_root_contents_and_tolerances)
        TOCSY_lines, TOCSY_sidechain_resonances_list = copy_aaindices_from_root_spectrum_2(nonoverlapping_root_contents_and_tolerances, args.TOCSY_fname, "TOCSY")
        NOESY_lines = copy_aaindices_from_root_spectrum_2(nonoverlapping_root_contents_and_tolerances, args.NOESY_FILE, "NOESY", TOCSY_sidechain_resonances_list)
        
        ## Load TOCSY file
        TOCSY_contents = []
        print("DEBUG: TOCSY_lines=", TOCSY_lines)
        for TOCSY_line in TOCSY_lines:
            TOCSY_words_list = [a for a in re.split('\s+', TOCSY_line)[:-1] if a!= ''] # discard the last element which is ''
            try:
                if len(TOCSY_words_list)==5 and float(TOCSY_words_list[1]) and float(TOCSY_words_list[2]) and float(TOCSY_words_list[3]) and float(TOCSY_words_list[4]):
                    TOCSY_contents.append(TOCSY_words_list)
                # otherwise check if the file contains intensities in the last column
                elif len(TOCSY_words_list)==6 and float(TOCSY_words_list[1]) and float(TOCSY_words_list[2]) and float(TOCSY_words_list[3]) and float(TOCSY_words_list[4]) and float(TOCSY_words_list[5]):
                    TOCSY_contents.append(TOCSY_words_list)
            except ValueError:
                print("Discarding the following line from TOCSY:")
                TOCSY_line
                continue
        TOCSY_contents.sort(key=itemgetter(0))  # sort TOCSY lines by the random AAIG (2nd column)
        
        #print "DEBUG: TOCSY_contents=",TOCSY_contents
        #print "DEBUG: NOESY_lines =",NOESY_lines
        # DEPRECATED
        #with open(args.NOESY_FILE, 'r') as f:
        #    NOESY_lines=f.readlines()
        
        ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
        ##                                          1st ROUND                                         ##
        ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
        
        ##
        ## GET POSSIBLE CONNECTIVITIES AND PROCESS THEM
        ##
        
        def create_connectivity_dictionaries():
            connectivities_mdict = get_possible_connectivities(TOCSY_contents,NOESY_lines)
            # connectivities_mdict is a multidimensional dictionary with structure:
            # TOCSY AAIG => NOESY AAIG => [number of matched TOCSY w2,w3 resonances in NOESY for that NOESY AAIG, total number of TOCSY w2,w3 resonances for that TOCSY AAIG]
            i_iminus1_dict = {} # a dictionary containing all possible connectivities of every TOCSY AAIG; has keys the TOCSY aa indices and values lists of triplets (tuples) \
                                # consisting of (NOESYaaindex, occupancy, numOfResonances)
            i_iminus1_complete_dict = {} # the same dictionary but with all possible connectivities, including those below args.RESONANCE_MATCH_CUTOFF
            for TOCSYaaindex in list(connectivities_mdict.keys()):
                NOESYaaindex_list = list(connectivities_mdict[TOCSYaaindex].keys())
                if TOCSYaaindex in NOESYaaindex_list:
                    NOESYaaindex_list.remove(TOCSYaaindex)  # remove residue i from the list of possible residues (i-1)
                # Put the NOESYaaindex, occupancy and total number of resonances for that particular NOESYaaindex into a tuple for convenience
                NOESYaaindex_occupancy_numOfResonances_tuple__list = [] 
                for NOESYaaindex in NOESYaaindex_list:
                    occupancy = connectivities_mdict[TOCSYaaindex][NOESYaaindex][0]
                    numOfResonances = connectivities_mdict[TOCSYaaindex][NOESYaaindex][1]
                    NOESYaaindex_occupancy_numOfResonances_tuple__list.append((NOESYaaindex, occupancy, numOfResonances))
                sorted_NOESYaaindex_occupancy_numOfResonances_tuple__list = sorted(NOESYaaindex_occupancy_numOfResonances_tuple__list, key=itemgetter(1, 2), reverse=True)
                #print "DEBUG:",numOfResonances, TOCSYaaindex, sorted_NOESYaaindex_occupancy_numOfResonances_tuple__list
                # Iterate over all NOESY aa indices and, if applicable, keep those that match in all resonance pairs w2,w3
                for triplet in sorted_NOESYaaindex_occupancy_numOfResonances_tuple__list:
                    if float(triplet[1])/numOfResonances >= float(args.RESONANCE_MATCH_CUTOFF):     # adjustable cutoff
                        try:
                            i_iminus1_dict[TOCSYaaindex].append(triplet)
                        except KeyError:
                            i_iminus1_dict[TOCSYaaindex] = [(triplet)]
                    try:
                        i_iminus1_complete_dict[TOCSYaaindex].append(triplet)
                    except KeyError:
                        i_iminus1_complete_dict[TOCSYaaindex] = [(triplet)]
            
            # FILTER i_iminus1_dict ACCORDING TO THE Z-SCORE CUTOFF
            new_i_iminus1_dict = {}
            for i in list(i_iminus1_dict.keys()):
                prob_list = [float(triplet[1])/triplet[2] for triplet in i_iminus1_dict[i]]
                if len(prob_list) == 1:     
                    zscore_array = np.array([10.00000])     # default Z-score in the case of just one connectivity
                elif len(set(prob_list)) == 1:  # if all connectivities have the same value keep them all
                    #print "DEBUG: all connectivities of ", i," have the same value:", prob_list
                    zscore_array = np.array([10.00000]*len(prob_list))
                    #print "DEBUG: zscore_array=", zscore_array
                elif len(prob_list) > 2:
                    zscore_array = zscore(prob_list)  # convert probabilities to Z-scores
                elif len(prob_list) == 2:   # if only 2 connectivities exist, set the 2nd to 0.000 to allow you to use high -zcutoff values (<= 0)
                    zscore_array = np.array([1.00000, 0.00000])
                elif len(prob_list) == 0:   # if no connectivity was made create an empty array
                    zscore_array = np.array([])
                
                new_i_iminus1_dict[i] = []  # add i to the new dictionary
                for triplet, Zscore in zip(i_iminus1_dict[i], zscore_array):
                    ## ONLY IF THE Z-SCORE OF THE CONNECTIVITY IS GREATER THAN THE CUTOFF
                    if Zscore >= args.ZSCORE_RESONANCE_MATCH_CUTOFF or approx_equal(Zscore, args.ZSCORE_RESONANCE_MATCH_CUTOFF):
                        new_i_iminus1_dict[i].append(triplet)
            
            del i_iminus1_dict
            i_iminus1_dict = new_i_iminus1_dict
            del new_i_iminus1_dict
            
            return i_iminus1_dict, i_iminus1_complete_dict
        
        
        if not args.POOL_CONNECTIVITIES_FILE:
            print("Calculating Connectivities...")
            i_iminus1_dict, i_iminus1_complete_dict = create_connectivity_dictionaries()
            
        else:
            print("Loading Connectivities from file", args.POOL_CONNECTIVITIES_FILE)
            i_iminus1_pool_dict = {} # same dictionary with the pool of possible connectivities remained after filtering using absolute consensus matches of a previous run
            i_iminus1_complete_dict = {} # the same dictionary but with all possible connectivities, including those below args.RESONANCE_MATCH_CUTOFF
        
            ## I THINK THIS IS NOT NECESSARY ##
            #with open(args.POOL_CONNECTIVITIES_FILE, 'r') as f:
            #    connectivities_file_contents = f.readlines()
            #for line in connectivities_file_contents[1:]:
            #    word_list = line.split()
            #    TOCSY_aaindex = word_list[0]
            #    i_iminus1_dict[TOCSY_aaindex] = []
            #    values_string = ''.join(word_list[1:])
            #    elements_string = re.sub('[\(\)\'\"]', '',  values_string ).split("),(")[0]
            #    elements_list = elements_string.split(",")
            #    for aa, occupancy, TOCSY_resonnum in zip(elements_list[0::3], elements_list[1::3], elements_list[2::3]):
            #        #print aa, occupancy, TOCSY_resonnum
            #        i_iminus1_dict[TOCSY_aaindex].append((aa, int(occupancy), int(TOCSY_resonnum)))
            ###################################
            
            # Now load the pool connectivities file to calculate correct probabilities
            with open(args.POOL_CONNECTIVITIES_FILE, 'r') as f:
                pool_connectivities_file_contents = f.readlines()
            for line in pool_connectivities_file_contents[1:]:
                word_list = line.split()
                TOCSY_aaindex = word_list[0]
                i_iminus1_pool_dict[TOCSY_aaindex] = []
                values_string = ''.join(word_list[1:])
                elements_string = re.sub('[\(\)\'\"]', '',  values_string ).split("),(")[0]
                elements_list = elements_string.split(",")
                for aa, occupancy, TOCSY_resonnum in zip(elements_list[0::3], elements_list[1::3], elements_list[2::3]):
                    #print aa, occupancy, TOCSY_resonnum
                    i_iminus1_pool_dict[TOCSY_aaindex].append((aa, int(occupancy), int(TOCSY_resonnum)))
            
            # Now load the complete connectivities file to calculate correct probabilities
            with open(args.COMPLETE_CONNECTIVITIES_FILE, 'r') as f:
                complete_connectivities_file_contents = f.readlines()
            for line in complete_connectivities_file_contents[1:]:
                word_list = line.split()
                TOCSY_aaindex = word_list[0]
                i_iminus1_complete_dict[TOCSY_aaindex] = []
                values_string = ''.join(word_list[1:])
                elements_string = re.sub('[\(\)\'\"]', '',  values_string ).split("),(")[0]
                elements_list = elements_string.split(",")
                for aa, occupancy, TOCSY_resonnum in zip(elements_list[0::3], elements_list[1::3], elements_list[2::3]):
                    #print aa, occupancy, TOCSY_resonnum
                    i_iminus1_complete_dict[TOCSY_aaindex].append((aa, int(occupancy), int(TOCSY_resonnum)))
            #for k,v in i_iminus1_dict.items():
            #        print k,v
            
            # IF THE USER HAS SPECIFIED TOLERANCES FOR MATCHING, USE THEM TO FIND NEW CONNECTIVITIES FOR THOSE AAIGs THAT WERE NOT CORRECTED
            if '-tolH' in sys.argv and '-tolC' in sys.argv:
                print("Calculating Connectivities using H tolerance", args.tolH, "and C torelance", args.tolC)
                new_i_iminus1_dict, new_i_iminus1_complete_dict = create_connectivity_dictionaries()
                for TOCSYaaindex in list(new_i_iminus1_complete_dict.keys()):
                    if TOCSYaaindex not in list(i_iminus1_complete_dict.keys()):  # if there was not connectivity for this TAAIG add it to the dictionaries               
                        i_iminus1_complete_dict[TOCSYaaindex] = new_i_iminus1_complete_dict[TOCSYaaindex]   # update the connectivities
                        i_iminus1_pool_dict[TOCSYaaindex] = new_i_iminus1_complete_dict[TOCSYaaindex]
                    elif len(new_i_iminus1_complete_dict[TOCSYaaindex]) > len(i_iminus1_complete_dict[TOCSYaaindex]): # if there were found extra connectivities using the new tolerances
                        if len(i_iminus1_pool_dict[TOCSYaaindex]) == 1 and len(i_iminus1_complete_dict[TOCSYaaindex]) > 1:
                            are_all_unique = True   # do all Tindices in the complete dict have a single connectivity in the pool
                            for tmp_triplet in i_iminus1_complete_dict[TOCSYaaindex]:
                                tmp_TOCSYaaindex = tmp_triplet[0]
                                if tmp_TOCSYaaindex in list(i_iminus1_pool_dict.keys()) and len(i_iminus1_pool_dict[tmp_TOCSYaaindex]) > 1:
                                    are_all_unique = False
                                    break
                            if are_all_unique == False: # this means that i_iminus1_pool_dict[TOCSYaaindex] has been modified and not cleaned from Tindices used elsewhere
                                continue    # leave the connectivity of this TOCSYaaindex unchanged
                            i_iminus1_complete_dict[TOCSYaaindex] = new_i_iminus1_complete_dict[TOCSYaaindex]   # update the connectivities
                            i_iminus1_pool_dict[TOCSYaaindex] = new_i_iminus1_complete_dict[TOCSYaaindex]
                        elif len(i_iminus1_pool_dict[TOCSYaaindex]) == 1 and len(i_iminus1_complete_dict[TOCSYaaindex]) == 1:
                            i_iminus1_complete_dict[TOCSYaaindex] = new_i_iminus1_complete_dict[TOCSYaaindex]   # update the connectivities
                            i_iminus1_pool_dict[TOCSYaaindex] = new_i_iminus1_complete_dict[TOCSYaaindex]
            
            # CALCULATING CONNECTIVITIES FROM THE POOL OF CONNECTIVITIES USING THE args.RESONANCE_MATCH_CUTOFF
            print("Finding Connectivities above the -mcutoff from the Pool of Connectivities...")
            # connectivities_mdict is a multidimensional dictionary with structure:
            # TOCSY AAIG => NOESY AAIG => [number of matched TOCSY w2,w3 resonances in NOESY for that NOESY AAIG, total number of TOCSY w2,w3 resonances for that TOCSY AAIG]
            i_iminus1_dict = {} # a dictionary containing all possible connectivities of every TOCSY AAIG above args.RESONANCE_MATCH_CUTOFF; has keys the TOCSY aa indices and values lists of triplets (tuples) \
                                # consisting of (NOESYaaindex, occupancy, numOfResonances)
            for TOCSYaaindex in list(i_iminus1_pool_dict.keys()):
                # Iterate over all NOESY aa indices and, if applicable, keep those that match in all resonance pairs w2,w3
                for triplet in i_iminus1_pool_dict[TOCSYaaindex]:
                    numOfResonances = triplet[2]
                    if float(triplet[1])/numOfResonances >= float(args.RESONANCE_MATCH_CUTOFF):     # adjustable cutoff
                        try:
                            i_iminus1_dict[TOCSYaaindex].append(triplet)
                        except KeyError:
                            i_iminus1_dict[TOCSYaaindex] = [(triplet)]
        
            # FILTER i_iminus1_dict ACCORDING TO THE Z-SCORE CUTOFF
            new_i_iminus1_dict = {}
            for i in list(i_iminus1_dict.keys()):
                prob_list = [float(triplet[1])/triplet[2] for triplet in i_iminus1_dict[i]]
                if len(prob_list) == 1:     
                    zscore_array = np.array([10.00000])     # default Z-score in the case of just one connectivity
                elif len(set(prob_list)) == 1:  # if all connectivities have the same value keep them all
                    zscore_array = np.array([10.00000]*len(prob_list))
                elif len(prob_list) > 2:
                    zscore_array = zscore(prob_list)  # convert probabilities to Z-scores
                elif len(prob_list) == 2:   # if only 2 connectivities exist, set the 2nd to 0.000 to allow you to use high -zcutoff values (<= 0)
                    zscore_array = np.array([1.00000, 0.00000])
                elif len(prob_list) == 0:   # if no connectivity was made create an empty array
                    zscore_array = np.array([])
                new_i_iminus1_dict[i] = []  # add i to the new dictionary
                for triplet, Zscore in zip(i_iminus1_dict[i], zscore_array):
                    ## ONLY IF THE Z-SCORE OF THE CONNECTIVITY IS GREATER THAN THE CUTOFF
                    if Zscore >= args.ZSCORE_RESONANCE_MATCH_CUTOFF or approx_equal(Zscore, args.ZSCORE_RESONANCE_MATCH_CUTOFF):
                        new_i_iminus1_dict[i].append(triplet)
            
            del i_iminus1_dict
            i_iminus1_dict = new_i_iminus1_dict
            del new_i_iminus1_dict
            
            
        ## Convert occupancies to probabilities
        i_iminus1_normProbabilities_dict = {}   # dict of the form TAAIG(i)->[(TAAIG(i-1), probability), (...), ...]; eventually will contain only the matches above the cutoff
        aaindex_weightSum_dict = {}
        for TOCSY_aaindex in list(i_iminus1_complete_dict.keys()):
            weight_sum = 0
            for triplet in i_iminus1_complete_dict[TOCSY_aaindex]:  # use all possible matches to calculate the weight_sum
                occupancy = triplet[1]
                TOCSY_resonnum = triplet[2]
                weight = float(occupancy)/TOCSY_resonnum   # occupancy/total number of TOCSY resonances
                weight_sum += weight
            aaindex_weightSum_dict[TOCSY_aaindex] = weight_sum
        
        for TOCSY_aaindex in list(i_iminus1_dict.keys()):     # use only the matches above the cutoff to save theirs probabilities
            #print "DEBUG: TOCSY_aaindex=", TOCSY_aaindex, "aaindex_weightSum_dict[TOCSY_aaindex] =",aaindex_weightSum_dict[TOCSY_aaindex]
            for triplet in i_iminus1_dict[TOCSY_aaindex]:
                aa = triplet[0]
                occupancy = triplet[1]
                TOCSY_resonnum = triplet[2]
                weight = float(occupancy)/TOCSY_resonnum   # occupancy/total number of TOCSY resonances
                duplet = (aa, weight/aaindex_weightSum_dict[TOCSY_aaindex])           # recover the probability by dividing by the weight_sum
                try:
                    i_iminus1_normProbabilities_dict[TOCSY_aaindex].append(duplet)
                except KeyError:
                    i_iminus1_normProbabilities_dict[TOCSY_aaindex] = [duplet]
        
         
        #for k,v in i_iminus1_normProbabilities_dict.items():
        #    print k,v
        #    prob_sum = 0
        #    for duplet in v:
        #        prob_sum += duplet[1]
        #    print prob_sum
        #
        #sys.exit(1)
        
        ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
        ##              FORM CHAINS FROM CONNECTIVITIES               ##
        ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
        
        if not args.NON_REDUNDANT_CHAINS_FILE:
            
            chain = Chain()
            with open("connectivities_cutoff_"+str(args.RESONANCE_MATCH_CUTOFF)+"_zmcutoff"+str(args.ZSCORE_RESONANCE_MATCH_CUTOFF), 'w') as f:
                f.write("i\tpossible i-1\n")
                for k,v in list(i_iminus1_dict.items()):
                    print(k,v)
                    f.write(k+"\t"+', '.join(str(e) for e in v)+"\n")
            
            with open("connectivities_all", 'w') as f:
                f.write("i\tpossible i-1\n")
                for k,v in list(i_iminus1_complete_dict.items()):
                    #print k,v
                    f.write(k+"\t"+', '.join(str(e) for e in v)+"\n")
            
            
            for i in list(i_iminus1_dict.keys()):
                chain.build_Chain_Tree(i, args.MAX_PEPTIDE_LENGTH, args.MIN_PEPTIDE_LENGTH)
            
            ##
            ## Write chains to files
            ##
            
            print("Saving all the chains to chains.list file ...")
            def write_chains_file(all_chainScore_set, outfname):
                """
                ARGUMENTS:
                all_chainScore_set: can be both a set or a list
                """
                all_chainScore_list = list(all_chainScore_set)
                all_chainScore_list.sort(key=itemgetter(-1), reverse=True)
                with open(outfname, 'w') as f:
                    for chainScore_list in all_chainScore_list:
                        chain = chainScore_list[0:-1]
                        score = chainScore_list[-1]
                        #print "Overall Score = ",score," chain = ", chain
                        f.write("Overall Score = "+str(score)+" chain = "+str(chain)+"\n")
            
            write_chains_file(all_chainScore_set, "chains.list")
            
            ##
            ##  Remove redudancy from the chain list
            ### Bash script to help in validating the redundancy removal
            #while read line;
            #    do chain=$(perl -pi -e "s/.* chain = \((.*)\)/\1/" <<< $line);
            #    score=$(awk '{print $4}' <<< $line);
            #    echo "---> redundant chain = $chain";
            #    grep "Score = $score chain = ($chain" chains.non-redundant.list | grep "$chain";
            #done < chains.redundant.list > redundant_lines.txt
            
            print("Removing redundant chains ...")
            all_chainScore_list = list(all_chainScore_set)  # convert the set to a list to keep track of the sequence that each score belongs to
            del all_chainScore_set  # delete it to save memory
            progbar = ProgressBar(100)
            all_chainScore_listList_list = chunkIt(all_chainScore_list, scoop.SIZE)
            try:
                shared.setConst(PROGBAR=progbar)
            except TypeError:
                pass
            shared.setConst(ALL_CHAINSCORE_LIST=all_chainScore_list)
            results = list(futures.map(remove_redundant_chains, all_chainScore_listList_list))
            # results is a list of sets with chains to be removed
            chainScores2remove_set = set()
            for s in results:
                chainScores2remove_set.union(s)
            
            # remove subchains to keep only the longest, unique, non-overlapping ones
            for chainScore in chainScores2remove_set:
                all_chainScore_list.remove(tuple(chainScore))
            
            write_chains_file(all_chainScore_list, "chains.non-redundant.list")
        
        else:   # otherwise load the specified non-redudant chains file
            all_chainScore_list = []
            with open(args.NON_REDUNDANT_CHAINS_FILE, 'r') as f:
                for line in f:
                    word_list = line.split()
                    chain_score = float(word_list[3])
                    values_string = ''.join(word_list[6:])
                    Tindices_string = re.sub('[\(\)\'\"]', '',  values_string ).split("),(")[0]
                    chain = Tindices_string.split(",")
                    chain.append(chain_score)   # append the chain score to the end of the chain list
                    all_chainScore_list.append(chain)
        
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        #   AUTOMATIC ASSIGNMENT OF AMINO ACID TYPES TO TOCSY RESONANCES    #
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        
        if not args.POOL_AA_TYPES_FILE and not args.COMPLETE_AA_TYPES_FILE:
            ## LOAD 1D BMRB HISTOGRAMS
            print("\nLoading BMRB chemical shift histograms...")
            aa_carbon_binDensityList_mdict = tree()   # multidict with amino acid name -> Carbon atom type -> [array of bin limits, array of probability density]
            aa_hydrogen_binDensityList_mdict = tree()   # multidict with amino acid name -> Hydrogen atom type -> [array of bin limits, array of probability density]
            fnames = os.listdir(HOME_DIR+"/BMRB_data/")
            fpattern = re.compile("[A-Z]{3}_[A-Z0-9]+_hist.txt$")
            hist_files_list = list(filter(fpattern.search, fnames))
            for hist_file in hist_files_list:
                aa = hist_file.split("_")[0]
                atom = hist_file.split("_")[1]
                #try:
                if atom in allowed_aa_atoms_dict[aa]:   # load histograms of the allowed atom types only
                    bin_list = []
                    density_list = []
                    with open(HOME_DIR+"/BMRB_data/"+hist_file, 'r') as f:
                        for line in f:
                            word_list = line.split()
                            bin_list.append(float(word_list[0]))
                            density_list.append(float(word_list[1]))
                            # print "DEBUG: saving bin ", word_list[0], "and density", word_list[1]
                    bin_array = np.array(bin_list)
                    density_array = np.array(density_list)/sum(density_list)
                    # use C & H resonances to assign aa type to the (i-1) residue 
                    if atom[0] == "C":
                        aa_carbon_binDensityList_mdict[aa][atom] = [bin_array, density_array]
                    elif atom[0] == "H":
                        aa_hydrogen_binDensityList_mdict[aa][atom] = [bin_array, density_array]
                #except KeyError:   # active this exception if use N, HN resonance to assign aa type to the ith residue
                #    if aa == "PRO": # Prolines cannot be detected in root spectrum (N-H HSQC) because they don't have "HN" atom
                #        continue
            
            
            if args.USE_2D_HISTOGRAMS == True:
                ## IF ASKED IN THE COMMAND LINE, LOAD 2D BMRB HISTOGRAMS, TOO
                print("\nLoading 2D BMRB chemical shift histograms...")
                aa_CHpair_binProbabilityList_mdict = tree()   # multidict with amino acid name -> Carbon atom type -> [array of bin limits, array of probability density]
                fnames = os.listdir(HOME_DIR+"/BMRB_data/")
                fpattern = re.compile("[A-Z]{3}_[A-Z0-9-]+_correlated_2Dhist.smoothed.txt$")
                hist_files_list = list(filter(fpattern.search, fnames))
                for hist_file in hist_files_list:
                    #print "Loading histogram", hist_file 
                    aa = hist_file.split("_")[0]
                    CH_pair = hist_file.split("_")[1]
                    #try:
                    if CH_pair in list(aa_CHpair_2Dhist_mdict[aa].keys()):   # load histograms of the allowed atom types only
                        x_bin_list = []     # for H
                        y_bin_list = []     # for C
                        probability_list = []
                        with open(HOME_DIR+"/BMRB_data/"+hist_file, 'r') as f:
                            for line in f:
                                word_list = line.split()
                                x_bin_list.append(float(word_list[0]))
                                y_bin_list.append(float(word_list[1]))
                                probability_list.append(float(word_list[2]))
                        x_bin_array = np.array(x_bin_list)
                        y_bin_array = np.array(y_bin_list)
                        probability_array = np.array(probability_list)
                        # use C & H resonances to assign aa type to the (i-1) residue
                        aa_CHpair_binProbabilityList_mdict[aa][CH_pair] = [x_bin_array, y_bin_array, probability_array]
                    #except KeyError:   # active this exception if use N, HN resonance to assign aa type to the ith residue
                    #    if aa == "PRO": # Prolines cannot be detected in root spectrum (N-H HSQC) because they don't have "HN" atom
                    #        continue
            
            
            iaaindex_iminus1aaTypesProbTupleList_dict = OrderedDict() # ordereddict with keys the AAIG of residue i and values lists of tuples of the form (residue i-1 matching aa type, average probability)
            Num_of_TOCSY_resonances = 0 # save here the number of TOCSY resonances for a particular AAIG
            previous_TOCSY_aaindex = None
            possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list = []
            carbon_groups_list = []  # list of tuples of the form (aaindex, H resonance, C resonance, N resonance, HN resonance, carbon group)
            aaindex_carbonGroupsTuple_dict = {} # dict of the form: aaindex -> list of the form (aaindex, H resonance, C resonance, N resonance, HN resonance, carbon group)
            counter = 0
            previous_aaindex = ""
            sorted_TOCSY_contents = sorted(TOCSY_contents, key=itemgetter(0))
            print("sorted_TOCSY_contents=", sorted_TOCSY_contents)
            #sys.exit(1)
            #sorted_TOCSY_contents= [
            #    ['A6', '4.148', '61.783', '131.690', '8.959'], ['A6', '1.829', '33.205', '131.701', '8.960'], ['A6', '1.038', '58.474', '131.692', '8.961'],
            #    ['A6', '0.756', '57.876', '131.676', '8.961']
            #    ]
            for TOCSY_words_list in sorted_TOCSY_contents:
                try:
                    # ignore 1st column
                    TOCSY_aaindex=TOCSY_words_list[0]    # residue i
                    print("DEBUG: TOCSY_aaindex=",TOCSY_aaindex,"previous_TOCSY_aaindex=",previous_TOCSY_aaindex)
                    if previous_TOCSY_aaindex != None and TOCSY_aaindex != previous_TOCSY_aaindex: # if we look at a different AAIG in TOCSY, print the matches and occupancies
                                                                                                   # of the previous TOCSY AAIG and clear matchingNOESYaaindex_occupancy_dict
                        
                        print("Assigning possible aa types to the aa upstream of AAIG",previous_TOCSY_aaindex)
                        #print "DEBUG: switched from AAIG", previous_TOCSY_aaindex, "to", TOCSY_aaindex
                        #print "DEBUG: possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list=",possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list
                        #print "DEBUG: Num_of_TOCSY_resonances=",Num_of_TOCSY_resonances
                        # GROUP THE CARBON RESONANCE OF EACH TAAIG
                        new_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list = group_carbons(possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list)
                        iaaindex_iminus1aaTypesProbTupleList_dict[previous_TOCSY_aaindex] = get_aatypes_from_all_H_C_resonpairs(new_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list, Num_of_TOCSY_resonances)
                        #print "DEBUG: previous_TOCSY_aaindex=",previous_TOCSY_aaindex,"iaaindex_iminus1aaTypesProbTupleList_dict[previous_TOCSY_aaindex]=",iaaindex_iminus1aaTypesProbTupleList_dict[previous_TOCSY_aaindex]
                        possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list = [] # list of tuples of the form (aa type, probability, H_resonance, C_resonance, TOCSY_reson_index) containing all matching sets of H,C resonances
                        aaindex_carbonGroupsTuple_dict[TOCSY_aaindex] = carbon_groups_list
                        #print "DEBUG: aaindex_carbonGroupsTuple_dict=", aaindex_carbonGroupsTuple_dict
                        carbon_groups_list = []  # list of tuples of the form (TOCSY_aaindex, H resonance, C resonance, N resonance, HN, resonance, carbon group)
                        
                        Num_of_TOCSY_resonances = 0
                        #counter += 1
                        #if counter == 10:
                        #    break
                    
                    TOCSY_H_resonance=float(TOCSY_words_list[1]) # aliphatic H resonance of residue i-1 
                    TOCSY_C_resonance=float(TOCSY_words_list[2]) # aliphatic C (Ca,Cb,Cc,Cg,Ce,etc.) resonance of residue i-1; this C is covalently bonded to the above H
                    TOCSY_N_resonance=float(TOCSY_words_list[3])
                    TOCSY_HN_resonance=float(TOCSY_words_list[4])
                    Num_of_TOCSY_resonances += 1
                    #print "DEBUG: TOCSY_C_resonance=",TOCSY_C_resonance,"TOCSY_H_resonance=",TOCSY_H_resonance, "Num_of_TOCSY_resonances", Num_of_TOCSY_resonances
                    # valid_matches_list: list of lists of the form [aa type, Carbon name, Hydrogen name, weighted average probability, TOCSY_reson_index, H resonance, C resonance, carbon group]
                    if args.USE_2D_HISTOGRAMS == True:
                        valid_matches_list = get_aatypes_from_H_C_resonpair_2Dhist(TOCSY_H_resonance, TOCSY_C_resonance, Num_of_TOCSY_resonances)
                    else:
                        valid_matches_list = get_aatypes_from_H_C_resonpair(TOCSY_H_resonance, TOCSY_C_resonance, Num_of_TOCSY_resonances)
                    #if TOCSY_aaindex == "Y130":
                    #print "DEBUG: valid_matches_list=",valid_matches_list
                    possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list.extend(valid_matches_list)
                    carbon_groups_list.append([TOCSY_aaindex, TOCSY_H_resonance, TOCSY_C_resonance, TOCSY_N_resonance, TOCSY_HN_resonance, None])
                
                    previous_TOCSY_aaindex = TOCSY_aaindex
                except (ValueError, IndexError):
                    print("WARNING: the 3rd and 4th elements of the following TOCSY file line are not numbers:")
                    print("TOCSY file line:", TOCSY_words_list)
                    continue
            # THIS IS FOR THE LAST TOCSY INDEX
            #if not counter == 10:
            print("Assigning possible aa types to the last aa, which is upstream of AAIG",previous_TOCSY_aaindex)
            #print "DEBUG: possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list=",possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list
            # GROUP THE CARBON RESONANCE OF EACH TAAIG
            new_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list = group_carbons(possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list)
            iaaindex_iminus1aaTypesProbTupleList_dict[previous_TOCSY_aaindex] = get_aatypes_from_all_H_C_resonpairs(new_possible_aatype_prob_C_H_resonpair_TOCSYindex_list_list, Num_of_TOCSY_resonances)
            #print "DEBUG: previous_TOCSY_aaindex=",previous_TOCSY_aaindex,"iaaindex_iminus1aaTypesProbTupleList_dict[previous_TOCSY_aaindex]=",iaaindex_iminus1aaTypesProbTupleList_dict[previous_TOCSY_aaindex]
            
            ##### NEW WAY: CALCULATE Z-SCORE FOR EACH TOCSY INDEX INDIVIDUAL #####
            ## Calculate Z-scores from aa type prediction probabilites and keep only those predictions above the cutoff
            iaaindex_iminus1aaTypesZscoreTupleList_dict = OrderedDict() # ordereddict with keys the AAIG of residue i and values lists of tuples of the form (residue i-1 matching aa type, Z-score)
            iaaindex_iminus1aaTypesCutoffProbTupleList_dict = OrderedDict() # contains only the tuples (residue i-1 matching aa type, probability) with Z-score above the cutoff
            for iaaindex, iminus1aaTypesProbTuple_list in list(iaaindex_iminus1aaTypesProbTupleList_dict.items()):
                # DECIDE THE PROBABILITY THRESHOLD TO DISCARD AA TYPE PREDICTIONS BELOW IT
                try:
                    max_prob = np.max([duplet[1] for duplet in iminus1aaTypesProbTuple_list])   # highest aa type probability
                except ValueError:  # if no aa type predictions exist for this TAAIG group
                    max_prob = 0
                if max_prob > 10e-10:   # Dedice the probability threshold
                    prob_threshold = 1000
                elif max_prob > 10e-20:
                    prob_threshold = 10000
                elif max_prob <= 10e-20:
                    prob_threshold = 100000
                print("DEBUG: max_prob", max_prob, "prob_threshold", prob_threshold, "args.LOG=", args.LOG)
                aatype_list, prob_list = [], []
                for duplet in iminus1aaTypesProbTuple_list:
                    if args.DELETE_AA_TYPE_PREDICTIONS == True:     # delete low probability aa type predictions
                        try:
                            if max_prob/float(duplet[1]) < prob_threshold:
                                aatype_list.append(duplet[0])
                                prob_list.append(duplet[1])  # unsorted list of probabilities
                        except ZeroDivisionError:
                            print("DEBUG: ZeroDivisionError ", max_prob, "/", float(duplet[1]))
                    elif args.DELETE_AA_TYPE_PREDICTIONS == False:
                        aatype_list.append(duplet[0])
                        prob_list.append(duplet[1])  # unsorted list of probabilities
                print("DEBUG: prob_list=", prob_list)
                if len(prob_list) == 1:     
                    zscore_array = np.array([10.00000])     # default Z-score in the case of just one prediction
                elif len(prob_list) > 2:
                    if args.LOG == True:
                        try:
                            ratio = float(np.max(prob_list))/np.min(prob_list)    # if the min probability is at least 3 orders of magnitude smaller, convert them to logarithmic scale
                        except ZeroDivisionError:
                            print("DEBUG: ZeroDivisionError ", np.max(prob_list), "/", np.min(prob_list)) 
                            ratio = -1.0
                        if ratio > 1000:
                            zscore_array = zscore(np.log(prob_list))  # convert probabilities the logarithms and then to Z-scores
                        else:
                            zscore_array = zscore(prob_list)  # convert probabilities to Z-scores
                    elif args.LOG == False:
                        zscore_array = zscore(prob_list)  # convert probabilities to Z-scores
                elif len(prob_list) == 2:   # if only 2 predictions exist, set the 2nd to 0.000 to allow you to use high -zcutoff values (<= 0)
                    zscore_array = np.array([1.00000, 0.00000]) 
                else:   # if no aa type prediction was made create an empty array
                    zscore_array = np.array([])
                print("DEBUG2: saving zscore_array = ", zscore_array)
                iminus1aaTypesZscoreTuple_list = []
                iminus1aaTypesCutoffProbTuple_list = []
                for aatype, Zscore, prob in zip(aatype_list, zscore_array, prob_list):
                    ## ONLY IF THE Z-SCORE OF THE AA TYPE PREDICTION IS GREATER THAN THE CUTOFF AND THE AA TYPE PREDICTIONS ARE AT LEAST args.MIN_NUM_OF_PREDICTIONS, THEN INCLUDE IT IN THE LIST
                    if Zscore > args.ZSCORE_ASSIGNMENT_CUTOFF and len(zscore_array) >= args.MIN_NUM_OF_PREDICTIONS:
                        aatypeProb_tuple = (aatype, prob)
                        aatypeZscore_tuple = (aatype, Zscore)  # replace the probability with the respective Z-score
                        iminus1aaTypesZscoreTuple_list.append(aatypeZscore_tuple)
                        iminus1aaTypesCutoffProbTuple_list.append(aatypeProb_tuple) # but also save the aatype probabilities that satisfy the Z-score cutoff
                    ## IF THE AA TYPE PREDICTIONS ARE LESS THAN args.MIN_NUM_OF_PREDICTIONS, THEN INCLUDED ALL OF THEM IN THE LIST
                    elif len(zscore_array) < args.MIN_NUM_OF_PREDICTIONS:
                        aatypeProb_tuple = (aatype, prob)
                        aatypeZscore_tuple = (aatype, Zscore)  # replace the probability with the respective Z-score
                        iminus1aaTypesZscoreTuple_list.append(aatypeZscore_tuple)
                        iminus1aaTypesCutoffProbTuple_list.append(aatypeProb_tuple) # but also save the aatype probabilities that satisfy the Z-score cutoff
                iaaindex_iminus1aaTypesZscoreTupleList_dict[iaaindex] = iminus1aaTypesZscoreTuple_list
                iaaindex_iminus1aaTypesCutoffProbTupleList_dict[iaaindex] = iminus1aaTypesCutoffProbTuple_list
            
            
            # Save ALL amino acid type predictions with the respective probabilities
            #with open("amino_acid_type_prediction_probabilities.mcutoff"+str(args.RESONANCE_MATCH_CUTOFF)+".acutoff"+str(args.ASSIGNMENT_CUTOFF)+".zacutoff"+str(args.ZSCORE_ASSIGNMENT_CUTOFF)
                      #+".wH"+str(args.H_weight)+".wC"+str(args.C_weight)+".tolH"+str(args.tolH)+".tolC"+str(args.tolC)+".rtolH"+str(args.rtolH)+".rtolN"+str(args.rtolN), 'w') as f:
            with open("amino_acid_type_prediction_probabilities", 'w') as f:
                f.write("i AAIG\tpossible i-1 aa types\n")
                print("Probabilities:")
                for i_aaindex in list(iaaindex_iminus1aaTypesProbTupleList_dict.keys()):
                    print(i_aaindex,"---> (i-1) aa type",sorted(iaaindex_iminus1aaTypesProbTupleList_dict[i_aaindex], key=itemgetter(1), reverse=True))
                    f.write(i_aaindex + " ---> (i-1) aa type " + str(sorted(iaaindex_iminus1aaTypesProbTupleList_dict[i_aaindex], key=itemgetter(1), reverse=True)) + "\n")
            
            # Save amino acid type predictions that are above the Z-Score cutoff, along with the respective probabilities
            #with open("amino_acid_type_prediction_probabilities.mcutoff"+str(args.RESONANCE_MATCH_CUTOFF)+".acutoff"+str(args.ASSIGNMENT_CUTOFF)+".zacutoff"+str(args.ZSCORE_ASSIGNMENT_CUTOFF)
                      #+".wH"+str(args.H_weight)+".wC"+str(args.C_weight)+".tolH"+str(args.tolH)+".tolC"+str(args.tolC)+".rtolH"+str(args.rtolH)+".rtolN"+str(args.rtolN), 'w') as f:
            with open("amino_acid_type_prediction_probabilities_above_cutoff.zacutoff"+str(args.ZSCORE_ASSIGNMENT_CUTOFF), 'w') as f:
                f.write("i AAIG\tpossible i-1 aa types\n")
                print("Probabilities:")
                for i_aaindex in list(iaaindex_iminus1aaTypesCutoffProbTupleList_dict.keys()):
                    print(i_aaindex,"---> (i-1) aa type",sorted(iaaindex_iminus1aaTypesCutoffProbTupleList_dict[i_aaindex], key=itemgetter(1), reverse=True))
                    f.write(i_aaindex + " ---> (i-1) aa type " + str(sorted(iaaindex_iminus1aaTypesCutoffProbTupleList_dict[i_aaindex], key=itemgetter(1), reverse=True)) + "\n")
            
            # Save amino acid type predictions with the respective Z-scores
            #with open("amino_acid_type_prediction_Z-scores.mcutoff"+str(args.RESONANCE_MATCH_CUTOFF)+".acutoff"+str(args.ASSIGNMENT_CUTOFF)+".zacutoff"+str(args.ZSCORE_ASSIGNMENT_CUTOFF)+
                      #".wH"+str(args.H_weight)+".wC"+str(args.C_weight)+".tolH"+str(args.tolH)+".tolC"+str(args.tolC)+".rtolH"+str(args.rtolH)+".rtolN"+str(args.rtolN), 'w') as f:
            with open("amino_acid_type_prediction_Z-scores.zacutoff"+str(args.ZSCORE_ASSIGNMENT_CUTOFF), 'w') as f:
                f.write("i AAIG\tpossible i-1 aa types\n")
                print("Z-scores:")
                for i_aaindex in list(iaaindex_iminus1aaTypesZscoreTupleList_dict.keys()):
                    print(i_aaindex,"---> (i-1) aa type",sorted(iaaindex_iminus1aaTypesZscoreTupleList_dict[i_aaindex], key=itemgetter(1), reverse=True))
                    f.write(i_aaindex + " ---> (i-1) aa type " + str(sorted(iaaindex_iminus1aaTypesZscoreTupleList_dict[i_aaindex], key=itemgetter(1), reverse=True)) + "\n")
        else:
            ##########################################################################################################################################################################
            ##  LOAD AMINO ACID TYPE PREDICTIONS FOR FASTER DEBUGGING
            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            print("")    # change line after the previous progress bar
            print("Loading amino acid prediction files ...")
            minimum_Zscore = 10000000
            iaaindex_iminus1aaTypesProbPoolTupleList_dict = OrderedDict() # ordereddict with keys the AAIG of residue i and values lists of tuples of the form
                                                                        # (residue i-1 matching aa type, Z-score)
            with open(args.POOL_AA_TYPES_FILE, 'r') as f:
                pool_aa_type_file_contents = f.readlines()
                for line in pool_aa_type_file_contents[1:]:
                    #print "DEBUG: line=",line
                    word_list = re.sub('---> \(i-1\) aa type', '', line).split()
                    key = word_list[0]
                    iaaindex_iminus1aaTypesProbPoolTupleList_dict[key] = []
                    values_string = ''.join(word_list[1:])
                    elements_string = re.sub('[\(\)\[\]\'\"]', '',  values_string ).split("),(")[0]
                    elements_list = elements_string.split(",")
                    for aa, Zscore in zip(elements_list[0::2], elements_list[1::2]):
                        #print "aa=",aa,"probability=",probability
                        duplet = (aa, float(Zscore))
                        iaaindex_iminus1aaTypesProbPoolTupleList_dict[key].append(duplet)
                        
            ## Now Calculate Z-scores from the pool of aa type prediction probabilites and keep only those predictions above the cutoff
            iaaindex_iminus1aaTypesZscoreTupleList_dict = OrderedDict() # ordereddict with keys the AAIG of residue i and values lists of tuples of the form (residue i-1 matching aa type, Z-score)
            iaaindex_iminus1aaTypesCutoffProbTupleList_dict = OrderedDict() # contains only the tuples (residue i-1 matching aa type, probability) with Z-score above the cutoff
            for iaaindex, iminus1aaTypesProbTuple_list in list(iaaindex_iminus1aaTypesProbPoolTupleList_dict.items()):
                # DECIDE THE PROBABILITY THRESHOLD TO DISCARD AA TYPE PREDICTIONS BELOW IT
                try:
                    max_prob = np.max([duplet[1] for duplet in iminus1aaTypesProbTuple_list])   # highest aa type probability
                except ValueError:  # if no aa type predictions exist for this TAAIG group
                    max_prob = 0
                if max_prob > 10e-10:   # Dedice the probability threshold
                    prob_threshold = 1000
                elif max_prob > 10e-20:
                    prob_threshold = 10000
                elif max_prob <= 10e-20:
                    prob_threshold = 100000
                print("DEBUG: max_prob", max_prob, "prob_threshold", prob_threshold, "args.LOG=", args.LOG)
                aatype_list, prob_list = [], []
                for duplet in iminus1aaTypesProbTuple_list:
                    if args.DELETE_AA_TYPE_PREDICTIONS == True:     # delete low probability aa type predictions
                        try:
                            if max_prob/float(duplet[1]) < prob_threshold:
                                aatype_list.append(duplet[0])
                                prob_list.append(duplet[1])  # unsorted list of probabilities
                        except ZeroDivisionError:
                            print("DEBUG: ZeroDivisionError ", max_prob, "/", float(duplet[1]))
                    elif args.DELETE_AA_TYPE_PREDICTIONS == False:
                        aatype_list.append(duplet[0])
                        prob_list.append(duplet[1])  # unsorted list of probabilities
                print("DEBUG: prob_list=", prob_list)
                if len(prob_list) == 1:     
                    zscore_array = np.array([10.00000])     # default Z-score in the case of just one prediction
                elif len(prob_list) > 1:
                    if args.LOG == True:
                        try:
                            ratio = float(np.max(prob_list))/np.min(prob_list)    # if the min probability is at least 3 orders of magnitude smaller, convert them to logarithmic scale
                        except ZeroDivisionError:
                            print("DEBUG: ZeroDivisionError ", np.max(prob_list), "/", np.min(prob_list)) 
                            ratio = -1.0
                        if ratio > 1000:
                            zscore_array = zscore(np.log(prob_list))  # convert probabilities the logarithms and then to Z-scores
                        else:
                            zscore_array = zscore(prob_list)  # convert probabilities to Z-scores
                    elif args.LOG == False:
                        zscore_array = zscore(prob_list)  # convert probabilities to Z-scores
                else:   # if no aa type prediction was made create an empty array
                    zscore_array = np.array([])
                print("DEBUG1: saving zscore_array = ", zscore_array)
                iminus1aaTypesZscoreTuple_list = []
                iminus1aaTypesCutoffProbTuple_list = []
                for aatype, Zscore, prob in zip(aatype_list, zscore_array, prob_list):
                    if Zscore < minimum_Zscore:
                        minimum_Zscore = Zscore
                    ## ONLY IF THE Z-SCORE OF THE AA TYPE PREDICTION IS GREATER THAN THE CUTOFF AND THE AA TYPE PREDICTIONS ARE AT LEAST args.MIN_NUM_OF_PREDICTIONS, THEN INCLUDE IT IN THE LIST
                    if Zscore > args.ZSCORE_ASSIGNMENT_CUTOFF and len(zscore_array) >= args.MIN_NUM_OF_PREDICTIONS:
                        aatypeProb_tuple = (aatype, prob)
                        aatypeZscore_tuple = (aatype, Zscore)  # replace the probability with the respective Z-score
                        iminus1aaTypesZscoreTuple_list.append(aatypeZscore_tuple)
                        iminus1aaTypesCutoffProbTuple_list.append(aatypeProb_tuple) # but also save the aatype probabilities that satisfy the Z-score cutoff
                    ## IF THE AA TYPE PREDICTIONS ARE LESS THAN args.MIN_NUM_OF_PREDICTIONS, THEN INCLUDED ALL OF THEM IN THE LIST
                    elif len(zscore_array) < args.MIN_NUM_OF_PREDICTIONS:
                        aatypeProb_tuple = (aatype, prob)
                        aatypeZscore_tuple = (aatype, Zscore)  # replace the probability with the respective Z-score
                        iminus1aaTypesZscoreTuple_list.append(aatypeZscore_tuple)
                        iminus1aaTypesCutoffProbTuple_list.append(aatypeProb_tuple) # but also save the aatype probabilities that satisfy the Z-score cutoff
                iaaindex_iminus1aaTypesZscoreTupleList_dict[iaaindex] = iminus1aaTypesZscoreTuple_list
                iaaindex_iminus1aaTypesCutoffProbTupleList_dict[iaaindex] = iminus1aaTypesCutoffProbTuple_list
            
            
            iaaindex_iminus1aaTypesProbTupleList_dict = OrderedDict()   # ordereddict with keys the AAIG of residue i and values lists of tuples of the form
                                                                        # (residue i-1 matching aa type, average probability)
            with open(args.COMPLETE_AA_TYPES_FILE, 'r') as f:
                complete_aa_type_file_contents = f.readlines()
                for line in complete_aa_type_file_contents[1:]:
                    #print "DEBUG: line=",line
                    word_list = re.sub('---> \(i-1\) aa type', '', line).split()
                    key = word_list[0]
                    iaaindex_iminus1aaTypesProbTupleList_dict[key] = []
                    values_string = ''.join(word_list[1:])
                    elements_string = re.sub('[\(\)\[\]\'\"]', '',  values_string ).split("),(")[0]
                    elements_list = elements_string.split(",")
                    for aa, hist_prob in zip(elements_list[0::2], elements_list[1::2]):
                        #print "aa=",aa,"probability=",probability
                        duplet = (aa, float(hist_prob))
                        iaaindex_iminus1aaTypesProbTupleList_dict[key].append(duplet)
            
            #for i_aaindex in iaaindex_iminus1aaTypesZscoreTupleList_dict.keys():
            #    print i_aaindex,"---> (i-1) aa type",sorted(iaaindex_iminus1aaTypesZscoreTupleList_dict[i_aaindex], key=itemgetter(1), reverse=True)
            ##########################################################################################################################################################################
            
            with open("amino_acid_type_prediction_Z-scores.zacutoff"+str(args.ZSCORE_ASSIGNMENT_CUTOFF), 'w') as f:
                f.write("i AAIG\tpossible i-1 aa types\n")
                print("Z-scores:")
                for i_aaindex in list(iaaindex_iminus1aaTypesZscoreTupleList_dict.keys()):
                    print(i_aaindex,"---> (i-1) aa type",sorted(iaaindex_iminus1aaTypesZscoreTupleList_dict[i_aaindex], key=itemgetter(1), reverse=True))
                    f.write(i_aaindex + " ---> (i-1) aa type " + str(sorted(iaaindex_iminus1aaTypesZscoreTupleList_dict[i_aaindex], key=itemgetter(1), reverse=True)) + "\n")
        
        #sys.exit(1)
        ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        ##              CALCULATE BAYESIAN STATISTICS                        #
        ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        
        ## CALCULATE THE PROBABILITY P[AATYPE(I-1)]
        ############################### NEW WAY #####################################
        template_sequence_string = ""
        with open(args.template_sequence_file, 'r') as f:
            for line in f:
                if not re.match("^>", line):
                    template_sequence_string += re.sub(r"[^A-Z]", "", line)
        template_sequence_list = list(template_sequence_string)
        
        aatype_P_dict = {}    # dict with the P[aatype(i-1)] for every aatype(i-1): the percentage of each aatype in the template sequence
        for aa_type in list(aatype_maxH_C_pairs_dict.keys()):
            aatype_P_dict[aa_type] = template_sequence_list.count(aa3to1_dict[aa_type]) / float(len(template_sequence_list))
        
        ## DO A MATHEMATIC TRANSFORM ON P
        # iaaindex_iminus1aaTypesProbTupleList_dict is an ordereddict with keys the AAIG of residue i and values lists of tuples of the form
        # (residue i-1 matching aa type, average probability)
        all_prob_list = []
        for iaaindex, iminus1aaTypesProbTuple_list in list(iaaindex_iminus1aaTypesProbTupleList_dict.items()):
            for duplet in iminus1aaTypesProbTuple_list:
                prob = duplet[1]
                all_prob_list.append(prob)
        
        all_prob_array = np.array(all_prob_list)
        transformed_all_prob_array = transform(all_prob_array)
        
        iaaindex_iminus1aaTypesTransformedProbTupleList_dict = OrderedDict()
        array_index = 0
        for iaaindex, iminus1aaTypesProbTuple_list in list(iaaindex_iminus1aaTypesProbTupleList_dict.items()):
            iaaindex_iminus1aaTypesTransformedProbTupleList_dict[iaaindex] = []
            for duplet in iminus1aaTypesProbTuple_list:
                aatype = duplet[0]
                trans_prob = transformed_all_prob_array[array_index]
                iaaindex_iminus1aaTypesTransformedProbTupleList_dict[iaaindex].append((aatype, trans_prob))
                array_index += 1
        
        
        #for k in iaaindex_iminus1aaTypesProbTupleList_dict.keys():
        #    print(k)
        #    print(iaaindex_iminus1aaTypesProbTupleList_dict[k])
        #    print(iaaindex_iminus1aaTypesTransformedProbTupleList_dict[k])
        #    
        #sys.exit(1)
        
        ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        ##       BUILD TREES OF POSSIBLE PEPTIDE SEQUENCES                   #
        ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        
        if os.path.exists('tmp_peptide_folder/') == True and args.RESUME == False:   # if the folder exists and this is not a resumed run, remove it and make a new one
            shutil.rmtree('tmp_peptide_folder/', ignore_errors=True)
            os.makedirs('tmp_peptide_folder/')
        elif os.path.exists('tmp_peptide_folder/') == False and args.RESUME == False: # if the folder does not exists, make a new one
            os.makedirs('tmp_peptide_folder/')
        elif os.path.exists('tmp_peptide_folder/') == False and args.RESUME == True:
            print("ERROR: folder tmp_peptide_folder/ with peptide sequence files does not exist! The process cannot be resumed!")
        # LOAD ALL PREVIOUSLY WRITEN PEPTIDE FILES
        fnames=os.listdir('tmp_peptide_folder/')
        fpattern = re.compile('chainIndex_[0-9]+.pickle.bz2')
        peptide_files_list = list(filter(fpattern.search, fnames))
        chainIndices2skip_set = set()
        for peptide_file in peptide_files_list:
            mo = re.search('chainIndex_([0-9]+).pickle.bz2', peptide_file)
            if mo:
                chainIndex = int(mo.group(1))
                chainIndices2skip_set.add(chainIndex)
        
        total_chain_number = len(all_chainScore_list)
        remaining_chainIndex_list = []
        remaining_chainScore_list = []
        # SET SHARED VARIABLES BETWEEN THREADS
        shared.setConst(TOTAL_CHAIN_NUMBER=total_chain_number)
        shared.setConst(ARGS=args)
        #print "DEBUG: = iaaindex_iminus1aaTypesCutoffProbTupleList_dict=", iaaindex_iminus1aaTypesCutoffProbTupleList_dict
        shared.setConst(IAAINDEX_IMINUS1AATYPESCUTOFFPROBTUPLELIST_DICT=iaaindex_iminus1aaTypesCutoffProbTupleList_dict)
        shared.setConst(AA3TO1_DICT=aa3to1_dict)
        shared.setConst(IAAINDEX_IMINUS1AATYPESTRANSFORMEDPROBTUPLELIST_DICT=iaaindex_iminus1aaTypesTransformedProbTupleList_dict)
        #shared.setConst(IAAINDEX_IMINUS1AATYPESPAATYPETUPLELIST_DICT=iaaindex_iminus1aaTypesPaatypeTupleList_dict)
        shared.setConst(AATYPE_P_DICT=aatype_P_dict)
        shared.setConst(I_IMINUS1_DICT=i_iminus1_dict)
        # OMMIT BUILDING PEPTIDE FOR THE CHAINS THAT ALREADY EXIST IN tmp_peptide_folder/
        for chainIndex, chainScore in enumerate(all_chainScore_list):
            if not chainIndex in chainIndices2skip_set:
                remaining_chainIndex_list.append(chainIndex)
                remaining_chainScore_list.append(chainScore)
        
        results = list(futures.map(build_Peptide_Tree, remaining_chainIndex_list, remaining_chainScore_list))   # build peptide tree from chain
        # results will be empty!!
        #for chainIndex, chainScore in zip(remaining_chainIndex_list, remaining_chainScore_list):
        #    build_Peptide_Tree(remaining_chainIndex_list, remaining_chainScore_list)
        
        ##
        ## Write peptides to files
        ##
        # LOAD ALL WRITEN PEPTIDE FILES
        progbar = ProgressBar(100)
        peptide_file_num = 0
        #while peptide_file_num < total_chain_number:    # CREATE MISSING PEPTIDE TREES
        fnames=os.listdir('tmp_peptide_folder/')
        fpattern = re.compile('chainIndex_[0-9]+.pickle.bz2')
        peptide_files_list = list(filter(fpattern.search, fnames))
        peptide_file_num = len(peptide_files_list)
        for chainIndex, chainScore in enumerate(all_chainScore_list):
            if not 'chainIndex_'+str(chainIndex)+'.pickle.bz2' in peptide_files_list:
                build_Peptide_Tree(chainIndex, chainScore)
        
        print("Saving all the peptide sequences to peptides.list file ...")
        peptide_filesList_list = chunkIt(peptide_files_list, scoop.SIZE)
        parts_list = [str(part) for part in range(1, scoop.SIZE +1)]
        try:
            shared.setConst(ARGS=args)
        except TypeError:
            pass
        try:
            shared.setConst(PROGBAR=progbar)
        except TypeError:
            pass
        shared.setConst(NEW_ALL_CHAINSCORE_LIST=all_chainScore_list)
        results = list(futures.map(save_peptides_to_file, peptide_filesList_list, parts_list))
        # results will be an empty list
        # Now concatenate all parts of peptides.list and peptides.fasta
        fasta_filehandler = open("peptides."+str(args.MAX_PEPTIDE_LENGTH)+"mers.fasta", 'w')
        f = open("peptides."+str(args.MAX_PEPTIDE_LENGTH)+"mers.list", 'w')
        index = 0
        for part in parts_list:
            with open("peptides."+str(args.MAX_PEPTIDE_LENGTH)+"mers.list.chunk"+part, 'r') as fin:
                for line in fin:
                    mo = re.search("^(.* peptide )P_[0-9]+_[0-9]+ = ([A-Z]+)( .*)$", line)
                    if mo:
                        line_part1 = mo.group(1)
                        peptide_seq = mo.group(2)
                        line_part2 = mo.group(3)
                        peptide_name = "P_" + str(len(peptide_seq)) + "_" + str(index)
                        #f.write("Chain probability = "+str(chain[-1])+" peptide P_"+str(len(peptide))+"_"+str(index+1) + " = "+''.join(peptide)+" "+'-'.join(chain[0:-1])+" "+','.join(map(str, peptideProbList[:-1]))+"\n")
                        #fasta_filehandler.write(">P_"+str(len(peptide))+"_"+str(index+1)+"\n"+''.join(peptide)+"\n")
                        f.write(line_part1 + peptide_name + " = " + peptide_seq + line_part2 + "\n")
                        fasta_filehandler.write(">" + peptide_name + "\n" + peptide_seq + "\n")
                        index += 1
            os.remove("peptides."+str(args.MAX_PEPTIDE_LENGTH)+"mers.list.chunk"+part)
        shutil.rmtree('tmp_peptide_folder/', ignore_errors=True)
    
    except:
        type, value, tb = sys.exc_info()
        lines = traceback.format_exception(type, value, tb)
        print((''.join(lines)))
        raise