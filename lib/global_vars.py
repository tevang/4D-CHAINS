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


from lib.global_func import tree

## Set global variables
code3 = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
code1 = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
aa1to3_dict = dict((c1,c3) for c1,c3 in zip(code1,code3))
aa3to1_dict = dict((c3,c1) for c3,c1 in zip(code3,code1))

# Excluding the aromatic C-H
allowed_aa_atoms_dict = {
"ALA" : ["HA", "HB", "CA", "CB", "N", "H"],
"ARG" : ["HA", "HB2", "HB3", "HG2", "HG3", "HD2", "HD3", "CA", "CB", "CG", "CD", "N", "H"],
"ASP" : ["HA", "HB2", "HB3", "CA", "CB", "N", "H"],
"ASN" : ["HA", "HB2", "HB3", "CA", "CB", "N", "H"],
"CYS" : ["HA", "HB2", "HB3", "CA", "CB", "N", "H"],
"GLU" : ["HA", "HB2", "HB3", "HG2", "HG3", "CA", "CB", "CG", "N", "H"],
"GLN" : ["HA", "HB2", "HB3", "HG2", "HG3", "CA", "CB", "CG", "N", "H"],
"GLY" : ["HA2", "HA3", "CA", "N", "H"],
"HIS" : ["HA", "HB2", "HB3", "CA", "CB", "N", "H"],
"ILE" : ["HA", "HB", "HG12", "HG13", "HG2", "HD1", "CA", "CB", "CG1", "CG2", "CD1", "N", "H"],
"LEU" : ["HA", "HB2", "HB3", "HG", "HD1", "HD2", "CA", "CB", "CG", "CD1", "CD2", "N", "H"],
"LYS" : ["HA", "HB2", "HB3", "HG2", "HG3", "HD2", "HD3", "HE2", "HE3", "CA", "CB", "CG", "CD", "CE", "N", "H"],
"MET" : ["HA", "HB2", "HB3", "HG2", "HG3", "CA", "CB", "CG", "N", "H"],
"PHE" : ["HA", "HB2", "HB3", "CA", "CB", "N", "H"],
"PRO" : ["HA", "HB2", "HB3", "HG2", "HG3", "HD2", "HD3", "CA", "CB", "CG", "CD", "N"], # Prolines are not detected by the method at position "i" due to lack of HN hydrogen
"SER" : ["HA", "HB2", "HB3", "CA", "CB", "N", "H"],
"THR" : ["HA", "HB", "HG2", "CA", "CB", "CG2", "N", "H"],
"TRP" : ["HA", "HB2", "HB3", "CA", "CB", "N", "H"],
"TYR" : ["HA", "HB2", "HB3", "CA", "CB", "N", "H"],
"VAL" : ["HA", "HB", "HG1", "HG2", "CA", "CB", "CG1", "CG2", "N", "H"]
}

aatype_maxH_C_pairs_dict = {
"ALA" : 2,
"ARG" : 7,
"ASP" : 3,
"ASN" : 3,
"CYS" : 3,
"GLU" : 5,
"GLN" : 5,
"GLY" : 2,
"HIS" : 3,
"ILE" : 6,
"LEU" : 6,
"LYS" : 9,
"MET" : 5,
"PHE" : 3,
"PRO" : 7, # Prolines are not detected by the method at position "i" due to lack of HN hydrogen
"SER" : 3,
"THR" : 3,
"TRP" : 3,
"TYR" : 3,
"VAL" : 4
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
    "CB" : ["HB2", "HB3"],
    "CG" : ["HG2", "HG3"],
    "CD" : ["HD2", "HD3"]
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
    "CB" : ["HB2", "HB3"],
    "CG" : ["HG2", "HG3"],
}
#"GLN" : ["HA", "HB2", "HB3", "HG2", "HG3", "CA", "CB", "CG"],
aa_carbonBondedHydrogensDict_dict["GLN"] = {
    "CA" : ["HA"],
    "CB" : ["HB2", "HB3"],
    "CG" : ["HG2", "HG3"]
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
    "CB" : ["HB"],
    "CG1" : ["HG12", "HG13"],
    "CG2" : ["HG2"],
    "CD1" : ["HD1"]
}
#"LEU" : ["HA", "HB2", "HB3", "HG", "HD1", "HD2", "CA", "CB", "CG", "CD1", "CD2"],
aa_carbonBondedHydrogensDict_dict["LEU"] = {
    "CA" : ["HA"],
    "CB" : ["HB2", "HB3"],
    "CG" : ["HG"],
    "CD1" : ["HD1"],
    "CD2" : ["HD2"]
}
#"LYS" : ["HA", "HB2", "HB3", "HG2", "HG3", "HD2", "HD3", "HE2", "HE3", "CA", "CB", "CG", "CD", "CE"],
aa_carbonBondedHydrogensDict_dict["LYS"] = {
    "CA" : ["HA"],
    "CB" : ["HB2", "HB3"],
    "CG" : ["HG2", "HG3"],
    "CD" : ["HD2", "HD3"],
    "CE" : ["HE2", "HE3"]
}
#"MET" : ["HA", "HB2", "HB3", "HG2", "HG3", "CA", "CB", "CG"],
aa_carbonBondedHydrogensDict_dict["MET"] = {
    "CA" : ["HA"],
    "CB" : ["HB2", "HB3"],
    "CG" : ["HG2", "HG3"],
    "CE" : ["HE1", "HE2", "HE3"]
}
#"PHE" : ["HA", "HB2", "HB3", "CA", "CB"],
aa_carbonBondedHydrogensDict_dict["PHE"] = {
    "CA" : ["HA"],
    "CB" : ["HB2", "HB3"],
    "CD1" : ["HD1"],
    "CE1" : ["HE1"],
    "CZ" : ["HZ"]
}
#"PRO" : ["HA", "HB2", "HB3", "HG2", "HG3", "HD2", "HD3", "CA", "CB", "CG", "CD"],
aa_carbonBondedHydrogensDict_dict["PRO"] = {
    "CA" : ["HA"],
    "CB" : ["HB2", "HB3"],
    "CG" : ["HG2", "HG3"],
    "CD" : ["HD2", "HD3"]
}
#"SER" : ["HA", "HB2", "HB3", "CA", "CB"],
aa_carbonBondedHydrogensDict_dict["SER"] = {
    "CA" : ["HA"],
    "CB" : ["HB2", "HB3"]
}
#"THR" : ["HA", "HB", "HG2", "CA", "CB", "CG2"],
aa_carbonBondedHydrogensDict_dict["THR"] = {
    "CA" : ["HA"],
    "CB" : ["HB"],
    "CG2" : ["HG2"]
}
#"TRP" : ["HA", "HB2", "HB3", "CA", "CB"],
aa_carbonBondedHydrogensDict_dict["TRP"] = {
    "CA" : ["HA"],
    "CB" : ["HB2", "HB3"]
}
#"TYR" : ["HA", "HB2", "HB3", "CA", "CB"],
aa_carbonBondedHydrogensDict_dict["TYR"] = {
    "CA" : ["HA"],
    "CB" : ["HB2", "HB3"],
    "CD1" : ["HD1"],
    "CE1" : ["HE1"],
}
#"VAL" : ["HA", "HB", "HG1", "HG2", "CA", "CB", "CG1", "CG2"]
aa_carbonBondedHydrogensDict_dict["VAL"] = {
    "CA" : ["HA"],
    "CB" : ["HB"],
    "CG1" : ["HG1"],
    "CG2" : ["HG2"]
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



# THIS DICT BELOW IS A TEMPORARY FIX. ILE CG1 is the problem, which according to BMRB has 2 protons, HG12 and HG13, are distinguishable, but in sometimes in real experiments they
# are degenerate.
# Degenerate protons (protons that can be never distinguished according to BMRB; usually methyl protons)
aa_CarbonListwithDegenerateH_dict = {}
aa_CarbonListwithDegenerateH_dict["ILE"] = ["CD1", "CG2"]
# aa_CarbonListwithDegenerateH_dict["LEU"] = ["CD1", "CD2"]
# aa_CarbonListwithDegenerateH_dict["VAL"] = ["CG1", "CG2"]
aa_CarbonListwithDegenerateH_dict["THR"] = ["CG2"]
aa_CarbonListwithDegenerateH_dict["MET"] = ["CE"]


aatype_carbon_nondegenerateHlist_mdict = tree() # names of non-degenerate methylene protons (in that case they are named HG2, HG3, etc.)
aatype_carbon_nondegenerateHlist_mdict['GLY']['CA'] = ['HA2', 'HA3']

aatype_carbon_nondegenerateHlist_mdict['ASP']['CB'] = ['HB2', 'HB3']
aatype_carbon_nondegenerateHlist_mdict['GLU']['CB'] = ['HB2', 'HB3']
aatype_carbon_nondegenerateHlist_mdict['HIS']['CB'] = ['HB2', 'HB3']
aatype_carbon_nondegenerateHlist_mdict['LYS']['CB'] = ['HB2', 'HB3']
aatype_carbon_nondegenerateHlist_mdict['LEU']['CB'] = ['HB2', 'HB3']
aatype_carbon_nondegenerateHlist_mdict['ASN']['CB'] = ['HB2', 'HB3']
aatype_carbon_nondegenerateHlist_mdict['PHE']['CB'] = ['HB2', 'HB3']
aatype_carbon_nondegenerateHlist_mdict['PRO']['CB'] = ['HB2', 'HB3']
aatype_carbon_nondegenerateHlist_mdict['GLN']['CB'] = ['HB2', 'HB3']
aatype_carbon_nondegenerateHlist_mdict['ARG']['CB'] = ['HB2', 'HB3']
aatype_carbon_nondegenerateHlist_mdict['SER']['CB'] = ['HB2', 'HB3']
aatype_carbon_nondegenerateHlist_mdict['CYS']['CB'] = ['HB2', 'HB3']
aatype_carbon_nondegenerateHlist_mdict['MET']['CB'] = ['HB2', 'HB3']
aatype_carbon_nondegenerateHlist_mdict['TRP']['CB'] = ['HB2', 'HB3']
aatype_carbon_nondegenerateHlist_mdict['TYR']['CB'] = ['HB2', 'HB3']

aatype_carbon_nondegenerateHlist_mdict['ARG']['CG'] = ['HG2', 'HG3']
aatype_carbon_nondegenerateHlist_mdict['GLN']['CG'] = ['HG2', 'HG3']
aatype_carbon_nondegenerateHlist_mdict['GLU']['CG'] = ['HG2', 'HG3']
aatype_carbon_nondegenerateHlist_mdict['ILE']['CG1'] = ['HG12', 'HG13']
aatype_carbon_nondegenerateHlist_mdict['LYS']['CG'] = ['HG2', 'HG3']
aatype_carbon_nondegenerateHlist_mdict['MET']['CG'] = ['HG2', 'HG3']
aatype_carbon_nondegenerateHlist_mdict['PRO']['CG'] = ['HG2', 'HG3']

aatype_carbon_nondegenerateHlist_mdict['ARG']['CD'] = ['HD2', 'HD3']
aatype_carbon_nondegenerateHlist_mdict['LYS']['CD'] = ['HD2', 'HD3']
aatype_carbon_nondegenerateHlist_mdict['PRO']['CD'] = ['HD2', 'HD3']

aatype_carbon_nondegenerateHlist_mdict['LYS']['CE'] = ['HE2', 'HE3']


aatype_carbon_degenerateH_mdict = tree() # names of degenerate methylene protons (in that case are named QB, QG, etc.)
aatype_carbon_degenerateH_mdict['GLY']['CA'] = 'QA'

aatype_carbon_degenerateH_mdict['ASP']['CB'] = 'QB'
aatype_carbon_degenerateH_mdict['GLU']['CB'] = 'QB'
aatype_carbon_degenerateH_mdict['HIS']['CB'] = 'QB'
aatype_carbon_degenerateH_mdict['LYS']['CB'] = 'QB'
aatype_carbon_degenerateH_mdict['LEU']['CB'] = 'QB'
aatype_carbon_degenerateH_mdict['ASN']['CB'] = 'QB'
aatype_carbon_degenerateH_mdict['PHE']['CB'] = 'QB'
aatype_carbon_degenerateH_mdict['PRO']['CB'] = 'QB'
aatype_carbon_degenerateH_mdict['GLN']['CB'] = 'QB'
aatype_carbon_degenerateH_mdict['ARG']['CB'] = 'QB'
aatype_carbon_degenerateH_mdict['SER']['CB'] = 'QB'
aatype_carbon_degenerateH_mdict['CYS']['CB'] = 'QB'
aatype_carbon_degenerateH_mdict['MET']['CB'] = 'QB'
aatype_carbon_degenerateH_mdict['TRP']['CB'] = 'QB'
aatype_carbon_degenerateH_mdict['TYR']['CB'] = 'QB'

aatype_carbon_degenerateH_mdict['ARG']['CG'] = 'QG'
aatype_carbon_degenerateH_mdict['GLN']['CG'] = 'QG'
aatype_carbon_degenerateH_mdict['GLU']['CG'] = 'QG'
aatype_carbon_degenerateH_mdict['ILE']['CG1'] = 'QG1'
aatype_carbon_degenerateH_mdict['LYS']['CG'] = 'QG'
aatype_carbon_degenerateH_mdict['MET']['CG'] = 'QG'
aatype_carbon_degenerateH_mdict['PRO']['CG'] = 'QG'

aatype_carbon_degenerateH_mdict['ARG']['CD'] = 'QD'
aatype_carbon_degenerateH_mdict['LYS']['CD'] = 'QD'
aatype_carbon_degenerateH_mdict['PRO']['CD'] = 'QD'

aatype_carbon_degenerateH_mdict['LYS']['CE'] = 'QE'


# Protons that cannot be grouped together
# (I added H* and M* alternative names to the first block to recognize alternatice user proton names and not raise an error)
aatype_carbon_nongeminalHname_mdict = tree()
aatype_carbon_nongeminalHname_mdict['ALA']['CB'] = ['QB', 'HB', 'MB']
aatype_carbon_nongeminalHname_mdict['ILE']['CG2'] = ['QG2', 'HG2', 'MG2']
aatype_carbon_nongeminalHname_mdict['ILE']['CD1'] = ['QD1', 'HD1', 'MD1']
aatype_carbon_nongeminalHname_mdict['LEU']['CD1'] = ['QD1', 'HD1', 'MD1']
aatype_carbon_nongeminalHname_mdict['LEU']['CD2'] = ['QD2', 'HD2', 'MD2']
aatype_carbon_nongeminalHname_mdict['THR']['CG2'] = ['QG2', 'HG2', 'MG2']
aatype_carbon_nongeminalHname_mdict['VAL']['CG1'] = ['QG1', 'HG1', 'MG1']
aatype_carbon_nongeminalHname_mdict['VAL']['CG2'] = ['QG2', 'HG2', 'MG2']
aatype_carbon_nongeminalHname_mdict['MET']['CE'] = ['QE', 'HE', 'ME']

aatype_carbon_nongeminalHname_mdict["ALA"]["CA"] = "HA"
aatype_carbon_nongeminalHname_mdict["ARG"]["CA"] = "HA"
aatype_carbon_nongeminalHname_mdict["ASP"]["CA"] = "HA"
aatype_carbon_nongeminalHname_mdict["ASN"]["CA"] = "HA"
aatype_carbon_nongeminalHname_mdict["CYS"]["CA"] = "HA"
aatype_carbon_nongeminalHname_mdict["GLU"]["CA"] = "HA"
aatype_carbon_nongeminalHname_mdict["GLN"]["CA"] = "HA"
aatype_carbon_nongeminalHname_mdict["HIS"]["CA"] = "HA"
aatype_carbon_nongeminalHname_mdict["ILE"]["CA"] = "HA"
aatype_carbon_nongeminalHname_mdict["ILE"]["CB"] = "HB"
aatype_carbon_nongeminalHname_mdict["LEU"]["CA"] = "HA"
aatype_carbon_nongeminalHname_mdict["LEU"]["CG"] = "HG"
aatype_carbon_nongeminalHname_mdict["LYS"]["CA"] = "HA"
aatype_carbon_nongeminalHname_mdict["MET"]["CA"] = "HA"
aatype_carbon_nongeminalHname_mdict["PHE"]["CA"] = "HA"
aatype_carbon_nongeminalHname_mdict["PRO"]["CA"] = "HA"
aatype_carbon_nongeminalHname_mdict["SER"]["CA"] = "HA"
aatype_carbon_nongeminalHname_mdict["THR"]["CA"] = "HA"
aatype_carbon_nongeminalHname_mdict["THR"]["CB"] = "HB"
aatype_carbon_nongeminalHname_mdict["TRP"]["CA"] = "HA"
aatype_carbon_nongeminalHname_mdict["TYR"]["CA"] = "HA"
aatype_carbon_nongeminalHname_mdict["VAL"]["CA"] = "HA"
aatype_carbon_nongeminalHname_mdict["VAL"]["CB"] = "HB"


aatype_carbon_possibleDegenerateHlist_mdict = tree() # protons that may be degenerate (in that case are named QB, QG, etc.) or not (in that case they are named HG1, HG2, etc.)
aatype_carbon_possibleDegenerateHlist_mdict['GLY']['CA'] = ['HA2', 'HA3']

aatype_carbon_possibleDegenerateHlist_mdict['ASP']['CB'] = ['HB2', 'HB3']
aatype_carbon_possibleDegenerateHlist_mdict['GLU']['CB'] = ['HB2', 'HB3']
aatype_carbon_possibleDegenerateHlist_mdict['HIS']['CB'] = ['HB2', 'HB3']
aatype_carbon_possibleDegenerateHlist_mdict['LYS']['CB'] = ['HB2', 'HB3']
aatype_carbon_possibleDegenerateHlist_mdict['LEU']['CB'] = ['HB2', 'HB3']
aatype_carbon_possibleDegenerateHlist_mdict['ASN']['CB'] = ['HB2', 'HB3']
aatype_carbon_possibleDegenerateHlist_mdict['PHE']['CB'] = ['HB2', 'HB3']
aatype_carbon_possibleDegenerateHlist_mdict['PRO']['CB'] = ['HB2', 'HB3']
aatype_carbon_possibleDegenerateHlist_mdict['GLN']['CB'] = ['HB2', 'HB3']
aatype_carbon_possibleDegenerateHlist_mdict['ARG']['CB'] = ['HB2', 'HB3']
aatype_carbon_possibleDegenerateHlist_mdict['SER']['CB'] = ['HB2', 'HB3']
aatype_carbon_possibleDegenerateHlist_mdict['CYS']['CB'] = ['HB2', 'HB3']
aatype_carbon_possibleDegenerateHlist_mdict['MET']['CB'] = ['HB2', 'HB3']
aatype_carbon_possibleDegenerateHlist_mdict['TRP']['CB'] = ['HB2', 'HB3']
aatype_carbon_possibleDegenerateHlist_mdict['TYR']['CB'] = ['HB2', 'HB3']

aatype_carbon_possibleDegenerateHlist_mdict['ARG']['CG'] = ['HG2', 'HG3']
aatype_carbon_possibleDegenerateHlist_mdict['GLN']['CG'] = ['HG2', 'HG3']
aatype_carbon_possibleDegenerateHlist_mdict['GLU']['CG'] = ['HG2', 'HG3']
aatype_carbon_possibleDegenerateHlist_mdict['ILE']['CG'] = ['HG12', 'HG13']
aatype_carbon_possibleDegenerateHlist_mdict['LYS']['CG'] = ['HG2', 'HG3']
aatype_carbon_possibleDegenerateHlist_mdict['MET']['CG'] = ['HG2', 'HG3']
aatype_carbon_possibleDegenerateHlist_mdict['PRO']['CG'] = ['HG2', 'HG3']

aatype_carbon_possibleDegenerateHlist_mdict['ARG']['CD'] = ['HD2', 'HD3']
aatype_carbon_possibleDegenerateHlist_mdict['LYS']['CD'] = ['HD2', 'HD3']
aatype_carbon_possibleDegenerateHlist_mdict['PRO']['CD'] = ['HD2', 'HD3']

aatype_carbon_possibleDegenerateHlist_mdict['LYS']['CE'] = ['HE2', 'HE3']


## PROBABILITIES OF THE PRESENCE OF A CHEMICAL SHIFT OF A PARTICULAR NUCLEAR TYPE
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
Prob_CS["MET"]["CE"] = 0.46229
Prob_CS["MET"]["N"] = 0.682714
Prob_CS["MET"]["H"] = 0.718516
Prob_CS["MET"]["HE1"] = 0.349773    # correct it!
Prob_CS["MET"]["HE2"] = 0.349773
Prob_CS["MET"]["HE3"] = 0.34986
Prob_CS["MET"]["HE"] = 0.34986
Prob_CS["PHE"]["HA"] = 0.628538
Prob_CS["PHE"]["HB2"] = 0.612522
Prob_CS["PHE"]["HB3"] = 0.574238
Prob_CS["PHE"]["CA"] = 0.718078
Prob_CS["PHE"]["CB"] = 0.677572
Prob_CS["PHE"]["HD1"] = 0.605456
# Prob_CS["PHE"]["HD2"] = 0.509747
Prob_CS["PHE"]["HE1"] = 0.532136
# Prob_CS["PHE"]["HE2"] = 0.45499
Prob_CS["PHE"]["HZ"] = 0.381171
Prob_CS["PHE"]["CD1"] = 0.349268
# Prob_CS["PHE"]["CD2"] = 0.254757
Prob_CS["PHE"]["CE1"] = 0.305589
# Prob_CS["PHE"]["CE2"] = 0.223253
Prob_CS["PHE"]["CZ"] = 0.235562
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
Prob_CS["TYR"]["HD1"] = 0.62963
# Prob_CS["TYR"]["HD2"] = 0.532103
Prob_CS["TYR"]["CD1"] = 0.369983
#Prob_CS["TYR"]["CD2"] = 0.254576
Prob_CS["TYR"]["HE1"] = 0.59921
# Prob_CS["TYR"]["HE2"] = 0.511475
Prob_CS["TYR"]["CE1"] = 0.367236
# Prob_CS["TYR"]["CE2"] = 0.252603
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


## AVERAGE HISTOGRAM PROBABILITIES OF A PARTICULAR NUCLEAR TYPE
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


aatype_carbon_methylHydrogens_mdict = tree()   # dictionary with methyl hydrogens (geminal protons; cannot be distinguished)

aatype_carbon_methylHydrogens_mdict['ALA']['CB'] = 'QB'
aatype_carbon_methylHydrogens_mdict['ILE']['CG2'] = 'QG2'
#aatype_carbon_methylHydrogens_mdict['ILE']['CG1'] = 'QG1'
aatype_carbon_methylHydrogens_mdict['ILE']['CD1'] = 'QD1'
aatype_carbon_methylHydrogens_mdict['LEU']['CD1'] = 'QD1'
aatype_carbon_methylHydrogens_mdict['LEU']['CD2'] = 'QD2'
aatype_carbon_methylHydrogens_mdict['THR']['CG2'] = 'QG2'
aatype_carbon_methylHydrogens_mdict['VAL']['CG1'] = 'QG1'
aatype_carbon_methylHydrogens_mdict['VAL']['CG2'] = 'QG2'

# a multidict with keys the amino acid --> carbon-hydrogen pair and values two arrays, the carbon chemical shifts and the respective hydrogen chemical shifts
# this multidict holds all the information in the input BMRB database and will be used to calculate the histograms
aa_CHpair_2Dhist_mdict = tree()
aa_CHpair_2Dhist_mdict["ALA"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_mdict["ALA"]["CB-HB"] = [[], []]

#"ARG" : ["HA", "HB2", "HB3", "HG2", "HG3", "HD2", "HD3", "CA", "CB", "CG", "CD"],
aa_CHpair_2Dhist_mdict["ARG"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_mdict["ARG"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_mdict["ARG"]["CB-HB3"] = [[], []]
aa_CHpair_2Dhist_mdict["ARG"]["CG-HG2"] = [[], []]
aa_CHpair_2Dhist_mdict["ARG"]["CG-HG3"] = [[], []]
aa_CHpair_2Dhist_mdict["ARG"]["CD-HD2"] = [[], []]
aa_CHpair_2Dhist_mdict["ARG"]["CD-HD3"] = [[], []]

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
aa_CHpair_2Dhist_mdict["GLU"]["CG-HG2"] = [[], []]
aa_CHpair_2Dhist_mdict["GLU"]["CG-HG3"] = [[], []]

#"GLN" : ["HA", "HB2", "HB3", "HG2", "HG3", "CA", "CB", "CG"],
aa_CHpair_2Dhist_mdict["GLN"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_mdict["GLN"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_mdict["GLN"]["CB-HB3"] = [[], []]
aa_CHpair_2Dhist_mdict["GLN"]["CG-HG2"] = [[], []]
aa_CHpair_2Dhist_mdict["GLN"]["CG-HG3"] = [[], []]

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
aa_CHpair_2Dhist_mdict["ILE"]["CG1-HG12"] = [[], []]
aa_CHpair_2Dhist_mdict["ILE"]["CG1-HG13"] = [[], []]
aa_CHpair_2Dhist_mdict["ILE"]["CG2-HG2"] = [[], []]
aa_CHpair_2Dhist_mdict["ILE"]["CD1-HD1"] = [[], []]

#"LEU" : ["HA", "HB2", "HB3", "HG", "HD1", "HD2", "CA", "CB", "CG", "CD1", "CD2"],
aa_CHpair_2Dhist_mdict["LEU"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_mdict["LEU"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_mdict["LEU"]["CB-HB3"] = [[], []]
aa_CHpair_2Dhist_mdict["LEU"]["CG-HG"] = [[], []]
aa_CHpair_2Dhist_mdict["LEU"]["CD1-HD1"] = [[], []]
aa_CHpair_2Dhist_mdict["LEU"]["CD2-HD2"] = [[], []]

#"LYS" : ["HA", "HB2", "HB3", "HG2", "HG3", "HD2", "HD3", "HE2", "HE3", "CA", "CB", "CG", "CD", "CE"],
aa_CHpair_2Dhist_mdict["LYS"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_mdict["LYS"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_mdict["LYS"]["CB-HB3"] = [[], []]
aa_CHpair_2Dhist_mdict["LYS"]["CG-HG2"] = [[], []]
aa_CHpair_2Dhist_mdict["LYS"]["CG-HG3"] = [[], []]
aa_CHpair_2Dhist_mdict["LYS"]["CD-HD2"] = [[], []]
aa_CHpair_2Dhist_mdict["LYS"]["CD-HD3"] = [[], []]
aa_CHpair_2Dhist_mdict["LYS"]["CE-HE2"] = [[], []]
aa_CHpair_2Dhist_mdict["LYS"]["CE-HE3"] = [[], []]

#"MET" : ["HA", "HB2", "HB3", "HG2", "HG3", "CA", "CB", "CG"],
aa_CHpair_2Dhist_mdict["MET"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_mdict["MET"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_mdict["MET"]["CB-HB3"] = [[], []]
aa_CHpair_2Dhist_mdict["MET"]["CG-HG2"] = [[], []]
aa_CHpair_2Dhist_mdict["MET"]["CG-HG3"] = [[], []]
aa_CHpair_2Dhist_mdict["MET"]["CE-HE"] = [[], []]

#"PHE" : ["HA", "HB2", "HB3", "CA", "CB"],
aa_CHpair_2Dhist_mdict["PHE"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_mdict["PHE"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_mdict["PHE"]["CB-HB3"] = [[], []]
aa_CHpair_2Dhist_mdict["PHE"]["CD1-HD1"] = [[], []]
aa_CHpair_2Dhist_mdict["PHE"]["CE1-HE1"] = [[], []]

#"PRO" : ["HA", "HB2", "HB3", "HG2", "HG3", "HD2", "HD3", "CA", "CB", "CG", "CD"],
aa_CHpair_2Dhist_mdict["PRO"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_mdict["PRO"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_mdict["PRO"]["CB-HB3"] = [[], []]
aa_CHpair_2Dhist_mdict["PRO"]["CG-HG2"] = [[], []]
aa_CHpair_2Dhist_mdict["PRO"]["CG-HG3"] = [[], []]
aa_CHpair_2Dhist_mdict["PRO"]["CD-HD2"] = [[], []]
aa_CHpair_2Dhist_mdict["PRO"]["CD-HD3"] = [[], []]

#"SER" : ["HA", "HB2", "HB3", "CA", "CB"],
aa_CHpair_2Dhist_mdict["SER"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_mdict["SER"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_mdict["SER"]["CB-HB3"] = [[], []]

#"THR" : ["HA", "HB", "HG2", "CA", "CB", "CG2"],
aa_CHpair_2Dhist_mdict["THR"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_mdict["THR"]["CB-HB"] = [[], []]
aa_CHpair_2Dhist_mdict["THR"]["CG2-HG2"] = [[], []]

#"TRP" : ["HA", "HB2", "HB3", "CA", "CB"],
aa_CHpair_2Dhist_mdict["TRP"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_mdict["TRP"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_mdict["TRP"]["CB-HB3"] = [[], []]

#"TYR" : ["HA", "HB2", "HB3", "CA", "CB"],
aa_CHpair_2Dhist_mdict["TYR"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_mdict["TYR"]["CB-HB2"] = [[], []]
aa_CHpair_2Dhist_mdict["TYR"]["CB-HB3"] = [[], []]
aa_CHpair_2Dhist_mdict["TYR"]["CD1-HD1"] = [[], []]
aa_CHpair_2Dhist_mdict["TYR"]["CE1-HE1"] = [[], []]

#"VAL" : ["HA", "HB", "HG1", "HG2", "CA", "CB", "CG1", "CG2"]
aa_CHpair_2Dhist_mdict["VAL"]["CA-HA"] = [[], []]
aa_CHpair_2Dhist_mdict["VAL"]["CB-HB"] = [[], []]
aa_CHpair_2Dhist_mdict["VAL"]["CG1-HG1"] = [[], []]
aa_CHpair_2Dhist_mdict["VAL"]["CG2-HG2"] = [[], []]


# Multidict with percentiles. A 0.9 percentile of value X means that the top 10% of the probability density values were above X. 
percentile_mdict = tree()
percentile_mdict["ALA"]["CA-HA"][0.9]=3.383211e-05 
percentile_mdict["ALA"]["CA-HA"][0.85]=8.607352e-06 
percentile_mdict["ALA"]["CA-HA"][0.80]=3.065289e-06 
percentile_mdict["ALA"]["CA-HA"][0.75]=1.225065e-06 
percentile_mdict["ALA"]["CA-HA"][0.0]=0.0
percentile_mdict["ALA"]["CB-HB"][0.9]=3.638728e-05 
percentile_mdict["ALA"]["CB-HB"][0.85]=1.198696e-05 
percentile_mdict["ALA"]["CB-HB"][0.80]=4.254142e-06 
percentile_mdict["ALA"]["CB-HB"][0.75]=1.401189e-06 
percentile_mdict["ALA"]["CB-HB"][0.0]=0.0
percentile_mdict["ARG"]["CA-HA"][0.9]=5.946551e-05 
percentile_mdict["ARG"]["CA-HA"][0.85]=1.767025e-05 
percentile_mdict["ARG"]["CA-HA"][0.80]=6.98521e-06 
percentile_mdict["ARG"]["CA-HA"][0.75]=2.949357e-06 
percentile_mdict["ARG"]["CA-HA"][0.0]=0.0
percentile_mdict["ARG"]["CB-HB2"][0.9]=3.446836e-05 
percentile_mdict["ARG"]["CB-HB2"][0.85]=1.165121e-05 
percentile_mdict["ARG"]["CB-HB2"][0.80]=5.137525e-06 
percentile_mdict["ARG"]["CB-HB2"][0.75]=2.595404e-06 
percentile_mdict["ARG"]["CB-HB2"][0.0]=0.0
percentile_mdict["ARG"]["CB-HB3"][0.9]=3.446836e-05 
percentile_mdict["ARG"]["CB-HB3"][0.85]=1.165121e-05 
percentile_mdict["ARG"]["CB-HB3"][0.80]=5.137525e-06 
percentile_mdict["ARG"]["CB-HB3"][0.75]=2.595404e-06 
percentile_mdict["ARG"]["CB-HB3"][0.0]=0.0
percentile_mdict["ARG"]["CD-HD2"][0.9]=2.888977e-05 
percentile_mdict["ARG"]["CD-HD2"][0.85]=1.189263e-05 
percentile_mdict["ARG"]["CD-HD2"][0.80]=5.836075e-06 
percentile_mdict["ARG"]["CD-HD2"][0.75]=2.592758e-06 
percentile_mdict["ARG"]["CD-HD2"][0.0]=0.0
percentile_mdict["ARG"]["CD-HD3"][0.9]=2.888977e-05 
percentile_mdict["ARG"]["CD-HD3"][0.85]=1.189263e-05 
percentile_mdict["ARG"]["CD-HD3"][0.80]=5.836075e-06 
percentile_mdict["ARG"]["CD-HD3"][0.75]=2.592758e-06 
percentile_mdict["ARG"]["CD-HD3"][0.0]=0.0
percentile_mdict["ARG"]["CG-HG2"][0.9]=5.934924e-05 
percentile_mdict["ARG"]["CG-HG2"][0.85]=2.220837e-05 
percentile_mdict["ARG"]["CG-HG2"][0.80]=1.165093e-05 
percentile_mdict["ARG"]["CG-HG2"][0.75]=6.073822e-06 
percentile_mdict["ARG"]["CG-HG2"][0.0]=0.0
percentile_mdict["ARG"]["CG-HG3"][0.9]=5.934924e-05 
percentile_mdict["ARG"]["CG-HG3"][0.85]=2.220837e-05 
percentile_mdict["ARG"]["CG-HG3"][0.80]=1.165093e-05 
percentile_mdict["ARG"]["CG-HG3"][0.75]=6.073822e-06 
percentile_mdict["ARG"]["CG-HG3"][0.0]=0.0
percentile_mdict["ASN"]["CA-HA"][0.9]=0.0001172794 
percentile_mdict["ASN"]["CA-HA"][0.85]=4.051621e-05 
percentile_mdict["ASN"]["CA-HA"][0.80]=1.599306e-05 
percentile_mdict["ASN"]["CA-HA"][0.75]=8.039179e-06 
percentile_mdict["ASN"]["CA-HA"][0.0]=0.0
percentile_mdict["ASN"]["CB-HB2"][0.9]=4.639792e-05 
percentile_mdict["ASN"]["CB-HB2"][0.85]=1.781918e-05 
percentile_mdict["ASN"]["CB-HB2"][0.80]=7.606467e-06 
percentile_mdict["ASN"]["CB-HB2"][0.75]=4.229953e-06 
percentile_mdict["ASN"]["CB-HB2"][0.0]=0.0
percentile_mdict["ASN"]["CB-HB3"][0.9]=4.639792e-05 
percentile_mdict["ASN"]["CB-HB3"][0.85]=1.781918e-05 
percentile_mdict["ASN"]["CB-HB3"][0.80]=7.606467e-06 
percentile_mdict["ASN"]["CB-HB3"][0.75]=4.229953e-06 
percentile_mdict["ASN"]["CB-HB3"][0.0]=0.0
percentile_mdict["ASP"]["CA-HA"][0.9]=9.333629e-05 
percentile_mdict["ASP"]["CA-HA"][0.85]=2.681896e-05 
percentile_mdict["ASP"]["CA-HA"][0.80]=1.045079e-05 
percentile_mdict["ASP"]["CA-HA"][0.75]=4.367671e-06 
percentile_mdict["ASP"]["CA-HA"][0.0]=0.0
percentile_mdict["ASP"]["CB-HB2"][0.9]=2.917872e-05 
percentile_mdict["ASP"]["CB-HB2"][0.85]=9.975445e-06 
percentile_mdict["ASP"]["CB-HB2"][0.80]=3.406245e-06 
percentile_mdict["ASP"]["CB-HB2"][0.75]=1.887798e-06 
percentile_mdict["ASP"]["CB-HB2"][0.0]=0.0
percentile_mdict["ASP"]["CB-HB3"][0.9]=2.917872e-05 
percentile_mdict["ASP"]["CB-HB3"][0.85]=9.975445e-06 
percentile_mdict["ASP"]["CB-HB3"][0.80]=3.406245e-06 
percentile_mdict["ASP"]["CB-HB3"][0.75]=1.887798e-06 
percentile_mdict["ASP"]["CB-HB3"][0.0]=0.0
percentile_mdict["CYS"]["CA-HA"][0.9]=0.0001552716 
percentile_mdict["CYS"]["CA-HA"][0.85]=6.162234e-05 
percentile_mdict["CYS"]["CA-HA"][0.80]=2.580545e-05 
percentile_mdict["CYS"]["CA-HA"][0.75]=1.215023e-05 
percentile_mdict["CYS"]["CA-HA"][0.0]=0.0
percentile_mdict["CYS"]["CB-HB2"][0.9]=2.399013e-05 
percentile_mdict["CYS"]["CB-HB2"][0.85]=8.012495e-06 
percentile_mdict["CYS"]["CB-HB2"][0.80]=3.381601e-06 
percentile_mdict["CYS"]["CB-HB2"][0.75]=1.623929e-06 
percentile_mdict["CYS"]["CB-HB2"][0.0]=0.0
percentile_mdict["CYS"]["CB-HB3"][0.9]=2.399013e-05 
percentile_mdict["CYS"]["CB-HB3"][0.85]=8.012495e-06 
percentile_mdict["CYS"]["CB-HB3"][0.80]=3.381601e-06 
percentile_mdict["CYS"]["CB-HB3"][0.75]=1.623929e-06 
percentile_mdict["CYS"]["CB-HB3"][0.0]=0.0
percentile_mdict["GLN"]["CA-HA"][0.9]=0.0001124071 
percentile_mdict["GLN"]["CA-HA"][0.85]=3.525829e-05 
percentile_mdict["GLN"]["CA-HA"][0.80]=1.359177e-05 
percentile_mdict["GLN"]["CA-HA"][0.75]=6.015156e-06 
percentile_mdict["GLN"]["CA-HA"][0.0]=0.0
percentile_mdict["GLN"]["CB-HB2"][0.9]=8.753691e-05 
percentile_mdict["GLN"]["CB-HB2"][0.85]=3.10272e-05 
percentile_mdict["GLN"]["CB-HB2"][0.80]=1.421064e-05 
percentile_mdict["GLN"]["CB-HB2"][0.75]=7.209848e-06 
percentile_mdict["GLN"]["CB-HB2"][0.0]=0.0
percentile_mdict["GLN"]["CB-HB3"][0.9]=8.753691e-05 
percentile_mdict["GLN"]["CB-HB3"][0.85]=3.10272e-05 
percentile_mdict["GLN"]["CB-HB3"][0.80]=1.421064e-05 
percentile_mdict["GLN"]["CB-HB3"][0.75]=7.209848e-06 
percentile_mdict["GLN"]["CB-HB3"][0.0]=0.0
percentile_mdict["GLN"]["CG-HG2"][0.9]=5.680763e-05 
percentile_mdict["GLN"]["CG-HG2"][0.85]=2.375643e-05 
percentile_mdict["GLN"]["CG-HG2"][0.80]=1.235801e-05 
percentile_mdict["GLN"]["CG-HG2"][0.75]=7.074422e-06 
percentile_mdict["GLN"]["CG-HG2"][0.0]=0.0
percentile_mdict["GLN"]["CG-HG3"][0.9]=5.680763e-05 
percentile_mdict["GLN"]["CG-HG3"][0.85]=2.375643e-05 
percentile_mdict["GLN"]["CG-HG3"][0.80]=1.235801e-05 
percentile_mdict["GLN"]["CG-HG3"][0.75]=7.074422e-06 
percentile_mdict["GLN"]["CG-HG3"][0.0]=0.0
percentile_mdict["GLU"]["CA-HA"][0.9]=1.231792e-05 
percentile_mdict["GLU"]["CA-HA"][0.85]=2.976818e-06 
percentile_mdict["GLU"]["CA-HA"][0.80]=6.918509e-07 
percentile_mdict["GLU"]["CA-HA"][0.75]=4.008398e-08 
percentile_mdict["GLU"]["CA-HA"][0.0]=0.0
percentile_mdict["GLU"]["CB-HB2"][0.9]=3.069784e-05 
percentile_mdict["GLU"]["CB-HB2"][0.85]=1.003659e-05 
percentile_mdict["GLU"]["CB-HB2"][0.80]=4.061968e-06 
percentile_mdict["GLU"]["CB-HB2"][0.75]=1.927056e-06 
percentile_mdict["GLU"]["CB-HB2"][0.0]=0.0
percentile_mdict["GLU"]["CB-HB3"][0.9]=3.069784e-05 
percentile_mdict["GLU"]["CB-HB3"][0.85]=1.003659e-05 
percentile_mdict["GLU"]["CB-HB3"][0.80]=4.061968e-06 
percentile_mdict["GLU"]["CB-HB3"][0.75]=1.927056e-06 
percentile_mdict["GLU"]["CB-HB3"][0.0]=0.0
percentile_mdict["GLU"]["CG-HG2"][0.9]=5.575773e-05 
percentile_mdict["GLU"]["CG-HG2"][0.85]=1.930426e-05 
percentile_mdict["GLU"]["CG-HG2"][0.80]=9.536563e-06 
percentile_mdict["GLU"]["CG-HG2"][0.75]=5.376753e-06 
percentile_mdict["GLU"]["CG-HG2"][0.0]=0.0
percentile_mdict["GLU"]["CG-HG3"][0.9]=5.575773e-05 
percentile_mdict["GLU"]["CG-HG3"][0.85]=1.930426e-05 
percentile_mdict["GLU"]["CG-HG3"][0.80]=9.536563e-06 
percentile_mdict["GLU"]["CG-HG3"][0.75]=5.376753e-06 
percentile_mdict["GLU"]["CG-HG3"][0.0]=0.0
percentile_mdict["GLY"]["CA-HA2"][0.9]=3.082408e-05 
percentile_mdict["GLY"]["CA-HA2"][0.85]=9.641981e-06 
percentile_mdict["GLY"]["CA-HA2"][0.80]=3.752225e-06 
percentile_mdict["GLY"]["CA-HA2"][0.75]=1.683045e-06 
percentile_mdict["GLY"]["CA-HA2"][0.0]=0.0
percentile_mdict["GLY"]["CA-HA3"][0.9]=3.082408e-05 
percentile_mdict["GLY"]["CA-HA3"][0.85]=9.641981e-06 
percentile_mdict["GLY"]["CA-HA3"][0.80]=3.752225e-06 
percentile_mdict["GLY"]["CA-HA3"][0.75]=1.683045e-06 
percentile_mdict["GLY"]["CA-HA3"][0.0]=0.0
percentile_mdict["HIS"]["CA-HA"][0.9]=7.246005e-05 
percentile_mdict["HIS"]["CA-HA"][0.85]=2.619243e-05 
percentile_mdict["HIS"]["CA-HA"][0.80]=1.215831e-05 
percentile_mdict["HIS"]["CA-HA"][0.75]=5.907301e-06 
percentile_mdict["HIS"]["CA-HA"][0.0]=0.0
percentile_mdict["HIS"]["CB-HB2"][0.9]=1.3788e-05 
percentile_mdict["HIS"]["CB-HB2"][0.85]=4.902179e-06 
percentile_mdict["HIS"]["CB-HB2"][0.80]=1.340235e-06 
percentile_mdict["HIS"]["CB-HB2"][0.75]=9.781946e-08 
percentile_mdict["HIS"]["CB-HB2"][0.0]=0.0
percentile_mdict["HIS"]["CB-HB3"][0.9]=1.3788e-05 
percentile_mdict["HIS"]["CB-HB3"][0.85]=4.902179e-06 
percentile_mdict["HIS"]["CB-HB3"][0.80]=1.340235e-06 
percentile_mdict["HIS"]["CB-HB3"][0.75]=9.781946e-08 
percentile_mdict["HIS"]["CB-HB3"][0.0]=0.0
percentile_mdict["ILE"]["CA-HA"][0.9]=9.534934e-05 
percentile_mdict["ILE"]["CA-HA"][0.85]=2.466723e-05 
percentile_mdict["ILE"]["CA-HA"][0.80]=7.983614e-06 
percentile_mdict["ILE"]["CA-HA"][0.75]=3.294043e-06 
percentile_mdict["ILE"]["CA-HA"][0.0]=0.0
percentile_mdict["ILE"]["CB-HB"][0.9]=6.020752e-05 
percentile_mdict["ILE"]["CB-HB"][0.85]=1.662912e-05 
percentile_mdict["ILE"]["CB-HB"][0.80]=6.640653e-06 
percentile_mdict["ILE"]["CB-HB"][0.75]=3.263402e-06 
percentile_mdict["ILE"]["CB-HB"][0.0]=0.0
percentile_mdict["ILE"]["CD1-HD1"][0.9]=2.902411e-05 
percentile_mdict["ILE"]["CD1-HD1"][0.85]=6.746207e-06 
percentile_mdict["ILE"]["CD1-HD1"][0.80]=2.208233e-06 
percentile_mdict["ILE"]["CD1-HD1"][0.75]=4.640454e-07 
percentile_mdict["ILE"]["CD1-HD1"][0.0]=0.0
percentile_mdict["ILE"]["CG1-HG12"][0.9]=6.095986e-05 
percentile_mdict["ILE"]["CG1-HG12"][0.85]=1.66325e-05 
percentile_mdict["ILE"]["CG1-HG12"][0.80]=7.616313e-06 
percentile_mdict["ILE"]["CG1-HG12"][0.75]=4.219306e-06 
percentile_mdict["ILE"]["CG1-HG12"][0.0]=0.0
percentile_mdict["ILE"]["CG1-HG13"][0.9]=6.095986e-05 
percentile_mdict["ILE"]["CG1-HG13"][0.85]=1.66325e-05 
percentile_mdict["ILE"]["CG1-HG13"][0.80]=7.616313e-06 
percentile_mdict["ILE"]["CG1-HG13"][0.75]=4.219306e-06 
percentile_mdict["ILE"]["CG1-HG13"][0.0]=0.0
percentile_mdict["ILE"]["CG2-HG2"][0.9]=8.220756e-05 
percentile_mdict["ILE"]["CG2-HG2"][0.85]=3.156711e-05 
percentile_mdict["ILE"]["CG2-HG2"][0.80]=1.577987e-05 
percentile_mdict["ILE"]["CG2-HG2"][0.75]=9.59985e-06 
percentile_mdict["ILE"]["CG2-HG2"][0.0]=0.0
percentile_mdict["LEU"]["CA-HA"][0.9]=2.380675e-05 
percentile_mdict["LEU"]["CA-HA"][0.85]=4.332169e-06 
percentile_mdict["LEU"]["CA-HA"][0.80]=1.472931e-06 
percentile_mdict["LEU"]["CA-HA"][0.75]=3.97495e-07 
percentile_mdict["LEU"]["CA-HA"][0.0]=0.0
percentile_mdict["LEU"]["CB-HB2"][0.9]=1.980984e-05 
percentile_mdict["LEU"]["CB-HB2"][0.85]=5.631579e-06 
percentile_mdict["LEU"]["CB-HB2"][0.80]=2.338489e-06 
percentile_mdict["LEU"]["CB-HB2"][0.75]=9.633315e-07 
percentile_mdict["LEU"]["CB-HB2"][0.0]=0.0
percentile_mdict["LEU"]["CB-HB3"][0.9]=1.980984e-05 
percentile_mdict["LEU"]["CB-HB3"][0.85]=5.631579e-06 
percentile_mdict["LEU"]["CB-HB3"][0.80]=2.338489e-06 
percentile_mdict["LEU"]["CB-HB3"][0.75]=9.633315e-07 
percentile_mdict["LEU"]["CB-HB3"][0.0]=0.0
percentile_mdict["LEU"]["CD1-HD1"][0.9]=1.149187e-05 
percentile_mdict["LEU"]["CD1-HD1"][0.85]=2.695695e-06 
percentile_mdict["LEU"]["CD1-HD1"][0.80]=7.594961e-07 
percentile_mdict["LEU"]["CD1-HD1"][0.75]=7.143093e-08 
percentile_mdict["LEU"]["CD1-HD1"][0.0]=0.0
percentile_mdict["LEU"]["CD2-HD2"][0.9]=1.149187e-05 
percentile_mdict["LEU"]["CD2-HD2"][0.85]=2.695695e-06 
percentile_mdict["LEU"]["CD2-HD2"][0.80]=7.594961e-07 
percentile_mdict["LEU"]["CD2-HD2"][0.75]=7.143093e-08 
percentile_mdict["LEU"]["CD2-HD2"][0.0]=0.0
percentile_mdict["LEU"]["CG-HG"][0.9]=4.21333e-05 
percentile_mdict["LEU"]["CG-HG"][0.85]=1.514254e-05 
percentile_mdict["LEU"]["CG-HG"][0.80]=7.091997e-06 
percentile_mdict["LEU"]["CG-HG"][0.75]=3.766562e-06 
percentile_mdict["LEU"]["CG-HG"][0.0]=0.0
percentile_mdict["LYS"]["CA-HA"][0.9]=2.643352e-05 
percentile_mdict["LYS"]["CA-HA"][0.85]=6.437593e-06 
percentile_mdict["LYS"]["CA-HA"][0.80]=2.276942e-06 
percentile_mdict["LYS"]["CA-HA"][0.75]=7.757168e-07 
percentile_mdict["LYS"]["CA-HA"][0.0]=0.0
percentile_mdict["LYS"]["CB-HB2"][0.9]=4.664275e-05 
percentile_mdict["LYS"]["CB-HB2"][0.85]=1.53263e-05 
percentile_mdict["LYS"]["CB-HB2"][0.80]=6.457721e-06 
percentile_mdict["LYS"]["CB-HB2"][0.75]=2.966778e-06 
percentile_mdict["LYS"]["CB-HB2"][0.0]=0.0
percentile_mdict["LYS"]["CB-HB3"][0.9]=4.664275e-05 
percentile_mdict["LYS"]["CB-HB3"][0.85]=1.53263e-05 
percentile_mdict["LYS"]["CB-HB3"][0.80]=6.457721e-06 
percentile_mdict["LYS"]["CB-HB3"][0.75]=2.966778e-06 
percentile_mdict["LYS"]["CB-HB3"][0.0]=0.0
percentile_mdict["LYS"]["CD-HD2"][0.9]=1.878168e-05 
percentile_mdict["LYS"]["CD-HD2"][0.85]=6.637545e-06 
percentile_mdict["LYS"]["CD-HD2"][0.80]=2.318034e-06 
percentile_mdict["LYS"]["CD-HD2"][0.75]=4.65699e-07 
percentile_mdict["LYS"]["CD-HD2"][0.0]=0.0
percentile_mdict["LYS"]["CD-HD3"][0.9]=1.878168e-05 
percentile_mdict["LYS"]["CD-HD3"][0.85]=6.637545e-06 
percentile_mdict["LYS"]["CD-HD3"][0.80]=2.318034e-06 
percentile_mdict["LYS"]["CD-HD3"][0.75]=4.65699e-07 
percentile_mdict["LYS"]["CD-HD3"][0.0]=0.0
percentile_mdict["LYS"]["CE-HE2"][0.9]=2.17513e-05 
percentile_mdict["LYS"]["CE-HE2"][0.85]=1.013308e-05 
percentile_mdict["LYS"]["CE-HE2"][0.80]=3.919163e-06 
percentile_mdict["LYS"]["CE-HE2"][0.75]=9.060029e-07 
percentile_mdict["LYS"]["CE-HE2"][0.0]=0.0
percentile_mdict["LYS"]["CE-HE3"][0.9]=2.17513e-05 
percentile_mdict["LYS"]["CE-HE3"][0.85]=1.013308e-05 
percentile_mdict["LYS"]["CE-HE3"][0.80]=3.919163e-06 
percentile_mdict["LYS"]["CE-HE3"][0.75]=9.060029e-07 
percentile_mdict["LYS"]["CE-HE3"][0.0]=0.0
percentile_mdict["LYS"]["CG-HG2"][0.9]=2.256604e-05 
percentile_mdict["LYS"]["CG-HG2"][0.85]=5.884895e-06 
percentile_mdict["LYS"]["CG-HG2"][0.80]=2.36138e-06 
percentile_mdict["LYS"]["CG-HG2"][0.75]=5.55815e-07 
percentile_mdict["LYS"]["CG-HG2"][0.0]=0.0
percentile_mdict["LYS"]["CG-HG3"][0.9]=2.256604e-05 
percentile_mdict["LYS"]["CG-HG3"][0.85]=5.884895e-06 
percentile_mdict["LYS"]["CG-HG3"][0.80]=2.36138e-06 
percentile_mdict["LYS"]["CG-HG3"][0.75]=5.55815e-07 
percentile_mdict["LYS"]["CG-HG3"][0.0]=0.0
percentile_mdict["MET"]["CA-HA"][0.9]=8.288802e-05 
percentile_mdict["MET"]["CA-HA"][0.85]=2.389277e-05 
percentile_mdict["MET"]["CA-HA"][0.80]=9.225335e-06 
percentile_mdict["MET"]["CA-HA"][0.75]=4.012606e-06 
percentile_mdict["MET"]["CA-HA"][0.0]=0.0
percentile_mdict["MET"]["CB-HB2"][0.9]=2.178438e-05 
percentile_mdict["MET"]["CB-HB2"][0.85]=7.360763e-06 
percentile_mdict["MET"]["CB-HB2"][0.80]=2.740282e-06 
percentile_mdict["MET"]["CB-HB2"][0.75]=1.010214e-06 
percentile_mdict["MET"]["CB-HB2"][0.0]=0.0
percentile_mdict["MET"]["CB-HB3"][0.9]=2.178438e-05 
percentile_mdict["MET"]["CB-HB3"][0.85]=7.360763e-06 
percentile_mdict["MET"]["CB-HB3"][0.80]=2.740282e-06 
percentile_mdict["MET"]["CB-HB3"][0.75]=1.010214e-06 
percentile_mdict["MET"]["CB-HB3"][0.0]=0.0
percentile_mdict["MET"]["CE-HE"][0.9]=1.087333e-05 
percentile_mdict["MET"]["CE-HE"][0.85]=4.045831e-06 
percentile_mdict["MET"]["CE-HE"][0.80]=6.60065e-07 
percentile_mdict["MET"]["CE-HE"][0.75]=3.457968e-08 
percentile_mdict["MET"]["CE-HE"][0.0]=0.0
percentile_mdict["MET"]["CG-HG2"][0.9]=2.140774e-05 
percentile_mdict["MET"]["CG-HG2"][0.85]=7.389893e-06 
percentile_mdict["MET"]["CG-HG2"][0.80]=2.85931e-06 
percentile_mdict["MET"]["CG-HG2"][0.75]=6.390732e-07 
percentile_mdict["MET"]["CG-HG2"][0.0]=0.0
percentile_mdict["MET"]["CG-HG3"][0.9]=2.140774e-05 
percentile_mdict["MET"]["CG-HG3"][0.85]=7.389893e-06 
percentile_mdict["MET"]["CG-HG3"][0.80]=2.85931e-06 
percentile_mdict["MET"]["CG-HG3"][0.75]=6.390732e-07 
percentile_mdict["MET"]["CG-HG3"][0.0]=0.0
percentile_mdict["PHE"]["CA-HA"][0.9]=0.0001295514 
percentile_mdict["PHE"]["CA-HA"][0.85]=4.904552e-05 
percentile_mdict["PHE"]["CA-HA"][0.80]=2.120015e-05 
percentile_mdict["PHE"]["CA-HA"][0.75]=8.560579e-06 
percentile_mdict["PHE"]["CA-HA"][0.0]=0.0
percentile_mdict["PHE"]["CB-HB2"][0.9]=0.0001203751 
percentile_mdict["PHE"]["CB-HB2"][0.85]=4.12753e-05 
percentile_mdict["PHE"]["CB-HB2"][0.80]=1.826156e-05 
percentile_mdict["PHE"]["CB-HB2"][0.75]=9.439478e-06 
percentile_mdict["PHE"]["CB-HB2"][0.0]=0.0
percentile_mdict["PHE"]["CB-HB3"][0.9]=0.0001203751 
percentile_mdict["PHE"]["CB-HB3"][0.85]=4.12753e-05 
percentile_mdict["PHE"]["CB-HB3"][0.80]=1.826156e-05 
percentile_mdict["PHE"]["CB-HB3"][0.75]=9.439478e-06 
percentile_mdict["PHE"]["CB-HB3"][0.0]=0.0
percentile_mdict["PHE"]["CD1-HD1"][0.9]=6.019081e-06 
percentile_mdict["PHE"]["CD1-HD1"][0.85]=2.660154e-07 
percentile_mdict["PHE"]["CD1-HD1"][0.80]=9.817369e-10 
percentile_mdict["PHE"]["CD1-HD1"][0.75]=3.022086e-13 
percentile_mdict["PHE"]["CD1-HD1"][0.0]=0.0
percentile_mdict["PHE"]["CD2-HD2"][0.9]=6.019081e-06 
percentile_mdict["PHE"]["CD2-HD2"][0.85]=2.660154e-07 
percentile_mdict["PHE"]["CD2-HD2"][0.80]=9.817369e-10 
percentile_mdict["PHE"]["CD2-HD2"][0.75]=3.022086e-13 
percentile_mdict["PHE"]["CD2-HD2"][0.0]=0.0
percentile_mdict["PHE"]["CE1-HE1"][0.9]=7.263719e-05 
percentile_mdict["PHE"]["CE1-HE1"][0.85]=2.079553e-05 
percentile_mdict["PHE"]["CE1-HE1"][0.80]=9.134318e-06 
percentile_mdict["PHE"]["CE1-HE1"][0.75]=5.056216e-06 
percentile_mdict["PHE"]["CE1-HE1"][0.0]=0.0
percentile_mdict["PHE"]["CE2-HE2"][0.9]=7.263719e-05 
percentile_mdict["PHE"]["CE2-HE2"][0.85]=2.079553e-05 
percentile_mdict["PHE"]["CE2-HE2"][0.80]=9.134318e-06 
percentile_mdict["PHE"]["CE2-HE2"][0.75]=5.056216e-06 
percentile_mdict["PHE"]["CE2-HE2"][0.0]=0.0
percentile_mdict["PHE"]["CZ-HZ"][0.9]=6.675602e-05 
percentile_mdict["PHE"]["CZ-HZ"][0.85]=1.897093e-05 
percentile_mdict["PHE"]["CZ-HZ"][0.80]=9.413795e-06 
percentile_mdict["PHE"]["CZ-HZ"][0.75]=5.637111e-06 
percentile_mdict["PHE"]["CZ-HZ"][0.0]=0.0
percentile_mdict["PRO"]["CA-HA"][0.9]=4.682164e-05 
percentile_mdict["PRO"]["CA-HA"][0.85]=1.526172e-05 
percentile_mdict["PRO"]["CA-HA"][0.80]=6.948323e-06 
percentile_mdict["PRO"]["CA-HA"][0.75]=3.702431e-06 
percentile_mdict["PRO"]["CA-HA"][0.0]=0.0
percentile_mdict["PRO"]["CB-HB2"][0.9]=1.513565e-05 
percentile_mdict["PRO"]["CB-HB2"][0.85]=3.97404e-06 
percentile_mdict["PRO"]["CB-HB2"][0.80]=1.588342e-06 
percentile_mdict["PRO"]["CB-HB2"][0.75]=3.316277e-07 
percentile_mdict["PRO"]["CB-HB2"][0.0]=0.0
percentile_mdict["PRO"]["CB-HB3"][0.9]=1.513565e-05 
percentile_mdict["PRO"]["CB-HB3"][0.85]=3.97404e-06 
percentile_mdict["PRO"]["CB-HB3"][0.80]=1.588342e-06 
percentile_mdict["PRO"]["CB-HB3"][0.75]=3.316277e-07 
percentile_mdict["PRO"]["CB-HB3"][0.0]=0.0
percentile_mdict["PRO"]["CD-HD2"][0.9]=1.009566e-05 
percentile_mdict["PRO"]["CD-HD2"][0.85]=3.636002e-06 
percentile_mdict["PRO"]["CD-HD2"][0.80]=9.729604e-07 
percentile_mdict["PRO"]["CD-HD2"][0.75]=1.20932e-07 
percentile_mdict["PRO"]["CD-HD2"][0.0]=0.0
percentile_mdict["PRO"]["CD-HD3"][0.9]=1.009566e-05 
percentile_mdict["PRO"]["CD-HD3"][0.85]=3.636002e-06 
percentile_mdict["PRO"]["CD-HD3"][0.80]=9.729604e-07 
percentile_mdict["PRO"]["CD-HD3"][0.75]=1.20932e-07 
percentile_mdict["PRO"]["CD-HD3"][0.0]=0.0
percentile_mdict["PRO"]["CG-HG2"][0.9]=3.031702e-05 
percentile_mdict["PRO"]["CG-HG2"][0.85]=1.009437e-05 
percentile_mdict["PRO"]["CG-HG2"][0.80]=4.534183e-06 
percentile_mdict["PRO"]["CG-HG2"][0.75]=1.819586e-06 
percentile_mdict["PRO"]["CG-HG2"][0.0]=0.0
percentile_mdict["PRO"]["CG-HG3"][0.9]=3.031702e-05 
percentile_mdict["PRO"]["CG-HG3"][0.85]=1.009437e-05 
percentile_mdict["PRO"]["CG-HG3"][0.80]=4.534183e-06 
percentile_mdict["PRO"]["CG-HG3"][0.75]=1.819586e-06 
percentile_mdict["PRO"]["CG-HG3"][0.0]=0.0
percentile_mdict["SER"]["CA-HA"][0.9]=8.379703e-05 
percentile_mdict["SER"]["CA-HA"][0.85]=2.83415e-05 
percentile_mdict["SER"]["CA-HA"][0.80]=1.186364e-05 
percentile_mdict["SER"]["CA-HA"][0.75]=5.433984e-06 
percentile_mdict["SER"]["CA-HA"][0.0]=0.0
percentile_mdict["SER"]["CB-HB2"][0.9]=6.764455e-05 
percentile_mdict["SER"]["CB-HB2"][0.85]=2.214627e-05 
percentile_mdict["SER"]["CB-HB2"][0.80]=1.046156e-05 
percentile_mdict["SER"]["CB-HB2"][0.75]=5.96634e-06 
percentile_mdict["SER"]["CB-HB2"][0.0]=0.0
percentile_mdict["SER"]["CB-HB3"][0.9]=6.764455e-05 
percentile_mdict["SER"]["CB-HB3"][0.85]=2.214627e-05 
percentile_mdict["SER"]["CB-HB3"][0.80]=1.046156e-05 
percentile_mdict["SER"]["CB-HB3"][0.75]=5.96634e-06 
percentile_mdict["SER"]["CB-HB3"][0.0]=0.0
percentile_mdict["THR"]["CA-HA"][0.9]=0.0001044913 
percentile_mdict["THR"]["CA-HA"][0.85]=3.080183e-05 
percentile_mdict["THR"]["CA-HA"][0.80]=1.125661e-05 
percentile_mdict["THR"]["CA-HA"][0.75]=4.798165e-06 
percentile_mdict["THR"]["CA-HA"][0.0]=0.0
percentile_mdict["THR"]["CB-HB"][0.9]=6.826693e-07 
percentile_mdict["THR"]["CB-HB"][0.85]=8.216753e-09 
percentile_mdict["THR"]["CB-HB"][0.80]=2.226973e-12 
percentile_mdict["THR"]["CB-HB"][0.75]=1.177657e-17 
percentile_mdict["THR"]["CB-HB"][0.0]=0.0
percentile_mdict["THR"]["CG2-HG2"][0.9]=3.784801e-05 
percentile_mdict["THR"]["CG2-HG2"][0.85]=1.2547e-05 
percentile_mdict["THR"]["CG2-HG2"][0.80]=4.593529e-06 
percentile_mdict["THR"]["CG2-HG2"][0.75]=9.198615e-07 
percentile_mdict["THR"]["CG2-HG2"][0.0]=0.0
percentile_mdict["TRP"]["CA-HA"][0.9]=0.0001443163 
percentile_mdict["TRP"]["CA-HA"][0.85]=5.697479e-05 
percentile_mdict["TRP"]["CA-HA"][0.80]=2.284171e-05 
percentile_mdict["TRP"]["CA-HA"][0.75]=8.139375e-06 
percentile_mdict["TRP"]["CA-HA"][0.0]=0.0
percentile_mdict["TRP"]["CB-HB2"][0.9]=0.0001373969 
percentile_mdict["TRP"]["CB-HB2"][0.85]=5.099959e-05 
percentile_mdict["TRP"]["CB-HB2"][0.80]=2.464142e-05 
percentile_mdict["TRP"]["CB-HB2"][0.75]=1.23831e-05 
percentile_mdict["TRP"]["CB-HB2"][0.0]=0.0
percentile_mdict["TRP"]["CB-HB3"][0.9]=0.0001373969 
percentile_mdict["TRP"]["CB-HB3"][0.85]=5.099959e-05 
percentile_mdict["TRP"]["CB-HB3"][0.80]=2.464142e-05 
percentile_mdict["TRP"]["CB-HB3"][0.75]=1.23831e-05 
percentile_mdict["TRP"]["CB-HB3"][0.0]=0.0
percentile_mdict["TYR"]["CA-HA"][0.9]=9.23798e-05 
percentile_mdict["TYR"]["CA-HA"][0.85]=3.358453e-05 
percentile_mdict["TYR"]["CA-HA"][0.80]=1.309448e-05 
percentile_mdict["TYR"]["CA-HA"][0.75]=5.689913e-06 
percentile_mdict["TYR"]["CA-HA"][0.0]=0.0
percentile_mdict["TYR"]["CB-HB2"][0.9]=0.0001800234 
percentile_mdict["TYR"]["CB-HB2"][0.85]=7.46402e-05 
percentile_mdict["TYR"]["CB-HB2"][0.80]=3.359517e-05 
percentile_mdict["TYR"]["CB-HB2"][0.75]=1.818607e-05 
percentile_mdict["TYR"]["CB-HB2"][0.0]=0.0
percentile_mdict["TYR"]["CB-HB3"][0.9]=0.0001800234 
percentile_mdict["TYR"]["CB-HB3"][0.85]=7.46402e-05 
percentile_mdict["TYR"]["CB-HB3"][0.80]=3.359517e-05 
percentile_mdict["TYR"]["CB-HB3"][0.75]=1.818607e-05 
percentile_mdict["TYR"]["CB-HB3"][0.0]=0.0
percentile_mdict["TYR"]["CD1-HD1"][0.9]=3.695459e-05 
percentile_mdict["TYR"]["CD1-HD1"][0.85]=1.12606e-05 
percentile_mdict["TYR"]["CD1-HD1"][0.80]=4.915059e-06 
percentile_mdict["TYR"]["CD1-HD1"][0.75]=1.396167e-06 
percentile_mdict["TYR"]["CD1-HD1"][0.0]=0.0
percentile_mdict["TYR"]["CD2-HD2"][0.9]=3.695459e-05 
percentile_mdict["TYR"]["CD2-HD2"][0.85]=1.12606e-05 
percentile_mdict["TYR"]["CD2-HD2"][0.80]=4.915059e-06 
percentile_mdict["TYR"]["CD2-HD2"][0.75]=1.396167e-06 
percentile_mdict["TYR"]["CD2-HD2"][0.0]=0.0
percentile_mdict["TYR"]["CE1-HE1"][0.9]=7.804912e-06 
percentile_mdict["TYR"]["CE1-HE1"][0.85]=1.936751e-06 
percentile_mdict["TYR"]["CE1-HE1"][0.80]=3.149482e-07 
percentile_mdict["TYR"]["CE1-HE1"][0.75]=2.825009e-08 
percentile_mdict["TYR"]["CE1-HE1"][0.0]=0.0
percentile_mdict["TYR"]["CE2-HE2"][0.9]=7.804912e-06 
percentile_mdict["TYR"]["CE2-HE2"][0.85]=1.936751e-06 
percentile_mdict["TYR"]["CE2-HE2"][0.80]=3.149482e-07 
percentile_mdict["TYR"]["CE2-HE2"][0.75]=2.825009e-08 
percentile_mdict["TYR"]["CE2-HE2"][0.0]=0.0
percentile_mdict["VAL"]["CA-HA"][0.9]=1.737019e-05 
percentile_mdict["VAL"]["CA-HA"][0.85]=3.981306e-06 
percentile_mdict["VAL"]["CA-HA"][0.80]=1.112401e-06 
percentile_mdict["VAL"]["CA-HA"][0.75]=2.931649e-07 
percentile_mdict["VAL"]["CA-HA"][0.0]=0.0
percentile_mdict["VAL"]["CB-HB"][0.9]=2.457079e-05 
percentile_mdict["VAL"]["CB-HB"][0.85]=8.136396e-06 
percentile_mdict["VAL"]["CB-HB"][0.80]=2.72143e-06 
percentile_mdict["VAL"]["CB-HB"][0.75]=1.044243e-06 
percentile_mdict["VAL"]["CB-HB"][0.0]=0.0
percentile_mdict["VAL"]["CG1-HG1"][0.9]=6.968235e-05 
percentile_mdict["VAL"]["CG1-HG1"][0.85]=1.862415e-05 
percentile_mdict["VAL"]["CG1-HG1"][0.80]=6.260681e-06 
percentile_mdict["VAL"]["CG1-HG1"][0.75]=2.630323e-06 
percentile_mdict["VAL"]["CG1-HG1"][0.0]=0.0
percentile_mdict["VAL"]["CG2-HG2"][0.9]=6.968235e-05 
percentile_mdict["VAL"]["CG2-HG2"][0.85]=1.862415e-05 
percentile_mdict["VAL"]["CG2-HG2"][0.80]=6.260681e-06 
percentile_mdict["VAL"]["CG2-HG2"][0.75]=2.630323e-06 
percentile_mdict["VAL"]["CG2-HG2"][0.0]=0.0


## carbons to exclude from i-i+1 NOESY peak matching and assignment
exclude_from_i_iplus_matching_dict = {}
exclude_from_i_iplus_matching_dict["MET"] = ['CE']
exclude_from_i_iplus_matching_dict["ILE"] = ['CG1']
exclude_from_i_iplus_matching_dict["LYS"] = ['CD', 'CE']
exclude_from_i_iplus_matching_dict["ARG"] = ['CD']
exclude_from_i_iplus_matching_dict["TYR"] = ['CD1', 'CD2', 'CE1', 'CE2']
exclude_from_i_iplus_matching_dict["PHE"] = ['CD1', 'CD2', 'CE1', 'CE2', 'CZ']
