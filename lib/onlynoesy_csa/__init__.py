"""
For an alternative design look at: http://www.qtrac.eu/pyclassmulti.html
"""


import copy

from lib.alignment import Alignment, read_NHassigned_spectrum_file
from lib.probhist import ProbHist_Loader
from lib.global_vars import *
from lib.global_func import *

class HCNH_CSA():

    def __init__(self, args, spectrum_combo):

        self.args = args

        # Allow CE-HE of MET in HCNH assignments (must not be allowed in TOCSY assignments)
        self.allowed_aa_atoms_dict = allowed_aa_atoms_dict
        self.allowed_aa_atoms_dict["MET"].extend(["HE", "CE"])

        # LOAD PROBABILITY HISTOGRAMS
        self.histload = ProbHist_Loader()
        self.histload.load_1Dhistograms()
        self.histload.load_2Dhistograms()

        # Create all alignment files
        self.ali = Alignment(args.NHmap_FILE, args.FIRST_RESIDUE_NUMBER)
        self.ali.create_absolute_AAIGmatches_alignment()
        assert self.ali.spectrum_combo == spectrum_combo, \
            ColorPrint("ERROR: you NHmap file %s starts with N/A at first residue"
                    " (implies TOCSY-HCNH) but you have specified that the spectrum combination is HCNH-HCNH."
                    " If the spectrum combination is indeed HCNH-HCNH, then replace 'N/A' with '-'." %
                       args.NHmap_FILE, "FAIL")
        self.ali.create_absolute_matches_alignment()
        print("DEBUG: self.ali.protein_alignment_list=", self.ali.protein_alignment_list)
        print("DEBUG: self.ali.absolute_AAIGmatches_alignment_list=", self.ali.absolute_AAIGmatches_alignment_list)
        print("DEBUG: self.ali.absolute_matches_alignment_list=", self.ali.absolute_matches_alignment_list)

        self.iter1_allowed_atoms_dict = {}   # resname->list of allowed carbons in iteration1
        self.iter1_allowed_atoms_dict["LEU"] = ["CD1", "CD2", "CG"]
        self.iter1_allowed_atoms_dict["VAL"] = ["CG1", "CG2"]
        self.iter1_allowed_atoms_dict["ILE"] = ["CD1", "CG2"]
        self.iter1_allowed_atoms_dict["THR"] = ["CG2"]
        self.iter1_allowed_atoms_dict["ALA"] = ["CB"]
        # Also methyl is Met CE-QE which must not be included in iteration 1 & 3

        # Create dictionary of intensity threshold for each aa type
        self.I_THRESHOLD = 1000000
        self.aa2ithres_dict = {}
        if type(args.ithresA) == float:    self.aa2ithres_dict["ALA"] = args.ithresA
        if type(args.ithresT) == float:    self.aa2ithres_dict["THR"] = args.ithresT
        if type(args.ithresV) == float:    self.aa2ithres_dict["VAL"] = args.ithresV
        if type(args.ithresI) == float:    self.aa2ithres_dict["ILE"] = args.ithresI
        if type(args.ithresL) == float:    self.aa2ithres_dict["LEU"] = args.ithresL
        # Create dictionary of percentile threshold for each aa type
        self.aa2pthres_dict = {}
        print("DEBUG: args.pthresA=", args.pthresA)
        if type(args.pthresA) == float:    self.aa2pthres_dict["ALA"] = args.pthresA
        if type(args.pthresT) == float:    self.aa2pthres_dict["THR"] = args.pthresT
        if type(args.pthresV) == float:    self.aa2pthres_dict["VAL"] = args.pthresV
        if type(args.pthresI) == float:    self.aa2pthres_dict["ILE"] = args.pthresI
        if type(args.pthresL) == float:    self.aa2pthres_dict["LEU"] = args.pthresL


        # LOAD ASSIGNED TOCSY CONTENTS (NOT APPLICABLE IN THIS SCRIPT)
        self.TOCSY_residue_assignments_dict, \
        self.TOCSY_residue_NHresonances_dict, \
        self.i_to_iminus1_dict, \
        self.patched_residues_list, \
        self.TOCSY_residue_peak_intensity_mdict = {}, {}, {}, {}, {}

        # Add the N-terminus residue in ali.absolute_matches_alignment_list, which does not have N-H
        # try:    # if the N-H of the 2nd residue have been assigned
        #     print "DEBUG: patching the N-terminus residue ", i_to_iminus1_dict[ali.absolute_matches_alignment_list[1]], " to ali.absolute_matches_alignment_list"
        #     ali.absolute_matches_alignment_list[0] = i_to_iminus1_dict[ali.absolute_matches_alignment_list[1]]
        # except IndexError:
        #     pass
        ## TEMPORARILY DEACTIVATE THE UPDATE OF ali.absolute_matches_alignment_list TO SEE IF ALL THE C-TERMINAL RESIDUES ARE FOUND
        # for i in i_to_iminus1_dict.keys():
        #     if i in ali.absolute_matches_alignment_list:
        #         i_index = ali.absolute_matches_alignment_list.index(i)
        #         if i_index > 0 and ali.absolute_matches_alignment_list[i_index-1] in ['-', 'N/A']:
        #             ali.absolute_matches_alignment_list[i_index-1] = i_to_iminus1_dict[i]
        # print "DEBUG: after patching ali.absolute_matches_alignment_list=", ali.absolute_matches_alignment_list


        # # LOAD HCNH FILE CONTENTS
        # HCNH_residue_assignments_dict, HCNH_residue_NHresonances_dict, original_HCNH_peaks_dict, HCNH_residue_peak_intensity_mdict = read_assigned_TOCSY_file(args.HCNH_FILE, 'HCNH', ali.absolute_AAIGmatches_alignment_list, ali.absolute_matches_alignment_list)

        # LOAD HCNH FILE CONTENTS
        self.HCNH_residue_assignments_dict, \
        self.HCNH_residue_NHresonances_dict, \
        original_HCNH_peaks_dict, \
        self.HCNH_residue_peak_intensity_mdict = \
            read_NHassigned_spectrum_file(args.HCNH_FILE, 'HCNH', ali=self.ali)
        # HCNH_residue_assignments_dict: never mutated in ONLYHCNH
        # self.HCNH_residue_NHresonances_dict: never mutated in ONLYHCNH
        # original_HCNH_peaks_dict: never used in ONLYHCNH
        # self.HCNH_residue_peak_intensity_mdict: used but never mutated in ONLYHCNH

        # This is to remember the C-terminal residues
        self.Cterm_residue_set = set()  # TODO: this is never populated but it is used as empty dict. It was
                                        # TODO: populated only in the OBSOLETE function write_Cterm_residues().

        # For proof-reading
        self.aa_equivalentCarbons_dict = {}
        self.aa_equivalentCarbons_dict["LEU"] = [["CD1", "CD2"]]
        self.aa_equivalentCarbons_dict["VAL"] = [["CG1", "CG2"]]

        # The following all the objects needed in all 5 iterations, which are saved into pickled files
        self.matched_HCNH_residue_assignments_dict = {}
        self.checked_residues_set = set()
        self.TOCSY_assigned_residue_list = []
        self.residues_with_written_NH_list = []
        self.clean_resid_assignments_dict = {}
        self.revordered_residue_keys_iter1 = []
        self.revordered_residue_keys_iter2 = []
        self.revordered_residue_keys_iter3 = []
        self.matched_i_to_iplus1_peaks_mdict = tree()
        self.proofread_clean_resid_assignments_dict = {}
        self.user_HCNH_resid_assignments_dict = {} # same as self.user_HCNH_resid_assignments_dict, but with keys resids instead of residues
        self.user_TOCSY_resid_assignments_dict = {}  # same as self.user_TOCSY_resid_assignments_dict, but with keys resids instead of residues
        self.user_VAL_LEU_HCNH_resid_assignments_dict = {}
        self.user_VAL_LEU_TOCSY_resid_assignments_dict = {}

        # KEEP ONLY SELECTED RESIDUES FOR DEBUGGING (IF REQUESTED)
        if args.DEBUG_RESIDUES:
            sel_residues_list = args.DEBUG_RESIDUES.split(",")
            self.keep_selected_residues(sel_residues_list=sel_residues_list)

    def keep_selected_residues(self, sel_residues_list=[]):
        """
        Only for DEBUGGING. Deletes all other residues but the specified ones.

        :param sel_residues_list: list of the form ['V124', 'V125']
        :return:
        """
        for i, residue in enumerate(self.ali.absolute_matches_alignment_list):
            if residue in sel_residues_list:
                continue
            self.ali.absolute_AAIGmatches_alignment_list[i] = "-"
            self.ali.absolute_matches_alignment_list[i] = "-"

        self.HCNH_residue_assignments_dict = {r: self.HCNH_residue_assignments_dict[r] for r in sel_residues_list}
        self.HCNH_residue_NHresonances_dict = {r: self.HCNH_residue_NHresonances_dict[r] for r in sel_residues_list}
        self.HCNH_residue_peak_intensity_mdict = {r: self.HCNH_residue_peak_intensity_mdict[r] for r in sel_residues_list}


    # Imported methods
    from ._onlynoesy_iter1 import onlynoesy_csa_iter1
    from ._onlynoesy_iter2 import onlynoesy_csa_iter2
    from ._proofread_iter12 import proofread_iter12
    from ._onlynoesy_iter3 import onlynoesy_csa_iter3
    from ._onlynoesy_iter4 import onlynoesy_csa_iter4
    from ._proofread_iter34 import proofread_iter34
    from ._onlynoesy_iter5 import onlynoesy_csa_iter5
    from ._onlynoesy_iter5 import proofread_iter5
