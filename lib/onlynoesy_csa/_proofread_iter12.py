import copy
import re

from lib.csa_for_noesy import *
from lib.global_func import *
from lib.global_vars import *


def proofread_iter12(self, iteration):
    """
    same as same as proofread_iter1()

    :param self:
    :return:
    """
    print("STARTING PROOF-READING HCNH CHEMICAL SHIFT ASSIGNMENTS")

    # OLD # LOAD USER  ASSIGNED TOCSY CONTENTS
    # OLD user_TOCSY_residue_assignments_dict, \
    # OLD user_TOCSY_residue_NHresonances_dict, \
    # OLD user_i_to_iminus1_dict, \
    # OLD user_patched_residues_list, \
    # OLD user_TOCSY_residue_peak_intensity_mdict = \
    # OLD    read_assigned_TOCSY_file(self.args.user_TOCSY_FILE,
    # OLD                             "TOCSY",
    # OLD                             self.ali.absolute_AAIGmatches_alignment_list,
    # OLD                             self.ali.absolute_matches_alignment_list)

    # LOAD USER ASSIGNED TOCSY CONTENTS
    user_TOCSY_residue_assignments_dict, \
    user_TOCSY_residue_NHresonances_dict, \
    user_i_to_iminus1_dict, \
    user_patched_residues_list, \
    user_TOCSY_residue_peak_intensity_mdict = \
        read_NHassigned_spectrum_file(self.args.user_TOCSY_FILE, "TOCSY", ali=self.ali)

    # print "DEBUG: i_to_iminus1_dict=", i_to_iminus1_dict
    print("DEBUG: user_TOCSY_residue_assignments_dict=", user_TOCSY_residue_assignments_dict)
    print("DEBUG: user_TOCSY_residue_NHresonances_dict=", user_TOCSY_residue_NHresonances_dict)
    for k, v in list(user_TOCSY_residue_assignments_dict.items()):
        print(k, "-->", v)

    # OLD # LOAD USER ASSIGNED HCNH FILE CONTENTS
    # OLD user_HCNH_residue_assignments_dict, user_HCNH_residue_NHresonances_dict, user_original_HCNH_peaks_dict, user_HCNH_residue_peak_instensity_dict = read_assigned_TOCSY_file(self.args.user_HCNH_FILE, 'HCNH', self.ali.absolute_AAIGmatches_alignment_list, self.ali.absolute_matches_alignment_list)

    # LOAD USER ASSIGNED HCNH FILE CONTENTS
    user_HCNH_residue_assignments_dict, \
    user_HCNH_residue_NHresonances_dict, \
    user_original_HCNH_peaks_dict, \
    user_HCNH_residue_peak_instensity_dict = \
        read_NHassigned_spectrum_file(self.args.user_HCNH_FILE,
                                      spectrum_combo='HCNH',
                                      ali=self.ali)

    print("DEBUG: user_HCNH_residue_assignments_dict=", user_HCNH_residue_assignments_dict)
    print("DEBUG: user_HCNH_residue_NHresonances_dict=", user_HCNH_residue_NHresonances_dict)
    print("DEBUG: user_original_HCNH_peaks_dict=", user_original_HCNH_peaks_dict)
    # rename TYR and PHE CD-QD>CD1-HD1 and CE-QE->CE1-HE1
    for residue in user_HCNH_residue_assignments_dict:
        if residue[0] in ['Y', 'F']:
            for assignment in user_HCNH_residue_assignments_dict[residue]:
                if assignment[1] == "CD":
                    assignment[1] = "CD1"
                elif assignment[1] == "QD":
                    assignment[1] = "HD1"
                elif assignment[1] == "CE":
                    assignment[1] = "CE1"
                elif assignment[1] == "QE":
                    assignment[1] = "HE1"
    for residue in user_TOCSY_residue_assignments_dict:
        if residue[0] in ['Y', 'F']:
            for assignment in user_TOCSY_residue_assignments_dict[residue]:
                if assignment[1] == "CD":
                    assignment[1] = "CD1"
                elif assignment[1] == "QD":
                    assignment[1] = "HD1"
                elif assignment[1] == "CE":
                    assignment[1] = "CE1"
                elif assignment[1] == "QE":
                    assignment[1] = "HE1"
    for k, v in list(user_HCNH_residue_assignments_dict.items()):
        print(k, "-->", v)

    # convert the keys from residues to resids to be compatible with self.clean_resid_assignments_dict
    for k, v in list(user_HCNH_residue_assignments_dict.items()):
        if not is_valid_residue(k, self.ali.protein_alignment_list, self.args, spectrum_combo='HCNH-HCNH'):
            print("WARNING: omitting non-valid user assigned residue ", k, "!")
            continue
        self.user_HCNH_resid_assignments_dict[int(k[1:])] = v
    for k, v in list(user_TOCSY_residue_assignments_dict.items()):
        if not is_valid_residue(k, self.ali.protein_alignment_list, self.args, spectrum_combo='HCNH-HCNH'):
            print("WARNING: omitting non-valid user assigned residue ", k, "!")
            continue
        self.user_TOCSY_resid_assignments_dict[int(k[1:])] = v

    print("DEBUG: self.user_HCNH_resid_assignments_dict=", self.user_HCNH_resid_assignments_dict)
    print("DEBUG: ITERATION %s point 2 self.clean_resid_assignments_dict= %s" %
          (iteration, self.clean_resid_assignments_dict))

    # CREATE ONLY ASSIGNED VAL AND LEU EQUIVALENT CARBONS DICTIONARIES
    self.user_VAL_LEU_HCNH_resid_assignments_dict, \
    self.user_VAL_LEU_TOCSY_resid_assignments_dict = \
        create_equivalentCarbons_assignments_dict(self.user_HCNH_resid_assignments_dict,
                                                  self.user_TOCSY_resid_assignments_dict,
                                                  self.ali.protein_alignment_list,
                                                  self.ali.absolute_matches_alignment_list,
                                                  self.aa_equivalentCarbons_dict,
                                                  self.args)
    print("DEBUG: self.user_VAL_LEU_HCNH_resid_assignments_dict=", self.user_VAL_LEU_HCNH_resid_assignments_dict)
    print("DEBUG: self.user_VAL_LEU_TOCSY_resid_assignments_dict=", self.user_VAL_LEU_TOCSY_resid_assignments_dict)

    ##
    ##  COMPARE YOUR CARBON ASSIGNMENTS WITH THE USER ASSIGNMENTS
    ##
    ## Remember that the comments can be the following:
    # HCNH
    # HCNH (C-term hanging residue)
    # HCNH + TOCSY-HCNH
    # TOCSY (average)
    # TOCSY-HCNH
    # TOCSY-HCNH + HCNH
    # TOCSY-HCNH + TOCSY (unmatched)
    # TOCSY (unmatched)
    # TOCSY (unmatched) + TOCSY-HCNH
    # Recall that if we have 2 peaks in the same C-group but only one comes from HCNH ('TOCSY (unmatched) + TOCSY-HCNH', or
    # 'TOCSY-HCNH + TOCSY (unmatched)'), then we save only the HCNH Carbon resonance. 'TOCSY (unmatched)' -> only one resonance available
    # only in TOCSY, 'TOCSY (average)' -> 2 resonances available only from TOCSY and they were averaged

    self.proofread_clean_resid_assignments_dict = copy.deepcopy(
        self.clean_resid_assignments_dict)  # same as self.clean_resid_assignments_dict but with comments CORRECT or WRONG at the carbon lines
    for resid in list(self.proofread_clean_resid_assignments_dict.keys()):
        print("DEBUG ITERATION%s: proof reading resid= %i" %(iteration, resid))
        resname = self.proofread_clean_resid_assignments_dict[resid][0][6]
        if resname in list(self.aa_equivalentCarbons_dict.keys()):  # If this is a LEU or VAL
            Carbon_pair = self.aa_equivalentCarbons_dict[resname][0]
            if len([a for a in self.proofread_clean_resid_assignments_dict[resid] if
                    a[3] in Carbon_pair]) > 0:  # AND if the methyl carbons have been assigned
                Proton_pair = ["H" + C[1:] for C in Carbon_pair]  # protons of the equivalent Carbons
                # get only the methyl Carbon assignments (LEU --> CD1,CD2 or VAL --> CG1,CG2) made by the program, to compare them with the user made assignments
                program_assignments_list = [assignment for assignment in
                                            self.proofread_clean_resid_assignments_dict[resid] if assignment[
                                                3] in Carbon_pair + Proton_pair]  # keep only the resonances assigned by the program to the equivalent methyl carbons and protons
                # proof_read_equivalentCarbons(program_assignments_list, resid, self.user_VAL_LEU_HCNH_resid_assignments_dict, self.user_VAL_LEU_TOCSY_resid_assignments_dict, Carbon_pair)   # OBSOLETE!!! it doesn't matter if 1 or 2 of the peaks have been proof-read in iteration 1 & 2
                proof_read_equivalentMethyls(program_assignments_list, resid,
                                             self.user_VAL_LEU_HCNH_resid_assignments_dict,
                                             self.user_VAL_LEU_TOCSY_resid_assignments_dict,
                                             Carbon_pair,
                                             Proton_pair)  # it doesn't matter if 1 or 2 of the peaks have been proof-read in iteration 1 & 2
        # THEN PROCEED NORMALLY
        for assignment in self.proofread_clean_resid_assignments_dict[resid]:
            print("DEBUG ITERATION%s: proof reading assignment=%s" % (iteration, assignment))
            spectrum_combo = re.sub(r' iter[1-9]$', '', assignment[
                7])  # remove the "iter[1-9]" suffix to be able to distinguish the type of assignment (HCNH, HCNH+TOCSY, etc.)
            comment = None
            aa_type = assignment[6]
            nucleus = assignment[3]  # the nucleus name
            if aa_type in list(self.aa_equivalentCarbons_dict.keys()) and nucleus in \
                    self.aa_equivalentCarbons_dict[aa_type][
                        0]:  # if this peak has been assigned by the program to LEU or VAL CD1, CD2 or CG1,CG2 respectively, skip it
                continue

            # IF THIS CARBON IS A METHYL OR HAS ONLY ONE PROTON, CONSIDE ALSO THE PROTON RESONANCE FOR PROOFREADING
            if aa_type in list(aatype_carbon_nongeminalHname_mdict.keys()) and nucleus in list(
                    aatype_carbon_nongeminalHname_mdict[aa_type].keys()):
                # get the proton assignment made by the program
                C_assignment = assignment  # the program-made carbon assignment
                H_assignment = None  # default value for the program-made proton assignment
                print("DEBUG: resid=", resid, "self.proofread_clean_resid_assignments_dict[resid]=",
                      self.proofread_clean_resid_assignments_dict[resid])
                for a in self.proofread_clean_resid_assignments_dict[resid]:
                    if ((type(aatype_carbon_nongeminalHname_mdict[aa_type][nucleus]) == list and a[3] in
                         aatype_carbon_nongeminalHname_mdict[aa_type][nucleus]) or
                            (type(aatype_carbon_nongeminalHname_mdict[aa_type][nucleus]) == str and a[3] ==
                             aatype_carbon_nongeminalHname_mdict[aa_type][nucleus])):
                        H_assignment = a
                        break
                print("DEBUG: C_assignment=", C_assignment, "H_assignment=", H_assignment)
                comment = ""
                if spectrum_combo in ['HCNH', 'HCNH (C-term hanging residue)', 'HCNH + TOCSY-HCNH', 'TOCSY-HCNH',
                                     'TOCSY-HCNH + HCNH', 'TOCSY (unmatched) + TOCSY-HCNH',
                                     'TOCSY-HCNH + TOCSY (unmatched)']:  # if this resonance comes only from HCNH
                    user_assignment = get_Carbon_resonance(resid,
                                                           nucleus,
                                                           spectrum_combo='HCNH',
                                                           user_HCNH_resid_assignments_dict=self.user_HCNH_resid_assignments_dict,
                                                           user_TOCSY_resid_assignments_dict=self.user_TOCSY_resid_assignments_dict,
                                                           aa_equivalentCarbons_dict=self.aa_equivalentCarbons_dict,
                                                           protein_alignment_list=self.ali.protein_alignment_list,
                                                           absolute_matches_alignment_list=self.ali.absolute_matches_alignment_list,
                                                           args=self.args,
                                                           get_assignment=True)
                    if user_assignment == False:  # if no user assignment was made for this Carbon, add no extra comment
                        print("WARNING: no user assignment made in HCNH for resid", resid, "and Carbon", nucleus,
                              ". Checking TOCSY...")
                        user_assignment = get_Carbon_resonance(resid, nucleus, 'TOCSY',
                                                               self.user_HCNH_resid_assignments_dict,
                                                               self.user_TOCSY_resid_assignments_dict,
                                                               self.aa_equivalentCarbons_dict,
                                                               self.ali.protein_alignment_list,
                                                               self.ali.absolute_matches_alignment_list, self.args,
                                                               get_assignment=True)
                        if user_assignment == False:
                            print("WARNING: no user assignment made also in TOCSY for resid", resid, "and Carbon",
                                  nucleus, ". Skipping nucleus.")
                            assignment[7] += " <NOT ASSIGNED>"
                            # self.proofread_clean_resid_assignments_dict[resid].append(assignment)
                            continue
                    user_Creson = user_assignment[2]
                    user_Hreson = user_assignment[4]
                    print("DEBUG: user_Creson=", user_Creson, "user_Hreson=", user_Hreson)
                    if approx_equal(C_assignment[1], user_Creson, 0.4) and approx_equal(H_assignment[1],
                                                                                        user_Hreson,
                                                                                        0.04):  # if the user-made assignment agrees with ours, add to the comment "CORRECT"
                        print("DEBUG: CORRECT C_assignment found:", C_assignment)
                        print("DEBUG: CORRECT H_assignment found:", H_assignment)
                        comment = " <CORRECT>"
                    else:
                        print("DEBUG: WRONG C_assignment found:", C_assignment)
                        print("DEBUG: WRONG H_assignment found:", H_assignment)
                        comment = " <WRONG>"
                elif spectrum_combo in ['TOCSY (average)', 'TOCSY (unmatched)']:
                    user_assignment = get_Carbon_resonance(resid, nucleus, 'TOCSY',
                                                           self.user_HCNH_resid_assignments_dict,
                                                           self.user_TOCSY_resid_assignments_dict,
                                                           self.aa_equivalentCarbons_dict,
                                                           self.ali.protein_alignment_list,
                                                           self.ali.absolute_matches_alignment_list, self.args,
                                                           get_assignment=True)
                    if user_assignment == False:  # if no user assignment was made for this Carbon, add no extra comment
                        print("WARNING: no user assignment made in TOCSY for resid", resid, "and Carbon", nucleus,
                              ". Checking HCNH...")
                        user_assignment = get_Carbon_resonance(resid,
                                                               nucleus,
                                                               spectrum_combo='HCNH',
                                                               user_HCNH_resid_assignments_dict=self.user_HCNH_resid_assignments_dict,
                                                               user_TOCSY_resid_assignments_dict=self.user_TOCSY_resid_assignments_dict,
                                                               aa_equivalentCarbons_dict=self.aa_equivalentCarbons_dict,
                                                               protein_alignment_list=self.ali.protein_alignment_list,
                                                               absolute_matches_alignment_list=self.ali.absolute_matches_alignment_list,
                                                               args=self.args,
                                                               get_assignment=True)
                        if user_assignment == False:
                            print("WARNING: no user assignment made also in HCNH for resid", resid, "and Carbon",
                                  nucleus, ". Skipping nucleus.")
                            assignment[7] += " <NOT ASSIGNED>"
                            # self.proofread_clean_resid_assignments_dict[resid].append(assignment)
                            continue
                    user_Creson = user_assignment[2]
                    user_Hreson = user_assignment[4]
                    if approx_equal(C_assignment[1], user_Creson, 0.4) and approx_equal(H_assignment[1],
                                                                                        user_Hreson,
                                                                                        0.04):  # if the user-made assignment agrees with ours, add to the comment "CORRECT"
                        print("DEBUG: CORRECT C_assignment found:", C_assignment)
                        print("DEBUG: CORRECT H_assignment found:", H_assignment)
                        comment = " <CORRECT>"
                    else:
                        print("DEBUG: WRONG C_assignment found:", C_assignment)
                        print("DEBUG: WRONG H_assignment found:", H_assignment)
                        comment = " <WRONG>"

                if comment == "":
                    comment = " <NOT ASSIGNED>"
                C_assignment[7] += comment
                H_assignment[7] += comment
                comment = None  # a trick to add no more comments to "assignment" at the end of this block of code

            elif nucleus[
                0] == 'C':  # only if the assignment line belongs to a Carbon then check if the user has made the same assignment
                comment = ""
                if spectrum_combo in ['HCNH', 'HCNH (C-term hanging residue)', 'HCNH + TOCSY-HCNH', 'TOCSY-HCNH',
                                     'TOCSY-HCNH + HCNH', 'TOCSY (unmatched) + TOCSY-HCNH',
                                     'TOCSY-HCNH + TOCSY (unmatched)']:  # if this resonance comes only from HCNH
                    user_assignment = get_Carbon_resonance(resid,
                                                           nucleus,
                                                           spectrum_combo='HCNH',
                                                           user_HCNH_resid_assignments_dict=self.user_HCNH_resid_assignments_dict,
                                                           user_TOCSY_resid_assignments_dict=self.user_TOCSY_resid_assignments_dict,
                                                           aa_equivalentCarbons_dict=self.aa_equivalentCarbons_dict,
                                                           protein_alignment_list=self.ali.protein_alignment_list,
                                                           absolute_matches_alignment_list=self.ali.absolute_matches_alignment_list,
                                                           args=self.args)
                    if user_assignment == False:  # if no user assignment was made for this Carbon, add no extra comment
                        print("WARNING: no user assignment made in HCNH for resid", resid, "and Carbon", nucleus,
                              ". Checking TOCSY...")
                        user_assignment = get_Carbon_resonance(resid, nucleus, 'TOCSY',
                                                               self.user_HCNH_resid_assignments_dict,
                                                               self.user_TOCSY_resid_assignments_dict,
                                                               self.aa_equivalentCarbons_dict,
                                                               self.ali.protein_alignment_list,
                                                               self.ali.absolute_matches_alignment_list, self.args)
                        if user_assignment == False:
                            print("WARNING: no user assignment made also in TOCSY for resid", resid, "and Carbon",
                                  nucleus, ". Skipping nucleus.")
                            assignment[7] += " <NOT ASSIGNED>"
                            # self.proofread_clean_resid_assignments_dict[resid].append(assignment)
                            continue
                    if approx_equal(assignment[1], user_assignment,
                                    0.4):  # if the user-made assignment agrees with ours, add to the comment "CORRECT"
                        print("DEBUG: CORRECT assignment found:", assignment)
                        comment = " <CORRECT>"
                    else:
                        print("DEBUG: WRONG assignment found:", assignment)
                        comment = " <WRONG>"
                elif spectrum_combo in ['TOCSY (average)', 'TOCSY (unmatched)']:
                    user_assignment = get_Carbon_resonance(resid, nucleus, 'TOCSY',
                                                           self.user_HCNH_resid_assignments_dict,
                                                           self.user_TOCSY_resid_assignments_dict,
                                                           self.aa_equivalentCarbons_dict,
                                                           self.ali.protein_alignment_list,
                                                           self.ali.absolute_matches_alignment_list, self.args)
                    if user_assignment == False:  # if no user assignment was made for this Carbon, add no extra comment
                        print("WARNING: no user assignment made in TOCSY for resid", resid, "and Carbon", nucleus,
                              ". Checking HCNH...")
                        user_assignment = get_Carbon_resonance(resid,
                                                               nucleus,
                                                               spectrum_combo='HCNH',
                                                               user_HCNH_resid_assignments_dict=self.user_HCNH_resid_assignments_dict,
                                                               user_TOCSY_resid_assignments_dict=self.user_TOCSY_resid_assignments_dict,
                                                               aa_equivalentCarbons_dict=self.aa_equivalentCarbons_dict,
                                                               protein_alignment_list=self.ali.protein_alignment_list,
                                                               absolute_matches_alignment_list=self.ali.absolute_matches_alignment_list,
                                                               args=self.args)
                        if user_assignment == False:
                            print("WARNING: no user assignment made also in HCNH for resid", resid, "and Carbon",
                                  nucleus, ". Skipping nucleus.")
                            assignment[7] += " <NOT ASSIGNED>"
                            # self.proofread_clean_resid_assignments_dict[resid].append(assignment)
                            continue
                    if approx_equal(assignment[1], user_assignment,
                                    0.4):  # if the user-made assignment agrees with ours, add to the comment "CORRECT"
                        print("DEBUG: CORRECT assignment found:", assignment)
                        comment = " <CORRECT>"
                    else:
                        print("DEBUG: WRONG assignment found:", assignment)
                        comment = " <WRONG>"

            if comment == None:  # this is a proton or a methyl Carbon or a Carbon with only one proton
                # assignment[7] += iteration_comment
                continue
            elif comment == "":
                comment = " <NOT ASSIGNED>"
            assignment[7] += comment

    print("DEBUG: self.proofread_clean_resid_assignments_dict=", self.proofread_clean_resid_assignments_dict)

    ## WRITE THE PROOF-READ XEASY FILE
    xeasy_out = open(self.args.NHmap_FILE + ".HCNH.proofread.iter%s.xeasy" % iteration, 'w')
    atom_index = 0  # reinitialize atom index
    for resid in list(self.proofread_clean_resid_assignments_dict.keys()):
        for peak in self.proofread_clean_resid_assignments_dict[resid]:
            atom_index += 1
            peak[0] = atom_index
            line = "\t".join([str(p) for p in peak])
            xeasy_out.write(line + "\n")
    xeasy_out.close()