from lib.csa_for_noesy import *
from lib.global_func import get_resid_from_residue


def onlynoesy_csa_iter2(self, iter1_pkl=None):

    if iter1_pkl:
        # LOAD THE NECESSARY DATA FROM ITERATION 1
        self.matched_HCNH_residue_assignments_dict, \
        self.checked_residues_set, \
        self.TOCSY_assigned_residue_list, \
        self.residues_with_written_NH_list, \
        self.clean_resid_assignments_dict, \
        self.revordered_residue_keys_iter1 = load_pickle(iter1_pkl)

    print("DEBUG: ENTERING ITERATION 2 self.matched_HCNH_residue_assignments_dict=",
          self.matched_HCNH_residue_assignments_dict)
    print("DEBUG: ENTERING ITERATION 2 self.clean_resid_assignments_dict=", self.clean_resid_assignments_dict)

    ## FIND MISSING CARBON ASSIGNMENTS AND TRY TO PREDICT THEM FROM HCNH
    # ATTENTION: recall that in TOCSY_residue_assignments_dict both the keys and the C,H assignments correspond to residue i-1.
    # But in matched_noesy_residue_assignments_dict the keys are residue i and the values are peaks from residues i (mostly the strongest), residue i-1 and maybe i+1 or
    # other neighboring residues.
    # So in order to match the HCNH peaks of residue i we need primarily the values of TOCSY_residue_assignments_dict keys i respectively.
    xeasy_fout = open(self.args.NHmap_FILE + ".HCNH.iter2.xeasy", 'w')
    sparky_fout = open(self.args.NHmap_FILE + ".HCNH.iter2.sparky", 'w')
    atom_index = 1
    # start the peak assignment from the last residue to the first
    self.revordered_residue_keys_iter2 = [(k[0], int(k[1:])) for k in
                                     list(self.matched_HCNH_residue_assignments_dict.keys()) if
                                     k in self.ali.absolute_matches_alignment_list]
    self.revordered_residue_keys_iter2.sort(key=itemgetter(1),
                                       reverse=True)  # reverse ordered keys, namely they start from the last residue in the protein sequence
    self.matched_i_to_iplus1_peaks_mdict = tree()  # residue -> i peak -> matched i+1 peak
    residues_with_written_prevAssigned_peaks = []

    for i_aa_resid in self.revordered_residue_keys_iter2:  # both the keys and the C,H assignments correspond to residue i-1
        i_residue = i_aa_resid[0] + str(i_aa_resid[1])
        # if not i_residue in ['L215']:
        #    continue
        # FIRST MATCH THE TOCSY ASSIGNMENTS OF RESIDUE i TO HCNH PEAKS
        # i_residue = get_i_residue_from_iminus1_residue(iminus1_residue)
        all_possible_assignments_list = []  # list with all the peaks for this residue (TOCSY-HCNH, HCNH) and the respective possible C-H assignments
        print("ITERATION 2: Matching residue i ", i_residue, " TOCSY PEAKS to residue's i HCNH peaks...")
        i_resid = get_resid_from_residue(i_residue)
        self.checked_residues_set.add(i_residue)
        aa_type = aa1to3_dict[get_aa_type_from_residue(i_residue)]
        all_allowed_carbons_list = [x for x in allowed_aa_atoms_dict[aa_type] if x[0] == "C" and not (
                    aa_type in list(self.iter1_allowed_atoms_dict.keys()) and x in self.iter1_allowed_atoms_dict[
                aa_type])]  # all allowed carbons for this aa type (exclude methyls)
        print("DEBUG: all_allowed_carbons_list=", all_allowed_carbons_list)
        matched_iplus1, \
        self.matched_HCNH_residue_assignments_dict, \
        self.matched_i_to_iplus1_peaks_mdict = \
            match_HCNH_i_iplus1_peaks(i_residue,
                                       self.ali.absolute_matches_alignment_list,
                                       self.matched_HCNH_residue_assignments_dict,
                                       self.matched_i_to_iplus1_peaks_mdict)  # add labels to the i->i+1 matched peaks. If i_residue is a C-term residue or presents a gap, skip it from this iteration
        if matched_iplus1 == False:
            continue
        else:  # otherwise proceed to assigned but not from the extreme carbons in exclude_from_i_iplus_matching_dict
            assigned_carbons_list = []  # list of carbons which have ALL their protons assigned in TOCSY
            partly_assigned_carbons_list = []  # list of methylene carbons for which only one proton frequency was found in TOCSY
            missing_carbons_list = []
            tmp_missing_carbons_list, tmp_partly_assigned_carbons_list = get_missing_carbons_from_xeasy_dict(aa_type,
                                                                                                             i_resid,
                                                                                                             self.clean_resid_assignments_dict,
                                                                                                             self.matched_HCNH_residue_assignments_dict,
                                                                                                             allowed_aa_atoms_dict)
            print("DEBUG: tmp_missing_carbons_list=", tmp_missing_carbons_list, "tmp_partly_assigned_carbons_list=",
                  tmp_partly_assigned_carbons_list)
            for carbon in tmp_missing_carbons_list:
                if aa_type in list(self.iter1_allowed_atoms_dict.keys()) and carbon in self.iter1_allowed_atoms_dict[
                    aa_type]:  # to exclude methyls from this iteration
                    continue
                if not aa_type in list(exclude_from_i_iplus_matching_dict.keys()) or (
                not carbon in exclude_from_i_iplus_matching_dict[aa_type]):
                    missing_carbons_list.append(carbon)
            for carbon in tmp_partly_assigned_carbons_list:
                if not aa_type in list(exclude_from_i_iplus_matching_dict.keys()) or (
                not carbon in exclude_from_i_iplus_matching_dict[aa_type]):
                    partly_assigned_carbons_list.append(carbon)
        print("DEBUG: dict key i_residue", i_residue, "missing_carbons_list=", missing_carbons_list)

        unassigned_HCNH_peaks_list = []  # ini
        if len(
                missing_carbons_list) != 0:  # IF THERE ARE MISSING CARBONS, TRY TO FIND THEM (EXCLUDE THE ASSIGNED OR MATCHED HCNH PEAKS)
            print("Residue ", i_residue,
                  " was not completely assigned in ITERATION 1. Using HCNH to assign resonances to the rest Carbon types.")
            try:  # DANGEROUS but works for the N-term residue (A200 in aLP)
                unassigned_HCNH_peaks_list = [peak for peak in self.matched_HCNH_residue_assignments_dict[i_residue]
                                               if peak[0] == "i+1" and peak[1] == "?" and peak[3] == "?"]
            except KeyError:
                continue
            print("DEBUG: unassigned_HCNH_peaks_list=", unassigned_HCNH_peaks_list)
            for HCNH_peak in unassigned_HCNH_peaks_list:  # make C-H type predictions for each unassigned HCNH peak
                HCNH_Creson = HCNH_peak[2]
                HCNH_Hreson = HCNH_peak[4]
                # OLD all_possible_assignments_list.extend(get_probabilities_from_H_C_resonpair_2Dhist(i_residue, HCNH_Hreson, HCNH_Creson, missing_carbons_list,
                # OLD                                                partly_assigned_carbons_list, self.args, self.aa2pthres_dict.aa_CHpair_binProbabilityList_mdict, self.aa2pthres_dict.aa_carbon_binDensityList_mdict,
                # OLD                                                self.HCNH_residue_peak_intensity_mdict, self.matched_i_to_iplus1_peaks_mdict, self.aa2pthres_dict,
                # OLD                                                self.aa2ithres_dict, self.ali.protein_alignment_list, self.ali.absolute_matches_alignment_list,
                # OLD                                                self.aa2pthres_dict.aa_hydrogen_binDensityList_mdict, iteration=2))
                all_possible_assignments_list.extend(get_probabilities_from_H_C_resonpair_2Dhist(i_residue,
                                                                                                 HCNH_Hreson,
                                                                                                 HCNH_Creson,
                                                                                                 missing_carbons_list,
                                                                                                 partly_assigned_carbons_list,
                                                                                                 self.args,
                                                                                                 self.HCNH_residue_peak_intensity_mdict,
                                                                                                 self.matched_i_to_iplus1_peaks_mdict,
                                                                                                 self.aa2pthres_dict,
                                                                                                 self.aa2ithres_dict,
                                                                                                 self.ali.protein_alignment_list,
                                                                                                 self.ali.absolute_matches_alignment_list,
                                                                                                 self.histload,
                                                                                                 iteration=1,
                                                                                                 I_THRESHOLD=self.I_THRESHOLD)
                                                     )
            [assignment.append('HCNH') for assignment in
             all_possible_assignments_list]  # we don't need to created a new list, the original all_possible_assignments_list will be updated
            print("DEBUG: i_residue=", i_residue, "all_possible_assignments_list=", all_possible_assignments_list)
            written_nucleiNames_list = []
            if i_resid in list(self.clean_resid_assignments_dict.keys()):
                written_nucleiNames_list = [line[3] for line in self.clean_resid_assignments_dict[
                    i_resid]]  # list of nuclei that have been written already, to avoid duplication
            if len(all_possible_assignments_list) > 0:  # IF WE HAVE FOUND POSSIBLE HCNH PEAK ASSIGNMENTS OF THE MISSING CARBONS
                # WRITE THE BEST HCNH PEAK COMBINATION TO COMPLETE THE ASSIGNMENT OF THE MISSING CARBONS
                extra_written_nucleiNames_list, \
                self.matched_HCNH_residue_assignments_dict, \
                self.TOCSY_assigned_residue_list = \
                    select_best_HCNH_peak_combination(i_residue,
                                                       all_possible_assignments_list,
                                                       partly_assigned_carbons_list,
                                                       atom_index,
                                                       xeasy_fout,
                                                       sparky_fout,
                                                       self.ali.protein_alignment_list,
                                                       self.ali.absolute_matches_alignment_list,
                                                       self.args,
                                                       self.HCNH_residue_NHresonances_dict,
                                                       self.matched_HCNH_residue_assignments_dict,
                                                       self.TOCSY_assigned_residue_list,
                                                       self.iter1_allowed_atoms_dict,
                                                       iteration=2)
                written_nucleiNames_list.extend(extra_written_nucleiNames_list)
            # WRITE ALSO THE TOCSY-HCNH MATCHED PEAKS OF THIS RESIDUE
            residues_with_written_prevAssigned_peaks.append(i_residue)
            # WRITE ALSO THE PEAKS ASSIGNED IN ITERATION 1
            write_peaks_assigned_in_previous_iteration(i_residue, self.clean_resid_assignments_dict, atom_index, xeasy_fout,
                                                       self.HCNH_residue_NHresonances_dict,
                                                       self.residues_with_written_NH_list)

        elif len(missing_carbons_list) == 0 and len(
                partly_assigned_carbons_list) == 0:  # NOTHING TO DO FOR THIS RESIDUE
            residues_with_written_prevAssigned_peaks.append(i_residue)
            # WRITE ALSO THE PEAKS ASSIGNED IN ITERATION 1 & 2
            write_peaks_assigned_in_previous_iteration(i_residue, self.clean_resid_assignments_dict, atom_index, xeasy_fout,
                                                       self.HCNH_residue_NHresonances_dict,
                                                       self.residues_with_written_NH_list)

        else:  # IF ALL THE CARBONS WERE ASSIGNED, FIND & SAVE PARTLY ASSIGNED CARBONS, SAVE THE MATCHED HCNH RESONANCES, AND IF NOT MATCHED, SAVE THE TOCSY RESONANCES
            print("Residue ", i_residue, " was completely assigned in ITERATION 1.")
            try:
                aa_type = aa1to3_dict[get_aa_type_from_residue(i_residue)]
            except (IndexError, KeyError):  # skip residue i
                print("EXITING: COULD NOT FIND AA_TYPE OF RESIDUE", i_residue)
                sys.exit(1)

            assigned_peaks_list = get_peaks_from_xeasy(i_resid, self.clean_resid_assignments_dict, atom_index, xeasy_fout,
                                                       self.HCNH_residue_NHresonances_dict,
                                                       self.residues_with_written_NH_list)  # get all assigned peaks of this residue from xeasy format
            for assigned_peak in assigned_peaks_list:
                spectrum_combo = assigned_peak[5]
                Cname = assigned_peak[1]
                Creson = assigned_peak[2]  # use the TOCSY resonance by default
                Hname = assigned_peak[3]
                Hreson = assigned_peak[4]
                all_possible_assignments_list.append([aa_type, Cname, Hname, 10000000.0, Hreson, Creson,
                                                      spectrum_combo])  # set the probability to 1 cause it has been already assigned in iteration 1
                print("DEBUG: i_residue=", i_residue, "aa_type=", aa_type, "all_possible_assignments_list=",
                      all_possible_assignments_list)
            written_nucleiNames_list = []
            if i_resid in list(self.clean_resid_assignments_dict.keys()):
                written_nucleiNames_list = [line[3] for line in self.clean_resid_assignments_dict[
                    i_resid]]  # list of nuclei that have been written already, to avoid duplication
            if len(all_possible_assignments_list) > 0:  # IF WE HAVE FOUND POSSIBLE HCNH PEAK ASSIGNMENTS OF THE MISSING CARBONS
                # WRITE THE BEST HCNH PEAK COMBINATION TO COMPLETE THE ASSIGNMENT OF THE MISSING CARBONS
                written_nucleiNames_list,\
                self.matched_HCNH_residue_assignments_dict, \
                self.TOCSY_assigned_residue_list = \
                    select_best_HCNH_peak_combination(i_residue,
                                                      all_possible_assignments_list,
                                                      partly_assigned_carbons_list,
                                                      atom_index,
                                                      xeasy_fout,
                                                      sparky_fout,
                                                      self.ali.protein_alignment_list,
                                                      self.ali.absolute_matches_alignment_list,
                                                      self.args,
                                                      self.HCNH_residue_NHresonances_dict,
                                                      self.matched_HCNH_residue_assignments_dict,
                                                      self.TOCSY_assigned_residue_list,
                                                      self.iter1_allowed_atoms_dict,
                                                      iteration=2)
        self.matched_HCNH_residue_assignments_dict = \
            remove_iplus1_labels(i_residue,
                                 self.matched_HCNH_residue_assignments_dict,
                                 self.ali.protein_alignment_list,
                                 self.ali.absolute_matches_alignment_list,
                                 self.args)  # remove all the 'i+1', 'i-1' labels from this residue and residue i+1 (if applicable)

    print("DEBUG iteration2 point2: self.matched_HCNH_residue_assignments_dict=",
          self.matched_HCNH_residue_assignments_dict)

    # WRITE ALSO THE PEAKS ASSIGNED IN PREVIOUS ITERATIONS OF ALL RESIDUES THAT WERE EXCLUDED FROM THIS ITERATION
    for i_aa_resid in set(self.revordered_residue_keys_iter1 + self.revordered_residue_keys_iter2):
        i_residue = i_aa_resid[0] + str(i_aa_resid[1])
        aa_type = aa1to3_dict[i_aa_resid[0]]
        if not i_residue in residues_with_written_prevAssigned_peaks:
            write_peaks_assigned_in_previous_iteration(i_residue, self.clean_resid_assignments_dict, atom_index, xeasy_fout,
                                                       self.HCNH_residue_NHresonances_dict,
                                                       self.residues_with_written_NH_list)

    # Now find if there are any residue in the alignment that haven't been checked:
    for residue in self.ali.absolute_matches_alignment_list:
        if not residue in self.checked_residues_set and residue != '-':
            print("WARNING: Residue ", residue, " was not checked!!!")

    # ##
    # ##  WRITE THE HN & N RESONANCES TO XEASY ONLY FOR THOSE RESIDUES THAT NO C,H ASSIGNMENTS COULD BE MADE BUT THEIR HN & H RESONANCES ARE KNOWN FROM TOCSY OR HCNH
    # ##
    write_NH_of_residues_without_CH(self.ali.absolute_matches_alignment_list,
                                    self.HCNH_residue_NHresonances_dict,
                                    {},  # TODO: it was TOCSY_residue_NHresonances_dict before
                                    self.residues_with_written_NH_list,
                                    atom_index,
                                    xeasy_fout,
                                    sparky_fout,
                                    self.ali.protein_alignment_list,
                                    self.args,
                                    self.Cterm_residue_set,
                                    iteration=2)

    xeasy_fout.close()
    sparky_fout.close()

    ##
    ##  REWRITE THE XEASY FILE BY REMOVING DUPLICATES OR AVERAGING, AND CORRECTING THE ORDER
    ##

    ## CLEAN THE CONTENTS OF XEASY FILE
    # self.clean_resid_assignments_dict contains the final, non-redundant assignments written to xeasy file. E.g.
    # 1: [[2336, 56.022, 0.2, 'CA', 1, '#', 'MET', 'TOCSY (unmatched) (iter2)'], [2338, 34.786, 0.2, 'CB', 1, '#', 'MET', 'TOCSY (unmatched) (iter2)'],
    # [2340, 30.946, 0.2, 'CG', 1, '#', 'MET', 'TOCSY (unmatched) (iter2)'], [2337, 4.194, 0.02, 'HA', 1, '#', 'MET', 'TOCSY (unmatched) (iter2)'],
    # [2339, 2.14, 0.02, 'QB', 1, '#', 'MET', 'TOCSY (unmatched) (iter2)'], [2341, 2.487, 0.02, 'QG', 1, '#', 'MET', 'TOCSY (unmatched) (iter2)']]
    self.clean_resid_assignments_dict = clean_xeasy_file(self.args.NHmap_FILE + ".HCNH.iter2.xeasy")
    # sort the peak assignments of each resid
    for resid, peak_list in list(self.clean_resid_assignments_dict.items()):
        peak_list.sort(key=itemgetter(3))
        self.clean_resid_assignments_dict[resid] = peak_list

    print("DEBUG: ITERATION 2 point 1 self.clean_resid_assignments_dict=", self.clean_resid_assignments_dict)
    ## WRITE THE CLEAN XEASY FILE
    xeasy_out = open(self.args.NHmap_FILE + ".HCNH.cleaned.iter2.xeasy", 'w')
    atom_index = 0  # reinitialize atom index
    for resid in list(self.clean_resid_assignments_dict.keys()):
        for peak in self.clean_resid_assignments_dict[resid]:
            atom_index += 1
            peak[0] = atom_index
            line = "\t".join([str(p) for p in peak])
            xeasy_out.write(line + "\n")
    xeasy_out.close()

    ##
    ##  IF PROVIDED, LOAD THE USER-ANNOTATED TOCSY AND HCNH FILES FOR PROOF-READING
    ##
    print("DEBUG: self.args.user_TOCSY_FILE=", self.args.user_TOCSY_FILE, "self.args.user_HCNH_FILE=", self.args.user_HCNH_FILE)
    self.user_HCNH_resid_assignments_dict = {}  # same as user_HCNH_residue_assignments_dict, but with keys resids instead of residues
    self.user_TOCSY_resid_assignments_dict = {}  # same as self.user_TOCSY_resid_assignments_dict, but with keys resids instead of residues
    self.user_VAL_LEU_HCNH_resid_assignments_dict = {}  # ONLY ASSIGNED VAL AND LEU EQUIVALENT CARBONS DICTIONARIES
    self.user_VAL_LEU_TOCSY_resid_assignments_dict = {}
    if self.args.user_TOCSY_FILE and self.args.user_HCNH_FILE:
        self.proofread_iter12(iteration="2")

    # Save the necessary objects for the next ITERATION
    save_pickle("iter2.pkl",
                self.matched_HCNH_residue_assignments_dict,
                self.checked_residues_set,
                self.TOCSY_assigned_residue_list,
                self.residues_with_written_NH_list,
                self.clean_resid_assignments_dict,
                self.revordered_residue_keys_iter1,
                self.revordered_residue_keys_iter2,
                self.matched_i_to_iplus1_peaks_mdict,
                self.proofread_clean_resid_assignments_dict,
                self.user_HCNH_resid_assignments_dict,
                self.user_TOCSY_resid_assignments_dict,
                self.user_VAL_LEU_HCNH_resid_assignments_dict,
                self.user_VAL_LEU_TOCSY_resid_assignments_dict
                )

#                   REDUNDANT TODO: more thorough testing needed.
# def proofread_iter2_short(self):
#     print("STARTING PROOF-READING HCNH CHEMICAL SHIFT ASSIGNMENTS")
#
#     # LOAD USER ASSIGNED TOCSY CONTENTS
#     user_TOCSY_residue_assignments_dict, \
#     user_TOCSY_residue_NHresonances_dict, \
#     user_i_to_iminus1_dict, \
#     user_patched_residues_list, \
#     user_TOCSY_residue_peak_intensity_mdict = \
#         read_NHassigned_spectrum_file(self.args.user_TOCSY_FILE,
#                                       "TOCSY",
#                                       self.ali.absolute_AAIGmatches_alignment_list,
#                                       self.ali.absolute_matches_alignment_list
#                                       )
#
#     # LOAD USER ASSIGNED HCNH FILE CONTENTS
#     user_HCNH_residue_assignments_dict, \
#     user_HCNH_residue_NHresonances_dict, \
#     user_original_HCNH_peaks_dict, \
#     user_HCNH_residue_peak_instensity_dict = \
#         read_NHassigned_spectrum_file(
#             self.args.user_HCNH_FILE,
#             'HCNH',
#             self.ali.absolute_AAIGmatches_alignment_list,
#             self.ali.absolute_matches_alignment_list
#         )
#
#     # rename TYR and PHE CD-QD>CD1-HD1 and CE-QE->CE1-HE1
#     for residue in user_HCNH_residue_assignments_dict:
#         if residue[0] in ['Y', 'F']:
#             for assignment in user_HCNH_residue_assignments_dict[residue]:
#                 if assignment[1] == "CD":
#                     assignment[1] = "CD1"
#                 elif assignment[1] == "QD":
#                     assignment[1] = "HD1"
#                 elif assignment[1] == "CE":
#                     assignment[1] = "CE1"
#                 elif assignment[1] == "QE":
#                     assignment[1] = "HE1"
#     for residue in user_TOCSY_residue_assignments_dict:
#         if residue[0] in ['Y', 'F']:
#             for assignment in user_TOCSY_residue_assignments_dict[residue]:
#                 if assignment[1] == "CD":
#                     assignment[1] = "CD1"
#                 elif assignment[1] == "QD":
#                     assignment[1] = "HD1"
#                 elif assignment[1] == "CE":
#                     assignment[1] = "CE1"
#                 elif assignment[1] == "QE":
#                     assignment[1] = "HE1"
#
#     # convert the keys from residues to resids to be compatible with self.clean_resid_assignments_dict
#     for k, v in list(user_HCNH_residue_assignments_dict.items()):
#         if not is_valid_residue(k, self.ali.protein_alignment_list, self.args):
#             print("WARNING: omitting non-valid user assigned residue ", k, "!")
#             continue
#         self.user_HCNH_resid_assignments_dict[int(k[1:])] = v
#     for k, v in list(user_TOCSY_residue_assignments_dict.items()):
#         if not is_valid_residue(k, self.ali.protein_alignment_list, self.args):
#             print("WARNING: omitting non-valid user assigned residue ", k, "!")
#             continue
#         self.user_TOCSY_resid_assignments_dict[int(k[1:])] = v
#
#     # CREATE ONLY ASSIGNED VAL AND LEU EQUIVALENT CARBONS DICTIONARIES
#     self.user_VAL_LEU_HCNH_resid_assignments_dict, \
#     self.user_VAL_LEU_TOCSY_resid_assignments_dict = \
#         create_equivalentCarbons_assignments_dict(self.user_HCNH_resid_assignments_dict,
#                                                   self.user_TOCSY_resid_assignments_dict,
#                                                   self.ali.protein_alignment_list,
#                                                   self.ali.absolute_matches_alignment_list,
#                                                   self.aa_equivalentCarbons_dict,
#                                                   self.args)
#
#     self.proofread_clean_resid_assignments_dict = copy.deepcopy(
#         self.clean_resid_assignments_dict)  # same as self.clean_resid_assignments_dict but with comments CORRECT or WRONG at the carbon lines
#     for resid in list(self.proofread_clean_resid_assignments_dict.keys()):
#         resname = self.proofread_clean_resid_assignments_dict[resid][0][6]
#         if resname in list(self.aa_equivalentCarbons_dict.keys()):  # If this is a LEU or VAL
#             Carbon_pair = self.aa_equivalentCarbons_dict[resname][0]
#             if len([a for a in self.proofread_clean_resid_assignments_dict[resid] if
#                     a[3] in Carbon_pair]) > 0:  # AND if the methyl carbons have been assigned
#                 Proton_pair = ["H" + C[1:] for C in Carbon_pair]  # protons of the equivalent Carbons
#                 # get only the methyl Carbon assignments (LEU --> CD1,CD2 or VAL --> CG1,CG2) made by the program, to compare them with the user made assignments
#                 program_assignments_list = [assignment for assignment in
#                                             self.proofread_clean_resid_assignments_dict[resid] if assignment[
#                                                 3] in Carbon_pair + Proton_pair]  # keep only the resonances assigned by the program to the equivalent methyl carbons and protons
#                 # proof_read_equivalentCarbons(program_assignments_list, resid, self.user_VAL_LEU_HCNH_resid_assignments_dict, self.user_VAL_LEU_TOCSY_resid_assignments_dict, Carbon_pair)   # OBSOLETE!!! it doesn't matter if 1 or 2 of the peaks have been proof-read in iteration 1 & 2
#                 proof_read_equivalentMethyls(program_assignments_list, resid,
#                                              self.user_VAL_LEU_HCNH_resid_assignments_dict,
#                                              self.user_VAL_LEU_TOCSY_resid_assignments_dict,
#                                              Carbon_pair,
#                                              Proton_pair)  # it doesn't matter if 1 or 2 of the peaks have been proof-read in iteration 1 & 2
#         # THEN PROCEED NORMALLY
#         for assignment in self.proofread_clean_resid_assignments_dict[resid]:
#             spectrum_combo = re.sub(r' iter[1-9]$', '', assignment[
#                 7])  # remove the "iter[1-9]" suffix to be able to distinguish the type of assignment (HCNH, HCNH+TOCSY, etc.)
#             comment = None
#             aa_type = assignment[6]
#             nucleus = assignment[3]  # the nucleus name
#             if aa_type in list(self.aa_equivalentCarbons_dict.keys()) and nucleus in \
#                     self.aa_equivalentCarbons_dict[aa_type][
#                         0]:  # if this peak has been assigned by the program to LEU or VAL CD1, CD2 or CG1,CG2 respectively, skip it
#                 continue
#
#             # IF THIS CARBON IS A METHYL OR HAS ONLY ONE PROTON, CONSIDE ALSO THE PROTON RESONANCE FOR PROOFREADING
#             if aa_type in list(aatype_carbon_nongeminalHname_mdict.keys()) and nucleus in list(
#                     aatype_carbon_nongeminalHname_mdict[aa_type].keys()):
#                 # get the proton assignment made by the program
#                 C_assignment = assignment  # the program-made carbon assignment
#                 H_assignment = None  # default value for the program-made proton assignment
#                 for a in self.proofread_clean_resid_assignments_dict[resid]:
#                     if ((type(aatype_carbon_nongeminalHname_mdict[aa_type][nucleus]) == list and a[3] in
#                          aatype_carbon_nongeminalHname_mdict[aa_type][nucleus]) or
#                             (type(aatype_carbon_nongeminalHname_mdict[aa_type][nucleus]) == str and a[3] ==
#                              aatype_carbon_nongeminalHname_mdict[aa_type][nucleus])):
#                         H_assignment = a
#                         break
#                 comment = ""
#                 if spectrum_combo in ['HCNH', 'HCNH (C-term hanging residue)', 'HCNH + TOCSY-HCNH', 'TOCSY-HCNH',
#                                      'TOCSY-HCNH + HCNH', 'TOCSY (unmatched) + TOCSY-HCNH',
#                                      'TOCSY-HCNH + TOCSY (unmatched)']:  # if this resonance comes only from HCNH
#                     user_assignment = get_Carbon_resonance(resid, nucleus, 'HCNH',
#                                                            self.user_HCNH_resid_assignments_dict,
#                                                            self.user_TOCSY_resid_assignments_dict,
#                                                            self.aa_equivalentCarbons_dict,
#                                                            self.ali.protein_alignment_list,
#                                                            self.ali.absolute_matches_alignment_list, self.args,
#                                                            get_assignment=True)
#                     if user_assignment == False:  # if no user assignment was made for this Carbon, add no extra comment
#                         print("WARNING: no user assignment made in HCNH for resid", resid, "and Carbon", nucleus,
#                               ". Checking TOCSY...")
#                         user_assignment = get_Carbon_resonance(resid, nucleus, 'TOCSY',
#                                                                self.user_HCNH_resid_assignments_dict,
#                                                                self.user_TOCSY_resid_assignments_dict,
#                                                                self.aa_equivalentCarbons_dict,
#                                                                self.ali.protein_alignment_list,
#                                                                self.ali.absolute_matches_alignment_list, self.args,
#                                                                get_assignment=True)
#                         if user_assignment == False:
#                             print("WARNING: no user assignment made also in TOCSY for resid", resid, "and Carbon",
#                                   nucleus, ". Skipping nucleus.")
#                             assignment[7] += " <NOT ASSIGNED>"
#                             # self.proofread_clean_resid_assignments_dict[resid].append(assignment)
#                             continue
#                     user_Creson = user_assignment[2]
#                     user_Hreson = user_assignment[4]
#                     if approx_equal(C_assignment[1], user_Creson, 0.4) and approx_equal(H_assignment[1], user_Hreson,
#                                                                                         0.04):  # if the user-made assignment agrees with ours, add to the comment "CORRECT"
#                         comment = " <CORRECT>"
#                     else:
#                         comment = " <WRONG>"
#                 elif spectrum_combo in ['TOCSY (average)', 'TOCSY (unmatched)']:
#                     user_assignment = get_Carbon_resonance(resid, nucleus, 'TOCSY',
#                                                            self.user_HCNH_resid_assignments_dict,
#                                                            self.user_TOCSY_resid_assignments_dict,
#                                                            self.aa_equivalentCarbons_dict,
#                                                            self.ali.protein_alignment_list,
#                                                            self.ali.absolute_matches_alignment_list, self.args,
#                                                            get_assignment=True)
#                     if user_assignment == False:  # if no user assignment was made for this Carbon, add no extra comment
#                         print("WARNING: no user assignment made in TOCSY for resid", resid, "and Carbon", nucleus,
#                               ". Checking HCNH...")
#                         user_assignment = get_Carbon_resonance(resid, nucleus, 'HCNH',
#                                                                self.user_HCNH_resid_assignments_dict,
#                                                                self.user_TOCSY_resid_assignments_dict,
#                                                                self.aa_equivalentCarbons_dict,
#                                                                self.ali.protein_alignment_list,
#                                                                self.ali.absolute_matches_alignment_list, self.args,
#                                                                get_assignment=True)
#                         if user_assignment == False:
#                             print("WARNING: no user assignment made also in HCNH for resid", resid, "and Carbon",
#                                   nucleus, ". Skipping nucleus.")
#                             assignment[7] += " <NOT ASSIGNED>"
#                             # self.proofread_clean_resid_assignments_dict[resid].append(assignment)
#                             continue
#                     user_Creson = user_assignment[2]
#                     user_Hreson = user_assignment[4]
#                     if approx_equal(C_assignment[1], user_Creson, 0.4) and approx_equal(H_assignment[1], user_Hreson,
#                                                                                         0.04):  # if the user-made assignment agrees with ours, add to the comment "CORRECT"
#                         comment = " <CORRECT>"
#                     else:
#                         comment = " <WRONG>"
#
#                 if comment == "":
#                     comment = " <NOT ASSIGNED>"
#                 C_assignment[7] += comment
#                 H_assignment[7] += comment
#                 comment = None  # a trick to add no more comments to "assignment" at the end of this block of code
#
#             elif nucleus[
#                 0] == 'C':  # only if the assignment line belongs to a Carbon then check if the user has made the same assignment
#                 comment = ""
#                 if spectrum_combo in ['HCNH', 'HCNH (C-term hanging residue)', 'HCNH + TOCSY-HCNH', 'TOCSY-HCNH',
#                                      'TOCSY-HCNH + HCNH', 'TOCSY (unmatched) + TOCSY-HCNH',
#                                      'TOCSY-HCNH + TOCSY (unmatched)']:  # if this resonance comes only from HCNH
#                     user_assignment = get_Carbon_resonance(resid, nucleus, 'HCNH',
#                                                            self.user_HCNH_resid_assignments_dict,
#                                                            self.user_TOCSY_resid_assignments_dict,
#                                                            self.aa_equivalentCarbons_dict,
#                                                            self.ali.protein_alignment_list,
#                                                            self.ali.absolute_matches_alignment_list, self.args)
#                     if user_assignment == False:  # if no user assignment was made for this Carbon, add no extra comment
#                         print("WARNING: no user assignment made in HCNH for resid", resid, "and Carbon", nucleus,
#                               ". Checking TOCSY...")
#                         user_assignment = get_Carbon_resonance(resid, nucleus, 'TOCSY',
#                                                                self.user_HCNH_resid_assignments_dict,
#                                                                self.user_TOCSY_resid_assignments_dict,
#                                                                self.aa_equivalentCarbons_dict,
#                                                                self.ali.protein_alignment_list,
#                                                                self.ali.absolute_matches_alignment_list, self.args)
#                         if user_assignment == False:
#                             print("WARNING: no user assignment made also in TOCSY for resid", resid, "and Carbon",
#                                   nucleus, ". Skipping nucleus.")
#                             assignment[7] += " <NOT ASSIGNED>"
#                             continue
#                     if approx_equal(assignment[1], user_assignment,
#                                     0.4):  # if the user-made assignment agrees with ours, add to the comment "CORRECT"
#                         comment = " <CORRECT>"
#                     else:
#                         comment = " <WRONG>"
#                 elif spectrum_combo in ['TOCSY (average)', 'TOCSY (unmatched)']:
#                     user_assignment = get_Carbon_resonance(resid, nucleus, 'TOCSY',
#                                                            self.user_HCNH_resid_assignments_dict,
#                                                            self.user_TOCSY_resid_assignments_dict,
#                                                            self.aa_equivalentCarbons_dict,
#                                                            self.ali.protein_alignment_list,
#                                                            self.ali.absolute_matches_alignment_list, self.args)
#                     if user_assignment == False:  # if no user assignment was made for this Carbon, add no extra comment
#                         print("WARNING: no user assignment made in TOCSY for resid", resid, "and Carbon", nucleus,
#                               ". Checking HCNH...")
#                         user_assignment = get_Carbon_resonance(resid, nucleus, 'HCNH',
#                                                                self.user_HCNH_resid_assignments_dict,
#                                                                self.user_TOCSY_resid_assignments_dict,
#                                                                self.aa_equivalentCarbons_dict,
#                                                                self.ali.protein_alignment_list,
#                                                                self.ali.absolute_matches_alignment_list, self.args)
#                         if user_assignment == False:
#                             print("WARNING: no user assignment made also in HCNH for resid", resid, "and Carbon",
#                                   nucleus, ". Skipping nucleus.")
#                             assignment[7] += " <NOT ASSIGNED>"
#                             # self.proofread_clean_resid_assignments_dict[resid].append(assignment)
#                             continue
#                     if approx_equal(assignment[1], user_assignment,
#                                     0.4):  # if the user-made assignment agrees with ours, add to the comment "CORRECT"
#                         comment = " <CORRECT>"
#                     else:
#                         comment = " <WRONG>"
#
#             if comment == None:  # this is a proton or a methyl Carbon or a Carbon with only one proton
#                 continue
#             elif comment == "":
#                 comment = " <NOT ASSIGNED>"
#             assignment[7] += comment
#
#     ## WRITE THE PROOF-READ XEASY FILE
#     xeasy_out = open(self.args.NHmap_FILE + ".HCNH.proofread.iter2.xeasy", 'w')
#     atom_index = 0  # reinitialize atom index
#     for resid in list(self.proofread_clean_resid_assignments_dict.keys()):
#         for peak in self.proofread_clean_resid_assignments_dict[resid]:
#             atom_index += 1
#             peak[0] = atom_index
#             line = "\t".join([str(p) for p in peak])
#             xeasy_out.write(line + "\n")
#     xeasy_out.close()
