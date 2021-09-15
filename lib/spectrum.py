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


import csv

from .alignment import *
from .aaig import *
import pandas as pd

class Spectrum(object):
    # Class to store all information about a NOESY of TOCSY spectrum (AAIGs, connectivities, etc.)

    def __init__(self, AAIGs_dict={}, spectrum_type=""):
        self.AAIGs_dict = AAIGs_dict
        self.spectrum_type = spectrum_type

        # The following dictionaries are used to assign real residue codes to AAIG signatures
        self.residue2AAIGsign_dict = OrderedDict()  # OrderedDict with keys: residue (resname{3 letters}+resid) -> assigned AAIG signature
        # the AAIG c_orresponds to the i+1 residue, as that residue "sees" the TOCSY peaks of i residue
        self.AAIGsign2residue_dict = OrderedDict()  # AAIG signature -> residue (resname{3 letters}+resid)
        self.residue2position_dict = OrderedDict()  # OrderedDict with keys: residue (resname+resid) -> position in protein
        # sequence list (starting from 0 not starting_resid!)

    def __add_AAIG__(self, signature, AAIG):
        """
        It is recommended to include the N-H type into the AAIG signature, e.g. signature="N207ND2-HD2', signature='N207N-H', or signature='N207' if
        the N-H type is uknown.
        :param signature:
        :param AAIG:
        :return:
        """
        self.AAIGs_dict[signature] = AAIG

    def __replace_AAIG__(self, AAIG):
        assert AAIG.signature in self.AAIGs_dict.keys(), \
            Debuginfo("ERROR: AAIG signature %s in not present in %s spectrum!"
                      " If you have changed it use the function '__replace_AAIG_by_signature__' instead by"
                      " supplying the old signature value." % (AAIG.signature, self.spectrum_type), fail=True)
        self.AAIGs_dict[AAIG.signature] = AAIG

    def __replace_AAIG_by_signature__(self, AAIG, AAIG_signature):
        """

        :param AAIG: new AAIG object with new signature
        :param AAIG_signature:  old signature
        :return:
        """
        assert AAIG_signature in self.AAIGs_dict.keys(), \
            Debuginfo("ERROR: AAIG signature %s in not present in %s spectrum!" % (AAIG_signature, self.spectrum_type),
                      fail=True)
        self.AAIGs_dict[AAIG_signature] = AAIG

    def __del_AAIG__(self, signature="", AAIG=None):

        if signature:    # delete AAIG by its signature
            del self.AAIGs_dict[signature]
        elif not signature and AAIG:   # delete by matching peaks
            for signature, aaig in list(self.AAIGs_dict.items()):
                if AAIG.does_AAIG_match(aaig):
                    del self.AAIGs_dict[signature]

    def __get_all_AAIG_signatures__(self):
        return list(self.AAIGs_dict.keys())

    def __get_AAIG_num__(self):
        return len(list(self.AAIGs_dict.keys()))

    def __exists__(self, AAIG_signature):
        return AAIG_signature in list(self.AAIGs_dict.keys())

    def __get_all_AAIGs__(self):
        return list(self.AAIGs_dict.values())

    def __get_all_nonoverlapping_AAIGs__(self):
        return [AAIG for AAIG in list(self.AAIGs_dict.values()) if len(AAIG.overlapping_AAIGs)==0]

    def __get_all_overlapping_AAIGs__(self):
        return [AAIG for AAIG in list(self.AAIGs_dict.values()) if len(AAIG.overlapping_AAIGs)>0]

    def __get_all_peaks__(self):
        return [peak for AAIG in list(self.AAIGs_dict.values()) for peak in AAIG.__get_all_peaks__()]

    def __get_AAIG__(self, AAIG_signature):
        """
        Get the AAIG object that matches the given signature.
        :param AAIG_signature:
        :return:
        """
        return self.AAIGs_dict[AAIG_signature]

    def __add_peak2AAIG__(self, peak, AAIG_signature):
        self.AAIGs_dict[AAIG_signature].__add_peak__(peak)

    def __is_equal__(self, spectrum):
        AAIGs1 = self.__get_all_AAIGs__()
        AAIGs2 = spectrum.__get_all_AAIGs__()
        matched_AAIGs1, matched_AAIGs2 = [], []
        for AAIG1 in AAIGs1:
            for AAIG2 in AAIGs2:
                if AAIG2 in matched_AAIGs2:
                    continue
                if AAIG1.__is_equal__(AAIG2):
                    matched_AAIGs1.append(AAIG1)
                    matched_AAIGs2.append(AAIG2)
                    break
        # If all the AAIGs between Spectrum1 and Spectrum2 match
        if len(matched_AAIGs1) == len(AAIGs1) and \
            len(matched_AAIGs1) == len(AAIGs2) and \
            len(AAIGs1) == len(AAIGs2):
            return True
        else:
            return False

    def __get_spectrum_info__(self):
        AAIG_num = len(self.__get_all_AAIGs__())
        peak_num = len([peak for aaig in self.__get_all_AAIGs__() for peak in aaig.__get_all_peaks__()])
        ave_peaknum_per_aaig = peak_num/float(AAIG_num)
        min_peaknum_per_aaig = np.min([len(aaig.__get_all_peaks__()) for aaig in self.__get_all_AAIGs__()])
        max_peaknum_per_aaig = np.max([len(aaig.__get_all_peaks__()) for aaig in self.__get_all_AAIGs__()])
        return AAIG_num, peak_num, ave_peaknum_per_aaig, min_peaknum_per_aaig, max_peaknum_per_aaig

    def __scale_intensities__(self, even_intensities=False, wrt_aaig=False):
        """

        :param even_intensities:
        :param wrt_aaig: scale peaks intensities separately for each AAIG by divide with max intensity of the current AAIG.
        :return:
        """
        intensities = np.array([p.Intensity for p in self.__get_all_peaks__() if p.Intensity])
        print(intensities.size)
        if intensities.size > 0 and intensities.max() != 1.0:  # if needed, scale intensities (range 0.0-1.0)
            max_intensity = float(intensities.max())    # for global scaling
            for aaig in self.__get_all_AAIGs__():
                if wrt_aaig:    # scale locally
                    max_intensity = np.max([p.Intensity for p in aaig.__get_all_peaks__() if p.Intensity])
                    # otherwise use the globally maximum intensity
                for peak in aaig.__get_all_peaks__():
                    if even_intensities:
                        peak.scaled_Intensity = 1.0
                    else:
                        peak.scaled_Intensity = peak.Intensity/max_intensity
                    aaig.__replace_peak__(peak)  # save the Peak with scaled_Intensity property
                self.__replace_AAIG__(aaig) # save the AAIG with Peaks that have scaled_Intensity property

        elif intensities.size == 0 and not even_intensities:
            ColorPrint("WARNING: spectrum " + self.spectrum_type + " has peaks without intensities!", "WARNING")
        else:   # for TOCSY
            for aaig in self.__get_all_AAIGs__():
                for peak in aaig.__get_all_peaks__():
                    peak.Intensity = 1.0
                    peak.scaled_Intensity = 1.0
                    aaig.__replace_peak__(peak)  # save the Peak with scaled_Intensity property
                self.__replace_AAIG__(aaig)  # save the AAIG with Peaks that have scaled_Intensity property

    def __get_min_scaled_Intensity__(self):
        """
        Usefull for peak augmentation before 2D-histogram creation.
        :return:
        """
        scaled_intensities = np.array([p.scaled_Intensity for p in self.__get_all_peaks__() if p.scaled_Intensity])
        return scaled_intensities.min()

    def __AAIG2hist__(self, AAIG_signature, even_intensities=False, scale_density=False, **hist_borders):
        """
        Calculates the 2D-histogram of the whole spin system of the specified AAIG.
        :param AAIG_signature:
        :param even_intensities:
        :param scale_density:
        :param hist_borders:
        :return:
        """
        self.AAIGs_dict[AAIG_signature].AAIG2hist(even_intensities=even_intensities, scale_density=scale_density,  **hist_borders)

    def __AAIG_Peaks2hist__(self, AAIG_signature, even_edges=True, even_intensities=False, scale_density=False):
        """
        Calculates the 2D-histograms of all the individual Peaks in the specified AAIG.
        """
        self.AAIGs_dict[AAIG_signature].Peaks2hist(even_edges=even_edges,
                                                   even_intensities=even_intensities,
                                                   scale_density=scale_density)

    def __delete_all_hist2D__(self):
        for signature, aaig in list(self.AAIGs_dict.items()):
            aaig.__delete_hist__()
            self.AAIGs_dict[signature] = aaig

    def uniquify_spectrum_lines(self, query_fname, original_query_contents=[], uniquify_lines=True):
        """
        Load the contents of HSQC, TOCSY or NOESY and remove replicate lines by comparing only resonances,
        not labels.
        :param query_fname: can be a CSV or a Sparky list file
        :return:
        query_lines:    list of unique lines
        patched_residues_list:  if TOCSY, this list contains the pathed AAIGs
        """
        if not original_query_contents:
            if query_fname.endswith(".csv"):
                original_query_contents = [l for l in csv.reader(open(query_fname))]
                if original_query_contents[0][0] == "peak_index":
                    original_query_contents = [o[1:] for o in original_query_contents]
            else:
                with open(query_fname, 'r') as f:
                    original_query_contents = [l.split() for l in f.readlines()]  # contents of original query_fname (4D TOCSY or 4D NOESY) in 5 column format (name H C N HN)

        query_lines = []
        query_CS_set = set()  # store the numbers here to discard replicate lines
        for word_list in original_query_contents:  # CLEANING THE SPECTRUM FROM IRRELEVANT AND REPLICATE LINES
            try:
                if len(word_list) >= 6:  # if there is intensity column check if it is a number (>= reads also HNNH with AAIG distance columns)
                    float(word_list[5])
                    word_list[5] = str(abs(float(word_list[5])))  # convert all intensities to positive numbers

                if len(word_list) < 3:
                    print(bcolors.WARNING + "WARNING: Discarding wrong line: \n" + "".join(word_list) + bcolors.ENDC)
                    continue
                CS_values = tuple([float(v) for v in word_list[1:5]])    # checking if the line contains numbers (resonances)
                if (not uniquify_lines) or (CS_values not in query_CS_set):
                    query_lines.append(" ".join(word_list) + "\n")
                    query_CS_set.add(CS_values)
                else:
                    print(bcolors.WARNING + "WARNING: Discarding replicate line: " + "".join(word_list) + bcolors.ENDC)
            except (IndexError, ValueError):
                print(bcolors.WARNING + "WARNING: Discarding invalid line: " + "".join(word_list) + bcolors.ENDC)
                print(word_list)

        # remove duplicate lines from query_fname (4D TOCSY or 4D NOESY)
        lines2remove_set = set()
        for qline1 in query_lines:
            try:
                counter = 0
                query1_words_list = qline1.split()
                for qline2 in query_lines:
                    query2_words_list = qline2.split()
                    # print "query1_words_list=", query1_words_list
                    # print "query2_words_list=", query2_words_list
                    matches = [approx_equal(float(q1_w), float(q2_w)) for q1_w, q2_w in
                     zip(query1_words_list[1:], query2_words_list[1:])]
                    if all(matches) and len(matches) > 0:
                        counter += 1
                        if counter > 1: # remove the replicate line without labels (if applicable)
                            if "?-?" in qline1 and "?-?" not in qline2:
                                lines2remove_set.add(qline1)
                            elif "?-?" not in qline1 and "?-?" in qline2:
                                lines2remove_set.add(qline2)
                            elif "?-?" in qline1 and "?-?" in qline2:
                                lines2remove_set.add(qline2)
            except (ValueError, IndexError):
                # print "WARNING: the 2nd and 3rd elements of the following root file line are not numbers:"
                # print "Root file line:", root_line
                continue

        # now remove the duplicate lines
        for qline in lines2remove_set:
            ColorPrint("WARNING: removing replicate line: %s\n" %qline, "WARNING")
            query_lines.remove(qline)

        # Read the patched N-terminal residues written as a comment assigned TOCSY file (in sparky format)
        patched_residues_list = []
        for line in reversed(original_query_contents):  # ATTENTION: read the original file contents, not query_contents!!!
            if line[:19] == '# PATCHED RESIDUES:':
                patched_residues_list = line[19:].split()
                break

        return query_lines, patched_residues_list   # patched_residues_list will be populated only in case of assigned TOCSY file

    def find_AAIG2residue_correspondaces(self, NHmap_file, starting_resid):
        """
        This method read the NHmap file and finds correspondances between residues and assigned AAIG signatures.
        E.g. ('LEU89', 'X2NXHX'), ('GLN90', 'X5NXHX'), etc.
        to.
        :param NHmap_file:
        :param starting_resid:
        :return: populates dictionaries self.residue2AAIGsign_dict and self.AAIGsign2residue_dict.
        """
        absolute_AAIGmatches_alignment_list, protein_alignment_list = \
            Alignment.read_NHmap_file(NHmap_file, get_protein_alignment=True)

        # Remove terminal 'N/A' otherwise they will cause you trouble later
        NA_indices = []
        if protein_alignment_list[0] == 'N/A' or absolute_AAIGmatches_alignment_list[0] == 'N/A':
            NA_indices.append(0)
        if protein_alignment_list[-1] == 'N/A' or absolute_AAIGmatches_alignment_list[-1] == 'N/A':
            NA_indices.append(len(protein_alignment_list) - 1)
        protein_alignment_list = [protein_alignment_list[i] for i in range(len(protein_alignment_list)) if
                                  not i in NA_indices]
        absolute_AAIGmatches_alignment_list = [absolute_AAIGmatches_alignment_list[i] for i in
                                               range(len(absolute_AAIGmatches_alignment_list)) if not i in NA_indices]
        print("DEBUG: protein_alignment_list=", protein_alignment_list)
        print("DEBUG: absolute_AAIGmatches_alignment_list=", absolute_AAIGmatches_alignment_list)

        for position in range(1, len(protein_alignment_list)):
            resid = position + starting_resid  # position+1 because we removed the first 'N/A' from protein_alignment_list
            resname = protein_alignment_list[position - 1]  # recall that the TOCSY resonances correspond to residue i-1
            AAIG_signature = absolute_AAIGmatches_alignment_list[position]
            if resname == 'N/A':
                continue
            residue = aa1to3_dict[resname] + str(resid)
            self.residue2position_dict[residue] = position
            if AAIG_signature == '-' or AAIG_signature == 'N/A':
                self.residue2AAIGsign_dict[residue] = None
            else:
                self.residue2AAIGsign_dict[residue] = AAIG_signature
                self.AAIGsign2residue_dict[AAIG_signature] = residue

    def rename_AAIGs_to_residues(self, NHmap_file, starting_resid):
        """
        This method read the NHmap file and renames the AAIG signatures using the real resnames and resids they correspond
        to.
        :param NHmap_file:
        :param starting_resid:
        :return: populates dictionaries self.residue2AAIGsign_dict and self.AAIGsign2residue_dict.
        """

        ## READ FILE WITH ABSOLUTE MATCHES
        # absolute_AAIGmatches_alignment_list contains the absolutely matched HSQC AAIGs, namely the names given to each
        # N,HN peak in the HSQC. It does not necessarily contain real residue names!

        ColorPrint("Renaming AAIG signatures in %s spectrum to real residue codenames according to the NHmap file.",
                   "BOLDGREEN")
        self.find_AAIG2residue_correspondaces(NHmap_file, starting_resid)
        # Assign real residue codenames to AAIG signatures in this spectrum
        for AAIG in self.__get_all_AAIGs__():
            if AAIG.signature in list(self.AAIGsign2residue_dict.keys()):
                old_AAIG_signature = AAIG.signature
                AAIG.signature = aa3to1_dict[self.AAIGsign2residue_dict[AAIG.signature][:3]] + \
                                 self.AAIGsign2residue_dict[AAIG.signature][3:] + "NH" # update the signature
                # print("DEBUG: renaming in %s spectrum %s to %s" %
                #       (self.spectrum_type, old_AAIG_signature, AAIG.signature))
                AAIG.NH_name = "NH"    # if it was mapped to the sequence then it's definitely N-H
                # Update all Peak properties accordingly
                for peak in AAIG.__get_all_peaks__():
                    peak.j_AAIG_signature = AAIG.signature
                    peak.j_AAIG_name = remove_NH_suffix(AAIG.signature)
                    peak.j_Nname = 'N'
                    peak.j_HNname = 'H'
                    AAIG.__replace_peak__(peak)
                # Once you finish all renamings, replace the AAIG in the current spectrum to save the changes
                self.__replace_AAIG_by_signature__(AAIG, old_AAIG_signature)
            else:
                continue

    def write_sparky_list(self, fname):
        """
        Method to write the whole spectrum info into a Sparky list file. The AAIGs will be sorted
        :return:
        """
        # TODO: think if you must make it work with spectra saved in DataFrames, too.

        # TODO: untested for 'HNNH' and 'HSQC'.

        # So the AAIGs by their resid (virtual or real).
        assigned_AAIG_info_list, unassigned_AAIG_info_list = [], []
        for AAIG_signature in self.__get_all_AAIG_signatures__():
            AAIG_name, NH_name = split_AAIG_signature(AAIG_signature)
            if re.match("^[ACDEFGHIKLMNPQRSTVWY][0-9]+$", AAIG_name):    # leave the "X[0-9]+" out
                assigned_AAIG_info_list.append( (AAIG_name[0], int(AAIG_name[1:]), NH_name) )
            else:
                unassigned_AAIG_info_list.append( (AAIG_signature, "", "") )
        assigned_AAIG_info_list.sort(key=itemgetter(1))

        # Now write the sorted AAIGs and their peaks into a file
        fout = open(fname, 'w')
        for info in assigned_AAIG_info_list + unassigned_AAIG_info_list:
            AAIG_signature = info[0] + str(info[1]) + info[2]
            aaig = self.__get_AAIG__(AAIG_signature)
            for peak in aaig.__get_all_peaks__():
                if self.spectrum_type in ["TOCSY", "HCNH"]:
                    if peak.i_AAIG_name in [None, '?']:
                        i_label = "?-?"
                    else:
                        i_label = "%s%s-%s" % (peak.i_AAIG_name, peak.Hname, peak.Cname)
                    fout.write("%s-%s\t\t%.3f\t%.3f\t%.3f\t%.3f\t%i\n" %
                               (i_label, dash_to_NH(peak.j_AAIG_signature),
                                peak.Hreson, peak.Creson, peak.j_Nreson, peak.j_HNreson, peak.Intensity))
                elif self.spectrum_type in ['HNNH']:
                    fout.write("%s%s-%s-%s\t\t%.3f\t%.3f\t%.3f\t%.3f\t%i\n" %
                               (peak.i_AAIG_name, peak.i_HNname, peak.i_Nname, dash_to_NH(peak.j_AAIG_signature),
                                peak.Hreson, peak.Creson, peak.j_Nreson, peak.j_HNreson, peak.Intensity))
                elif self.spectrum_type in ["HSQC"]:
                    fout.write("%s%s-%s\t\t%.3f\t%.3f\n" %
                               (peak.j_AAIG_name, peak.j_HNname, peak.j_Nname,
                                peak.j_Nreson, peak.j_HNreson))
        fout.close()