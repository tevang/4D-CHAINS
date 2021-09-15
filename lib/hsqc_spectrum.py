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


import re
from operator import itemgetter
import regex
from collections import defaultdict

from .csa_for_tocsy import get_protein_sequence
from .global_func import *
from .spectrum import *
from .peak import *

class HSQC_spectrum(Spectrum):

    def __init__(self,
                 HSQC_FILE,
                 fasta,
                 H_bin_length=0.01,
                 C_bin_length=0.1,
                 bandwidth=0.0025):

        super(HSQC_spectrum, self).__init__(spectrum_type = "HSQC")
        self.HSQC_FILE = HSQC_FILE
        self.fasta = fasta
        self.H_bin_length = H_bin_length
        self.C_bin_length = C_bin_length
        self.bandwidth = bandwidth
        self.numlist_files = {} # 4D spectrum type -> 4D spectrum file with labels copied from HSQC

        HSQC_spectrum.annotate_HSQC(HSQC_FILE, fasta, OUT_fname=HSQC_FILE + "_annotated")
        self.HSQC_FILE = HSQC_FILE + "_annotated"

        self.uniquify_HSQC_AAIGs()

    @staticmethod
    def annotate_HSQC(HSQC_FILE, FASTA_fname, STARTING_RESID=None, ABSOLUTE_MATCHES_FILE=None, OUT_fname=None):
        """
        This method reads in the {N-H}-HSQC file in sparky list format. If there are labels inside that comply
        with the pattern '[A-Z][0-9]+N-H' (e.g. A13N-H), the script will check if these labels are valid given
        the protein sequence. If yes, it will keep them and rename all other peaks as 'X[0-9]+'. If not, then
        it will ignore all existing labels and relabel all peaks.

        :param HSQC_FILE:
        :param STARTING_RESID:
        :param FASTA_fname:
        :param ABSOLUTE_MATCHES_FILE:
        :param OUT_fname:
        :return:
        """
        if ABSOLUTE_MATCHES_FILE:  # TRY TO GUESS STARTING RESID FROM AAIG LABELS IN HSQC
            if STARTING_RESID == None:
                print(bcolors.WARNING + "WARNING: you have not specified the starting residue number with the option -rstart. I will try \
        to guess it from the labels." + bcolors.ENDC)
                ## READ FILE WITH ABSOLUTE MATCHES
            absolute_AAIGmatches_alignment_list, protein_sequence_list = \
                Alignment.read_NHmap_file(ABSOLUTE_MATCHES_FILE, get_protein_alignment=True)

            TIG2residue_dict = {}
            for i, aa, TIG in zip(list(range(len(protein_sequence_list))), protein_sequence_list,
                                  absolute_AAIGmatches_alignment_list):
                if not TIG in ['N/A', '-']:
                    TIG2residue_dict[TIG] = aa + str(i + STARTING_RESID)
            # print "DEBUG: TIG2residue_dict=", TIG2residue_dict

            with open(HSQC_FILE, 'r') as f:
                contents = f.readlines()
                index = 1
                for i in range(len(contents)):
                    line = contents[i]
                    mo = re.match(".*\s+([0-9.]+)\s+([0-9.]+)\s*", line)  # the last 2 columns must be numbers
                    if mo:
                        Nreson = mo.group(1)
                        Hreson = mo.group(2)
                        mp = re.match("^\s*([A-Za-z-]+[0-9]+)+N-H\s+[0-9.]+\s+[0-9.]+\s*",
                                      line)  # consider only N-H mapped AAIGs, ignore NX-HX, etc.
                        if mp:
                            TIG = mp.group(1)
                            if TIG in list(TIG2residue_dict.keys()):
                                residue = TIG2residue_dict[TIG]
                            else:
                                residue = TIG
                            contents[i] = "\t%s\t%s\t%s\n" % (str(residue) + "N-H", Nreson, Hreson)
                            index += 1

        else:
            used_labels = []  # AAIG labels that have been assigned to HSQC peaks already
            if not FASTA_fname:
                print(
                    bcolors.FAIL + "You must also provide a fasta sequence file using the option -fasta !" + bcolors.ENDC)
                sys.exit(1)
            protein_sequence_list = get_protein_sequence(FASTA_fname)
            # print "DEBUG: protein_sequence_list=", protein_sequence_list
            STARTING_RESID = HSQC_spectrum.guess_starting_resid(HSQC_FILE=HSQC_FILE, fasta=FASTA_fname, NHmap=None)
            # If no NH-mapping table was given, consider labels meeting the conditions as real labels from the user and relabel the rest of the peaks
            with open(HSQC_FILE, 'r') as f:
                contents = f.readlines()
                index = 1
                for i in range(len(contents)):
                    line = contents[i]
                    mo = re.match(".*\s+([0-9.]+)\s+([0-9.]+)\s*", line)  # the last 2 columns must be numbers
                    if mo:
                        AAIG_signature = line.split()[0]
                        Nreson = mo.group(1)
                        Hreson = mo.group(2)
                        mp = re.match("^\s*([A-Z])([0-9]+)NX?-HX?\s+[0-9.]+\s+[0-9.]+\s*$", line)
                        mq = re.match("^\s*([0-9]+)NX?-HX?\s+[0-9.]+\s+[0-9.]+\s*$", line)
                        if mp:  # if the label matches the conditions, keep it
                            aa_type = mp.group(1)
                            resid = int(mp.group(2))
                            try:
                                if not STARTING_RESID: raise ValueError
                                # print "DEBUG: resid=", resid, "STARTING_RESID=", STARTING_RESID
                                if protein_sequence_list[resid - STARTING_RESID] == aa_type:
                                    continue
                            except ValueError:
                                contents[i] = "\t%s\t%s\t%s\n" % ("X" + str(index) + "NX-HX", Nreson, Hreson)
                                index += 1
                            except IndexError:  # resid too high or too small, probably a sidechain, leave it as it is
                                continue
                        # TODO: if it is side chain let it in without controling if it agrees with the protein sequence
                        elif is_sidechain(AAIG_signature):
                            pass
                        elif mq:  # if the label is of the form "[0-9]+NX?-HX?", then keep it to ease the user in the re-assignment
                            user_index = mq.group(1)
                            if "X" + user_index in used_labels:
                                raise KeyError(
                                    "Label " + user_index + " has been used already for another HSQC peak! Go back to the HSQC file"
                                                            "and change the label in the following line:\n" + line)
                            contents[i] = "\t%s\t%s\t%s\n" % ("X" + user_index + "NX-HX", Nreson, Hreson)
                        else:
                            contents[i] = "\t%s\t%s\t%s\n" % ("X" + str(index) + "NX-HX", Nreson, Hreson)
                            index += 1

        # Now write the modified contents to a new file
        if OUT_fname:
            out_fname = OUT_fname
        else:
            out_fname = HSQC_FILE + "_annotated"
        with open(out_fname, 'w') as f:
            for line in contents:
                f.write(line)

    def find_nonoverlapping_groups_in_HSQC(self, remaining_HSQC_lines, rtolH, rtolN, keep_closest=False,
                                           get_only_AAIG_overlaps=False, skip_sidechains=False):
        """
        Method used only by self.uniquify_HSQC_AAIGs.
        :param rtolH:
        :param rtolN:
        :param keep_closest:
        :return:
        nonoverlapping_groups_lines_and_tolerances:    list with the line of each non-overlapping group of the HSQC spectrum along with the rtolH and rtolN
        remaining_groups_lines:                        list with the lines of the Root spectrum that contain overlapping groups with the current torelances
        """

        print("DEBUG: remaining_HSQC_contents=", remaining_HSQC_lines)
        nonoverlapping_groups_lines_and_tolerances = []
        remaining_groups_lines = [] # just for debugging
        AAIG_overlaps_dict = defaultdict(list)  # AAIG signature -> list of overlapping AAIG signatures according to the given tolH, tolN

        if keep_closest == False:
            for HSQC1_line in remaining_HSQC_lines:
                HSQC1_words_list = HSQC1_line.split()
                found_overlap = False
                try:
                    if skip_sidechains and HSQC1_words_list[0][-3:] != "N-H":  # if it's not a backbone amide, skip it
                        continue
                    HSQC1_AAIG_signature = HSQC1_words_list[0].replace("-", "")
                    HSQC1_Nreson = float(HSQC1_words_list[1])
                    HSQC1_HNreson = float(HSQC1_words_list[2])
                    for HSQC2_line in remaining_HSQC_lines:
                        if HSQC1_line == HSQC2_line:  # in case we are looking at the same line
                            continue
                        # print "DEBUG: HSQC2_line=",HSQC2_line
                        HSQC2_words_list = HSQC2_line.split()
                        try:
                            HSQC2_AAIG_signature = HSQC2_words_list[0].replace("-", "")
                            HSQC2_Nreson = float(HSQC2_words_list[1])
                            HSQC2_HNreson = float(HSQC2_words_list[2])
                            if abs(HSQC1_HNreson-HSQC2_HNreson) <= rtolH and abs(HSQC1_Nreson-HSQC2_Nreson) <= rtolN:
                                remaining_groups_lines.append(HSQC1_line)
                                found_overlap = True
                                AAIG_overlaps_dict[HSQC1_AAIG_signature].append(HSQC2_AAIG_signature)
                                break
                        except (ValueError, IndexError):
                            # print "WARNING: the 4th and 5th elements of the following HSQC file line are not numbers:"
                            # print "Root file line:", HSQC_line
                            continue
                    if found_overlap == False:  # if no noverlapping line was found for this line then save it
                        nonoverlapping_groups_lines_and_tolerances.append((HSQC1_line, rtolH, rtolN))
                except (ValueError, IndexError):
                    # print "WARNING: the 2nd and 3rd elements of the following HSQC file line are not numbers:"
                    # print "Root file line:", HSQC_line
                    continue

        if get_only_AAIG_overlaps:
            return AAIG_overlaps_dict
        # print "DEBUG: len(nonoverlapping_groups_lines_and_tolerances)=", len(nonoverlapping_groups_lines_and_tolerances)
        print("DEBUG: nonoverlapping_groups_lines_and_tolerances =", nonoverlapping_groups_lines_and_tolerances)
        print("DEBUG: remaining_groups_lines = ", remaining_groups_lines)
        return nonoverlapping_groups_lines_and_tolerances

    def uniquify_HSQC_AAIGs(self):
        """
        Method that gradually reduces N and HN tolerances to find more non-overlapping groups in HSQC spectrum
        and assign a unique label to each peak (aka AAIG, line) of the HSQC spectrum. This method doesn't deal
        with unlabeled peaks, namely '?-?'. Use instead annotate_HSQC().
        :return:
        """
        HSQC_lines, patched_residues = self.uniquify_spectrum_lines(self.HSQC_FILE)    # ignore patched_residues here
        remaining_HSQC_contents = []  # list of Root spectrum lines with overlapping groups
        nonoverlapping_HSQC_contents_and_tolerances = []  # list with the lines of the Root spectrum that contain non-overlapping groups and the associated rtoH, rtolN
        AAIG_overlaps_dict = self.find_nonoverlapping_groups_in_HSQC(HSQC_lines, rtolH=0.02, rtolN=0.2,
                                                                     get_only_AAIG_overlaps=True)  # get overlaps between HSQC AAIGs
        for rtolH, rtolN in zip([0.02, 0.01, 0.01, 0.005], [0.2, 0.1, 0.05, 0.015]):
            # for rtolH, rtolN in zip([0.02, 0.01, 0.01], [0.2, 0.1, 0.05]):
            # print "DEBUG: len(HSQC_lines)=", len(HSQC_lines)
            # FIND NON-OVERLAPPING RESONANCE GROUPS IN THE HSQC SPECTRUM
            nonoverlapping_HSQC_contents_and_tolerances.extend(
                self.find_nonoverlapping_groups_in_HSQC(HSQC_lines, rtolH, rtolN))
            print("DEBUG: rtolH=", rtolH, "rtolN=", rtolN, "nonoverlapping_HSQC_contents_and_tolerances=", nonoverlapping_HSQC_contents_and_tolerances)
            for triplet in nonoverlapping_HSQC_contents_and_tolerances:
                line = triplet[0]
                try:
                    HSQC_lines.remove(line)
                except ValueError:
                    continue

        # NOTE: up to this point HSQC labels are not changed
        # Load all non-overlapping AAIGs into AAIG objects
        renamed_AAIG_labels = {}    # to keep track of the invalid AAIG names
        name_index = 1000  # designates 'X1000', number starts from 1000 to avoid conflicts with manually assigned AAIG
        for line, rtolH, rtolN  in nonoverlapping_HSQC_contents_and_tolerances:
            try:
                label, Nreson, HNreson = line.split()
            except ValueError:
                raise Exception("Line %s failed does not have laber, Nreson, or HNreson." % line)
            # ATTENTION: I cannot have '' in Nnames because it will always prefer that. Therefore this is not valid label 'N107-?', but this is 'N107N-H'
            Nnames = ['N', 'ND2', 'ND2', 'NE2', 'NE2', 'NE', 'NE1', '?', 'NX']          # all alternative valid N names
            HNnames = ['H', 'HD21', 'HD22', 'HE21', 'HE22', 'HE', 'HE1', '?', 'HX']    # and the respective HN names
            p = regex.compile(r"^([^\s]+)(\L<Nname>-\L<HNname>)$", regex.X, Nname=Nnames, HNname=HNnames)
            m = p.search(label)
            if m:
                AAIG_name = m.group(1)
                NH_name = m.group(2)
                # print "DEBUG: AAIG_name=%s , NH_name=%s" % (AAIG_name, NH_name)
                if NH_name == "-?": # special case, e.g. label = 'A107-?'
                    NH_name = "?-?"
            else:
                ColorPrint("WARNING: HSQC line does not have a valid amide label:\t%s" % line, "WARNING")
                ColorPrint("Renaming AAIG and amide to '?'", "WARNING")
                AAIG_name, NH_name = '?', '?-?'

            # The following block renames only the AAIG_name, the NH_name is preserved
            if not re.match("^[ACDEFGHIKLMNPQRSTVWYX][0-9]+$", AAIG_name):
                ColorPrint("WARNING: HSQC line does not have a valid AAIG name:\t" + line, "WARNING")
                if AAIG_name in list(renamed_AAIG_labels.keys()): # if a valid AAIG name has been already assigned to this invalid AAIG name
                    AAIG_name = renamed_AAIG_labels[AAIG_name]
                else:
                    if re.match("^[0-9]+$", AAIG_name) and "X"+AAIG_name not in list(renamed_AAIG_labels.values()):   # special case
                        renamed_AAIG_labels[AAIG_name] = "X"+AAIG_name
                        ColorPrint("Renaming %s to %s" % (AAIG_name, renamed_AAIG_labels[AAIG_name]), "WARNING")
                        AAIG_name = renamed_AAIG_labels[AAIG_name]  # replace the invalid AAIG name with the valid
                    else:
                        while "X%i" % name_index in list(renamed_AAIG_labels.values()):
                            name_index += 1
                        renamed_AAIG_labels[AAIG_name] = "X%i" % name_index
                        ColorPrint("Renaming %s to %s" % (AAIG_name, renamed_AAIG_labels[AAIG_name]), "WARNING")
                        AAIG_name = renamed_AAIG_labels[AAIG_name]  # replace the invalid AAIG name with the valid
                        name_index += 1
                if NH_name == "?-?":    # finally rename the unknown N-H names
                    NH_name = "NX-HX"

            # Add one AAIG object for each line
            aaig = AAIG(signature=AAIG_name + NH_name.replace('-', ''),
                        NH_name=NH_name.replace('-', ''),
                        overlapping_AAIGs=AAIG_overlaps_dict[AAIG_name + NH_name.replace('-', '')],
                        user_assignment=label,
                        H_bin_length=self.H_bin_length,
                        C_bin_length=self.C_bin_length,
                        bandwidth=self.bandwidth)
            peak = Peak(j_AAIG_signature=aaig.signature, j_AAIG_name=remove_NH_suffix(aaig.signature),
                        j_Nreson=round(float(Nreson), 3), j_Nname=NH_name.split('-')[0],
                        j_HNreson=round(float(HNreson), 3), j_HNname=NH_name.split('-')[1])
            aaig.__add_peak__(peak)
            self.__add_AAIG__(signature=aaig.signature, AAIG=aaig)


    def match_j_AAIG_resonances(self,
                                query_lines,
                                tolHN,
                                tolN,
                                query_fname,
                                spectrum_type,
                                keep_CH_assignments=False,
                                debug=False,
                                exclude_alternatives=False,
                                get_traindata=False):
        """
        Auxiliary method for copy_AAIGnames_from_HSQC_spectrum(). It allows the latter to work also for HNNH spectra.

        :param query_lines:
        :param tolHN:
        :param tolN:
        :param query_fname:
        :return:
        """

        if keep_CH_assignments: # Save the C-H labels to add them later
            if keep_CH_assignments:
                CH_labels_dict = {}  # C-H label -> all numbers in the line in str format
                for i in range(len(query_lines)):
                    words = query_lines[i].split()
                    CH_labels_dict[tuple(words[1:])] = "-".join(words[0].split('-')[:2])
                    query_lines[i] = "\t".join(['?-?-?-?'] + words[1:])  # replace the label with questionmarks to be readable by match_j_AAIG_resonances()
        if debug:
            for i in range(len(query_lines)):
                words = query_lines[i].split()
                words.append('nan') # add an extra column for the distance of multiple alternative matched j_AAIG from the HSQC AAIG
                query_lines[i] = "\t".join(words)

        for AAIG in self.__get_all_AAIGs__():
        # for AAIG in self.__get_all_nonoverlapping_AAIGs__():  # TEMPORARILY DEACTIVATED
            peak = AAIG.__get_all_peaks__()[0]     # by convention HSQC has only one peak per AAIG
            try:
                NH_AAIG_signature = peak.j_AAIG_signature     # signature is f.e. N107ND2-HD21, or N107N-H
                # print("DEBUG: NH_AAIG_signature =%s" % NH_AAIG_signature)
                HSQC_Nreson = peak.j_Nreson
                HSQC_HNreson = peak.j_HNreson
                for q_index in range(0, len(query_lines)):
                    query_line = query_lines[q_index]
                    # print("DEBUG: query_line=",query_line)
                    try:
                        query_words_list = query_line.split()
                        # create one peak object for the current query spectrum line to facilitate comparison operations
                        query_peak = Peak(i_Nreson=float(query_words_list[2]), i_HNreson=float(query_words_list[1]),
                                          j_Nreson=float(query_words_list[3]), j_HNreson=float(query_words_list[4]))
                        # First remove the garbage that are far away from any peak, and then keep the closest match
                        if peak.does_peak_match(query_peak, tolH=tolHN, tolN=tolN, attypes=['N', 'HN']):
                            # print "DEBUG: ?-?-"+NH_AAIG_signature+"\t"+query_words_list[1]+"\t"+query_words_list[2]+"\t"+query_words_list[3]+"\t"+query_words_list[4]
                            #  If the line does not have an AAIG already, assign the currect AAIG to it
                            if re.search("^\s*\?-\?-\?-\?\s+", query_lines[q_index]):
                                if debug:
                                    if get_traindata:
                                        delta_HSQC = "+".join([str(d) for d in peak.distance(query_peak, attypes=['N', 'HN'], get_distlist=True)])
                                    else:
                                        delta_HSQC = peak.distance(query_peak, attypes=['N', 'HN'])
                                    query_words_list[6] = str(delta_HSQC)  # chage only the distance, replace nan
                                    query_lines[q_index] = " ".join(query_words_list + ["\n"])
                                # now substitude ? for the right label
                                query_lines[q_index] = re.sub("\?-\?-\?-\?", NH_AAIG_signature,
                                                                 query_lines[q_index]).replace("\\", "")
                            else:   # Else if it has already been assigned an AAIG, check if the currect AAIG has closer
                                    # N & HN resonances
                                if exclude_alternatives:  # this line has alternative assignments, hence exclude it!
                                    query_words_list[0] = "EXCLUDE"  # instead of distance, mark the line to be excluded
                                    query_lines[q_index] = " ".join(query_words_list + ["\n"])
                                    continue
                                # print "DEBUG: Already assigned AAIG. Line: ", query_lines[q_index]
                                mo = re.search('^\s*([A-Za-z0-9]+N[DE12X]*H[DE12X]*)\s+', query_lines[q_index])
                                if mo:
                                    previous_NH_AAIG_signature = mo.group(1)  # the AAIG that is already assigned to this peak from a previous round
                                    # print("DEBUG: previous_NH_AAIG_signature = %s" % previous_NH_AAIG_signature)
                                    for tmp_AAIG in self.__get_all_AAIGs__(): # iterate over all AAIGs until you find
                                                                            # previous_AAIG, in order to retrieve its N & HN resonances
                                    # for tmp_AAIG in self.__get_all_nonoverlapping_AAIGs__(): # iterate over all AAIGs until you find
                                                                            # previous_AAIG, in order to retrieve its N & HN resonances
                                        tmp_peak = tmp_AAIG.__get_all_peaks__()[0]  # by convention HSQC has only one peak per AAIG
                                        try:
                                            tmp_NH_AAIG_signature = tmp_peak.j_AAIG_signature
                                            # print "DEBUG: tmp_NH_AAIG_signature=", tmp_NH_AAIG_signature, "previous_NH_AAIG_signature=", previous_NH_AAIG_signature
                                            if tmp_NH_AAIG_signature == previous_NH_AAIG_signature:  # we found the N & HN resonances of the assigned AAIG
                                                # print "DEBUG: we found the N & HN resonances of the assigned AAIG:", tmp_HSQC_words_list
                                                previous_HSQC_Nreson = tmp_peak.j_Nreson
                                                previous_HSQC_HNreson = tmp_peak.j_HNreson
                                                # if the current AAIG has N & HN resonances closer to the query line than the previously
                                                # assigned N & HN resonances, change the AAIG
                                                # print "DEBUG: query_Nreson=", query_Nreson, "previous_HSQC_Nreson=", previous_HSQC_Nreson, "HSQC_Nreson=", HSQC_Nreson
                                                # print "DEBUG: query_HNreson=", query_HNreson, "previous_HSQC_HNreson=", previous_HSQC_HNreson, "HSQC_HNreson=", HSQC_HNreson
                                                delta_previous_HSQC = tmp_peak.distance(query_peak, attypes=['N', 'HN'])
                                                delta_HSQC = peak.distance(query_peak, attypes=['N', 'HN'])
                                                # print "DEBUG: delta_previous_HSQC=", delta_previous_HSQC, "delta_HSQC=", delta_HSQC
                                                if delta_previous_HSQC > delta_HSQC:
                                                    # print "DEBUG: substituting previous_NH_AAIG_signature", previous_NH_AAIG_signature, " for ", NH_AAIG_signature, " in ", query_lines[q_index]
                                                    if debug:   # I presume that the file contains intensity column
                                                        if get_traindata:
                                                            delta_HSQC = "+".join([str(d) for d in peak.distance(query_peak, attypes=['N', 'HN'], get_distlist=True)])
                                                        query_lines.append(" ".join(
                                                            [NH_AAIG_signature, query_words_list[1], query_words_list[2],
                                                             query_words_list[3], query_words_list[4],
                                                             query_words_list[5], str(delta_HSQC), "\n"]
                                                        ))
                                                    elif len(query_words_list) >= 6:  # >= to catch debug, too (extra word=distance)
                                                        query_lines[q_index] = " ".join(
                                                            [NH_AAIG_signature, query_words_list[1],
                                                             query_words_list[2],
                                                             query_words_list[3], query_words_list[4],
                                                             query_words_list[5], "\n"]
                                                        )
                                                    else:
                                                        query_lines[q_index] = " ".join(
                                                            [NH_AAIG_signature, query_words_list[1],
                                                             query_words_list[2],
                                                             query_words_list[3], query_words_list[4], "\n"]
                                                        )
                                                    # query_lines[q_index] = re.sub(previous_NH_AAIG_signature, NH_AAIG_signature, query_lines[q_index])
                                                # otherwise keep the previous AAIG
                                                # ATTENTION: THIS IS NOT ENTIRELY CORRECT, BECAUSE THE N RESONANCE MAY BE CLOSER BUT THE H NOT!
                                                break
                                        except (ValueError, IndexError):
                                            # print "WARNING: the 2nd and 3rd elements of the following HSQC file line are not numbers:"
                                            # print "Root file line:", HSQC_line
                                            continue
                                else:
                                    ColorPrint("ERROR: wrong line (modified or not) in %s" % query_fname, "FAIL")
                                    print(query_lines[q_index])
                                    sys.exit(1)

                    except (ValueError, IndexError):
                        # print "WARNING: the 4th and 5th elements of the following query file line are not numbers:"
                        # print "Query file line:", query_line
                        continue
            except (ValueError, IndexError):
                # print "WARNING: the 2nd and 3rd elements of the following HSQC file line are not numbers:"
                # print "Root file line:", HSQC_line
                continue

        # print "DEBUG: ", spectrum_type, " contents: "
        # print "".join(query_lines)
        if spectrum_type == "HNNH":
            new_query_lines = []
        else:
            new_query_lines = query_lines
        # In case of HNNH the num.list file will be written twice, but this doesn't hurt anybody :)
        with open(os.path.splitext(query_fname)[0] + "num.list", 'w') as f:
            for line in query_lines:
                try:
                    if exclude_alternatives and line.startswith("EXCLUDE"):  # remove peaks that were marked to be excluded
                        continue
                    # print "DEBUG: saving line:", line
                    words = line.split()
                    if debug:
                        appendix = "\t" + "\t".join(words[5:])  # an extra word (distance) exists, already converted to string
                    elif len(words) == 6:     # if the file has intensity column
                        appendix = "\t" + words[5]
                    else:
                        appendix = ""

                    if keep_CH_assignments:
                        CH_label = CH_labels_dict[tuple(words[1:])]
                    else:
                        CH_label = "?-?"
                    # print("DEBUG: words[0]=", words[0])
                    if not "?" in words[0]: # if an j_AAIG_signature has been assigned
                        new_line = "?-?-%s\t\t%s\t%s\t%s\t%s%s\n" % \
                                   (words[0], words[1], words[2], words[3], words[4], appendix)
                        new_fline = "%s-%s\t\t%s\t%s\t%s\t%s%s\n" % \
                                (CH_label, dash_to_NH(words[0]), words[1], words[2], words[3], words[4], appendix)
                    else:
                        new_line = "?-?-?-?\t\t%s\t%s\t%s\t%s%s\n" % \
                                   (words[1], words[2], words[3], words[4], appendix)
                        new_fline = "%s-?-?\t\t%s\t%s\t%s\t%s%s\n" % \
                                (CH_label, words[1], words[2], words[3], words[4], appendix)
                    f.write(new_fline)
                    if spectrum_type == "HNNH":
                        new_query_lines.append(new_line)
                except IndexError:
                    ColorPrint("WARNING: Discarding " + spectrum_type + " line:\n" + line, "WARNING")
                    # print "DEBUG: line=", line
                    # print "DEBUG: copy_aaindices_from_HSQC_spectrum point 2 words=", words
                    continue

        return new_query_lines

    def match_i_AAIG_resonances(self, query_lines, tolHN, tolN, query_fname, spectrum_type, debug=False,
                                exclude_alternatives=True, get_traindata=False):
        """
        Auxiliary method for copy_AAIGnames_from_HSQC_spectrum(). It allows the latter to work also for HNNH spectra,
        in which case this method is called after match_NH_resonances() and instead of C,H it uses the first HN, N
        resonances of NHHN file.

        :param query_lines:
        :param tolHN:
        :param tolN:
        :param query_fname:
        :param spectrum_type:
        :param debug:   if True, all HCNH matching lines will be saved, not only the closest to an AAIG in HSQC
        :return:
        """

        if debug:
            for i in range(len(query_lines)):
                words = query_lines[i].split()
                # I presume that already contains the distance column between j_AAIG and HSQC AAIG.
                words = words[:-1] + ['nan', words[-1]] # add an extra column for the distance of multiple alternative matched i_AAIG from the HSQC AAIG
                query_lines[i] = "\t".join(words)

        # TODO: not sure if this works correctly with signatures.
        j_AAIG_labels = [l.split('-')[-1]
                         for l in [ql.split()[0]
                                   for ql in query_lines]]  # save the N-HN labels to add them later
        for AAIG in self.__get_all_AAIGs__():
        # for AAIG in self.__get_all_nonoverlapping_AAIGs__():  # TEMPORARILY DEACTIVATE
            peak = AAIG.__get_all_peaks__()[0]     # by convention HSQC has only one peak per AAIG
            try:
                # HSQC has only one pair of HN, N resonances, therefore i refers to HNNH
                i_AAIG_signature = peak.j_AAIG_signature     # signature is f.e. N107ND2-HD21, or N107N-H
                print("DEBUG: i_AAIG_signature =%s" % i_AAIG_signature)
                for q_index in range(0, len(query_lines)):
                    query_line = query_lines[q_index]
                    try:
                        query_words_list = query_lines[q_index].split()
                        # create one peak object for the current query spectrum line to facilitate comparison operations
                        query_peak = Peak(i_Nreson=float(query_words_list[2]), i_HNreson=float(query_words_list[1]))
                        # First remove the garbage that are far away from any peak
                        if peak.does_peak_match(query_peak, tolHN=tolHN, tolN=tolN, attypes=['Njtoi', 'HNjtoi']):
                            #  If the line does not have a i_AAIG_signature already, assign the currect i_AAIG_signature to it
                            label = query_words_list[0] # the whole label ".*-.*-.*-.*"
                            if label[:4] == '?-?-':
                                if debug:
                                    if get_traindata:
                                        delta_HSQC = "+".join([str(d) for d in peak.distance(query_peak, attypes=['Njtoi', 'HNjtoi'], get_distlist=True)])
                                    else:
                                        delta_HSQC = peak.distance(query_peak, attypes=['Njtoi', 'HNjtoi'])
                                    query_words_list = query_lines[q_index].split()
                                    query_lines[q_index] = " ".join(
                                        [query_words_list[0], query_words_list[1], query_words_list[2],
                                         query_words_list[3], query_words_list[4],
                                         query_words_list[5], str(delta_HSQC), query_words_list[7],
                                         "\n"])
                                # now substitude ? for the right label
                                query_lines[q_index] = re.sub(label.replace('?', '\?'), i_AAIG_signature,
                                                                 query_lines[q_index]).replace("\\", "")

                            else:   # Else if it has already been assigned an AAIG, check if the currect AAIG has closer
                                    # H & C resonances
                                if exclude_alternatives:  # this line has alternative assignments, hence exclude it!
                                    query_words_list[0] = "EXCLUDE"  # instead of distance, mark the line to be excluded
                                    query_lines[q_index] = " ".join(query_words_list + ["\n"])
                                    continue
                                mo = re.search('^([A-Za-z0-9]+N[DE12X]*H[DE12X]*)$', label)
                                if mo:
                                    previous_i_AAIG_signature = mo.group(1)  # the AAIG that is already assigned to this peak from a previous iteration
                                    print("DEBUG: previous_i_AAIG_signature = %s" % previous_i_AAIG_signature)
                                    for tmp_AAIG in self.__get_all_AAIGs__(): # iterate over all AAIGs until you find
                                                                            # previous_AAIG, in order to retrieve its N & HN resonances
                                    # for tmp_AAIG in self.__get_all_nonoverlapping_AAIGs__(): # TEMPORARILY DEACTIVATED. iterate over all AAIGs until you find
                                    #                                     # previous_AAIG, in order to retrieve its N & HN resonances
                                        tmp_peak = tmp_AAIG.__get_all_peaks__()[0]  # by convention HSQC has only one peak per AAIG
                                        try:
                                            tmp_i_AAIG_signature = tmp_peak.j_AAIG_signature
                                            if tmp_i_AAIG_signature == previous_i_AAIG_signature:  # we found the N & HN resonances of the assigned AAIG
                                                # If the current AAIG has N & HN resonances closer to the query line than the previously
                                                # assigned N & HN resonances, change the AAIG
                                                delta_previous_HSQC = tmp_peak.distance(query_peak, attypes=['Njtoi', 'HNjtoi'])
                                                delta_HSQC = peak.distance(query_peak, attypes=['Njtoi', 'HNjtoi'])
                                                if delta_previous_HSQC > delta_HSQC:
                                                    if debug:   # I presume that contains intensity column
                                                        if get_traindata:
                                                            delta_HSQC = "+".join([str(d) for d in peak.distance(query_peak, attypes=['Njtoi', 'HNjtoi'], get_distlist=True)])
                                                        j_AAIG_labels.append(j_AAIG_labels[q_index])
                                                        query_lines.append(" ".join(
                                                            [i_AAIG_signature, query_words_list[1], query_words_list[2],
                                                             query_words_list[3], query_words_list[4],
                                                             query_words_list[5], str(delta_HSQC), query_words_list[7], "\n"]
                                                        ))
                                                    elif len(query_words_list) == 6:
                                                        query_lines[q_index] = " ".join(
                                                            [i_AAIG_signature, query_words_list[1],
                                                             query_words_list[2],
                                                             query_words_list[3], query_words_list[4],
                                                             query_words_list[5], "\n"]
                                                        )
                                                    else:
                                                        query_lines[q_index] = " ".join(
                                                            [i_AAIG_signature, query_words_list[1],
                                                             query_words_list[2],
                                                             query_words_list[3], query_words_list[4], "\n"]
                                                        )
                                                    print("DEBUG: changed label of line to:",query_lines[q_index])
                                                    # query_lines[q_index] = re.sub(previous_i_AAIG_signature, i_AAIG_signature, query_lines[q_index])
                                                # otherwise keep the previous AAIG
                                                # ATTENTION: THIS IS NOT ENTIRELY CORRECT, BECAUSE THE N RESONANCE MAY BE CLOSER BUT THE H NOT!
                                                break
                                        except (ValueError, IndexError):
                                            # print "WARNING: the 2nd and 3rd elements of the following HSQC file line are not numbers:"
                                            # print "Root file line:", HSQC_line
                                            continue
                                else:
                                    ColorPrint("ERROR: wrong line (modified or not) in %s" % query_fname, "FAIL")
                                    print(query_lines[q_index])
                                    sys.exit(1)

                    except (ValueError, IndexError):
                        # print "WARNING: the 4th and 5th elements of the following query file line are not numbers:"
                        # print "Query file line:", query_line
                        continue
            except (ValueError, IndexError):
                # print "WARNING: the 2nd and 3rd elements of the following HSQC file line are not numbers:"
                # print "Root file line:", HSQC_line
                continue

        # print "DEBUG: ", spectrum_type, " contents: "
        # print "".join(query_lines)
        # sys.exit(1)
        new_query_lines = []
        with open(os.path.splitext(query_fname)[0] + "num.list", 'w') as f:
            for q_index in range(len(query_lines)):
                try:
                    line = query_lines[q_index]
                    if exclude_alternatives and line.startswith("EXCLUDE"):  # remove peaks that were marked to be excluded
                        continue
                    j_AAIG = j_AAIG_labels[q_index]
                    # print "DEBUG: saving line:", line
                    word_list = line.split()
                    if debug:   # I presume that the file has intensity column
                        appendix = "\t" + "\t".join(word_list[5:])  # add intensity and all distances
                    elif len(word_list) == 6:     # if the file has intensity column
                        appendix = "\t" + word_list[5]
                    else:
                        appendix = ""

                    if "?" not in word_list[0]: # if an i_AAIG_signature has been assigned
                        new_line = "%s-%s\t\t%s\t%s\t%s\t%s%s\n" % \
                                (word_list[0], j_AAIG, word_list[1], word_list[2], word_list[3], word_list[4], appendix)
                        new_fline = "%s-%s\t\t%s\t%s\t%s\t%s%s\n" % \
                                (dash_to_NH(word_list[0]), dash_to_NH(j_AAIG), word_list[1], word_list[2], word_list[3], word_list[4], appendix)
                    else:
                        new_line = "?-?-%s\t\t%s\t%s\t%s\t%s%s\n" % \
                                (j_AAIG, word_list[1], word_list[2], word_list[3], word_list[4], appendix)
                        new_fline = "?-?-%s\t\t%s\t%s\t%s\t%s%s\n" % \
                                (dash_to_NH(j_AAIG), word_list[1], word_list[2], word_list[3], word_list[4], appendix)
                    f.write(new_fline)
                    new_query_lines.append(new_line)
                except IndexError:
                    ColorPrint("WARNING: Discarding " + spectrum_type + " line:\n" + line, "WARNING")
                    # print "DEBUG: line=", line
                    # print "DEBUG: copy_aaindices_from_HSQC_spectrum point 2 word_list=", word_list
                    continue

        return new_query_lines

    def copy_AAIGnames_from_HSQC_spectrum(self,
                                          spectrum4D_fname,
                                          spectrum4D_type,
                                          tolHN=0.04,
                                          tolN=0.4,
                                          TOCSY_sidechain_resonances_list=None,
                                          keep_CH_assignments=False,
                                          add_dash_to_NH=False,
                                          debug=False,
                                          exclude_alternatives=False,
                                          get_traindata=False,
                                          uniquify_lines=True):
        """
        Method that finds the closest TOCSY/NOESY peaks to each HSQC spectrum group (namely, each HSQC line). If the lowest
        distance of a TOCSY/NOESY peak from any HSQC group is greater than (0.04, 0.4) then it will be considered as garbage
        and will be discarded from TOCSYnum.list file. It writes the fully labeled spectrum into a file with the
        suffix "num.list", e.g. if spectrum4D_fname="hIL9_HCNH.list" then "hIL9_HCNHnum.list".

        :param spectrum4D_fname: a sparky list file (4D-TOCSY, 4D-NOESY, HNNH) to which to copy the HSQC AAIG names
        :param spectrum4D_type:   can be 'TOCSY', 'HCNH', 'HNNH'
        :param tolHN:
        :param tolN:
        :param TOCSY_sidechain_resonances_list: this is useful to discard the side-chain AAIGs from NOESY
        :param keep_CH_assignments:
        :param add_dash_to_NH:
        :param debug:
        :param exclude_alternatives:
        :param get_traindata:
        :return query_lines:    list of the lines of query_fname with N-H labels (w3,w4) copied from HSQC.
        :return TOCSY_sidechain_resonances_list: this list is returned only if spectrum4D_type == 'TOCSY'. It contains
                                                 unassigned TOCSY lines that correspond to sidechains. The resonances
                                                 must be later removed from NOESY.
        """

        query_spec = Spectrum(AAIGs_dict={}, spectrum_type=spectrum4D_type)
        if uniquify_lines:
            query_lines, patched_residues = \
                query_spec.uniquify_spectrum_lines(spectrum4D_fname) # ignore patched residues here
        else:   # Assuming that the spectrum is uniquified; used in nmr_pipeline (NEW CODE).
            if spectrum4D_fname.endswith(".csv"):
                original_query_contents = [l for l in csv.reader(open(spectrum4D_fname))]
                if original_query_contents[0][0] == "peak_index":
                    original_query_contents = [o[1:] for o in original_query_contents]
                if original_query_contents[0][-1] == "scaled_Intensity":
                    original_query_contents = [o[:-1] for o in original_query_contents]
                query_lines = [" ".join(l) for l in original_query_contents[1:]]
            else:
                with open(spectrum4D_fname, 'r') as f:
                    query_lines = f.readlines()  # contents of original spectrum4D_fname (4D TOCSY or 4D NOESY) in 5 column format (name H C N HN)

        # Check if the query fname has any label
        if keep_CH_assignments == False:
            for line in query_lines:
                assert line.split()[0] == '?-?-?-?', \
                    ColorPrint("ERROR: the input %s file %s had labeled peaks. You must provide "
                                "a file without labels in order to do HSQC-%s matching!"
                               % (spectrum4D_fname, spectrum4D_type, spectrum4D_type), "FAIL")

        query_lines = self.match_j_AAIG_resonances(query_lines,
                                                   tolHN,
                                                   tolN,
                                                   spectrum4D_fname,
                                                   spectrum4D_type,
                                                   keep_CH_assignments=keep_CH_assignments,
                                                   debug=debug,
                                                   exclude_alternatives=exclude_alternatives,
                                                   get_traindata=get_traindata)

        if spectrum4D_type == "HNNH": # this time match the first HN,N resonances (w1,w2) in the HNNH file
            query_lines = self.match_i_AAIG_resonances(query_lines, tolHN, tolN, spectrum4D_fname, spectrum4D_type,
                                                       debug=debug, exclude_alternatives=exclude_alternatives,
                                                       get_traindata=get_traindata)
                                                       # keep_CH_assignments=keep_CH_assignments)   # NOT APPLICABLE???
        self.numlist_files[spectrum4D_type] = os.path.abspath(os.path.splitext(spectrum4D_fname)[0] + "num.list")

        if spectrum4D_type == "TOCSY":
            # save unassigned TOCSY lines that correspond to sidechains. The resonances must be later removed from NOESY
            TOCSY_sidechain_resonances_list = []
            for qline in query_lines:
                try:
                    # print "DEBUG: qline.split()[0]:", qline.split()[0]
                    if qline.split()[0] == "?-?-?-?":
                        TOCSY_sidechain_resonances_list.append(qline)
                except IndexError:
                    continue
            for TOCSY_line in TOCSY_sidechain_resonances_list:  # TEMPORARILY REMOVE ALL "?-?-?-?" LINES (NOT ALL OF THEM ARE SIDE CHAIN
                                                                # RESONANCES FROM TOCSY FILE CONTENTS)
                query_lines.remove(TOCSY_line)   # now you won't see TOCSY index groups like "?-?-?-?" in the connectivities
                                                    # and aa type prediction files
            return query_lines, TOCSY_sidechain_resonances_list
        elif spectrum4D_type == "HCNH":   # the 4D-NOESY
            lines2remove_set = set()
            for qline in query_lines:  # FIND ALL LINES WITH "?-?-?-?"
                try:
                    # print "DEBUG: qline.split()[0]:",qline.split()[0]
                    if qline.split()[0] == "?-?-?-?":
                        lines2remove_set.add(qline)
                except IndexError:
                    continue
            # remove the sidechain resonances from NOESY
            # for TOCSY_line in TOCSY_sidechain_resonances_list:
            for qline in lines2remove_set:  # REMOVE ALL LINES WITH "?-?-?-?"
                # print "DEBUG: removing side chain line from NOESY: ",TOCSY_line
                try:
                    # query_lines.remove(TOCSY_line)
                    query_lines.remove(qline)
                except ValueError:
                    # print "DEBUG: side chain line not present in NOESY!"
                    continue
            return query_lines, []
        elif spectrum4D_type == "HNNH":
            if add_dash_to_NH:  # convert E45NH-F46NH to E45N-H-F46N-H
                for i in range(len(query_lines)):
                    words = query_lines[i].split()
                    label = words[0]
                    i_AAIG_signature, j_AAIG_signature = label.split('-')[0], label.split('-')[-1]
                    if i_AAIG_signature == '?':
                        i_AAIG_signature = "?-?"
                    if j_AAIG_signature == '?':
                        j_AAIG_signature = "?-?"
                    words[0] = "%s-%s\t" % (dash_to_NH(i_AAIG_signature), dash_to_NH(j_AAIG_signature))    # replace the 1st word
                    query_lines[i] = "\t".join(words) + "\n"
            # TODO: think if removing side chains like in NOESY above, applies in HNNH, too.
            return query_lines, []

    @staticmethod
    def guess_starting_resid(HSQC_FILE, fasta, NHmap=None):
        """
            FUNCTION that tries to guess the starting resid from the labels present in the HSQC HSQC file. The function will return
        """
        
        if not NHmap:
            protein_sequence_list = get_protein_sequence(fasta)
        else:
            absolute_AAIGmatches_alignment_list, protein_sequence_list = \
                Alignment.read_NHmap_file(NHmap, get_protein_alignment=True)
        labeled_residues_list = []  # list of the form: [(resname, resid), ...]
        with open(HSQC_FILE, 'r') as f:
            contents = f.readlines()
            index = 1
            for i in range(len(contents)):
                line = contents[i]
                mo = re.match(".*\s+([0-9.]+)\s+([0-9.]+)\s*", line) # the last 2 columns must be numbers
                if mo:
                    Nreson = mo.group(1)
                    Hreson = mo.group(2)
                    mp = re.match("^\s*([ARNDCEQGHILKMFPSTWYV])([0-9]+)N-H\s+[0-9.]+\s+[0-9.]+\s*", line)    # consider only N-H mapped AAIG, ignore NX-HX, etc.
                    if mp: # if the label matches the conditions, save it
                        aa_type = mp.group(1)
                        resid = int(mp.group(2))
                        labeled_residues_list.append((aa_type, int(resid)))
        
        if len(labeled_residues_list) == 0:
            print(bcolors.BOLDBLUE + "I could not guess the starting residue number from the labels in the file. Therefore I will relabel", end=' ')
            print("everything starting from 1. If you want to keep your labels please check again for consistency with the", end=' ')
            print("protein sequence." + bcolors.ENDBOLD)
            return None
        labeled_residues_list.sort(key=itemgetter(1))
        existing_resnames = [d[0] for d in labeled_residues_list]
        first_resid = labeled_residues_list[0][1]   # the number in the label with the smallest resid
        existing_resids = [d[1]-first_resid for d in labeled_residues_list] # subtract from the labeled resids the first so that they start from 0
        possible_starting_positions = []
        # print "DEBUG: existing_resids=", existing_resids
        # print "DEBUG: existing_resnames=", existing_resnames
        # print "DEBUG: protein_sequence_list=", protein_sequence_list
        for start in range(len(protein_sequence_list)):
            # print "DEBUG: Checking starting position", start
            try:
                subsequence = protein_sequence_list[start:]
                # matches = 0
                unmatched = 0
                # print "DEBUG: subsequence=", subsequence
                for pos,resname in zip(existing_resids, existing_resnames):
                    if pos >= len(subsequence):
                        break
                    # print "DEBUG: Checking if pos=", pos, subsequence[pos], "==", resname, subsequence[pos] == resname
                    if subsequence[pos] != resname:
                        unmatched += 1
                        # matches += 1
                # if matches == len(existing_resids[:existing_resids.index(pos)]):
                if unmatched == 0:
                    # print "DEBUG: Found one possible starting position."
                    possible_starting_positions.append(start)
            except IndexError as e:
                print(e)
                break
        
        if len(possible_starting_positions) == 0:
            print(bcolors.BOLDBLUE + "I could not guess the starting residue number from the labels in the file. There must be a conflict in the residue \
    names you placed and the protein sequence, therefore I will relabel everything starting from 1. If you want to keep your labels please check again \
    for consistency with the protein sequence." + bcolors.ENDBOLD)
            return None
        elif len(possible_starting_positions) == 1:
            print(bcolors.BOLDBLUE + "I found one possible starting position in the sequence ("+str(possible_starting_positions[0])+") from", end=' ')
            print("the labels in the file. Residue numbering starts from "+str(first_resid - possible_starting_positions[0])+". Keeping the", end=' ')
            print("current labels and relabeling all other peaks starting from 'X1N-H', 'X2N-H', etc." + bcolors.ENDBOLD)
            return first_resid - possible_starting_positions[0]
        elif len(possible_starting_positions) > 1:
            print(bcolors.BOLDBLUE + "ERROR: I found multiple possible starting positions in the sequence", end=' ')
            print("("+",".join([str(p) for p in possible_starting_positions])+"). I will keep the first, which covers more sequence.", end=' ')
            print("Residue numbering starts from "+str(first_resid - possible_starting_positions[0])+". Keeping the", end=' ')
            print("current labels and relabeling all other peaks starting from 'X1N-H', 'X2N-H', etc." + bcolors.ENDBOLD)
            return first_resid - possible_starting_positions[0]

    def rewrite_HSQC_spectrum(self):
        """
        Adds comments to lines with overlapping AAIGs. The output file with comments is named HSQC_FILE+".comments".
        :return:
        """
        fout = open(self.HSQC_FILE+".comments", 'w')
        for AAIG in self.__get_all_AAIGs__():
            comment = ""
            if len(AAIG.overlapping_AAIGs) > 0:
                comment = "# " + ",".join([dash_to_NH(a) for a in AAIG.overlapping_AAIGs])
            peak = AAIG.__get_all_peaks__()[0]  # HSQC has only one peak per AAIG
            fout.write("%s\t%.3f\t%.3f\t%s\n" % (dash_to_NH(AAIG.signature), peak.j_HNreson, peak.j_Nreson, comment))