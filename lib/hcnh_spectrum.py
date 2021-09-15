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


from .spectrum import *
from .peak import *

class HCNH_Spectrum(Spectrum):

    def __init__(self,
                 HCNH_FILE=None,
                 HCNH_spectrum_csv=None,
                 Alignment=None,
                 H_bin_length=0.01,
                 C_bin_length=0.1,
                 bandwidth=0.0025):
        """
        If you provide an Alignment object then the constructor will load the assigned NOESY to a DataFrame.
        If not, then the constructor loads the NOESY file into AAIGs that can be used to identify connectivities
        (TOCSY-HCNH or HCNH-HCNH).


        :param HCNH_FILE:
        :param Alignment:   Alignment object
        :param H_bin_length:
        :param C_bin_length:
        :param bandwidth:
        """

        super(HCNH_Spectrum, self).__init__(AAIGs_dict={}, spectrum_type="HCNH")

        # i_residue -> gives the Creson, Hreson
        # j_residue -> gives the Nreson, HNreson
        colnames = ['NresonIndex', 'i_residue', 'j_residue', 'Cname', 'Creson', 'Hname', 'Hreson', 'Nreson',
                    'HNreson', 'aveNreson', 'aveHNreson', 'stdevNreson', 'stdevHNreson', 'Intensity',
                    'scaled_Intensity']
        self.spectrum_df = pd.DataFrame([], columns=colnames)   # appart from AAIGs and peaks, store the whole spectrum information
                                                                # into a DataFrame
        if HCNH_FILE and Alignment:
            self.load_HCNH_to_DataFrame(query_fname=HCNH_FILE, Alignment=Alignment)
        elif HCNH_FILE:
            self.load_HCNH_to_AAIGs(HCNH_FILE, H_bin_length=H_bin_length,
                                    C_bin_length=C_bin_length,
                                    bandwidth=bandwidth)
        elif HCNH_spectrum_csv: # only for nmr_pipeline
            self.load_HCNH_csv_to_AAIGs(HCNH_spectrum_csv=HCNH_spectrum_csv,
                                        H_bin_length=H_bin_length,
                                        C_bin_length=C_bin_length,
                                        bandwidth=bandwidth)


    def load_HCNH_to_AAIGs(self, fname, H_bin_length=0.01, C_bin_length=0.1, bandwidth=0.0025):
        """
        Loads each line of the NOESY sparky list file to Peak objects and subsequently - if N-H labels are present- to AAIGs.

        :param fname:
        :param H_bin_length:
        :param C_bin_length:
        :param bandwidth:
        :return:
        """
        print(bcolors.BOLDBLUE + "Loading NOESY spectrum lines with assigned N-H to AAIG objects." + bcolors.ENDBOLD)
        NOESY_lines, patched_residues_list = self.uniquify_spectrum_lines(query_fname=fname)

        # continue by processing the file lines but use the AAIG of the N,HN peaks (i). The i-1 will be derived from the matching
        # NOESY AAIGs ...
        for line in NOESY_lines:
            words = line.split()
            # print("DEBUG: words=", words)
            j_AAIG_signature = "%s%s" % (words[0].split('-')[-2], words[0].split('-')[-1])  # e.g. 'X127NXHX'
            j_AAIG_name = remove_NH_suffix(j_AAIG_signature)
            if j_AAIG_name == '':  # it is WRONG to have an AAIG without name, because it will contain all unassigned peaks "?-?-?-?"!
                continue
            label = words[0].split('-')[0]
            # m = re.search("^([A-Z][0-9]+)[CMQH].*", label)    # old but worked
            m = re.search('([A-Za-z0-9][0-9]+)[HQM][ABGDEHZ][1-3]{0,2}.*', label) # more complete (I think)
            if m:
                i_AAIG_name = m.group(1)  # iminus1 because it is NOESY
            elif label == "?":
                i_AAIG_name = "?"
            else:
                raise NameError(bcolors.FAIL + "ERROR: wrong label in NOESY line:\n" + line + bcolors.ENDC)
            Hreson = round(float(words[1]), 3)
            Creson = round(float(words[2]), 3)
            Nreson = round(float(words[3]), 3)
            HNreson = round(float(words[4]), 3)
            try:
                intensity = float(words[5])
            except IndexError:
                intensity = None
            # Now save the labes in the input file
            components = remove_NH_suffix(words[0]).split('-')  # this will give something like ['A297HA', 'CA', 'A298']
            f_Cname = components[1]
            f_Hname = "?"
            m = re.search('([A-Za-z0-9][0-9]+)([HQM][ABGDEHZ][1-3]{0,2})', components[0])
            if m:
                f_Hname = m.group(2)
            # NOTE: we cannot know iminus1_AAIG_signature from TOCSY because the N-H belongs to resid i
            peak = Peak(i_AAIG_name=i_AAIG_name,
                        j_AAIG_signature=j_AAIG_signature,
                        j_AAIG_name=j_AAIG_name,
                        file_label=words[0],
                        Creson=Creson,
                        f_Cname=f_Cname,
                        Hreson=Hreson,
                        f_Hname=f_Hname,
                        j_Nreson=Nreson,
                        j_HNreson=HNreson,
                        Intensity=intensity,
                        H_bin_length=H_bin_length,
                        C_bin_length=C_bin_length,
                        bandwidth=bandwidth)
            if not self.__exists__(j_AAIG_signature):
                # by convention name the AAIGs by their NH peak label
                aaig = AAIG(signature=j_AAIG_signature,
                            NH_name=j_AAIG_name,
                            user_assignment=dash_to_NH(j_AAIG_signature),
                            H_bin_length=H_bin_length,
                            C_bin_length=C_bin_length,
                            bandwidth=bandwidth)
                self.__add_AAIG__(signature=aaig.signature, AAIG=aaig)
            self.__add_peak2AAIG__(peak=peak, AAIG_signature=j_AAIG_signature)


    def load_HCNH_to_DataFrame(self,
                               query_fname,
                               Alignment):
        """
        Method to read the NOESY file with the assignments (just the N-H [e.g. *num.list] or all
        [as produced by cs_assignment.py script]) in SPARKY format. Unassigned "?-?-?-?" lines will
        also be loaded to the DataFrame.

        We don't use the notation 'AAIG' here, but 'residue' because we know the assignments.
        Also, we don't use 'i' and 'iplus1' because in NOESY it can be any residue. Instead of 'iplus1' we use 'j'.
        i_residue -> gives the Creson, Hreson
        j_residue -> gives the Nreson, HNreson

        :param query_fname:
        :param absolute_AAIGmatches_alignment_list:
        :param absolute_matches_alignment_list:
        :return:
        """

        spectrum_type = "HCNH"
        print("Loading " + spectrum_type + " file " + query_fname)
        NOESY_lines, patched_residues_list = self.uniquify_spectrum_lines(query_fname=query_fname)

        # populate a dictionary with keys the resids and values the list of the respective C-H & N-H resonances from the
        # spectrum file
        query_lineLists_list = []  # list of the lines of query frame in list form not in string
        for qline in NOESY_lines:
            qline_list = qline.split()
            query_lineLists_list.append(qline_list)
        sorted_query_lineLists_list = sorted(query_lineLists_list,
                                             key=itemgetter(0))  # sort the spectrum lines by the assignment

        for qline_list in sorted_query_lineLists_list:
            print("DEBUG: ------------------------->")
            print("DEBUG: reading line:", qline_list)
            if qline_list[0] == '?-?-?-?':  # if this is an non-labeled peak (?-?-?-?), skipt it
                print("WARNING: the following line will be omitted because it is not assigned:", qline_list)
                continue
            elif qline_list[0][-3:] != 'N-H':   # only N-H mapped must pass, not NX-HX or sidechains
                print("WARNING: the following line will be omitted because the amide is not of the backbone's:", qline_list)
                continue
            components = qline_list[0].replace('N-H', '').split('-')  # this will give something like ['A297HA', 'CA', 'A298'] (only N-H mapped have passed)
            # print "DEBUG: j_residue=", j_residue
            i_Cname = components[1]
            if components[0] == '?' and components[1] == '?':
                i_Hname = "?"
                i_residue = "?"
            else:  # only if this is a user-provided annotated NOESY file, it will contain the i residue label (in addition to the i-1)
                mo = re.search('([A-Za-z0-9][0-9]+)([HQM][ABGDEHZ][1-3]{0,2})',
                               components[0])  # read only lines with atom type assignments
                if mo:
                    i_residue = mo.group(1)
                    i_Hname = mo.group(2)
                else:
                    print("WARNING: the following line will be ommited because Hname and TAAIG could not be found:", \
                    components[0])
                    # sys.exit(1)
                    continue

            # Rename VAL QG1-->HG1, QG2-->HG2 and LEU QD1-->HD1, QD2-->HD2 for compatibility with the code and the histogram file names
            if i_residue[0] == 'L' and i_Hname == "QD1":
                i_Hname = "HD1"
            elif i_residue[0] == 'L' and i_Hname == "QD2":
                i_Hname = "HD2"
            if i_residue[0] == 'V' and i_Hname == "QG1":
                i_Hname = "HG1"
            elif i_residue[0] == 'V' and i_Hname == "QG2":
                i_Hname = "HG2"

            if components[2] == '':
                j_residue = i_residue
            elif Alignment == 'ignore_alignment':   # special case for NOESY featvec creation
                j_residue = components[2]
            else:  # if this is a NOESY *num.list file or user-provided annotated NOESY file, it will contain the i residue label
                # (apart from the i-1)
                j_residue = get_residue_from_AAIGsignature(components[2], Alignment.absolute_AAIGmatches_alignment_list,
                                                           Alignment.absolute_matches_alignment_list)  # rename the RIG to the residue actual name for consistency with TOCSY
            print("DEBUG: i_residue=", i_residue, "j_residue=", j_residue, "spectrum_type=", spectrum_type, \
                "components=", components)

            row_dict = {}  # stores values from this line for the self.spectrum_df
            row_dict['j_residue'] = j_residue
            row_dict['NresonIndex'] = self.spectrum_df[self.spectrum_df.j_residue == j_residue].shape[0]
            row_dict['i_residue'], row_dict['Cname'], row_dict['Creson'], row_dict['Hname'], row_dict['Hreson'] = \
                rename_lone_methylene_proton([i_residue, i_Cname, float(qline_list[2]), i_Hname, float(qline_list[1])],
                                             aatype_carbon_degenerateH_mdict,
                                             aatype_carbon_nondegenerateHlist_mdict)
            row_dict['Intensity'] = float(qline_list[5])
            row_dict['Nreson'], row_dict['HNreson'] = float(qline_list[3]), float(qline_list[4])
            self.spectrum_df = self.spectrum_df.append(row_dict, ignore_index=True)  # save this line to the dataframe

        # Finally average and save the N and HN resonances for each residue, and normalize Intensities
        for j_residue in set(self.spectrum_df['j_residue']):
            print("DEBUG: averaging N,HN and normalizing intensity of j_residue %s" % j_residue)
            Nresons = self.spectrum_df.loc[(self.spectrum_df['j_residue'] == j_residue)].get('Nreson').dropna(
                axis=0, how='any')
            aveNreson = Nresons.mean()
            stdevNreson = Nresons.std()
            HNresons = self.spectrum_df.loc[(self.spectrum_df['j_residue'] == j_residue)].get('HNreson').dropna(
                axis=0, how='any')
            aveHNreson = HNresons.mean()
            stdevHNreson = HNresons.std()
            self.spectrum_df.aveNreson[self.spectrum_df.j_residue == j_residue] = aveNreson
            self.spectrum_df.stdevNreson[self.spectrum_df.j_residue == j_residue] = stdevNreson
            self.spectrum_df.aveHNreson[self.spectrum_df.j_residue == j_residue] = aveHNreson
            self.spectrum_df.stdevHNreson[self.spectrum_df.j_residue == j_residue] = stdevHNreson
            # Normalize Intensities (within AAIG)
            max_Intensity = float(self.spectrum_df.Intensity[self.spectrum_df.j_residue == j_residue].max())
            # TODO: use indices to speed up the following.
            for Intensity in self.spectrum_df.Intensity[self.spectrum_df.j_residue == j_residue]:
                self.spectrum_df.scaled_Intensity[(self.spectrum_df.j_residue == j_residue) &
                                                (self.spectrum_df.Intensity == Intensity)] = Intensity / max_Intensity

    def load_HCNH_csv_to_DataFrame(self, HCNH_frames_csv):
        """
        UNTESTED (needs scaled_intensities)
        ** SPECIAL METHOD FOR NMR_PIPELINE. **

        Method to read the NOESY file with the assignments (just the N-H [e.g. *num.list] or all
        [as produced by cs_assignment.py script]) in SPARKY format. Unassigned "?-?-?-?" lines will
        also be loaded to the DataFrame.

        We don't use the notation 'AAIG' here, but 'residue' because we know the assignments.
        Also, we don't use 'i' and 'iplus1' because in NOESY it can be any residue. Instead of 'iplus1' we use 'j'.
        i_residue -> gives the Creson, Hreson
        j_residue -> gives the Nreson, HNreson

        :param HCNH_frames_csv:
        :return:
        """

        # TODO: create Alignment() object and this function is ready!

        colnames = ['NresonIndex', 'i_residue', 'j_residue', 'Cname', 'Creson', 'Hname', 'Hreson', 'Nreson',
                    'HNreson', 'aveNreson', 'aveHNreson', 'stdevNreson', 'stdevHNreson', 'Intensity',
                    'scaled_Intensity']
        self_spectrum_df = pd.DataFrame([], columns=colnames)   # appart from AAIGs and peaks, store the whole spectrum information
                                                                # into a DataFrame

        HCNH_spectrum_df = pd.read_csv(HCNH_frames_csv)
        # sort the spectrum lines by the assignment
        HCNH_spectrum_df.sort_values(by="assignment", inplace=True)

        # Adds extra column with scaled peak intensity with respect to the maximum intensity of the
        # whole spectrum.
        HCNH_peaks_df["scaled_intensity"] = HCNH_peaks_df["intensity"]/HCNH_peaks_df["intensity"].max()

        for idx, row in HCNH_spectrum_df.iterrows():
            print("DEBUG: ------------------------->")
            print("DEBUG: reading row:", row)
            if row["assignment"] == '?-?-?-?':  # if this is an non-labeled peak (?-?-?-?), skipt it
                print("WARNING: the following line will be omitted because it is not assigned:", row)
                continue
            elif row["assignment"][-3:] != 'N-H':   # only N-H mapped must pass, not NX-HX or sidechains
                print("WARNING: the following line will be omitted because the amide is not of the backbone's:", row)
                continue
            components = row["assignment"].replace('N-H', '').split('-')  # this will give something like ['A297HA', 'CA', 'A298'] (only N-H mapped have passed)
            # print "DEBUG: j_residue=", j_residue
            i_Cname = components[1]
            if components[0] == '?' and components[1] == '?':
                i_Hname = "?"
                i_residue = "?"
            else:  # only if this is a user-provided annotated NOESY file, it will contain the i residue label (in addition to the i-1)
                mo = re.search('([A-Za-z0-9][0-9]+)([HQM][ABGDEHZ][1-3]{0,2})',
                               components[0])  # read only lines with atom type assignments
                if mo:
                    i_residue = mo.group(1)
                    i_Hname = mo.group(2)
                else:
                    print("WARNING: the following line will be ommited because Hname and TAAIG could not be found:", \
                    components[0])
                    # sys.exit(1)
                    continue

            # Rename VAL QG1-->HG1, QG2-->HG2 and LEU QD1-->HD1, QD2-->HD2 for compatibility with the code and the histogram file names
            if i_residue[0] == 'L' and i_Hname == "QD1":
                i_Hname = "HD1"
            elif i_residue[0] == 'L' and i_Hname == "QD2":
                i_Hname = "HD2"
            if i_residue[0] == 'V' and i_Hname == "QG1":
                i_Hname = "HG1"
            elif i_residue[0] == 'V' and i_Hname == "QG2":
                i_Hname = "HG2"

            if components[2] == '':
                j_residue = i_residue
            elif Alignment == 'ignore_alignment':   # special case for NOESY featvec creation
                j_residue = components[2]
            else:  # if this is a NOESY *num.list file or user-provided annotated NOESY file, it will contain the i residue label
                # (apart from the i-1)
                j_residue = get_residue_from_AAIGsignature(components[2], Alignment.absolute_AAIGmatches_alignment_list,
                                                           Alignment.absolute_matches_alignment_list)  # rename the RIG to the residue actual name for consistency with TOCSY
            print("DEBUG: i_residue=", i_residue, "j_residue=", j_residue, "spectrum_type= HCNH NOESY", \
                "components=", components)

            row_dict = {}  # stores values from this line for the self_spectrum_df
            row_dict['j_residue'] = j_residue
            row_dict['NresonIndex'] = self_spectrum_df[self_spectrum_df.j_residue == j_residue].shape[0]
            row_dict['i_residue'], row_dict['Cname'], row_dict['Creson'], row_dict['Hname'], row_dict['Hreson'] = \
                rename_lone_methylene_proton([i_residue, i_Cname, row["C"], i_Hname, row["HC"]],
                                             aatype_carbon_degenerateH_mdict,
                                             aatype_carbon_nondegenerateHlist_mdict)

            row_dict['Intensity'], row_dict['Nreson'], row_dict['HNreson'] = \
                row["intensity"], row["N"], row["HN"]
            self_spectrum_df = self_spectrum_df.append(row_dict, ignore_index=True)  # save this line to the dataframe

        # Finally average and save the N and HN resonances for each residue, and normalize Intensities
        for j_residue in set(self_spectrum_df['j_residue']):
            print("DEBUG: averaging N,HN and normalizing intensity of j_residue %s" % j_residue)
            Nresons = self_spectrum_df.loc[(self_spectrum_df['j_residue'] == j_residue)].get('Nreson').dropna(
                axis=0, how='any')
            aveNreson = Nresons.mean()
            stdevNreson = Nresons.std()
            HNresons = self_spectrum_df.loc[(self_spectrum_df['j_residue'] == j_residue)].get('HNreson').dropna(
                axis=0, how='any')
            aveHNreson = HNresons.mean()
            stdevHNreson = HNresons.std()
            self_spectrum_df.aveNreson[self_spectrum_df.j_residue == j_residue] = aveNreson
            self_spectrum_df.stdevNreson[self_spectrum_df.j_residue == j_residue] = stdevNreson
            self_spectrum_df.aveHNreson[self_spectrum_df.j_residue == j_residue] = aveHNreson
            self_spectrum_df.stdevHNreson[self_spectrum_df.j_residue == j_residue] = stdevHNreson
            # Normalize Intensities (within AAIG)
            max_Intensity = float(self_spectrum_df.Intensity[self_spectrum_df.j_residue == j_residue].max())
            # TODO: use indices to speed up the following.
            for Intensity in self_spectrum_df.Intensity[self_spectrum_df.j_residue == j_residue]:
                self_spectrum_df.scaled_Intensity[(self_spectrum_df.j_residue == j_residue) &
                                                (self_spectrum_df.Intensity == Intensity)] = Intensity / max_Intensity

        return self_spectrum_df


    def load_HCNH_csv_to_AAIGs(self, HCNH_spectrum_csv, H_bin_length=0.01, C_bin_length=0.1, bandwidth=0.0025):
        """
        ** SPECIAL METHOD FOR NMR PIPELINE **

        Loads each line of the NOESY sparky list file to Peak objects and subsequently - if N-H labels are
        present- to AAIGs.

        :param fname:
        :param H_bin_length:
        :param C_bin_length:
        :param bandwidth:
        :return:
        """
        ColorPrint("Loading NOESY spectrum lines with assigned N-H to AAIG objects, "
                   "from file %s." % HCNH_spectrum_csv, "BOLDBLUE")
        HCNH_spectrum_df = pd.read_csv(HCNH_spectrum_csv)
        HCNH_spectrum_df.fillna("", inplace=True)
        # sort the spectrum lines by the assignment
        HCNH_spectrum_df.sort_values(by="assignment", inplace=True)

        # continue by processing the file lines but use the AAIG of the N,HN peaks (i). The i-1 will be
        # derived from the matching NOESY AAIGs ...
        for idx, row in HCNH_spectrum_df.iterrows():
            # j_AAIG_signature = "%s%s" % (row["assignment"].split('-')[-2], row["assignment"].split('-')[-1])  # e.g. 'X127NXHX'
            j_AAIG_signature = row["frame"]
            print("DEBUG: j_AAIG_signature=", j_AAIG_signature)
            j_AAIG_name = remove_NH_suffix(j_AAIG_signature)
            if j_AAIG_name == '':  # it is WRONG to have an AAIG without name, because it will contain all unassigned peaks "?-?-?-?"!
                continue
            label = row["assignment"].split('-')[0]
            # m = re.search("^([A-Z][0-9]+)[CMQH].*", label)    # old but worked
            m = re.search('([A-Za-z0-9][0-9]+)[HQM][ABGDEHZ][1-3]{0,2}.*', label) # more complete (I think)
            if m:
                i_AAIG_name = m.group(1)  # iminus1 because it is NOESY
            elif label == "?":
                i_AAIG_name = "?"
            else:
                raise NameError(bcolors.FAIL + "ERROR: wrong label in NOESY line:\n" + line + bcolors.ENDC)
            Hreson = round(row["HC"], 3)
            Creson = round(row["C"], 3)
            Nreson = round(row["N"], 3)
            HNreson = round(row["HN"], 3)
            try:
                intensity = row["intensity"]
            except IndexError:
                intensity = None
            # Now save the labes in the input file
            components = remove_NH_suffix(row["assignment"]).split('-')  # this will give something like ['A297HA', 'CA', 'A298']
            f_Cname = components[1]
            f_Hname = "?"
            m = re.search('([A-Za-z0-9][0-9]+)([HQM][ABGDEHZ][1-3]{0,2})', components[0])
            if m:
                f_Hname = m.group(2)
            # NOTE: we cannot know iminus1_AAIG_signature from TOCSY because the N-H belongs to resid i
            peak = Peak(i_AAIG_name=i_AAIG_name,
                        j_AAIG_signature=j_AAIG_signature,
                        j_AAIG_name=j_AAIG_name,
                        file_label=row["assignment"],
                        Creson=Creson,
                        f_Cname=f_Cname,
                        Hreson=Hreson,
                        f_Hname=f_Hname,
                        j_Nreson=Nreson,
                        j_HNreson=HNreson,
                        Intensity=intensity,
                        H_bin_length=H_bin_length,
                        C_bin_length=C_bin_length,
                        bandwidth=bandwidth)
            if not self.__exists__(j_AAIG_signature):
                # by convention name the AAIGs by their NH peak label
                aaig = AAIG(signature=j_AAIG_signature,
                            NH_name=j_AAIG_name,
                            user_assignment=dash_to_NH(j_AAIG_signature),
                            H_bin_length=H_bin_length,
                            C_bin_length=C_bin_length,
                            bandwidth=bandwidth)
                self.__add_AAIG__(signature=aaig.signature, AAIG=aaig)
            self.__add_peak2AAIG__(peak=peak, AAIG_signature=j_AAIG_signature)


################################################ PRINT FUNCTIONS #############################################

    def __print_contents__(self):

        for aaig in self.__get_all_AAIGs__():
            print("\nContents of AAIG %s:" % aaig.signature)
            aaig.__print_peaks__(only_CHresons=True)