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

class TOCSY_Spectrum(Spectrum):

    def __init__(self,
                 TOCSY_fname,
                 Alignment=None,
                 H_bin_length=0.01,
                 C_bin_length=0.1,
                 bandwidth=0.0025):
        """
        If you provide an Alignment object then the constructor will load the assigned NOESY to a DataFrame.
        If not, then the constructor loads the NOESY file into AAIGs that can be used to identify connectivities
        (TOCSY-HCNH or HCNH-HCNH).

        :param TOCSY_fname:
        :param Alignment:
        :param H_bin_length:
        :param C_bin_length:
        :param bandwidth:
        """

        super(TOCSY_Spectrum, self).__init__(AAIGs_dict={}, spectrum_type="")

        colnames = ['i_residue', 'iplus1_residue', 'Cname', 'Creson', 'Hname', 'Hreson', 'Nreson', 'HNreson',
                    'aveNreson', 'aveHNreson', 'stdevNreson', 'stdevHNreson', 'Intensity', 'scaled_Intensity']
        self.spectrum_df = pd.DataFrame([], columns=colnames)   # appart from AAIGs and peaks, store the whole spectrum information
                                                                # into a DataFrame

        if Alignment:
            self.load_TOCSY_to_DataFrame(query_fname=TOCSY_fname,
                                         Alignment=Alignment)
        else:
            self.load_TOCSY_to_AAIGs(TOCSY_fname,
                                     H_bin_length=H_bin_length,
                                     C_bin_length=C_bin_length,
                                     bandwidth=bandwidth)

    def load_TOCSY_to_AAIGs(self, fname, H_bin_length=0.01, C_bin_length=0.1, bandwidth=0.0025):
        """
        Loads each line of the TOCSY sparky list file to Peak objects and subsequently - if N-H labels are present- to AAIGs.

        :param fname:
        :param H_bin_length:
        :param C_bin_length:
        :param bandwidth:
        :return:
        """
        ColorPrint("Loading TOCSY spectrum lines with assigned N-H to AAIG objects.", "BOLDBLUE")
        TOCSY_lines, patched_residues_list = self.uniquify_spectrum_lines(query_fname=fname)

        # continue by processing the file lines but use the AAIG of the N,HN peaks (i). The i-1 will be derived from the matching
        # NOESY AAIGs ...
        for line in TOCSY_lines:
            words = line.split()
            # print "DEBUG: words=", words
            i_AAIG_signature = "%s%s" % (words[0].split('-')[-2], words[0].split('-')[-1])  # e.g. 'X127NX-HX'
            i_AAIG_name = remove_NH_suffix(i_AAIG_signature)
            if i_AAIG_name == '':  # it is WRONG to have an AAIG without name, because it will contain all unassigned peaks "?-?-?-?"!
                continue
            label = words[0].split('-')[0]
            # m = re.search("^([A-Z][0-9]+)[CMQH].*", label)    # old but worked
            m = re.search('([A-Za-z0-9][0-9]+)[HQM][ABGDEHZ][1-3]{0,2}.*', label) # more complete (I think)
            if m:
                iminus1_AAIG_name = m.group(1)  # iminus1 because it is TOCSY
            elif label == "?":
                iminus1_AAIG_name = "?"
            else:
                raise NameError(bcolors.FAIL + "ERROR: wrong label in TOCSY line:\n" + line + bcolors.ENDC)
            Hreson = round(float(words[1]), 3)
            Creson = round(float(words[2]), 3)
            Nreson = round(float(words[3]), 3)
            HNreson = round(float(words[4]), 3)
            try:
                intensity = float(words[5])
            except IndexError:
                intensity = None
            # NOTE: we cannot know iminus1_AAIG_signature from TOCSY because the N-H belongs to resid i
            peak = Peak(i_AAIG_name=iminus1_AAIG_name, j_AAIG_signature=i_AAIG_signature, j_AAIG_name=i_AAIG_name,
                        file_label=words[0], Creson=Creson, Hreson=Hreson, j_Nreson=Nreson, j_HNreson=HNreson,
                        Intensity=intensity, H_bin_length=H_bin_length, C_bin_length=C_bin_length, bandwidth=bandwidth)
            if not self.__exists__(i_AAIG_signature):
                # by convention name the AAIGs by their NH peak label
                aaig = AAIG(signature=i_AAIG_signature,
                            CH_name=iminus1_AAIG_name,
                            user_assignment=dash_to_NH(i_AAIG_signature),
                            H_bin_length=H_bin_length,
                            C_bin_length=C_bin_length,
                            bandwidth=bandwidth)
                self.__add_AAIG__(signature=aaig.signature, AAIG=aaig)
            self.__add_peak2AAIG__(peak=peak, AAIG_signature=i_AAIG_signature)  # name this AAIG by the N-H resonance source

    def load_TOCSY_to_DataFrame(self,
                                query_fname,
                                Alignment):
        """
        Method to read the TOCSY file with the assignments (just the N-H [e.g. *num.list] or all [as produced by cs_assignment.py script])
        in SPARKY format.

        :param query_fname:
        :return:  i_residue & iplus1_residue in this function are RIG not valid residue names!
        :return:  residue_assignments_dict:    i_residue ->   [ (iminus1_residue, Cname, Cresonance, Hname, Hresonance),
                                                                (iminus1_residue, Cname, Cresonance, Hname, Hresonance), ... ]
        :return:  residue_NHresonances_dict:   i_residue -> (average_Nresonance, average_Hresonance, stdev_N_resonance, stdev_H_resonance)
        """

        spectrum_type = "TOCSY"
        print("Loading " + spectrum_type + " file " + query_fname)
        TOCSY_lines, patched_residues_list = self.uniquify_spectrum_lines(query_fname=query_fname)

        query_lineLists_list = []  # list of the lines of query frame in list form not in string
        for qline in TOCSY_lines:
            qline_list = qline.split()
            query_lineLists_list.append(qline_list)
        sorted_query_lineLists_list = sorted(query_lineLists_list, key=itemgetter(0))  # sort the spectrum lines by the assignment

        i_to_iminus1_dict = {}  # dictionary to be populated from the assigned TOCSY file
        for qline_list in sorted_query_lineLists_list:
            print("DEBUG: ------------------------->")
            print("DEBUG: reading line:", qline_list)
            if qline_list[0] == '?-?-?-?':  # if this is an non-labeled peak (?-?-?-?), skipt it
                print("WARNING: the following line will be omitted because it is not assigned:", qline_list)
                continue
            elif not qline_list[0][-3:] in ['N-H', 'NX-HX']:   # only N-H mapped must pass, not NX-HX or sidechains
                print("WARNING: the following line will be omitted because the amide is not of the backbone's:", qline_list)
                continue
            components = remove_NH_suffix(qline_list[0]).split('-')  # this will give something like ['A297HA', 'CA', 'A298'] (only N-H mapped have passed)
            # print "DEBUG: iplus1_residue=", iplus1_residue
            i_Cname = components[1]
            mo = re.search('([A-Za-z0-9][0-9]+)([HQM][ABGDEHZ][1-3]{0,2})', components[0])  # read only lines with atom type assignments
            if mo:
                i_residue = mo.group(1)
                i_Hname = mo.group(2)
            else:
                raise  Exception(bcolors.FAIL + "ERROR: cannot find Hname and AAIG in string: " + " ".join(components[
                    0]) + ". \nPlease correct the following line in TOCSY:\n" + qline_list + bcolors.ENDC)

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
                iplus1_residue = i_residue
            else:   # if this is a TOCSY *num.list file or user-provided annotated TOCSY file, it will contain the i residue label
                    # (apart from the i-1)
                iplus1_residue = get_residue_from_AAIGsignature(components[2],
                                                                Alignment.absolute_AAIGmatches_alignment_list,
                                                                Alignment.absolute_matches_alignment_list)  # rename the RIG to the residue actual name for consistency with TOCSY
            print("DEBUG: i_residue=", i_residue, "iplus1_residue=", iplus1_residue, \
                "spectrum_type=", spectrum_type, "components=", components)
            # print "DEBUG: point 0, iplus1_residue=", iplus1_residue

            row_dict = {}   # stores values from this line for the self.spectrum_df
            row_dict['iplus1_residue'] = iplus1_residue
            row_dict['i_residue'], row_dict['Cname'], row_dict['Creson'], row_dict['Hname'], row_dict['Hreson'] = \
                rename_lone_methylene_proton([i_residue, i_Cname, float(qline_list[2]), i_Hname, float(qline_list[1])],
                aatype_carbon_degenerateH_mdict,
                aatype_carbon_nondegenerateHlist_mdict)
            if len(qline_list) == 6:    # if the TOCSY has intensities
                row_dict['Intensity'] = float(qline_list[5])
            i_to_iminus1_dict[iplus1_residue] = i_residue
            row_dict['Nreson'], row_dict['HNreson'] = float(qline_list[3]), float(qline_list[4])
            self.spectrum_df = self.spectrum_df.append(row_dict, ignore_index=True)   # save this line to the dataframe

        # Finaly average and save the N and HN resonances for each residue, and normalize Intensities
        for iplus1_residue in set(self.spectrum_df['iplus1_residue']):
            Nresons = self.spectrum_df.loc[(self.spectrum_df['iplus1_residue'] == iplus1_residue)].get(
                'Nreson').dropna(axis=0, how='any')
            aveNreson = Nresons.mean()
            stdevNreson = Nresons.std()
            HNresons = self.spectrum_df.loc[(self.spectrum_df['iplus1_residue'] == iplus1_residue)].get(
                'HNreson').dropna(axis=0, how='any')
            aveHNreson = HNresons.mean()
            stdevHNreson = HNresons.std()
            self.spectrum_df.aveNreson[self.spectrum_df.iplus1_residue == iplus1_residue] = aveNreson
            self.spectrum_df.stdevNreson[self.spectrum_df.iplus1_residue == iplus1_residue] = stdevNreson
            self.spectrum_df.aveHNreson[self.spectrum_df.iplus1_residue == iplus1_residue] = aveHNreson
            self.spectrum_df.stdevHNreson[self.spectrum_df.iplus1_residue == iplus1_residue] = stdevHNreson
            # Normalize Intensities
            # TODO: use indices to speed up the following.
            if not self.spectrum_df.Intensity.isnull().all():
                max_Intensity = float(
                    self.spectrum_df.Intensity[self.spectrum_df.iplus1_residue == iplus1_residue].max())
                for Intensity in self.spectrum_df.Intensity[self.spectrum_df.iplus1_residue == iplus1_residue]:
                    self.spectrum_df.scaled_Intensity[(self.spectrum_df.iplus1_residue == iplus1_residue) &
                                                      (
                                                                  self.spectrum_df.Intensity == Intensity)] = Intensity / max_Intensity

        return i_to_iminus1_dict, patched_residues_list