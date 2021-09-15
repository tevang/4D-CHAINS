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


import pickle
from collections import defaultdict

from .spectrum import *


class HNNH_Spectrum(Spectrum):

    def __init__(self, HNNH_FILE, H_bin_length=0.01, C_bin_length=0.1, bandwidth=0.0025, uniquify_lines=True):
        super(HNNH_Spectrum, self).__init__(AAIGs_dict={}, spectrum_type="HNNH")

        # i_AAIG -> gives the i_Nreson, i_HNreson
        # j_AAIG -> gives the j_Nreson, j_HNreson
        colnames = ['i_AAIG', 'j_AAIG', 'i_Nreson', 'i_HNreson', 'j_Nreson', 'j_HNreson',
                    'i_aveNreson', 'i_aveHNreson', 'i_stdevNreson', 'i_stdevHNreson',
                    'j_aveNreson', 'j_aveHNreson', 'j_stdevNreson', 'j_stdevHNreson',
                    'Intensity', 'scaled_Intensity', 'distance', 'prob_ASN', 'prob_GLN', 'prob_UNK',
                    'i_HSQCdist', 'j_HSQCdist', 'test_assignment', 'i_agreement', 'j_agreement']
        # NOTE: by defining the dtype you will be able to retrieve the value of a column by simply doing row[colname]
        self.spectrum_df = pd.DataFrame([],
                                        columns=colnames, dtype='float64')  # appart from AAIGs and peaks, store the whole spectrum information
        self.bb_Wsc_spectrum_df = pd.DataFrame([], columns=colnames, dtype='float64')    # only backbone N,HN pairs
        self.NQsc_spectrum_df = pd.DataFrame([], columns=colnames, dtype='float64')    # only sidechain N,HN pairs
        # Set dtype of non-numeric columns explicitly
        self.spectrum_df = self.spectrum_df.astype({'i_AAIG':'str', 'j_AAIG':'str', 'test_assignment':'str',
                                                                    'i_agreement':'bool', 'j_agreement':'bool'})    # they will appear as 'object' not 'str'
        self.bb_Wsc_spectrum_df = self.bb_Wsc_spectrum_df.astype({'i_AAIG': 'str', 'j_AAIG': 'str', 'test_assignment':'str',
                                                                    'i_agreement':'bool', 'j_agreement':'bool'})    # they will appear as 'object' not 'str'
        self.NQsc_spectrum_df = self.NQsc_spectrum_df.astype({'i_AAIG': 'str', 'j_AAIG': 'str', 'test_assignment':'str',
                                                                    'i_agreement':'bool', 'j_agreement':'bool'})    # they will appear as 'object' not 'str'


        # into a DataFrame
        self.HNNH_FILE = HNNH_FILE
        self.load_HNNH_to_DataFrame(query_fname=HNNH_FILE, uniquify_lines=uniquify_lines)

        # TODO: see in the future if you need to compare AAIG histograms
        # else:
        #     self.load_HNNH_to_AAIGs(HNNH_FILE,
        #                              H_bin_length=H_bin_length,
        #                              C_bin_length=C_bin_length,
        #                              bandwidth=bandwidth)

    def __scale_intensities(self):
        """
        OVERRIDES
        Method used only by HNNH-connectivities. In that spectrum we don't have AAIGs, therefore a separate function
        was writen that scales all intensities in the DataFrame with respect to the maximum intensity of the HNNH spectrum.
        :return:
        """
        self.spectrum_df['scaled_Intensity'] = self.spectrum_df['Intensity'].div(self.spectrum_df['Intensity'].max())

    def distance_Series(self, row1, row2):
        """
        Special distance function for HNNH to identify the ASN,GLN sidechain N-H.
        row1 and row2 must be in Series format.
        :param row1:
        :param row2:
        :return:
        """
        distance = np.sqrt(np.sum(
            (row1['i_HNreson'] - row2['j_HNreson']) ** 2 + \
            ((row1['i_Nreson'] - row2['i_Nreson']) / 6) ** 2 + \
            (row1['j_Nreson'] - row2['j_Nreson']) ** 2 + \
            ((row1['j_HNreson'] - row2['i_HNreson']) / 6) ** 2
        ))
        return distance

    def distance_DataFrame(self, row1, row2):
        """
        Special distance function for HNNH to identify the ASN,GLN sidechain N-H.
        Each HNNH lines contains the N and HN resonances of i_AAIG and j_AAIG. The distance is
        defined as: (i_HNreson1-j_HNreson2)**2 + ((i_Nreson1-i_Nreson2)**2)/6 +
                    (j_HNreson1-i_HNreson2)**2 + ((j_Nreson1-j_Nreson2)**2)/6
        # TODO: why i_Nreson1-i_Nreson2 and not i_Nreson1-j_Nreson2?

        row1 and row2 must be in DataFrame format.
        :param row1:
        :param row2:
        :return distance:
        """
        row1 = row1.reset_index()
        row2 = row2.reset_index()
        distance = np.sqrt(np.sum(
            (row1.at[0, 'i_HNreson'] - row2.at[0, 'j_HNreson']) ** 2 + \
            ((row1.at[0, 'i_Nreson'] - row2.at[0, 'i_Nreson']) / 6) ** 2 + \
            (row1.at[0, 'j_Nreson'] - row2.at[0, 'j_Nreson']) ** 2 + \
            ((row1.at[0, 'j_HNreson'] - row2.at[0, 'i_HNreson']) / 6) ** 2
        ))
        return distance

    def load_HNNH_to_DataFrame(self, query_fname, uniquify_lines=True):
        """
        Method to read the HNNH file with the assignments (just the N-H [e.g. *num.list] or all
        [as produced by cs_assignment.py script])
        in SPARKY format. Unassigned "?-?-?-?" lines will also be loaded to the DataFrame.

        We don't use the notation 'AAIG' here, but 'residue' because we know the assignments.
        Also, we don't use 'i' and 'iplus1' because in HNNH it can be any residue. Instead of 'iplus1' we use 'j'.
        i_AAIG -> gives the Creson, Hreson
        j_AAIG -> gives the Nreson, HNreson

        :param query_fname: can be a CSV or a Sparky list file (HNNH num.list file)
        :return:
        """

        spectrum_type = "HNNH"
        print("Loading " + spectrum_type + " file " + query_fname)
        HNNH_lines, patched_residues_list = \
            self.uniquify_spectrum_lines(query_fname=query_fname, uniquify_lines=uniquify_lines)

        # populate a dictionary with keys the resids and values the list of the respective C-H & N-H
        # resonances from the spectrum file
        query_lineLists_list = []  # list of the lines of query frame in list form not in string
        for qline in HNNH_lines:
            qline_list = qline.split()
            query_lineLists_list.append(qline_list)
        sorted_query_lineLists_list = sorted(query_lineLists_list,
                                             key=itemgetter(0))  # sort the spectrum lines by the assignment
        residue_assignments_dict = {}  # i_AAIG -> [ (i_AAIG, Cname, Cresonance, Hname, Hresonance), (i_AAIG, Cname, Cresonance, Hname, Hresonance), ... ]
        residue_NHresonances_dict = {}  # j_AAIG -> [(j_AAIG, Nresonance, Hresonance), (j_AAIG, Nresonance, Hresonance), ...]
        original_HNNH_peaks_dict = {}  # residue -> [Hreson, Creson, Nreson, HNreson]
        residue_peak_intensity_mdict = tree()  # peak (e.g. (23.45, 0,823)) -> intensity

        for qline_list in sorted_query_lineLists_list:
            print("DEBUG: ------------------------->")
            print("DEBUG: reading line:", qline_list)
            if qline_list[0] == '?-?-?-?':  # if this is an non-labeled peak (?-?-?-?), skipt it
                print("WARNING: the following line will be omitted because it is not assigned:", qline_list)
                continue
            # elif qline_list[0][-3:] != 'N-H': # OBSOLETE
            #     print "WARNING: the following line will be ommited because the amide is not of the backbone's:", qline_list
            #     continue
            words = qline_list[0].split('-')
            components = "%s-%s" % (words[0],words[1]), "%s-%s" % (words[2],words[3])  # this will give something like
            # ['N140ND2-HD21', 'E139N-H'], and like ['??', 'I108N-H'] for partly assigned lines
            # print "DEBUG: j_AAIG=", j_AAIG
            i_AAIG, j_AAIG = components
            if i_AAIG == "??" or j_AAIG == "??":
                ColorPrint("WARNING: the following line will be omitted because AAIGs could not be found:", "WARNING")
                "\t".join(qline_list)
                # sys.exit(1)
                continue
            if len(remove_NH_suffix(j_AAIG)) == 0 and len(remove_NH_suffix(reverse_NH_suffix(i_AAIG))) > 0:
                # Corrects cases like "Q133HE22-NE2-NE2-HE21" to this "Q133HE22-NE2-Q133NE2-HE21"
                j_AAIG = remove_NH_suffix(reverse_NH_suffix(i_AAIG)) + j_AAIG

            row_dict = {}  # stores values from this line for the self.spectrum_df
            row_dict['j_AAIG'] = j_AAIG
            row_dict['i_AAIG'], row_dict['j_AAIG'], row_dict['i_HNreson'], \
            row_dict['i_Nreson'], row_dict['j_Nreson'], row_dict['j_HNreson'] = \
                i_AAIG, j_AAIG, float(qline_list[1]), float(qline_list[2]), float(qline_list[3]), float(qline_list[4])
            row_dict['Intensity'] = float(qline_list[5])
            if len(qline_list) >= 8:    # if the file contains distance columns, save them too
                row_dict['i_HSQCdist'] = qline_list[6]
                row_dict['j_HSQCdist'] = qline_list[7]
            if len(qline_list) == 11:    # if the file contains user assignment and proof reading agreement
                row_dict['test_assignment'] = qline_list[8]
                row_dict['i_agreement'] = qline_list[9]
                row_dict['j_agreement'] = qline_list[10]
            self.spectrum_df = self.spectrum_df.append(row_dict, ignore_index=True)  # save this line to the dataframe

        # Finaly average and save the N and HN resonances for each residue, and normalize Intensities
        i_AAIGs = set(self.spectrum_df['i_AAIG'])
        j_AAIGs = set(self.spectrum_df['j_AAIG'])
        all_AAIGs = i_AAIGs.union(j_AAIGs)
        for AAIG in all_AAIGs:
            print("DEBUG: averaging N,HN of AAIG %s" % AAIG)
            Nresons = self.spectrum_df.loc[(self.spectrum_df['j_AAIG'] == AAIG)].get('j_Nreson').dropna(axis=0, how='any')
            Nresons = pd.concat([Nresons, self.spectrum_df.loc[(self.spectrum_df['i_AAIG'] == AAIG)].get('i_Nreson').dropna(axis=0, how='any')])
            aveNreson = Nresons.mean()
            stdevNreson = Nresons.std()
            HNresons = self.spectrum_df.loc[(self.spectrum_df['j_AAIG'] == AAIG)].get('j_HNreson').dropna(axis=0, how='any')
            Nresons = pd.concat([Nresons, self.spectrum_df.loc[(self.spectrum_df['i_AAIG'] == AAIG)].get('i_HNreson').dropna(axis=0, how='any')])
            aveHNreson = HNresons.mean()
            stdevHNreson = HNresons.std()
            self.spectrum_df.j_aveNreson[(self.spectrum_df.j_AAIG == AAIG)] = aveNreson
            self.spectrum_df.i_aveNreson[(self.spectrum_df.i_AAIG == AAIG)] = aveNreson
            self.spectrum_df.j_stdevNreson[(self.spectrum_df.j_AAIG == AAIG)] = stdevNreson
            self.spectrum_df.i_stdevNreson[(self.spectrum_df.i_AAIG == AAIG)] = stdevNreson
            self.spectrum_df.j_aveHNreson[(self.spectrum_df.j_AAIG == AAIG)] = aveHNreson
            self.spectrum_df.i_aveHNreson[(self.spectrum_df.i_AAIG == AAIG)] = aveHNreson
            self.spectrum_df.j_stdevHNreson[(self.spectrum_df.j_AAIG == AAIG)] = stdevHNreson
            self.spectrum_df.i_stdevHNreson[(self.spectrum_df.i_AAIG == AAIG)] = stdevHNreson
            # TODO: I am not sure if this makes sense for HNNH.
            # # Normalize Intensities (within AAIG)
            # max_Intensity = float(self.spectrum_df.Intensity[self.spectrum_df.j_AAIG == j_AAIG].max())
            # for Intensity in self.spectrum_df.Intensity[self.spectrum_df.j_AAIG == j_AAIG]:
            #     self.spectrum_df.scaled_Intensity[(self.spectrum_df.j_AAIG == j_AAIG) &
            #                                     (self.spectrum_df.Intensity == Intensity)] = Intensity / max_Intensity
        # By default, std() returns NaN for single values, but we want it to return 0.
        self.spectrum_df.fillna(
            value={'i_stdevNreson': 0, 'i_stdevHNreson': 0, 'j_stdevNreson': 0, 'j_stdevHNreson': 0}, inplace=True)

        # Finally, scale the intensities between 0 and 1.
        self.__scale_intensities()

    def find_sidechains(self, SC_NQ_num, tolHN=0.03, tolN=0.3, max_distratio=5.0):
        """
        Method to find ASN and GLN side chain AAIGs.

        :param SC_NQ_num:   total number of ASN and GLN sidechain AAIGs according to the fasta sequence.
        :param tolHN:
        :param tolN:
        :param max_distratio:
        :return:
        """
        # Load pre-trained Random Forest models
        CHAINS_LIB_DIR = os.path.dirname(os.path.realpath(__file__))
        with open(CHAINS_LIB_DIR + '/../ML_models_DELETE/asn_gln_model.pickle', 'rb') as infile: # NEW CLASSIFIER
            rf = pickle.load(infile)

        sc_row_indices = []
        for index, row in self.spectrum_df.iterrows():
            if approx_equal(row['i_Nreson'], row['j_Nreson'], tolN):
                sc_row_indices.append(index)

        sc_pairs = []
        for i in range(len(sc_row_indices)):
            index1 = sc_row_indices[i]
            for j in range(i+1, len(sc_row_indices)):
                index2 = sc_row_indices[j]
                row1 = self.spectrum_df.iloc[index1]
                row2 = self.spectrum_df.iloc[index2]
                # Name and Resonance condition. Example:
                #        54H-N-6N-H      7.658    113.821    113.877      6.701    113763480
                #        6H-N-54N-H      6.701    113.819    113.878      7.658    116716272
                if row1[['i_AAIG', 'j_AAIG']].tolist() == list(reversed(row2[['i_AAIG', 'j_AAIG']])) and \
                    approx_equal(row1['i_HNreson'], row2['j_HNreson'], tolHN) and \
                    approx_equal(row1['i_Nreson'], row2['i_Nreson'], tolN) and \
                    approx_equal(row1['j_Nreson'], row2['j_Nreson'], tolN) and \
                    approx_equal(row1['j_HNreson'], row2['i_HNreson'], tolHN):

                    # Check also the highest Intensity condition. The sidechain should be the strongest peak for the given AAIG.
                    AAIG1, AAIG2 = row1[['i_AAIG', 'j_AAIG']].tolist()
                    df1 = self.spectrum_df[(self.spectrum_df['i_AAIG'] == AAIG1) | (self.spectrum_df['j_AAIG'] == AAIG1)]
                    sorted_df1 = df1.sort_values(by=['Intensity'], ascending=False).reset_index()
                    df2 = self.spectrum_df[(self.spectrum_df['i_AAIG'] == AAIG2) | (self.spectrum_df['j_AAIG'] == AAIG2)]
                    sorted_df2 = df2.sort_values(by=['Intensity'], ascending=False).reset_index()

                    if set(sorted_df1.iloc[0][['i_AAIG', 'j_AAIG']]) == {AAIG1, AAIG2} and \
                        set(sorted_df2.iloc[0][['i_AAIG', 'j_AAIG']]) == {AAIG1, AAIG2}:
                        distance = self.distance_Series(row1, row2)
                        sc_pairs.append((index1, index2, distance))
                    else:
                        ColorPrint("WARNING: The following 2 HNNH peaks passed the Name and Resonance conditions but failed "
                                   "in the Intensity conditions:", "WARNING")
                        print(row1)
                        print(row2)

        # Sort by distance, and make sure that each row (peak) is "paired" only once!
        sc_pairs.sort(key=itemgetter(2))
        unique_sc_pairs = []
        used_indices = set()
        for index1,index2,distance in sc_pairs:
            if not index1 in used_indices and not index2 in used_indices:
                unique_sc_pairs.append((index1,index2,distance))
                used_indices.add(index1)
                used_indices.add(index2)
        del sc_pairs    # delete to avoid confusion. Only unique_sc_pairs must be used.

        # Populate the backbone spectrum Dataframe
        all_sc_indices = list(flatten(unique_sc_pairs))
        for index in range(self.spectrum_df.shape[0]):
            if not index in all_sc_indices:
                row = self.spectrum_df.iloc[[index]]
                self.bb_Wsc_spectrum_df = self.bb_Wsc_spectrum_df.append(row, ignore_index=True)

        # Populate the side-chain spectrum DataFrame
        for k in range(len(unique_sc_pairs)):
            i, j, d = unique_sc_pairs[k]
            row1 = self.spectrum_df.iloc[[i]].reset_index()
            row2 = self.spectrum_df.iloc[[j]].reset_index()
            distance = self.distance_DataFrame(row1, row2)
            row1.at[0, 'distance'] = distance
            row2.at[0, 'distance'] = distance
            self.NQsc_spectrum_df = self.NQsc_spectrum_df.append(row1, ignore_index=True)
            self.NQsc_spectrum_df = self.NQsc_spectrum_df.append(row2, ignore_index=True)

        # OBSOLETE: Sort by distance and keep thos above the mindist_ratio. The remaining transfer them
        #            to bb_spectrum_df
        sorted_sc_spectrum_df = self.NQsc_spectrum_df.sort_values(by=['distance']).reset_index()
        # Reset NQsc_spectrum_df
        self.NQsc_spectrum_df = pd.DataFrame([], columns=sorted_sc_spectrum_df.columns, dtype='float64')  # only sidechain N,HN pairs
        self.NQsc_spectrum_df = self.NQsc_spectrum_df.astype({'i_AAIG': 'str', 'j_AAIG': 'str'})
        sc_nq_num = 0   # the actual number of ASN+GLN SC peaks found [<= 2*(#ASN+#GLN)]
        # Apply the max_distratio condition (OPTIONAL). To deactivate set max_distratio=0
        # previous_distance = 0.0
        for index, row in sorted_sc_spectrum_df.iterrows():    # SC_NQ_num is the maximum number of ASN+GLN peaks
            pred_label = rf.predict(np.array(sorted(row[['i_HNreson', 'j_HNreson']]) +
                [np.mean(row[['i_aveNreson', 'j_aveNreson']])]).reshape(1, -1))[0]
            pred_proba = rf.predict_proba(np.array(sorted(row[['i_HNreson', 'j_HNreson']]) +
                [np.mean(row[['i_aveNreson', 'j_aveNreson']])]).reshape(1, -1))[0]
            row['prob_ASN'], row['prob_GLN'], row['prob_UNK'] = pred_proba

            if row['prob_UNK'] > 0.4:   # <== CHANGE ME
                # Namely, if the probability of this pair of HNNH lines NOT BEING ASN OR GLN SIDECHAIN.
                self.bb_Wsc_spectrum_df = self.bb_Wsc_spectrum_df.append(row, ignore_index=True)
            else:
                # if previous_distance > 0 and row['distance'] / previous_distance > max_distratio:
                #     break
                sc_nq_num += 1
                self.NQsc_spectrum_df = self.NQsc_spectrum_df.append(row, ignore_index=True)
                # previous_distance = row['distance']
        # self.bb_Wsc_spectrum_df = self.bb_Wsc_spectrum_df.append(sorted_sc_spectrum_df.iloc[sc_nq_num:],    # append the peaks that did not
        #                                                         # satisfy the max_distratio criterion to the bb spectrum.
        #                                                          ignore_index=True)

    def proof_read_line(self, row):
        """
        Special method just for proofreading assigned HNNH peak lists.
        :return:
        """
        # i_AAIG_name, i_NH_name = split_AAIG_signature(row['i_AAIG'])
        i_AAIG_name_N_name = row['i_AAIG'].split('-')[0]
        # j_AAIG_name, j_NH_name = split_AAIG_signature(row['j_AAIG'])
        j_AAIG_name_N_name = row['j_AAIG'].split('-')[0]
        if i_AAIG_name_N_name[0] in ['N', 'Q'] and i_AAIG_name_N_name == j_AAIG_name_N_name:
            return True
        else:
            return False

    def write_sc_bb_sparky_lists(self):

        with open(self.HNNH_FILE.replace("num.list", "_bbnum.list"), 'w') as f:
            for index, row in self.bb_Wsc_spectrum_df.iterrows():
                # f.write("%s-%s\t\t%f\t%f\t%f\t%f\t%f\n" % (row['i_AAIG'], row['j_AAIG'], row['i_HNreson'], row['i_Nreson'],
                #                                                row['j_Nreson'], row['j_HNreson'], row['Intensity'])
                #         )
                # Include Random Forest probabilities
                f.write("%s-%s\t\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\n" % (row['i_AAIG'], row['j_AAIG'], row['i_HNreson'], row['i_Nreson'],
                                                                       row['j_Nreson'], row['j_HNreson'], row['Intensity'],
                                                                       row['prob_ASN'], row['prob_GLN'], row['prob_UNK'],
                                                                        not self.proof_read_line(row))
                        )

        with open(self.HNNH_FILE.replace("num.list", "_scnum.list"), 'w') as f:
            for index, row in self.NQsc_spectrum_df.iterrows():
        #         f.write("%s-%s\t\t%f\t%f\t%f\t%f\t%f\n" % (row['i_AAIG'], row['j_AAIG'], row['i_HNreson'], row['i_Nreson'],
        #                                                        row['j_Nreson'], row['j_HNreson'], row['Intensity'])
        #                 )
                f.write("%s-%s\t\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\n" % (row['i_AAIG'], row['j_AAIG'], row['i_HNreson'], row['i_Nreson'],
                                                                       row['j_Nreson'], row['j_HNreson'], row['Intensity'],
                                                                       row['prob_ASN'], row['prob_GLN'], row['prob_UNK'],
                                                                       self.proof_read_line(row))
                        )

    def get_NQsc_AAIGs(self):
        """
        Each sidechain N of ASN and GLN has two H. Therefore, the expected number of NQsc_AAIGs is 2*(ASN+GLN).
        :return:
        """
        i_AAIGs = set(self.NQsc_spectrum_df['i_AAIG'])
        j_AAIGs = set(self.NQsc_spectrum_df['j_AAIG'])
        return i_AAIGs.union(j_AAIGs)

    def get_bb_Wsc_AAIGs(self):
        NQsc_AAIGs = self.get_NQsc_AAIGs()
        i_AAIGs = set(self.bb_Wsc_spectrum_df['i_AAIG'])
        j_AAIGs = set(self.bb_Wsc_spectrum_df['j_AAIG'])
        all_AAIGs = i_AAIGs.union(j_AAIGs)
        return all_AAIGs.difference(NQsc_AAIGs)

    def get_replicate_lines(self, fields=['i_HSQCdist', 'j_HSQCdist', 'test_assignment', 'i_agreement', 'j_agreement']):
        """
        Especially for training HSQC-HNNH matching RF.
        :return:
        """
        all_peaks = []  # these peaks are just tuples, not Peak objects
        for index, row in self.spectrum_df.iterrows():
            all_peaks.append((index, row['i_HNreson'], row['i_Nreson'],
                                row['j_Nreson'], row['j_HNreson'], row['Intensity']))
        all_peaks.sort(key=itemgetter(1,2,3,4,5))
        replicate_peak_indices = []   # the indices in the self.spectrum_df of the replicate lines
        peak_group = set()
        print(all_peaks)
        for i in range(len(all_peaks)-1):
            peak_group.add(all_peaks[i][0])
            if all_peaks[i][1:] == all_peaks[i+1][1:]:
                peak_group.add(all_peaks[i+1][0])
            else:
                replicate_peak_indices.append(peak_group)
                peak_group = set()

        replicate_lines = defaultdict(list) # group index -> [1st replicate line, 2nd replicate line, ...]
        print(replicate_peak_indices)
        for grp, group in enumerate(replicate_peak_indices):
            for index in group:
                row = self.spectrum_df.iloc[index]
                replicate_lines[grp].append([row[f] for f in fields])    # prepend the index of the group
        return replicate_lines  # as dict of word lists