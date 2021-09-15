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


from lib.connectivities.__init__ import *
from lib.proofread.print_functions import *
from copy import deepcopy


class HNNH_Connectivities(Connectivities):
    """
    TERM DEFINITIONS:
    "connectivity": the existence of one or more common peaks between two AAIGs, namely AAIG1 and AAIG2.
    "connection": a pair of peaks with matching C and H resonances, one in AAIG1 and one in AAIG2.
    Thereby, two AAIGs can have only one "connectivity" but many "connections".


    This Class is meant to find connectivities within the NOESY HCNH spectrum. The NOESY HNNH spectrum can be OPTIONALLY
    used to filter out connectivities between non-contiguous residues.

    Below is described the general strategy for finding connectivities using only 4D NOESY spectra (HCNH, HNNH):
    1. Match HSQC-HNNH 3rd and 4th column
    2. Match HSQC-HNNH 1st and 2nd column (lower resolution => more errors)
    3. Identify side-chain N-H (ASN,GLN). These have the same values in columns 1,3 and 2,4 and very high intensity.
       Use the intensity to identify the intra-residue peaks.
    4. Remove the side-chains from all the spectra (HSQC, HNNH, HCNH).
    5. * Redo HSQC-HNNH matching and HSQC-HCHN matching, this time without the side-chains.
    6. Find HNNH-HNNH connectivities.
    7. Find HCNH-HCNH connectivities and use the HNNH-HNNH connectivities to filter them.
    """

    def __init__(self, tolH=0.04, tolC=0.4, H_bin_length=0.0005, C_bin_length=0.005, bandwidth=0.2):

        super(HNNH_Connectivities, self).__init__()

        self.HNNH_peak_intensities_mdict = tree()  # AAIG(1) -> matching AAIG(2) -> peak intensity

    def clean_connectivities(self):
        ColorPrint("Cleaning all saved connectivities.", "BOLDGREEN")

        self.AAIG_connectivities_complete_dict = {}  # AAIG(i) -> [(NOESY AAIG(!=i), 1, tot. num. connectivities of AAIG(i),
                                                         # peak intensity, peak intensity), ...]

    def write_spectrum_info(self):
        if not self.HNNH_spec:
            raise Exception(bcolors.FAIL + "ERROR: you must load the NOESY spectrum in order to write statistics "
                                           "in file!" + bcolors.ENDC)
        with open("spectrum_statistics.txt", 'w') as f:
            AAIG_num, peak_num, ave_peaknum_per_aaig, min_peaknum_per_aaig, max_peaknum_per_aaig = \
                self.HNNH_spec.__get_spectrum_info__()
            f.write("\nNOESY:\n")
            f.write("AAIG number: %g\n" % AAIG_num)
            f.write("Peak number: %g\n" % peak_num)
            f.write("Average peak number per AAIG: %f\n" % ave_peaknum_per_aaig)
            f.write("Minimum peak number per AAIG: %i\n" % min_peaknum_per_aaig)
            f.write("Maximum peak number per AAIG: %i\n" % max_peaknum_per_aaig)

#####*  I moved these methods from __init__.py to here as they are relevent only to HCNH-HCNH.

    def write_HNNH_connectivities_to_txt(self,
                                         fname,
                                         pickle_fname="",
                                         selected_AAIG_connectivities_dict={}):
        """
        Method to write the connectivity present in self.HNNH_peak_connectivities_mdict into a file.
        If you want to pass connectivities from an external dictionary use method write_HNNH_connectivities_from_dict().
        If the selected_AAIG_connectivities_dict is provided, then the method will still write all connectivities in
        self.HNNH_peak_connectivities_mdict but when an AAIG is a key in both dicts then only the connectivities
        in selected_AAIG_connectivities_dict will be written. Therefore, if self.HNNH_peak_connectivities_mdict
        is empty or contains fewer connectivities than selected_AAIG_connectivities_dict, then this method is useless and
        dangerous!!!

        :param fname:
        :param pickle_fname:
        :param oexp:    3.0448217240563293 is globally optimized value in the set of 9 globular proteins for global intensity scaling
                        and peak-peak 2D-hist connectivities.
        :param selected_AAIG_connectivities_dict:   dict with selected connectivities of every NOESY AAIG1; has keys the
                                                    AAIG signatures and values lists of quintuplets (tuples) consisting of
                                                    (NOESY AAIG2 signature, # matching peaks, total # peaks, intersection,
                                                    multi score). If this dict is provided then this method will write
                                                    to files only the connectivites in this dict.

        :return:
        """
        if not pickle_fname:
            pickle_fname = os.path.splitext(fname)[0] + ".pkl"
        AAIG_connectivities_complete_dict = {}    # AAIG1 signature -> [(connected AAIG2 signature, 1,
                                        # number of AAIG1 connectivities, HNNH peak intensity,
                                        # peak intensity), ...]

        with open(fname, 'w') as f:
            f.write("# AAIG1 signature --> (connected AAIG2 signature, 1, total number of AAIG1 connectivities"
                    ", HNNH peak intensity, 'multi' score [combines occupancy and "
                    "intensity]), ...\n")
            HNNH_peak_intensities_mdict = deepcopy(self.HNNH_peak_intensities_mdict)
            for AAIG1_signature in list(self.HNNH_peak_intensities_mdict.keys()):
                if selected_AAIG_connectivities_dict and not AAIG1_signature in selected_AAIG_connectivities_dict.keys():
                    continue
                # NOTE: MULTI score is pointless in HNNH because the occupancy of all connectivities of AAIG1 is the same.
                # NOTE: for compatibility with other types of connectivities I use the intensity instead.
                connectivities_list = [(AAIG2_signature,
                                        1,
                                        len(self.HNNH_peak_intensities_mdict[AAIG1_signature].keys()),
                                        self.HNNH_peak_intensities_mdict[AAIG1_signature][AAIG2_signature],
                                        self.HNNH_peak_intensities_mdict[AAIG1_signature][AAIG2_signature])
                                       for AAIG2_signature in list(self.HNNH_peak_intensities_mdict[AAIG1_signature].keys())]
                connectivities_list.sort(key=itemgetter(3), reverse=True)   # sort by peak intensity
                connectivities_list = [(c[0], c[1], c[2], c[3], c[4]) for c in connectivities_list]
                if selected_AAIG_connectivities_dict:   # FILTER: keep only the selected connectivities
                    connAAIG2_list = [conn[0] for conn in selected_AAIG_connectivities_dict[AAIG1_signature]]
                    connectivities_list = [conn for conn in connectivities_list if conn[0] in connAAIG2_list]
                    for AAIG2_signature in self.HNNH_peak_intensities_mdict[AAIG1_signature].keys():
                        if not AAIG2_signature in connAAIG2_list:
                            del HNNH_peak_intensities_mdict[AAIG1_signature][AAIG2_signature]
                if len(connectivities_list) > 0:    # don't write AAIGs without valid connectivities
                    f.write(AAIG1_signature + "\t" + ', '.join(str(c) for c in connectivities_list) + "\n")
                AAIG_connectivities_complete_dict[AAIG1_signature] = connectivities_list
        # Save the same information in a pickle file, along withthe dictionarry with the all matching peaks between pairs
        # of AAIGs.
        save_pickle(pickle_fname, AAIG_connectivities_complete_dict,
                    HNNH_peak_intensities_mdict)  # save the filtered versions of both dicts

    def __get_HNNH_connectivities_dataframe__(self, selected_AAIG_connectivities_dict={}):
        """
        New method for nmr_pipeline.

        :param selected_AAIG_connectivities_dict:   dict with selected connectivities of every NOESY AAIG1; has keys the
                                            AAIG signatures and values lists of quintuplets (tuples) consisting of
                                            (NOESY AAIG2 signature, # matching peaks, total # peaks, intersection,
                                            multi score). If this dict is provided then this method will write
                                            to files only the connectivites in this dict.
        :return:
        """
        conn_data = []  # [ ["frame1", "frame2", "number of common peaks", "occupancy", "2D-hist intersection"], ...]
        # self.HNNH_peak_intensities_mdict: AAIG signature (1) -> matching AAIG signature (2) -> peak intensity
        for AAIG1_signature in list(self.HNNH_peak_intensities_mdict.keys()):
            if selected_AAIG_connectivities_dict and AAIG1_signature not in selected_AAIG_connectivities_dict.keys():
                continue
            # NOTE: for compatibility with other types of connectivities I use the intensity instead.
            # frame1, frame2, 1, total num. of frame1's peaks, scaled peak's intensity
            total_peaknum = len(self.HNNH_peak_intensities_mdict[AAIG1_signature].keys())
            conn_list =  [
                            [
                              AAIG1_signature,
                              AAIG2_signature,
                              1,
                              total_peaknum,
                              self.HNNH_peak_intensities_mdict[AAIG1_signature][AAIG2_signature]
                            ]
                            for AAIG2_signature in
                            list(self.HNNH_peak_intensities_mdict[AAIG1_signature].keys())
                        ]
            conn_list.sort(key=itemgetter(4), reverse=True)  # sort frame1's connectivities by scaled peaks' intensity
            if selected_AAIG_connectivities_dict:  # FILTER: keep only the selected connectivities
                connAAIG2_list = [conn[0] for conn in selected_AAIG_connectivities_dict[AAIG1_signature]]   # list only the AAIG2 signatures
                conn_list = [conn for conn in conn_list if conn[1] in connAAIG2_list]   # conn[1] because we save the AAIG1 signature first (see above)
            conn_data.extend(conn_list)

        conn_df = pd.DataFrame(data=conn_data,
                               columns=["source_frame", "target_frame", "num_common_peaks",
                                        "total_num_peaks", "scaled_intensity"])
        conn_df = conn_df.round(4)
        return conn_df


####### *

    @staticmethod
    def write_HNNH_connectivities_from_dict(AAIG_connectivities_dict, fname):
        """
        Method to write the contents of a connectivity dict to a file.

        :param AAIG_connectivities_dict:    connectivity dict to be written to the file.
        :param fname:
        :return:
        """

        with open(fname, 'w') as f:
            f.write("# AAIG1 signature --> (connected AAIG2 signature, always 1.0, total number of "
                    "AAIG1 connectivites, AAIG1-AAIG2 peak scaled intensity, 'multi' score [combines occupancy and "
                    "peak intensity]), ...\n")    # header
            for k, v in list(AAIG_connectivities_dict.items()):
                # print k,v
                f.write(k + "\t" + ', '.join(map(str, v)) + "\n")

    def find_HNNH_connectivities(self, selected_AAIG_connectivities_dict=None):
        """
        Method to find connectivities from HNNH peak labels.

        :return: it populates self.HNNH_peak_intensities_mdict
        """
        # Clean previously stored connectivities
        self.clean_connectivities()

        # TODO: integrate selected_AAIG_connectivities_dict into the code.
        # Instead of 2D-hist intersections, save the HNNH peak intensities
        for i, row in self.HNNH_spec.spectrum_df.iterrows():
            AAIG1_signature = row['i_AAIG'].replace('-','')
            AAIG2_signature = row['j_AAIG'].replace('-','')
            if AAIG1_signature == AAIG2_signature:
                continue
            self.HNNH_peak_intensities_mdict[AAIG1_signature][AAIG2_signature] = row['scaled_Intensity']
        # Check if the reverse connectivities exist, and if not, save them with the same intensity
        for i, row in self.HNNH_spec.spectrum_df.iterrows():
            AAIG1_signature = row['i_AAIG'].replace('-', '')
            AAIG2_signature = row['j_AAIG'].replace('-', '')
            if AAIG1_signature == AAIG2_signature:
                continue
            if AAIG2_signature not in self.HNNH_peak_intensities_mdict.keys() \
                    or AAIG1_signature not in self.HNNH_peak_intensities_mdict[AAIG2_signature].keys():
                self.HNNH_peak_intensities_mdict[AAIG2_signature][AAIG1_signature] = row['scaled_Intensity']

        # DEACTIVATE IT FOR NOW.
        # # Finally remove connectivities with 0 intensity
        # self.deconvolute_HNNH_connectivities()

    def deconvolute_HNNH_connectivities(self):
        """
        Remove connectivities with 0 intensity.
        :return:
        """
        # TODO: if a connectivity peak has intensity 0, delete also the reverse connectivity???
        for AAIG1_signature in list(self.HNNH_peak_intensities_mdict.keys()):
            # if len(self.NOESY_2Dhist_intersections_mdict[AAIG1_signature].values()) == 0:
            #     continue
            # max_intersection = np.max(self.NOESY_2Dhist_intersections_mdict[AAIG1_signature].values())
            for AAIG2_signature in list(self.HNNH_peak_intensities_mdict[AAIG1_signature].keys()):
                if self.HNNH_peak_intensities_mdict[AAIG1_signature][AAIG2_signature] == 0: #/max_intersection < 0.00001:   # <== CHANGE ME
                    del self.HNNH_peak_intensities_mdict[AAIG1_signature][AAIG2_signature]
            if len(list(self.HNNH_peak_intensities_mdict[AAIG1_signature].keys())) == 0:    # remove the NOESY AAIG if no connectivities
                                                                                            # were left
                del self.HNNH_peak_intensities_mdict[AAIG1_signature]


    def get_best_NOESY_connectivity(self, connections_dict, AAIG1_AAIG2, parallel=False):
        """
        Method to find the maximum possible number of matching peaks between NOESY AAIG1 and NOESY AAIG2.
        :param connections_dict:    contains as keys all the peaks of AAIG1 and values list of
                                    (matching AAIG2 peak, C-H distance)
        :param AAIG1_AAIG2: string that contains the signatures of AAIG1 and AAIG2 which have matching peaks
        :param parallel:    highly parallel more (worthy only if AAIG1 has more that 4 matching peaks with AAIG2)
        :return: AAIG1_signature:
        :return: AAIG2_signature:
        :return: best_chain: the optimum list of [(AAIG1 Peak object, matching AAIG2 Peak object), ...] (without the distance)
        """
        # TODO: adapt it for HNNH or remove it!
        AAIG1_signature, AAIG2_signature = AAIG1_AAIG2.split("_")
        con_tree = Connectivities_Tree(connections_dict)
        if parallel:
            # matching peaks of AAIG1 I mean.
            ColorPrint("Uniquifying big connectivities (%i matching peaks) between NOESY AAIG %s and NOESY AAIG %s " %
                       (len(connections_dict), AAIG1_signature, AAIG2_signature), "OKBLUE")
            con_tree.form_chains_parallel(length='max')
            ColorPrint("Finished big connectivities between NOESY AAIG " + AAIG1_signature + \
                       " and NOESY AAIG " + AAIG2_signature, "OKBLUE")
        else:
            ColorPrint("Uniquifying small connectivities (%i matching peaks) between NOESY AAIG %s and NOESY AAIG %s " %
                       (len(connections_dict), AAIG1_signature, AAIG2_signature), "OKBLUE")
            con_tree.form_chains(length='max')
            ColorPrint("Finished small connectivities between NOESY AAIG " + AAIG1_signature + \
                  " and NOESY AAIG " + AAIG2_signature, "OKBLUE")
        return AAIG1_signature, AAIG2_signature, con_tree.get_best_chain(length='max')


######################################### START OF GETTER METHODS ########################################

    def __get_HNNH_connectivities_dict__(self):
        """
        Getter method that returns all possible connectivities of every NOESY AAIG with other NOESY AAIGs in a dictionary form.
        It presumes that match_NOESY_to_NOESY_Peaks() has been already called, otherwise it will return an empty
        dictionary.

        :return: AAIG_connectivities_complete_dict: dictionary of the form: NOESY AAIG1 signature --> [(NOESY AAIG2 signature,
                                            occupancy, numOfResonances, intersection, multi score), ...]
        """

        for AAIG1_signature in list(self.HNNH_peak_intensities_mdict.keys()):
            # NOTE: MULTI score is pointless in HNNH because the occupancy of all connectivities of AAIG1 is the same.
            # NOTE: for compatibility with other types of connectivities I use the intensity instead.
            try:
                connectivities_list = [(AAIG2_signature,
                                        1,
                                        len(self.HNNH_peak_intensities_mdict[AAIG1_signature].keys()),
                                        self.HNNH_peak_intensities_mdict[AAIG1_signature][AAIG2_signature],
                                        self.HNNH_peak_intensities_mdict[AAIG1_signature][AAIG2_signature])
                                       for AAIG2_signature in list(self.HNNH_peak_intensities_mdict[AAIG1_signature].keys())]

            except IndexError:
                # if AAIG2_signature == "" or len(connectivities_list) == 0:
                #     pass
                # else:
                print("DEBUG: AAIG1_signature=", AAIG1_signature)
                print("DEBUG: AAIG2_signature=", AAIG2_signature)
                print("DEBUG: type(AAIG2_signature)=", type(AAIG2_signature))
                print("DEBUG: connectivities_list=", connectivities_list)
                raise IndexError()
            connectivities_list.sort(key=itemgetter(3), reverse=True)   # sort by peak intensity
            connectivities_list = [(c[0], c[1], c[2], c[3], c[4]) for c in connectivities_list]
            self.AAIG_connectivities_complete_dict[AAIG1_signature] = connectivities_list
        return self.AAIG_connectivities_complete_dict

    def __get_filtered_HNNH_connectivities_dict__(self,
                                                 INTENSITY_RATIO=10,
                                                 reverse_conn=False,
                                                 proofread=False,
                                                 pool_AAIG_connectivities_dict=None):
        """
        Filter the self.AAIG_connectivities_complete_dict, by applying the INTENSITY_RATIO criterion. Also, for any two AAIGs,
        AAIG1 and AAIG2, look if both connectivities AAIG1-AAIG2 and AAIG2-AAIG1 exist in the produced chains under the
        given INTENSITY_RATIO. If not, then remove the one that exists (either AAIG1-AAIG2, or AAIG2-AAIG1).

        :param INTENSITY_RATIO: keep all connected AAIG2s that have 'intensity' up to INTENSITY_RATIO times lower than the maximum
                            'intensity' of the connectivities of AAIG1. The higher the INTENSITY_RATIO, the more connectivities
                            are retained.
        :param reverse_conn:    if NOT both AAIG1-AAIG2 and AAIG2-AAIG1 exist, then remove any of them from the i_iminus1_dict
        :param proofread:   print ERROR whenever a wrong connectivity is removed by the REVERSE CONNECTIVITY criterion.
        :param pool_AAIG_connectivities_dict: if this connectivities dict is provided, then it will be used instead
                                              self.AAIG_connectivities_complete_dict and be filtered by the INTENSITY_RATIO &
                                              reverse_conn criteria.
        :return i_iminus1_dict: dict with selected connectivities of every NOESY AAIG1; has keys the
                AAIG signatures and values lists of quintuplets (tuples) consisting of (NOESY AAIG2 signature, # matching peaks,
                total # peaks, intersection, intensity).
        """
        ColorPrint("Filtering the HCNH-HCNH connectivities by applying the INTENSITY_RATIO criterion.", "BOLDGREEN")
        if pool_AAIG_connectivities_dict == None:   # if no refined connectivity pool is provided, use all possible connectivities
            pool_AAIG_connectivities_dict = self.__get_HNNH_connectivities_dict__()

        num_before = len(list(flatten(pool_AAIG_connectivities_dict.values())))
        # Iterate over all NOESY aa indices and, if applicable, keep those that exceed the intensity threshold
        i_iminus1_dict = {}
        AAIG_AAIG_connBool_mdict = tree()  # AAIG1->AAIG2->True/False depending on the existence of this connectivity
        for i_AAIG_signature, connectivities in list(pool_AAIG_connectivities_dict.items()):
            filtered_connectivities = []
            max_intensity = np.max([float(q[4]) for q in connectivities])  # max intensity for this AAIG i
            for quintuplet in connectivities:
                j_AAIG_signature = quintuplet[0]
                AAIG_AAIG_connBool_mdict[i_AAIG_signature][j_AAIG_signature] = False  # initialize this dict. It will be updated later
                AAIG_AAIG_connBool_mdict[j_AAIG_signature][i_AAIG_signature] = False  # initialize this dict. It will be updated later
                intensity = float(quintuplet[3])    # multi in HNNH is the intensity
                if max_intensity / intensity <= INTENSITY_RATIO and intensity > 10e-6:  # <== CHANGE ME
                    filtered_connectivities.append(quintuplet)
                elif proofread and is_valid_signature(i_AAIG_signature) and is_valid_signature(j_AAIG_signature) \
                        and abs(get_resid(i_AAIG_signature) - get_resid(j_AAIG_signature)) == 1 and not \
                        is_sidechain(i_AAIG_signature) and not is_sidechain(j_AAIG_signature):
                    ColorPrint("ERROR: correct connectivity %s-%s was removed by INTENSITY_RATIO criterion!" %
                               (i_AAIG_signature, j_AAIG_signature),
                               "OKRED")
            i_iminus1_dict[i_AAIG_signature] = filtered_connectivities
        num_after = len(list(flatten(i_iminus1_dict.values())))
        ColorPrint("The number of connectivities was reduced from %i to %i using INTENSITY_RATIO=%.1f, "
                   "namely by %.3f %%." %
                   (num_before, num_after, INTENSITY_RATIO, 100 * (num_before - num_after) / num_before),
                   "OKGREEN")

        if reverse_conn:
            ColorPrint("Filtering the HCNH-HCNH connectivities by applying the REVERSE CONNECTIVITY criterion.", "BOLDGREEN")
            # Make sure that  both connectivities AAIG1-AAIG2 and AAIG2-AAIG1 passed the INTENSITY_RATIO criterion
            for i_AAIG_signature, connectivities in i_iminus1_dict.items():
                for quintuplet in connectivities:
                    j_AAIG_signature = quintuplet[0]
                    AAIG_AAIG_connBool_mdict[i_AAIG_signature][j_AAIG_signature] = True     # set right values to the dict
            i_iminus1_dict_copy = deepcopy(i_iminus1_dict)
            for i_AAIG_signature, connectivities in i_iminus1_dict_copy.items():
                filtered_connectivities = []
                for quintuplet in connectivities:
                    j_AAIG_signature = quintuplet[0]
                    if AAIG_AAIG_connBool_mdict[i_AAIG_signature][j_AAIG_signature] and \
                        AAIG_AAIG_connBool_mdict[j_AAIG_signature][i_AAIG_signature]:    # both AAIG1-AAIG2 and AAIG2-AAIG1 must exist
                        filtered_connectivities.append(quintuplet)
                    elif proofread and is_valid_signature(i_AAIG_signature) and is_valid_signature(j_AAIG_signature) \
                            and abs(get_resid(i_AAIG_signature)-get_resid(j_AAIG_signature)) == 1 and not \
                            is_sidechain(i_AAIG_signature) and not is_sidechain(j_AAIG_signature):
                        ColorPrint("ERROR: correct connectivity %s-%s was removed by REVERSE CONNECTIVITY criterion!" %
                                   (i_AAIG_signature, j_AAIG_signature),
                                   "OKRED")
                # TODO: check if you must delete the key if the values are an empty list.
                i_iminus1_dict[i_AAIG_signature] = filtered_connectivities


        num_after2 = len(list(flatten(i_iminus1_dict.values())))
        ColorPrint("The number of connectivities was further reduced to %i using reverse_conn=%s, "
                   "namely by %.3f %%." %
                   (num_after2, reverse_conn, 100*(num_after-num_after2)/num_before), "OKGREEN")
        return i_iminus1_dict

############################################# START OF PROOF-READING METHODS ####################################################

    def __get_intensity_position_correlation__(self):
        """
        TODO: method that counts how many times the intensity of i->i+1 connectivity is higher than the intensity of i->i-1
        TODO: and the opposite.
        :return:
        """
        num_iplus1_gt_iminus1 = 0    # intensity(i->i+1)
        num_iplus1_lt_iminus1 = 0

        for AAIG1_signature in self.HNNH_peak_intensities_mdict.keys():
            if not is_valid_signature(AAIG1_signature) or is_sidechain(AAIG1_signature):
                continue
            # AAIG1_signature=i, now search the AAIG signatures of i-1 and i+1
            iplus1_signature, iminus1_signature = None, None
            for AAIG2_signature in self.HNNH_peak_intensities_mdict[AAIG1_signature].keys():
                if not is_valid_signature(AAIG2_signature) or is_sidechain(AAIG2_signature):
                     continue
                if get_resid(AAIG1_signature) == get_resid(AAIG2_signature) -1: # AAIG1=i, AAIG2=i+1
                    iplus1_signature = AAIG2_signature
                elif get_resid(AAIG1_signature) == get_resid(AAIG2_signature) +1: # # AAIG1=i, AAIG2=i-1
                    iminus1_signature = AAIG2_signature
            # If connectivities for both i-1 and i+1 exist, then compare their intensities
            if iplus1_signature and iminus1_signature:
                if self.HNNH_peak_intensities_mdict[AAIG1_signature][iplus1_signature] > self.HNNH_peak_intensities_mdict[AAIG1_signature][iminus1_signature]:
                    num_iplus1_gt_iminus1 += 1
                elif self.HNNH_peak_intensities_mdict[AAIG1_signature][iplus1_signature] < self.HNNH_peak_intensities_mdict[AAIG1_signature][iminus1_signature]:
                    num_iplus1_lt_iminus1 += 1
        ColorPrint("The HNNH connectivity i->i+1 had greater intensity than the connectivity i->i-1 %i times, while lower %i times." %
                   (num_iplus1_gt_iminus1, num_iplus1_lt_iminus1), "BOLDGREEN")

# ## JUST FOR TESTING OF THE CLASS ##
# if __name__ == '__main__':
#
#     INTENSITY_RATIO = 5
#     HNNH_con = HNNH_Connectivities()
#     HNNH_con.load_HNNH_spectrum("/home2/thomas/Documents/4D-CHAINS_regtests/FD3A/ASN_GLN_sc_assignment/FD3A_HNNHnum.list")
#     HNNH_con.find_HNNH_connectivities()
#     AAIG_connectivities_dict = HNNH_con.__get_filtered_HNNH_connectivities_dict__(INTENSITY_RATIO=INTENSITY_RATIO,
#                                                                                  reverse_conn=True, # DANGEROUS: removes correct connectivities
#                                                                                  proofread=True)
#     HNNH_con.write_HNNH_connectivities(fname="/home2/thomas/Documents/4D-CHAINS_regtests/FD3A/ASN_GLN_sc_assignment/HNNH_connectivities_intratio%.1f" % INTENSITY_RATIO,
#                                        selected_AAIG_connectivities_dict=AAIG_connectivities_dict)




