from lib.connectivities import *
from collections import defaultdict
from lib.proofread import *
from lib.global_func import *
import pandas as pd
from scoop import futures, shared


class HCNH_Connectivities(Connectivities):
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

        super(HCNH_Connectivities, self).__init__()

        self.NOESY_NOESY_peak_connectivities_mdict = tree()  # AAIG(1) signature -> matching AAIG(2) signature ->
                                                            # [(matching NOESY peak1, matching NOESY peak 2, distance), ...]
        self.NOESY_NOESY_2Dhist_intersections_mdict = tree()  # AAIG(1) -> matching AAIG(2) -> 2D-hist intersection
        self.AAIG_recurrent_peaks_mdict = tree()    # AAIG_signature -> Peak -> frequency of recurrence

    def clean_connectivities(self):
        ColorPrint("Cleaning all saved connectivities.", "BOLDGREEN")
        self.NOESY_NOESY_peak_connectivities_mdict = tree()  # AAIG(1) signature -> matching AAIG(2) signature ->
        # [(matching NOESY peak1, matching NOESY peak 2, distance), ...]
        self.AAIG_connectivities_complete_dict = {}  # NOESY AAIG(i) -> [(NOESY AAIG(!=i), occupancy, tot. num. peaks, intersection, multi score),
        # ...]

    def clean_intersections(self):
        print(bcolors.BOLDGREEN + "Cleaning all saved intersections." + bcolors.ENDBOLD)
        self.NOESY_2Dhist_connectivities_mdict = tree()  # NOESY AAIG(1) -> matching NOESY AAIG(2) -> 2D-hist intersection

    def __scale_intensities__(self, wrt_aaig=False):
        """
        Method used only by NOESY-connectivities.
        :param wrt_aaig:    if True scale the intensities wrt the AAIG, otherwise wrt the whole spectrum maximum.
        :return:
        """
        self.HCNH_spec.__scale_intensities__(even_intensities=False, wrt_aaig=wrt_aaig)

    def scale_NOESY_NOESY_2Dhist_connectivities(self, AAIG1_signature):
        """
        Method to scale the intersections of the connectivities of NOESY AAIG1 to be between 0 and 1.0.
        :param AAIG1:
        :return:
        """
        if len(list(self.NOESY_NOESY_2Dhist_intersections_mdict[AAIG1_signature].values())) == 0: # if no connectivities were found
            return
        max_intersection = np.max(list(self.NOESY_NOESY_2Dhist_intersections_mdict[AAIG1_signature].values()))
        if max_intersection == 0:   # if all connectivities have intersection 0, leave them as they are (they will be removed)
            return
        for AAIG2_signature in list(self.NOESY_NOESY_2Dhist_intersections_mdict[AAIG1_signature].keys()):
            self.NOESY_NOESY_2Dhist_intersections_mdict[AAIG1_signature][AAIG2_signature] /= max_intersection

    def save_NOESY_2Dhist_connectivity(self, AAIG1, AAIG2):
        """
        Method to calculate and save the intersection (value 0-1) between AAIG1 & AAIG2 expressed as
        2D density histograms.
        :param AAIG1:
        :param AAIG2:
        :return:
        """
        intersection = AAIG1.calc_intersection(AAIG2)
        self.NOESY_NOESY_2Dhist_intersections_mdict[AAIG1.signature][AAIG2.signature] = intersection

#####*  I moved these methods from __init__.py to here as they are relevant only to HCNH-HCNH.

    ####### *

    def save_NOESY_NOESY_peak_connectivity(self, AAIG1_signature, peak1, AAIG2_signature, peak2):
        """
        Method to save all possible connectivities and their distances into dataframes. You must call
        uniquify_NOESY_NOESY_peak_connectivities() afterwards to find the maximum possible
        number of peak pairs!

        :param AAIG1_signature:  the NOESY AAIG name
        :param peak1:
        :param AAIG2_signature:  the NOESY AAIG name
        :param peak2:
        :return:
        """
        distance = peak1.distance(peak2, attypes=['C', 'H'])
        # print("DEBUG: saving peak1", peak1.__get_CHresons__(), "peak2", peak2.__get_CHresons__(), "distance", distance)

        try:
            df = self.NOESY_NOESY_peak_connectivities_mdict[AAIG1_signature][AAIG2_signature]
            df = df.append({'NOESY_peak1': peak1, 'NOESY_peak2': peak2, 'distance': distance}, ignore_index=True)
            self.NOESY_NOESY_peak_connectivities_mdict[AAIG1_signature][AAIG2_signature] = df  # save changes
            # NOTE: I think that adding the reverse connectivity (peak2-peak1) is not necessary, as it will be addeded
            # NOTE: since match_one_NOESY_AAIG_against_all_NOESY_Peaks() iterates over all NOESY peaks.
        except (AttributeError, KeyError):
            # Create the dataframe if it doesn't exist.
            self.NOESY_NOESY_peak_connectivities_mdict[AAIG1_signature][AAIG2_signature] = \
                pd.DataFrame([(peak1, peak2, distance)], columns=("NOESY_peak1", "NOESY_peak2", "distance" ))

    def match_one_NOESY_AAIG_against_all_NOESY_Peaks(self, AAIG1, tolH=0.04, tolC=0.4, selected_AAIG_connectivities_dict=None):
        """

        :param AAIG1: NOESY AAIG
        :param tolH:
        :param tolC:
        :param selected_AAIG_connectivities_dict:
        :return:
        """
        AAIG1_signature = AAIG1.signature
        if selected_AAIG_connectivities_dict and not AAIG1_signature in selected_AAIG_connectivities_dict:
            return
        ColorPrint("Searching for the HCNH-HCNH connectivities of AAIG " + AAIG1_signature, "OKBLUE")
        for peak1 in AAIG1.__get_all_peaks__():
            for AAIG2_signature in self.HCNH_spec.__get_all_AAIG_signatures__():
                if AAIG1_signature == AAIG2_signature:  continue
                if selected_AAIG_connectivities_dict and AAIG2_signature not in [q[0] for q in selected_AAIG_connectivities_dict[AAIG1_signature]]:
                    continue
                AAIG2 = self.HCNH_spec.__get_AAIG__(AAIG2_signature)
                for peak2 in AAIG2.__get_all_peaks__():
                    # if peak2.Hreson == 0.000 and peak1.does_peak_match(peak2, tolH=tolH, tolC=tolC,
                    #                                 attypes=['C', 'H']):
                    #     print("DEBUG: peaks (%.3f, %.3f) and (%.3f, %.3f) match? " %
                    #           (peak1.Creson, peak1.Hreson, peak2.Creson, peak2.Hreson),
                    #           peak1.does_peak_match(peak2, tolH=tolH, tolC=tolC,
                    #                                 attypes=['C', 'H']))
                    #     peak1.__print_peak__()
                    #     peak2.__print_peak__()
                    if peak1.does_peak_match(peak2, tolH=tolH, tolC=tolC, attypes=['C', 'H']):
                        self.save_NOESY_NOESY_peak_connectivity(AAIG1_signature, peak1, AAIG2_signature, peak2)  # store this connectivity

    def measure_all_NOESY_NOESY_Peak_2Dhist_intersections(self, tolH=0.04, tolC=0.4, selected_AAIG_connectivities_dict=None, set_all_to_1=False):
        """
        For each NOESY AAIG1, look at the Peak connectivities and calculate the overall 2D-histogram
        intersection for each connectivity.

        :param tolH:
        :param tolC:
        :return:
        """

        for AAIG1_signature in list(self.NOESY_NOESY_peak_connectivities_mdict.keys()):
            if selected_AAIG_connectivities_dict and not AAIG1_signature in selected_AAIG_connectivities_dict:
                return
            AAIG1 = self.HCNH_spec.__get_AAIG__(AAIG1_signature)
            ColorPrint("Measuring Peak-Peak 2D-histogram intersection for the connectivities of NOESY AAIG " \
                  + AAIG1_signature, "OKBLUE")
            for AAIG2_signature in list(self.NOESY_NOESY_peak_connectivities_mdict[AAIG1_signature].keys()):
                if AAIG1_signature == AAIG2_signature:  continue
                if selected_AAIG_connectivities_dict and AAIG2_signature not in [q[0] for q in selected_AAIG_connectivities_dict[AAIG1_signature]]:
                    continue
                AAIG2 = self.HCNH_spec.__get_AAIG__(AAIG2_signature)
                if set_all_to_1:
                    self.NOESY_NOESY_2Dhist_intersections_mdict[AAIG1_signature][AAIG2_signature] = 1.0
                else:
                    total_intersection = 0.0
                    try:
                        for peak1, peak2 in self.NOESY_NOESY_peak_connectivities_mdict[AAIG1_signature][AAIG2_signature]:
                            updated_peak1 = AAIG1.__get_peak__(peak1)   # peaks1 and peaks2 are older version that don't have hist2D
                            updated_peak2 = AAIG2.__get_peak__(peak2)
                            total_intersection += updated_peak1.calc_intersection(updated_peak2)
                    except TypeError:
                        raise TypeError(Debuginfo("%s and %s had None connectivities in NOESY_NOESY_peak_connectivities_mdict!" %
                                                  (AAIG1_signature, AAIG2_signature) ,fail=True))
                    # Save the total intersection between matching Peaks of AAIG1 and AAIG2
                    self.NOESY_NOESY_2Dhist_intersections_mdict[AAIG1_signature][AAIG2_signature] = total_intersection

            self.scale_NOESY_NOESY_2Dhist_connectivities(AAIG1_signature)


    def match_NOESY_to_NOESY_Peaks(self,
                                   tolH=0.04,
                                   tolC=0.4,
                                   calc_intersection=True,
                                   uniquify=True,
                                   selected_AAIG_connectivities_dict=None,
                                   remove_matched_AAIGs_without_CA=True,
                                   min_peak_num=1,
                                   debug=False):
        """
        Method to find connectivities of each AAIG in NOESY with all other AAIGs in NOESY by comparing the 2D density
        histograms of their matching __Peaks__. The HCNH-HCNH connectivities are stored in multidict
        self.NOESY_peak_connectivities_mdict. To expedite the connectivities calculation, you can provide a
        dict (selected_AAIG_connectivities_dict) with a subset of connectivities that you want to calculate including the intersection
        and multi values.

        THEORY:
        each NOESY 4D peak has C,H of i,i+1,i-1,etc. and N,HN of i.
        The NOEYSY(i-1)-NOESY(i) & NOESY(i)-NOESY(i+1) matching is based on the fact that each NOESY AAIG i contains peaks
        of i-1, and i+1, too. As such, we iterate over all NOESY AAIGs and label them according to their N,HN label (i)
        and compare them with all NOESY AAIGs which we label according to their N,HN label. By definition the
        NOESY AAIG i will match with the NOESY AAIG i-1 and i+1. If we remove the label of AAIG i from the NOESY matches,
        we have the i-1, the i+1 and the rest. These i-1 and i+1 label corresponds to the NOESY C,H shifts, and hence in the
        output connectivities file we have 2 columns:
            i       possible i-1 and i+1

        :return:
        """
        # Clean previously stored connectivities
        self.clean_connectivities()

        # First, find connectivities without 2D-histogram comparison
        for AAIG1 in self.HCNH_spec.__get_all_AAIGs__():
            self.match_one_NOESY_AAIG_against_all_NOESY_Peaks(AAIG1, tolH=tolH, tolC=tolC,
                                                              selected_AAIG_connectivities_dict=selected_AAIG_connectivities_dict)

        if uniquify:
            # remove redundancy
            self.uniquify_NOESY_NOESY_peak_connectivities(selected_AAIG_connectivities_dict=selected_AAIG_connectivities_dict)
        else:   # NECESSARY FOR LATER STEPS!
            self.convert_NOESY_NOESY_peak_connectivities_mdict_values_to_lists()

        # Compare the 2D-histograms of each Peak pair to find the overall intersection of each connectivity
        if calc_intersection:
            # Clean previously stored intersections
            self.clean_intersections()
            # Scale intensities
            # IMPORTANT: self.scale_intensities() must be called outside of this function for the updates to be applied to
            # IMPORTANT: self.noesy_spec under parallel execution!

            # First generate the 2D-histograms for all AAIGs in the spectrum
            if selected_AAIG_connectivities_dict:
                NOESY_AAIG_signatures = [q[0] for v in list(selected_AAIG_connectivities_dict.values()) for q in v] + list(selected_AAIG_connectivities_dict.keys())
                NOESY_AAIG_signatures = list(set(NOESY_AAIG_signatures))
            else:
                NOESY_AAIG_signatures = self.HCNH_spec.__get_all_AAIG_signatures__()
            # if debug == False:
            #     try:  # PARALLEL EXECUTION WITH SCOOP
            #         shared.setConst(EVEN_EDGES=True,
            #                         SCALE_DENSITY=False)
            #     except TypeError:   # "TypeError: This constant already exists: EVEN_EDGES." (means we already invoked this method before)
            #         pass
            ColorPrint("Generating the 2D-histograms of all NOESY Peaks that participate in connectivities.", "BOLDBLUE")
            if debug:
                ## Serial NOESY Execution for debugging
                for AAIG_signature in NOESY_AAIG_signatures:
                    new_AAIG = self.NOESY_Peaks2hist(AAIG_signature, even_edges=True, scale_density=False)
                    self.noesy_spec.__replace_AAIG__(new_AAIG)
            else:
                ## Parallel NOESY Execution
                new_AAIGs = list(futures.map(self.NOESY_Peaks2hist,
                                             NOESY_AAIG_signatures,
                                             [True] * len(NOESY_AAIG_signatures),
                                             [False] * len(NOESY_AAIG_signatures)
                                             )
                                 )
                for aaig in new_AAIGs:  self.HCNH_spec.__replace_AAIG__(aaig)

            # Now that you have the Peak 2D-histograms measure all the intersections!
            self.measure_all_NOESY_NOESY_Peak_2Dhist_intersections(tolH=tolH, tolC=tolC,
                                                                   selected_AAIG_connectivities_dict=selected_AAIG_connectivities_dict)

            # # Release memory
            # self.noesy_spec.__delete_all_hist2D__()
            # self.noesy_spec.__delete_all_hist2D__()
            # gc.collect()
        else:   # set all intersections to 1.0
            self.measure_all_NOESY_NOESY_Peak_2Dhist_intersections(tolH=tolH, tolC=tolC,
                                                                   selected_AAIG_connectivities_dict=selected_AAIG_connectivities_dict,
                                                                   set_all_to_1=True)

        # deconvolute connectivities
        if calc_intersection:
            self.deconvolute_NOESY_NOESY_connectivities()

        # remove AAIG1->AAIG2 connections without CA atom, to reduce the number of connectivities
        if remove_matched_AAIGs_without_CA:
            self.remove_matched_AAIGs_without_CA(min_peak_num=min_peak_num)

        # populate self.AAIG_connectivities_complete_dict and self.i_iminus1_complete_dict (the are the same!)
        # TODO: control oexp parameter
        self.AAIG_connectivities_complete_dict = self.__get_NOESY_NOESY_connectivities_dict__()

    def deconvolute_NOESY_NOESY_connectivities(self):
        """
        Remove connectivities with very weak intensity intersection. DANGEROUS for small bandwidth values!!!
        :return:
        """
        for AAIG1_signature in list(self.NOESY_NOESY_2Dhist_intersections_mdict.keys()):
            # if len(self.NOESY_2Dhist_intersections_mdict[AAIG1_signature].values()) == 0:
            #     continue
            # max_intersection = np.max(self.NOESY_2Dhist_intersections_mdict[AAIG1_signature].values())
            for AAIG2 in list(self.NOESY_NOESY_2Dhist_intersections_mdict[AAIG1_signature].keys()):
                if self.NOESY_NOESY_2Dhist_intersections_mdict[AAIG1_signature][AAIG2] == 0: #/max_intersection < 0.00001:   # <== CHANGE ME
                    del self.NOESY_NOESY_2Dhist_intersections_mdict[AAIG1_signature][AAIG2]
                    del self.NOESY_NOESY_peak_connectivities_mdict[AAIG1_signature][AAIG2]  # AAIG1_signature, AAIG2 are names here
            if len(list(self.NOESY_NOESY_2Dhist_intersections_mdict[AAIG1_signature].keys())) == 0:    # remove the NOESY AAIG if no connectivities
                                                                                            # were left
                del self.NOESY_NOESY_2Dhist_intersections_mdict[AAIG1_signature]
                if AAIG1_signature in list(self.NOESY_NOESY_peak_connectivities_mdict.keys()):
                    del self.NOESY_NOESY_peak_connectivities_mdict[AAIG1_signature]


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

    def convert_NOESY_NOESY_peak_connectivities_mdict_values_to_lists(self):
        """
        Helper method to be called after match_NOESY_to_NOESY_Peaks(uniquify=False) in order to convert the values of
        self.NOESY_NOESY_peak_connectivities_mdict from Pandas DataFrames to lists of tuples (peak1, peak2, distance)
        :return:
        """
        NOESY_NOESY_peak_connectivities_mdict = \
            copy.deepcopy(self.NOESY_NOESY_peak_connectivities_mdict)
        for AAIG1_signature in list(NOESY_NOESY_peak_connectivities_mdict.keys()):
            for AAIG2_signature in list(NOESY_NOESY_peak_connectivities_mdict[AAIG1_signature].keys()):
                peakPeakDist_list = []
                df = self.NOESY_NOESY_peak_connectivities_mdict[AAIG1_signature][AAIG2_signature]
                df = df.sort_values(by="distance", ascending=True)
                for i, row in df.iterrows():
                    # print("DEBUG: row=", row)
                    peak1, peak2, distance = row
                    peakPeakDist_list.append( (peak1, peak2) )
                self.NOESY_NOESY_peak_connectivities_mdict[AAIG1_signature][AAIG2_signature] = peakPeakDist_list

    def uniquify_NOESY_NOESY_peak_connectivities(self, selected_AAIG_connectivities_dict=None):
        """
        ATTENTION: calling this method will convert the values of self.NOESY_NOESY_peak_connectivities_mdict from
        Pandas DataFrame to list of (peak,peak,distance) tuples! If you don't call it then self.NOESY_NOESY_peak_connectivities_mdict
        will contain Pandas DataFrames! This is importance for later processing steps!

        This method will find and save the maximum possible number of matching peaks between AAIG1 & AAIG2.
        NEW VERSION THAT WORKS WITH CHUNKS OF connections_dict TO REDUCE THE MEMORY OVERLOAD.

        There is the scenario where:
        p11 -> p21 = 0.2
        p12 -> p21 = 0.3
        p11 -> p22 = 0.4
        If we consider only the minimum distance, at end the method will keep only the peak pair p11 -> p21, but
        p12 and p22 will have no matches. However, we want to maximize the number of matched peaks in the whole
        spectrum, therefore we want this matching:  p11 -> p22, p12 -> p21. Now all 4 peaks are matched.

        Example:
        NOESY:
        ?-?-K60N-H		4.370	55.138	120.816	7.780	3514959.0
        ?-?-K60N-H		3.532	65.361	120.825	7.779	9851030.0
        ?-?-K60N-H		1.567	18.082	120.776	7.780	6530732.0
        ?-?-K60N-H		2.315	36.340	120.799	7.777	6136392.0
        ?-?-K60N-H		2.165	29.731	120.810	7.778	22834484.0
        ?-?-K60N-H		4.090	59.555	120.807	7.778	29782930.0
        ?-?-K60N-H		1.590	25.192	120.789	7.780	9520304.0
        ?-?-K60N-H		1.444	25.210	120.792	7.780	6372204.0
        ?-?-K60N-H		1.953	32.612	120.806	7.778	39343272.0
        ?-?-K60N-H		1.726	29.560	120.800	7.778	8433470.0
        ?-?-K60N-H		3.031	42.272	120.797	7.782	5141307.0
        TOCSY:
        K56HA-CA-I57N-H      4.106     59.644    122.370      8.069     19138414
        K56HB2-CB-I57N-H      2.045     32.276    122.373      8.068      9431224
        K56HB3-CB-I57N-H      1.899     32.312    122.357      8.069      9890725
        K56HG2-CG-I57N-H      1.697     25.610    122.368      8.066      5929990
        K56HG3-CG-I57N-H      1.522     25.639    122.354      8.071      5407076
        K56QD-CD-I57N-H      1.715     29.077    122.367      8.070      7114867
        the algorithm matches only 4 out of 6 peaks!

        :updates:   self.NOESY_NOESY_peak_connectivities_mdict: AAIG1_signature -> AAIG2_signature ->
                                                        [(AAIG1 Peak Object, matching AAIG2 Peak Object, C-H resonance distance), ...]
                    BUT IT REMOVES THE 'C-H resonance distance'!

        """
        # print("DEBUG: BEFORE UNIQUIFICATION self.NOESY_NOESY_peak_connectivities_mdict=")
        # print_dict_contents(self.NOESY_NOESY_peak_connectivities_mdict)

        small_connections_dict_args, big_connections_dict_args = [], [] # for aliphatic C-H
        small_AAIG1_AAIG2_args, big_AAIG1_AAIG2_args = [], []   # for aliphatic C-H
        for AAIG1_signature in list(self.NOESY_NOESY_peak_connectivities_mdict.keys()):
            for AAIG2_signature in list(self.NOESY_NOESY_peak_connectivities_mdict[AAIG1_signature].keys()):
                # The connectivities of AAIG2-AAIG1 are the same as AAIG1-AAIG2, but with the peak pairs reversed.
                # To save time, we skip them now and create them later. APPLICABLE ONLY IN HCNH-HCNH!
                if AAIG2_signature + "_" + AAIG1_signature in small_AAIG1_AAIG2_args + big_AAIG1_AAIG2_args:
                    continue

                connections_dict = {}  # NOESY aliphatic AAIG1 peak -> (NOESY aliphatic AAIG2 peak, C-H resonance distance)
                df = self.NOESY_NOESY_peak_connectivities_mdict[AAIG1_signature][AAIG2_signature]
                # distance: is the C,H resonance distance, not 2Dhist intersection!
                df = df.sort_values(by="distance", ascending=True)  # sort by ascending distance???
                # print("DEBUG: %s %s df=%s" % (AAIG1_signature, AAIG2_signature, df))
                for i, row in df.iterrows():
                    # print("DEBUG: row=", row)
                    peak1, peak2, distance = row
                    try:
                        connections_dict[peak1].append((peak2, distance))
                    except KeyError:
                        connections_dict[peak1] = [(peak2, distance)]

                # Chunk the connections_dict into smaller parts to reduce the memory overhead
                for connections_chunk_dict in Connectivities.chunk_connections_dict(connections_dict, tolC=0.2):
                    if 4 >= len(connections_chunk_dict.keys()) > 0:  # AAIG1 & AAIG2 share <= 4 common aliphatic peaks
                        small_connections_dict_args.append(connections_chunk_dict)
                        small_AAIG1_AAIG2_args.append(AAIG1_signature + "_" + AAIG2_signature)
                    elif len(connections_chunk_dict.keys()) > 4:  # AAIG1 & AAIG2 share more than 4 common aliphatic peaks
                        big_connections_dict_args.append(connections_chunk_dict)
                        big_AAIG1_AAIG2_args.append(AAIG1_signature + "_" + AAIG2_signature)

        ColorPrint("CALCULATING PEAK CONNECTIVITIES", "BOLDGREEN")
        # Parallel execution of small aliphatic connections
        # print("DEBUG: small_connections_dict_args=")
        # for args in small_connections_dict_args:
        #     print_dict_contents(args)
        #     print("\n")
        small_results = list(futures.map(self.get_best_NOESY_connectivity,
                                         small_connections_dict_args,
                                         small_AAIG1_AAIG2_args))
        # print("DEBUG: small_results=")
        # for r in small_results:
        #     print(r[0], r[1], [(unfold_peak(d[0]), unfold_peak(d[1]))for d in r[2]])
        # Parallel execution of big aliphatic connections
        big_results = []
        for connections_dict, AAIG1_AAIG2 in zip(big_connections_dict_args, big_AAIG1_AAIG2_args):
            big_results.append(self.get_best_NOESY_connectivity(connections_dict, AAIG1_AAIG2, parallel=True)
                               )

        ColorPrint("Combining the individual connectivity results.", "OKGREEN")
        # Because self.NOESY_NOESY_peak_connectivities_mdict is not empty, but __UPDATED__, we need to clean it up first
        all_AAIG1_AAIG2_set = set(small_AAIG1_AAIG2_args + big_AAIG1_AAIG2_args)
        for AAIG1_AAIG2 in all_AAIG1_AAIG2_set:
            AAIG1_signature, AAIG2_signature = AAIG1_AAIG2.split('_')
            self.NOESY_NOESY_peak_connectivities_mdict[AAIG1_signature][AAIG2_signature] = None
            self.NOESY_NOESY_peak_connectivities_mdict[AAIG2_signature][AAIG1_signature] = None # clean the reverse, too

        assert len(small_results + big_results) == len(small_AAIG1_AAIG2_args + big_AAIG1_AAIG2_args), \
            Debuginfo("ERROR: not all connections_chunk_dict produced results!", fail=True)

        # best_chain: the optimum list of [(AAIG1 Peak object, matching AAIG2 Peak object), ...] (without the distance)
        for AAIG1_signature, AAIG2_signature, best_chain in small_results + big_results:
            # print("DEBUG: AAIG1=%s AAIG2=%s best_chain=%s" % (AAIG1_signature, AAIG2_signature, best_chain))
            # for pp in best_chain:
            #     print((unfold_peak(pp[0]), unfold_peak(pp[1])))
            try:
                self.NOESY_NOESY_peak_connectivities_mdict[AAIG1_signature][AAIG2_signature].extend(best_chain)
            except AttributeError:
                self.NOESY_NOESY_peak_connectivities_mdict[AAIG1_signature][AAIG2_signature] = list(best_chain)

        # Copy the connectivities of AAIG1-AAIG2 to AAIG2-AAIG1, but reverse the peak pairs in the chains.
        # APPLICABLE ONLY IN HCNH-HCNH!
        for AAIG1_signature, AAIG2_signature, best_chain in small_results + big_results:
            self.NOESY_NOESY_peak_connectivities_mdict[AAIG2_signature][AAIG1_signature] = \
                [tuple(reversed(pp)) for pp in self.NOESY_NOESY_peak_connectivities_mdict[AAIG1_signature][AAIG2_signature]]

        # print("DEBUG: AFTER UNIQUIFICATION self.NOESY_NOESY_peak_connectivities_mdict=")
        # print(self.NOESY_NOESY_peak_connectivities_mdict)
        # print_dict_contents(self.NOESY_NOESY_peak_connectivities_mdict)

    def remove_matched_AAIGs_without_CA(self, OFFSET=0.4999, min_peak_num=1):
        """
        Method to remove AAIG1->AAIG2 connections that have common peaks (any number) but none of them
        corresponds to Ca atom.
        :param OFFSET: determines the borders of the probability twilight zone, within the answer of the CA-HA classifier is
                        uncertain.
        :param min_peak_num: remove only connectivities between AAIGs that have at least 'min_peak_num' matching peaks
                                the most. The lower the more connectivities it removes?
        :return:
        """

        def within_twilight(prob_series, OFFSET):
            """
            Returns True if any of the probabilities is within the twilight zone defined by OFFSET.
            :param prob_series:
            :param OFFSET:
            :return:
            """
            answers = []
            for prob in prob_series:
                if (prob < 0.5 - OFFSET) or (prob > 0.5 + OFFSET):
                    answers.append(False)
                else:
                    answers.append(True)
            return any(answers)

        PredCA = pickle.load(open(CHAINS_LIB_DIR[:-3] + "ML_models_DELETE/ca_rf.pickle", "rb"))
        NOESY_NOESY_peak_connectivities_mdict = \
            copy.deepcopy(self.NOESY_NOESY_peak_connectivities_mdict)
        for AAIG1_signature in list(NOESY_NOESY_peak_connectivities_mdict.keys()):
            for AAIG2_signature in list(NOESY_NOESY_peak_connectivities_mdict[AAIG1_signature].keys()):
                # Combine it with an occupancy threshold, otherwise it will also remove correct connections.
                if len(NOESY_NOESY_peak_connectivities_mdict[AAIG1_signature][AAIG2_signature]) < min_peak_num:    # only connectivities with 1 common peak will be considered
                    continue

                found_CA = False
                for peak1, peak2 in NOESY_NOESY_peak_connectivities_mdict[AAIG1_signature][AAIG2_signature]:
                    df = pd.DataFrame([(peak1.Creson, peak1.Hreson), (peak2.Creson, peak2.Hreson)],
                                  columns=["Creson", "Hreson"])
                    preds = PredCA.predict(df)

                    if within_twilight(preds[True], OFFSET):      # these are uncertain peaks, therefore consider them as CA-HA
                        found_CA = True                           # and keep the connectivity
                        break

                    if (preds[True] >= 0.5).any():              # if any of peak1 and peak2 is CA-HA and none is in twilight zone,
                                                                # keep the connectivity
                        found_CA = True
                        break
                if not found_CA:
                    ColorPrint("Removing connectivity %s -> %s due to lack of matched CA-HA peak." %
                               (AAIG1_signature, AAIG2_signature), "OKBLUE")
                    del self.NOESY_NOESY_peak_connectivities_mdict[AAIG1_signature][AAIG2_signature]
                    del self.NOESY_NOESY_2Dhist_intersections_mdict[AAIG1_signature][AAIG2_signature]

            # Delete AAIGs without connectivities
            if len(self.NOESY_NOESY_2Dhist_intersections_mdict[AAIG1_signature]) == 0:
                ColorPrint("WARNING: none of the connections of %s contained CA-HA peak and thus this AAIG"
                           " has no connectivities!" % AAIG1_signature, "WARNING")
                del self.NOESY_NOESY_2Dhist_intersections_mdict[AAIG1_signature]
                del self.NOESY_NOESY_peak_connectivities_mdict[AAIG1_signature]

    def find_recurrent_peak(self, cutoff=1):
        """
        Method to iterate over all AAIG(i)->AAIG(?) connectivities and find those AAIG(i) peaks that tend to match
        most frequently with peaks from other AAIGs, aka "recurrent peaks". These are more likely to correspond to
        atoms of AAIG(i).
        :param cutoff: minimum number of peaks of AAIGs other that (i) that a peak of AAIG(i) must match in
                        order to be considered as recurrent.
        :return: it just populates self.AAIG_recurrent_peaks_mdict.
        """
        ColorPrint("Searching for recurrent peaks in the AAIGs of NOESY spectrum.", "BOLDGREEN")
        for AAIG1_signature in self.AAIG_connectivities_complete_dict.keys():
            peakID_frequency_dict = defaultdict(int)    # dict of zeros
            for quintuplet in self.AAIG_connectivities_complete_dict[AAIG1_signature]:
                AAIG2_signature = quintuplet[0]
                for peak1, peak2 in self.NOESY_NOESY_peak_connectivities_mdict[AAIG1_signature][AAIG2_signature]:
                    peakID_frequency_dict[peak1.ID] += 1
            for peakID, freq in peakID_frequency_dict.items():
                if freq >= cutoff:
                    aaig = self.HCNH_spec.__get_AAIG__(AAIG1_signature)
                    peak = aaig.__get_peak_by_ID__(peakID)
                    self.AAIG_recurrent_peaks_mdict[AAIG1_signature][peak] = freq

    # Import all the sub-modules
    from ._write_methods import write_NOESY_NOESY_connectivities, write_NOESY_NOESY_connectivities_from_dict, write_spectrum_info
    from ._getter_methods import __get_NOESY_NOESY_connectivities_dict__, __get_filtered_NOESY_NOESY_connectivities_dict__, __get_NOESY_NOESY_connectivities_dataframe__, __get_NOESY_NOESY_peak_connectivities_mdict__, __get_NOESY_NOESY_common_peak_dataframe__, _scale_multiscores
    from ._proofread import proofread_connectivities



