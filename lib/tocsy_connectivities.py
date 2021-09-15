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


import os
from operator import itemgetter

from lib.aaig import AAIG
from lib.connectivities import Connectivities, Connectivities_Tree, HCNH_Spectrum
from lib.global_func import tree, bcolors, ColorPrint, Debuginfo, save_pickle, get_NH_name, dash_to_NH
from lib.tocsy_spectrum import TOCSY_Spectrum
import pandas as pd
import numpy as np
from scoop import futures, shared


class TOCSY_Connectivities(Connectivities):
    """
    This Class is meant to find connectivities between TOCSY AAIGs and NOESY AAIGs.
    """

    def __init__(self,
                 tolH=0.04,
                 tolC=0.4,
                 H_bin_length=0.0005,
                 C_bin_length=0.005,
                 bandwidth=0.2):

        super(TOCSY_Connectivities, self).__init__()

        self.TOCSY_NOESY_peak_connectivities_mdict = tree()  # TOCSY AAIG(i) -> matching NOESY AAIG(i-1) ->
        # [(matching TOCSY peak1, matching NOESY peak2, distance), ...]
        self.TOCSY_NOESY_2Dhist_intersections_mdict = tree()  # TOCSY AAIG(i) -> matching NOESY AAIG(i-1) -> 2D-hist intersection
        # Make TOCSY-relevant copies for compatibility with the old code. They will be automatically updated.
        self.i_iminus1_complete_dict = self.AAIG_connectivities_complete_dict

    def clean_connectivities(self):
        print(bcolors.BOLDGREEN + "Cleaning all saved connectivities." + bcolors.ENDBOLD)
        self.TOCSY_NOESY_peak_connectivities_mdict = tree()  # TOCSY AAIG(i) -> matching NOESY AAIG(i-1) ->
        # [(matching TOCSY peak1, matching NOESY peak2, distance), ...]
        # The following clean self.AAIG_connectivities_complete_dict, too.
        self.i_iminus1_complete_dict = {}  # TOCSY AAIG(i) -> [(NOESY AAIG(i-1), occupancy, tot. num. peaks, intersection, multi score),
        # ...]

    def clean_intersections(self):
        print(bcolors.BOLDGREEN + "Cleaning all saved intersections." + bcolors.ENDBOLD)
        self.TOCSY_NOESY_2Dhist_intersections_mdict = tree()  # TOCSY AAIG(i) -> matching NOESY AAIG(i-1) -> 2D-hist intersection

    def load_TOCSY_spectrum(self, TOCSY_fname):
        self.TOCSY_spec = TOCSY_Spectrum(TOCSY_fname=TOCSY_fname,
                                         H_bin_length=self.H_bin_length,
                                         C_bin_length=self.C_bin_length,
                                         bandwidth=self.bandwidth)
        # if the NOESY spectrum has been also loaded, check if it is identical to the TOCSY spectrum
        if self.HCNH_spec and self.TOCSY_spec.__is_equal__(self.HCNH_spec):
            raise Exception(
                bcolors.FAIL + "ERROR: TOCSY and NOESY spectra are identical! Check your input files!" + bcolors.ENDC)

    def write_spectrum_info(self):
        if not self.TOCSY_spec or not self.HCNH_spec:
            raise Exception(bcolors.FAIL + "ERROR: you must load both TOCSY and NOESY spectra in order to write statistics "
                                           "in file!" + bcolors.ENDC)
        with open("spectrum_statistics.txt", 'w') as f:
            AAIG_num, peak_num, ave_peaknum_per_aaig, min_peaknum_per_aaig, max_peaknum_per_aaig = \
                self.TOCSY_spec.__get_spectrum_info__()
            f.write("\nTOCSY:\n")
            f.write("AAIG number: %g\n" % AAIG_num)
            f.write("Peak number: %g\n" % peak_num)
            f.write("Average peak number per AAIG: %f\n" % ave_peaknum_per_aaig)
            f.write("Minimum peak number per AAIG: %i\n" % min_peaknum_per_aaig)
            f.write("Maximum peak number per AAIG: %i\n" % max_peaknum_per_aaig)
            AAIG_num, peak_num, ave_peaknum_per_aaig, min_peaknum_per_aaig, max_peaknum_per_aaig = \
                self.HCNH_spec.__get_spectrum_info__()
            f.write("\nNOESY:\n")
            f.write("AAIG number: %g\n" % AAIG_num)
            f.write("Peak number: %g\n" % peak_num)
            f.write("Average peak number per AAIG: %f\n" % ave_peaknum_per_aaig)
            f.write("Minimum peak number per AAIG: %i\n" % min_peaknum_per_aaig)
            f.write("Maximum peak number per AAIG: %i\n" % max_peaknum_per_aaig)

    def save_TOCSY_NOESY_peak_connectivity(self, AAIG1_signature, peak1, AAIG2_signature, peak2):
        """
        Method to save all possible connectivities and their distances into dataframes. You must call
        uniquify_TOCSY_peak_connectivities() afterwards to remove redundancy!

        :param AAIG1_signature:  the TOCSY AAIG name (residue i)
        :param peak1:
        :param AAIG2_signature:  the NOESY AAIG name (putative residue i-1)
        :param peak2:
        :return:
        """
        distance = peak1.distance(peak2, attypes=['C', 'H'])
        try:
            df = self.TOCSY_NOESY_peak_connectivities_mdict[AAIG1_signature][AAIG2_signature]
            df = df.append({'TOCSY_peak': peak1, 'NOESY_peak': peak2, 'distance': distance}, ignore_index=True)
            self.TOCSY_NOESY_peak_connectivities_mdict[AAIG1_signature][AAIG2_signature] = df  # save changes
        except (AttributeError, KeyError):
            self.TOCSY_NOESY_peak_connectivities_mdict[AAIG1_signature][AAIG2_signature] = \
                pd.DataFrame([(peak1, peak2, distance)], columns=("TOCSY_peak", "NOESY_peak", "distance" ))

    def get_best_TOCSY_connectivity(self, connections_dict, AAIG1_AAIG2, parallel=False):
        """
        Method to find the maximum possible number of matching peaks between AAIG1 and AAIG2.
        :param parallel:    highly parallel more (worthy only if AAIG1 has more that 4 matching peaks with AAIG2)
        :return:    best_chain: the optimum list of [(AAIG1 Peak object, matching AAIG2 Peak object), ...] (without the distance)
        """
        AAIG1_signature, AAIG2_signature = AAIG1_AAIG2.split("_")
        con_tree = Connectivities_Tree(connections_dict)
        if parallel:
            print(bcolors.OKBLUE + "Uniquifying big connectivities between NOESY AAIG " + AAIG1_signature + \
                  " and TOCSY AAIG " + AAIG2_signature + bcolors.ENDC)
            con_tree.form_chains_parallel(length='max')
        else:
            print(bcolors.OKBLUE + "Uniquifying small connectivities between NOESY AAIG " + AAIG1_signature + \
                  " and TOCSY AAIG " + AAIG2_signature + bcolors.ENDC)
            con_tree.form_chains(length='max')
        return AAIG1_signature, AAIG2_signature, con_tree.get_best_chain(length='max')

    # TEMPORARILY DEACTIVATED TO TEST THE NEW VERSION THAT WORKS WITH CHUNKS OF connections_dict
    # def uniquify_TOCSY_NOESY_peak_connectivities(self, selected_i_iminus1_dict=None):
    #     """
    #     This method will find and save the maximum possible number of matching peaks between AAIG1 & AAIG2.
    #
    #     There is the scenario where:
    #     Tp1 -> Np1 = 0.2
    #     Tp2 -> Np1 = 0.3
    #     Tp1 -> Np2 = 0.4
    #     So at the end the method will keep only:
    #     Tp1 -> Np1
    #     Tp2 will have no match!
    #     Example:
    #     NOESY:
    #     ?-?-K60N-H		4.370	55.138	120.816	7.780	3514959.0
    #     ?-?-K60N-H		3.532	65.361	120.825	7.779	9851030.0
    #     ?-?-K60N-H		1.567	18.082	120.776	7.780	6530732.0
    #     ?-?-K60N-H		2.315	36.340	120.799	7.777	6136392.0
    #     ?-?-K60N-H		2.165	29.731	120.810	7.778	22834484.0
    #     ?-?-K60N-H		4.090	59.555	120.807	7.778	29782930.0
    #     ?-?-K60N-H		1.590	25.192	120.789	7.780	9520304.0
    #     ?-?-K60N-H		1.444	25.210	120.792	7.780	6372204.0
    #     ?-?-K60N-H		1.953	32.612	120.806	7.778	39343272.0
    #     ?-?-K60N-H		1.726	29.560	120.800	7.778	8433470.0
    #     ?-?-K60N-H		3.031	42.272	120.797	7.782	5141307.0
    #     TOCSY:
    #     K56HA-CA-I57N-H      4.106     59.644    122.370      8.069     19138414
    #     K56HB2-CB-I57N-H      2.045     32.276    122.373      8.068      9431224
    #     K56HB3-CB-I57N-H      1.899     32.312    122.357      8.069      9890725
    #     K56HG2-CG-I57N-H      1.697     25.610    122.368      8.066      5929990
    #     K56HG3-CG-I57N-H      1.522     25.639    122.354      8.071      5407076
    #     K56QD-CD-I57N-H      1.715     29.077    122.367      8.070      7114867
    #     the algorithm matches only 4 out of 6 peaks!
    #
    #     :updates:   self.TOCSY_NOESY_peak_connectivities_mdict: AAIG1_signature -> AAIG2_signature ->
    #                                                                 [(AAIG1 Peak Object, matching AAIG2 Peak Object), ...]
    #
    #     """
    #
    #     small_connections_dict_args, big_connections_dict_args = [], []
    #     small_AAIG1_AAIG2_args, big_AAIG1_AAIG2_args = [], []
    #     for AAIG1_signature in list(self.TOCSY_NOESY_peak_connectivities_mdict.keys()):
    #         for AAIG2_signature in list(self.TOCSY_NOESY_peak_connectivities_mdict[AAIG1_signature].keys()):
    #             # The connectivities of AAIG2-AAIG1 are the same as AAIG1-AAIG2, but with the peak pairs reversed.
    #             # To save time, we skip them now and create them later.
    #             if AAIG2_signature + "_" + AAIG1_signature in small_AAIG1_AAIG2_args + big_AAIG1_AAIG2_args:
    #                 continue
    #             connections_dict = {}   # TOCSY peak -> (NOESY peak, distance)
    #             df = self.TOCSY_NOESY_peak_connectivities_mdict[AAIG1_signature][AAIG2_signature]
    #             df = df.sort_values(by="distance", ascending=False)  # sort by ascending distance
    #             for i, row in df.iterrows():
    #                 peak1, peak2, distance = row
    #                 try:
    #                     connections_dict[peak1].append((peak2, distance))
    #                 except KeyError:
    #                     connections_dict[peak1] = [(peak2, distance)]
    #             if 4 >= len(connections_dict) > 0:
    #                 small_connections_dict_args.append(connections_dict)
    #                 small_AAIG1_AAIG2_args.append(AAIG1_signature + "_" + AAIG2_signature)
    #             elif len(connections_dict) > 4: # AAIG1 & AAIG2 share more that 4 common peaks
    #                 big_connections_dict_args.append(connections_dict)
    #                 big_AAIG1_AAIG2_args.append(AAIG1_signature + "_" + AAIG2_signature)
    #     # Parallel execution of small connections
    #     small_results = list(futures.map(self.get_best_TOCSY_connectivity,
    #                                      small_connections_dict_args,
    #                                      small_AAIG1_AAIG2_args))
    #     # Parallel execution of big connections
    #     big_results = []
    #     for connections_dict, AAIG1_AAIG2 in zip(big_connections_dict_args, big_AAIG1_AAIG2_args):
    #         big_results.append(self.get_best_TOCSY_connectivity(
    #                                         connections_dict,
    #                                         AAIG1_AAIG2,
    #                                         parallel=True)
    #                            )
    #     for AAIG1_signature, AAIG2_signature, best_chain in small_results + big_results:
    #         self.TOCSY_NOESY_peak_connectivities_mdict[AAIG1_signature][AAIG2_signature] = best_chain
    #
    #     # Copy the connectivities of AAIG1-AAIG2 to AAIG2-AAIG1, but reverse the peak pairs in the chains.
    #     for AAIG1_signature, AAIG2_signature, best_chain in small_results + big_results:
    #         self.TOCSY_NOESY_peak_connectivities_mdict[AAIG2_signature][AAIG1_signature] = \
    #             [tuple(reversed(pp)) for pp in
    #              self.TOCSY_NOESY_peak_connectivities_mdict[AAIG1_signature][AAIG2_signature]]
    
    def uniquify_TOCSY_NOESY_peak_connectivities(self, selected_i_iminus1_dict=None):
        """
        This method will find and save the maximum possible number of matching peaks between AAIG1 & AAIG2.
        NEW VERSIONS THAT WORKS WITH CHUNKS OF connections_dict TO REDUCE THE MEMORY OVERLOAD.

        There is the scenario where:
        Tp1 -> Np1 = 0.2
        Tp2 -> Np1 = 0.3
        Tp1 -> Np2 = 0.4
        So at the end the method will keep only:
        Tp1 -> Np1
        Tp2 will have no match!
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

        :updates:   self.TOCSY_NOESY_peak_connectivities_mdict: AAIG1_signature -> AAIG2_signature ->
                                                                    [(AAIG1 Peak Object, matching AAIG2 Peak Object), ...]

        """

        small_connections_dict_args, big_connections_dict_args = [], [] # for aliphatic C-H
        small_AAIG1_AAIG2_args, big_AAIG1_AAIG2_args = [], []   # for aliphatic C-H
        for AAIG1_signature in list(self.TOCSY_NOESY_peak_connectivities_mdict.keys()):
            for AAIG2_signature in list(self.TOCSY_NOESY_peak_connectivities_mdict[AAIG1_signature].keys()):
                connections_dict = {}  # TOCSY aliphatic AAIG1 peak -> (NOESY aliphatic AAIG2 peak, distance)
                df = self.TOCSY_NOESY_peak_connectivities_mdict[AAIG1_signature][AAIG2_signature]
                df = df.sort_values(by="distance", ascending=False)  # sort by ascending distance???
                # print("DEBUG: %s %s df=%s" % (AAIG1_signature, AAIG2_signature, df))
                for i, row in df.iterrows():
                    peak1, peak2, distance = row
                    try:
                        connections_dict[peak1].append((peak2, distance))
                    except KeyError:
                        connections_dict[peak1] = [(peak2, distance)]

                # Chunk the connections_dict into smaller parts to reduce the memory consumption
                for connections_chunk_dict in Connectivities.chunk_connections_dict(connections_dict, tolC=0.2):
                    if 4 >= len(connections_chunk_dict) > 0:  # AAIG1 & AAIG2 share <= 4 common aliphatic peaks
                        small_connections_dict_args.append(connections_chunk_dict)
                        small_AAIG1_AAIG2_args.append(AAIG1_signature + "_" + AAIG2_signature)
                    elif len(connections_chunk_dict) > 4:  # AAIG1 & AAIG2 share more than 4 common aliphatic peaks
                        big_connections_dict_args.append(connections_chunk_dict)
                        big_AAIG1_AAIG2_args.append(AAIG1_signature + "_" + AAIG2_signature)

        ColorPrint("CALCULATING PEAK CONNECTIVITIES", "BOLDGREEN")
        # Parallel execution of small aliphatic connections
        small_results = list(futures.map(self.get_best_TOCSY_connectivity,
                                         small_connections_dict_args,
                                         small_AAIG1_AAIG2_args))
        # Parallel execution of big aliphatic connections
        big_results = []
        for connections_dict, AAIG1_AAIG2 in zip(big_connections_dict_args, big_AAIG1_AAIG2_args):
            big_results.append(self.get_best_TOCSY_connectivity(connections_dict, AAIG1_AAIG2, parallel=True)
                               )

        ColorPrint("Combining the individual connectivity results.", "OKGREEN")
        # Because self.TOCSY_NOESY_peak_connectivities_mdict is not empty, but __UPDATED__, we need to clean it up first
        all_AAIG1_AAIG2_set = set(small_AAIG1_AAIG2_args + big_AAIG1_AAIG2_args)
        for AAIG1_AAIG2 in all_AAIG1_AAIG2_set:
            AAIG1_signature, AAIG2_signature = AAIG1_AAIG2.split('_')
            self.TOCSY_NOESY_peak_connectivities_mdict[AAIG1_signature][AAIG2_signature] = None

        assert len(small_results + big_results) == len(small_AAIG1_AAIG2_args + big_AAIG1_AAIG2_args), \
            Debuginfo("ERROR: not all connections_chunk_dict produced results!", fail=True)

        # best_chain: the optimum list of [(AAIG1 Peak object, matching AAIG2 Peak object), ...] (without the distance)
        for AAIG1_signature, AAIG2_signature, best_chain in small_results + big_results:
            try:
                self.TOCSY_NOESY_peak_connectivities_mdict[AAIG1_signature][AAIG2_signature].extend(best_chain)
            except AttributeError:
                self.TOCSY_NOESY_peak_connectivities_mdict[AAIG1_signature][AAIG2_signature] = list(best_chain)

    
    def deconvolute_TOCSY_NOESY_connectivities(self):
        """
        Remove connectivities with very weak intensity intersection. DANGEROUS for small bandwidth values!!!
        :return:
        """
        for AAIG1 in list(self.TOCSY_NOESY_2Dhist_intersections_mdict.keys()):
            # if len(self.TOCSY_NOESY_2Dhist_intersections_mdict[AAIG1].values()) == 0:
            #     continue
            # max_intersection = np.max(self.TOCSY_NOESY_2Dhist_intersections_mdict[AAIG1].values())
            for AAIG2 in list(self.TOCSY_NOESY_2Dhist_intersections_mdict[AAIG1].keys()):
                if self.TOCSY_NOESY_2Dhist_intersections_mdict[AAIG1][AAIG2] == 0: #/max_intersection < 0.00001:   # <== CHANGE ME
                    del self.TOCSY_NOESY_2Dhist_intersections_mdict[AAIG1][AAIG2]
                    del self.TOCSY_NOESY_peak_connectivities_mdict[AAIG1][AAIG2]  # AAIG1, AAIG2 are names here
            if len(list(self.TOCSY_NOESY_2Dhist_intersections_mdict[AAIG1].keys())) == 0:    # remove the TOCSY AAIG if no connectivities
                                                                                            # were left
                del self.TOCSY_NOESY_2Dhist_intersections_mdict[AAIG1]
                if AAIG1 in list(self.TOCSY_NOESY_peak_connectivities_mdict.keys()):
                    del self.TOCSY_NOESY_peak_connectivities_mdict[AAIG1]

    def save_TOCSY_NOESY_2Dhist_AAIG_connectivity(self, AAIG1, AAIG2, use_Peak_2Dhistograms=False):
        """
        Method to save the intersection between AAIG1 & AAIG2 expressed as 2D density histograms.

        :param AAIG1: TOCSY AAIG (has the C,H or residue i-1)
        :param AAIG2: NOESY AAIG (putative residue i-1)
        :param use_Peak_2Dhistograms: use a separate 2D-histogram for each peak of AAIG1 & AAIG2 (NOT FUNCTIONAL!)
        :return:
        """
        if use_Peak_2Dhistograms:
            pass
        else:
            intersection = AAIG1.calc_intersection(AAIG2)  # it applies only to C,H resonances

        self.TOCSY_NOESY_2Dhist_intersections_mdict[AAIG1.signature][AAIG2.signature] = intersection

    def __scale_intensities__(self, wrt_aaig=False):
        """
        Method used only by TOCSY- connectivities.
        :param wrt_aaig:    if True scale the intensities wrt the AAIG, otherwise wrt the whole spectrum maximum.
        :return:
        """
        self.HCNH_spec.__scale_intensities__(even_intensities=False, wrt_aaig=wrt_aaig)
        self.TOCSY_spec.__scale_intensities__(even_intensities=True, wrt_aaig=wrt_aaig)

    def scale_TOCSY_NOESY_2Dhist_connectivities(self, AAIG1_signature):
        """
        Method to scale the intersections of the connectivities of TOCSY AAIG1 to be between 0 and 1.0.
        :param AAIG1:
        :return:
        """
        if len(list(self.TOCSY_NOESY_2Dhist_intersections_mdict[AAIG1_signature].values())) == 0: # if no connectivities were found
            return
        max_intersection = np.max(list(self.TOCSY_NOESY_2Dhist_intersections_mdict[AAIG1_signature].values()))
        if max_intersection == 0:   # if all connectivities have intersection 0, leave them as they are (they will be removed)
            return
        for AAIG2_signature in list(self.TOCSY_NOESY_2Dhist_intersections_mdict[AAIG1_signature].keys()):
            self.TOCSY_NOESY_2Dhist_intersections_mdict[AAIG1_signature][AAIG2_signature] /= max_intersection

    def TOCSY_AAIG2hist(self, AAIG_signature):
        """
       Create the 2D-histogram of this TOCSY AAIG.
       :param AAIG_signature:
       :return:
       """
        self.TOCSY_spec.__AAIG2hist__(AAIG_signature, even_intensities=True, scale_density=False, **self.hist_borders)
        return self.TOCSY_spec.__get_AAIG__(AAIG_signature)

    def TOCSY_Peaks2hist(self, AAIG_signature):
        """
        Convenience method for parallel execution, creates the 2D-histogram of all individual peaks of this TOCSY AAIG.
       :param AAIG_signature:
       :return:
       """
        even_edges = shared.getConst('EVEN_EDGES')
        scale_density = shared.getConst('SCALE_DENSITY')
        self.TOCSY_spec.__AAIG_Peaks2hist__(AAIG_signature, even_edges=even_edges, even_intensities=True,
                                            scale_density=scale_density)
        return self.TOCSY_spec.__get_AAIG__(AAIG_signature)

    # def match_TOCSY_AAIG_against_all_NOESY_AAIGs(self,
    #                                              AAIG1,
    #                                              tolH=0.04,
    #                                              tolC=0.4,
    #                                              calc_intersection=True,
    #                                              selected_i_iminus1_dict=None):
    #     """
    #     NOT USED! USE INSTEAD match_TOCSY_AAIG_against_all_NOESY_Peaks()
    #     :param AAIG1: TOCSY AAIG, contains C,H resonances of i-1 and N,H of i, and has the name of residue i.
    #     :param tolH:
    #     :param tolC:
    #     :return:
    #     """
    #     AAIG1_signature = AAIG1.signature
    #     if selected_i_iminus1_dict and not AAIG1_signature in selected_i_iminus1_dict:
    #         return
    #     print(bcolors.OKBLUE + "Searching for the TCOSY-NOESY connectivities of AAIG " + AAIG1_signature + bcolors.ENDC)
    #     for peak1 in AAIG1.__get_all_peaks__():  # the C,H come from i-1 but the N,HN resonances from i
    #         for AAIG2_signature in self.noesy_spec.__get_all_AAIG_signatures__():
    #             if AAIG1_signature == AAIG2_signature:  continue  # recall that AAIG1_signature is i but AAIG1 contains the C,H of i-1!
    #             if selected_i_iminus1_dict and not AAIG2_signature in [q[0] for q in selected_i_iminus1_dict[AAIG1_signature]]:
    #                 continue
    #             AAIG2 = self.noesy_spec.__get_AAIG__(AAIG2_signature)
    #             for peak2 in AAIG2.__get_all_peaks__():
    #                 if peak1.does_peak_match(peak2, tolH=tolH, tolC=tolC, attypes=['C', 'H']):
    #                     self.save_TOCSY_NOESY_peak_connectivity(AAIG1_signature, peak1, AAIG2_signature,
    #                                                             peak2)  # store this connectivity
    #                     # if not already present, calculate and save the 2D-hist intersection between the two AAIGs
    #                     if not AAIG1_signature in list(self.TOCSY_NOESY_2Dhist_intersections_mdict.keys()) or \
    #                             not AAIG2_signature in list(self.TOCSY_NOESY_2Dhist_intersections_mdict[AAIG1.signature].keys()):
    #                         if calc_intersection:   # save AAIG1-AAIG2 2D density histogram intersection
    #                             self.save_TOCSY_NOESY_2Dhist_AAIG_connectivity(AAIG1, AAIG2)
    #                         else: # if you don't want to calculate intersection, save the default value 1.0
    #                             self.TOCSY_NOESY_2Dhist_intersections_mdict[AAIG1.signature][AAIG2.signature] = 1.0  # <== CHANGE ME
    #     self.scale_TOCSY_NOESY_2Dhist_connectivities(AAIG1_signature)  # scale all intersections for this TOCSY AAIG to be between 0 and 1

    def match_TOCSY_AAIG_against_all_NOESY_Peaks(self, AAIG1, tolH=0.04, tolC=0.4, selected_i_iminus1_dict=None):
        """

        :param AAIG1: TOCSY AAIG, contains C,H resonances of i-1 and N,H of i, and has the name of residue i.
        :param tolH:
        :param tolC:
        :return:
        """
        AAIG1_signature = AAIG1.signature
        if selected_i_iminus1_dict and not AAIG1_signature in selected_i_iminus1_dict:
            return
        print(bcolors.OKBLUE + "Searching for the TCOSY-NOESY connectivities of AAIG " + AAIG1_signature + bcolors.ENDC)
        for peak1 in AAIG1.__get_all_peaks__():  # the C,H come from i-1 but the N,HN resonances from i
            for AAIG2_signature in self.HCNH_spec.__get_all_AAIG_signatures__():
                if AAIG1_signature == AAIG2_signature:  continue  # recall that AAIG1_signature is i but AAIG1 contains the C,H of i-1!
                if selected_i_iminus1_dict and not AAIG2_signature in [q[0] for q in selected_i_iminus1_dict[AAIG1_signature]]:
                    continue
                AAIG2 = self.HCNH_spec.__get_AAIG__(AAIG2_signature)
                for peak2 in AAIG2.__get_all_peaks__():
                    if peak1.does_peak_match(peak2, tolH=tolH, tolC=tolC, attypes=['C', 'H']):
                        self.save_TOCSY_NOESY_peak_connectivity(AAIG1_signature, peak1, AAIG2_signature, peak2)  # store this connectivity

    def measure_all_TOCSY_NOESY_Peak_2Dhist_intersections(self, tolH=0.04, tolC=0.4, selected_i_iminus1_dict=None, set_all_to_1=False):
        """
        For each TOCSY AAIG1, look at the Peak connectivities and calculate the overall 2D-histogram
        intersection for each connectivity.

        :param tolH:
        :param tolC:
        :return:
        """
        # save_pickle("TOCSY_NOESY_peak_connectivities_mdict.pkl", self.TOCSY_NOESY_peak_connectivities_mdict)  # for DEBUGGING
        # save_pickle("TOCSY_spec.pkl", self.TOCSY_spec)    # FOR DEBUGGING
        # save_pickle("noesy_spec.pkl", self.noesy_spec)    # FOR DEBUGGING
        for AAIG1_signature in list(self.TOCSY_NOESY_peak_connectivities_mdict.keys()):
            if selected_i_iminus1_dict and not AAIG1_signature in selected_i_iminus1_dict:
                return
            AAIG1 = self.TOCSY_spec.__get_AAIG__(AAIG1_signature)
            print(bcolors.OKBLUE + "Measuring Peak-Peak 2D-histogram intersection for the connectivities of TOCSY AAIG " \
                  + AAIG1_signature + bcolors.ENDC)
            for AAIG2_signature in list(self.TOCSY_NOESY_peak_connectivities_mdict[AAIG1_signature].keys()):
                if AAIG1_signature == AAIG2_signature:  continue  # recall that AAIG1_signature is i but AAIG1 contains the C,H of i-1!
                if selected_i_iminus1_dict and not AAIG2_signature in [q[0] for q in selected_i_iminus1_dict[AAIG1_signature]]:
                    continue
                AAIG2 = self.HCNH_spec.__get_AAIG__(AAIG2_signature)
                if set_all_to_1:
                    self.TOCSY_NOESY_2Dhist_intersections_mdict[AAIG1_signature][AAIG2_signature] = 1.0
                else:
                    total_intersection = 0.0
                    # print("DEBUG: AAIG1_signature = %s ; AAIG2_signature = %s" % (AAIG1_signature, AAIG2_signature))
                    for (peak1, peak2) in self.TOCSY_NOESY_peak_connectivities_mdict[AAIG1_signature][AAIG2_signature]:
                        updated_peak1 = AAIG1.__get_peak__(peak1)   # peaks1 and peaks2 are older version that don't have hist2D
                        updated_peak2 = AAIG2.__get_peak__(peak2)
                        total_intersection += updated_peak1.calc_intersection(updated_peak2)
                    # Save the total intersection between matching Peaks of AAIG1 and AAIG2
                    self.TOCSY_NOESY_2Dhist_intersections_mdict[AAIG1_signature][AAIG2_signature] = total_intersection

            self.scale_TOCSY_NOESY_2Dhist_connectivities(AAIG1_signature)


    def match_TOCSY_to_NOESY_Peaks(self,
                                   tolH=0.04, tolC=0.4,
                                   calc_intersection=True,
                                   uniquify=True,
                                   selected_i_iminus1_dict=None,
                                   clean_native_peaks=False):
        """
        Method to find connectivities of each AAIG in TOCSY with all other AAIGs in NOESY by comparing the 2D density
        histograms of their matching __Peaks__. The NOESY-TOCSY connectivities are stored in multidict
        self.TOCSY_NOESY_peak_connectivities_mdict. To expedite the connectivities calculation, you can provide a
        dict (selected_i_iminus1_dict) with a subset of connectivities that you want to calculate including the intersection
        and multi values.

        THEORY:
        each TOCSY 4D peak has C,H of i and N,HN of i+1
        each NOESY 4D peak has C,H of i,i+1,i-1,etc. and N,HN of i.
        The TOCSY(i-1)-NOESY(i) matching is based on the fact that each NOESY AAIG i contains peaks of i-1, too.
        As such, we iterate over all TOCSY AAIGs and label them according to their N,HN label (i) and compare them with
        all NOESY AAIGs which we label according to their N,HN label. By definition the TOCSY AAIG i will match with the NOESY
        AAIG i but also with the NOESY AAIG i-1. If we remove the label of AAIG i from the NOESY matches, we have the i-1 and the rest.
        That i-1 label corresponds to the TOCSY C,H shifts, and hence in the connectivities file we have 2 columns:
        i       possible i-1

        :return:
        """
        # Clean previously stored connectivities
        self.clean_connectivities()

        # Clean native peaks before you calculate 2D-histograms and intersections!!!
        if clean_native_peaks:
            self.remove_native_peaks_from_NOESY()

        # Now find connectivities without 2D-histogram comparison
        for AAIG1 in self.TOCSY_spec.__get_all_AAIGs__():   # contains the N,HN or residue i and the C,H or residue i-1
            self.match_TOCSY_AAIG_against_all_NOESY_Peaks(AAIG1, tolH=tolH, tolC=tolC,
                                                          selected_i_iminus1_dict=selected_i_iminus1_dict)

        if uniquify:
            # remove redundancy
            self.uniquify_TOCSY_NOESY_peak_connectivities(selected_i_iminus1_dict=selected_i_iminus1_dict)

        # Compare the 2D-histograms of each Peak pair to find the overall intersection of each connectivity
        if calc_intersection:
            # Clean previously stored intersections
            self.clean_intersections()
            # Scale intensities
            # IMPORTANT: self.__scale_intensities__() must be called outside of this function in order updated to be applied to
            # IMPORTANT: self.noesy_spec and self.TOCSY_spect under parallel execution!

            # First generate the 2D-histograms for all AAIGs in the spectrum
            if selected_i_iminus1_dict:
                NOESY_AAIG_signatures = [q[0] for v in list(selected_i_iminus1_dict.values()) for q in v]
                TOCSY_AAIG_signatures = list(selected_i_iminus1_dict.keys())
            else:
                NOESY_AAIG_signatures = self.HCNH_spec.__get_all_AAIG_signatures__()
                TOCSY_AAIG_signatures = self.TOCSY_spec.__get_all_AAIG_signatures__()
            try:
                shared.setConst(EVEN_EDGES=True,
                                SCALE_DENSITY=False)
            except TypeError:   # "TypeError: This constant already exists: EVEN_EDGES." (means we already invoked this method before)
                pass
            print(bcolors.BOLDBLUE + "Generating the 2D-histograms of all NOESY Peaks that participate in connectivities." + bcolors.ENDBOLD)
            # ## Serial NOESY Execution for debugging
            # for AAIG_signature in NOESY_AAIG_signatures:
            #     new_AAIG = self.NOESY_Peaks2hist(AAIG_signature, even_edges=True, scale_density=False)
            #     self.HCNH_spec.__replace_AAIG__(new_AAIG)
            # Parallel NOESY Execution
            new_AAIGs = list(futures.map(self.NOESY_Peaks2hist,
                                         NOESY_AAIG_signatures,
                                         [True] * len(NOESY_AAIG_signatures),
                                         [False] * len(NOESY_AAIG_signatures)
                                         )
                             )
            for aaig in new_AAIGs:  self.HCNH_spec.__replace_AAIG__(aaig)

            print(bcolors.BOLDBLUE + "Generating the 2D-histograms of all TOCSY Peaks that participate in connectivities." + bcolors.ENDBOLD)
            # ## Serial TOCSY Execution for debugging
            # for AAIG_signature in TOCSY_AAIG_signatures:
            #     new_AAIG = self.TOCSY_Peaks2hist(AAIG_signature)
            #     self.TOCSY_spec.__replace_AAIG__(new_AAIG)
            ## Parallel TOCSY Execution
            new_AAIGs = list(futures.map(self.TOCSY_Peaks2hist, TOCSY_AAIG_signatures))
            for aaig in new_AAIGs:  self.TOCSY_spec.__replace_AAIG__(aaig)

            # Now that you have the Peak 2D-histograms measure all the intersections!
            self.measure_all_TOCSY_NOESY_Peak_2Dhist_intersections(tolH=tolH, tolC=tolC,
                                                                   selected_i_iminus1_dict=selected_i_iminus1_dict)

            # # Release memory
            # self.noesy_spec.__delete_all_hist2D__()
            # self.TOCSY_spec.__delete_all_hist2D__()
            # gc.collect()
        else:   # set all intersections to 1.0
            self.measure_all_TOCSY_NOESY_Peak_2Dhist_intersections(tolH=tolH, tolC=tolC,
                                                                   selected_i_iminus1_dict=selected_i_iminus1_dict,
                                                                   set_all_to_1=True)

        # deconvolute connectivities
        if calc_intersection:
            self.deconvolute_TOCSY_NOESY_connectivities()

        # populate self.AAIG_connectivities_complete_dict and self.i_iminus1_complete_dict (the are the same!)
        # TODO: control oexp parameter
        self.AAIG_connectivities_complete_dict = self.get_TOCSY_NOESY_connectivities_dict()

    def pairwise_TOCSY_to_NOESY_match_AAIGs(self, TOCSY_AAIG_signature, NOESY_AAIG_signature, tolH=0.04, tolC=0.4):
        """
        Method to find connectivities of each AAIG in NOESY (i) with all other AAIGs in TOCSY (i-1).
        The NOESY-TOCSY connectivities are stored in multidict self.TOCSY_NOESY_peak_connectivities_mdict.
        :return:
        """
        # First generate the 2D-histograms for all AAIGs in the spectrum
        self.find_2Dhist_boundaries()
        aaig1 = self.NOESY_AAIG2hist(NOESY_AAIG_signature)
        aaig1.__print_peaks__()
        aaig2 = self.TOCSY_AAIG2hist(TOCSY_AAIG_signature)
        aaig2.__print_peaks__()
        intersection = aaig1.calc_intersection(aaig2)
        print("Intersection between NOESY AAIG " + NOESY_AAIG_signature + " and TOCSY AAIG " + TOCSY_AAIG_signature + " = " + str(intersection))

    def write_TOCSY_NOESY_connectivities(self, fname, pickle_fname="", oexp=3.0448217240563293):
        """

        :param fname:
        :param oexp:    3.0448217240563293 is globally optimized value in the set of 9 globular proteins for global intensity scaling
                        and peak-peak 2D-hist connectivities.
        :return:
        """
        if not pickle_fname:
            pickle_fname = os.path.splitext(fname)[0] + ".pkl"
        with open(fname, 'w') as f:
            f.write("i\tpossible i-1\n")
            for AAIG1_signature in list(self.i_iminus1_complete_dict.keys()):    # i -> i-1
                # NOTE: the following commented code is not necessary since self.i_iminus1_complete_dict is
                # NOTE: already populated and is the same as self.AAIG_connectivities_complete_dict.

                # # NOTE: Before the introduction of 'multi', conectivities were sorted byt occupancy and secondarily by resid.
                # #       After the introduction of 'multi', they are sorted by 'multi' only, because it is unique for each connection.
                # #       Moreover, replacement of AAIG names with signatures (e.g. R34N-H) renders sorting by resid more complicated,
                # #       therefore I decided to completely remove it.
                # connectivities_list = [(AAIG2_signature,
                #                         len(self.TOCSY_NOESY_peak_connectivities_mdict[AAIG1_signature][AAIG2_signature]),
                #                         len(self.TOCSY_spec.__get_AAIG__(AAIG1_signature).__get_all_peaks__()),
                #                         self.TOCSY_NOESY_2Dhist_intersections_mdict[AAIG1_signature][AAIG2_signature],
                #                         self.TOCSY_NOESY_2Dhist_intersections_mdict[AAIG1_signature][AAIG2_signature] * \
                #                         (float(len(self.TOCSY_NOESY_peak_connectivities_mdict[AAIG1_signature][AAIG2_signature])) / \
                #                          len(self.TOCSY_spec.__get_AAIG__(AAIG1_signature).__get_all_peaks__())) ** oexp)
                #                        for AAIG2_signature in list(self.TOCSY_NOESY_peak_connectivities_mdict[AAIG1_signature].keys())]
                # connectivities_list.sort(key=itemgetter(4), reverse=True)   # sort by 'multi'
                # connectivities_list = [(c[0], c[1], c[2], c[3], c[4]) for c in connectivities_list]

                connectivities_list = self.i_iminus1_complete_dict[AAIG1_signature]
                if len(connectivities_list) > 0:
                    f.write(AAIG1_signature + "\t" + ', '.join(str(c) for c in connectivities_list) + "\n")
        # Save the same information in a pickle file, along withthe dictionarry with the all matching peaks between pairs
        # of AAIGs.
        save_pickle(pickle_fname, self.AAIG_connectivities_complete_dict,
                    self.TOCSY_NOESY_peak_connectivities_mdict)  # save the filtered versions of both dicts

    def remove_native_peaks_from_NOESY(self):
        """
        APPLICABLE ONLY IN TOCSY->NOESY MATCHING.
        This method must be called after you load the TOCSY and NOESY files to AAIGS.
        When you match TOCSY peaks of AAIG(i) to NOESY peaks in order to find the AAIG(i-1),
        you don't need to use the AAIG(i) in NOESY (these are common between NOESY and TOCSY).
        This function iterates over all TOCSY AAIG(i), finds the NOESY AAIG(i) and removes the peaks that
        match from NOESY_contents. There is a problem with the removal of peaks from the same AAIG. Look
        at this example from nEIt. There are 3 consecutive E in the sequence, E37,E38,E39.

        In [62]: [w for w in NOESY_contents if w[0]=='E38']
        Out[62]:
        [['E38', '4.684', '55.456', '121.017', '8.064', '10353967.0'],
         ['E38', '1.022', '22.189', '120.987', '8.066', '3772312.0'],
         ['E38', '3.054', '42.354', '120.975', '8.069', '3508862.0'],
         ['E38', '4.011', '60.287', '121.009', '8.066', '8719190.0'],
         ['E38', '2.136', '28.789', '120.997', '8.067', '32965242.0'],
         ['E38', '4.207', '59.453', '120.970', '8.066', '41168996.0'],
         ['E38', '2.432', '36.962', '121.004', '8.066', '22151098.0'],
         ['E38', '2.337', '36.957', '121.002', '8.067', '16331055.0'],
         ['E38', '2.251', '29.850', '121.003', '8.067', '60914216.0']]

        In [63]: [w for w in TOCSY if w[0] in ['E38']]
        Out[63]:
        [['E38', '4.010', '60.370', '120.991', '8.063'],
         ['E38', '2.429', '36.937', '120.991', '8.063'],
         ['E38', '2.343', '36.943', '120.986', '8.063'],
         ['E38', '2.133', '28.864', '120.995', '8.061']]

        In [64]: [w for w in TOCSY if w[0] in ['E39']]
        Out[64]:
        [['E39', '4.207', '59.379', '121.895', '8.254'],
         ['E39', '2.430', '36.954', '121.892', '8.253'],
         ['E39', '2.328', '36.986', '121.886', '8.253'],
         ['E39', '2.250', '29.810', '121.908', '8.254']]

        As you may notice, there is only one pair of carbons 36 ppm (only one methylene) in NOESY. Normally, it should be two at least because both E38 and E39 have this type of methylene (E37 too). If I remove the TOCSY peaks of E38 from the NOESY list I end up with this:

        In [61]: [w for w in clean_NOESY_contents if w[0]=='E38']
        Out[61]:
        [['E38', '4.684', '55.456', '121.017', '8.064', '10353967.0'],
         ['E38', '1.022', '22.189', '120.987', '8.066', '3772312.0'],
         ['E38', '3.054', '42.354', '120.975', '8.069', '3508862.0'],
         ['E38', '4.207', '59.453', '120.970', '8.066', '41168996.0'],
         ['E38', '2.251', '29.850', '121.003', '8.067', '60914216.0']]

        As you can see there are no carbons at 36 ppm to match with E39. This results in lower occupancy ('E38', 2, 4)<--E39.

        This was just one scenario (when you have consecutive residues of the same type you lose connectivities). Although the connectivities become fewer, the NH-mapping deteriorates when I remove NOESY peaks from the same AAIG when I do the TOCY->NOESY matching. Therefore I suggest to leave it as it is.
        """
        print("Removing native peaks from NOESY AAIGs.")
        NOESY_contents2remove = []  # the NOESY lines that contains peaks of AAIG(i)
        for TOCSY_AAIG in self.TOCSY_spec.__get_all_AAIGs__():  # residue i
            for TOCSY_peak in TOCSY_AAIG.__get_all_peaks__():
                try:
                    NOESY_AAIG = self.HCNH_spec.__get_AAIG__(TOCSY_AAIG.signature)  # save the NOESY AAIG to modify it
                except KeyError:    # if this TOCSY AAIG label does not exist in NOESY, move to the next
                    continue
                for NOESY_peak in NOESY_AAIG.__get_all_peaks__():      # residue i
                    if NOESY_peak.does_peak_match(tolH=self.tolH,    # remove i-1 C,H resonances from i in NOESY
                                                  tolC=self.tolC,
                                                  peak2=TOCSY_peak,
                                                  attypes=['C', 'H']):
                        NOESY_AAIG.__del_peak__(peak=NOESY_peak)
                        break   # delete only the first peak that matches (multiple can match but not all of them come from i-1)
                self.HCNH_spec.__replace_AAIG__(NOESY_AAIG)    # replace with the modified AAIG in the NOESY spectrum

    def get_TOCSY_NOESY_connectivities_dict(self, oexp=3.0448217240563293):
        """
        Method that returns all possible connectivities of every TOCSY AAIG with NOESY AAIGs in a dictionary form.
        It presumes that match_TOCSY_to_NOESY_Peaks() has been already called, otherwise it will return an empty
        dictionary.

        :param oexp:    3.0448217240563293 is globally optimized value in the set of 9 globular proteins for global
                        intensity scaling and peak-peak 2D-hist connectivities.
        :return: i_iminus1_complete_dict: dictionary of the form: TOCSY AAIG i signature --> [(NOESY AAIG i-1 signature,
                                            occupancy, numOfResonances, intersection, multi score), ...]
        """

        for AAIG1_signature in list(self.TOCSY_NOESY_peak_connectivities_mdict.keys()):  # i -> i-1
            # NOTE: Before the introduction of 'multi', conectivities were sorted byt occupancy and secondarily by resid.
            #       After the introduction of 'multi', they are sorted by 'multi' only, because it is unique for each connection.
            #       Moreover, replacement of AAIG names with signatures (e.g. R34N-H) renders sorting by resid more complicated,
            #       therefore I decided to completely remove it.
            try:
                connectivities_list = [(AAIG2_signature,
                                        len(self.TOCSY_NOESY_peak_connectivities_mdict[AAIG1_signature][AAIG2_signature]),
                                        len(self.TOCSY_spec.__get_AAIG__(AAIG1_signature).__get_all_peaks__()),
                                        self.TOCSY_NOESY_2Dhist_intersections_mdict[AAIG1_signature][AAIG2_signature],
                                        self.TOCSY_NOESY_2Dhist_intersections_mdict[AAIG1_signature][AAIG2_signature] * \
                                        (float(len(self.TOCSY_NOESY_peak_connectivities_mdict[AAIG1_signature][AAIG2_signature])) / \
                                         len(self.TOCSY_spec.__get_AAIG__(AAIG1_signature).__get_all_peaks__())) ** oexp)
                                       for AAIG2_signature in list(self.TOCSY_NOESY_peak_connectivities_mdict[AAIG1_signature].keys())]
            except IndexError:
                # if AAIG2_signature == "" or len(connectivities_list) == 0:
                #     pass
                # else:
                print("DEBUG: AAIG1_signature=", AAIG1_signature)
                print("DEBUG: AAIG2_signature=", AAIG2_signature)
                print("DEBUG: type(AAIG2_signature)=", type(AAIG2_signature))
                print("DEBUG: connectivities_list=", connectivities_list)
                raise IndexError()
            connectivities_list.sort(key=itemgetter(4), reverse=True)  # sort by 'multi'
            connectivities_list = [(c[0], c[1], c[2], c[3], c[4]) for c in connectivities_list]
            self.i_iminus1_complete_dict[AAIG1_signature] = connectivities_list
        return self.i_iminus1_complete_dict

    def write_TOCSY_NOESY_common_peak_list(self, out_fname, tolH=0.04, tolC=0.4):
        """
        Method to rewrite the HCNH NOESY file but only with the peaks that belong to the native AAIG.
        Recall that
        :return:
        """
        # TODO: uniquify peaks and make it more sophisticated.
        noesy_spec = HCNH_Spectrum()
        for tAAIG in self.TOCSY_spec.__get_all_AAIGs__():  # contains the N,HN or residue i and the C,H or residue i-1
            tAAIG_signature = tAAIG.signature
            for peak1 in tAAIG.__get_all_peaks__():  # the C,H come from i-1 but the N,HN resonances from i
                try:
                    nAAIG = self.HCNH_spec.__get_AAIG__(tAAIG_signature)
                except KeyError:
                    ColorPrint("WARNING: the TOCSY AAIG %s does not exist in HCNH NOESY." % tAAIG_signature, "WARNING")
                    continue
                for peak2 in nAAIG.__get_all_peaks__():
                    if peak1.does_peak_match(peak2, tolH=tolH, tolC=tolC, attypes=['C', 'H']):
                        # print("DEBUG: %f == %f && %f == %f ?" % (peak1.Creson, peak2.Creson, peak1.Hreson, peak2.Hreson))
                        if not noesy_spec.__exists__(tAAIG_signature):
                            # by convention name the AAIGs by their NH peak label
                            aaig = AAIG(signature=tAAIG_signature,
                                        NH_name=get_NH_name(tAAIG_signature),
                                        user_assignment=dash_to_NH(tAAIG_signature))
                            noesy_spec.__add_AAIG__(signature=aaig.signature, AAIG=aaig)
                        noesy_spec.__add_peak2AAIG__(peak=peak2, AAIG_signature=tAAIG_signature)
        noesy_spec.write_sparky_list(out_fname) # write a sparky list with the common peaks

# if __name__ == "__main__":
#
#     try:
#         # CONSUMES TOO MUCH MEMORY WITH SMALL BIN LENGTHS!!!
#         H_bin_length = 0.0005
#         C_bin_length = 0.005
#         # For Peak 2D-histograms try bandwidths 0.1-0.3
#         bandwidth = 0.2
#         tolH = 0.04
#         tolC = 0.4
#
#         ## HCNH-HCNH Connectivities example
#         """
#         NOESY_con = Connectivities("/home2/thomas/Documents/4D-CHAINS_regtests/nEIt/onlyNOESY_assignment/NHmap",
#                                          FIRST_RESIDUE_NUMBER=1,
#                                          tolH=tolH,
#                                          tolC=tolC,
#                                          H_bin_length=H_bin_length,
#                                          C_bin_length=C_bin_length,
#                                          bandwidth=bandwidth)
#         NOESY_con.load_NOESY_spectrum("/home2/thomas/Documents/4D-CHAINS_regtests/nEIt/onlyNOESY_assignment/nEIt_INTnoesy.23.9.2016num.list")
#         NOESY_con.load_TOCSY_spectrum("/home2/thomas/Documents/4D-CHAINS_regtests/nEIt/onlyNOESY_assignment/nE1t_INTtocsyHEIGHT.list")
#         # NOESY_con.pairwise_TOCSY_to_NOESY_match("G88", "V89")
#         # NOESY_con.match_all_against_all()
#         NOESY_con.match_TOCSY_to_NOESY(tolH=tolH, tolC=tolC)
#         # NOESY_con.write_NOESY_NOESY_connectivities("/home2/thomas/Documents/4D-CHAINS_regtests/nEIt/onlyNOESY_assignment/"
#         #                                "NOESY_connectivitis.tolH0.02_tolC0.2_Hbin"+str(H_bin_length)+"_Cbin"
#         #                                +str(C_bin_lengt
#         NOESY_con.write_TOCSY_NOESY_connectivities("/home2/thomas/Documents/4D-CHAINS_regtests/nEIt/onlyNOESY_assignment/"
#                                            "TOCSY_connectivities.tolH"+str(tolH)+"_tolC"+str(tolC)+"_Hbin"+str(H_bin_length)+"_Cbin"
#                                            +str(C_bin_length)+"_bw"+str(bandwidth)+".new.list")
#         """
#         ## NOESY-TOCSY Connectivities example
#         con = Connectivities(tolH=tolH,
#                                    tolC=tolC,
#                                    H_bin_length=H_bin_length,
#                                    C_bin_length=C_bin_length,
#                                    bandwidth=bandwidth)
#         # con.load_NOESY_spectrum("/home2/thomas/Documents/4D-CHAINS_regtests/aLP/aLP_devtest/aLP_INTnoesy.23.9.2016num.list")
#         # con.load_TOCSY_spectrum("/home2/thomas/Documents/4D-CHAINS_regtests/aLP/aLP_devtest/aLP_INTtocsy.23.9.2016num.list")
#         con.load_NOESY_spectrum("/home2/thomas/Documents/4D-CHAINS_regtests/aLP/aLP_devtest/sample_NOESY.list")
#         con.load_TOCSY_spectrum("/home2/thomas/Documents/4D-CHAINS_regtests/aLP/aLP_devtest/sample_TOCSY.list")
#         con.__scale_intensities__(wrt_aaig=False)
#         con.match_TOCSY_to_NOESY_Peaks(tolH=tolH, tolC=tolC)
#         i_iminus1_dict = con.get_TOCSY_NOESY_connectivities()
#         # print "DEBUG: i_iminus1_dict="
#         # print_dict_contents(i_iminus1_dict)
#         con.write_TOCSY_NOESY_connectivities(
#             "TOCSY-HCNH_connectivities.tolH" + str(tolH) + "_tolC" + str(tolC) + "_Hbin" + str(H_bin_length) + "_Cbin"
#             + str(C_bin_length) + "_bw" + str(bandwidth) + ".specIntensities.list")
#         del con
#         gc.collect()
#
#
#     except:
#         type, value, tb = sys.exc_info()
#         lines = traceback.format_exception(type, value, tb)
#         print(''.join(lines))
#         raise