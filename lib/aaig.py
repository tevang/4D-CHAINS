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


import numpy as np
from sklearn.neighbors import KernelDensity
import matplotlib.pyplot as plt
from .global_func import *
from .inthist import *
from scipy.sparse import csr_matrix
import uuid
import pandas as pd
from decimal import *
getcontext().prec = 7   # 5 decimal points

class AAIG():
    # Class to store all information regarding an Amino Acid Index Group (aka spin system, or set of peaks).

    def __init__(self,
                 signature="",
                 CH_name="",
                 NH_name="",
                 overlapping_AAIGs=[],
                 user_assignment="",
                 H_bin_length=Dec(0.01),
                 C_bin_length=Dec(0.1),
                 al_H_minmax=np.array([Dec(-0.42), Dec(6.496)]),
                 al_C_minmax=np.array([Dec(6.762), Dec(75.934)]),
                 ar_H_minmax=np.array([Dec(5.831), Dec(8.021)]),
                 ar_C_minmax=np.array([Dec(113.675), Dec(135.581)]), bandwidth=0.0025):
        """
        2D-histograms have H value range [-3.78, 10.22] and C value range [1.90, 143.90].
        :param H_bin_length:
        :param C_bin_length:
        :param H_minmax:
        :param C_minmax:
        :param bandwidth:
        """
        self.ID = uuid.uuid4()  # unique identifier for this AAIG object
        self.signature = signature  # always refers to NH_name, e.g. X40NXHX or F45NH
        self.CH_name = CH_name
        self.NH_name = NH_name
        self.overlapping_AAIGs = overlapping_AAIGs
        self.user_assignment = user_assignment
        self.H_bin_length = Dec(H_bin_length)
        self.C_bin_length = Dec(C_bin_length)
        self.al_H_minmax = np.array([Dec(al_H_minmax[0]), Dec(al_H_minmax[1])])
        self.al_C_minmax = np.array([Dec(al_C_minmax[0]), Dec(al_C_minmax[1])])
        self.ar_H_minmax = np.array([Dec(ar_H_minmax[0]), Dec(ar_H_minmax[1])])
        self.ar_C_minmax = np.array([Dec(ar_C_minmax[0]), Dec(ar_C_minmax[1])])
        self.bandwidth = bandwidth
        self.peaks, self.protons, self.carbons, self.intensities = [], [], [], []
        self.hist2D = None

    def __add_peak__(self, peak):
        self.peaks.append(peak)

    def __del_peak__(self, peak):
        for i, p in enumerate(self.peaks):
            if peak.does_peak_match(p, tolC=0, tolH=0, tolHN=0, tolN=0):    # if not tol=0 it will delete wrong peak!
                del self.peaks[i]
                break

    def __replace_peak__(self, peak):
        """
        Replaces the existing Peak object with a new one that updated properties. The resonances must match!
        :param peak:
        :return:
        """
        for i, p in enumerate(self.peaks):
            if peak.does_peak_match(p, tolC=0, tolH=0, tolHN=0, tolN=0):    # if not tol=0 it will replace wrong peak!
                self.peaks[i] = peak    # replace this Peak with the updated version
                break

    def __get_peak__(self, peak):
        """
        Gets the Peak object (usually with updated properties) that matches with the given one. The resonances must match!
        :param peak:
        :return:
        """
        for i, p in enumerate(self.peaks):
            if peak.does_peak_match(p, tolC=0, tolH=0, tolHN=0, tolN=0):    # if not tol=0 it will return wrong peak!
                return self.peaks[i]

    def __get_peak_by_ID__(self, peakID):
        """
        Gets the Peak object (usually with updated properties) that has the given peak ID.
        :param peak:
        :return:
        """
        for i, p in enumerate(self.peaks):
            if p.ID == peakID:
                return self.peaks[i]

    def __get_all_peaks__(self):
        return self.peaks

    def __get_CH_dataframe_for_AAtype_prediction__(self):
        """

        :return: a Pandas dataframe in the format required by the NOESY AA-type Classifier.
        """
        df = pd.DataFrame([], columns=("atom1", "atom2"))
        for peak in self.__get_all_peaks__():
            df = df.append({'atom1': peak.Creson, 'atom2': peak.Hreson}, ignore_index=True)
        return df

    def __print_peaks__(self, only_CHresons=False):
        if only_CHresons:
            CHresons = [p.__get_CHresons__() for p in self.peaks]
            CHresons.sort(key=itemgetter(1))
            for CH in CHresons:
                print("Hreson=", CH[0], "Creson=", CH[1])
        else:
            for peak in self.peaks:
                peak.__print_peak__()

    def __is_equal__(self, aaig):
        peaks1 = self.__get_all_peaks__()
        peaks2 = aaig.__get_all_peaks__()
        matched_peaks1, matched_peaks2 = [], []
        for peak1 in peaks1:
            for peak2 in peaks2:
                if peak2 in matched_peaks2:
                    continue
                if peak1.does_peak_match(peak2, tolH=0.0, tolC=0., tolHN=0.0, tolN=0.0):
                    matched_peaks1.append(peak1)
                    matched_peaks2.append(peak2)
                    break
        # If all the peaks between AAIG1 and AAIG2 match
        if len(matched_peaks1) == len(peaks1) and \
            len(matched_peaks1) == len(peaks2) and \
            len(peaks1) == len(peaks2):
            return True
        else:
            return False


    def __delete_hist__(self):
        del self.hist2D
        self.hist2D = None

    def __scale_intensities__(self, even_intensities=False):
        """
        OBSOLE: the same operations are conducted from within Spectrum() class.
        Scale the intensities of all Peaks of this AAIG to 0-1.0 and save the value into Peak.scaled_Intensity.
        :return:
        """
        intensities = np.array([p.Intensity for p in self.__get_all_peaks__() if p.Intensity])
        if intensities.size > 0 and intensities.max() != 1.0:  # if needed, scale intensities (range 0.0-1.0)
            for peak in self.__get_all_peaks__():
                if even_intensities:
                    peak.scaled_Intensity = 1.0
                else:
                    peak.scaled_Intensity = peak.Intensity/float(intensities.max())
                self.__replace_peak__(peak)  # save the Peak with scaled_Intensity property

        elif intensities.size == 0 and not even_intensities:
            print(bcolors.WARNING + "WARNING: AAIG " + self.signature + " has peaks without intensities!" + bcolors.ENDC)
        else:
            for peak in self.__get_all_peaks__():
                peak.Intensity = 1.0
                peak.scaled_Intensity = 1.0

    def __set_hist_borders__(self,
                             alHmin=Dec(-0.5),
                             alHmax=Dec(6.5),
                             alCmin=Dec(8.8),
                             alCmax=Dec(74.0),
                             arHmin=Dec(6.3),
                             arHmax=Dec(77),
                             arCmin=Dec(113.7),
                             arCmax=Dec(135.6)):
        self.al_H_minmax = np.array([Dec(alHmin), Dec(alHmax)])
        self.al_C_minmax = np.array([Dec(alCmin), Dec(alCmax)])
        self.ar_H_minmax = np.array([Dec(arHmin), Dec(arHmax)])
        self.ar_C_minmax = np.array([Dec(arCmin), Dec(arCmax)])

    def even_edges(self, H_minmax, C_minmax):
        """
        Adjusts the extrema of the H or the C axis in order the have the same number of bins.
        This method is called before smoothing.
        :return:
        """
        assert H_minmax[0] < H_minmax[1], Debuginfo("ERROR: the H axis does not have right extrema (%f < %f)" %
                                                    (H_minmax[0], H_minmax[1]), fail=True)
        assert C_minmax[0] < C_minmax[1], Debuginfo("ERROR: the C axis does not have right extrema (%f < %f)" %
                                                    (C_minmax[0], C_minmax[1]), fail=True)
        # Convert all numbers to Decimal objects with fixed decimal points
        H_minmax[0], H_minmax[1] = Dec(H_minmax[0]), Dec(H_minmax[1])
        C_minmax[0], C_minmax[1] = Dec(C_minmax[0]), Dec(C_minmax[1])
        self.H_bin_length = Dec(self.H_bin_length)
        self.C_bin_length = Dec(self.C_bin_length)
        print(bcolors.BOLDBLUE + "Making H and C edges even." + bcolors.ENDBOLD)
        bindiff_tol = 1
        binNumH = int((H_minmax[1] - H_minmax[0]) / self.H_bin_length) - 1
        binNumC = int((C_minmax[1] - C_minmax[0]) / self.C_bin_length) - 1
        while abs(binNumH-binNumC) > bindiff_tol:
            if binNumH-binNumC > bindiff_tol:
                C_minmax[0] -= self.C_bin_length
                C_minmax[1] += self.C_bin_length
            elif binNumC-binNumH > bindiff_tol:
                H_minmax[0] -= self.H_bin_length
                H_minmax[1] += self.H_bin_length
            binNumH = int((H_minmax[1] - H_minmax[0]) / self.H_bin_length) - 1
            binNumC = int((C_minmax[1] - C_minmax[0]) / self.C_bin_length) - 1
        # Finally pad only one end to even out the difference in the bin number
        if binNumH - binNumC > 0:
            C_minmax[0] -= self.C_bin_length
        elif binNumC - binNumH > 0:
            H_minmax[0] -= self.H_bin_length

        return H_minmax, C_minmax

    def create_2Dhist(self, augm_protons, augm_carbons, H_minmax, C_minmax, even_edges=True):
        """

        :param augm_protons:
        :param augm_carbons:
        :param H_minmax:
        :param C_minmax:
        :param even_edges:
        :return:
        """
        if augm_protons.size == 0 and augm_carbons.size == 0:   # if empty lists were given (e.g. no aromatics)
            return [], [], np.array([])
        if even_edges:
            # VERY IMPORTANT: the two axes must have the same number of bins in order smoothing to work correctly!!!
            H_minmax, C_minmax = self.even_edges(H_minmax, C_minmax)
        binNumH = int((H_minmax[1] - H_minmax[0]) / self.H_bin_length) # - 1 + 0j  # convert the bin number in x-axis to complex number to emulate e.g. "100j"
        binNumC = int((C_minmax[1] - C_minmax[0]) / self.C_bin_length) # - 1 + 0j  # convert the bin number in y-axis to complex number to emulate e.g. "100j"
        if binNumH != binNumC:
            raise Exception(bcolors.FAIL + "ERROR: the number of bins in H and in C axes are not equal! " +
                            str(binNumH)+" != "+str(binNumC) +
                            " H_minmax= " + str(H_minmax[0]) + ", " + str(H_minmax[1]) +
                            " C_minmax= " + str(C_minmax[0]) + ", " + str(C_minmax[1]) +
                            "\n" + "Use aaig.AAIG2hist(even_edges=True) or set "
                            "the parameters H_bin_length, C_bin_length, al_H_minmax, al_C_minmax, ar_H_minmax, "
                            "ar_C_minmax appropriately." + bcolors.ENDC)

        print(bcolors.BOLDGREEN + "H bin length = " + str(self.H_bin_length) + " ppm." + bcolors.ENDBOLD)
        print(bcolors.BOLDGREEN + "C bin length = " + str(self.C_bin_length) + " ppm." + bcolors.ENDBOLD)
        print(bcolors.BOLDGREEN + "2D-histogram dimension = " + str(binNumH) + " x " + str(binNumC) + bcolors.ENDBOLD)

        H_edges, C_edges, hist2D = IntHistogram().kde2D(augm_protons,
                                                        augm_carbons,
                                                        bandwidth=self.bandwidth,
                                                        xbins=binNumH,
                                                        ybins=binNumC,
                                                        x_minmax=(0., 1.),
                                                        y_minmax=(0., 1.),  # Carbon values in ppm have 10x larger scale than Proton
                                                        kernel='cosine',
                                                        metric='euclidean')
        return H_edges, C_edges, hist2D

    def AAIG2hist(self, even_edges = False, even_intensities=False, scale_density=False, **hist_borders):
        """
        Method to convert a list of NOESY peaks of a specific AAIG, to a 2D-histogram where the x and y axes
        are the proton and carbon resonances, while z axis will be the intensity of the peak.
        :param peaks: list [[Hreson, Creson, norm. intensity], ...]
        :param H_bin_length:
        :param C_bin_length:
        :param H_minmax:
        :param C_minmax:
        :return:
        """
        print(bcolors.OKBLUE + "Generating the 2D-histogram of AAIG " + self.signature + bcolors.ENDC)

        if not hist_borders:    # if not provided, use the default values
            hist_borders = {
            'arHmax': Dec(7.646),
            'alCmax': Dec(73.984),
            'alHmin': Dec(-0.42),
            'arCmax': Dec(135.581),
            'alHmax': Dec(6.496),
            'alCmin': Dec(8.762),
            'arCmin': Dec(113.675),
            'arHmin': Dec(6.21),
            }

        self.al_H_minmax =  np.array([Dec(hist_borders['alHmin']), Dec(hist_borders['alHmax'])])
        self.al_C_minmax =  np.array([Dec(hist_borders['alCmin']), Dec(hist_borders['alCmax'])])
        self.ar_H_minmax =  np.array([Dec(hist_borders['arHmin']), Dec(hist_borders['arHmax'])])
        self.ar_C_minmax =  np.array([Dec(hist_borders['arCmin']), Dec(hist_borders['arCmax'])])

        if even_intensities:
            intensities = [1.0 for p in self.peaks]
        else:
            intensities = [p.scaled_Intensity for p in self.peaks]
        if None in self.intensities:
            raise Warning(
                "AAIG " + self.signature + " could not be converted to 2D-histogram because there aren't intensity values for all its peaks!")
            return
        intensities = np.array(intensities)
        if intensities.size > 0 and intensities.max() > 1.0:    # valid for both local and global intensity scaling
            raise Exception(bcolors.FAIL + "ERROR: intensities are not scaled in order to calculate the 2D-histogram of AAIG " + self.signature +
                            "!" + bcolors.ENDC)
        else:
            scaled_intensities = intensities
        # Add an extra property to each peak, the "scaled_Intensity"
        for i in range(len(self.peaks)):
            self.peaks[i].scaled_Intensity = scaled_intensities[i]

        al_protons, ar_protons, al_carbons, ar_carbons, al_intensities, ar_intensities = [], [], [], [], [], []
        for p in self.peaks:
            if p.Creson < 100:
                al_carbons.append(p.Creson)
                al_protons.append(p.Hreson)
                if even_intensities:
                    al_intensities.append(1.0)
                else:
                    al_intensities.append(p.scaled_Intensity)
            elif p.Creson >= 100:
                ar_carbons.append(p.Creson)
                ar_protons.append(p.Hreson)
                if even_intensities:
                    ar_intensities.append(1.0)
                else:
                    ar_intensities.append(p.scaled_Intensity)

        self.al_protons = np.array(al_protons)
        self.al_carbons = np.array(al_carbons)
        self.al_intensities = np.array(al_intensities)
        self.ar_protons = np.array(ar_protons)
        self.ar_carbons = np.array(ar_carbons)
        self.ar_intensities = np.array(ar_intensities)

        self.al_H_edges, self.ar_H_edges, self.al_C_edges, self.ar_C_edges, self.al_hist2D, self.ar_hist2D = \
            None, None, None, None, None, None
        # include the intensity implicitly by augmenting the respective H and C resonances accordingly
        augm_al_protons, augm_al_carbons = IntHistogram().augment_peaks_by_intensity(self.al_protons,
                                                                                     self.al_carbons,
                                                                                     self.al_intensities)
        augm_ar_protons, augm_ar_carbons = IntHistogram().augment_peaks_by_intensity(self.ar_protons,
                                                                                     self.ar_carbons,
                                                                                     self.ar_intensities)
        # IMPORTANT: Scale the proton and carbon resonances between 0. and 1., otherwise smoothing won't work equally in both dimensions!!!
        augm_al_protons = minmax_scale(augm_al_protons, minmax=[self.al_H_minmax[0], self.al_H_minmax[1]])
        augm_ar_protons = minmax_scale(augm_ar_protons, minmax=[self.ar_H_minmax[0], self.ar_H_minmax[1]])
        self.scaled_al_protons = minmax_scale(al_protons, minmax=[self.al_H_minmax[0], self.al_H_minmax[1]])
        self.scaled_ar_protons = minmax_scale(ar_protons, minmax=[self.ar_H_minmax[0], self.ar_H_minmax[1]])
        augm_al_carbons = minmax_scale(augm_al_carbons, minmax=[self.al_C_minmax[0], self.al_C_minmax[1]], final_range=[0., 1.])
        augm_ar_carbons = minmax_scale(augm_ar_carbons, minmax=[self.ar_C_minmax[0], self.ar_C_minmax[1]], final_range=[0., 1.])
        self.scaled_al_carbons = minmax_scale(al_carbons, minmax=[self.al_C_minmax[0], self.al_C_minmax[1]], final_range=[0., 1.])
        self.scaled_ar_carbons = minmax_scale(ar_carbons, minmax=[self.ar_C_minmax[0], self.ar_C_minmax[1]], final_range=[0., 1.])

        ##~~~ ALIPHATIC 2D-hist ~~~##
        self.al_H_edges, self.al_C_edges, al_hist2D = self.create_2Dhist(augm_al_protons, augm_al_carbons, self.al_H_minmax,
                                                                                   self.al_C_minmax, even_edges=even_edges)
        self.al_hist2D = csr_matrix(al_hist2D)
        ##~~~ AROMATIC 2D-hist ~~~##
        self.ar_H_edges, self.ar_C_edges, ar_hist2D = self.create_2Dhist(augm_ar_protons, augm_ar_carbons, self.ar_H_minmax,
                                                                         self.ar_C_minmax, even_edges=even_edges)
        self.ar_hist2D = csr_matrix(ar_hist2D)

        # Scale the densities to sum to 1.0 in both aliphatic and aromatic histograms
        if (self.al_hist2D.sum() + self.ar_hist2D.sum()) == 0.0:    # instead of scaling densities use this
            raise Exception(bcolors.FAIL + "ERROR: the bandwidth value (" + str(self.bandwidth) +" is too low to create density for"
                                        " AAIG " + self.signature + " ! Try to increase it and rerun the program." + bcolors.ENDC)
        if scale_density:
            if self.ar_hist2D.size > 0: # if there are aromatic peaks in this AAIG
                total_density = float(self.al_hist2D.sum()+self.ar_hist2D.sum())
                self.al_hist2D = self.al_hist2D/total_density
                self.ar_hist2D = self.ar_hist2D/total_density
            else:
                self.al_hist2D = self.al_hist2D/float(self.al_hist2D.sum())

    def Peaks2hist(self,
                   even_edges=True,
                   even_intensities=False,
                   scale_density=False,):
        """
        Calculates the 2D-histograms of all the peaks in this AAIG.
        :return:
        """
        print(bcolors.OKBLUE + "Generating 2D-histograms from all peaks of AAIG " + self.signature + bcolors.ENDC)

        for Peak in self.__get_all_peaks__():
            # if Peak.i_AAIG_signature + "-" + Peak.j_AAIG_signature == '?-I206':   # DEBUGGING
            # Peak.__print_peak__()
            Peak.peak2hist(even_edges=even_edges, scale_density=scale_density)
            # if Peak.i_AAIG_signature + "-" + Peak.j_AAIG_signature == '?-G204':   # DEBUGGING
            # Peak.plot_Peak()  # plot the 2D-hist of this Peak
            self.__replace_peak__(Peak)  # replaces the existing Peak with the one that has the 2D-histogram

    def plot_aliphatic_AAIG(self):
        plt.pcolormesh(self.al_H_edges, self.al_C_edges, np.array(self.al_hist2D.todense()))
        # for c,h in zip(self.scaled_carbons, self.scaled_protons):
        #     print c,h
        plt.scatter(self.scaled_al_protons, self.scaled_al_carbons, s=2, facecolor='white')
        plt.xlabel('scaled aliphatic H chemical shift', fontsize=24);
        plt.ylabel('scaled aliphatic C chemical shift', fontsize=24);
        plt.show()

    def plot_aromatic_AAIG(self):
        if self.ar_hist2D.size == 0:
            print(bcolors.WARNING + "AAIG " + self.signature + " does not contain aromatic peaks!" + bcolors.ENDC)
            return
        plt.pcolormesh(self.ar_H_edges, self.ar_C_edges, np.array(self.ar_hist2D.todense()))
        # for c,h in zip(self.scaled_carbons, self.scaled_protons):
        #     print c,h
        plt.scatter(self.scaled_ar_protons, self.scaled_ar_carbons, s=2, facecolor='white')
        plt.xlabel('scaled aromatic H chemical shift', fontsize=24);
        plt.ylabel('scaled aromatic C chemical shift', fontsize=24);
        plt.show()

    def calc_intersection(self, AAIG2):
        """
        Method to calculate the intersection between the current AAIG and AAIG2 expressed as 2D density histograms.
        :param AAIG2:
        :return: intersection:  a value in the range [0.0, 1.0]
        """
        if not np.array_equal(self.al_H_minmax, AAIG2.al_H_minmax) or \
            not np.array_equal(self.al_C_minmax, AAIG2.al_C_minmax) or \
            not np.array_equal(self.ar_H_minmax, AAIG2.ar_H_minmax) or \
            not np.array_equal(self.ar_C_minmax, AAIG2.ar_C_minmax):
            raise Exception("AAIG1 " + self.signature + " and AAIG2 " + AAIG2.name + " do not have 2D-histograms with "
                                                                        "equal dimensions!\n" +
            " ,".join([str(i) for i in [self.al_H_minmax.tolist() +
                                        self.al_C_minmax.tolist() +
                                        self.ar_H_minmax.tolist() +
                                        self.ar_C_minmax.tolist()]]) +
            " != " + " ,".join([str(j) for j in [AAIG2.al_H_minmax.tolist() +
                                                 AAIG2.al_C_minmax.tolist() +
                                                 AAIG2.ar_H_minmax.tolist() +
                                                 AAIG2.ar_C_minmax.tolist()]]))
        if self.al_hist2D.shape != AAIG2.al_hist2D.shape:
            raise Exception("AAIG1 " + self.signature + " and AAIG2 " + AAIG2.name + " do not have aliphatic 2D-histograms with "
                                                                        "equal dimensions!\n" +
                            " ,".join([str(i) for i in self.al_hist2D.shape]) +
                            " != " + " ,".join([str(j) for j in AAIG2.al_hist2D.shape]))
        if self.ar_hist2D.shape != AAIG2.ar_hist2D.shape and \
                not 0 in self.ar_hist2D.shape + AAIG2.ar_hist2D.shape:
            raise Exception(
                "AAIG1 " + self.signature + " and AAIG2 " + AAIG2.name + " do not have aromatic 2D-histograms with "
                                                                    "equal dimensions!\n" +
                " ,".join([str(i) for i in self.ar_hist2D.shape]) +
                " != " + " ,".join([str(j) for j in AAIG2.ar_hist2D.shape]))

        try:
            nonzero_al_voxels_set1 = set([(i, j) for i, j in zip(*self.al_hist2D.nonzero())])
        except AttributeError:
            raise Exception("AAIG1 " + self.signature + " instance has no attribute 'al_hist2D'")
        try:
            nonzero_al_voxels_set2 = set([(i, j) for i, j in zip(*AAIG2.al_hist2D.nonzero())])
        except AttributeError:
            raise Exception("AAIG2 "+AAIG2.name+" instance has no attribute 'al_hist2D'")
        common_nonzero_al_voxels = nonzero_al_voxels_set1.intersection(nonzero_al_voxels_set2)
        intersection = 0.0
        for i,j in common_nonzero_al_voxels:    # sum up the intersection between densities of aliphatic peaks
            intersection += min(self.al_hist2D[i,j], AAIG2.al_hist2D[i,j])

        if self.ar_hist2D.size > 0 and AAIG2.ar_hist2D.size > 0:    # If at least one AAIG1 contains aromatic peaks
            try:
                nonzero_ar_voxels_set1 = set([(i, j) for i, j in zip(*self.ar_hist2D.nonzero())])
            except AttributeError:
                raise Exception("AAIG1 " + self.signature + " instance has no attribute 'ar_hist2D'")
            try:
                nonzero_ar_voxels_set2 = set([(i, j) for i, j in zip(*AAIG2.ar_hist2D.nonzero())])
            except AttributeError:
                raise Exception("AAIG2 " + AAIG2.name + " instance has no attribute 'ar_hist2D'")
            common_nonzero_ar_voxels = nonzero_ar_voxels_set1.intersection(nonzero_ar_voxels_set2)
            for i, j in common_nonzero_ar_voxels:   # add the intersection between densities of aromatic peaks
                intersection += min(self.ar_hist2D[i, j], AAIG2.ar_hist2D[i, j])

        return intersection



# if __name__ == "__main__":
#
# # EXAMPLE:
# from peak import *
#
#
# N78 = [
# [3.968, 71.015, 1.0],
# [1.067, 21.741, 1.0]
# ]
#
# T77 = [
# [1.494,	17.201,	4041598.0],
# [1.065,	21.548,	16629995.0],
# [1.156,	22.121,	4526034.0],
# [2.375,	31.518,	6262566.0],
# [1.999,	33.474,	5397962.0],
# [4.241,	56.259,	14267978.0],
# [3.955,	71.009,	6147898.0]
# ]
#
# AAIG_dict = {}
# # for peaks, name in zip([Y122, L123, F65, I13, K15, G14, A16, V218], ['Y122', 'L123', 'F65', 'I13', 'K15', 'G14', 'A16', 'V218']):
# for peaks, name in zip([N78, T77], ['N78', 'T77']):
#     print name, peaks
#     aaig = AAIG(name=name,
#                 H_bin_length=0.001,
#                 C_bin_length=0.01,
#                 bandwidth=0.2)
#     for p in peaks:
#         peak = Peak(Creson=p[0], Hreson=p[1], Intensity=p[2])
#         aaig.__add_peak__(peak)
#     aaig.AAIG2hist(even_edges=True)
#     AAIG_dict[name] = aaig
#     print "DEBUG: aaig.al_H_edges.shape=", aaig.al_H_edges.shape
#     print "DEBUG: aaig.al_C_edges.shape=", aaig.al_C_edges.shape
#     print "DEBUG: aaig.al_hist2D.sum()=", aaig.al_hist2D.sum()
#     print "DEBUG: aaig.ar_hist2D.sum()=", aaig.ar_hist2D.sum()
#     aaig.plot_aliphatic_AAIG()
#     # aaig.plot_aromatic_AAIG()
#     # del aaig
#
# # # for name in ['Y122', 'L123', 'F65', 'I13', 'G14', 'K15', 'A16', 'V218']:
# # for name in ['Y122', 'L123']:
# #     print "Intesection between Y122 and "+name+" is", AAIG_dict['Y122'].calc_intersection(AAIG_dict[name]), \
# #         AAIG_dict[name].calc_intersection(AAIG_dict['Y122'])

