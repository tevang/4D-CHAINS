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
import sys
from .global_func import *
from .aaig import *
from .inthist import *
from sklearn.neighbors import KernelDensity
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
import uuid
from decimal import *
getcontext().prec = 7   # 5 decimal points, it applies only for operations between Decimals, also for my Dec objects
# # Example:
# getcontext().prec = 7
# v1 = Decimal(1.456)
# v2 = Decimal(2.621)
# v1
# Decimal('1.455999999999999960920149533194489777088165283203125')
# v2
# Decimal('2.620999999999999996447286321199499070644378662109375')
# v1*v2
# Decimal('3.8162')


class Peak():

    def is_valid_bin_length(self, bin_length):
        """
        The bin_length (either H or C) must be a divisor of 1.0!!!
        :param bin_length:
        :return:
        """
        alpha = multiplier_for_int(bin_length)
        if alpha % int(alpha*bin_length) != 0.0:
            raise Exception(bcolors.FAIL + "Error: bin length value " + str(bin_length) + " is invalid as it is not a divisor of 1.0!"
                            + bcolors.ENDC)

    def __init__(self,
                 i_AAIG_signature=None,
                 i_AAIG_name=None,
                 j_AAIG_signature=None,
                 j_AAIG_name=None,
                 file_label=None,
                 Cname=None,
                 f_Cname=None,  # Cname according to the input file
                 Creson=None,
                 Hname=None,
                 f_Hname=None,
                 Hreson=None,
                 i_Nreson=None,
                 i_Nname=None,
                 fi_Nname=None,
                 i_HNreson=None,
                 i_HNname=None,
                 fi_HNname=None,
                 j_Nreson=None,
                 j_Nname=None,
                 fj_Nname=None,
                 j_HNreson=None,
                 j_HNname=None,
                 fj_HNname=None,
                 Intensity=None,
                 scaled_Intensity=None,
                 H_bin_length=Dec(0.001),
                 C_bin_length=Dec(0.01),
                 bandwidth=0.2):
        self.ID = uuid.uuid4()  # unique identifier for this Peak object
        # Check is bin lengths are valid for 2D-histogram construction
        self.is_valid_bin_length(C_bin_length)
        self.is_valid_bin_length(H_bin_length)

        self.i_AAIG_signature = i_AAIG_signature
        self.i_AAIG_name = i_AAIG_name
        self.j_AAIG_signature = j_AAIG_signature
        self.j_AAIG_name = j_AAIG_name
        self.Cname = Cname
        self.Creson = Creson
        self.Hname = Hname
        self.Hreson = Hreson
        self.i_Nreson = i_Nreson
        self.i_Nname = i_Nname
        self.i_HNreson = i_HNreson
        self.i_HNname = i_HNname
        self.j_Nreson = j_Nreson
        self.j_Nname = j_Nname
        self.j_HNreson = j_HNreson
        self.j_HNname = j_HNname
        self.file_label = file_label    # the full label in the input file
        self.f_Cname = f_Cname  # the atom names according to the input file
        self.f_Hname = f_Hname
        self.fi_Nname = fi_Nname
        self.fi_HNname = fi_HNname
        self.fj_Nname = fj_Nname
        self.fj_HNname = fj_HNname
        self.Intensity = Intensity
        self.scaled_Intensity = scaled_Intensity
        self.H_bin_length = Dec(H_bin_length)
        self.C_bin_length = Dec(C_bin_length)
        self.bandwidth = bandwidth
        self.hist2D = None
        self.H_values = []
        self.C_values = []

        if Intensity == 0:
            print(bcolors.FAIL + "ERROR: the following peak has 0 intensity! Please remove it and run " \
                                           "4D-CHAINS again." + bcolors.ENDC)
            self.__print_peak__()
            sys.exit(1)

    def __print_peak__(self, print_2Dhist=False):
        print("i_AAIG_signature=", self.i_AAIG_signature)
        print("i_AAIG_name=", self.i_AAIG_name)
        print("j_AAIG_signature=", self.j_AAIG_signature)
        print("j_AAIG_name=", self.j_AAIG_name)
        print("Cname=", self.Cname)
        print("Creson=", self.Creson)
        print("Hname=", self.Hname)
        print("Hreson=", self.Hreson)
        print("i_Nreson=", self.i_Nreson)
        print("i_HNreson=", self.i_HNreson)
        print("j_Nreson=", self.j_Nreson)
        print("j_HNreson=", self.j_HNreson)
        print("Intensity=", self.Intensity)
        print("scaled_Intensity=", self.scaled_Intensity)
        print("H_bin_length=", self.H_bin_length)
        print("C_bin_length=", self.C_bin_length)
        print("bandwidth=", self.bandwidth)
        # if type(self.hist2D) is not type(None):     # this works for numpy.array, too
        if print_2Dhist:
            print("hist2D=", self.hist2D.tolist())
            print("H_values=", self.H_values)
            print("C_values=", self.C_values)

    def __print_CHresons__(self):
        print("Hreson=", self.Hreson, "Creson=", self.Creson)

    def __get_CHresons__(self):
        return (self.Hreson, self.Creson)

    def does_peak_match(self, peak2, tolH=0.02, tolC=0.2, tolHN=0.04, tolN=0.4, attypes=['C', 'H', 'N', 'HN', 'Njtoi', 'HNjtoi']):
        """
        Method to tell you if the current Peak object matches with another Peak object given the current tolerances.
        :param peak2:
        :param tolH:
        :param tolC:
        :param tolHN:
        :param tolN:
        :param attypes: 'Njtoi', 'HNjtoi' are special cases for HSQC->HNNH matching. The difference between N from HSQC peak and C from HNNH peak
        as well as between HN from HSQC and H from HNNH peak are calculated. In this scenario the 'self' must always be the HSQC peak!
        :return:
        """
        for a in attypes:
            if a not in ['C', 'H', 'N', 'HN', 'Njtoi', 'HNjtoi']:
                raise ValueError(bcolors.FAIL + "ERROR: atom type " + a + " is not supported!" + bcolors.ENDC)

        match_list = []
        # WARNING: You have to explicitly ask for Hreson!=None because if Hreson==0.0 and
        # WARNING: "if peak2.Hreson and self.Heron and 'H' in aatypes" will never be satisfied even if 'H' in attypes!
        # WARNING: Likewise for the other nuclei.
        if not peak2.Hreson==None and not self.Hreson==None and 'H' in attypes:
            match_list.append( (peak2.Hreson - tolH) <= self.Hreson <= (peak2.Hreson + tolH) )
        if not peak2.Creson==None and not self.Creson==None and 'C' in attypes:
            match_list.append( (peak2.Creson - tolC) <= self.Creson <= (peak2.Creson + tolC) )
        if not peak2.i_HNreson==None and not self.i_HNreson==None and 'HN' in attypes:
            match_list.append((peak2.i_HNreson - tolHN) <= self.i_HNreson <= (peak2.i_HNreson + tolHN))
        if not peak2.i_Nreson==None and not self.i_Nreson==None and 'N' in attypes:
            match_list.append((peak2.i_Nreson - tolN) <= self.i_Nreson <= (peak2.i_Nreson + tolN))
        if not peak2.j_HNreson==None and not self.j_HNreson==None and 'HN' in attypes:
            match_list.append((peak2.j_HNreson - tolHN) <= self.j_HNreson <= (peak2.j_HNreson + tolHN))
        if not peak2.j_Nreson==None and not self.j_Nreson==None and 'N' in attypes:
            match_list.append((peak2.j_Nreson - tolN) <= self.j_Nreson <= (peak2.j_Nreson + tolN))

        if not self.j_HNreson==None and not peak2.i_HNreson==None and 'HNjtoi' in attypes:
            match_list.append((peak2.i_HNreson - tolHN) <= self.j_HNreson <= (peak2.i_HNreson + tolHN))
        if not self.j_Nreson==None and not peak2.i_Nreson==None and 'Njtoi' in attypes:
            match_list.append((peak2.i_Nreson - tolN) <= self.j_Nreson <= (peak2.i_Nreson + tolN))

        return all(match_list) and len(match_list) > 0  # if all conditions are met

    def distance(self, peak2, attypes=['C', 'H', 'N', 'HN', 'Njtoi', 'HNjtoi'], get_distlist=False):
        """
        Method to measure the Euclidian distance from another peak.
        :param peak2:
        :param attypes: 'Njtoi', 'HNjtoi' are special cases for HSQC->HNNH matching. The difference between N from
                        HSQC peak and w1 from HNNH peak as well as between HN from HSQC and w2 from HNNH peak are
                        calculated. In this scenario the 'self' must always be the HSQC peak!
        :param get_distlist:    instead of a single number get a list of individual squared distances of the specified
                                attypes.
        :return:
        """
        dist = []
        if peak2.Hreson and self.Hreson and 'H' in attypes:
            dist.append( (self.Hreson - peak2.Hreson)**2 )
        if peak2.Creson and self.Creson and 'C' in attypes:
            dist.append( ((self.Creson - peak2.Creson) / 6)**2 )
        if peak2.i_HNreson and self.i_HNreson and 'HN' in attypes:
            # print "DEBUG: HN distance: ", self.i_HNreson, '-', peak2.i_HNreson, "=", (self.i_HNreson - peak2.i_HNreson)**2
            dist.append((self.i_HNreson - peak2.i_HNreson) ** 2)
        if peak2.i_Nreson and self.i_Nreson and 'N' in attypes:
            # print "DEBUG: N distance: ", self.i_Nreson, '-', peak2.i_Nreson, "=", ((self.i_Nreson - peak2.i_Nreson) / 6)**2
            dist.append(((self.i_Nreson - peak2.i_Nreson) / 6) ** 2)
        if peak2.j_HNreson and self.j_HNreson and 'HN' in attypes:
            # print "DEBUG: HN distance: ", self.j_HNreson, '-', peak2.j_HNreson, "=", (self.j_HNreson - peak2.j_HNreson)**2
            dist.append((self.j_HNreson - peak2.j_HNreson) ** 2)
        if peak2.j_Nreson and self.j_Nreson and 'N' in attypes:
            # print "DEBUG: N distance: ", self.j_Nreson, '-', peak2.j_Nreson, "=", ((self.j_Nreson - peak2.j_Nreson) / 6)**2
            dist.append(((self.j_Nreson - peak2.j_Nreson) / 6) ** 2)

        if self.j_HNreson and peak2.i_HNreson and 'HNjtoi' in attypes:
            # print "DEBUG: HN distance: ", self.j_HNreson, '-', peak2.j_HNreson, "=", (self.j_HNreson - peak2.j_HNreson)**2
            dist.append((self.j_HNreson - peak2.i_HNreson) ** 2)
        if self.j_Nreson and peak2.i_Nreson and 'Njtoi' in attypes:
            # print "DEBUG: N distance: ", self.j_Nreson, '-', peak2.j_Nreson, "=", ((self.j_Nreson - peak2.j_Nreson) / 6)**2
            dist.append(((self.j_Nreson - peak2.i_Nreson) / 6) ** 2)

        if get_distlist:
            return [np.sqrt(d) for d in dist]

        return np.sqrt(np.sum(dist))

    def even_edges(self, H_minmax, C_minmax):
        """
        Adjusts the extrema of the H or C axis in order the difference of their bins to be 0.
        This method is called before smoothing.
        :return:
        """
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
        binNumH = int((H_minmax[1] - H_minmax[0]) / self.H_bin_length) # + 0j  # convert the bin number in x-axis to complex number to emulate e.g. "100j"
        binNumC = int((C_minmax[1] - C_minmax[0]) / self.C_bin_length) # + 0j  # convert the bin number in y-axis to complex number to emulate e.g. "100j"
        if binNumH != binNumC:
            raise Exception(bcolors.FAIL + "ERROR: the number of bins in H and in C axes are not equal! " +
                            str(binNumH)+" != "+str(binNumC) +
                            " H_minmax= " + str(H_minmax[0]) + ", " + str(H_minmax[1]) +
                            " C_minmax= " + str(C_minmax[0]) + ", " + str(C_minmax[1]) +
                            "\n" + "Use aaig.AAIG2hist(even_edges=True) or set "
                            "the parameters H_bin_length, C_bin_length, H_minmax, C_minmax, ar_H_minmax, "
                            "ar_C_minmax appropriately." + bcolors.ENDC)

        print(bcolors.BOLDGREEN + "H bin length = " + str(self.H_bin_length) + " ppm." + bcolors.ENDBOLD)
        print(bcolors.BOLDGREEN + "C bin length = " + str(self.C_bin_length) + " ppm." + bcolors.ENDBOLD)
        print(bcolors.BOLDGREEN + "2D-histogram dimension = " + str(binNumH) + " x " + str(binNumC) + bcolors.ENDBOLD)

        # hist2D is Nh x Nc array (Proton values in the rows and Carbon values in the columns)
        scaled_H_edges, scaled_C_edges, hist2D = IntHistogram().kde2D(augm_protons,
                                                                      augm_carbons,
                                                                      bandwidth=self.bandwidth,
                                                                      xbins=binNumH,
                                                                      ybins=binNumC,
                                                                      x_minmax=(0., 1.),
                                                                      y_minmax=(0., 1.),  # Carbon values in ppm have 10x larger scale than Proton
                                                                      kernel='cosine',
                                                                      metric='euclidean')
        # Recover the original edges in ppm
        # H_edges = np.array([[H_minmax[0] + i * self.H_bin_length] * binNumH for i in range(binNumC)])
        # C_edges = np.array([np.arange(C_minmax[0], C_minmax[1], self.C_bin_length)] * binNumH)
        # ATTENTION: I encountered one peak in MS6282 where Decimal failed! It produced -0e-60 instead of 0.0 and I couldn't
        # ATTENTION: align its 2D-histogram! It was ?-S35 (23.812, 0.038).
        # print("DEBUG: H_minmax=", H_minmax, "self.H_bin_length=", self.H_bin_length, "binNumH=", binNumH)
        # print("DEBUG: C_minmax=", C_minmax, "self.C_bin_length=", self.C_bin_length, "binNumC=", binNumC)
        H_values = [H_minmax[0] + i * self.H_bin_length for i in range(binNumH)]
        C_values = [C_minmax[0] + i * self.C_bin_length for i in range(binNumC)]
        # H_values = [round(h, 4) for h in H_values]  # don't use Decimal!! # TEMPORARILY DEACTIVATED
        # C_values = [round(c, 4) for c in C_values]    # TEMPORARILY DEACTIVATED

        return scaled_H_edges, scaled_C_edges, hist2D, H_values, C_values

    def peak2hist(self, even_edges=True, scale_density=False):
        """
        Generated the 2D density histogram of this Peak.

        :param even_edges:
        :param scale_density:
        :return:
        """

        ColorPrint("Generating the 2D-histogram of peak " + self.i_AAIG_name + "-" + self.j_AAIG_signature + \
              " (" + str(self.Creson) + ", " + str(self.Hreson) + ").", "OKBLUE")
        # Create a local copy of Creson and Hreson and convert them to Decimals for consistence with the calculations
        Creson = Dec(self.Creson)
        Hreson = Dec(self.Hreson)

        # Find the minimum and maximum H, C values of the 2D-histogram
        if Hreson >= 0.0:
            grid = np.arange(int(Hreson), int(Hreson) + 1 + self.H_bin_length, self.H_bin_length)
        elif Hreson < 0.0:
            grid = np.arange(int(Hreson) - 1, int(Hreson) + self.H_bin_length, self.H_bin_length)
        grid = np.array([Dec(v) for v in grid], dtype=Decimal)   # convert all values to fixed decimal length floats! (NECESSARY))
        ind = find_nearest_index(grid, Dec(Hreson))    # find the edge in the axis that is nearest to the Hreson
        edge = Dec(grid[ind])
        H_subgrid = np.arange(edge - Dec(0.04), edge + Dec(0.04) + self.H_bin_length, self.H_bin_length)
        # H_subgrid = [round(v, 4) for v in H_subgrid]  # convert all values to fixed decimal length floats! (NECESSARY) # TEMPORARILY DEACTIVATED
        H_minmax = [Dec(H_subgrid[0]), Dec(H_subgrid[-1])]

        if Creson >= 0.0:
            grid = Decarange(int(Creson), Dec(int(Creson) + 1) + Dec(self.C_bin_length), self.C_bin_length)
        elif Creson < 0.0:
            grid = np.arange(int(Creson) - 1, int(Creson) + self.C_bin_length, self.C_bin_length)
        grid = np.array([Dec(v) for v in grid], dtype=Decimal)  # convert all values to fixed decimal length floats! (NECESSARY)))
        ind = find_nearest_index(grid, Dec(Creson))
        edge = Dec(grid[ind])
        C_subgrid = np.arange(edge - Dec(0.4), edge + Dec(0.4) + self.C_bin_length, self.C_bin_length)
        # C_subgrid = [round(v, 4) for v in C_subgrid]  # convert all values to fixed decimal length floats! (NECESSARY) # TEMPORARILY DEACTIVATED
        C_minmax = [Dec(C_subgrid[0]), Dec(C_subgrid[-1])]

        if even_edges:
            aaig = AAIG(H_bin_length=self.H_bin_length, C_bin_length=self.C_bin_length)
            H_minmax, C_minmax = aaig.even_edges(H_minmax=H_minmax, C_minmax=C_minmax)

        # IMPORTANT: Scale the proton and carbon resonances between 0. and 1., otherwise smoothing won't work equally in both dimensions!!!
        self.scaled_Creson = minmax_scale(Creson, minmax=C_minmax)
        self.scaled_Hreson = minmax_scale(Hreson, minmax=H_minmax)
        # include the intensity implicitly by augmenting the respective H and C resonances accordingly
        augm_protons, augm_carbons = IntHistogram().augment_single_peak(self.scaled_Hreson, self.scaled_Creson)
                                                                            # np.array([self.scaled_Intensity]))

        self.scaled_H_edges, \
        self.scaled_C_edges, \
        self.hist2D, \
        self.H_values, \
        self.C_values = self.create_2Dhist(augm_protons,
                                          augm_carbons,
                                          H_minmax,
                                          C_minmax,
                                          even_edges=even_edges)


        # print "DEBUG: Peak (%f, %f)\n" % (Hreson, Creson)
        # print "DEBUG: hist2D=", self.hist2D.tolist()
        # print "DEBUG: H_values=", self.H_values
        # print "DEBUG: C_values=", self.C_values
        # self.hist2D = csr_matrix(hist2D)

        # Scale the densities to sum to 1.0 in the 2D-histogram
        if self.hist2D.sum() == 0.0:    # instead of scaling densities use this
            raise Exception(bcolors.FAIL + "ERROR: the bandwidth value (" + str(self.bandwidth) +" is too low to create density for"
                                        " Peak " + self.i_AAIG_name + "-" + self.j_AAIG_signature +
                            " ! Try to increase it and rerun the program." + bcolors.ENDC)
        # print "DEBUG: self.hist2D.sum()=", self.hist2D.sum()
        # IMPORTANT: smoothing seems to produce the same density irrespective of the number of augmented peak list.
        # IMPORTANT: therefore since we have only one peak in the 2D-histogram and solution would be to scale the density to 0-1.0
        # IMPORTANT: and multiply it by the scaled peak intensity.
        if not self.scaled_Intensity:
            raise Exception(bcolors.FAIL + "ERROR: you forgot to scale the intensities with "
                                           "Connectivities().__scale_intensities__(wrt_aaig=False) function!" + bcolors.ENDC)
        self.hist2D = self.scaled_Intensity * self.hist2D/float(self.hist2D.sum())

    def plot_Peak(self):
        plt.pcolormesh(self.scaled_H_edges, self.scaled_C_edges, np.array(self.hist2D))
        # for c,h in zip(self.scaled_carbons, self.scaled_protons):
        #     print c,h
        plt.scatter([self.scaled_Hreson], [self.scaled_Creson], s=2, facecolor='white')
        plt.xlabel('scaled aliphatic H chemical shift', fontsize=24);
        plt.ylabel('scaled aliphatic C chemical shift', fontsize=24);
        plt.show()

    def calc_intersection(self, peak2):
        """
        Method to calculate the intersection between two Peaks expressed as 2D density histograms.

        :param peak2:
        :return: intersection:  a value in the range [0.0, 1.0]
        """
        name1 = self.i_AAIG_name + "-" + self.j_AAIG_signature + " (" + str(self.Creson) + ", " + str(self.Hreson) + ")"
        name2 = peak2.i_AAIG_name + "-" + peak2.j_AAIG_signature + " (" + str(peak2.Creson) + ", " + str(peak2.Hreson) + ")"

        try:
            H_min1, H_max1 = self.H_values[0], self.H_values[-1]
        except AttributeError:
            raise AttributeError(bcolors.FAIL + "Peak " + name1 + " has not attribute 'H_values'!" + bcolors.ENDC)
        try:
            H_min2, H_max2 = peak2.H_values[0], peak2.H_values[-1]
        except AttributeError:
            raise AttributeError(bcolors.FAIL + "Peak " + name2 + " has not attribute 'H_values'!" + bcolors.ENDC)
        C_min1, C_max1 = self.C_values[0], self.C_values[-1]
        C_min2, C_max2 = peak2.C_values[0], peak2.C_values[-1]

        try:
            if H_min1 <= H_min2:
                row_start1 = self.H_values.index(H_min2)
                row_start2 = 0
                row_end1 = len(self.H_values)               # the ending index must be +1 to return also the last element
                row_end2 = peak2.H_values.index(H_max1) + 1 # the ending index must be +1 to return also the last element
            else:
                row_start1 = 0
                row_start2 = peak2.H_values.index(H_min1)
                row_end1 = self.H_values.index(H_max2) + 1  # the ending index must be +1 to return also the last element
                row_end2 = len(self.H_values)               # the ending index must be +1 to return also the last element

            if C_min1 <= C_min2:
                col_start1 = self.C_values.index(C_min2)
                col_start2 = 0
                col_end1 = len(self.C_values)               # the ending index must be +1 to return also the last element
                col_end2 = peak2.C_values.index(C_max1) + 1 # the ending index must be +1 to return also the last element
            else:
                col_start1 = 0
                col_start2 = peak2.C_values.index(C_min1)
                col_end1 = self.C_values.index(C_max2) + 1  # the ending index must be +1 to return also the last element
                col_end2 = len(self.C_values)               # the ending index must be +1 to return also the last element
        except ValueError:
            raise ValueError(bcolors.FAIL + "ERROR: the 2D-histograms of the peaks " + name1 + " and " + name2 +
                             " do not overlap!!!\n" +
                             "H_values1=" + ",".join([str(v) for v in self.H_values]) + "\n" +
                             "H_values2=" + ",".join([str(v) for v in peak2.H_values]) + "\n" +
                             "C_values1=" + ",".join([str(v) for v in self.C_values]) + "\n" +
                             "C_values2=" + ",".join([str(v) for v in peak2.C_values]) + "\n" + bcolors.ENDC
                             )

        try:
            intersection = np.minimum(self.hist2D[row_start1:row_end1, col_start1:col_end1],
                       peak2.hist2D[row_start2:row_end2, col_start2:col_end2]).sum()
        except AttributeError:
            raise Exception("Peak " + name1 + " or " + name2 + " instance has no attribute 'hist2D'")
        return intersection


# if __name__ == "__main__":
#     N78 = [
#         [3.968, 71.015, 1.0],
#         [1.067, 21.741, 1.0]
#     ]
#
#     T77 = [
#         [3.955, 71.009, 6147898.0/832898816.0],
#         [1.065, 21.548, 16629995.0/832898816.0]
#     ]
#
#     for i in range(2):
#         peak1 = Peak(i_AAIG_name="N78", j_AAIG_signature="N78", Hreson=N78[i][0], Creson=N78[i][1], scaled_Intensity=N78[i][2],
#                      H_bin_length=0.005, C_bin_length=0.05)
#         peak2 = Peak(i_AAIG_name="T77", j_AAIG_signature="T77", Hreson=T77[i][0], Creson=T77[i][1], scaled_Intensity=T77[i][2],
#                      H_bin_length=0.005, C_bin_length=0.05)
#         peak1.peak2hist()
#         peak2.peak2hist()
#         peak1.__print_peak__()
#         peak2.__print_peak__()
#         print "Intersection =", peak1.calc_intersection(peak2)
#         break