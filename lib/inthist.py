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


import os, re, scipy
import numpy as np
from sklearn.neighbors import KernelDensity
from .global_vars import *
from .global_func import *

class IntHistogram():

    def __init__(self):
        pass

    def kde2D(self,
              x,
              y,
              bandwidth,
              xbins=100,
              ybins=100,
              x_minmax=(0., 1.),
              y_minmax=(0., 1.),
              **kwargs):
        """
         Build 2D kernel density estimate (KDE).

        :param x:
        :param y:
        :param bandwidth:
        :param xbins:
        :param ybins:
        :param x_minmax:
        :param y_minmax:
        :param kwargs:
        :return: zz:    has the same dimensions as xx (Nx X Ny).
        """

        # # create grid of sample locations (default: 100x100)
        # xx, yy = np.mgrid[x_minmax[0]:x_minmax[1]:xbins, y_minmax[0]:y_minmax[1]:ybins]

        x_bin_length = (x_minmax[1] - x_minmax[0]) / float(xbins)
        y_bin_length = (y_minmax[1] - y_minmax[0]) / float(ybins)

        xx = np.array([[x_minmax[0] + i * x_bin_length] * xbins for i in range(ybins)])
        yy = np.array([np.arange(y_minmax[0], y_minmax[1], y_bin_length)] * xbins)

        xy_sample = np.vstack([yy.ravel(), xx.ravel()]).T
        xy_train = np.vstack([y, x]).T

        kde_skl = KernelDensity(bandwidth=bandwidth, **kwargs)
        kde_skl.fit(xy_train)

        # score_samples() returns the log-likelihood of the samples
        z = np.exp(kde_skl.score_samples(xy_sample))
        zz = np.reshape(z, xx.shape)    # the smoothed 2D histogram of peak intensity

        return xx, yy, zz

    def do_kdtree(self, combined_x_y_arrays, points):
        """
        Find position of points in flatten meshgrid.
        :param combined_x_y_arrays:
        :param points:
        :return:
        """
        mytree = scipy.spatial.cKDTree(combined_x_y_arrays)
        dist, indexes = mytree.query(points)
        return indexes

    def augment_single_peak(self, Hreson, Creson):
        """
        No intensity is needed. Works only for 2D-histograms for single peaks.
        :param Hreson:
        :param Creson:
        :return:
        """
        alpha = 1000    # it doesn't make a difference if you increase it
        return np.array([Hreson]*alpha), np.array([Creson]*alpha)

    def augment_peaks_by_intensity(self, protons, carbons, intensities):
        """
        Converts the scaled intensity of each peak to an integer I and then replicates (augments) each H,C resonance
        pairs of each peak xI times in order that point on the 2D-hist to have higher cosine hill.
        :param protons:
        :param carbons:
        :param intensities:
        :return:
        """
        if intensities.size == 0:
            return np.array([]), np.array([])
        # if not min_intensity:
        #     min_intensity = min(intensities)
        # FIX THE DECIMAL POINTS
        decimals = 4    # according to 9 proteins, the highest ratio max/min intensity is 8048. With decimals=8 it runs out of memory!
        alpha = 10 ** decimals
        # while int(alpha*min_intensity) == 0.0:    # if don't want 0 intensity because that peak will have no density in the 2D-hist!
        #     decimals += 1
        #     alpha = 10**decimals
        #     if decimals > 5:
        #         i = intensities.tolist().index(0.0)
        #         raise Exception("There was one peak with intensity "+str(intensities[i])+"! Please remove this peak "
        #                         "and run the program again!\n"
        #                         "peak = ("+str(carbons[i])+", "+str(protons[i])+", "+str(intensities[i])+")")
        # print "DEBUG: intensities=", intensities.tolist()
        # print "DEBUG: protons=", protons
        # print "DEBUG: carbons=", carbons
        int_intensities = [int(alpha*i) for i in intensities]
        augm_protons, augm_carbons = [], []
        for i, h, c in zip(int_intensities, protons, carbons):
            augm_protons.extend([h] * i)
            augm_carbons.extend([c] * i)
        return np.array(augm_protons), np.array(augm_carbons)

def get_peak_position_in_hist(H_resonance, C_resonance, x_bin_array, y_bin_array):
    # If out of borders of the 2D-hist, return 0 probability
    if (H_resonance > x_bin_array.max() or H_resonance < x_bin_array.min() or
            C_resonance > y_bin_array.max() or C_resonance < y_bin_array.min()):
        raise Exception

    x_min = round(x_bin_array[0], 2)
    x_offset = x_bin_array[x_bin_array == x_min].shape[0]
    x_bin_list = list(set(x_bin_array[0:2 * x_offset]))
    x_bin_list.sort()
    x_bin_length = round(x_bin_list[1] - x_bin_list[0], 2)

    y_bin_length = round(y_bin_array[1] - y_bin_array[0], 2)

    # Now predict approximately the location of the peak in the 2D hist (both x_bin_array & y_bin_array have the same length)
    start = int((round(H_resonance - x_min, 2) / x_bin_length) * x_offset)

    # print "DEBUG: H_resonance=", H_resonance, "C_resonance=", C_resonance
    # print "DEBUG: x_bin_array=", x_bin_array.tolist()
    # print "DEBUG: y_bin_array=", y_bin_array.tolist()
    # print "DEBUG: probability_array=", probability_array.tolist()
    index = start - x_offset  # start measuring from the starting index
    # print "DEBUG: start=", start, "offset=", x_offset, "index=", index, x_bin_array.shape, y_bin_array.shape
    probability = 0.0
    for x_hist_bin, y_hist_bin in zip(np.nditer(x_bin_array[start - x_offset:], order='K'),
                                      np.nditer(y_bin_array[start - x_offset:], order='K')):
        # print "DEBUG: x_hist_bin=", x_hist_bin, "H_resonance=", H_resonance,"x_hist_bin+x_bin_length=",x_hist_bin+x_bin_length
        # print "DEBUG: y_hist_bin=", y_hist_bin, "C_resonance=", C_resonance,"y_hist_binyx_bin_length=",y_hist_bin+y_bin_length
        if x_hist_bin <= H_resonance < x_hist_bin + x_bin_length and y_hist_bin <= C_resonance < y_hist_bin + y_bin_length:
            break
        index += 1

    return index

