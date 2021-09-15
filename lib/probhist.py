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
from decimal import Decimal, getcontext
# getcontext().prec = 7   # 6 decimal points
CHAINS_BIN_LIB = os.path.dirname(os.path.realpath(__file__))
from .global_vars import *
from .global_func import *
from collections import Counter
from collections import Mapping

class ProbHist_Loader():

    def __init__(self, PROBHIST_DIR=""):
        if not PROBHIST_DIR:
            self.PROBHIST_DIR = CHAINS_BIN_LIB + "/../databases/histograms/"
        else:
            self.PROBHIST_DIR = PROBHIST_DIR
        self.aa_carbon_binDensityList_mdict = tree()  # multidict with amino acid name -> Carbon atom type -> [array of bin limits, array of probability density]
        self.aa_hydrogen_binDensityList_mdict = tree()  # multidict with amino acid name -> Hydrogen atom type -> [array of bin limits, array of probability density]

        self.aa_CHpair_binProbabilityList_mdict = tree()    # multidict with amino acid name -> Carbon nad Proton atom type ->
                                                                # [array of x bin limits, array of y bin limits, array of probability density]
    def load_1Dhistograms(self):
        ## LOAD 1D HISTOGRAMS
        print(bcolors.BOLDBLUE +  "\nLoading chemical shift 1D-histograms..." + bcolors.ENDBOLD)
        fnames = os.listdir(self.PROBHIST_DIR)
        fpattern = re.compile("[A-Z]{3}_[A-Z0-9]+_hist.txt$")
        hist_files_list = list(filter(fpattern.search, fnames))
        for hist_file in hist_files_list:
            aa = hist_file.split("_")[0]
            atom = hist_file.split("_")[1]
            # try:
            atoms2load = [a for a in allowed_aa_atoms_dict[aa]]
            if aa == "PHE":
                atoms2load.extend(["HD1", "CD1"])
            elif aa == "TYR":
                atoms2load.extend(["HD1", "HE1", "CD1", "CE1"])
            atoms2load.extend([])
            if atom in atoms2load:  # load histograms of the allowed atom types only
                bin_list = []
                density_list = []
                with open(self.PROBHIST_DIR + hist_file, 'r') as f:
                    for line in f:
                        word_list = line.split()
                        bin_list.append(float(word_list[0]))
                        density_list.append(float(word_list[1]))
                bin_array = np.array(bin_list)
                density_array = np.array(density_list) / sum(density_list)
                # use C & H resonances to assign aa type to the (i-1) residue
                if atom[0] == "C":
                    self.aa_carbon_binDensityList_mdict[aa][atom] = [bin_array, density_array]
                elif atom[0] == "H":
                    self.aa_hydrogen_binDensityList_mdict[aa][atom] = [bin_array, density_array]
            # except KeyError:   # active this exception if use N, HN resonance to assign aa type to the ith residue
            #    if aa == "PRO": # Prolines cannot be detected in root spectrum (N-H HSQC) because they don't have "HN" atom
            #        continue

    def get_step(self, bin_array):
        srt = list(set(bin_array))
        srt.sort()
        return Dec(srt[1] - srt[0])
        # return float(srt[1] - srt[0])

    def get_row_length(self, bin_array):
        counts_dict = Counter(bin_array)
        return next(iter(list(counts_dict.values())))

    def map_nested_dicts(self, ob, func):
        if isinstance(ob, Mapping):
            return {k: self.map_nested_dicts(v, func) for k, v in list(ob.items())}
        else:
            return func(ob)

    def load_2Dhistograms(self):
        """
        In general, when loading 2D-histograms all calculations are significantly slown down!
        :return:
        """

        ## LOAD 2D HISTOGRAMS
        print(bcolors.BOLDBLUE + "\nLoading chemical shift 2D-histograms..." + bcolors.ENDBOLD)
        fnames = os.listdir(self.PROBHIST_DIR)
        fpattern = re.compile("[A-Z]{3}_[A-Z0-9-]+_correlated_2Dhist.smoothed.txt$")
        hist_files_list = list(filter(fpattern.search, fnames))
        for hist_file in hist_files_list:
            # print "Loading histogram", hist_file
            aa = hist_file.split("_")[0]
            CH_pair = hist_file.split("_")[1]
            # try:
            if CH_pair in list(aa_CHpair_2Dhist_mdict[aa].keys()):  # load histograms of the allowed atom types only
                x_bin_list = []  # for H
                y_bin_list = []  # for C
                probability_list = []
                with open(self.PROBHIST_DIR + hist_file, 'r') as f:
                    for line in f:
                        word_list = line.split()
                        # x_bin_list.append(Dec(word_list[0]))  # for the new 2Dhist prob function
                        # y_bin_list.append(Dec(word_list[1]))  # for the new 2Dhist prob function
                        x_bin_list.append(float(word_list[0]))  # for the old 2Dhist prob function
                        y_bin_list.append(float(word_list[1]))  # for the old 2Dhist prob function
                        probability_list.append(float(word_list[2]))
                x_bin_array = np.array(x_bin_list)
                y_bin_array = np.array(y_bin_list)
                probability_array = np.array(probability_list)
                # use C & H resonances to assign aa type to the (i-1) residue
                self.aa_CHpair_binProbabilityList_mdict[aa][CH_pair] = [x_bin_array, y_bin_array, probability_array]
            # except KeyError:   # active this exception if use N, HN resonance to assign aa type to the ith residue
            #    if aa == "PRO": # Prolines cannot be detected in root spectrum (N-H HSQC) because they don't have "HN" atom
            #        continue

        # # The multidicts below are used only by the NEW FUNCTION to find the probability from a 2D histogram
        # self.minmax_mdict = self.map_nested_dicts(self.aa_CHpair_binProbabilityList_mdict,
        #                                     lambda v: (v[0].min(), v[0].max(), v[1].min(), v[1].max()))
        #
        # self.step_mdict = self.map_nested_dicts(self.aa_CHpair_binProbabilityList_mdict,
        #                                   lambda v: (self.get_step(v[0]), self.get_step(v[1])))
        #
        # self.row_length_mdict = self.map_nested_dicts(self.aa_CHpair_binProbabilityList_mdict, lambda v: self.get_row_length(v[0]))


class ProbHist_Calculator():

    def __init__(self):
        pass

    def fix_decimal_points(self, prob, decimal_places=8):
        if prob < 10e-40:
            return 0.0
        product = 0.0
        exponent = 1
        while product < 1.0:
            product = prob * 10 ** exponent
            exponent += 1

        return float(("%." + str(decimal_places) + "fe-" + str(exponent)) % (prob * 10 ** exponent))

    # def get_probability_from_2Dhistogram(self,
    #                                      CH_pair,
    #                                      aa,
    #                                      Hreson,
    #                                      Creson):
    #     """
    #         NEW FUNCTION to find the probability from a 2D histogram. New faster version, knows where to start searching.
    #     """
    #     _, _, probability_array = self.aa_CHpair_binProbabilityList_mdict[aa][CH_pair]
    #     x_min, x_max, y_min, y_max = self.minmax_mdict[aa][CH_pair]
    #     x_step, y_step = self.step_mdict[aa][CH_pair]
    #     row_length = self.row_length_mdict[aa][CH_pair]
    #
    #     Hreson = Dec(Hreson)
    #     Creson = Dec(Creson)
    #     # Hreson = float(Hreson)
    #     # Creson = float(Creson)))
    #
    #     # If out of borders of the 2D-hist, return 0 probability
    #     if (Hreson > x_max or Hreson < x_min or
    #             Creson > y_max or Creson < y_min):
    #         return 0.0
    #
    #     x = int((Hreson - x_min) / x_step)
    #     y = int((Creson - y_min) / y_step)
    #
    #     index = x * row_length + y
    #
    #     try:
    #         probability = self.fix_decimal_points(float(probability_array[index]))  # I checked it and works correctly !
    #         # print "Converted ", probability_array[index], " to ", probability
    #     except OverflowError:
    #         # print "OverflowError when converting ", probability_array[index], "to float!"
    #         sys.exit(1)
    #
    #     return probability

    def get_probability_from_histogram(self, resonance, bin_array, density_array):
        """
        FUNCTION to get the probability of a resonance from the respective 1D histogram.

    	EXAMPLE USAGE:
    	probablity = get_probability_from_histogram(56.4, aa_carbon_binDensityList_mdict["MET"]["CA"][0], aa_carbon_binDensityList_mdict["MET"]["CA"][1])
        """
        previous_hist_bin = bin_array[0]
        resonance = float(resonance)
        probability = 0.0  # in case the resonance falls completely out of the borders of the histogram
        for index, hist_bin in enumerate(np.nditer(bin_array[1:],
                                                   order='K')):  # start counting from the 2nd element because we have appended "0" at the beginning of bin_list
            # print "DEBUG: previous_hist_bin=", previous_hist_bin, "resonance=", resonance,"hist_bin=",hist_bin
            if previous_hist_bin <= resonance and resonance < hist_bin:
                # print "DEBUG: index=",index, "previous_hist_bin=",previous_hist_bin,"resonance=",resonance,"hist_bin=",hist_bin
                probability = float(density_array[index])  # checked it and works correctly !
                break
            previous_hist_bin = hist_bin
        return probability

    def get_H_probability_from_1Dhistogram(self, aa, Hname, Hreson, histload):
        bin_array, density_array = histload.aa_hydrogen_binDensityList_mdict[aa][Hname]
        return self.get_probability_from_histogram(Hreson, bin_array, density_array)

    def get_C_probability_from_1Dhistogram(self, aa, Cname, Creson, histload):
        bin_array, density_array = histload.aa_carbon_binDensityList_mdict[aa][Cname]
        return self.get_probability_from_histogram(Creson, bin_array, density_array)

    def get_CH_probability_from_1Dhistograms(self,
                                             aa,
                                             Hname,
                                             Hreson,
                                             Cname,
                                             Creson,
                                             histload,
                                             PROBABILITY_MODEL,
                                             H_weight,
                                             C_weight):
        """
        FUNCTION to find the probability from H and C 1D-histograms. This functions implements the get_probability_from_histogram()
            to load 1D histograms and provides more functionalities. It is used only once within get_probabilities_from_H_C_resonpair_2Dhist().

        :param aa:
        :param Hname:
        :param Hreson:
        :param Cname:
        :param Creson:
        :param PROBABILITY_MODEL:
        :param H_weight:
        :param C_weight:
        :return:
        """

        C_bin_array, C_density_array = histload.aa_carbon_binDensityList_mdict[aa][Cname]
        H_bin_array, H_density_array = histload.aa_hydrogen_binDensityList_mdict[aa][Hname]

        H_prob = self.get_probability_from_histogram(Hreson, H_bin_array, H_density_array)
        C_prob = self.get_probability_from_histogram(Creson, C_bin_array, C_density_array)

        if H_prob > 0.0:
            if PROBABILITY_MODEL == 1:  # treat the H and C probabilities as dependent events
                # print "DEBUG: treating the H and C probabilities as dependent events."
                if aa == "LEU":
                    probability = (1.0 * H_prob + 0.1 * C_prob) / float((1.0 + 0.1))
                else:
                    probability = (H_weight * H_prob + C_weight * C_prob) / float((H_weight + C_weight))
            elif PROBABILITY_MODEL == 2:
                # print  "DEBUG: treating the H and C probabilities as independent events."
                probability = H_prob * C_prob
        elif H_prob == 0.0:  # if hist(H)=0, consider only hist(C) but make it negative to distiguish it from weighted average probabilities
            probability = -1 * C_prob  # temporarily deactivate to see the outcome

        probability = self.fix_decimal_points(probability)
        # print "DEBUG: 1D hist probability=", probability

        return probability

    def get_CH_probability_from_2Dhistogram(self, CH_pair, aa, Hreson, Creson, histload):
        """
                OLD FUNCTION to find the probability from a 2D histogram. New faster version, knows where to start searching.

        If you want to validate this function:
        awk '{if(NR%500==0){print $1+0.01,$2+0.1, $3}}' 2D_hist/LEU_CB-HB2_correlated_2Dhist.smoothed.txt > test_LEU_CB_HB2.txt
        # and then in ipython
        from data_augmentation import aa_CHpair_binProbabilityList_mdict, get_probability_from_2Dhistogram
        fname="test_LEU_CB_HB2.txt"
        trash, aa, C, H = fname.replace(".txt", "").split('_')
        x_bin_array, y_bin_array, probability_array = aa_CHpair_binProbabilityList_mdict[aa][C+"-"+H]
        with open(fname, 'r') as f:
            for line in f:
                Hreson = float(line.split()[0])
                Creson = float(line.split()[1])
                corr_prob = float(line.split()[2])
                prob = get_probability_from_2Dhistogram(Hreson, Creson, x_bin_array, y_bin_array, probability_array)
                print prob, "=", corr_prob
        :param CH_pair:
        :param aa:
        :param Hreson:
        :param Creson:
        :param histload:
        :return:
        """
        x_bin_array, y_bin_array, probability_array = histload.aa_CHpair_binProbabilityList_mdict[aa][CH_pair]

        # If out of borders of the 2D-hist, return 0 probability
        if (Hreson > x_bin_array.max() or Hreson < x_bin_array.min() or
                Creson > y_bin_array.max() or Creson < y_bin_array.min()):
            return 0.0

        x_min = round(x_bin_array[0], 2)
        x_offset = x_bin_array[x_bin_array == x_min].shape[0]
        x_bin_list = list(set(x_bin_array[0:2 * x_offset]))
        x_bin_list.sort()
        x_bin_length = round(x_bin_list[1] - x_bin_list[0], 2)

        y_bin_length = round(y_bin_array[1] - y_bin_array[0], 2)

        # Now predict approximately the location of the peak in the 2D hist (both x_bin_array & y_bin_array have the same length)
        start = int((round(Hreson - x_min, 2) / x_bin_length) * x_offset)

        # print "DEBUG: Hreson=", Hreson, "Creson=", Creson
        # print "DEBUG: x_bin_array=", x_bin_array.tolist()
        # print "DEBUG: y_bin_array=", y_bin_array.tolist()
        # print "DEBUG: probability_array=", probability_array.tolist()
        index = start - x_offset  # start measuring from the starting index
        # print "DEBUG: start=", start, "offset=", x_offset, "index=", index, x_bin_array.shape, y_bin_array.shape
        probability = 0.0
        for x_hist_bin, y_hist_bin in zip(np.nditer(x_bin_array[start - x_offset:], order='K'),
                                          np.nditer(y_bin_array[start - x_offset:], order='K')):
            # print "DEBUG: x_hist_bin=", x_hist_bin, "Hreson=", Hreson,"x_hist_bin+x_bin_length=",x_hist_bin+x_bin_length
            # print "DEBUG: y_hist_bin=", y_hist_bin, "Creson=", Creson,"y_hist_binyx_bin_length=",y_hist_bin+y_bin_length
            if x_hist_bin <= Hreson and Hreson < x_hist_bin + x_bin_length and y_hist_bin <= Creson and Creson < y_hist_bin + y_bin_length:
                # print "DEBUG: index=",index, "x_hist_bin=", x_hist_bin, "Hreson=", Hreson,"x_hist_bin+x_bin_length=",x_hist_bin+x_bin_length
                # print "DEBUG: index=",index, "y_hist_bin=", y_hist_bin, "Creson=", Creson,"y_hist_bin+y_bin_length=",y_hist_bin+y_bin_length
                try:
                    probability = self.fix_decimal_points(float(probability_array[index]))  # I checked it and works correctly !
                    # print "Converted ", probability_array[index], " to ", probability
                except OverflowError:
                    # print "OverflowError when converting ", probability_array[index], "to float!"
                    sys.exit(1)
                break
            index += 1
        # print "DEBUG: 2D hist probability=", probability
        return probability

    def smoothen_distribution(self, proton, carbon, Hedges, Cedges):
        """
        proton: hydrogen chemical shift array
        carbon: carbon chemical shift array. Both proton & carbon must have the same length.
        RETURN:
        Z0:     array with the smoothed probability distribution.
        """
        xmin = proton.min()
        xmax = proton.max()
        ymin = carbon.min()
        ymax = carbon.max()
        # X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
        binNox = len(Hedges) - 1 + 0j  # convert the bin number in x-axis to complex number to emulate e.g. "100j"
        binNoy = len(Cedges) - 1 + 0j  # convert the bin number in y-axis to complex number to emulate e.g. "100j"
        X, Y = np.mgrid[xmin:xmax:binNox, ymin:ymax:binNoy]
        positions = np.vstack(
            [X.ravel(), Y.ravel()])  # stack the two arrays vertically (one in each row); creates a 2x10000 array
        values = np.vstack([proton,
                            carbon])  # place the proton observations on the 1st row and the carbon observations of the 2nd; creates a 2xN array (N is the number of observations)
        kernel = stats.gaussian_kde(values)  # apply
        Z = np.reshape(kernel(positions).T,
                       X.shape)  # convert the 1D array with 10,000 values to a 2D array with 100x100 values
        Z0 = Z / Z.sum()  # normalize the distribution (sum of all values must be 1, like in a 2D probability distribution)

        return Z0

        # # NOW PLOT IT
        # import matplotlib.pyplot as plt
        #
        # fig = plt.figure()
        # ax = fig.add_subplot(111)
        # ax.imshow(np.rot90(Z0), cmap=plt.cm.gist_earth_r, extent=[xmin, xmax, ymin, ymax])
        # ax.plot(proton, carbon, 'k.', markersize=2)
        # ax.set_xlim([xmin, xmax])
        # ax.set_ylim([ymin, ymax])
        # plt.show()