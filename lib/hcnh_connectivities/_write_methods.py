import os
from copy import deepcopy
from operator import itemgetter

from lib.global_func import bcolors, save_pickle
from lib.open_func import bcolors


def write_spectrum_info(self):
    if not self.HCNH_spec:
        raise Exception(bcolors.FAIL + "ERROR: you must load the NOESY spectrum in order to write statistics "
                                       "in file!" + bcolors.ENDC)
    with open("spectrum_statistics.txt", 'w') as f:
        AAIG_num, peak_num, ave_peaknum_per_aaig, min_peaknum_per_aaig, max_peaknum_per_aaig = \
            self.HCNH_spec.__get_spectrum_info__()
        f.write("\nNOESY:\n")
        f.write("AAIG number: %g\n" % AAIG_num)
        f.write("Peak number: %g\n" % peak_num)
        f.write("Average peak number per AAIG: %f\n" % ave_peaknum_per_aaig)
        f.write("Minimum peak number per AAIG: %i\n" % min_peaknum_per_aaig)
        f.write("Maximum peak number per AAIG: %i\n" % max_peaknum_per_aaig)


def write_NOESY_NOESY_connectivities(self,
                                     fname,
                                     pickle_fname="",
                                     oexp=1.6497056483276888,
                                     selected_AAIG_connectivities_dict={}):
    """
    OLD oexp=3.0448217240563293, OPTIMIZED FOR TOCSY-NOESY.

    Method to write the connectivity present in self.NOESY_NOESY_peak_connectivities_mdict into a file.
    If you want to pass connectivities from an external dictionary use method write_HNNH_connectivities_from_dict().
    If the selected_AAIG_connectivities_dict is provided, then the method will still write all connectivities in
    self.NOESY_NOESY_peak_connectivities_mdict but when an AAIG is a key in both dicts then only the connectivities
    in selected_AAIG_connectivities_dict will be written. Therefore, if self.NOESY_NOESY_peak_connectivities_mdict
    is empty or contains fewer connectivities than selected_AAIG_connectivities_dict, then this method is useless and
    dangerous!!!

    :param fname:
    :param pickle_fname:
    :param oexp:    1.6497056483276888 is the globally optimized value for HCNH-HCNH in the set of 13 globular
                    proteins for global intensity scaling and peak-peak 2D-hist connectivities.
    :param selected_AAIG_connectivities_dict:   dict with selected connectivities of every NOESY AAIG1; has keys the
                                                AAIG signatures and values lists of quintuplets (tuples) consisting of
                                                (NOESY AAIG2 signature, # matching peaks, total # peaks, intersection,
                                                multi score). If this dict is provided then this method will write
                                                to files only the connectivites in this dict.

    :return:
    """
    if not pickle_fname:
        pickle_fname = os.path.splitext(fname)[0] + ".pkl"
    AAIG_connectivities_complete_dict = {}  # AAIG1 signature -> [(connected AAIG2 signature, number of matched AAIG1 peaks,
    # total number of AAIG1 peaks, AAIG1-AAIG2 2D density histogram intersection,
    # 'multi' score [combines occupancy and intersection]), ...]

    with open(fname, 'w') as f:
        f.write("# AAIG1 signature --> (connected AAIG2 signature, number of matched AAIG1 peaks, total number of "
                "AAIG1 peaks, AAIG1-AAIG2 2D density histogram intersection, 'multi' score [combines occupancy and "
                "intersection]), ...\n")
        NOESY_NOESY_peak_connectivities_mdict = deepcopy(self.NOESY_NOESY_peak_connectivities_mdict)
        for AAIG1_signature in list(self.NOESY_NOESY_peak_connectivities_mdict.keys()):
            if selected_AAIG_connectivities_dict and not AAIG1_signature in selected_AAIG_connectivities_dict.keys():
                continue

            # NOTE: the following commented code is not necessary, as these information are already stored in
            # NOTE: self.AAIG_connectivities_complete_dict.

            # # NOTE: Before the introduction of 'multi', conectivities were sorted byt occupancy and secondarily by resid.
            # #       After the introduction of 'multi', they are sorted by 'multi' only, because it is unique for each connection.
            # #       Moreover, replacement of AAIG names with signatures (e.g. R34N-H) renders sorting by resid more complicated,
            # #       therefore I decided to completely remove it.
            # connectivities_list = [(AAIG2_signature,
            #                         len(self.NOESY_NOESY_peak_connectivities_mdict[AAIG1_signature][AAIG2_signature]),
            #                         len(self.HCNH_spec.__get_AAIG__(AAIG1_signature).__get_all_peaks__()),
            #                         self.NOESY_NOESY_2Dhist_intersections_mdict[AAIG1_signature][AAIG2_signature],
            #                         self.NOESY_NOESY_2Dhist_intersections_mdict[AAIG1_signature][AAIG2_signature] * \
            #                         (float(len(self.NOESY_NOESY_peak_connectivities_mdict[AAIG1_signature][AAIG2_signature])) / \
            #                          len(self.HCNH_spec.__get_AAIG__(AAIG1_signature).__get_all_peaks__())) ** oexp)
            #                        for AAIG2_signature in list(self.NOESY_NOESY_peak_connectivities_mdict[AAIG1_signature].keys())]

            connectivities_list = self.AAIG_connectivities_complete_dict[AAIG1_signature]

            connectivities_list.sort(key=itemgetter(4), reverse=True)  # sort by 'multi'
            connectivities_list = [(c[0], c[1], c[2], c[3], c[4]) for c in connectivities_list]
            if selected_AAIG_connectivities_dict:  # FILTER: keep only the selected connectivities
                connAAIG2_list = [conn[0] for conn in selected_AAIG_connectivities_dict[AAIG1_signature]]
                connectivities_list = [conn for conn in connectivities_list if conn[0] in connAAIG2_list]
                for AAIG2_signature in self.NOESY_NOESY_peak_connectivities_mdict[AAIG1_signature].keys():
                    if not AAIG2_signature in connAAIG2_list:
                        del NOESY_NOESY_peak_connectivities_mdict[AAIG1_signature][AAIG2_signature]
            if len(connectivities_list) > 0:
                f.write(AAIG1_signature + "\t" + ', '.join(str(c) for c in connectivities_list) + "\n")
            AAIG_connectivities_complete_dict[AAIG1_signature] = connectivities_list
    # Save the same information in a pickle file, along withthe dictionarry with the all matching peaks between pairs
    # of AAIGs.
    save_pickle(pickle_fname, AAIG_connectivities_complete_dict,
                NOESY_NOESY_peak_connectivities_mdict)  # save the filtered versions of both dicts


@staticmethod
def write_NOESY_NOESY_connectivities_from_dict(AAIG_connectivities_dict, fname):
    """
    Method to write the contents of a connectivity dict to a file.

    :param AAIG_connectivities_dict:    connectivity dict to be written to the file. Must have the form:
                                    AAIG1 signature -> [(connected AAIG2 signature, number of matched AAIG1 peaks,
                                    total number of AAIG1 peaks, AAIG1-AAIG2 2D density histogram intersection,
                                    'multi' score [combines occupancy and intersection]), ...]
    :param fname:
    :return:
    """

    with open(fname, 'w') as f:
        f.write("# AAIG1 signature --> (connected AAIG2 signature, number of matched AAIG1 peaks, total number of "
                "AAIG1 peaks, AAIG1-AAIG2 2D density histogram intersection, 'multi' score [combines occupancy and "
                "intersection]), ...\n")    # header
        for k, v in list(AAIG_connectivities_dict.items()):
            # print k,v
            f.write(k + "\t" + ', '.join(map(str, v)) + "\n")


