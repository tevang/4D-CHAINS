######################################### START OF GETTER METHODS ########################################
from operator import itemgetter

from lib.global_func import *
from copy import deepcopy
from sklearn.preprocessing import MinMaxScaler
import pandas as pd

def __get_NOESY_NOESY_peak_connectivities_mdict__(self):
    """
    Getter method that returns the connectivities in 4D-CHAINS object's form.

    :return NOESY_NOESY_peak_connectivities_mdict: 3D dictionary of the form:
    AAIG(1) signature -> matching AAIG(2) signature -> [(matching NOESY peak1, matching NOESY peak2, distance), ...]
    """
    return self.NOESY_NOESY_peak_connectivities_mdict

def __get_NOESY_NOESY_connectivities_dict__(self, oexp=1.6497056483276888):
    """
    OLD oexp=3.0448217240563293, OPTIMIZED FOR TOCSY-NOESY.

    Getter method that returns all possible connectivities of every NOESY AAIG with other NOESY AAIGs in a dictionary form.
    It presumes that match_NOESY_to_NOESY_Peaks() has been already called, otherwise it will return an empty
    dictionary.

    :param oexp:    1.6497056483276888 is the globally optimized value in the set of 13 globular proteins for global
                    intensity scaling and peak-peak 2D-hist connectivities.
    :return: AAIG_connectivities_complete_dict: dictionary of the form: NOESY AAIG1 signature --> [(NOESY AAIG2 signature,
                                        occupancy, numOfResonances, intersection, multi score), ...]
    """

    for AAIG1_signature in list(self.NOESY_NOESY_peak_connectivities_mdict.keys()):
        # NOTE: Before the introduction of 'multi', conectivities were sorted by occupancy and secondarily by resid.
        #       After the introduction of 'multi', they are sorted by 'multi' __ONLY__, because it is unique for each connection.
        #       Moreover, replacement of AAIG names with signatures (e.g. R34N-H) renders sorting by resid more complicated,
        #       therefore I decided to completely remove it.
        try:
            connectivities_list = [(AAIG2_signature,
                                    len(self.NOESY_NOESY_peak_connectivities_mdict[AAIG1_signature][AAIG2_signature]),
                                    len(self.HCNH_spec.__get_AAIG__(AAIG1_signature).__get_all_peaks__()),
                                    self.NOESY_NOESY_2Dhist_intersections_mdict[AAIG1_signature][AAIG2_signature],
                                    self.NOESY_NOESY_2Dhist_intersections_mdict[AAIG1_signature][AAIG2_signature] * \
                                    (float(len(self.NOESY_NOESY_peak_connectivities_mdict[AAIG1_signature][AAIG2_signature])) / \
                                     len(self.HCNH_spec.__get_AAIG__(AAIG1_signature).__get_all_peaks__())) ** oexp)
                                   for AAIG2_signature in list(self.NOESY_NOESY_peak_connectivities_mdict[AAIG1_signature].keys())]
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
        connectivities_list = [[c[0], c[1], c[2], c[3], c[4]] for c in connectivities_list]
        self.AAIG_connectivities_complete_dict[AAIG1_signature] = connectivities_list
    self._scale_multiscores(scale_intersections=True)
    return self.AAIG_connectivities_complete_dict

def __get_NOESY_NOESY_connectivities_dataframe__(self, scale=True):
    """
    OLD oexp=3.0448217240563293, OPTIMIZED FOR TOCSY-NOESY.

    Getter method that returns all possible connectivities of every NOESY AAIG with other NOESY AAIGs in a dictionary form.
    It presumes that match_NOESY_to_NOESY_Peaks() has been already called, otherwise it will return an empty
    dictionary.

    :param oexp:    1.6497056483276888 is the globally optimized value in the set of 13 globular proteins for global
                    intensity scaling and peak-peak 2D-hist connectivities.
    :return: AAIG_connectivities_complete_dict: dictionary of the form: NOESY AAIG1 signature --> [(NOESY AAIG2 signature,
                                        occupancy, numOfResonances, intersection, multi score), ...]
    """
    conn_data = []  # [ ["frame1", "frame2", "number of common peaks", "occupancy", "2D-hist intersection"], ...]
    for AAIG1_signature in list(self.NOESY_NOESY_peak_connectivities_mdict.keys()):
        # NOTE: Before the introduction of 'multi', conectivities were sorted by occupancy and secondarily by resid.
        #       After the introduction of 'multi', they are sorted by 'multi' __ONLY__, because it is unique for each connection.
        #       Moreover, replacement of AAIG names with signatures (e.g. R34N-H) renders sorting by resid more complicated,
        #       therefore I decided to completely remove it.
        conn_list = []
        total_peaknum = len(self.HCNH_spec.__get_AAIG__(AAIG1_signature).__get_all_peaks__())   # of AAIG1
        try:
            for AAIG2_signature in self.NOESY_NOESY_peak_connectivities_mdict[AAIG1_signature].keys():
                conn_list.append( [AAIG1_signature,
                                   AAIG2_signature,
                                   len(self.NOESY_NOESY_peak_connectivities_mdict[AAIG1_signature][AAIG2_signature]),
                                   len(self.NOESY_NOESY_peak_connectivities_mdict[AAIG1_signature][AAIG2_signature]),
                                   len(self.NOESY_NOESY_peak_connectivities_mdict[AAIG1_signature][AAIG2_signature])/total_peaknum,
                                   self.NOESY_NOESY_2Dhist_intersections_mdict[AAIG1_signature][AAIG2_signature]]
                                  )
        except IndexError:
            # if AAIG2_signature == "" or len(connectivities_list) == 0:
            #     pass
            # else:
            print("DEBUG: AAIG1_signature=", AAIG1_signature)
            print("DEBUG: AAIG2_signature=", AAIG2_signature)
            print("DEBUG: type(AAIG2_signature)=", type(AAIG2_signature))
            print("DEBUG: conn_list=", conn_list)
            raise IndexError()
        max_common_peaknum = np.max([c[2] for c in conn_list])  # the maximum common peak number for AAIG1
        for c in conn_list: c[3] /= max_common_peaknum     # relative occupancy
        conn_list.sort(key=itemgetter(5), reverse=True)  # sort by 'intersection'
        conn_data.extend(conn_list)
    conn_df = pd.DataFrame(data=conn_data,
                           columns=["source_frame", "target_frame", "num_common_peaks",
                                    "relative_occupancy", "occupancy", "intersection"])

    if scale == True:
        mscaler = MinMaxScaler(feature_range=(conn_df.relative_occupancy.min(), 100))
        conn_df["scaled_relative_occupancy"] = [0] * conn_df["relative_occupancy"].size
        conn_df["scaled_occupancy"] = [0] * conn_df["occupancy"].size
        conn_df["scaled_intersection"] = [0] * conn_df["occupancy"].size
        conn_df[["scaled_relative_occupancy", "scaled_occupancy", "scaled_intersection"]] = \
            mscaler.fit_transform(conn_df[["relative_occupancy", "occupancy", "intersection"]])

    conn_df = conn_df.round(4)
    df_data = [pp for g in
               conn_df.apply(lambda r: [[r["source_frame"], r["target_frame"],
                                         p1.Hreson, p1.Creson, p1.j_HNreson, p1.j_Nreson,
                                         p1.Intensity, p2.Hreson, p2.Creson, p2.j_HNreson,
                                         p2.j_Nreson,  p2.Intensity] for p1, p2 in
                                        self.NOESY_NOESY_peak_connectivities_mdict[
                                            r["source_frame"]][
                                            r["target_frame"]]],
                             axis=1 )
               for pp in g]
    common_peak_df = pd.DataFrame(data=df_data, columns=["source_frame", "target_frame", "source_HC",
                                                         "source_C", "source_HN", "source_N",
                                                         "source_intensity", "target_HC",
                                                         "target_C", "target_HN", "target_N",
                                                         "target_intensity"])
    # Add dashes to frame signatures before returning the DataFrames
    conn_df["source_frame"] = conn_df["source_frame"].apply(dash_to_NH)
    conn_df["target_frame"] = conn_df["target_frame"].apply(dash_to_NH)
    common_peak_df["source_frame"] = common_peak_df["source_frame"].apply(dash_to_NH)
    common_peak_df["target_frame"] = common_peak_df["target_frame"].apply(dash_to_NH)

    return conn_df, common_peak_df.astype({"source_intensity": int, "target_intensity": int})

def __get_NOESY_NOESY_common_peak_dataframe__(self):

    # AAIG(1) signature -> matching AAIG(2) signature ->
    # [(matching NOESY peak1, matching NOESY peak 2), ...]
    common_peak_data = []
    columns = ["source_frame", "source_C", "source_H", "source_N", "source_HN",
               "target_frame", "target_C", "target_H", "target_N", "target_HN"]
    for AAIG1_sign in self.NOESY_NOESY_peak_connectivities_mdict.keys():
        for AAIG2_sign in self.NOESY_NOESY_peak_connectivities_mdict[AAIG1_sign].keys():
            for p1, p2 in self.NOESY_NOESY_peak_connectivities_mdict[AAIG1_sign][AAIG2_sign]:
                common_peak_data.append( [p1.j_AAIG_signature, p1.Creson, p1.Hreson, p1.j_Nreson, p1.j_HNreson,
                                          p2.j_AAIG_signature, p2.Creson, p2.Hreson, p2.j_Nreson, p2.j_HNreson]
                                         )
    common_peak_df = pd.DataFrame(data=common_peak_data, columns=columns)
    return common_peak_df

def _scale_multiscores(self, scale_intersections=False):
    """
    Internal method that transforms all MULTISCORE values to the range [min(occupancy),100].
    For the Graph algorithm, it makes sense to scale the values globally, namely wrt the whole spectrum.

    :param self:
    :return: updates self.AAIG_connectivities_complete_dict
    """

    all_occup, all_multi = [], []
    for connectivities_list in self.AAIG_connectivities_complete_dict.values():
        for c in connectivities_list:
            all_occup.append(c[1])
            all_multi.append(c[4])
    all_multi = np.array(all_multi).reshape(-1, 1)
    mscaler = MinMaxScaler(feature_range=(min(all_occup), 100))
    mscaler.fit(all_multi)
    if scale_intersections:
        all_intersections = [c[3] for connectivities_list in self.AAIG_connectivities_complete_dict.values()
                     for c in connectivities_list]
        all_intersections = np.array(all_intersections).reshape(-1, 1)
        iscaler = MinMaxScaler(feature_range=(min(all_occup), 100))

        iscaler.fit(all_intersections)

    for AAIGsign, connectivities_list in self.AAIG_connectivities_complete_dict.items():
        for c in connectivities_list:
            c[4] = mscaler.transform([[c[4]]])[0][0]   # AAIG_connectivities_complete_dict is updated automatically
            if scale_intersections:
                c[3] = iscaler.transform([[c[3]]])[0][0]  # AAIG_connectivities_complete_dict is updated automatically


def __get_filtered_NOESY_NOESY_connectivities_dict__(self,
                                                     MULTI_RATIO,
                                                     reverse_conn=False,
                                                     proofread=False,
                                                     pool_AAIG_connectivities_dict=None):
    """
    Filter the self.AAIG_connectivities_complete_dict, by applying the MULTI_RATIO criterion. Also, for any two AAIGs,
    AAIG1 and AAIG2, look if both connectivities AAIG1-AAIG2 and AAIG2-AAIG1 exist in the produced chains under the
    given MULTI_RATIO. If not, then remove the one that exists (either AAIG1-AAIG2, or AAIG2-AAIG1).

    :param MULTI_RATIO: keep all connected AAIG2s that have 'multi' score up to MULTI_RATIO times lower than the maximum
                        'multi' score of the connectivities of AAIG1. The higher the MULTI_RATIO, the more connectivities
                        are retained.
    :param reverse_conn:    if not both AAIG1-AAIG2 and AAIG2-AAIG1 exist, then remove any of them from the i_iminus1_dict
    :param proofread:   print ERROR whenever a wrong connectivity is removed by the REVERSE CONNECTIVITY criterion.
    :param pool_AAIG_connectivities_dict: if this connectivities dict is provided, then it will be used instead
                                          self.AAIG_connectivities_complete_dict and be filtered by the MULTI_RATIO &
                                          reverse_conn criteria.
    :return i_iminus1_dict: dict with selected connectivities of every NOESY AAIG1; has keys the
            AAIG signatures and values lists of quintuplets (tuples) consisting of (NOESY AAIG2 signature, # matching peaks,
            total # peaks, intersection, multi score).
    """
    ColorPrint("Filtering the HCNH-HCNH connectivities by applying the MULTI_RATIO criterion.", "BOLDGREEN")
    if pool_AAIG_connectivities_dict == None:   # if no refined connectivity pool is provided, use all possible connectivities
        pool_AAIG_connectivities_dict = self.AAIG_connectivities_complete_dict

    num_before = len(list(flatten(pool_AAIG_connectivities_dict.values())))
    # Iterate over all NOESY aa indices and, if applicable, keep those that exceed the multi threshold
    i_iminus1_dict = {}
    AAIG_AAIG_connBool_mdict = tree()  # AAIG1->AAIG2->True/False depending on the existence of this connectivity
    for i_AAIG_signature, connectivities in list(pool_AAIG_connectivities_dict.items()):
        filtered_connectivities = []
        max_multi = np.max([float(q[4]) for q in connectivities])  # max multi for this AAIG i
        for quintuplet in connectivities:
            j_AAIG_signature = quintuplet[0]
            AAIG_AAIG_connBool_mdict[i_AAIG_signature][j_AAIG_signature] = False  # initialize this dict. It will be updated later
            AAIG_AAIG_connBool_mdict[j_AAIG_signature][i_AAIG_signature] = False  # initialize this dict. It will be updated later
            multi = float(quintuplet[4])
            if max_multi / multi <= MULTI_RATIO and multi > 10e-6:  # <== CHANGE ME
                filtered_connectivities.append(quintuplet)
            elif proofread and is_valid_signature(i_AAIG_signature) and is_valid_signature(j_AAIG_signature) \
                    and abs(get_resid(i_AAIG_signature) - get_resid(j_AAIG_signature)) == 1 and not \
                    is_sidechain(i_AAIG_signature) and not is_sidechain(j_AAIG_signature):
                ColorPrint("ERROR: correct connectivity %s-%s was removed by MULTI_RATIO criterion!" %
                           (i_AAIG_signature, j_AAIG_signature),
                           "OKRED")
        i_iminus1_dict[i_AAIG_signature] = filtered_connectivities
    num_after = len(list(flatten(i_iminus1_dict.values())))
    ColorPrint("The number of connectivities was reduced from %i to %i using MULTI_RATIO=%.1f, "
               "namely by %.3f %%." %
               (num_before, num_after, MULTI_RATIO, 100 * (num_before - num_after) / num_before),
               "OKGREEN")

    if reverse_conn:
        ColorPrint("Filtering the HCNH-HCNH connectivities by applying the REVERSE CONNECTIVITY criterion.", "BOLDGREEN")
        # Make sure that  both connectivities AAIG1-AAIG2 and AAIG2-AAIG1 passed the MULTI_RATIO criterion
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

