from lib.global_func import *


############################################# START OF PROOF-READING METHODS ####################################################

def proofread_connectivities(self):
    # Proof read the connectivities and print statistics
    iplus1, iminus1, total_NH, total_N, total_connectivities = 0, 0, 0, 0, 0
    AAIGsign_iplus1_dict, AAIGsign_iplus1_CAHA_dict, AAIGsign_iminus1_dict, \
    AAIGsign_iminus1_CAHA_dict, AAIGsign_conn_dict, AAIGsign_conn_CAHA_dict = {}, {}, {}, {}, {}, {}     # values are 'True' if any of i->i+1 or i->i-1 connectivities was found, and 'False' if none was found.
    results_dict = {}
    for AAIG1sign in self.NOESY_NOESY_peak_connectivities_mdict.keys():
        total_connectivities += len(list(self.NOESY_NOESY_peak_connectivities_mdict[AAIG1sign].keys()))
        NH1_name = get_NH_name(AAIG1sign)
        total_N += 1
        if NH1_name != "NH":
            continue
        total_NH += 1
        AAIGsign_iplus1_dict[AAIG1sign] = False
        AAIGsign_iplus1_CAHA_dict[AAIG1sign] = False
        AAIGsign_iminus1_dict[AAIG1sign] = False
        AAIGsign_iminus1_CAHA_dict[AAIG1sign] = False
        AAIGsign_conn_dict[AAIG1sign] = False
        AAIGsign_conn_CAHA_dict[AAIG1sign] = False
        for AAIG2sign in self.NOESY_NOESY_peak_connectivities_mdict[AAIG1sign].keys():
            NH2_name = get_NH_name(AAIG2sign)
            if NH2_name != "NH":
                continue
            resid1 = get_resid(AAIG1sign)
            resid2 = get_resid(AAIG2sign)
            if resid2 == resid1-1:
                iminus1 += 1
                AAIGsign_iminus1_dict[AAIG1sign] = True
                AAIGsign_conn_dict[AAIG1sign] = True
                # Find out if the i->(i-1) connectivity contains a CA-HA connection
                AAIGsign_iminus1_CAHA_dict[AAIG1sign] = \
                    any([pp[0].f_Cname=="CA" or pp[1].f_Cname=="CA" for pp in self.NOESY_NOESY_peak_connectivities_mdict[AAIG1sign][AAIG2sign]])
                AAIGsign_conn_CAHA_dict[AAIG1sign] = AAIGsign_conn_CAHA_dict[AAIG1sign] or \
                                                     AAIGsign_iminus1_CAHA_dict[AAIG1sign]
            elif resid2 == resid1+1:
                iplus1 += 1
                AAIGsign_iplus1_dict[AAIG1sign] = True
                AAIGsign_conn_dict[AAIG1sign] = True
                # Find out if the i->(i+1) connectivity contains a CA-HA connection
                AAIGsign_iplus1_CAHA_dict[AAIG1sign] = \
                    any([pp[0].f_Cname=="CA" or pp[1].f_Cname=="CA" for pp in
                         self.NOESY_NOESY_peak_connectivities_mdict[AAIG1sign][AAIG2sign]])
                AAIGsign_conn_CAHA_dict[AAIG1sign] = AAIGsign_conn_CAHA_dict[AAIG1sign] or \
                                                     AAIGsign_iplus1_CAHA_dict[AAIG1sign]

    results_dict["total_connectivities"] = total_connectivities # including the sidechains
    results_dict["mean_connectivities_per_AAIG"] = total_connectivities/total_N # including sidechains
    results_dict["total_iplus1"] = np.sum(list(AAIGsign_iplus1_dict.values()))
    results_dict["total_iplus1_CAHA"] = np.sum(list(AAIGsign_iplus1_CAHA_dict.values()))
    results_dict["total_iminus1"] = np.sum(list(AAIGsign_iminus1_dict.values()))
    results_dict["total_iminus1_CAHA"] = np.sum(list(AAIGsign_iminus1_CAHA_dict.values()))
    results_dict["total_NH"] = total_NH
    results_dict["total_iplus1_or_iminus1"] = np.sum(list(AAIGsign_conn_dict.values())) # number of AAIGs that had connectivity with
                                                                                        # at least one of these two.
    results_dict["total_iplus1_or_iminus1_CAHA"] = np.sum(list(AAIGsign_conn_CAHA_dict.values())) # number of AAIGs that had connectivity with
                                                                                        # at least one of these two.
    return results_dict