from lib.hnnh_spectrum import *


class ProofreadHNNH:

    def __init__(self, TEST_HNNH, REF_HNNH, OUT_NAME, WORKDIR="."):

        test_HNNH_spec = HNNH_Spectrum(TEST_HNNH, uniquify_lines=False).spectrum_df
        save_pickle("%s/test_HNNH_spectrum.pkl" % WORKDIR, test_HNNH_spec)
        ref_HNNH_spec = HNNH_Spectrum(REF_HNNH, uniquify_lines=False).spectrum_df
        save_pickle("%s/ref_HNNH_spectrum.pkl" % WORKDIR, ref_HNNH_spec)
        proofread_test_HNNH_spec = pd.DataFrame([],
                                                columns=ref_HNNH_spec.columns.tolist(),
                                                dtype='float64')
        proofread_test_HNNH_spec = proofread_test_HNNH_spec.astype({'i_AAIG': 'str', 'j_AAIG': 'str', 'test_assignment':'str',
                                                                    'i_agreement': 'bool', 'j_agreement': 'bool'})
        # Write proof-read user file
        true_labels, false_labels =set(), set()
        for ui, urow in test_HNNH_spec.iterrows():
            for ri, rrow in ref_HNNH_spec.iterrows():
                if self.__resonances_match(urow, rrow):
                    # THE FIRST if-else CLAUSE IS JUST FOR STATISTICS
                    # One ?-? is overseen if the other AAIG is correct, but two ?-? are considered False.
                    if reverse_NH_suffix(urow['i_AAIG']) == rrow['i_AAIG'] and urow['j_AAIG'] == rrow['j_AAIG']: # or \
                            # reverse_NH_suffix(urow['i_AAIG']) == "?-?" and urow['j_AAIG'] == rrow['j_AAIG'] or \
                            # reverse_NH_suffix(urow['i_AAIG']) == rrow['i_AAIG'] and urow['j_AAIG'] == "?-?":
                        if 'H-N' not in rrow['i_AAIG'] and 'N-H' not in rrow['j_AAIG']:
                            true_labels.add("%s-%s" % (rrow['i_AAIG'], rrow['j_AAIG']))
                    else:
                        if 'H-N' not in rrow['i_AAIG'] and 'N-H' not in rrow['j_AAIG']:
                            false_labels.add("%s-%s" % (rrow['i_AAIG'], rrow['j_AAIG']))

                    if reverse_NH_suffix(urow['i_AAIG']) == rrow['i_AAIG']:
                        i_agree = True
                    elif reverse_NH_suffix(urow['i_AAIG']) == "?-?":
                        i_agree = np.nan
                    else:
                        i_agree = False

                    if urow['j_AAIG'] == rrow['j_AAIG']:
                        j_agree = True
                    elif urow['j_AAIG'] == "?-?":
                        j_agree = np.nan
                    else:
                        j_agree = False

                    # a = pd.Series(["%s-%s" % (rrow['i_AAIG'], rrow['j_AAIG']), i_agree, j_agree],
                    #               index=['test_assignment', 'i_agreement', 'j_agreement'])
                    # urow = urow.append(a) # OBSOLETE because these fields already exist by default in hnnh_spec!
                    urow['test_assignment'] = "%s-%s" % (rrow['i_AAIG'], rrow['j_AAIG'])
                    urow['i_agreement'] = i_agree
                    urow['j_agreement'] = j_agree
                    proofread_test_HNNH_spec = proofread_test_HNNH_spec.append(urow, ignore_index=True)
                    break

        # because test_HNNN_spec may contain replicate lines (but why?)
        total_true = len(true_labels)
        total_false = len(false_labels)
        all_labels = true_labels.union(false_labels)
        wrong_labels = list(all_labels.difference(true_labels))
        with open("%s/%s" % (WORKDIR, OUT_NAME), 'w') as f:
            for index, row in proofread_test_HNNH_spec.iterrows():
                f.write("%s-%s\t\t%f\t%f\t%f\t%f\t%f\t%s\t%s\t%s\t%s\t%s\n" % (
                    row['i_AAIG'], row['j_AAIG'], row['i_HNreson'], row['i_Nreson'],
                    row['j_Nreson'], row['j_HNreson'], row['Intensity'],
                    row['i_HSQCdist'], row['j_HSQCdist'],
                    row['test_assignment'], row['i_agreement'], row['j_agreement'])
                        )

        save_pickle("%s/proofread_test_HNNH_spec.pkl" % WORKDIR, proofread_test_HNNH_spec)

        missing_ref_lines = []
        for ri, rrow in ref_HNNH_spec.iterrows():
            if 'H-N' in rrow['i_AAIG'] or 'N-H' in rrow['j_AAIG']:  # we want only side-chains
                continue
            found = False
            for ui, urow in proofread_test_HNNH_spec.iterrows():
                if self.__resonances_match(urow, rrow):
                    # ref2test_matches += 1
                    found = True
                    break
            if not found:
                missing_ref_lines.append("%s-%s\t\t%f\t%f\t%f\t%f\t%f\n" % (rrow['i_AAIG'], rrow['j_AAIG'], rrow['i_HNreson'], rrow['i_Nreson'],
                    rrow['j_Nreson'], rrow['j_HNreson'], rrow['Intensity']))

        # Print statistics
        test_labels = set()  # assigned side chain HNNH peak labels by the user
        for ri, rrow in ref_HNNH_spec.iterrows():
            if 'H-N' not in rrow['i_AAIG'] and 'N-H' not in rrow['j_AAIG']:
                test_labels.add("%s-%s" % (rrow['i_AAIG'], rrow['j_AAIG']))
        ColorPrint("Total number of assigned side chain AAIGs in the reference HNNH file: %i" % len(test_labels), "BOLDGREEN")
        ColorPrint("Number of True side chain AAIG assignments: %i" % total_true, "BOLDGREEN")
        ColorPrint("Number of False side chain AAIG assignments: %i" % total_false, "BOLDGREEN")
        ColorPrint("The following labels were never assigned correctly: %s" % wrong_labels, "BOLDGREEN")
        if len(missing_ref_lines) > 0:
            ColorPrint("The following side-chain lines in ref spectrum were not assigned at all in the user spectrum:", "BOLDGREEN")
            for l in missing_ref_lines:
                print(l)
        with open("%s/HNNH_proofreading_stats.txt" % WORKDIR, 'w') as f:
            f.write("Total number of assigned side chain AAIGs in the reference HNNH file: %i\n" % len(test_labels))
            f.write("Number of True side chain AAIG assignments: %i\n" % total_true)
            f.write("Number of False side chain AAIG assignments: %i\n" % total_false)
            f.write("The following labels were never assigned correctly: %s\n" % wrong_labels)
            if len(missing_ref_lines) > 0:
                f.write("The following side-chain lines in ref spectrum were not assigned at all in the user spectrum:\n")
                for l in missing_ref_lines:
                    f.write(l + "\n")

    def __resonances_match(self, urow, rrow):
        return urow['i_HNreson'] == rrow['i_HNreson'] and \
               urow['i_Nreson'] == rrow['i_Nreson'] and \
               urow['j_HNreson'] == rrow['j_HNreson'] and \
               urow['j_Nreson'] == rrow['j_Nreson']