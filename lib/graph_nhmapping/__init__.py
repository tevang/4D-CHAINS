from lib.alignment import Alignment
from lib.fasta import FASTA
from lib.global_vars import aa1to3_dict
from lib.global_func import *

class Graph():

    def __init__(self):
        pass

    @staticmethod
    def write_Graph_algorthm_input(fasta_file, hsqc_spec, AAIG_connectivities_complete_dict):
        """
        Method to write the input files for Graph Algorithm.

        sequence.csv
    - sequence of the protein, for each amino acid, the frame is assigned
    - CSV format, the first row contains amino acids, the second row probability of gap (missing frame). The default probability is set to 0.1, if any part of sequence is known to have no frames, it can be rised up to 1.0 (which means "this is a gap, don't put frame there")

    frames.csv
    - a list of all frames, the number of frames must be the same or of greater as the size of the sequence
    - CSV format with just one row (the frames)

    adjacency-2d-nums.csv
    - a list of common peaks for each combination of frames
    - CSV format, each row is one pair of frames, three columns: frame1, frame2, number of peaks

    adjacency-3d-nums.csv
    - a list of common peaks classified as peaks of an amino acid
    - CSV format, each row is one triplet of two frames and one amino acid, columns: frame1, frame2, acid, number of  peaks

    assignment.csv
    - assignment of frames to amino acids
    - CSV format, each row is combination of frame and amino acid, columns: frame, acid, score

        :param AAIG_connectivities_complete_dict:
        :return:
        """
        fasta = FASTA(fasta_file)

        with open("sequence.csv", 'w') as f:
            for aa in fasta.protseq_list:
                f.write(aa1to3_dict[aa]+",")
            f.write("\n")
            f.write( ",".join(['0.1']*len(fasta.protseq_list)) + "\n")

        with open("frames.csv", 'w') as f:
            AAIG_signatures = ",".join(hsqc_spec.__get_all_AAIG_signatures__()) # without the PRO
            for i in range(fasta.protseq_list.count("P")):  # append the PRO
                AAIG_signatures += ",PRO%i" % i
            f.write(AAIG_signatures)

        with open("adjacency-2d-nums.csv", 'w') as f:
            for AAIG1sign, connectivities in AAIG_connectivities_complete_dict.items():
                for AAIG2sign, common_peak_num, total_peak_num, intersection, muliscore in connectivities:
                    f.write("%s,%s,%s\n" % (AAIG1sign, AAIG2sign, common_peak_num))

        # adjacency-3d-nums.csv needs the aa-type classifier.
        # with open("adjacency-3d-nums.csv", 'w') as f:
        #     valid_aas = list(aa1to3_dict.values())
        #     valid_aas.remove("PRO")
        #     for AAIG1sign, connectivities in AAIG_connectivities_complete_dict.items():
        #         for AAIG2sign, common_peak_num, total_peak_num, intersection, muliscore in connectivities:
        #             for aa in valid_aas:
        #                 f.write("%s,%s,%s\n" % (AAIG1sign, AAIG2sign, common_peak_num))


    @staticmethod
    def replace_peaknum_in_adjacency_2d_nums(adjacency_2d_file, AAIG_connectivities_complete_dict, metric="multi"):
        """
        Method to replace the number of common peaks in an existing adjacency-2d-nums.csv file with the
        specified metric from the input AAIG_connectivities_complete_dict.

        :param AAIG_connectivities_complete_dict:
        :param metric: can be 'interscection' for 2D-hist intersection, or 'multi' for MULTISCORE
        :return:
        """
        if metric == 'multi':
            i = 4
        elif metric == 'intersection':
            i = 3

        # Reformat AAIG_connectivities_complete_dict
        AAIG1_AAIG2_connvalue_dict = tree()
        for AAIG1sign, connectivities_list in AAIG_connectivities_complete_dict.items():
            for q in connectivities_list:
                AAIG2sign = q[0]
                AAIG1_AAIG2_connvalue_dict[AAIG1sign][AAIG2sign] = q[i]

        fname, ftype = os.path.splitext(adjacency_2d_file)
        fout = open(fname + "_mod" + ftype, 'w')
        with open(adjacency_2d_file, 'r') as fin:
            for line in fin:
                AAIG1sign, AAIG2sign, common_peak_num = line.split(',')
                try:
                    fout.write("%s,%s,%f\n" % (AAIG1sign, AAIG2sign, AAIG1_AAIG2_connvalue_dict[AAIG1sign][AAIG2sign]))
                except TypeError:
                    # NOTE: adjacency_2d_file contains all possible pairwise combinations of AAIGs, even if they don't
                    # NOTE: have any connectivity.
                    if int(common_peak_num) > 0:
                        ColorPrint("WARNING: there is no connectivity between %s and %s. Setting it to 0." %
                                   (AAIG1sign, AAIG2sign), "WARNING")
                    fout.write("%s,%s,%f\n" % (AAIG1sign, AAIG2sign, 0.0))
                fout.flush()
        fout.close()
        print("Created file %s with replaced values." % (fname + "_mod" + ftype))


"""
EXAMPLE: replace the number of common peaks in an existing adjacency-2d-nums.csv file with multi score.

from lib.graph_nhmapping import *
from lib.global_func import *

adjacency_2d_file = "/home2/thomas/Documents/4d-chains_dev/Graph_algorithm/FD3A/adjacency-2d-nums.csv"
AAIG_connectivities_complete_dict, NOESY_NOESY_peak_connectivities_mdict = \
load_pickle("/home2/thomas/Documents/4D-CHAINS_regtests/FD3A/NOESY_connectivities/FD3A_HCNH-HCNH_connectivities.pkl")

Graph.replace_peaknum_in_adjacency_2d_nums(adjacency_2d_file, AAIG_connectivities_complete_dict, metric="multi")
"""