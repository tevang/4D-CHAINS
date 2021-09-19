import copy

from cluster import HierarchicalClustering
from lib.hcnh_spectrum import *
from lib.hnnh_spectrum import *
from lib.trees.connectivities_tree import *


class Connectivities(object):

    def __init__(self,
                 # ABSOLUTE_MATCHES_FILE,
                 # FIRST_RESIDUE_NUMBER,
                 tolH=0.04,
                 tolC = 0.4,
                 H_bin_length=0.0005,
                 C_bin_length=0.005,
                 bandwidth=0.2
                 ):
        """

        :param tolH: always use 0.04
        :param tolC: always use 0.4
        :param H_bin_length: for AAIG connectivities use 0.01, for peak connectivities use
        :param C_bin_length: for AAIG connectivities use 0.1, for peak connectivities use
        :param bandwidth: for AAIG connectivities use 0.0025, for peak connectivities use
        """
        # self.ABSOLUTE_MATCHES_FILE, self.FIRST_RESIDUE_NUMBER = ABSOLUTE_MATCHES_FILE, FIRST_RESIDUE_NUMBER
        self.tolH, self.tolC = tolH, tolC
        self.H_bin_length, self.C_bin_length, self.bandwidth = H_bin_length, C_bin_length, bandwidth
        self.hist_borders = {}

        self.AAIG_connectivities_complete_dict = {}
        # ** In case of TOCSY-HCNH:
        #   TOCSY AAIG(i) -> [(NOESY AAIG(i-1), occupancy, tot. num. peaks, intersection, multi score), ...]
        # ** In case of HCNH-HCNH:
        #   NOESY AAIG(i) -> [(NOESY AAIG(j), occupancy, tot. num. peaks, intersection, multi score), ...]
        # ** In case of HNNH:
        #   NOESY AAIG(i) -> [(NOESY AAIG(j), 1, tot. num. connectivities of AAIG(i), peak intensity, peak intensity), ...]

        self.TOCSY_spec = None
        self.HCNH_spec = None
        self.HNNH_spec = None

        # # Create all alignment files # ZATIM JE ZBYTECNY
        # self.ali = Alignment(self.ABSOLUTE_MATCHES_FILE, self.FIRST_RESIDUE_NUMBER)
        # self.ali.create_absolute_AAIGmatches_alignment()
        # self.ali.create_absolute_matches_alignment()

    def load_HCNH_spectrum(self, HCNH_FILE):
        """
        Method used by both TOCSY- and NOESY- connectivities.
        :param HCNH_FILE:
        :return:
        """
        self.HCNH_spec = HCNH_Spectrum(HCNH_FILE=HCNH_FILE,
                                       H_bin_length=self.H_bin_length,
                                       C_bin_length=self.C_bin_length,
                                       bandwidth=self.bandwidth)
        # if the TOCSY spectrum has been also loaded, check if it is identical to the NOESY spectrum
        if self.TOCSY_spec and self.HCNH_spec.__is_equal__(self.TOCSY_spec):
            raise Exception(
                bcolors.FAIL + "ERROR: NOESY and TOCSY spectra are identical! Check your input files!" + bcolors.ENDC)

    def load_NMR_PIPELINE_HCNH_spectrum(self, HCNH_spectrum_csv):
        """
        Method used by both TOCSY- and NOESY- connectivities.
        :param HCNH_FILE:
        :return:
        """
        self.HCNH_spec = HCNH_Spectrum(HCNH_spectrum_csv=HCNH_spectrum_csv,
                                       H_bin_length=self.H_bin_length,
                                       C_bin_length=self.C_bin_length,
                                       bandwidth=self.bandwidth)
        # if the TOCSY spectrum has been also loaded, check if it is identical to the NOESY spectrum
        if self.TOCSY_spec and self.HCNH_spec.__is_equal__(self.TOCSY_spec):
            raise Exception(
                bcolors.FAIL + "ERROR: NOESY and TOCSY spectra are identical! Check your input files!" + bcolors.ENDC)

    def load_HNNH_spectrum(self, HNNH_FILE):
        """
        Method used by both TOCSY- and NOESY- connectivities.
        :param HNNH_FILE: can be a CSV or a Sparky list file
        :return:
        """
        self.HNNH_spec = HNNH_Spectrum(HNNH_FILE=HNNH_FILE, H_bin_length=self.H_bin_length,
                                       C_bin_length=self.C_bin_length, bandwidth=self.bandwidth)
        # Intensities will be scaled when you will load the HNNH spectrum to a DataFrame using HNNH_Spectrum.load_HNNH_to_DataFrame()


    def find_2Dhist_boundaries(self):
        """
         Basically split the 2D-hist representation into aliphatic and into aromatic region to save memory
         and time and find their boundaries from the spectrum chemical shifts. Method used by both TOCSY-
         and NOESY-connectivities.
        :return:
        """
        aliphatic_carbons = []
        aliphatic_protons = []
        aromatic_carbons = []
        aromatic_protons = []
        aromatic_protons = []
        for peak in self.HCNH_spec.__get_all_peaks__():
            Creson, Hreson = peak.Creson, peak.Hreson
            if Creson < 100:    # aliphatic carbon
                aliphatic_carbons.append(Creson)
                aliphatic_protons.append(Hreson)
            elif Creson >= 100: # aromatic carbon
                aromatic_carbons.append(Creson)
                aromatic_protons.append(Hreson)
        if self.TOCSY_spec:
            for peak in self.TOCSY_spec.__get_all_peaks__():
                Creson, Hreson = peak.Creson, peak.Hreson
                if Creson < 100:    # aliphatic carbon
                    aliphatic_carbons.append(Creson)
                    aliphatic_protons.append(Hreson)
                elif Creson >= 100: # aromatic carbon
                    aromatic_carbons.append(Creson)
                    aromatic_protons.append(Hreson)
        aaig = AAIG(H_bin_length=self.H_bin_length, C_bin_length=self.C_bin_length)
        al_H_minmax, al_C_minmax = aaig.even_edges(
            H_minmax=[min(aliphatic_protons) - 0.2, max(aliphatic_protons) + 0.2],
            C_minmax=[min(aliphatic_carbons) - 2.0, max(aliphatic_carbons) + 2.0])
        if len(aromatic_protons) > 0:
            ar_H_minmax, ar_C_minmax = aaig.even_edges(
                H_minmax=[min(aromatic_protons) - 0.2, max(aromatic_protons) + 0.2],
                C_minmax=[min(aromatic_carbons) - 2.0, max(aromatic_carbons) + 2.0])
        else:   # if not aromatic peaks are in the spectrum, use the default values
            ar_H_minmax, ar_C_minmax = [Dec(6.0), Dec(7.7)], [Dec(113.0), Dec(136.0)]

        self.hist_borders = {   'alHmin': al_H_minmax[0],
                                'alHmax': al_H_minmax[1],
                                'alCmin': al_C_minmax[0],
                                'alCmax': al_C_minmax[1],
                                'arHmin': ar_H_minmax[0],
                                'arHmax': ar_H_minmax[1],
                                'arCmin': ar_C_minmax[0],
                                'arCmax': ar_C_minmax[1]
                                }

    def NOESY_AAIG2hist(self, AAIG_name):
        """
        Create the 2D-histogram of this NOESY AAIG.
        Method used by both TOCSY- and NOESY- connectivities.
        :param AAIG_name:
        :return:
        """
        self.HCNH_spec.__AAIG2hist__(AAIG_name, even_intensities=False, scale_density=False, **self.hist_borders)
        return self.HCNH_spec.__get_AAIG__(AAIG_name)

    def NOESY_Peaks2hist(self, AAIG_name, even_edges, scale_density):
        """
        Convenience method for parallel execution with SCOOP, creates the 2D-histogram of all individual peaks of
        this NOESY AAIG. Method used by both TOCSY- and NOESY- connectivities.
        :param AAIG_name:
        :return:
        """
        self.HCNH_spec.__AAIG_Peaks2hist__(AAIG_name,
                                           even_edges=even_edges,
                                           even_intensities=False,
                                           scale_density=scale_density)
        return self.HCNH_spec.__get_AAIG__(AAIG_name)

    def intersect_connectivities(self, ref_conn):
        """
        Method to filter self.AAIG_connectivities_complete_dict by keeping only those connectivities that are common with
        the input ref_AAIG_connectivities_complete_dict, basically finding the intersection between these two.
        This method is used to deconvolute the HCNH NOESY connectivities using the HCNH TOCSY connectivities or the
        HNNH connectivities.
        :param ref_AAIG_connectivities_complete_dict:
        :return:
        """
        AAIG_connectivities_complete_dict = copy.deepcopy(self.AAIG_connectivities_complete_dict)
        for AAIG1_signature in list(AAIG_connectivities_complete_dict.keys()):
            connectivities_list = self.AAIG_connectivities_complete_dict[AAIG1_signature]
            try:
                ref_connectivities_list = ref_conn.AAIG_connectivities_complete_dict[AAIG1_signature]
            except KeyError:
                continue
            AAIG2_signature_set = {q[0] for q in connectivities_list}
            ref_AAIG2_signature_set = {q[0] for q in ref_connectivities_list}
            common_AAIG2_signatures_set = AAIG2_signature_set.intersection(ref_AAIG2_signature_set)
            self.AAIG_connectivities_complete_dict[AAIG1_signature] = [q for q in connectivities_list if q[0] in common_AAIG2_signatures_set]    # keep only the common connectivities

    @staticmethod
    def load_connectivities_from_file(CON_FILE):
        """
        Works with old (3 elements) and new (5 elements) connectivity formats.
        Method used by both TOCSY- and NOESY- connectivities.
        :param CON_FILE:
        :return:
        """
        AAIG_connectivities_dict = defaultdict(list)
        with open(CON_FILE, 'r') as f:
            for line in f:
                m = re.search("^([A-Za-z0-9]+N[DE12X]*H[DE12X]*)\s+(\(.*)$", line)
                if m:
                    i_AAIG_name = m.group(1)
                    s = "connectivities = " + m.group(2) + ","  ; # comma to support the single connectivities, too.
                    namespace = {}  # in Python 3 you must define a new namespace for exec()
                    exec(s, namespace)
                    connectivities = list(namespace['connectivities'])
                    AAIG_connectivities_dict[i_AAIG_name] = connectivities
        return AAIG_connectivities_dict     # the i_iminus1_dict

    @staticmethod
    def chunk_connections_dict(connections_dict, tolC=0.2):
        """
        Method to chunk the connections_dict into smaller dicts that contain as keys peaks with Carbons close to each other.
        This can reduce the memory overhead during search for the maximum possible number of peak pairs between two AAIGs.

        :param connections_dict:
        :param tolC:
        :return:
        """
        # print("\nDEBUG: initial connections_dict=")
        # for k,v in connections_dict.items():
        #     print(k.__get_CHresons__(), "-->", [p[0].__get_CHresons__() for p in v])

        # connections_dict[peak1] = [(peak2A, distance), (peak2B, distance), ...], where peak1 belongs to AAIG1 and peak2[AB] to AAIG2
        # ATTENTION: do not use Creson dictionary because the AAIG may have multiple peaks with the same carbon resonance!
        Creson_list = [AAIG1peak.Creson for AAIG1peak in connections_dict.keys()]
        # print("DEBUG: Creson_list=", Creson_list)
        cl = HierarchicalClustering(Creson_list, lambda x, y: float(abs(x - y)))
        # ATTENTIONS: YOU MIGHT NEED TO REDUCE IT TO SOMETHING LIKE 4*tolC IN CASE OF MEMORY OVERLOAD.
        cutoff = 4*tolC  # for security, although 2*tolC may have been enough
        cluster_list = cl.getlevel(cutoff)
        # print("DEBUG: cluster_list=", cluster_list)
        connections_dict_list = []
        for cluster in cluster_list:
            if type(cluster) == float:  # if for some reason there is only one Creson value in the connections_dict
                                        # but not necessarily only one peak.
                cluster = [cluster]
            # print("DEBUG: connections_dict=", connections_dict)
            # print("DEBUG: cluster=", cluster)
            connections_chunk_dict = {k:v for k,v in connections_dict.items() if k.Creson in cluster}
            connections_dict_list.append(connections_chunk_dict)

        # VERIFICATION
        assert len([k for d in connections_dict_list for k in d.keys()]) == len(connections_dict.keys()), \
            Debuginfo("ERROR: the chunks of connections_dict do not have all the keys of the initial connections_dict!", fail=True)

        assert len([v for d in connections_dict_list for vl in d.values() for v in vl]) == \
               len([v for vl in connections_dict.values() for v in vl]), \
            Debuginfo("ERROR: the chunks of connections_dict do not have all the values of the initial connections_dict!", fail=True)

        # print("\nDEBUG: list of connections_dict chunks=")
        # for d in connections_dict_list:
        #     for k,v in d.items():
        #         print(k.__get_CHresons__(), "-->", [p[0].__get_CHresons__() for p in v])

        return connections_dict_list

