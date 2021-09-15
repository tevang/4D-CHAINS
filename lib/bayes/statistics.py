import pickle
import re
from collections import OrderedDict, defaultdict
from operator import itemgetter
import numpy as np
from lib.fasta import FASTA
from lib.global_func import tree, ColorPrint
from scipy import stats
from scipy.stats import zscore


class Probability():

    def __init__(self, fasta_file, classifier_file=None):

        fasta = FASTA(fasta_file=fasta_file, full_seq=True)
        self.aatype_P_dict = fasta.count_AA_frequencies()
        self.classifier_file = classifier_file

        # The goal of this Class is to populate the following dicts.
        self.AAIG_aaTypesProbTupleList_dict = OrderedDict()   # dict with ALL AAIGs, their predicted aa-types and associated probabilities

        self.AAIG_aaType_Condprob_mdict = tree()  # multidict of the form AAIG of residue i -->
                                                      # residue i matching aa type --> conditional probability
        self.AAIG_aaTypesTransformedProbTupleList_dict = OrderedDict()  # AAIG --> [ (aatype, transformed probability), ... ]
        self.AAIG_aaTypesProbPoolTupleList_dict = OrderedDict()         # same as the above, but with only the possible aa-types
                                                                        # that were not filtered out from previous cycles.

        # Make TOCSY-relevant copies for compatibility with the old code. They will be automatically updated.
        self.iAAIG_iminus1aaType_Condprob_mdict = self.AAIG_aaType_Condprob_mdict  # multidict of the form AAIG of residue i -->
                                                              # residue i-1 matching aa type --> conditional probability
        self.iAAIG_iminus1aaTypesTransformedProbTupleList_dict = self.AAIG_aaTypesTransformedProbTupleList_dict  # AAIG --> [ (aatype, transformed probability), ...]
        self.iAAIG_iminus1aaTypesProbPoolTupleList_dict = self.AAIG_aaTypesProbPoolTupleList_dict     # same as the above, but with only the possible aa-types
                                                                            # that were not filtered out from previous cycles.

    @staticmethod
    def transform(probability_array, transform_type):
        """
            FUNCTION to do a mathematical transform on a value
        """
        if transform_type == "None":
            return probability_array
        elif transform_type == "log":
            return np.log(probability_array)
        elif transform_type == "log10":
            return np.log10(probability_array)
        elif transform_type == "boxcox_pearson":
            lmax_pearsonr = stats.boxcox_normmax(probability_array)
            prob_pearson = stats.boxcox(probability_array, lmbda=lmax_pearsonr)
            return prob_pearson
        elif transform_type == "boxcox_mle":
            prob_mle, lmax_mle = stats.boxcox(probability_array)
            return prob_mle

    def calc_TOCSY_aa_type_conditional_probs(self, iAAIG_iminus1aaTypesProbTupleList_dict, transform_type):
        """
        This method was designed to work with the old TOCSY assignment code, not with TOCSY Spectrum objects.

        :param iAAIG_iminus1aaTypesProbTupleList_dict: raw probabilities either calculated or loaded from a file
        :param transform_type:
        :return:
        """
        ## CALCULATE THE PROBABILITY P[AATYPE(I-1)]
        ## DO A MATHEMATICAL TRANSFORM ON P
        # iAAIG_iminus1aaTypesProbTupleList_dict is an ordereddict with keys the AAIG of residue i and values lists of tuples of the form
        # (residue i-1 matching aa type, average probability)
        all_prob_list = []
        for iAAIG, iminus1aaTypesProbTuple_list in list(iAAIG_iminus1aaTypesProbTupleList_dict.items()):
            for duplet in iminus1aaTypesProbTuple_list:
                prob = duplet[1]
                all_prob_list.append(prob)

        all_prob_array = np.array(all_prob_list)
        transformed_all_prob_array = Probability.transform(all_prob_array, transform_type=transform_type)

        array_index = 0
        for iAAIG, iminus1aaTypesProbTuple_list in list(iAAIG_iminus1aaTypesProbTupleList_dict.items()):
            self.iAAIG_iminus1aaTypesTransformedProbTupleList_dict[iAAIG] = []
            for duplet in iminus1aaTypesProbTuple_list:
                aatype = duplet[0]
                trans_prob = transformed_all_prob_array[array_index]
                self.iAAIG_iminus1aaTypesTransformedProbTupleList_dict[iAAIG].append((aatype, trans_prob))
                array_index += 1

        # Calculate and save the Conditional Probabilities for every aa type prediction
        for i_AAIG in list(iAAIG_iminus1aaTypesProbTupleList_dict.keys()):
            aatypeCondprob_list = []  # list of the form [(aa type, conditional probability), (aa type, conditional probability), ...]
            for aa_type, prob in sorted(iAAIG_iminus1aaTypesProbTupleList_dict[i_AAIG], key=itemgetter(1), reverse=True):
                condprob = self.__calculate_TOCSY_Paaiminus1_TAAIGi(i_AAIG, aa_type)
                self.iAAIG_iminus1aaType_Condprob_mdict[i_AAIG][aa_type] = condprob

    def __populate_NOESY_AAIG_aaTypesProbTupleList_dict(self, noesy_spec):
        """
        This internal method invokes the NOESY AA-type classifier and gets all possible aa-type predictions for each
        NOESY AAIG of the spectrum along with their probabilities.
        :param noesy_spec:
        :return:
        """
        ColorPrint("Loading the AA-type Classifier and predicting probabilities for all NOESY AAIGs.", "BOLDGREEN")
        # load the NOESY Classifier here. TODO: FIX THE PATH!!!
        # with open("/home2/thomas/Documents/4D-CHAINS_NOESY_Classifier/single_frame_noesy_classifier_12_2019.pickle", 'rb') as infile:
        # with open("/home2/thomas/Documents/4D-CHAINS_NOESY_Classifier/single_frame_noesy_classifier_12_2019_lite.pickle", 'rb') as infile:
        with open(self.classifier_file, 'rb') as infile:
            model = pickle.load(infile)
        AAIG_aaTypesProbTupleList_dict = {}
        for aaig in noesy_spec.__get_all_AAIGs__():
            aaTypesProbTupleList = model.predict( aaig.__get_CH_dataframe_for_AAtype_prediction__() )
            AAIG_aaTypesProbTupleList_dict[aaig.signature] = [d for d in aaTypesProbTupleList if d[1]>0]    # keep only non-zero prob aa-types
        return AAIG_aaTypesProbTupleList_dict

    def calc_NOESY_aa_type_conditional_probs(self, noesy_spec, COMPLETE_AA_TYPES_FILE=None, transform_type="None"):
        """
        Method to either calculate AA-type probabilities from scratch with the ML Classifier, or to load them from
        a file. It then calculates Conditional AA-type probabilities to be used for peptide formation.

        :param noesy_spec:      NOESY spectrum object.
        :param transform_type:AAIG_aaTypes
        :return: populates self.AAIG_aaTypesProbTupleList_dict and self.AAIG_aaType_Condprob_mdict.
        """
        # TODO: make a separate method for noesy_spec and for aa-type probabilities file.
        if COMPLETE_AA_TYPES_FILE:
        # If a file is provided, load ALL AA-type prediction probabilities for speed
            self.AAIG_aaTypesProbTupleList_dict = \
                Probability.load_aatype_probs_from_file(fname=COMPLETE_AA_TYPES_FILE)
        else:
        # Otherwise calculate all probabilities P[AATYPE(I)] using the ML classifier
            self.AAIG_aaTypesProbTupleList_dict = self.__populate_NOESY_AAIG_aaTypesProbTupleList_dict(noesy_spec)

        ## DO A MATHEMATICAL TRANSFORM ON P (CURRENTLY DEACTIVATED IN NOESY)
        # iAAIG_iminus1aaTypesProbTupleList_dict is an ordereddict with keys the AAIG of residue i and values lists of tuples of the form
        # (residue i-1 matching aa type, average probability)

        # Collect all aa-type probabilities in one list and transform them
        all_prob_list = []
        for AAIG, aaTypesProbTuple_list in list(self.AAIG_aaTypesProbTupleList_dict.items()):
            for duplet in aaTypesProbTuple_list:
                prob = duplet[1]
                all_prob_list.append(prob)
        all_prob_array = np.array(all_prob_list)
        transformed_all_prob_array = Probability.transform(all_prob_array, transform_type=transform_type)

        # Populate the dict with transformed probabilities
        array_index = 0
        for AAIG, aaTypesProbTuple_list in list(self.AAIG_aaTypesProbTupleList_dict.items()):
            self.AAIG_aaTypesTransformedProbTupleList_dict[AAIG] = []
            for duplet in aaTypesProbTuple_list:
                aatype = duplet[0]
                trans_prob = transformed_all_prob_array[array_index]
                self.AAIG_aaTypesTransformedProbTupleList_dict[AAIG].append((aatype, trans_prob))
                array_index += 1

        # Calculate and save the Conditional Probabilities for every aa type prediction
        for AAIG in list(self.AAIG_aaTypesProbTupleList_dict.keys()):
            aatypeCondprob_list = []  # list of the form [(aa type, conditional probability), (aa type, conditional probability), ...]
            for aa_type, prob in sorted(self.AAIG_aaTypesProbTupleList_dict[AAIG], key=itemgetter(1), reverse=True):
                condprob = self.__calculate_Paa_AAIG(AAIG, aa_type)
                self.AAIG_aaType_Condprob_mdict[AAIG][aa_type] = condprob

    @staticmethod
    def save_aa_type_predictions_to_file(AAIG_aaType_prob_dict,
                                         fname="amino_acid_type_prediction_conditional_probabilities",
                                         spectrum_type="TOCSY"):

        if spectrum_type == "TOCSY":
            header = "i AAIG\tpossible i-1 aa types\n"
        elif spectrum_type == "NOESY":
            header = "i AAIG\tpossible i aa types\n"

        with open(fname, 'w') as f:
            f.write(header)
            # print(header)

            # Save ALL amino acid type predictions with the respective conditional probabilities
            if type(AAIG_aaType_prob_dict) in [ OrderedDict, dict, defaultdict ]:
                for i_AAIG in list(AAIG_aaType_prob_dict.keys()):
                    # print(i_AAIG,"---> (i-1) aa type",sorted(iAAIG_iminus1aaTypesZscoreTupleList_dict[i_AAIG], key=itemgetter(1), reverse=True))
                    f.write(i_AAIG + " ---> (i-1) aa type " + str(sorted(AAIG_aaType_prob_dict[i_AAIG],
                                                                     key=itemgetter(1),
                                                                     reverse=True)) + "\n")
            elif type(AAIG_aaType_prob_dict) == tree:   # its the self.AAIG_aaType_Condprob_mdict
                for i_AAIG in list(AAIG_aaType_prob_dict.keys()):
                    aa_type_condprob_list = []  # list of the form [(aa type, conditional probability), (aa type, conditional probability), ...]
                    for aa_type in list(AAIG_aaType_prob_dict[i_AAIG].keys()):
                        aa_type_condprob_list.append((aa_type, AAIG_aaType_prob_dict[i_AAIG][aa_type]))
                    # print(i_AAIG, "---> (i-1) aa type", sorted(aa_type_condprob_list, key=itemgetter(1), reverse=True))
                    f.write(i_AAIG + " ---> (i-1) aa type " + str(sorted(aa_type_condprob_list,
                                                                       key=itemgetter(1),
                                                                       reverse=True)) + "\n")

    def get_pruned_condprob_dict(self, prob_cutoff=0.1):
        AAIG_aaTypesCutoffProbTupleList_dict = OrderedDict()
        for AAIGsign, probTuple_list in self.AAIG_aaTypesTransformedProbTupleList_dict.items():
            AAIG_aaTypesCutoffProbTupleList_dict[AAIGsign] = [d for d in probTuple_list if d[1]>=prob_cutoff]
        return AAIG_aaTypesCutoffProbTupleList_dict

    @staticmethod
    def get_pruned_prob_dicts_by_zscore(AAIG_aaTypesProbPoolTupleList_dict,
                                        ZSCORE_ASSIGNMENT_CUTOFF=-1.0,
                                        DELETE_AA_TYPE_PREDICTIONS=False,
                                        LOG_TRANSFORM=False,
                                        MIN_NUM_OF_PREDICTIONS=4):

        # # TODO: fix it, THIS IS NOESY-specify
        # if not AAIG_aaTypesProbPoolTupleList_dict:
        #     # In this case I will assume that the whole set of AA-type predictions and probabilities should be used.
        #     AAIG_aaTypesProbPoolTupleList_dict = self.AAIG_aaTypesTransformedProbTupleList_dict

        ## Now Calculate Z-scores from the pool of aa type prediction probabilities and keep only those predictions above the cutoff
        minimum_Zscore = 10000000
        AAIG_aaTypesZscoreTupleList_dict = OrderedDict()  # ordereddict with keys the AAIG of residue i and values lists of tuples of the form (residue i-1 matching aa type, Z-score)
        AAIG_aaTypesCutoffProbTupleList_dict = OrderedDict()  # contains only the tuples (residue i-1 matching aa type, probability) with Z-score above the cutoff
        for iAAIG, iminus1aaTypesProbTuple_list in list(AAIG_aaTypesProbPoolTupleList_dict.items()):
            # DECIDE THE PROBABILITY THRESHOLD TO DISCARD AA TYPE PREDICTIONS BELOW IT
            try:
                max_prob = np.max([duplet[1] for duplet in iminus1aaTypesProbTuple_list])  # highest aa type probability
            except ValueError:  # if no aa type predictions exist for this TAAIG group
                max_prob = 0
            if max_prob > 10e-10:  # Dedice the probability threshold
                prob_threshold = 1000
            elif max_prob > 10e-20:
                prob_threshold = 10000
            elif max_prob <= 10e-20:
                prob_threshold = 100000
            print("DEBUG: max_prob", max_prob, "prob_threshold", prob_threshold, "LOG_TRANSFORM=", LOG_TRANSFORM)
            aatype_list, prob_list = [], []
            for duplet in iminus1aaTypesProbTuple_list:
                if DELETE_AA_TYPE_PREDICTIONS == True:  # delete low probability aa type predictions
                    try:
                        if max_prob / float(duplet[1]) < prob_threshold:
                            aatype_list.append(duplet[0])
                            prob_list.append(duplet[1])  # unsorted list of probabilities
                    except ZeroDivisionError:
                        print("DEBUG: ZeroDivisionError ", max_prob, "/", float(duplet[1]))
                elif DELETE_AA_TYPE_PREDICTIONS == False:
                    aatype_list.append(duplet[0])
                    prob_list.append(duplet[1])  # unsorted list of probabilities
            print("DEBUG: prob_list=", prob_list)
            if len(prob_list) == 1:
                zscore_array = np.array([10.00000])  # default Z-score in the case of just one prediction
            elif len(prob_list) > 1:
                if LOG_TRANSFORM == True:
                    try:
                        ratio = float(np.max(prob_list)) / np.min(
                            prob_list)  # if the min probability is at least 3 orders of magnitude smaller, convert them to logarithmic scale
                    except ZeroDivisionError:
                        print("DEBUG: ZeroDivisionError ", np.max(prob_list), "/", np.min(prob_list))
                        ratio = -1.0
                    if ratio > 1000:
                        zscore_array = zscore(np.log(prob_list))  # convert probabilities the logarithms and then to Z-scores
                    else:
                        zscore_array = zscore(prob_list)  # convert probabilities to Z-scores
                elif LOG_TRANSFORM == False:
                    zscore_array = zscore(prob_list)  # convert probabilities to Z-scores
            else:  # if no aa type prediction was made create an empty array
                zscore_array = np.array([])
            # print "DEBUG1: saving zscore_array = ", zscore_array
            iminus1aaTypesZscoreTuple_list = []
            iminus1aaTypesCutoffProbTuple_list = []
            for aatype, Zscore, prob in zip(aatype_list, zscore_array, prob_list):
                if Zscore < minimum_Zscore:
                    minimum_Zscore = Zscore
                ## ONLY IF THE Z-SCORE OF THE AA TYPE PREDICTION IS GREATER THAN THE CUTOFF AND THE AA TYPE PREDICTIONS ARE AT LEAST args.MIN_NUM_OF_PREDICTIONS, THEN INCLUDE IT IN THE LIST
                if Zscore > ZSCORE_ASSIGNMENT_CUTOFF and len(zscore_array) >= MIN_NUM_OF_PREDICTIONS:
                    aatypeProb_tuple = (aatype, prob)
                    aatypeZscore_tuple = (aatype, Zscore)  # replace the probability with the respective Z-score
                    iminus1aaTypesZscoreTuple_list.append(aatypeZscore_tuple)
                    iminus1aaTypesCutoffProbTuple_list.append(aatypeProb_tuple)  # but also save the aatype probabilities that satisfy the Z-score cutoff
                ## IF THE AA TYPE PREDICTIONS ARE LESS THAN MIN_NUM_OF_PREDICTIONS, THEN INCLUDED ALL OF THEM IN THE LIST
                elif len(zscore_array) < MIN_NUM_OF_PREDICTIONS:
                    aatypeProb_tuple = (aatype, prob)
                    aatypeZscore_tuple = (aatype, Zscore)  # replace the probability with the respective Z-score
                    iminus1aaTypesZscoreTuple_list.append(aatypeZscore_tuple)
                    iminus1aaTypesCutoffProbTuple_list.append(aatypeProb_tuple)  # but also save the aatype probabilities that satisfy the Z-score cutoff
            AAIG_aaTypesZscoreTupleList_dict[iAAIG] = iminus1aaTypesZscoreTuple_list
            AAIG_aaTypesCutoffProbTupleList_dict[iAAIG] = iminus1aaTypesCutoffProbTuple_list

        return AAIG_aaTypesCutoffProbTupleList_dict, AAIG_aaTypesZscoreTupleList_dict

    # Private method to the class (it is used only internally)
    def __calculate_TOCSY_Paaiminus1_TAAIGi(self, TAAIGi, aaiminus1):
        """
            This method was initially writen for TOCSY-HCNH assignment, hence the variable names TAAIG, aaiminu1.
            Calculate the conditional probability P(AA^{u} | AAIG^{CS}), which quantifies how probable is to get an
            amino-acid of type u given a set of chemical shifts AAIG^{CS}.
            ARGS:
            TAAIGi:    the AAIG or TOCSY Index Group aligned to position i
            aaiminus1:  the aa type if the residue in position i-1
        """
        # aatype_P_dict: probability P[aatype(i-1)], dict with the P[aatype(i-1)] for every aatype(i-1)
        # aatype_P_dict = shared.getConst('AATYPE_P_DICT')
        # iAAIG_iminus1aaTypesTransformedProbTupleList_dict = shared.getConst('IAAINDEX_IMINUS1AATYPESTRANSFORMEDPROBTUPLELIST_DICT')

        PTAAIGi = 0    # P[TAAIG(i)] = Sum{P[TAAIG(i)|aatype(i-1)]}
        for duplet in self.AAIG_aaTypesTransformedProbTupleList_dict[TAAIGi]:
            #print "DEBUG: duplet=", duplet
            aatype = duplet[0]
            PTAAIGi += duplet[1]
            if aatype == aaiminus1:
                PTAAIGi_aaiminus1 = duplet[1]  # this is P[TAAIG(i)|aatype(i-1)]

        #print "DEBUG: aaiminus1=",aaiminus1
        #print "DEBUG: self.aatype_P_dict=", self.aatype_P_dict
        Paaiminus1 = self.aatype_P_dict[aaiminus1]  # get the frequency of this aa-type in the protein
        #print "DEBUG: PTAAIGi_aaiminus1=",PTAAIGi_aaiminus1
        #print "DEBUG: Paaiminus1=",Paaiminus1
        #print "DEBUG: PTAAIGi=", PTAAIGi
        Paaiminus1_TAAIGi = PTAAIGi_aaiminus1 * Paaiminus1 / PTAAIGi
        #print "DEBUG: Paaiminus1_TAAIGi=", Paaiminus1_TAAIGi
        return Paaiminus1_TAAIGi

    # Private method to the class (it is used only internally)
    def __calculate_NOESY_Paa_AAIG(self, AAIG, aa):
        """
        Calculate the conditional probability P(AA^{u} | AAIG^{CS}) for NOESY, which quantifies how probable is to get an
        amino-acid of type u given a set of chemical shifts AAIG^{CS}.

        The NOESY AA-type classifier returns the joint probability P(AA, AAIG) = P(AAIG|AA) * P(AA), unlike the Histograms,
        which return P(AAIG|AA). Therefore, we need to divide the transformed probabilities by P(AA) to get P(AAIG|AA).

        :param AAIG: the AAIG or TOCSY Index Group aligned to position i
        :param aa: the aa type if the residue in position i-1
        :return:
        """
        # aatype_P_dict: probability P[aatype(i-1)], dict with the P[aatype(i-1)] for every aatype(i-1)
        # aatype_P_dict = shared.getConst('AATYPE_P_DICT')
        # AAIG_aaTypesTransformedProbTupleList_dict = shared.getConst('AAIG_AATYPESTRANSFORMEDPROBTUPLELIST_DICT')

        PAAIG = 0    # P[TAAIG(i)] = Sum{P[TAAIG(i)|aatype(i-1)]}
        for duplet in self.AAIG_aaTypesTransformedProbTupleList_dict[AAIG]:
            aatype = duplet[0]
            Paa = self.aatype_P_dict[aatype]  # get the frequency of aatype in the protein
            print("DEBUG: AA=%s , P(AA, AAIG)/P(AA) =" % aatype, duplet[1], "/", Paa)
            if Paa == 0:    # TEMPORARY FIX
                PAAIG += 0  # duplet[1] is P(AA, AAIG) = P(AAIG|AA) * P(AA), that's why we divide by Paa
                if aatype == aa:
                    PAAIG_aa = 0  # this is P[AAIG(i)|aatype(i)]
            else:
                PAAIG += duplet[1]/Paa  # duplet[1] is P(AA, AAIG) = P(AAIG|AA) * P(AA), that's why we divide by Paa
                if aatype == aa:
                    PAAIG_aa = duplet[1]/Paa  # this is P[AAIG(i)|aatype(i)]

        Paa = self.aatype_P_dict[aa]  # get the frequency of this aa-type in the protein
        Paa_AAIG = PAAIG_aa * Paa / PAAIG
        return Paa_AAIG

    # For NOESY name compatibility, make a private copy of the above method.
    def __calculate_Paa_AAIG(self, AAIG, aa):
        return self.__calculate_NOESY_Paa_AAIG(AAIG, aa)

    @staticmethod
    def load_aatype_probs_from_file(fname):

        print("")  # change line after the previous progress bar
        print("Loading amino acid prediction files ...")
        AAIG_aaTypesProbTupleList_dict = OrderedDict()  # ordereddict with keys the AAIG of residue i and values lists of tuples of the form

        # (residue i-1 matching aa type, Z-score)
        with open(fname, 'r') as f:
            aa_type_file_contents = f.readlines()
            for line in aa_type_file_contents[1:]:
                # print "DEBUG: line=",line
                word_list = re.sub('---> \(i-1\) aa type', '', line).split()
                key = word_list[0]
                AAIG_aaTypesProbTupleList_dict[key] = []
                values_string = ''.join(word_list[1:])
                elements_string = re.sub('[\(\)\[\]\'\"]', '', values_string).split("),(")[0]
                elements_list = elements_string.split(",")
                for aa, Zscore in zip(elements_list[0::2], elements_list[1::2]):
                    # print "aa=",aa,"probability=",probability
                    duplet = (aa, float(Zscore))
                    AAIG_aaTypesProbTupleList_dict[key].append(duplet)
        return AAIG_aaTypesProbTupleList_dict
