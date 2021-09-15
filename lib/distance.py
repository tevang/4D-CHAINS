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


from lib.global_func import *
from sklearn.preprocessing import normalize, StandardScaler
from scipy.spatial.distance import cdist

# Scalar vector distance types (without silly ones, namely 'chebyshev'):
SCIPY_SCALAR_DISTANCE_TYPES = ['braycurtis', 'canberra', 'cityblock', 'correlation', 'cosine', 'euclidean',
'jensenshannon', 'minkowski', 'seuclidean', 'sokalsneath', 'sqeuclidean'] # 'mahalanobis' is very slow

def calc_norm_dist_matrix(mat1, mat2, disttype):
    """
    Method that takes two 2D matrices, calculates all pairwise scalar vector distances and scales
    them to [0,1].
    :param mat1: N samples x M features (by default the crossval set). Can be a list of a numpy array.
    :param mat2: K samples x M features (by default the xtest set). Can be a list of a numpy array.
    :param disttype: distance type
    :return ndistmat:   N x K normalized distance matrix
    """
    assert disttype in SCIPY_SCALAR_DISTANCE_TYPES, ColorPrint("ERROR: invalid distance type: %s" %
                                                               disttype, "FAIL")
    distmat = cdist(mat1, mat2, disttype)  # N x K distance matrix
    ndistmat = normalize(distmat)   # scale between [0,1]
    return ndistmat

def ensemble_maxvariance_sim_matrix(mat1, mat2, is_distance=False):
    """
    Method that takes two 2D matrices, calculates multiple types of pairwise scalar vector distances,
    scales them to [0,1], finds the variance of each distance type, and returns the -1*z-scores
    of the distance type with the highest variance.

    RECOMMENDED FOR FINDING WHICH TRAINING COMPOUNDS ARE SIMILAR TO THE TEST COMPOUND (~DOMAIN OF APPLICABILITY).

    :param mat1: N samples x M features (by default the crossval set). Can be a list of a numpy array.
    :param mat2: K samples x M features (by default the xtest set). Can be a list of a numpy array.
    :param is_distance: whether the returned matrix will contain distances of similarities.
    :return scaled_distmat: N x K z-score matrix of distances or similarities of the distance type with
                            the highest variance
    """
    disttype_ndistmat_dict = {}
    mat1 = np.array(mat1)
    mat2 = np.array(mat2)
    assert mat1.shape[1] == mat2.shape[1], \
        Debuginfo("ERROR: mat1 and mat2 do not have the same number of columns (%i vs %i)!" %
                                               (mat1.shape[1], mat2.shape[1]), fail=True)

    # Create the N x K distance matrix
    disttype_variance_list = []
    for disttype in SCIPY_SCALAR_DISTANCE_TYPES:
        ColorPrint("Calculating %s distances." % disttype, "OKBLUE")
        try:
            ndistmat = calc_norm_dist_matrix(mat1, mat2, disttype)  # N x K normalized distance matrix
        except Exception as e:
            print(repr(e))
            ColorPrint("Distance type %s won't be considered." % disttype, "OKRED")
            continue
        except np.linalg.LinAlgError:
            ColorPrint("Singular matrix. %s distances won't be calculated." % disttype, "OKRED")
            continue
        variance = ndistmat.var()
        ColorPrint("\t\t\t\tvariance=%f." % (variance), "OKBLUE")
        disttype_ndistmat_dict[disttype] = ndistmat
        disttype_variance_list.append( (disttype, variance) )    # the overall variance

    disttype_variance_list.sort(key=itemgetter(1), reverse=True)    # sort in descending order
    best_disttype = disttype_variance_list[0][0]    # keep the distance type with the maximum variance
    max_variance = disttype_variance_list[0][1]
    ColorPrint("The distance type with the max variance (%f) was %s" % (max_variance, best_disttype), "BOLDBLUE")

    scaler = StandardScaler()
    scaler.fit(disttype_ndistmat_dict[best_disttype].flatten().reshape(-1, 1))

    if is_distance:
        scaled_distmat = scaler.transform(disttype_ndistmat_dict[best_disttype])  # return z-score distances
    else:
        scaled_distmat = -1 * scaler.transform(disttype_ndistmat_dict[best_disttype])  # return similarities (-1*z-score distances)

    return scaled_distmat