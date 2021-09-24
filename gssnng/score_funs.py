import numpy as np
from scipy import sparse
import pandas as pd
import gssnng.util as si
import statsmodels.robust.scale


def summed_up(su):
    """
    Just sum up the counts

    :param x: the pandas data frame of ranks, all genes
    :param su: the ranked list of genes *IN* the gene set
    :param sig_len_up: the number of expressed genes matched in the set
    :param norm_method: 'standard or theoretical' # from singscore
    :param score_up: is the rank up or down?  True or False
    """
    # normalise the score for the number of genes in the signature
    mad_su = statsmodels.robust.scale.mad(su)
    score = np.sum(su)
    return((score, mad_su))


def average(su):
    """
    Average Z score

    :param x: the pandas data frame of ranks, all genes
    :param su: the ranked list of genes *IN* the gene set
    :param sig_len_up: the number of expressed genes matched in the set
    :param norm_method: 'standard or theoretical' # from singscore
    :param score_up: is the rank up or down?  True or False
    """
    # normalise the score for the number of genes in the signature
    cnts_mean = np.mean(su)
    std_su = np.std(su)
    score = cnts_mean / len(su)
    return((score, std_su))


def mean_z(su):
    """
    Average Z score

    :param x: the pandas data frame of ranks, all genes
    :param su: the ranked list of genes *IN* the gene set
    :param sig_len_up: the number of expressed genes matched in the set
    :param norm_method: 'standard or theoretical' # from singscore
    :param score_up: is the rank up or down?  True or False
    """
    # normalise the score for the number of genes in the signature
    cnts_mean = np.mean(su)
    std_su = np.std(su)
    centered_cnts = [ ((x - cnts_mean) / std_su) for x in su ]
    score = np.mean(centered_cnts)
    return((score, std_su))


def robust_std(su):
    """
    Median of median standardized counts

    :param x: the pandas data frame of ranks, all genes
    :param su: the ranked list of genes *IN* the gene set
    :param sig_len_up: the number of expressed genes matched in the set
    :param norm_method: 'standard or theoretical' # from singscore
    :param score_up: is the rank up or down?  True or False
    """
    # normalise the score for the number of genes in the signature
    cnts_med = np.median(su)
    mad_su = statsmodels.robust.scale.mad(su)
    centered_cnts = [ ((x - cnts_med) / mad_su) for x in su ]
    score = np.median(centered_cnts)
    return((score, mad_su))


def singscore(x, su, sig_len_up, norm_method):
    """
    The singscore method

    :param x: the pandas data frame of ranks, all genes
    :param su: the ranked list of genes *IN* the gene set
    :param sig_len_up: the number of expressed genes matched in the set
    :param norm_method: 'standard or theoretical' # from singscore
    :param score_up: is the rank up or down?  True or False
    """
    # normalise the score for the number of genes in the signature
    score_up = np.mean(su)
    norm_up = si.normalisation(norm_method=norm_method,
                               library_len=len(x.index),
                               score_list=su,
                               score=score_up,
                               sig_len=sig_len_up)
    norm_up = norm_up - 0.5
    mad_up = statsmodels.robust.scale.mad(su)
    return((norm_up, mad_up))



def scorefun(x, su, sig_len_up, norm_method, score_up, method='singscore'):
    """
    given a ranked list, produce a score

    :param x: the pandas data frame of ranks, all genes
    :param su: the ranked list of genes *IN* the gene set
    :param sig_len_up: the number of expressed genes matched in the set
    :param norm_method: 'standard or theoretical' # from singscore
    :param score_up: is the rank up or down?  True or False

    :return a tuple (a,b)  score is 'a', extra info is in 'b'
    """
    if method == 'singscore':
        res0 = singscore(x, su, sig_len_up, norm_method)

    if method == 'robust_std':
        res0 = robust_std(su)

    if method == 'summed_up':
        res0 = summed_up(su)

    if method == 'average':
        res0 = average(su)

    if method == 'meanz':
        res0 = mean_z(su)

    return(res0)


def _ms_sing(geneset: list,
             x: pd.Series,
             score_method: str,
             norm_method: str,
             rankup: bool,
             dorank: bool,
             ) -> dict:
    """
    bare bones version of scsing scoring. Their function (see scsingscore.py)
    does a ton of stuff, here's the essentials

    :param genest: Geneset to score against
    :param x: pd.Series with the gene expression of a single sample. One gene per row
    :param norm_method: how to normalize the scores
    :param rankup: direction of ranking, up: True, down: False
    """

    sig_len_up = len(geneset)
    assert isinstance(x, pd.Series)
    if dorank:
        up_sort = x.rank(method='min', ascending=rankup)  #
    else:
        up_sort = x
    su = []

    # for every gene in the list gene get the value at that
    # index/rowname (the gene) and the sample that is equal to i
    if True:
        for j in geneset:
            if j in up_sort.index:
                su.append(up_sort[j])
            else:
                sig_len_up = sig_len_up - 1
    else:
        # dict acces would be faster, but dict generation takes too loading
        # damn
        d = up_sort.to_dict()
        for g in geneset:
            if g in d:
                su.append(d[g])
            else:
                sig_len_up = sig_len_up - 1

    total_score, variance = scorefun(x, su, sig_len_up, norm_method, rankup, score_method)
    return dict(score=total_score, var=variance)
