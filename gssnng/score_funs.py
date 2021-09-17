import numpy as np
from scipy import sparse
import pandas as pd
import gssnng.util as si
import statsmodels.robust.scale


def scorefun(x, su, sig_len_up, norm_method, score_up, ):
    """
    given a ranked list, produce a score
    TODO add more functions like median etc

    :param x: the pandas data frame of ranks
    :param su: the ranked list
    :param sig_len_up: the number of genes matched in the set
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


def _ms_sing(geneset: list, x: pd.Series, norm_method: str, rankup: bool, dorank: bool) -> dict:
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

    total_score, mad_up = scorefun(x, su, sig_len_up, norm_method, rankup)
    return dict(total_score=total_score, mad_up=mad_up)
