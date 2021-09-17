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
