import numpy as np
from scipy import sparse
import pandas as pd
import gssnng.util as si
import statsmodels.robust.scale
from gssnng.genesets import geneset

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


def median_score(su):
    """
    Average Z score

    :param x: the pandas data frame of ranks, all genes
    :param su: the ranked list of genes *IN* the gene set
    :param sig_len_up: the number of expressed genes matched in the set
    :param norm_method: 'standard or theoretical' # from singscore
    :param score_up: is the rank up or down?  True or False
    """
    # normalise the score for the number of genes in the signature
    cnts_med = np.median(su)
    mad_su = statsmodels.robust.scale.mad(su)
    return((cnts_med, mad_su))


def average_score(su):
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
    return((cnts_mean, std_su))


def mean_z(df, su):
    """
    Average Z score

    :param x: the pandas data frame of ranks, all genes
    :param su: the ranked list of genes *IN* the gene set
    :param sig_len_up: the number of expressed genes matched in the set
    :param norm_method: 'standard or theoretical' # from singscore
    :param score_up: is the rank up or down?  True or False
    """
    # normalise the score for the number of genes in the signature
    cnts_mean = np.mean(df)
    std_su = np.std(df)
    centered_cnts = [ ((x - cnts_mean) / std_su) for x in su ]
    score = np.mean(centered_cnts)
    return((score, std_su))


def robust_std(df, su):
    """
    Median of median standardized counts

    :param x: the pandas data frame of ranks, all genes
    :param su: the ranked list of genes *IN* the gene set
    :param sig_len_up: the number of expressed genes matched in the set
    :param norm_method: 'standard or theoretical' # from singscore
    :param score_up: is the rank up or down?  True or False
    """
    # normalise the score for the number of genes in the signature
    cnts_med = np.median(df)
    mad_su = statsmodels.robust.scale.mad(df)
    centered_cnts = [ ((x - cnts_med) / mad_su) for x in su ]
    score = np.median(centered_cnts)
    return((score, mad_su))


def rank_biased_overlap(x, su, gs, rankup, limit=100):
    """
    Rank biased overlap method

    :param x: the pandas data frame of ranks, all genes
    :param su: the ranked list of genes *IN* the gene set
    :param gs: the gene set
    """

    rbo_score = 0.0

    sortorder = not rankup  # if it's rank up, we want ascending false

    y = x.sort_values(ascending=sortorder)

    for i,gi in enumerate(gs):
        subrank = y[0:(i+1)]
        subset = set(subrank.index)
        rbo_score += len(subset.intersection(gs))
        if i > limit:
            break

    mad_up = statsmodels.robust.scale.mad(su)
    return((rbo_score, mad_up))


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



def scorefun(gs, x, su, sig_len_up, norm_method, score_up, method, rbo_depth):
    """
    given a ranked list, produce a score

    :param gs: the gene set
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
        res0 = robust_std(x, su)

    if method == 'summed_up':
        res0 = summed_up(su)

    if method == 'median_score':
        res0 = median_score(x, su)

    if method == 'average_score':
        res0 = average_score(x, su)

    if method == 'mean_z':
        res0 = mean_z(x, su)

    if method == 'rank_biased_overlap':
        res0 = rank_biased_overlap(x, su, gs, score_up, rbo_depth)

    return(res0)


def _ms_sing(geneset: geneset,
             x: pd.Series,
             score_method: str,
             norm_method: str,
             rankup: bool,
             dorank: bool,
             rbo_depth: int,
             ) -> dict:
    """
    bare bones version of scsing scoring. Their function (see scsingscore.py)
    does a ton of stuff, here's the essentials

    :param genest: Geneset to score against
    :param x: pd.Series with the gene expression of a single sample. One gene per row
    :param norm_method: how to normalize the scores
    :param rankup: direction of ranking, up: True, down: False
    """

    geneset_genes = geneset.genes_up

    sig_len_up = len(geneset_genes)

    exprdat = x.uprank

    assert isinstance(x, pd.DataFrame)

    su = []

    # for every gene in the list gene get the value at that
    # index/rowname (the gene) and the sample that is equal to i
    for j in geneset_genes:
        if j in exprdat.index:
            su.append(exprdat[j])
        else:
            sig_len_up = sig_len_up - 1

    total_score, variance = scorefun(geneset_genes, exprdat, su, sig_len_up, norm_method, rankup, score_method, rbo_depth)
    return dict(score=total_score, var=variance)
