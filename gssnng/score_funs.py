import numpy as np
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
    score = np.sum(su)
    return(score)


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
    return(cnts_med)


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
    return(cnts_mean)


def mean_z(allexprvals, genesetvals):
    """
    Average Z score

    :param exprdat: the pandas data frame of ranks, all genes
    :param su: the ranked list of genes *IN* the gene set
    :param sig_len_up: the number of expressed genes matched in the set
    :param norm_method: 'standard or theoretical' # from singscore
    :param score_up: is the rank up or down?  True or False
    """
    # normalise the score for the number of genes in the signature
    vals_mean = np.mean(allexprvals)
    vals_std = np.std(allexprvals)
    centered = [ (np.abs(x - vals_mean) / vals_std) for x in genesetvals ]
    score = np.mean(centered)
    return(score)


def robust_std(exprdat, su):
    """
    Median of median standardized counts

    :param exprdat: the pandas data frame of ranks, all genes
    :param su: the ranked list of genes *IN* the gene set
    :param sig_len_up: the number of expressed genes matched in the set
    :param norm_method: 'standard or theoretical' # from singscore
    :param score_up: is the rank up or down?  True or False
    """
    # normalise the score for the number of genes in the signature
    cnts_med = np.median(exprdat)
    mad_su = statsmodels.robust.scale.mad(exprdat)
    centered_cnts = [ (np.abs(x - cnts_med) / mad_su) for x in su ]
    score = np.median(centered_cnts)
    return(score)


def rank_biased_overlap(x, exprcol, gs, geneset_genes, limit):
    """
    Rank biased overlap method

    :param x: the pandas data frame of ranks, all genes
    :param su: the ranked list of genes *IN* the gene set
    :param gs: the gene set
    """
    rbo_score = 0.0

    # if undirected, then sort based on centered values
    if gs.mode == '?':
        # center & absolute value ranks
        maxN = np.ceil(len(x.index)/2.0)
        x['undir'] = [ np.abs(xi - maxN) for xi in x[exprcol]]
        exprcol = 'undir'

    # get sorted dataframe
    x_sorted = x.sort_values(by=exprcol, ascending=False)
    y = x_sorted[exprcol]

    for i in range(limit):
        subset = set(y.index[0:(i+1)])
        rbo_score += len(subset.intersection(geneset_genes))

    return( rbo_score )


def singscore(x, su, sig_len, norm_method, gs):
    """
    The singscore method

    :param x: the pandas data frame of ranks, all genes
    :param su: the ranked list of genes *IN* the gene set
    :param sig_len_up: the number of expressed genes matched in the set
    :param norm_method: 'standard or theoretical' # from singscore
    :param score_up: is the rank up or down?  True or False
    :param gs: gene set object
    """
    # normalise the score for the number of genes in the signature
    if gs.mode == '?':
        # center & absolute value ranks
        maxN = np.ceil(len(x.index)/2.0)
        su = [ np.abs(xi - maxN) for xi in su]

    mean_rank = np.mean(su)
    norm_up = si.normalisation(norm_method=norm_method,
                               gs_mode=gs.mode,
                               score=mean_rank,
                               library_len=len(x.index),
                               sig_len=sig_len)

    if gs.mode != '?':
        norm_up = norm_up - 0.5

    return(norm_up)


def ssgsea(x, su, sig_len, omega, gs):
    """
    The ssGSEA method

    :param x: the pandas data frame of ranks, all genes
    :param su: the ranked list of genes *IN* the gene set
    :param sig_len_up: the number of expressed genes matched in the set
    :param norm_method: 'standard or theoretical' # from singscore
    :param score_up: is the rank up or down?  True or False
    :param gene_set: gene set object

    """
    gene_set = set(gs)

    #first sort by absolute expression value, starting with the highest expressed genes first
    xsorted = x.sort_values(axis=0, ascending=False, inplace=False)
    keys_sorted = xsorted.index.tolist()

    #values representing the ECDF of genes in the geneset
    P_GW_numerator = 0
    P_GW_denominator = 0

    #determining denominator value
    i = 1 #current rank stepping through listing of sorted genes
    for gene in keys_sorted:
        if gene in gene_set:
            P_GW_denominator += i ** omega
        i += 1

    P_GW = lambda : P_GW_numerator / P_GW_denominator

    #values representing the ECDF of genes not in the geneset
    P_NG_numerator = 0
    P_NG_denominator = len(x) - len(gene_set)
    P_NG = lambda : P_NG_numerator / P_NG_denominator

    #integrate different in P_GW and P_NG
    i = 1 #current index in the traversal of sorted genes
    scores = []
    for gene in keys_sorted:
        if gene in gene_set:
            P_GW_numerator += i ** omega
        else:
            P_NG_numerator += 1
        scores.append(P_GW() - P_NG())
        i += 1

    return sum(scores)



def expr_format(x, exprcol, geneset_genes):
    sig_len_up = len(geneset_genes)
    su = []
    for j in geneset_genes:
        if j in x.index:
            su.append(x[exprcol][j])
        else:
            sig_len_up = sig_len_up - 1
    return( (su, sig_len_up) )



def method_selector(gs, x, exprcol, geneset_genes, method, method_params):
    """
    :param gs: the gene set
    :param x: the gene expr data frame
    :param exprcol: the column containing values we'll compute on
    :param geneset_genes: genes in the gene set
    :param method: the method we'll call
    :param method_params: dictionary of method parameters
    :param barcode: cell barcode

    :return: dictionary of results
    """

    (su, sig_len) = expr_format(x, exprcol, geneset_genes)
    exprdat = x[exprcol]

    if method == 'singscore':
        res0 = singscore(exprdat, su, sig_len, method_params['normalization'], gs)

    elif method == 'robust_std':
        res0 = robust_std(exprdat, su)

    elif method == 'summed_up':
        res0 = summed_up(su)

    elif method == 'median_score':
        res0 = median_score(su)

    elif method == 'average_score':
        res0 = average_score(su)

    elif method == 'mean_z':
        res0 = mean_z(exprdat, su)

    elif method == 'rank_biased_overlap':
        res0 = rank_biased_overlap(x, exprcol, gs, geneset_genes, method_params['rbo_depth'])

    elif method == 'ssgsea':
        #x, su, sig_len, omega, gene_set
        res0 = ssgsea(exprdat, su, sig_len, method_params['omega'], geneset_genes)

    else:
        return(np.nan)

    return(res0)


def scorefun(gs,
             x,
             method,
             method_params,
             ranked):
    """
    given a ranked list, produce a score

    :param gs: the gene set
    :param x: the pandas data frame of ranks, all genes
    :param method: the method we'll call
    :param method_params: dictionary of method parameters
    :param ranked: ranked data? True | False

    :return a score
    """

    try:
        if (gs.mode == 'UP') and (ranked == False):
            res0 = method_selector(gs, x, 'counts', gs.genes_up, method, method_params)

        elif (gs.mode == 'DN') and (ranked == False):
            res0 = method_selector(gs, x, 'counts', gs.genes_dn, method, method_params)

        elif (gs.mode == 'BOTH') and (ranked == False):
            res0_up = method_selector(gs, x, 'counts', gs.genes_up, method, method_params)
            res0_dn = method_selector(gs, x, 'counts', gs.genes_dn, method, method_params)
            res0 = (res0_up + res0_dn)

        elif (gs.mode == '?') and (ranked == False):
            res0 = method_selector(gs, x, 'counts', gs.genes_up, method, method_params)

        elif (gs.mode == 'UP') and (ranked == True):
            res0 = method_selector(gs, x, 'uprank', gs.genes_up, method, method_params)

        elif (gs.mode == 'DN') and (ranked == True):
            res0 = method_selector(gs, x, 'dnrank', gs.genes_dn, method, method_params)

        elif (gs.mode == 'BOTH') and (ranked == True):
            res0_up = method_selector(gs, x, 'uprank', gs.genes_up , method, method_params)
            res0_dn = method_selector(gs, x, 'dnrank', gs.genes_dn, method, method_params)
            res0 = (res0_up + res0_dn)

        elif (gs.mode == '?') and (ranked == True):
            res0 = method_selector(gs, x, 'uprank', gs.genes_up, method, method_params)

    except ():
        #res1 = dict(barcode = barcode, name=gs.name, mode=gs.mode, score=np.nan, var=np.nan)
        res0 = np.nan

    return(res0)

