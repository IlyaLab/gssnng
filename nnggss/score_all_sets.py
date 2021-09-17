"""
MS-DOS version
"""
import numpy as np
from scipy import sparse
import pandas as pd
import nnggss.util as si
import tqdm
import statsmodels.robust.scale
from anndata import AnnData
from nnggss.smoothing import nn_smoothing
from nnggss.util import read_gene_sets, error_checking


def score_cells_all_sets_up(
        adata=None,
        gene_set_file=None,
        direction='Up',
        key_added='GeneSetNames',
        samp_neighbors=5,
        noise_trials=0,
        mode='average'
):

    """
    gene set scoring (all gene sets in file) with nearest neighbor smoothing of the expression matrix

    Improved single cell scoring by:
    - smoothing the data matrix
        - adding noise to the nearest neighbor smoothing via `samp_neighbors`
    - adding noise to the expression data itself (via noise_trials)

    :param adata: anndata.AnnData containing the cells to be scored
    :param gene_set_file: the gene set file with list of gene sets, gmt, one per line
    :param direction: gene sets will be scored in the 'up' direction
    :param key_added: name given to the new entry of adata.obs['key_added'] -- taken from file
    :param samp_neighbors: number of neighbors to sample
    :param noise_trials: number of noisy samples to create, integer
    :param mode: average or theoretical normalization of scores

    :returns: sparse matrix of scores, one per gene set and per cell in adata
    """
    if error_checking(adata, samp_neighbors) == 'ERROR':
        return()

    gs = read_gene_sets('data/h.all.v7.2.symbols.gmt')

    # NOTE: this is cells x genes
    smoothed_matrix = nn_smoothing(adata.X, adata, 'connectivity', samp_neighbors)
    # for easier handling with gene names
    smoothed_adata = AnnData(smoothed_matrix, obs=adata.obs, var=adata.var)
    """
    since we're doing all cells at the same time now,
    the following gets probelmatic (the df kills the sparsity)
    Ideas:
    - loop over cells, do the scoring
    - batch the cells, i.e. create a df with 100 cells (puling them into mem) and score those in one go
    """
    all_scores = _score_all_cells_all_sets(gene_set_dict=gs, smoothed_adata=smoothed_adata,
                                           noise_trials=noise_trials, mode=mode)

    for gs_name in gs.keys():
        gs_scores = [x[gs_name]['total_score'] for x in all_scores]
        adata.obs[gs_name] = gs_scores
    return()




def _ms_sing_norank(geneset: list, up_sort: pd.Series, norm_method: str, rankup: bool) -> dict:
    """
    bare bones version of scsing scoring. Their function (see scsingscore.py)
    does a ton of stuff, here's the essentials

    :param genest: Geneset to score against
    :param x: pd.Series with the gene expression of a single sample. One gene per row
    :param norm_method: how to normalize the scores
    :param rankup: direction of ranking, up: True, down: False
    """

    sig_len_up = len(geneset)
    assert isinstance(up_sort, pd.Series)
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

    # normalise the score for the number of genes in the signature
    score_up = np.mean(su)
    norm_up = si.normalisation(norm_method=norm_method,
                               library_len=len(up_sort.index),
                               score_list=su,
                               score=score_up,
                               sig_len=sig_len_up)
    norm_up = norm_up - 0.5
    mad_up = statsmodels.robust.scale.mad(su)
    total_score = norm_up
    return dict(total_score=total_score, mad_up=mad_up)


def get_ranked_cells(smoothed_adata, cell_ix, noise_trials, mode):
    """
    rank the expression data for each cell and return a list of data frames

    :param smoothed_adata:
    :param noise_trials:
    :return: a list of data frames, one per cell
    """
    # for each cell, rank the expression
    gene_mat = smoothed_adata.X[cell_ix]
    # then we subset it to only the genes with counts
    _, gdx, _ = sparse.find(gene_mat)
    # TODO we could do a dict instead of the df, that would be faster in _mssing too
    if gene_mat.ndim == 2:
        df = pd.DataFrame(gene_mat[:, gdx].A.flatten(), index=smoothed_adata.var.index[gdx])  ## ????
    else:
        df = pd.DataFrame(gene_mat[gdx],
                          index=smoothed_adata.var.index[gdx])  ## not sure why it's coming off as an array
    df.columns = ['gene_counts']
    # TODO not clear that adding noise really helps.
    if mode == 'average' and noise_trials > 0:
        # add some noise to gene counts.. create a n numbers of examples
        raise ValueError('not implemented')
        df_noise = si.add_noise(df, noise_trials, 0.01, 0.99)  ## slow part .. fixed
    else:
        df_noise = df
    up_sort = df_noise.rank(method='min', ascending=True)  #
    return(up_sort)


def _score_all_cells_all_sets(gene_set_dict=None, smoothed_adata=None, noise_trials=0, mode='average'):
    """
    Want to rank each cell once, but score all sets, and return a sparse matrix of scores.
    """

    results_list = []  # one per cell
    for cell_ix in tqdm.trange(smoothed_adata.shape[0]):
        results = dict()  # one score per cell in a list
        df_noise = get_ranked_cells(smoothed_adata, cell_ix, noise_trials, mode)
        for gs_i in gene_set_dict.keys():
            gene_list = gene_set_dict[gs_i]
            s = _ms_sing_norank(gene_list, df_noise['gene_counts'], norm_method='standard', rankup=True)
            s['CB'] = smoothed_adata.obs.index[cell_ix]
            results[gs_i] = s
        results_list.append( results )

    return(results_list)
