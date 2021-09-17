"""
MS-DOS version
"""
import numpy as np
from scipy import sparse
import pandas as pd
import gssnng.util as si
import tqdm
import statsmodels.robust.scale
from anndata import AnnData
from gssnng.smoothing import nn_smoothing
from gssnng.util import read_gene_sets, error_checking
from gssnng.score_funs import _ms_sing


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
            s = _ms_sing(gene_list, df_noise['gene_counts'], norm_method='standard', rankup=True, dorank=False)
            s['CB'] = smoothed_adata.obs.index[cell_ix]
            results[gs_i] = s
        results_list.append( results )

    return(results_list)
