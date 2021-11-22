import anndata
import numpy as np
from scipy import sparse
import pandas as pd
import tqdm
from anndata import AnnData
from gssnng.smoothing import nn_smoothing
from gssnng.util import error_checking
from gssnng.score_funs import scorefun
from gssnng.genesets import genesets


def score_cells_all_sets(
        adata: anndata.AnnData,
        gene_set_file: str,
        score_method: str,
        method_params: dict,
        samp_neighbors: int,
        noise_trials: int,
        keys_added: list
    ) -> anndata.AnnData:

    """
    gene set scoring (all gene sets in file) with nearest neighbor smoothing of the expression matrix

    Improved single cell scoring by:
    - smoothing the data matrix
        - adding noise to the nearest neighbor smoothing via `samp_neighbors`
    - adding noise to the expression data itself (via noise_trials)

    :param adata: anndata.AnnData containing the cells to be scored
    :param gene_set_file: the gene set file with list of gene sets, gmt, one per line
    :param samp_neighbors: number of neighbors to sample
    :param noise_trials: number of noisy samples to create, integer
    :param method_params: specific params for each method.

    :returns: adata with gene set scores in .obs
    """

    if error_checking(adata, samp_neighbors) == 'ERROR':
        return()

    gs_obj = genesets(gene_set_file)  # for scoring one set, make unit list, and proceed.

    # NOTE: this is cells x genes
    smoothed_matrix = nn_smoothing(adata.X, adata, 'connectivity', samp_neighbors)
    # for easier handling with gene names
    smoothed_adata = AnnData(smoothed_matrix, obs=adata.obs, var=adata.var)

    all_scores = _score_all_cells_all_sets(smoothed_adata=smoothed_adata, gene_set_obj=gs_obj,
                                           score_method=score_method, method_params=method_params,
                                           noise_trials=noise_trials)

    for gs in gs_obj.set_list:
        gs_name = gs.name
        gs_scores = [x[gs_name]['score'] for x in all_scores]  # for each cell, pull out gs_score gs_name
        adata.obs[gs_name] = gs_scores

    return(adata)


def _get_cell_data(
        smoothed_adata: anndata.AnnData,
        cell_ix: int,
        noise_trials: int,
        method_params: dict
) -> pd.DataFrame:
    """
    the processed expression data for each cell

    :param smoothed_adata: anndata.AnnData containing the cells to be scored
    :param cell_ix: index of the cell in adata
    :param noise_trials: number of noisy samples to create, integer
    :param method_params: specific params for each method.

    :return: a data frame
    """
    # for each cell, rank the expression
    gene_mat = smoothed_adata.X[cell_ix]
    # then we subset it to only the genes with counts
    _, gdx, _ = sparse.find(gene_mat)

    if gene_mat.ndim == 2:
        df = pd.DataFrame(gene_mat[:, gdx].A.flatten(), index=smoothed_adata.var.index[gdx])  ## ????
    else:
        df = pd.DataFrame(gene_mat[gdx],
                          index=smoothed_adata.var.index[gdx])  ## not sure why it's coming off as an array
    df.columns = ['counts']

    if ('normalization' in method_params) and (method_params['normalization'] == 'average') and (noise_trials > 0):
        # add some noise to gene counts.. create a n numbers of examples
        raise ValueError('not implemented')
        #df_noise = si.add_noise(df, noise_trials, 0.01, 0.99)  ## slow part .. fixed
    else:
        df_noise = df

    # for right now always ranking genes up #
    df_noise['uprank'] = df_noise.iloc[:,0].rank(method='min', ascending=True)  # up or down
    df_noise['dnrank'] = np.max(df['uprank']) - df['uprank']

    return(df_noise)


def _score_all_cells_all_sets(
        smoothed_adata: anndata.AnnData,
        gene_set_obj: genesets,
        score_method: str,
        method_params: dict,
        noise_trials: int
        ) -> list:
    """
    Process cells and score each with a list of gene sets and a method

    :param smoothed_adata: anndata.AnnData containing the cells to be scored
    :param gene_set_obj: list of geneset objects
    :param score_method: what method we'll be calling
    :param method_params: specific params for each method.
    :param noise_trials: number of noisy samples to create, integer

    :return: list of list of gene set score dictionaries
    """
    results_list = []  # one per cell
    for cell_ix in tqdm.trange(smoothed_adata.shape[0]):  # for each cell ID
        results = dict()                                  #   we will have one score per cell
        df_cell = _get_cell_data(smoothed_adata, cell_ix, noise_trials, method_params)  # process the cell's data
        for gs_i in gene_set_obj.set_list:                #   for each gene set
            res0 = scorefun(gs_i, df_cell, score_method, method_params, smoothed_adata.obs.index[cell_ix])
            results[gs_i.name] = res0
        results_list.append( results )
    return(results_list)
