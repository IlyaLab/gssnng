import anndata
import numpy as np
from scipy import sparse
import pandas as pd
import scanpy as sc
import tqdm
from anndata import AnnData
from gssnng.smoothing import nn_smoothing
from gssnng.util import error_checking
from gssnng.score_funs import scorefun
from gssnng.genesets import genesets
from typing import Union
from multiprocessing import Pool

def with_gene_sets(
        adata: anndata.AnnData,
        gene_set_file: str,
        groupby: Union[str,dict],
        recompute_neighbors: int,
        score_method: str,
        method_params: dict,
        samp_neighbors: int,
        noise_trials: int,
        ranked: bool,
        threads: int
    ) -> anndata.AnnData:

    """
    gene set scoring (all gene sets in file) with nearest neighbor smoothing of the expression matrix

    Improved single cell scoring by:
    - smoothing the data matrix
        - adding noise to the nearest neighbor smoothing via `samp_neighbors`
    - adding noise to the expression data itself (via noise_trials)

    :param adata: anndata.AnnData containing the cells to be scored
    :param gene_set_file: the gene set file with list of gene sets, gmt, one per line
    :param groupby: either a column label in adata.obs, and all categories taken, or a dict specifies one group.
    :param recompute_neighbors: should neighbors be recomputed within each group, 0 for no, >0 for yes and specifies N
    :param score_method: which scoring method to use
    :param method_params: specific params for each method.
    :param samp_neighbors: number of neighbors to sample
    :param noise_trials: number of noisy samples to create, integer
    :param ranked: whether the gene expression counts should be rank ordered
    :param threads: number of parallel processes to work through groupby groups

    :returns: adata with gene set scores in .obs
    """

    if error_checking(adata, samp_neighbors) == 'ERROR':
        return()

    gs_obj = genesets(gene_set_file)  # for scoring one set, make unit list, and proceed.

    all_scores = _proc_data(adata, gs_obj, groupby, recompute_neighbors,
                                  score_method, method_params, samp_neighbors,
                                  noise_trials, ranked, threads)

    for gs in gs_obj.set_list:
        gs_name = gs.name
        gs_scores = [x[gs_name]['score'] for x in all_scores]  # for each cell, pull out gs_score gs_name
        adata.obs[gs_name] = gs_scores

    return(adata)


def _smooth_out(adata, samp_neighbors):
    # NOTE: this is cells x genes
    smoothed_matrix = nn_smoothing(adata.X, adata, 'connectivity', samp_neighbors)
    # for easier handling with gene names
    smoothed_adata = AnnData(smoothed_matrix, obs=adata.obs, var=adata.var)
    return(smoothed_adata)


def _proc_data(
        adata: anndata.AnnData,
        gs_obj: genesets,
        groupby: Union[str,dict],
        recompute_neighbors: bool,
        score_method: str,
        method_params: dict,
        samp_neighbors: int,
        noise_trials: int,
        ranked: bool,
        threads: int
                     ):
    """
    In many cases, the neighbors should be defined.  If you have mixed clinical endpoints,
    you might want scores for only the early points together, and scores for later points.
    Or T-cell neighbors should maybe only be other T-cells. etc. Or by leiden cluster label.
    By building out the list of groups, it can also be run in parallel.

    :param adata: anndata.AnnData containing the cells to be scored
    :param gs_obs: the gene sets class object
    :param groupby: either a column label in adata.obs, and all categories taken, or a dict specifies one group.
    :param recompute_neighbors: should the neighbors be recomputed within each group?
    :param score_method: which scoring method to use
    :param method_params: specific params for each method.
    :param samp_neighbors: number of neighbors to sample
    :param noise_trials: number of noisy samples to create, integer
    :param ranked: whether the gene expression counts should be rank ordered
    :param threads: number of parallel processes to work through groupby groups

    :returns: scores in a dict for each cell in a list.
    """
    data_list = []  # list of dicts

    if groupby is None:  # take all cells
        if recompute_neighbors > 0:
            sc.pp.neighbors(adata, n_neighbors=recompute_neighbors)
        smoothed_adata =  _smooth_out(adata, samp_neighbors)  # list of 1 item
        params = {'smoothed_adata': smoothed_adata, 'gene_set_obj': gs_obj,
                  'score_method': score_method, 'method_params': method_params,
                  'noise_trials': noise_trials, 'ranked': ranked}
        data_list.append(params)

    elif type(groupby) == type("string"):
        cats = set(adata.obs[groupby])
        for ci in cats:
            ### then for each group
            qi = adata[adata.obs[groupby] == ci]
            if recompute_neighbors > 0:
                sc.pp.neighbors(qi, n_neighbors=recompute_neighbors)
            ### smooth and make an adata, and add to the list
            qi_smoothed = _smooth_out(qi, samp_neighbors)
            params = { 'smoothed_adata':qi_smoothed, 'gene_set_obj':gs_obj,
                        'score_method':score_method, 'method_params':method_params,
                        'noise_trials':noise_trials, 'ranked':ranked, 'group_name':ci }
            data_list.append( params )
            ### send off for scores

    with Pool(processes=threads) as pool:
        res0 = pool.map(_score_all_cells_all_sets, data_list)
    return(res0)


def _get_cell_data(
        smoothed_adata: anndata.AnnData,
        cell_ix: int,
        noise_trials: int,
        method_params: dict,
        ranked: bool
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

    if ranked:
        # for right now always ranking genes up #
        df_noise['uprank'] = df_noise.iloc[:,0].rank(method='min', ascending=True)  # up or down
        df_noise['dnrank'] = np.max(df['uprank']) - df['uprank']

    return(df_noise)


def _score_all_cells_all_sets(
        params
        ) -> list:
    """
    Process cells and score each with a list of gene sets and a method

    :param params: list of parameters from _build_data_list

    :return: list of list of gene set score dictionaries
    """
    smoothed_adata = params['smoothed_adata']
    gene_set_obj = params['gene_set_obj']
    score_method = params['score_method']
    method_params = params['method_params']
    noise_trials = params['noise_trials']
    ranked = params['ranked']
    group_name = params['group_name']

    print("running " + group_name)

    results_list = []  # one entry per cell
    for cell_ix in tqdm.trange(smoothed_adata.shape[0]):  # for each cell ID
        results = dict()                                  #   we will have one score per cell
        df_cell = _get_cell_data(smoothed_adata, cell_ix, noise_trials, method_params, ranked)  # process the cell's data
        for gs_i in gene_set_obj.set_list:                #   for each gene set
            res0 = scorefun(gs_i, df_cell, score_method, method_params, smoothed_adata.obs.index[cell_ix], ranked)
            results[gs_i.name] = res0
        results_list.append( results )
    return(results_list)
