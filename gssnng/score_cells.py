import anndata
import numpy as np
from scipy import sparse
import pandas as pd
import scanpy as sc
from anndata import AnnData
from gssnng.smoothing import nn_smoothing
from gssnng.util import error_checking
from gssnng.score_funs import scorefun
from gssnng.gene_sets import genesets_from_gmt, Genesets, genesets_from_decoupler_model
from typing import Union
from multiprocessing import Pool


def run_gssnng(
    mat, net, source, target, weight, 
    groupby: Union[str, list, dict],
    smooth_mode: str,
    recompute_neighbors: int,
    score_method: str,
    method_params: dict,
    #samp_neighbors: int,
    ranked: bool,
    cores: int,
    verbose=False,
    use_raw=True
    ):
    
    """
    GSSNNG: gene set scoring (all gene sets in file) with nearest neighbor smoothing of the expression matrix

    Improved single cell scoring by:
    - smoothing the data matrix
        - adding noise to the nearest neighbor smoothing via `samp_neighbors`
    - adding noise to the expression data itself (via noise_trials)


    Parameters
    ----------
    mat : list, DataFrame or AnnData
        List of [features, matrix], dataframe (samples x features) or an AnnData
        instance.
    net : DataFrame
        Network in long format.
    source : str
        Column name in net with source nodes.
    target : str
        Column name in net with target nodes.
    weight : str
        Column name in net with weights.
    groupby : str
        Either a column label in adata.obs, and all categories taken, or a dict specifies one group.
    smooth_mode : ['adjacency','connectivity']
        Which representation of the neighborhood graph to use.
        `adjacency` weights all neighbors equally, `connectivity` weights close neighbors more
    recompute_neighbors : bool
        should neighbors be recomputed within each group, 0 for no, >0 for yes and specifies N.
    score_method : str
        Which scoring method to use.
    method_params: dict
        specific params for each method.
    ranked: int
        Whether the gene expression counts should be rank ordered
    cores: int
        number of parallel processes to work through groupby groups

    Returns
    -------
    if mat is an AnnData object:
        gssnng scores in `mat.obsm['gssnng_estimate']` 
    else:
        pd.DataFrame with the scores
    """

    # This whole function is basically a wrapper around the old gssnng-API
    # but working with decoupler-input and -output
    noise_trials = 0  ### !TODO! not used currently

    # our gene set data object list
    gs_obj = genesets_from_decoupler_model(net, source, target , weight)

    samp_neighbors = None
    error_checking(mat, samp_neighbors, recompute_neighbors,
                   gs_obj, score_method, ranked, method_params)

    if method_params == None:
        method_params = dict()

    # score each cell with the list of gene sets
    all_scores = _proc_data(mat, gs_obj, groupby, smooth_mode, recompute_neighbors,
                                  score_method, method_params, samp_neighbors,
                                  noise_trials, ranked, cores)

    # warning: the all_scores rows might have a diferent order!
    # make sure to resort them according to the mat.obs.index
    all_scores = all_scores.loc[mat.obs.index]

    # AnnData support
    if isinstance(mat, AnnData):
        # Update obsm AnnData object
        mat.obsm['gssnng_estimate'] = all_scores
    else:
        return all_scores


def with_gene_sets(
        adata: anndata.AnnData,
        gene_set_file: str,
        groupby: Union[str, list, dict],
        smooth_mode: str,
        recompute_neighbors: int,
        score_method: str,
        method_params: dict,
        #samp_neighbors: int,
        ranked: bool,
        cores: int
    ) -> anndata.AnnData:

    """
    gene set scoring (all gene sets in file) with nearest neighbor smoothing of the expression matrix

    Improved single cell scoring by:
    - smoothing the data matrix
        - adding noise to the nearest neighbor smoothing via `samp_neighbors`
    - adding noise to the expression data itself (via noise_trials)

    :param adata
        anndata.AnnData containing the cells to be scored
    :param gene_set_file
        the gene set file with list of gene sets, gmt, one per line
    :param groupby
        either a column label in adata.obs, and all categories taken, or a dict specifies one group.
    :param smooth_mode
        `adjacency` or `connectivity`, which representation of the neighborhood graph to use.
        `adjacency` weights all neighbors equally, `connectivity` weights close neighbors more
    :param recompute_neighbors
        should neighbors be recomputed within each group, 0 for no, >0 for yes and specifies N
    :param score_method
        which scoring method to use
    :param method_params
        specific params for each method.
    :param ranked
        whether the gene expression counts should be rank ordered
    :param cores
        number of parallel processes to work through groupby groups

    :returns: adata with gene set scores in .obs
    """

    noise_trials = 0  ### not used currently

    # our gene set data object list
    gs_obj = genesets_from_gmt(gene_set_file)

    samp_neighbors = None
    error_checking(adata, samp_neighbors, recompute_neighbors,
                   gs_obj, score_method, ranked, method_params)

    if method_params == None:
        method_params = dict()

    # score each cell with the list of gene sets
    all_scores = _proc_data(adata, gs_obj, groupby, smooth_mode, recompute_neighbors,
                                  score_method, method_params, samp_neighbors,
                                  noise_trials, ranked, cores)
    ## join in new results
    adata.obs = adata.obs.join(all_scores, how='left')

    print("**done**")
    return(adata)


def _smooth_out(adata, samp_neighbors, smooth_mode):
    """
    calculates a neighourhood-smoothed signal of adata.X.
    :param samp_neighbors: how many NN to subsample (hence introducing some noise)
    :param smooth_mode: 'connectivity', 'adjacency', see nn_smoothing()

    Returns a complete AnnData with smooeth adata.X
    """
    if smooth_mode not in ['connectivity', 'adjacency', 'off']:
        raise Exception("ERROR:  please use smooth mode: `adjacency`, `connectivity`, or `off`.")

    if smooth_mode == 'off':
        adata.obsm['X_smooth'] = sparse.csr_matrix(adata.X)
    else:
        smoothed_matrix = nn_smoothing(adata.X, adata, smooth_mode, samp_neighbors)
        # for easier handling with gene names
        # smoothed_adata = AnnData(smoothed_matrix, obs=adata.obs, var=adata.var)
        adata.obsm['X_smooth'] = smoothed_matrix
    return(adata)


def _build_data_list(
        adata,
        groupby,
        cats,
        recompute_neighbors,
        samp_neighbors,
        smooth_mode,
):
    """
    creates the smoothed adata.
    Smoothing proceeds within splits via 'groupby' (i.e. a subset of the true NN graph)
    """
    data_list = []  # list of dicts
    for ci in cats:
        # then for each group
        qi = adata[adata.obs['gssnng_groupby'] == ci]
        if recompute_neighbors > 0:
            print("RECOMPUTING KNN")
            sc.pp.neighbors(qi, n_neighbors=recompute_neighbors)
        # smooth and make an adata, and add to the list
        qi_smoothed = _smooth_out(qi, samp_neighbors, smooth_mode)
        params = (qi_smoothed, ci)
        data_list.append(params)
    # send off for scores
    return(data_list)


def _proc_data(
        adata: anndata.AnnData,
        gs_obj: Genesets,
        groupby: Union[str, list, dict],
        smooth_mode: str,
        recompute_neighbors: int,
        score_method: str,
        method_params: dict,
        samp_neighbors: int,
        noise_trials: int,
        ranked: bool,
        cores: int
                     ):
    """
    In many cases, the neighbors should be defined.  If you have mixed clinical endpoints,
    you might want scores for only the early points together, and scores for later points.
    Or T-cell neighbors should maybe only be other T-cells. etc. Or by leiden cluster label.
    By building out the list of groups, it can also be run in parallel.

    :param adata: anndata.AnnData containing the cells to be scored
    :param gs_obj: the gene sets class object
    :param groupby: either a column label in adata.obs, and all categories taken, or a dict specifies one group.
    :param recompute_neighbors: should the neighbors be recomputed within each group?  SLOW SLOW SLOW
    :param score_method: which scoring method to use
    :param method_params: specific params for each method.
    :param samp_neighbors: number of neighbors to sample
    :param noise_trials: number of noisy samples to create, integer
    :param ranked: whether the gene expression counts should be rank ordered
    :param cores: number of parallel processes to work through groupby groups

    :returns: scores in a dict for each cell in a list.
    """
    data_list = []  # list of dicts

    if groupby is None:  # take all cells
        cats = ['cat1']
        adata.obs['gssnng_groupby'] = 'cat1'  # altered to gssnng_groupby

    # here we group by one category
    elif isinstance(groupby, str):
        adata.obs['gssnng_groupby'] = adata.obs[groupby]  #
        cats = set(adata.obs[groupby])  # categories

    # here we groupby an intersection of categories, up to 4
    elif isinstance(groupby, list):
        # probably more readable:
        # glist = [adata.obs[g] for g in groupby ]
        # zipped_cols = zip(*glist)
        if len(groupby) == 2:
            zipped_cols = list(zip(adata.obs[groupby[0]], adata.obs[groupby[1]]))
        elif len(groupby) == 3:
            zipped_cols = list(zip(adata.obs[groupby[0]], adata.obs[groupby[1]], adata.obs[groupby[2]]))
        elif len(groupby) == 4:
            zipped_cols = list(zip(adata.obs[groupby[0]], adata.obs[groupby[1]], adata.obs[groupby[2]], adata.obs[groupby[3]]))
        else:
            raise Exception("ERROR: Number of groups must be 4 or less")
        adata.obs['gssnng_groupby'] = zipped_cols  #
        cats = set(zipped_cols)  # categories

    elif isinstance(groupby, dict):
        # then we want to only score
        return("ERROR: dict not implemented")

    data_list = _build_data_list(adata, groupby, cats, recompute_neighbors, samp_neighbors, smooth_mode)
    # then we can start scoring cells #

    # building up the argument list for the parallel call of _score_all_cells_all_sets
    arglist = []
    for smoothed_adata, groupname in data_list:
        arglist.append(
            (smoothed_adata, gs_obj, score_method, method_params, noise_trials, ranked, groupname)
        )
    # how about lambda x: _score_all_cells_all_sets(x, gs_obj, score_method, method_params, noise_trials, ranked)
    # not sure if that works with parallel, due to pickling. might have to be a proper function
    with Pool(processes=cores) as pool:
        res0 = pool.starmap_async(_score_all_cells_all_sets, arglist).get()
    # we column-bind the vectors into a pandas table
    res1 = pd.concat(res0, axis=0)
    return(res1)


def _get_cell_data(
        smoothed_adata: anndata.AnnData,
        cell_ix: int,
        noise_trials: int,
        method_params: dict,
        ranked: bool
) -> pd.DataFrame:
    """
    the processed expression data for each cell: Pulls out the cells expression vector,
    optinoally adding some noise

    :param smoothed_adata: anndata.AnnData containing the cells to be scored
    :param cell_ix: index of the cell in adata
    :param noise_trials: number of noisy samples to create, integer
    :param method_params: specific params for each method.

    :return: a data frame
    """
    # for each cell, rank the expression
    #gene_mat = smoothed_adata.X[cell_ix]
    gene_mat = (smoothed_adata.obsm['X_smooth'])[cell_ix]
    # then we subset it to only the genes with counts
    _, gdx, _ = sparse.find(gene_mat)

    if gene_mat.ndim == 2:
        df = pd.DataFrame(gene_mat[:, gdx].A.flatten(),
                          index=smoothed_adata.var.index[gdx])  ## ????
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
        smoothed_adata: AnnData,
        gene_set_obj: Genesets,
        score_method: str,
        method_params: dict,
        noise_trials: int,
        ranked: bool,
        group_name: str
        ) -> pd.DataFrame:
    """
    Process cells and score each with a list of gene sets and a method

    :param smoothed_adata: anndata.AnnData containing the cells to be scored
    :param gene_set_obj: the gene sets class object
    :param score_method: which scoring method to use
    :param method_params: specific params for each method.
    :param noise_trials: number of noisy samples to create, integer
    :param ranked: whether the gene expression counts should be rank ordered
    :param group_name: group of cells currently being processed

    :return: list of list of gene set score dictionaries
    """
    print("running " + str(group_name))

    results_list = []
    barcodes = []
    for cell_ix in range(smoothed_adata.shape[0]):  # tqdm.trange(smoothed_adata.shape[0]):  # for each cell ID
        results = dict()
        df_cell = _get_cell_data(smoothed_adata, cell_ix, noise_trials, method_params, ranked)
        barcodes.append(smoothed_adata.obs.index[cell_ix])
        for gs_i in gene_set_obj.set_list:  # for each gene set
            results[gs_i.name] = scorefun(gs_i, df_cell, score_method, method_params, ranked)
        results_list.append(results)
    results_df = pd.DataFrame(results_list, index=barcodes)
    return(results_df)
