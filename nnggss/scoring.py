"""
MS version
"""
import numpy as np
from scipy import sparse
import pandas as pd
import nnggss.util as si
import tqdm
import statsmodels.robust.scale
from anndata import AnnData
from nnggss.smoothing import nn_smoothing


def sc_score(
        adata,
        noise_trials,
        samp_neighbors,
        gene_set,
        mode='average',
        ):

    """
    gene set scoring with nearest neighbor smoothing of the expression matrix

    Improves upon the usual scsing scoring by:
    - smoothing the data matrix
        - adding noise to the nearest neighbor smoothing via `samp_neighbors`
    - adding noise to the expression data itself (via noise_trials)

    :param adata: anndata.AnnData containing the cells to be scored
    :param noise_trials: number of noisy samples to create, integer
    :param samp_neighbors: number of neighbors to sample
    :param gene_set: the gene set of interest
    :param mode: average or theoretical normalization of scores

    :returns: np.array of scores, one per cell in adata
    """

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
    # # then we subset it to only the genes with counts
    # df = pd.DataFrame(gene_mat.T.A, index=adata.var.index, columns=adata.obs.index)
    # gdx = df.sum(1) == 0
    # df = df.iloc[gdx,:]
    # df.columns = ['gene_counts']
    all_scores = _score_one_by_one(gene_set, smoothed_adata, noise_trials, mode='average')
    return all_scores


def sc_score_ms(
        adata,
        noise_trials,
        samp_neighbors,
        gene_set_up,
        gene_set_down,
        mode='average'
        ):

    """
    gene set scoring with nearest neighbor smoothing of the expression matrix

    Improves upon the usual scsing scoring by:
    - smoothing the data matrix
        - adding noise to the nearest neighbor smoothing via `samp_neighbors`
    - adding noise to the expression data itself (via noise_trials)

    :param adata: anndata.AnnData containing the cells to be scored
    :param noise_trials: number of noisy samples to create, integer
    :param samp_neighbors: number of neighbors to sample
    :param gene_set_up: the gene set of interest for up expressed genes
    :param gene_set_down: the gene set of interest for down expressed genes
    :param mode: average or theoretical normalization of scores

    :returns: np.array of scores, one per cell in adata
    """

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
    all_scores = _score_all_at_once(gene_set_up=gene_set_up, gene_set_down=gene_set_down,
                                    smoothed_adata=smoothed_adata, noise_trials=noise_trials, mode=mode)
    return all_scores



def _score_one_by_one(gene_set, smoothed_adata, noise_trials, mode='average', ):
    scores_across_cells = []
    for cell_ix in tqdm.trange(smoothed_adata.shape[0]):
        gene_mat = smoothed_adata.X[cell_ix]
        # then we subset it to only the genes with counts
        _, gdx, _ = sparse.find(gene_mat)
        df = pd.DataFrame(gene_mat[:, gdx].A.flatten(), index=smoothed_adata.var.index[gdx])
        df.columns = ['gene_counts']

        # FROM HERE: parallel computing of scores across gene sets.
        # or parallelize on each column of dataframe

        if mode == 'average' and noise_trials > 0:
            # add some noise to gene counts.. create a n numbers of examples
            df_noise = si.add_noise(df, noise_trials, 0.01, 0.99)  # slow part .. fixed
        else:
            df_noise = df
        # score the neighborhoods
        s = si.score(up_gene=gene_set, sample=df_noise, norm_method='standard', full_data=False) # standard workin gbetter here than theoretical

        avg_score = s.mean()
        scores_across_cells.append(avg_score)

    return np.array(scores_across_cells).flatten()


def _ms_sing(geneset: list, x: pd.Series, norm_method: str, rankup: bool) -> dict:
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
    up_sort = x.rank(method='min', ascending=rankup)  #
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
                               library_len=len(x.index),
                               score_list=su,
                               score=score_up,
                               sig_len=sig_len_up)
    norm_up = norm_up - 0.5
    mad_up = statsmodels.robust.scale.mad(su)
    total_score = norm_up
    return dict(total_score=total_score, mad_up=mad_up)


def _score_all_at_once(gene_set_up=False, gene_set_down=False, smoothed_adata=None, noise_trials=0, mode='average'):
    """
    not really, but at least call `si.score` only once
    """
    results = []
    for cell_ix in tqdm.trange(smoothed_adata.shape[0]):
        gene_mat = smoothed_adata.X[cell_ix]
        # then we subset it to only the genes with counts
        _, gdx, _ = sparse.find(gene_mat)
        # TODO we could do a dict instead of the df, that would be faster in _mssing too
        df = pd.DataFrame(gene_mat[:, gdx].A.flatten(), index=smoothed_adata.var.index[gdx])
        df.columns = ['gene_counts']

        if mode == 'average' and noise_trials > 0:
            # add some noise to gene counts.. create a n numbers of examples
            raise ValueError('not implemented')
            df_noise = si.add_noise(df, noise_trials, 0.01, 0.99) ## slow part .. fixed
        else:
            df_noise = df

        # Handle the cases of up vs down gene sets #
        if gene_set_up and (not gene_set_down):
            s = _ms_sing(gene_set_up, df_noise['gene_counts'], norm_method='standard', rankup=True)
        elif (not gene_set_up) and gene_set_down:
            s = _ms_sing(gene_set_down, df_noise['gene_counts'], norm_method='standard', rankup=False)
            s['mad_down'] = s.pop('mad_up')
        else: # both gene sets
            s_up = _ms_sing(gene_set_up, df_noise['gene_counts'], norm_method='standard', rankup=True)
            s_down = _ms_sing(gene_set_down, df_noise['gene_counts'], norm_method='standard', rankup=False)
            s = dict(total_score=(s_up['total_score']+s_down['total_score']),
                     mad_up=s_up['mad_up'],
                     mad_down=s_down['mad_up'],
                     up_score=s_up['total_score'],
                     dn_score=s_down['total_score']
                     )
        s['CB'] = smoothed_adata.obs.index[cell_ix]
        results.append(s)
    return results
