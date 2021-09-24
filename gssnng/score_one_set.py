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


def score_cells(
        adata,
        gene_set_up,
        gene_set_down,
        score_method,
        key_added,
        samp_neighbors,
        noise_trials,
        mode='average'
):

    """
    gene set scoring (One gene set) with nearest neighbor smoothing of the expression matrix

    Improved single cell scoring by:
    - smoothing the data matrix
        - adding noise to the nearest neighbor smoothing via `samp_neighbors`
    - adding noise to the expression data itself (via noise_trials)

    :param adata: anndata.AnnData containing the cells to be scored
    :param gene_set_up: the gene set of interest for up expressed genes
    :param gene_set_down: the gene set of interest for down expressed genes
    :param key_added: name given to the new entry of adata.obs['key_added']
    :param samp_neighbors: number of neighbors to sample
    :param noise_trials: number of noisy samples to create, integer
    :param mode: average or theoretical normalization of scores

    :returns: np.array of scores, one per cell in adata
    """
    if error_checking(adata, samp_neighbors) == 'ERROR':
        return()

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
    all_scores = _score_all_cells_at_once(gene_set_up=gene_set_up, gene_set_down=gene_set_down,
                                          smoothed_adata=smoothed_adata, noise_trials=noise_trials,
                                          mode=mode, score_method=score_method)

    adata.obs[key_added] = [x['score'] for x in all_scores]
    return(all_scores)



def _score_all_cells_at_once(gene_set_up=None, gene_set_down=None,
                             smoothed_adata=None, noise_trials=0,
                             mode='average', score_method='singscore'):
    """
    not really, but at least call `si.score` only once
    """
    results = []
    for cell_ix in tqdm.trange(smoothed_adata.shape[0]):
        gene_mat = smoothed_adata.X[cell_ix]
        # then we subset it to only the genes with counts
        _, gdx, _ = sparse.find(gene_mat)
        # TODO we could do a dict instead of the df, that would be faster in _mssing too
        if gene_mat.ndim == 2:
            df = pd.DataFrame(gene_mat[:, gdx].A.flatten(), index=smoothed_adata.var.index[gdx]) ## ????
        else:
            df = pd.DataFrame(gene_mat[gdx], index=smoothed_adata.var.index[gdx]) ## not sure why it's coming off as an array
        df.columns = ['gene_counts']

        if mode == 'average' and noise_trials > 0:
            # add some noise to gene counts.. create a n numbers of examples
            raise ValueError('not implemented')
            df_noise = si.add_noise(df, noise_trials, 0.01, 0.99) ## slow part .. fixed
        else:
            df_noise = df

        # Handle the cases of up vs down gene sets #
        if (gene_set_up != None) and (gene_set_down == None):
            s = _ms_sing(gene_set_up, df_noise['gene_counts'], score_method, norm_method='standard', rankup=True, dorank=True)
        elif (gene_set_up == None) and (gene_set_down != None):
            s = _ms_sing(gene_set_down, df_noise['gene_counts'], score_method, norm_method='standard', rankup=False, dorank=True)
            s['mad_down'] = s.pop('mad_up')
        else: # both gene sets
            s_up = _ms_sing(gene_set_up, df_noise['gene_counts'], score_method, norm_method='standard', rankup=True, dorank=True)
            s_down = _ms_sing(gene_set_down, df_noise['gene_counts'], score_method, norm_method='standard', rankup=False, dorank=True)
            s = dict(total_score=(s_up['total_score']+s_down['total_score']),
                     mad_up=s_up['mad_up'],
                     mad_down=s_down['mad_up'],
                     up_score=s_up['total_score'],
                     dn_score=s_down['total_score']
                     )
        s['CB'] = smoothed_adata.obs.index[cell_ix]
        results.append(s)
    return results

