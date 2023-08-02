import numpy as np
from scipy import sparse
import logging


NN_DISTANCE_KEY = 'distances'  # scanpy names in .obsp
NN_CONN_KEY = 'connectivities'

# TODO test: should always sum to 1
# multiplying should leave a "one-vector" still sum to one


def get_smoothing_matrix(adata, mode, add_diag):
    """
    using the nearest neighbor graph in adata.obsp, calculate the smoothing
    matrix S such that S @ X smoothes the signal X over neighbors

    :param adata: anndata.AnnData gene expression matrix to smooth
    :param mode: `adjacency` or `connectivity`, which representation of the neighborhood graph to use.
        `adjacency` weights all neighbors equally, `connectivity` weights close neighbors more
    :param add_diag: should the datapoint itself be consider in the smoothing. This should ==True in almost all cases!

    :return: scipy.sparse sparse matrix with the smoothed expression
    """

    if mode == 'adjacency':
        A = (adata.obsp[NN_DISTANCE_KEY] > 0).astype(int)

        # add the diagnoal, ie. the datapoint itself should be represented
        # in the smoothing!
        assert np.all(A.diagonal() == 0), "diagonal of distance matrix not 0!!"
        if add_diag:
            A = A + sparse.diags(np.ones(A.shape[0]))

        # normalize to sum=1 per  row
        row_scaler = 1 / A.sum(axis=1).A.flatten()

        # multiplying with a diag matrix d on the left multiplies each row i by d_i
        normA = sparse.diags(row_scaler) @ A
        return normA

    # actually works exactly the same; could just switch out the obsp key
    elif mode == 'connectivity':
        A = adata.obsp[NN_CONN_KEY] ## neighbor graph
        # add the diagnoal, ie. the datapoint itself should be represented
        # in the smoothing! Note that the max connectivity == 1
        assert np.all(A.diagonal() == 0), "diagonal of connectivity matrix not 0!!"
        if add_diag:
            A = A + sparse.diags(np.ones(A.shape[0]))
        # normalize to sum=1 per  row
        row_scaler = 1 / A.sum(axis=1).A.flatten()
        normA = sparse.diags(row_scaler) @ A
        return normA
    else:
        raise ValueError(f'unknown mode {mode}')


def random_mask_a_nn_matrix(X, nn_to_keep):
    """
    for a nearest neighbour matrix (i.e. each row has N entries)
    subsample the neighours (setting some entries per row to 0)

    this is pretty slow, maybe there's a better way...

    :param X: nearest neighbor matrix (sparse)
    :nn_to_keep: out of the N nearest neighbors how many to keep (the others will be set to zero)

    :return: sparse.csr_matrix with randomly subsampled nearest neighbors
    """
    # TODO if X is in csr format, we can quickly subsample!
    # just set some elements of data[intptr[row]:intptr[row+1]] to zero
    # and eliminate zeros

    newrows = []
    newcols = []
    newvals = []
    for i in range(X.shape[0]):
        _, col_ix, vals = sparse.find(X[i])
        col_ix_len = len(col_ix)
        if nn_to_keep > col_ix_len:
            samp_n = col_ix_len
        else:
            samp_n = nn_to_keep  # sometimes there are not exactly the expected number of values.
        rand_ix = np.random.choice(col_ix_len, size=samp_n, replace=False)
        newrows.extend([i]*samp_n)
        newcols.extend(col_ix[rand_ix])
        newvals.extend(vals[rand_ix])

    # TODO choose matrix format as X
    return sparse.csr_matrix((newvals, (newrows, newcols)))


def nn_smoothing(X, adata, mode, samp_neighbors, add_diag=True):
    """
    smooth the expression matrix X (cells x genes) using the neighborhood
    graph stored in adata.

    :param X: data to smooth. (cell x feature) matrix
    :param adata: sc.AnnData, containing the neighborhood graph
    :param samp_neighbors: consider all neighbours or a subsample of size samp_neighbors
    :param add_diag: consider the datapoint itself in smoothing. Should be True in almost all cases

    :return: scipy.sparse matrix with the smoothed signals
    """
    assert X.shape[0] == adata.shape[0], "cell number mismatch"
    logging.info("creating smoothing matrix")
    smoothing_mat = get_smoothing_matrix(adata, mode, add_diag)

    if (samp_neighbors is not None) and (samp_neighbors > 0):
        # randomly set some edges/neighbours to zero
        # experimental and pretty slow!
        logging.info("creating random mask")
        smoothing_mat = random_mask_a_nn_matrix(smoothing_mat, samp_neighbors)

    smooth_signals = smoothing_mat @ X

    return(smooth_signals)
