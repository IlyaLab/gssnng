from anndata import AnnData
import numpy as np
from nnggss.smoothing import get_smoothing_matrix
from scipy import sparse

def test_get_smoothing_matrix():
    """
    assert the the rows of the smoothing matrix are normalized to 1
    """
    ncells = 4
    ngenes = 3

    adjacency = np.array([
        [0, 0.1, 0, 0],
        [0, 0, 0, 0],  # weird edge case where one cell doesnt have any neighbors!
        [0, 0, 0, 1],
        [0, 0, 1, 0]])

    adata = AnnData(np.random.rand(ncells,ngenes))
    adata.obsp['distances'] = sparse.csr_matrix(adjacency)

    # adjacency has to ignore the weightings of the neighbors
    S = get_smoothing_matrix(adata, mode='adjacency', add_diag=False)
    np.testing.assert_allclose(S.A.sum(1), [1, 0, 1, 1])


    adata.obsp['connectivities'] = sparse.csr_matrix(adjacency)
    S = get_smoothing_matrix(adata, mode='connectivity', add_diag=False)
    np.testing.assert_allclose(S.A.sum(1), [1, 0, 1, 1])


    # if there's diagnoal elements the results are differnt
    S = get_smoothing_matrix(adata, mode='adjacency', add_diag=True)
    np.testing.assert_allclose(S.A.sum(1), [1, 1, 1, 1])


    adata.obsp['connectivities'] = sparse.csr_matrix(adjacency)
    S = get_smoothing_matrix(adata, mode='connectivity', add_diag=True)
    np.testing.assert_allclose(S.A.sum(1), [1, 1, 1, 1])

test_get_smoothing_matrix()
print('test done')