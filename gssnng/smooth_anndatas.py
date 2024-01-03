import anndata
from gssnng.score_cells import _proc_data
#from gssnng.util import error_checking
from typing import Union

def smooth_anndata(
        adata: anndata.AnnData,
        groupby: Union[str, list, dict],
        smooth_mode: str,
        recompute_neighbors: int,
        method_params: dict,
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
    :param groupby
        either a column label in adata.obs, and all categories taken, or a dict specifies one group.
    :param smooth_mode
        `adjacency` or `connectivity`, which representation of the neighborhood graph to use.
        `adjacency` weights all neighbors equally, `connectivity` weights close neighbors more
    :param recompute_neighbors
        should neighbors be recomputed within each group, 0 for no, >0 for yes and specifies N
    :param method_params
        specific params for each method.
    :param cores
        number of parallel processes to work through groupby groups

    :returns: adata with gene set scores in .obs
    """

    return_data = 1
    noise_trials = 0  ### not used currently
    samp_neighbors = None

    #error_checking2(adata, recompute_neighbors, method_params)  # UPDATE

    if method_params == None:
        method_params = dict()

    # score each cell with the list of gene sets
    data_list = _proc_data(adata, None, groupby, smooth_mode, recompute_neighbors,
                                  None, method_params, None,
                                  noise_trials, None, cores, return_data)

    print("**done**")
    return(data_list)
