import anndata
from gssnng.score_cells import _proc_data
from gssnng.util import error_checking
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
    nearest neighbor smoothing of the expression matrix

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

    :returns: a list of adatas with smoothed data
    """

    return_data = 1
    noise_trials = 0  ### not used currently
    samp_neighbors = None ### also not used
    just_smoothing=1

    error_checking(adata, samp_neighbors, recompute_neighbors,
                   None, None, None, method_params, just_smoothing)

    if method_params == None:
        method_params = dict()

    # score each cell with the list of gene sets
    data_list = _proc_data(adata, None, groupby, smooth_mode, recompute_neighbors,
                                  None, method_params, samp_neighbors,
                                  noise_trials, None, cores, return_data)

    print("**done**")
    return(data_list)
