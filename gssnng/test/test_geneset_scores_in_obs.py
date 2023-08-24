
# here we test the that the number of scores in the obs table is correct.
# see the gene set ingestion test
# There are 22 signatures in the LM22, 1 gssnng column, and 4 others for 27 columns

# doubles:
# B.cells.naive.up and .dn  ## in order and next to one another
# T.cells.CD4.memory.resting_dn  ### are out of order
# T_cells_gamma_delta_up ### the two are separated
# Dendritic.cells.activated.up ### doesn't have a pair.
# Macrophages.M1.dn ### no match but is down

### MSIGDB is labeled as UP and DN ###


# rootdir: /home/runner/work/gssnng/gssnng

import scanpy as sc
from gssnng import score_cells
def test_get_number_of_columns_in_obs(adata, genesets):
    res0 = score_cells.with_gene_sets(adata=adata,
                          gene_set_file=genesets,
                          groupby='louvain',
                          smooth_mode='connectivity',
                          recompute_neighbors=32,
                          score_method='median_score',
                          method_params={},
                          ranked=False,
                          cores=4)
    return (res0.obs.shape[1])


def test_number_of_genesets():
    q = sc.datasets.pbmc3k_processed()
    gs = 'gssnng/test/data/cibersort_lm22.gmt'  # 'data/gene_set_test.gmt' #'data/cibersort_lm22.gmt'  #
    assert test_get_number_of_columns_in_obs(q,gs) == 27  ## The up and dn sets should be combined into one.

