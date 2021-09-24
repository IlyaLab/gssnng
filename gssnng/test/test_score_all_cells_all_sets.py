

import scanpy as sc
from gssnng.score_all_sets import score_cells_all_sets


def test_score_all_sets_fun(adata, genesets):
    res0 = score_cells_all_sets(
        adata=adata,
        gene_set_file=genesets,
        score_method='summed_up',
        set_direction='up',
        key_added='GeneSetNames', # only option right now
        samp_neighbors=8,
        noise_trials=0,
        mode='average')
    return(res0)


def test_score_all_sets():
    q = sc.read_h5ad('data/pbmc3k_processed.h5ad')
    q2 = q[q.obs['louvain'] == 'NK cells']
    gs = 'data/h.all.v7.2.symbols.gmt'
    score_list = test_score_all_sets_fun(q2, gs)
    print('*************')
    print('post function')
    print(q2.obs[0:5])


test_score_all_sets()
print('test score_all_sets done')

