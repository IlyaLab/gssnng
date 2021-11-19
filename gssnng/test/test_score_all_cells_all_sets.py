

import scanpy as sc
from gssnng.score_all_sets import score_cells_all_sets


def test_score_all_sets_fun(adata, genesets):
    res0 = score_cells_all_sets(
        adata=adata,
        gene_set_file=genesets,
        score_method='robust_std',
        set_direction='up',
        key_added='GeneSetNames', # only option right now
        samp_neighbors=27,
        noise_trials=0,
        mode='average')
    return(res0)


def test_score_all_sets():
    q = sc.read_h5ad('data/pbmc3k_processed.h5ad')
    q2 = q[q.obs.louvain == 'CD4 T cells']
    print(q2)
    gs = 'data/cibersort_lm22.gmt'  #'data/gene_set_test.gmt'
    print("computing knn...")
    sc.pp.neighbors(q2, n_neighbors=32)
    print('scoring...')
    score_list = test_score_all_sets_fun(q2, gs)
    print('******DONE*******')
    print(q2.obs.head())
    print(q2.obs.columns)
    #q.write_h5ad('data/pbmc3k_lm22_scores.h5ad')

test_score_all_sets()
print('test score_all_sets done')

