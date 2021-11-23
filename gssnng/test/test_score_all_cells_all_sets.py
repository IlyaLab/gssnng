

import scanpy as sc
from gssnng.score_cells import with_gene_sets


def test_score_all_sets_fun(adata, genesets):
    res0 = with_gene_sets(adata=adata, gene_set_file=genesets, score_method='mean_z', method_params=dict(),
                          samp_neighbors=27, noise_trials=0, ranked=True)
    return(res0)


def test_score_all_sets():
    q = sc.read_h5ad('data/pbmc3k_processed.h5ad')
    q2 = q[q.obs.louvain == 'CD4 T cells']
    gs = 'data/gene_set_test.gmt' #'data/cibersort_lm22.gmt'  #
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

