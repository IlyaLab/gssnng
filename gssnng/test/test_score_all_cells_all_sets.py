

import scanpy as sc
from gssnng.score_cells import with_gene_sets
import time

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
    t0 = time.time()
    print('start time: ' + str(t0))
    score_list = test_score_all_sets_fun(q2, gs)
    print('******DONE*******')
    t1 = time.time()
    print('end time: ' + str(t1))
    print('TOTAL TIME: ' + str(t1-t0))
    print(q2.obs.head())
    #q.write_h5ad('data/pbmc3k_lm22_scores.h5ad')

test_score_all_sets()
print('test score_all_sets done')
