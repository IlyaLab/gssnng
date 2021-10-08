

import scanpy as sc
from gssnng.score_all_sets import score_cells_all_sets


def test_score_all_sets_fun(adata, genesets):
    res0 = score_cells_all_sets(
        adata=adata,
        gene_set_file=genesets,
        score_method='rank_biased_overlap',
        set_direction='up',
        key_added='GeneSetNames', # only option right now
        samp_neighbors=27,
        noise_trials=0,
        mode='average',
        rbo_depth=33)
    return(res0)


def test_score_all_sets():
    q = sc.read_h5ad('data/pbmc3k_processed.h5ad')
    gs = 'data/celltypes_bone_and_blood_msigdb.gmt'
    print("computing knn...")
    sc.pp.neighbors(q, n_neighbors=32)
    print('scoring...')
    score_list = test_score_all_sets_fun(q, gs)
    print('******DONE*******')

test_score_all_sets()
print('test score_all_sets done')

