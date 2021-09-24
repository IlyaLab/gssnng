

import scanpy as sc

from gssnng.score_one_set import score_cells
from gssnng.util import read_gene_sets


def test_score_all_fun(adata, genesets):
    key1 = list(genesets.keys())[0]
    print('scoring: ' + key1)
    geneset1 = genesets[key1]
    res0 = score_cells(
        adata=adata,
        gene_set_up=geneset1,
        gene_set_down=None,
        score_method='summed_up',
        key_added=key1,
        samp_neighbors=8,
        noise_trials=0,
        mode='average')
    return(res0)


def test_score_all():
    q = sc.read_h5ad('data/pbmc3k_processed.h5ad')
    gs = read_gene_sets('data/h.all.v7.2.symbols.gmt')  #returns a dict
    q2 = q[q.obs['louvain'] == 'NK cells']
    print(q2)
    score_list = test_score_all_fun(q2, gs)
    print('*************')
    print('post function')
    print(q2.obs[0:5])


test_score_all()
print('test score_all done')

