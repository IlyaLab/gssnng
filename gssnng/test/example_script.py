from gssnng import score_cells
import scanpy as sc
import time

print("reading data")
q = sc.read_h5ad('gssnng/gssnng/test/data/pbmc3k_processed.h5ad')

print("computing neighborhood")
sc.pp.neighbors(q, n_neighbors=32)

t0 = time.time()
print('start time: ' + str(t0))

print("scoring cells")
score_cells.with_gene_sets(
        adata=q,
        gene_set_file="gssnng/gssnng/test/data/cibersort_lm22.gmt",
        groupby='louvain',
        recompute_neighbors=0,
        score_method="singscore",
        method_params={'normalization':'theoretical'},
        samp_neighbors=27,
        ranked=True,
        cores=6
    )

t1 = time.time()
print('end time: ' + str(t1))
print('TOTAL TIME: ' + str(t1-t0))
print("done")
