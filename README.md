# gssnng

https://colab.research.google.com/github/jakevdp/PythonDataScienceHandbook/blob/master/notebooks/Index.ipynb

Gene Set Scoring on the Nearest Neighbor Graph (gssnng) for Single Cell RNA-seq (scRNA-seq)

Works with AnnData objects stored as h5ad files. Takes values from adata.X.

Scoring functions:
```
    singscore:  mean(ranks) / n, where n is length of gene set
    
    robust_std:   med (x-med / mad), median of robust standardized ranks.
    
    mean_z:  mean( (x - mean)/stddv ), average z score.
    
    rank_biased_overlap: weighted average of agreement across depths, repeated intersection of set with ranked order.
    
    average_score:  sum(ranks) / n, average ranks.     
    
    median_score:  median score: med()
    
    summed_up: just sum up the ranks or counts.
```

# Instructions (see notebooks)

1. Read in an AnnData object using scanpy (an h5ad file).

2. Get gene sets formatted as a .gmt file.

3. Score cells, each gene set will show up as a column in adata.obs.

```
from gssnng import score_cells

q = sc.read_h5ad('gssnng/gssnng/test/data/pbmc3k_processed.h5ad')

sc.pp.neighbors(q, n_neighbors=32)

score_cells.with_gene_sets(
        adata=q,
        gene_set_file="gssnng/gssnng/test/data/cibersort_lm22.gmt",
        score_method="robust_std",
        method_params=dict(),
        samp_neighbors=29,
        noise_trials=0,
        keys_added=None
    )
    

sc.pl.umap(q, color=['louvain','T.cells.CD8'], wspace=0.35)
```

## References

rank biased overlap:  https://arxiv.org/pdf/1408.3587.pdf


