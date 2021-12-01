# gssnng

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Gibbsdavidl/gssnng/blob/main/notebooks/gssnng_quick_start.ipynb)


Gene Set Scoring on the Nearest Neighbor Graph (gssnng) for Single Cell RNA-seq (scRNA-seq)

Works with AnnData objects stored as h5ad files. Takes values from adata.X.

Scoring functions, works with ranked or unranked data (**"your mileage may vary"**):
```
    singscore:  mean(ranks) / n, where n is length of gene set
    
    robust_std:   med (x-med / mad), median of robust standardized values (ranks).
    
    mean_z:  mean( (x - mean)/stddv ), average z score. (recommend unranked)
    
    rank_biased_overlap: weighted average of agreement across depths, repeated intersection of set with ranked order.
    
    average_score:  sum(ranks) / n, average ranks.     
    
    median_score:  median score: med()
    
    summed_up: just sum up the ranks or counts.
```

## Installation

git clone https://github.com/Gibbsdavidl/gssnng
pip install -e gssnng

## Instructions (see notebooks for examples)

1. Read in an AnnData object using scanpy (an h5ad file).

2. Get gene sets formatted as a .gmt file. (default is _UP, can take _DN, and split gene sets _UP+_DN)

3. Score cells, each gene set will show up as a column in adata.obs.

```
from gssnng import score_cells

q = sc.read_h5ad('gssnng/gssnng/test/data/pbmc3k_processed.h5ad')

sc.pp.neighbors(q, n_neighbors=32)

scores_cells.with_gene_sets(adata=q,                            # AnnData object
                            gene_set_file='cibersort_lm22.gmt', # File path of gene sets
                            groupby='louvain',                  # Will sample neighbors within this group
                            recompute_neighbors=0,              # Rebuild nearest neighbor graph with groups, 0 turns off function
                            score_method='mean_z',              # Method of scoring
                            method_params=dict(),               # Special parameters for some methods 
                            samp_neighbors=27,                  # Number of sampled neighbors for pseudobulk
                            noise_trials=0,                     # ***not used currently***
                            ranked=False,                       # Use ranked data, True or False
                            cores=8)                            # Groups are scored in parallel.
    

sc.pl.umap(q, color=['louvain','T.cells.CD8'], wspace=0.35)
```

## Gene sets

We are following the mSigDB format, where gene sets default to up regulated, which can be marked with suffix _UP
(example: CD8_signature_UP).  In this case, if data is ranked, higher expressed genes have larger ranks. If the 
gene set has suffix _DN (example: CD8_signature_DN), then lowest expressed genes have largest ranks. In the 
of a Z score, the Zs are based on absolute values, so either direction will result in a large Z.

## References

rank biased overlap:  https://arxiv.org/pdf/1408.3587.pdf
singscore:  https://pubmed.ncbi.nlm.nih.gov/30400809/

