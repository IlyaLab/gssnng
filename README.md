# gssnng

Try it out!  ===>>>  [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Gibbsdavidl/gssnng/blob/main/notebooks/gssnng_quick_start.ipynb)


Gene Set Scoring on the Nearest Neighbor Graph (gssnng) for Single Cell RNA-seq (scRNA-seq)

Works with AnnData objects stored as h5ad files. Expression values are taken from adata.X.

Scoring functions, works with ranked or unranked data (**"your mileage may vary"**):
```
    singscore:   Normalised mean (median centered) ranks (requires ranked data)
        
    rank_biased_overlap:  Weighted average of agreement between sorted ranks and gene set.

    robust_std:  Med(x-med / mad), median of robust standardized values (recommend unranked).
    
    mean_z:      Mean( (x - mean)/stddv ), average z score. (recommend unranked).
    
    average_score:        Mean ranks or counts     
    
    median_score:         Median of counts or ranks
    
    summed_up:            Sum up the ranks or counts.
```

## Installation

```
git clone https://github.com/Gibbsdavidl/gssnng

pip install -e gssnng
```

## Instructions (see notebooks for examples)

1. Read in an AnnData object using scanpy (an h5ad file).

2. Get gene sets formatted as a .gmt file. (default is undirected, can take _UP,  _DN, and split gene sets _UP+_DN)

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

We are following the mSigDB nomenclature, where gene sets default to undirected, but can be marked with the suffix "_UP"
(example: CD8_signature_UP or CD8.signature.up).  In this case, when data is ranked, genes with higher expression have larger ranks. If the 
gene set has suffix "_DN" (example: CD8_signature_DN or CD8.signature.dn), then lowest expressed genes have largest ranks. In the 
of singscore or a Z score, the undirected case is based on absolute values, so either direction, in the extreme, will result in a large Z.

## References

rank biased overlap:  https://arxiv.org/pdf/1408.3587.pdf

singscore:  https://pubmed.ncbi.nlm.nih.gov/30400809/

