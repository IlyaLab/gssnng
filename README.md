# gssnng
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

# Instructions

1. Read in an AnnData object using scanpy (an h5ad file).

2. Get gene sets formatted as a .gmt file.

3. Score cells, each gene set will show up as a column in adata.obs.

```
    # the scores are written to adata.obs columns #
    score_cells_all_sets(
        adata=adata,               ## The AnnData
        gene_set_file=genesets,    ## The file name containing genes sets in .gmt format.
        score_method='summed_up',  ## The scoring method name 
        set_direction='up',        ## Gene sets direction,
        key_added='GeneSetNames',  ## Only option right now, takes names from gene_set_file
        samp_neighbors=8,          ## Number of neighbors to sample
        noise_trials=0,            ## Noise injection rounds
        mode='average')            ## Method to combine noisy examples
```

## More options

## Score all cells for one gene set (up, down, or up & down)
```
    res0 = score_cells(
        adata=adata,                ## The AnnData
        gene_set_up=geneset1,       ## List of gene symbols, use read_gene_sets(filepath) to get dict of gene sets
        gene_set_down=None,         ##
        score_method='summed_up',   ## Scoring method (see above)
        key_added=key1,             ## Key to add to AnnData, column in adata.obs['key_added'] = gene_set_scores
        samp_neighbors=8,           ## Number of neighbor cells to sample
        noise_trials=0,             ## Noise trials, injects some noise
        mode='average')             ## Method to combine noisy examples
```

## code example ## data included in gssnng/test/data ##
```
  import scanpy as sc
  from gssnng.score_all_sets import score_cells_all_sets

  def test_score_all_sets():
      q = sc.read_h5ad('data/pbmc3k_processed.h5ad')
      q2 = q[q.obs['louvain'] == 'NK cells']
      genesets = 'data/h.all.v7.2.symbols.gmt'
      res0 = score_cells_all_sets(
        adata=adata,
        gene_set_file=genesets,
        score_method='robust_std',
        set_direction='up',
        key_added='GeneSetNames', # only option right now
        samp_neighbors=8,
        noise_trials=0,
        mode='average')
      print('*************')
      print('post function')
      print(q2.obs[0:5])

  test_score_all_sets()
```

## References

rank biased overlap:  https://arxiv.org/pdf/1408.3587.pdf


