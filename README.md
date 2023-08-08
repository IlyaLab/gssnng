# gssnng

**Try it out!  ===>>>**  [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/IlyaLab/gssnng/blob/main/notebooks/gssnng_quick_start.ipynb)

**See the paper ===>>>** [bioRxiv](https://www.biorxiv.org/content/10.1101/2022.11.29.518384v1)

Gene Set Scoring on the Nearest Neighbor Graph (gssnng) for Single Cell RNA-seq (scRNA-seq).

The problem:  The sparsity of scRNA-seq data creates a poor overlap with many gene sets, 
which in turn makes gene set scoring difficult. 

The GSSNNG method is based on using the nearest neighbor graph of cells for data smoothing. This essentially creates 
mini-pseudobulk expression profiles for each cell, which can be scored by using single sample gene set scoring 
methods often associated with bulk RNA-seq. 

Nearest neighbor graphs (NNG) are constructed based on user defined groups (see the 'groupby' parameter below). 
The defined groups can be processed in parallel, speeding up the calculations. For example, a NNG could be 
constructed within each cluster or jointly by cluster *and* sample. Smoothing can be performed using either the 
adjacency matrix (all 1s) or the weighted graph to give less weight to more distant cells.

This package works with AnnData objects stored as h5ad files. Expression values are taken from adata.X.
For creating groups, up to four categorical variables can be used, which are found in the adata.obs table. 
Gene sets can be provided by using .gmt files or through the OmniPath API (see below).

Scoring functions work with ranked or unranked data (**"your mileage may vary"**):

Method references (singscore, RBO) are below. 

Some methods have additional parameters, see below!

```
    geneset_overlap: Number (or fraction) of genes, past an expression threshold, that overlap with each geneset.

    singscore:       Normalised mean (median centered) ranks (requires ranked data)
    
    ssgsea:          The well known single sample GSEA.
        
    rank_biased_overlap:  RBO, Weighted average of agreement between sorted ranks and gene set.

    robust_std:      Med(x-med / mad), median of robust standardized values (recommend unranked).
    
    mean_z:          Mean( (x - mean)/stddv ), average z score. (recommend unranked).
    
    average_score:   Mean ranks or counts     
    
    median_score:    Median of counts or ranks
    
    summed_up:       Sum up the ranks or counts.
    
```

## Installation from PyPI
```
pip3 install gssnng
```


## Installation from GitHub

```
# also gets you the demo data and some gene sets.
git clone https://github.com/IlyaLab/gssnng

pip install -e gssnng
```

## Example script

Copy the script out from the cloned repo and run, check the paths if you get an error.

```
 cp gssnng/gssnng/test/example_script.py  .
 
 python3 test_gssnng.py
```


## Usage 

See gssnng/notebooks for examples on all methods

1. Read in an AnnData object using scanpy (an h5ad file).

2. Get your gene sets formatted as a .gmt file. (default is undirected, can take _UP,  _DN, and split gene sets _UP+_DN)

3. Score cells, each gene set will show up as a column in adata.obs.

```
from gssnng import score_cells

q = sc.datasets.pbmc3k_processed()

score_cells.with_gene_sets(adata=adata,                         # AnnData object
                            gene_set_file='cibersort_lm22.gmt', # File path of gene sets
                            groupby='louvain',                  # Will sample neighbors within this group
                            smooth_mode='connectivity',         # Smooths matrix using distance weights from NN graph.
                            recompute_neighbors=32,             # Rebuild nearest neighbor graph with groups, 0 turns off function
                            score_method='singscore',           # Method of scoring
                            method_params={'normalization':'theoretical'},  # Special parameters for some methods 
                            ranked=True,                        # Use ranked data, True or False
                            cores=8)                            # Groups are scored in parallel.
    

sc.pl.umap(q, color=['louvain','T.cells.CD8.up'], wspace=0.35)
```

## Parameters

    adata:  AnnData object from scanpy.read_*
    AnnData containing the cells to be scored

    gene_set_file: str[path]
    The gene set file with list of gene sets, gmt, one per line. See `this definition <https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29>`_ .

    groupby: [str, list, dict]
    either a column label in adata.obs, and all categories taken, or a dict specifies one group.
    SEE DESCRIPTION BELOW

    smooth_mode: "adjacency", "connectivity", or "off"
    Dictates how to use the neighborhood graph.
    `adjacency` weights all neighbors equally, `connectivity` weights close neighbors more, `off` does no smoothing.

    recompute_neighbors: int
    should neighbors be recomputed within each group, 0 for no, >0 for yes and specifies N

    score_method: str
    which scoring method to use

    method_params: dict
    python dict with XGBoost params.

    ranked: bool
    whether the gene expression counts should be rank ordered

    cores: int
    number of parallel processes to work through groupby groups

## Groupby

The group of cells available for constructing the neighborhood graph can be controlled by using the groupby parameter. In the example
above, by setting groupby='louvain', only cells within a louvain cluster will be considered as being part of the
neighborhood and will available for sampling.

Groupby specifies a column name (or list of columns) that's found in the AnnData.obs table.
In that case, cells will be grouped as the intersection of categories. For example, using groupby=['louvain','phenotype']
will take cells that are first in a given louvain cluster and then also in a given phenotype group. By also setting
the recompute_neighbors, the nearest neighbor graph is recomputed within this subset of cells. Controlling the
neighborhood leads to more controlled smoothing of the count matrix and is more suitable for downstream comparisons.

## Gene sets

We are following the MSigDB nomenclature, where gene sets default to undirected, but can be marked with the suffix "_UP"
(example: CD8_signature_UP or CD8.signature.up).  In this case, when data is ranked, genes with higher expression have larger ranks. If the 
gene set has suffix "_DN" (example: CD8_signature_DN or CD8.signature.dn), then lowest expressed genes have largest ranks. In the 
use of singscore or Z scores, the undirected case is based on absolute values, so either direction, in the extreme, will result in a large score.

To access OmniPath signatures, see the ["decoupler-style"](https://github.com/IlyaLab/gssnng/blob/main/notebooks/gssnng_decoupler.ipynb) notebook.

## Method options

Some methods have some additional options. They are passed as a dictionary, method_params={param_name: param_value}.

    geneset_overlap: {'threshold': 0, 'fraction': False}
    
The geneset overlap compares the smoothed value or rank of geneset genes to the threshold. If ranked=False, then ranks are used.
If smooth_mode is set to "off", this method can be used with untransformed counts.  If 'percent' is set to True, then 
the overlap will be given as a percentage of the total number of genes in the gene set.

    singscore:  {'normalization': 'theoretical'}, {'normalization': 'standard'}

The singscore manuscript describes the theoretical method of standarization which involves determining the theoretical max and minimum ranks for the given gene set.

    ssGSEA: {'omega': 0.75}
    
The ssGSEA method uses this parameter as a exponent to the ranks. It has been strongly suggested to use 0.75.

    rank_biased_overlap:  {'rbo_depth': n}  (n: int)

Here, n is the depth that we decend down the ranks, where at each step, the overlap with the gene set is measured and added to the score.


*The following methods do not have additional options.*

    robust_std
    mean_z
    average_score
    median_score
    summed_up

## References

rank biased overlap:  https://arxiv.org/pdf/1408.3587.pdf

singscore:  https://pubmed.ncbi.nlm.nih.gov/30400809/

ssGSEA: https://gsea-msigdb.github.io/ssGSEA-gpmodule/v10/index.html

anndata: https://anndata.readthedocs.io/en/latest/

MSigDB: https://www.gsea-msigdb.org/gsea/msigdb/


