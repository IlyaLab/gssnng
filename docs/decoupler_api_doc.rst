.. GSSNNG documentation master file, created by
sphinx-quickstart on Wed Apr 27 09:20:15 2022.
You can adapt this file completely to your liking, but it should at least
contain the root `toctree` directive.

gssnng using the decoupler-omnipath api
==================================

Gene Set Scoring on the Nearest Neighbor Graph (gssnng) for Single Cell RNA-seq (scRNA-seq).

..
    .. toctree::
       :caption: Table of Contents
       :maxdepth: 2

       ___Usage


`**Notebook using gmt files**  <https://colab.research.google.com/github/IlyaLab/gssnng/blob/main/notebooks/gssnng_quick_start.ipynb>`_

`**Notebook using Decoupler/Omnipath style API** <https://colab.research.google.com/github/IlyaLab/gssnng/blob/main/notebooks/Scoring_PBMC_data_with_the_GSSNNG_decoupleR_API.ipynb>`_

`**Notebook for creating smoothed count matrices**<https://www.google.com>`_

`**See the paper** <https://academic.oup.com/bioinformaticsadvances/article/3/1/vbad150/7321111?login=false>`_

This package works with AnnData objects stored as h5ad files. Expression values are taken from adata.X.
For creating groups, up to four categorical variables can be used, which are found in the adata.obs table.
Gene sets can be provided by using .gmt files or through the OmniPath API (see below).

Scoring functions work with ranked or unranked data (**"your mileage may vary"**):

Method references (singscore, RBO) are below.

Some methods have additional parameters, see below!


Installation
============

Install the package using the following commands::

    python3 -m pip install gssnng

    # or to from github
    python3 -m pip install git+https://github.com/IlyaLab/gssnng



Example script
==============

Copy the script out from the cloned repo and run, check the paths if you get an error.

::

 cp gssnng/gssnng/test/example_decoupler_omnipath_api.py  .

 python3.10 example_decoupler_omnipath_api.py

Usage
======

See gssnng/notebooks for examples on all methods.

1. Read in an AnnData object using scanpy (an h5ad file).

2. Get the model from omnipath via the decoupler API.  You may want to filter out genes negatively associated with the pathway, see the example.

3. Score cells, each gene set will show up as a column in adata.obsm['gssnng_estimate'].

::

   from gssnng import score_cells

    q = sc.datasets.pbmc3k_processed()

    # OmniPath Model #
    model = dc.get_progeny().query('weight>0')

    score_cells.run_gssnng(
        adata, model,
        source='source',target='target', weight='weight',
        groupby="louvain",
        smooth_mode='connectivity',
        recompute_neighbors=32,
        score_method="mean_z",
        method_params={},
        ranked=False,
        cores=6
    )

    # Extracts activities as AnnData object.
    acts_gss = dc.get_acts(adata, obsm_key='gssnng_estimate')

    # Now we can plot the gene set scores
    sc.pl.umap(acts_gss, color=sorted(acts_gss.var_names), cmap='coolwarm')




Scoring Functions
=================

The list of scoring functions::

    geneset_overlap: For each geneset, number (or fraction) of genes expressed past a given threshold.

    singscore:      Normalised mean (median centered) ranks (requires ranked data)

    ssGSEA:         Single sample GSEA based on ranked data.

    rank_biased_overlap:  RBO, Weighted average of agreement between sorted ranks and gene set.

    robust_std:     Med(x-med / mad), median of robust standardized values (recommend unranked).

    mean_z:         Mean( (x - mean)/stddv ), average z score. (recommend unranked).

    average_score:  Mean ranks or counts

    median_score:   Median of counts or ranks

    summed_up:      Sum up the ranks or counts.






Parameters
==========

These parameters are used with the "scores_cells.with_gene_sets" function.::

    adata:  AnnData object from scanpy.read_*
    AnnData containing the cells to be scored

    model: str
    The decoupler gene set model. See Omnipath Wrappers (https://saezlab.github.io/decoupleR/reference/index.html#omnipath-wrappers).

    source: str
    weight: str
    target: str
    Each pathway in OmniPath is a collection of *target* genes from a *source* (i.e. pathway), where each has an interaction *weight*.

    groupby: [str, list, dict]
    either a column label in adata.obs, and all categories taken, or a dict specifies one group.
    SEE DESCRIPTION BELOW

    smooth_mode: "adjacency",  "connectivity", or "off"
    Dictates how to use the neighborhood graph.
    `adjacency` weights all neighbors equally, `connectivity` weights close neighbors more

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


Groupby
=======

The specific neighborhood for each cell can be controlled by using the groupby parameter. In the example
above, by setting groupby='louvain', only cells within a louvain cluster will be considered as being part of the
neighborhood and will available for sampling.

Groupby specifies a column name that's found in the AnnData.obs table, and it can also take a list of column names.
In that case, cells will be grouped as the intersection of categories. For example, using groupby=['louvain','phenotype']
will take cells that are first in a given louvain cluster and then also in a given phenotype group. By also setting
the recompute_neighbors, the nearest neighbor graph is recomputed within this subset of cells. Controlling the
neighborhood leads to more controlled smoothing of the count matrix and is more suitable for downstream comparisons.


Gene sets
=========

We are following the Omnipath wrapper APIs supplied by Decoupler.
See: https://saezlab.github.io/decoupleR/reference/index.html#omnipath-wrappers  for available gene sets.

## Method options

Some methods have some additional options. They are passed as a dictionary, method_params={param_name, param_value}.::

    singscore:  {'normalization', 'theoretical'}, {'normalization', 'standard'}

The singscore manuscript describes the theoretical method of standarization which involves determining the theoretical max and minimum ranks for the given gene set.::

    rank_biased_overlap:  {'rbo_depth', n}  (n: int)

Here, n is the depth that is decended down the ranks, where at each step, the overlap with the gene set is measured and added to the score.::

    ssGSEA: {'omega': 0.75}

The ssGSEA method uses this parameter as a exponent to the ranks. It has been strongly suggested to use 0.75.

*The following methods do not have additional options.*

    robust_std
    mean_z
    average_score
    median_score
    summed_up

References
==========

rank biased overlap:  https://arxiv.org/pdf/1408.3587.pdf

singscore:  https://pubmed.ncbi.nlm.nih.gov/30400809/

anndata: https://anndata.readthedocs.io/en/latest/

MSigDB: https://www.gsea-msigdb.org/gsea/msigdb/

ssGSEA: https://gsea-msigdb.github.io/ssGSEA-gpmodule/v10/index.html

decoupler: https://academic.oup.com/bioinformaticsadvances/article/2/1/vbac016/6544613

omnipath: https://omnipathdb.org/