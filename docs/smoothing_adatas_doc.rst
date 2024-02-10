.. GSSNNG documentation master file, created by
sphinx-quickstart on Wed Apr 27 09:20:15 2022.
You can adapt this file completely to your liking, but it should at least
contain the root `toctree` directive.

gssnng to make smoothed count matrices
======================================

Gene Set Scoring on the Nearest Neighbor Graph (gssnng) for Single Cell RNA-seq (scRNA-seq).

..
    .. toctree::
       :caption: Table of Contents
       :maxdepth: 1


`**Notebook using gmt files**  <https://colab.research.google.com/github/IlyaLab/gssnng/blob/main/notebooks/gssnng_quick_start.ipynb>`_

`**Notebook using Decoupler/Omnipath style API** <https://colab.research.google.com/github/IlyaLab/gssnng/blob/main/notebooks/Scoring_PBMC_data_with_the_GSSNNG_decoupleR_API.ipynb>`_

`**Notebook for creating smoothed count matrices**<https://www.google.com>`_

`**See the paper** <https://academic.oup.com/bioinformaticsadvances/article/3/1/vbad150/7321111?login=false>`_


This package works with AnnData objects stored as h5ad files. Expression values are taken from adata.X.
For creating groups, up to four categorical variables can be used, which are found in the adata.obs table.


Installation
------------

Install the package using the following commands::

    python3 -m pip install gssnng

    # or to from github
    python3 -m pip install git+https://github.com/IlyaLab/gssnng



Example script
--------------

Copy the script out from the cloned repo and run, check the paths if you get an error.

::

 cp gssnng/gssnng/test/example_smoothing_counts.py  .

 python3.10 example_smoothing_counts.py


Usage
-----

See gssnng/notebooks for examples on all methods.

1. Read in an AnnData object using scanpy (an h5ad file).

2. Get gene sets formatted as a .gmt file. (default is UP, also uses _UP,  _DN, and split gene sets _UP+_DN), see below for more details.

3. Score cells, each gene set will show up as a column in adata.obs.

::

   from gssnng import smoothing

    q = sc.datasets.pbmc3k_processed()

    q_list = smoothing.smooth_adata(adata=q,                    # AnnData object
                                       groupby='louvain',          # Will sample neighbors within this group, can take a list
                                       smooth_mode='connectivity', # Smooths matrix using distance weights from NN graph.
                                       recompute_neighbors=11,     # Rebuild nearest neighbor graph with groups, 0 turns off function
                                       cores=4)                    # Smoothed in parallel.



Parameters
----------

These parameters are used with the "scores_cells.with_gene_sets" function.::

    adata:  AnnData object from scanpy.read_*
    AnnData containing the cells to be scored

    groupby: [str, list, dict]
    either a column label in adata.obs, and all categories taken, or a dict specifies one group.
    SEE DESCRIPTION BELOW

    smooth_mode: "adjacency", "connectivity", or "off"
    Dictates how to use the neighborhood graph.
    `adjacency` weights all neighbors equally, `connectivity` weights close neighbors more

    recompute_neighbors: int
    should neighbors be recomputed within each group, 0 for no, >0 for yes and specifies N

    cores: int
    number of parallel processes to work through groupby groups


Groupby
-------

The specific neighborhood for each cell can be controlled by using the groupby parameter. In the example
above, by setting groupby='louvain', only cells within a louvain cluster will be considered as being part of the
neighborhood and will available for sampling.

Groupby specifies a column name that's found in the AnnData.obs table, and it can also take a list of column names.
In that case, cells will be grouped as the intersection of categories. For example, using groupby=['louvain','phenotype']
will take cells that are first in a given louvain cluster and then also in a given phenotype group. By also setting
the recompute_neighbors, the nearest neighbor graph is recomputed within this subset of cells. Controlling the
neighborhood leads to more controlled smoothing of the count matrix and is more suitable for downstream comparisons.


References
----------

rank biased overlap:  https://arxiv.org/pdf/1408.3587.pdf

singscore:  https://pubmed.ncbi.nlm.nih.gov/30400809/

anndata: https://anndata.readthedocs.io/en/latest/

MSigDB: https://www.gsea-msigdb.org/gsea/msigdb/

ssGSEA: https://gsea-msigdb.github.io/ssGSEA-gpmodule/v10/index.html

decoupler: https://academic.oup.com/bioinformaticsadvances/article/2/1/vbac016/6544613

omnipath: https://omnipathdb.org/
