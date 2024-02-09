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
       :maxdepth: 2

       Installation
       Scoring Functions
       Example script
       Usage
       Parameters
       Groupby
       Gene sets
       References


`**Notebook using gmt files**  <https://colab.research.google.com/github/IlyaLab/gssnng/blob/main/notebooks/gssnng_quick_start.ipynb>`_

`**Notebook using Decoupler/Omnipath style API** <https://colab.research.google.com/github/IlyaLab/gssnng/blob/main/notebooks/Scoring_PBMC_data_with_the_GSSNNG_decoupleR_API.ipynb>`_

`**Notebook for creating smoothed count matrices**<https://www.google.com>`_

`**See the paper** <https://academic.oup.com/bioinformaticsadvances/article/3/1/vbad150/7321111?login=false>`_


This package works with AnnData objects stored as h5ad files. Expression values are taken from adata.X.
For creating groups, up to four categorical variables can be used, which are found in the adata.obs table.


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

 cp gssnng/gssnng/test/example_smoothing_counts.py  .

 python3.10 example_smoothing_counts.py


Usage
======

See gssnng/notebooks for examples on all methods.

1. Read in an AnnData object using scanpy (an h5ad file).

2. Get gene sets formatted as a .gmt file. (default is UP, also uses _UP,  _DN, and split gene sets _UP+_DN), see below for more details.

3. Score cells, each gene set will show up as a column in adata.obs.

::

   from gssnng import nnsmooth

    q = sc.datasets.pbmc3k_processed()

    q_list = nnsmooth.smooth_adata(adata=q,                    # AnnData object
                                       groupby='louvain',          # Will sample neighbors within this group, can take a list
                                       smooth_mode='connectivity', # Smooths matrix using distance weights from NN graph.
                                       recompute_neighbors=32,     # Rebuild nearest neighbor graph with groups, 0 turns off function
                                       cores=4)                    # Smoothed in parallel.

