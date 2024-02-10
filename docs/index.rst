gssnng
======

Gene Set Scoring on the Nearest Neighbor Graph (gssnng) for Single Cell RNA-seq (scRNA-seq).

Notebook using gmt gene set files ===>>> `Open In Colab <https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/IlyaLab/gssnng/blob/main/notebooks/gssnng_quick_start.ipynb>`_

Notebook using the Decoupler/Omnipath style API ===>>> `Open In Colab <https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/IlyaLab/gssnng/blob/main/notebooks/Scoring_PBMC_data_with_the_GSSNNG_decoupleR_API.ipynb>`_

See the paper ===>>> `gssnng <https://academic.oup.com/bioinformaticsadvances/article/3/1/vbad150/7321111?login=false>`_

Contents
--------

.. toctree::
   gssnng with gene set gmt files <gmt_files_doc>
   gssnng with a decoupler/omnipath style <decoupler_api_doc>
   gssnng to smooth count matrices <smoothing_adatas_doc>

The problem: The sparsity of scRNA-seq data creates a poor overlap with many gene sets, which in turn makes gene set scoring difficult.
The GSSNNG method is based on using the nearest neighbor graph of cells for data smoothing. This essentially creates mini-pseudobulk expression profiles for each cell, which can be scored by using single sample gene set scoring methods often associated with bulk RNA-seq.
Nearest neighbor graphs (NNG) are constructed based on user defined groups (see the 'groupby' parameter below). The defined groups can be processed in parallel, speeding up the calculations. For example, a NNG could be constructed within each cluster or jointly by cluster and sample. Smoothing can be performed using either the adjacency matrix (all 1s) or the weighted graph to give less weight to more distant cells.

This package works with AnnData objects stored as h5ad files. Expression values are taken from adata.X. For creating groups, up to four categorical variables can be used, which are found in the adata.obs table. Gene sets can be provided by using .gmt files or through the OmniPath API (see below).

Scoring functions work with ranked or unranked data ("your mileage may vary")

.. note::

   This project is under active development. Please consider using a named release.

