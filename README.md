# gssnng
Gene Set Scoring on the Nearest Neighbor Graph (gssnng) for Single Cell RNA-seq (scRNA-seq).


This package works with Scanpy AnnData objects stored as h5ad files.


  * **Try it out!  ===>>>**  [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/IlyaLab/gssnng/blob/main/notebooks/gssnng_quick_start.ipynb)

  * **Decoupler/Omnipath style API ===>>>** [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/IlyaLab/gssnng/blob/main/notebooks/Scoring_PBMC_data_with_the_GSSNNG_decoupleR_API.ipynb)

  * **See the paper ===>>>** [gssnng](https://academic.oup.com/bioinformaticsadvances/article/3/1/vbad150/7321111?login=false)

  * and finally, [Read the Docs!](https://gssnng.readthedocs.io/en/latest/)



The GSSNNG method is based on using the nearest neighbor graph of cells for data smoothing. This essentially creates 
mini-pseudobulk expression profiles for each cell, which can be scored by using single sample gene set scoring 
methods often associated with bulk RNA-seq. 

Nearest neighbor graphs (NNG) are constructed based on user defined groups (see the 'groupby' parameter below). 
The defined groups can be processed in parallel, speeding up the calculations. For example, a NNG could be 
constructed within each cluster or jointly by cluster *and* sample. Smoothing can be performed using either the 
adjacency matrix (all 1s) or the weighted graph to give less weight to more distant cells.
