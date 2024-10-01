# gssnng
Gene Set Scoring on the Nearest Neighbor Graph (gssnng) for Single Cell RNA-seq (scRNA-seq).


This package is part of the [scverse ecosystem](https://scverse.org/packages/#ecosystem) and works with Scanpy AnnData objects stored as h5ad files.

  *  **[Read the Docs!](https://gssnng.readthedocs.io/en/latest/)**

  * **Create an AnnData objects with smoothed counts. New Sept. 2024 ===>>>**  [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/IlyaLab/gssnng/blob/main/notebooks/gssnng_data_smoothing.ipynb)

  * **Score cells using .gmt gene set files ===>>>**  [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/IlyaLab/gssnng/blob/main/notebooks/gssnng_quick_start.ipynb)

  * **Score cells using the Decoupler/Omnipath style API ===>>>** [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/IlyaLab/gssnng/blob/main/notebooks/Scoring_PBMC_data_with_the_GSSNNG_decoupleR_API.ipynb)

  * **Read the paper ===>>>** [gssnng](https://academic.oup.com/bioinformaticsadvances/article/3/1/vbad150/7321111?login=false)

The GSSNNG method is based on using the nearest neighbor graph of cells for data smoothing. This essentially creates 
mini-pseudobulk expression profiles for each cell, which can be scored by using single sample gene set scoring 
methods often associated with bulk RNA-seq. 

Nearest neighbor graphs (NNG) are constructed based on user defined groups (see the 'groupby' parameter below). 
The defined groups can be processed in parallel, speeding up the calculations. For example, a NNG could be 
constructed within each cluster or jointly by cluster *and* sample. Smoothing can be performed using either the 
adjacency matrix (all 1s) or the weighted graph to give less weight to more distant cells.


## Scoring Functions

The list of scoring functions:

**geneset_overlap**: For each geneset, number (or fraction) of genes expressed past a given threshold.

**singscore**:      Normalised mean (median centered) ranks (requires ranked data)

**ssGSEA**:         Single sample GSEA based on ranked data.

**rank_biased_overlap**:  RBO, Weighted average of agreement between sorted **ranks** and gene set.

**robust_std**:     Med(x-med / mad), median of robust standardized values (recommend unranked).

**mean_z**:         Mean( (x - mean)/stddv ), average z score. (recommend unranked).

**average_score**:  Mean ranks or counts

**median_score**:   Median of counts or ranks

**summed_up**:      Sum up the ranks or counts.


## Parameters

These parameters are used with the “scores_cells.with_gene_sets” function.:

**adata:**  AnnData object from scanpy.read_*<br/>
AnnData containing the cells to be scored

**gene_set_file:** str[path]<br/>
The gene set file with list of gene sets, gmt, one per line. See `this definition <https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29>`_ .

**groupby:** [str, list, dict]<br/>
either a column label in adata.obs, and all categories taken, or a dict specifies one group.
SEE DESCRIPTION BELOW

**smooth_mode:** "adjacency", "connectivity", or "off"<br/>
Dictates how to use the neighborhood graph.
`adjacency` weights all neighbors equally, `connectivity` weights close neighbors more

**recompute_neighbors:** int<br/>
should neighbors be recomputed within each group, 0 for no, >0 for yes and specifies N

**score_method:** str<br/>
which scoring method to use

**method_params:** dict<br/>
python dict with XGBoost params.

**ranked:** bool<br/>
whether the gene expression counts should be rank ordered

**cores:** int<br/>
number of parallel processes to work through groupby groups


## Method parameters

Some methods have some additional options. They are passed as a dictionary, method_params={param_name, param_value}.:<br/>

singscore:  {'normalization', 'theoretical'}, {'normalization', 'standard'}<br/>
The singscore manuscript describes the theoretical method of standardization which involves determining the theoretical max and minimum ranks for the given gene set.:

rank_biased_overlap:  {'rbo_depth', n}  (n: int)<br/>
Here, n is the depth that is decended down the ranks, where at each step, the overlap with the gene set is measured and added to the score.:

ssGSEA: {'omega': 0.75}<br/>
