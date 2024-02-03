from gssnng import score_cells, gene_sets
import decoupler as dc
import scanpy as sc
import pandas as pd
import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = [4.0, 3.0]

# load the data
adata = sc.datasets.pbmc3k_processed()
adata.X = sparse.csr_matrix(adata.X)

#Get the genesets
#This is a simple DataFrame, encoding the bipartite graph geneset-name-> genename.
#Note: These genesets contain genes negatively associated with the signature
# (i.e. low expression of a gene indicates the presence of a signature).
# We filter those out here as gssnng doesn't take into account the (negative) weight.
model = dc.get_progeny().query('weight>0')


score_cells.run_gssnng(
    adata, model,
    source='source',target='target', weight='weight',
    groupby="louvain", # None
    smooth_mode='connectivity',
    recompute_neighbors=32,
    score_method="mean_z",
    method_params={}, # 'normalization':'standard'
    ranked=False,
    cores=6
)

acts_gss = dc.get_acts(adata, obsm_key='gssnng_estimate')

sc.pl.umap(acts_gss, color=sorted(acts_gss.var_names), cmap='coolwarm')

