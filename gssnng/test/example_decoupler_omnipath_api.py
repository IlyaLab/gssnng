
import time
from gssnng import score_cells, gene_sets
import decoupler as dc
import omnipath
import scanpy as sc
#from scipy import sparse


def decoupler_test():

    # load the data
    adata = sc.datasets.pbmc3k_processed()
    #adata.X = sparse.csr_matrix(adata.X) #

    #Get the genesets
    #This is a simple DataFrame, encoding the bipartite graph geneset-name-> genename.
    #Note: These genesets contain genes negatively associated with the signature
    # (i.e. low expression of a gene indicates the presence of a signature).
    # We filter those out here as gssnng doesn't take into account the (negative) weight.
    print("getting omnipath model")
    model = dc.get_progeny().query('weight>0')

    print(model.head())

    t0 = time.time()
    print('gssnng start time: ' + str(t0))

    score_cells.run_gssnng(
        adata, model,
        source='source',target='target', weight='weight',
        groupby="louvain",
        smooth_mode='connectivity',
        recompute_neighbors=12,
        score_method="summed_up",
        method_params={},
        ranked=False,
        cores=4
    )

    # get the results back as an adata
    acts_gss = dc.get_acts(adata, obsm_key='gssnng_estimate')

    # and we can visualize it.
    #sc.pl.umap(acts_gss, color=sorted(acts_gss.var_names), cmap='coolwarm')
    print(acts_gss.obsm['gssnng_estimate'].head())

    t1 = time.time()
    print('end time: ' + str(t1))
    print('TOTAL TIME: ' + str(t1 - t0))
    print("done with decoupler api")

#if __name__ == '__main__':
decoupler_test()