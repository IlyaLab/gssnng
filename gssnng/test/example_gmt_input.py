from gssnng import score_cells
import scanpy as sc
import time

def run_gmt_input():

        print("reading data")
        q = sc.datasets.pbmc3k_processed()

        t0 = time.time()
        print('start time: ' + str(t0))

        print("scoring cells")
        score_cells.with_gene_sets(
                adata=q,
                gene_set_file="data/cibersort_lm22.gmt",
                groupby='louvain',
                smooth_mode='connectivity',
                recompute_neighbors=13,
                score_method="geneset_overlap",
                method_params= {'threshold': 0, 'fraction': False},
                ranked=True,
                cores=4
            )

        t1 = time.time()

        print("MEAN GENESET OVERLAP for the T.cells.CD8.up signature")
        print(q.obs.groupby(['louvain'])['T.cells.CD8.up'].mean().reset_index())

        print('end time: ' + str(t1))
        print('TOTAL TIME: ' + str(t1-t0))
        print("done")

#if __name__ == '__main__':
run_gmt_input()