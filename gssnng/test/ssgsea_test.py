if __name__ == '__main__':

        from gssnng import score_cells
        import scanpy as sc
        import time

        print("reading data")
        q = sc.datasets.pbmc3k_processed()

        print("computing neighborhood")
        sc.pp.neighbors(q, n_neighbors=64)

        t0 = time.time()
        print('start time: ' + str(t0))

        print("scoring cells")
        score_cells.with_gene_sets(
                adata=q,
                gene_set_file="data/cibersort_lm22.gmt",
                groupby='louvain',
                smooth_mode='connectivity',
                recompute_neighbors=32,
                score_method="ssgsea",
                method_params={'omega':0.25},
                samp_neighbors=27,
                ranked=True,
                cores=8
            )

        t1 = time.time()

        print("MEAN SCORES for the T.cells.CD8.up signature")
        print(q.obs.groupby(['louvain'])['T.cells.CD8.up'].mean().reset_index())

        print('end time: ' + str(t1))
        print('TOTAL TIME: ' + str(t1-t0))
        print("done")
