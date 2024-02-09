if __name__ == '__main__':

        from gssnng import nnsmooth
        import scanpy as sc
        import time

        print("reading data")
        q = sc.datasets.pbmc3k_processed()

        t0 = time.time()
        print('start time: ' + str(t0))

        print("scoring cells")
        q_list = smooth_anndatas.smooth_anndata(
                adata=q,
                groupby='louvain',
                smooth_mode='connectivity',
                recompute_neighbors=0,
                cores=8
            )

        t1 = time.time()

        print("Adata List with SMooTHed counts.")
        print(len(q_list))

        print('end time: ' + str(t1))
        print('TOTAL TIME: ' + str(t1-t0))
        print("done")
