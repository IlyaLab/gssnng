if __name__ == '__main__':

    import scanpy as sc
    from gssnng.smooth_anndatas import smooth_anndata
    import time

    def test_return_smoothed(adata):
        res0 = smooth_anndata(adata=adata,
                              groupby='louvain',
                              smooth_mode='adjacency',
                              recompute_neighbors=32,
                              method_params={},
                              cores=4)
        return(res0)


    def test_score_all_sets():
        q = sc.datasets.pbmc3k_processed()
        t0 = time.time()
        print('start time: ' + str(t0))
        data_list = test_return_smoothed(q)
        print('******DONE*******')
        t1 = time.time()
        print('end time: ' + str(t1))
        print('TOTAL TIME: ' + str(t1-t0))
        print(len(data_list))

    test_score_all_sets()
    print('test done')
