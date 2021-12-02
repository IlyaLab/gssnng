if __name__ == '__main__':

    import scanpy as sc
    from gssnng.score_cells import with_gene_sets
    import time

    def test_score_all_sets_fun(adata, genesets):
        res0 = with_gene_sets(adata=adata,
                              gene_set_file=genesets,
                              groupby='louvain',
                              recompute_neighbors=0,
                              score_method='singscore',  #'rank_biased_overlap',   #,
                              method_params={'normalization':'theoretical'},  #{'rbo_depth':50},   #{'normalization':'theoretical'},
                              samp_neighbors=27,
                              ranked=True,
                              cores=8)
        return(res0)


    def test_score_all_sets():
        q = sc.read_h5ad('data/pbmc3k_processed.h5ad')
        gs = 'data/cibersort_lm22.gmt' # 'data/gene_set_test.gmt' #'data/cibersort_lm22.gmt'  #
        print("computing knn...")
        sc.pp.neighbors(q, n_neighbors=32)
        print('scoring...')
        t0 = time.time()
        print('start time: ' + str(t0))
        score_list = test_score_all_sets_fun(q, gs)
        print('******DONE*******')
        t1 = time.time()
        print('end time: ' + str(t1))
        print('TOTAL TIME: ' + str(t1-t0))
        print(q.obs.head())
        print(q.obs.groupby(['louvain'])['T.cells.CD8.up'].mean().reset_index())
        #q.write_h5ad('data/pbmc3k_lm22_scores.h5ad')

    test_score_all_sets()
    print('test score_all_sets done')
