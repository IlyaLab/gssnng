# nnggss
Nearest Neighbor Graph Gene Set Scoring (nnggss) for Single Cell RNA-seq (scRNA-seq)

Works with AnnDatas processed with scanpy.

# Score all cells for one gene set (up, down, or up & down)
    res0 = score_cells(
        adata=adata,
        gene_set_up=geneset1,
        gene_set_down=None,
        key_added=key1,
        samp_neighbors=8,
        noise_trials=0,
        mode='average')

# Score all cells & all gene sets (up only):
   res0 = score_cells_all_sets_up(
        adata=adata,
        gene_set_file=geneset_filename, 
        samp_neighbors=8,
        noise_trials=0,
        mode='average')


# example

  import scanpy as sc
  from nnggss.score_all_sets import score_cells_all_sets_up


  def test_score_all_sets_fun(adata, genesets):
      res0 = score_cells_all_sets_up(
          adata=adata,
          gene_set_file=genesets,
          samp_neighbors=8,
          noise_trials=0,
          mode='average')
      return(res0)


  def test_score_all_sets():
      q = sc.read_h5ad('data/pbmc3k_processed.h5ad')
      q2 = q[q.obs['louvain'] == 'NK cells']
      gs = 'data/h.all.v7.2.symbols.gmt'
      score_list = test_score_all_sets_fun(q2, gs)
      print('*************')
      print('post function')
      print(q2.obs[0:5])


  test_score_all_sets()
  print('test score_all_sets done')
