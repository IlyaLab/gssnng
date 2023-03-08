
# here we test the ingestion of gene set files gmts

# doubles:
# B.cells.naive.up and .dn  ## in order and next to one another
# T.cells.CD4.memory.resting_dn  ### are out of order
# T_cells_gamma_delta_up ### the two are separated
# Dendritic.cells.activated.up ### doesn't have a pair.
# Macrophages.M1.dn ### no match but is down

### MSIGDB is labeled as UP and DN ###


# rootdir: /home/runner/work/gssnng/gssnng

from gssnng.gene_sets import genesets
def get_number_of_genesets():
    gmt_file = 'test/data/gene_set_test.gmt'
    gslist = genesets(gmt_file)
    return(gslist.num_genesets())

def test_number_of_genesets():
    assert get_number_of_genesets() == 22
