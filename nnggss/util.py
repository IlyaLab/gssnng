import matplotlib
matplotlib.use("Agg")

import numpy
import pandas

import logging

logger = logging.getLogger('singscore')


def read_gene_sets(filepath):
    txt = open(filepath).read().split('\n')
    gd = dict()
    for line in txt:
        bits = line.split('\t')
        if len(bits) >= 3:
            gd[bits[0]] = bits[2:]
    return(gd)


def add_noise(df, n, noise_low, noise_high):
    #df2 = df.copy()
    # draw all random noise in one call
    dfs = pandas.Series(df.gene_counts)
    randmat = pandas.DataFrame(numpy.random.uniform(noise_low, noise_high, (dfs.size, n) ) )
    randmat.index = df.index
    randmat = randmat.apply(lambda x: x + dfs)
    return(randmat)


def get_conn_dist(q, celli, nn):
    # q is an adata
    # celli is cell i in q
    # nn is the number of neighbors
    jdx = numpy.where(q.obsp['distances'][celli].todense() > 0)[1]
    y = [q.obsp['distances'][celli].todense().tolist()[0][ i ] for i in jdx]
    x = [q.obsp['connectivities'][celli].todense().tolist()[0][ i ] for i in jdx]
    df = pandas.DataFrame({'idx':jdx, 'conn':x, 'dist':y})
    sumconn = sum(df['conn'])
    df['prob'] = [x/sumconn for x in df['conn']]
    return(df)


def to_dense_transpose_list(gene_counts):
    gene_mat = gene_counts.todense().transpose().sum(axis=1)
    gdx = numpy.argwhere(gene_mat > 0)
    return( (gene_mat,  [x[0] for x in gdx] ) )
