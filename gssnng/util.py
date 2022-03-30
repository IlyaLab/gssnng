import matplotlib
import numpy as np

matplotlib.use("Agg")
import numpy
import pandas
import logging
logger = logging.getLogger('gssnng')

def error_checking(
        adata,
        samp_neighbors,
        gs_obj,
        score_method,
        ranked
):
    """
    QC on the adata. Need to make sure there's enough neighbors available given the sampling size.

    :param adata: the AnnData object
    :param samp_neighbors: integer, number of neighbors to sample
    """
    if ranked == False and score_method == 'singscore':
        return('ERROR: singscore requires ranked data, set ranked parameter to True')

    n_neighbors = adata.uns['neighbors']['params']['n_neighbors'] #[0] #in older AnnData versions need this??
    if n_neighbors < samp_neighbors:
        print('*******')
        print('WARNING: Number of neighbors too low for sampling parameter!')
        print('Please reduce number of neighbor samples or recompute neighbor graph.')
        return('ERROR')
    else:
        return('OK')


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


## From SingScore
def normalisation(norm_method,
                  gs_mode,
                  score,
                  library_len,
                  sig_len):

    """

    :param norm_method: method of normalisation, standard or theoretical
    :param score_list: list of scores will each be normalised
    :param score: average score (average of list)
    :param library_len: length of the library (int)
    :param sig_len: length of the signature used to generate the scores

    :return: a tuple, containing the normalised score (float) and an array of
    each genes score normalised
    """
    try:
        if norm_method == 'standard':
            norm = score / library_len
        elif norm_method == 'theoretical':
            if gs_mode != '?':
                low_bound = (sig_len + 1) / 2
                upper_bound = ((2*library_len) - sig_len + 1)/2
                norm = (score - low_bound) / (upper_bound - low_bound)
            else:
                low_bound = (np.ceil(sig_len/2) + 1) / 2
                upper_bound = (library_len - np.ceil(sig_len/2) + 1) / 2
                norm = (score - low_bound) / (upper_bound - low_bound)

        return norm

    except:
        logger.exception('Normalisation method must be standard or theoretical.')


def normalisation_rank(norm_method, ranks, library_len, sig_len):
    """
    :param norm_method: method of normalisation, standard or theoretical
    :param ranks: a dataframe of ranks
    :param library_len: length of library (int)
    :param sig_len: length of the signature used to generate the dataframe

    :return: a dataframe with normalised ranks for each gene in each sample
    """
    try:
        if norm_method == 'standard':
            ranks = ranks/library_len
        elif norm_method == 'theoretical':
            low_bound = (sig_len + 1) / 2
            upper_bound = ((2 * library_len) - sig_len + 1) / 2
            ranks = (ranks- low_bound)/(upper_bound-low_bound)

        return ranks

    except:
        logger.exception('Normalisation method must be standard or theoretical.')
