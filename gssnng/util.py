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
        recompute_neighbors,
        gs_obj,
        score_method,
        ranked,
        method_params
):
    """
    QC on the adata. Need to make sure there's enough neighbors available given the sampling size.

    :param adata: the AnnData object
    :param samp_neighbors: integer, number of neighbors to sample
    """

    if type(method_params) != type(dict()):
        raise Exception('ERROR: please use a dictionary to pass method params')

    if any([xi in adata.obs.columns for xi in gs_obj.get_gs_names()]):
        #raise Exception('ERROR: gene set names in columns of adata.obs, please drop.')
        print("Warning! Dropping gene set names from obs!")
        genesetlist = [x.name for x in gs_obj.set_list]
        for gsi in genesetlist:
            print('dropping: ' + gsi)
            adata.obs.drop(columns=[gsi], inplace=True)

    if 'gssnng_groupby' in adata.obs.columns:
        adata.obs.drop(columns='gssnng_groupby', inplace=True)
        #raise Exception("Error: please drop 'gssnng_groupby' as a column name.")
        print('... and dropping gssnng_groupby column...')

    if ranked == False and score_method == 'singscore':
        raise Exception('ERROR: singscore requires ranked data, set ranked parameter to True')

    if (recompute_neighbors == None) or (recompute_neighbors == 0):
        n_neighbors = adata.uns['neighbors']['params']['n_neighbors'] #[0]# in older AnnData versions need this??
    else:
        n_neighbors = recompute_neighbors

    #if n_neighbors < samp_neighbors:
    #    print('*******')
    #    print('WARNING: Number of neighbors too low for sampling parameter!')
    #    print('Please reduce number of neighbor samples or recompute neighbor graph.')
    #    raise Exception('Samp Neighbors Error')
    #else:
    return(True)


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


def genesets_long_to_gmt(file_name_in, gmt_file_out, gs_mode, set_name_column, symbol_column, desc_column=None, header=True):
    """
    :param file_name_in: file path to long formatted gene set file
    :param gmt_file_out: the path and filename for the newly formatted gmt file
    :param gs_mode: gene set scoring mode, 'up', 'dn', or 'both'
    :param set_name_column: integer, indicating the column containing the get set name
    :param symbol_column: integer, indicates column containing the gene symbol or ID
    :param desc_column: integer, indicates column containing the gene set description, if none, gene set name used.
    :param header: boolean, indicates whether the first row in file_name_in contains column names (a header)
    :return: None, no return
    """
    # if you have a two column file with
    # gene set in set_name_column and gene symbol in symbol_column
    # write a new gmt formatted file.
    txt = open(file_name_in, 'r').read().split('\n')
    outtxt = open(gmt_file_out,'w')
    gd = dict()

    if header: # drop the first row
        txt = txt[1:]

    if desc_column is None:
        desc_column = set_name_column

    i = set_name_column
    j = symbol_column
    k = desc_column

    for line in txt:
        if '.tsv' in file_name_in:
            bits = line.strip().split('\t')
        elif '.csv' in file_name_in:
            bits = line.strip().split(',')
        else:
            print('ERROR, only .csv or .tsv files allowed.')
            return()
        # collect the next gene
        if bits[i] in gd:
            gd[ bits[i] ] = gd[ bits[i] ] + [ bits[j] ]   # past items + next gene
        else:
            gd[ bits[i] ] = [ bits[k] ] + [ bits[j] ]  # desc + first gene
    #
    # now we walk the dictionary
    for ki in gd.keys():
        if gs_mode == 'up':
            outtxt.write(ki+'.up'+'\t'+'\t'.join(gd[ki]) +'\n')
        elif gs_mode == 'both':
            outtxt.write(ki+'.both'+'\t'+'\t'.join(gd[ki]) +'\n')
        else:
            outtxt.write(ki+'.dn'+'\t'+'\t'.join(gd[ki]) +'\n')
    return()
