import pandas as pd


class Geneset:
    """ class object that holds a gene set, can be combination of up and dn"""

    def __init__(self, name: str, info: str, gs_up: list, gs_dn: list, mode='NA'):
        """
        :param name: gene set name
        :param info: gene set info
        :param genes_up: list of up regulated genes
        :param genes_dn:  list of down regulated genes
        :param mode:  what type of geneset, up, dn, or both
        """
        self.name = name
        self.info = info
        if mode == 'UP':
            self.genes_up = gs_up
            self.mode = 'UP'
        elif mode == 'DN':
            self.genes_dn = gs_dn
            self.mode = 'DN'
        elif mode == '?':
            self.genes_up = gs_up
            self.mode = '?'
        else:
            self.genes_up = gs_up
            self.genes_dn = gs_dn
            self.mode = 'BOTH'

    def get_name(self):
        return(self.name)

    def up_gene_size(self):
        return (len(self.genes_up))

    def check_direction(self):
        if self.mode == "MATCHED":
            return("MATCHED")
        if 'UP' in self.name.upper():
            self.mode = 'UP'
            return('UP')
        elif 'DN' in self.name.upper():
            self.mode = 'DN'
            return ('DN')
        else:
            self.mode = 'UP'  # default
            return('default')

    def __repr__(self):
        if self.mode == "BOTH":
            return f'Geneset {self.name}\n{self.info}\nUP:' + ",".join(self.genes_up) +  "\nDOWN: " + ",".join(self.genes_dn) 
        if self.mode == "UP":
            return f'Geneset {self.name}\n{self.info}\nUP:' + ",".join(self.genes_up)
        if self.mode == "DN":
            return f'Geneset {self.name}\n{self.info}\nDN:' + ",".join(self.genes_dn)
        if self.mode == "?":
            return f'Geneset {self.name}\n{self.info}\n?:' + ",".join(self.genes_up)


def genesets_from_decoupler_model(df_model: pd.DataFrame, source:str, target: str, weight:str):
    """
    decoupler stores genesets as pd.DataFrames, with a source column (setname) and target column (genenames)
    and an optional weight column (importance + directionality of the gene)

    This function converts it into a Geneset object, adjusting the genes_up/genes_dn based on weight
    """
    geneset_list = []
    for gs_name, df in df_model.groupby(source):
        if weight is not None:
            genes_down = df[df[weight]<0][target].values.tolist()
            genes_up   = df[df[weight]>0][target].values.tolist()
            if len(genes_down) > 0 and len(genes_up) > 0:
                mode='BOTH'
            elif len(genes_down) > 0 and len(genes_up) == 0:
                mode='DN'
            elif len(genes_down) == 0 and len(genes_up) > 0:
                mode="UP"
            else:
                raise ValueError('empty geneset')
        # unweighted geneset, assume those are postivley weighted    
        else:
            genes_up   = df[target].values.tolist()
            genes_down = []
        geneset_list.append(
            Geneset(name=gs_name, info='', gs_up=genes_up, gs_dn=genes_down, mode='?')
        )
    return Genesets(geneset_list)
    

def genesets_from_gmt(gmt_file: str):

    # TODO: maybe use decoupler.read_gmt() -> pd.DataFrame

    # initially all gene sets are entered as UP
    # then in the cleaning step, they are assigned as UP,DN,or MATCHED    
    set_list = []  # list of geneset obj
    with open(gmt_file) as fh:
        for line in fh:
            bits = line.split('\t')
            if len(bits) > 3:
                set_list.append(Geneset(name=bits[0], info=bits[1], gs_up=[x.strip() for x in bits[2:]], gs_dn=[], mode='?'))
            else:
                raise Exception("ERROR: problem with the geneset definitions, make sure it's only tabs.")

    return Genesets(clean_sets(set_list))

def trim_gs_name(gs: Geneset):
    """
    genesets often have names like XYZ.up or XYZ.dn
    this prunes away the UP/DN qualifier
    """
    #TODO: pretty rough, regexp instead?!
    name = gs.name[0: (len(gs.name) - 2)]
    if name[(len(name)-1)] == '.':
        name = name[:name.rfind(".")]
    if name[(len(name)-1)] == '_':
        name = name[:name.rfind("_")]
    return(name)

def clean_sets(set_list):
    """
    when loading from a gmt file, it often contains "matched" genesets, i.e. a signature where some genes are up, some are down
    this is split into two lines/entries in the GMT, one named `somename.up`, the other `somename.down`
    This function pairs up entries in `set_list`

    look for up and dn signifiers
    and join up-dn pairs.
    undefined sets will remain "?", could be up or down
    """
    cleaned_sets = []
    for gs1 in set_list:
        di = gs1.check_direction()
        if gs1.mode != 'MATCHED' and di in ['UP', 'DN']: ### may need to be matched
            # then we need to check for a pair
            match_flag = 0
            trimmed_name = trim_gs_name(gs1)
            for gs2 in set_list:
                trimmed_name2 = trim_gs_name(gs2)
                if ( (trimmed_name == trimmed_name2) and (gs1.name != gs2.name) and
                        (gs1.mode != "MATCHED") and (gs2.mode != "MATCHED") ):
                    gs1.mode = "MATCHED"
                    gs2.mode = "MATCHED"
                    match_flag = 1
                    if di == 'UP': ### gs1 is UP and gs2 is DN
                        cleaned_sets.append(Geneset(name=trimmed_name, info=gs1.info,
                                                            gs_up=gs1.genes_up, gs_dn=gs2.genes_up, mode='BOTH'))
                    else: ### gs1 is down and gs2 is up
                        cleaned_sets.append(Geneset(name=trimmed_name, info=gs1.info,
                                                            gs_up=gs2.genes_up, gs_dn=gs1.genes_up, mode='BOTH'))
                    break
            if match_flag == 0: ### didn't match, update direction
                if di == 'UP':
                    cleaned_sets.append(Geneset(name=gs1.name, info=gs1.info, gs_up=gs1.genes_up, gs_dn=[], mode='UP'))
                else:
                    cleaned_sets.append(Geneset(name=gs1.name, info=gs1.info, gs_up=[], gs_dn=gs1.genes_up, mode='DN'))
        # then it's default and added as ?
        elif gs1.mode != "MATCHED":
            cleaned_sets.append(Geneset(name=gs1.name, info=gs1.info, gs_up=gs1.genes_up, gs_dn=[], mode='?'))
    return cleaned_sets


class Genesets:

    def __init__(self, set_list):
        self.set_list = set_list

    def num_genesets(self):
        return(len(self.set_list))

    def get_gs_names(self):
        gs_names = [gs.get_name() for gs in self.set_list]
        return(gs_names)