
class geneset:
    """ class object that holds a gene set, can be combination of up and dn"""
    name = ''
    info = ''
    genes_up = []  # list of up regulated genes
    genes_dn = []  # list of down regulated genes
    mode = 'NA'    # what type of geneset, up, dn, or both

    def __init__(self, name, info, gs_up, gs_dn, mode):
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


class genesets:
    set_list = []  # list of geneset obj

    def __init__(self, gmt_file):
        ### initially all gene sets are entered as UP
        ### then in the cleaning step, they are assigned as UP,DN,or MATCHED
        txt = open(gmt_file).read().split('\n')
        gd = dict()
        for line in txt:
            bits = line.split('\t')
            if len(bits) >= 3:
                self.set_list.append( geneset(name=bits[0], info=bits[1], gs_up=bits[2:], gs_dn=[], mode='?') )
        self.set_list = self.clean_sets()
        return(None)


    def clean_sets(self):
        ### look for up and dn signifiers
        ### and join up-dn pairs.
        cleaned_sets = []
        for gs1 in self.set_list:
            di = gs1.check_direction()
            if gs1.mode != 'MATCHED' and di in ['UP', 'DN']: ### may need to be matched
                # then we need to check for a pair
                match_flag = 0
                trimmed_name = gs1.name[0: (len(gs1.name)-2)]
                for gs2 in self.set_list:
                    trimmed_name2 = gs2.name[0: (len(gs2.name) - 2)]
                    if ( (trimmed_name == trimmed_name2) and (gs1.name != gs2.name) and
                            (gs1.mode != "MATCHED") and (gs2.mode != "MATCHED") ):
                        gs1.mode = "MATCHED"
                        gs2.mode = "MATCHED"
                        match_flag = 1
                        if di == 'UP': ### gs1 is UP and gs2 is DN
                            cleaned_sets.append(geneset(name=trimmed_name, info=gs1.info,
                                                             gs_up=gs1.genes_up, gs_dn=gs2.genes_up, mode='BOTH'))
                        else: ### gs1 is down and gs2 is up
                            cleaned_sets.append(geneset(name=trimmed_name, info=gs1.info,
                                                             gs_up=gs2.genes_up, gs_dn=gs1.genes_up, mode='BOTH'))
                        break
                if match_flag == 0: ### didn't match, update direction
                    if di == 'UP':
                        cleaned_sets.append(geneset(name=gs1.name, info=gs1.info, gs_up=gs1.genes_up, gs_dn=[], mode='UP'))
                    else:
                        cleaned_sets.append(geneset(name=gs1.name, info=gs1.info, gs_up=[], gs_dn=gs1.genes_up, mode='DN'))
            # then it's default and added as UP
            elif gs1.mode != "MATCHED":
                cleaned_sets.append(geneset(name=gs1.name, info=gs1.info, gs_up=gs1.genes_up, gs_dn=[], mode='UP'))
        return(cleaned_sets)


    def num_genesets(self):
        return(len(self.set_list))

