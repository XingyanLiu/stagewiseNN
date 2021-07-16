# -*- coding: utf-8 -*-
"""
Created on Mon Dec 23 10:40:54 2019

@author: xyliu
"""

import os
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse, io
from matplotlib import pyplot as plt
import seaborn as sns

WORKDIR = Path(r'D:\Users\xyliu\003')
os.chdir(WORKDIR)

import funx as fx
import Build as B



# In[]


DATADIR = Path(r'E:\Users\xyliu\data003\amph')
DATADIR_RAW = DATADIR / 'RAW'
FORMALDIR = DATADIR / 'afterQC_formal'

FORMALDIR_NEW = DATADIR / '_formal_data'
FORMALDIR_NEW_ENB = FORMALDIR_NEW / 'embryos'
FORMALDIR_NEW_ADT = FORMALDIR_NEW / 'adults'

LINEDIR = DATADIR / 'lineage' 
DATADIR_AD = DATADIR / 'Adult'
GENEDIE = DATADIR / 'genes'

CMAP_CLASS = sc.pl.palettes.zeileis_26
CMAP_STAGE0 = 'Spectral'
CMAP_STAGE = 'plasma_r' #'Spectral'
CMAP_VALUE = fx.diy_cmap_grey_bg()

Stages0 = ('E1', 'E2', 'E3', 'E4', 'E5', 'E6', 'E7', 'E8', 'E9', 
           'E10', 'E11', 'E12', 'E13', 'E14', 'E15')
#Snames = [stg for stg in Stages0 

usePolyT = ('E1', 'E2', 'E3', 'E4', 'E5', 'E6', 'E7', 'E8', 'E10',)

Stages = ('E6', 'E7', 'E8',  'E9', 
           'E10', 'E11', 'E12', 'E13', 'E14')
StageNames = ("B", "G3", "G4", "G5", "G6", "N0", "N1", "N3", "L0" )
StageNamesAll = ("2cell", "4cell", "8cell", "32cell", "256cell",
                 "B", "G3", "G4", "G5", "G6", "N0", "N1", "N3", "L0" )


Tissues = ('T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7', 
           'T8', 'T9', 'T10', 'T11', 'T12', 'T13', 'T14')

Primers = ('random', 'polyT')

StageNameDict = dict(
        E6 = "B",
        E7 = "G3",
        E8 = "G4",
        E9 = "G5",
        E10 = "G6",
        E11 = "N0",
        E12 = "N1",
        E13 = "N3",
        E14 = "L0"
        )
TissueNameDict = {
         'T1': 'Tentacles',
         'T2': 'Neural Tube',
         'T3': 'Male Gonad',
         'T4': 'Endostyle',
         'T5': 'Epidermis',
         'T6': 'Hepatic diverticulum',
         'T7': 'Fore End',
         'T8': 'Hind End',
         'T9': 'Gill Branchia',
         'T10': 'End Gut',
         'T11': 'Rostral Side',
         'T12': 'Female Gonad',
         'T13': 'Muscle',
         'T14': 'Notochord'
         }
TissueNamesUse = [
     'Female Gonad',
     'Hind End',
     'End Gut',
     'Fore End',
     'Rostral Side',
     'Gill Branchia',
     'Male Gonad',
     'Notochord',
     'Epidermis',
     'Neural Tube',
     'Tentacles',
     'Endostyle',
#     'Muscle',
#     'Hepatic diverticulum'
 ]

LineageOrdAll = (
        "B_0",
        "B_1",
        "B_2",
        "Primordial germ cells",
        "Epithelial ectoderm",
        "Neural ectoderm",
        "Notochord",
        "Mesoderm",
        "Unassigned",
        "Endoderm",
        )

LineageOrd = LineageOrdAll[-7:]
LineageAnno = dict(
  E7_0 = "Epithelial ectoderm",
  E7_1 = "Endoderm",
  E7_2 = "Mesoderm",
  E7_3 = "Notochord",
  E7_4 = "Neural ectoderm",
  E7_5 = "Unassigned",
  E7_6 = "Primordial germ cells"
)


LineageAnnoAll = dict(
        E6_0 = 'B_0',
        E6_1 = 'B_1',
        E6_2 = 'B_2',
        E6_3 = "Primordial germ cells",
        E7_0 = "Epithelial ectoderm",
        E7_1 = "Endoderm",
        E7_2 = "Mesoderm",
        E7_3 = "Notochord",
        E7_4 = "Neural ectoderm",
        E7_5 = "Unassigned",
        E7_6 = "Primordial germ cells"
)


MarkerDictShow = {
#        'bf_00013097': 'DNMT1',
        'bf_00019348': 'Hu_Elav', # neural ectoderm
#        'bf_00016583': 'Tbx2',
#        'bf_00004829': 'Nanos',
#        'bf_00008183': 'LRIG3',
        'bf_00002546': 'Sox17', # endoderm
        'bf_00006571': 'Wnt8', # Mesoderm
        'bf_00016125': 'FoxJ1', # Epithelial ectoderm
        'bf_00023325': 'Keratin', # Epithelial ectoderm
        'bf_00020985': 'CHRD', # notochord
        'bf_00012984': 'Cdx',
        }

# In[]

_candi0 = ['#596e79', '#b3b3b3', '#c7b198'] + \
           ['#40bad5', #'#377eb8', #'#66c2a5', 
           '#984ea3', #'#c060a1', #'#bc658d', # purple
#           '#f57b51',
           '#36622b', #'#00bd56', #'#009975', #'#007065', # green
#             '#91bd3a', #'#a6d854', #yellow-green
#           '#50d890', # green
           '#035aa6', #'#4f98ca', #'#40bad5', #'#2fc4b2', ## blue
           
           '#fcbf1e', #'#ffd800', #'#ffd31d', #  #yellow
           '#af0404', 
#            '#e78ac3', # pink
#           '#fa9191', #'#fc8d62', # pink-orenge
           '#dd7631', # pink-brown
#                 '#e41a1c', '#8da0cb', 
                  ][: 7]
_candi1 = ['#596e79', '#b3b3b3', '#c7b198'] +\
           ['#40bad5', '#984ea3', '#36622b', #'#00bd56', #'#009975', #'#007065', # green
           '#035aa6', '#fcbf1e', '#af0404', 
           '#dd7631']


LineageColorCandi = pd.DataFrame({
        'candi0': _candi0,
        'candi1': _candi1,
        }, \
    index = LineageOrdAll)

LineageColorFormal = pd.Series(
        _candi1 + ['#40bad5'], index = list(LineageOrdAll) + ['B_3'])
           
TissueColors0 = [
         '#8c564b',
         '#ff7f0e',
         '#2ca02c',
         '#d62728',
         '#9467bd',
         '#c49c94', #'Hepatic diverticulum'
         '#e377c2',
         '#bcbd22',
         '#17becf',
         '#aec7e8',
         '#ffbb78',
         '#98df8a',
         '#ff9896', #'Muscle'
         '#1f77b4',
#         '#c5b0d5',
#         '#f7b6d2',
#         '#c7c7c7',
#         '#dbdb8d',
#         '#9edae5'
#         '#7f7f7f',
         ]
TissueColors = [
         '#8c564b',
         '#ff7f0e',
         '#2ca02c',
         '#d62728',
         '#9467bd',
#         '#c49c94', #'Hepatic diverticulum'
         '#e377c2',
         '#bcbd22',
         '#17becf',
         '#aec7e8',
         '#ffbb78',
         '#98df8a',
#         '#ff9896', #'Muscle'
         '#1f77b4',
         ]






# In[]
#==========================================================

def sampleNames(stages = Stages0, polyT = usePolyT):
    
    snames = list(stages)
    for i, sn in enumerate(snames):
        if sn in polyT:
            snames[i] = sn + '_polyT'
    return snames


def rename_genes(gene_ids, mapping=None, tolist=False):
    
    print('the old gene names:', gene_ids[:5])
    mapping = fx.gene_id_name_mapping('fish') if mapping is None else mapping
    gnames = pd.Series(gene_ids).apply(lambda x: fx.rename_gene_fish(x, mapping))
    print('the new gene names:', gnames[:5])
    return gnames.to_list() if tolist else gnames

def make_gene_dict(s: str, n2id=True) -> dict:
    s = s.split()
    if not n2id:
        s = s[::-1]
    gene_dict = {}
    for k, v in zip(s[1::2], s[::2]):
        if k in gene_dict.keys():
            gene_dict[k] += f',{v}'
        else:
            gene_dict[k] = v
    return gene_dict
# In[]
'''
==========================================
        0-0: Getting datasets
==========================================
'''


def GeneAnnotations():
    df =  _load_data(DATADIR, 'gene_annos_merged_srt_bed.tsv',
                      pd.read_csv, sep='\t')
    df_note =  _load_data(DATADIR, 'gene_annos_full.csv',
                  pd.read_csv, index_col=0)
    df.set_index('ID', drop=False, inplace=True)
    df['Note'] = df_note['Note']
    return df


def SaveGenesWithAnnos(genes, resdir, name='tmp', cols=None):
    '''
    genes: list-like
    '''
    geneAnno = GeneAnnotations()

    if cols is None:
        cols = ['name', 'anno_human', 'Note', 'chr'] 
    geneAnno.loc[genes, cols].to_excel(resdir / f'{name}_genes_info.xlsx', 
                header=True)
    print(f'save genes with annotaions into:\n\t{resdir}')
    

#def GeneAnnotations1():
#    df =  _load_data(DATADIR, 'gene_annos_full.csv',
#                      pd.read_csv,)
#    df.set_index('ID', drop=False, inplace=True)
#    return df


def DataEmb(primer='random', qc=100, tail='', datadir=None):
    '''
    qc=100 (150 actually) or 200
    tail = Union['_gHVGs', '', ]
    '''
    # E:\Users\xyliu\data\amphioxus\MTX\whole_emb_mtx\nCcut_100\RNAcap
    datadir = DATADIR / f'nCcut_{qc}' / f'RNAcap{tail}' if datadir is None else datadir
    fname = f'{primer}.h5ad'
    return _load_data(datadir, fname)


def data_emb_at(stage = 'E1', primer='random', qc=100, tail='', datadir=None):
    '''
    primer = Union['random', 'polyT']
    tail = Union['', '_custom', ]
    '''
    # E:\Users\xyliu\data\amphioxus\MTX\whole_emb_mtx\nCcut_100\batch_random
    datadir = DATADIR / f'nCcut_{qc}' / f'batch_{primer}{tail}' if datadir is None else datadir
    if isinstance(stage, int): stage = f'E{stage}'
    fname = f'{stage}.h5ad'
    return _load_data(datadir, fname)

def DataEmbAt(stages=Stages, primer='random', qc=100, tail='', datadir=None):

    datadir = DATADIR / f'nCcut_{qc}' / f'batch_{primer}{tail}' if datadir is None else datadir
    adatas = [data_emb_at(stg, datadir=datadir) for stg in stages]
    return adatas

def RawDataEmb(qc = 100):
    
    datadir = DATADIR / f'nCcut_{qc}'
    fname = 'labeled_afterQC.h5ad'
    return _load_data(datadir, fname)


def HVGsEmbAt():
    
    pass

def HVGsEmb(primer='random', qc=100, tail='_gHVGs'):
    # E:\Users\xyliu\data\amphioxus\MTX\whole_emb_mtx\nCcut_100\RNAcap_gHVGs
    datadir = DATADIR / f'nCcut_{qc}' / f'RNAcap{tail}'
    fname = f'{primer}_hvgs_id.csv'
    return _load_data(datadir, fname, pd.read_csv, header=None).iloc[:, 0].to_list()


def _load_data(datadir, fname, _method=sc.read_h5ad, **kwds):
    fpath = datadir / fname
    print(f'Loading data from file:\n {fpath}')
    adata = _method(fpath, **kwds)
    return adata


# In[]
'''
        0-1: Manuplating datasets
==========================================
'''







def extract_meta(adatas, keys = 'batch', tags = None, add_key='tag'):
    '''
    tags: a list of snames of each dataset
    add_key: another column will be added for identifying data origins
    '''
    if isinstance(keys, str):
        keys = [keys]
    metas = [adt.obs[keys] for adt in adatas]
    meta_all = pd.concat(metas, axis=0)
    if tags is not None:
        tag = []
        for i, adt in enumerate(adatas):
            tag += [tags[i]] * adt.shape[0]
        meta_all[add_key] = tag
    n = meta_all.shape[0]
    print(f'Concatenated meta data: {n} observations with columns:\n{meta_all.columns}')
    return meta_all



  

# In[]
'''
==========================================
        0-2: Adjusting clusters
==========================================
'''
from sklearn.neighbors import KNeighborsClassifier

'''
>>> X = [[0], [1], [2], [3]]
>>> y = [0, 0, 1, 1]
>>> neigh = KNeighborsClassifier(n_neighbors=3)
>>> neigh.fit(X, y)
KNeighborsClassifier(...)
>>> print(neigh.predict([[1.1]]))
[0]
'''
def labelByKNN(adt, is_query, key_y, #='leiden', 
               key_x='X_pca', 
               inplace=False, rename_cat=False,
               k = 10, **kwds):
    '''
    adt: sc.AnnData
    is_query: pd.Series; bool; index must be the same as adt.obs_names
    '''
    X = adt.obsm[key_x]
    y = adt.obs[key_y] if inplace else adt.obs[key_y].copy() 
    knnCl = KNeighborsClassifier(n_neighbors=k, **kwds)
    knnCl.fit(X[~ is_query.values, :], y[~ is_query])
    y_que = knnCl.predict(X[is_query.values, :])
    y[is_query] = y_que#knnCl.classes_[y_que]
    y = y.astype(str).astype('category')
    if rename_cat:
        n_groups = len(y.cat.categories)
        new_cat = list(map(str, range(n_groups)))
        y.cat.categories = new_cat
        print(f'categories are renamed as:\n{new_cat}')
    if inplace:
        adt.obs[key_y] = y
    else:
        return y#, knnCl
    










# In[]

'''
==========================================
        1: Simple Tree from clusters
==========================================
'''
import networkx as nx

def TreeByStages(df, stages=None, ranks=None, cross_stg=True):
    '''
    df: dataframe with values as a symetric matrix, represents the correlations
        between clusters
        - columns (e.g.):
          "E1_0", "E2_0", "E2_1", ..., "E13_0", ..., "E13_12"  
    '''
    node_names = df.columns
    stages = []
    labels = []
    for nn in node_names:    
        stg, lb = nn.split('_')
        stages += [int(stg.strip('E'))]
        labels += [int(lb)]
    if ranks is None:
        ranks = rank_of_stages(stages)
    attrs = pd.DataFrame(dict(rank = ranks, stage = stages, label=labels), 
                         index=node_names)
    # initiating the graph
    G = nx.Graph()
    ## adding nodes
    G.add_nodes_from(node_names)
    nx.set_node_attributes(G, attrs.T.to_dict())
        
    ## adding edges
    for nd, rk in G.nodes.data('rank'):
        if rk == 0: continue
        # only consider nodes from the former stage
        pa = df[nd][attrs['rank'] == rk - 1].idxmax()
        w = df[nd][pa]
        G.nodes[nd]['parent'] = pa # set `parent` attribute for `nd` as `pa`
        G.add_edge(nd, pa, weight = w, diff = 1) # diff: rank difference
        if cross_stg:
            # the node most correlated with `nd` may be from earlier stages
            pa0 = df[nd][attrs['rank'] < rk].idxmax()
            if pa0 != pa:
                G.add_edge(nd, pa0, weight=df[nd][pa0], diff = rk - attrs['rank'][pa0])
    return G

def dfFromG(G, default_root=np.nan):
    node_pa = G.nodes.data('parent', default=default_root)
    df_tree = pd.DataFrame(node_pa, columns=['node', 'parent'])
    df_tree['label'] = df_tree['node'].copy()
    df_tree['stage'] = df_tree['node'].apply(lambda x: x.split('_')[0])
    df_tree['stage_int'] = df_tree['stage'].apply(lambda x: int(x.strip('E')))
    return df_tree



def plot_tree_by_stages(G, mode='ud c', figsize=(10, 15), 
                        node_color_by = 'rank',
                        save=None, show=True,
                        **kwds):
    '''
    stage-by-stage tree layout
    '''
    pos = pos_tree_by_stage(G, mode=mode)
    
    node_color = list(nx.get_node_attributes(G, node_color_by).values())
    # diff: rank difference between 2 nodes
    skip_edges = [(u, v) for u, v, diff in G.edges.data('diff') if diff > 1]
    near_edges = [(u, v) for u, v, diff in G.edges.data('diff') if diff == 1]
    
    plt.subplots(figsize=figsize)
    nx.draw_networkx_nodes(G, pos, with_labels=True, 
                            node_color=node_color, 
                            cmap = 'Spectral', **kwds)
    if len(skip_edges) >= 1:
        nx.draw_networkx_edges(G, pos, skip_edges, edge_color='grey',
                            width=1.2, alpha = 0.5, **kwds)
    nx.draw_networkx_edges(G, pos, near_edges, **kwds)#edge_color = 'grey', )
    nx.draw_networkx_labels(G, pos, font_size=10)
    fig = plt.gcf()
    plt.axis('off')
    
#    plt.tight_layout()
    if save: # BUG!!!
        fname = save if isinstance(save, (Path, str)) else 'lineage_tree_by_stage_temp.pdf' 
        fig.savefig(fname, bbox_inches='tight')
        print(f'figure has been saved into {fname}')
    if show: plt.show()


def plot_tree_by_express(G, mode='ud c', figsize=(10, 15), 
                        node_color_by = 'rank',
                        save=None, show=True,
                        **kwds):
    '''
    stage-by-stage tree layout
    '''
    pos = pos_tree_by_stage(G, mode=mode)
    
    node_color = list(nx.get_node_attributes(G, node_color_by).values())
    # diff: rank difference between 2 nodes
    skip_edges = [(u, v) for u, v, diff in G.edges.data('diff') if diff > 1]
    near_edges = [(u, v) for u, v, diff in G.edges.data('diff') if diff == 1]
    
    plt.subplots(figsize=figsize)
    nx.draw_networkx_nodes(G, pos, with_labels=True, 
                            node_color=node_color, 
                            cmap = 'Spectral', **kwds)
    nx.draw_networkx_edges(G, pos, skip_edges, edge_color='grey',
                            width=1.2, alpha = 0.5, **kwds)
    nx.draw_networkx_edges(G, pos, near_edges, **kwds)#edge_color = 'grey', )
    nx.draw_networkx_labels(G, pos, font_size=10)
    fig = plt.gcf()
    plt.axis('off')
    
#    plt.tight_layout()
    if save: # BUG!!!
        fname = save if isinstance(save, (Path, str)) else 'lineage_tree_by_stage_temp.pdf' 
        fig.savefig(fname, bbox_inches='tight')
        print(f'figure has been saved into {fname}')
    if show: plt.show()
    
    
def plot_tree_auto(G, pos = None, init = False, layout='kk', **kwds):
    
    if layout == 'kk':
        pos0 = pos_tree_by_stage(G) if init else None
        pos = nx.drawing.kamada_kawai_layout(G, pos=pos0)
#    pos = nx.drawing.spring_layout(G, pos=pos0) # by default
#    pos = nx.drawing.spectral_layout(G) # bad layout!
    near_edges = [(u, v) for u, v, diff in G.edges.data('diff') if diff == 1]
    node_color = list(nx.get_node_attributes(G, 'rank').values())
    
    plt.subplots(figsize=(12, 12))
    nx.draw(G, pos, with_labels=True, edge_list = near_edges,
            node_color=node_color, 
            font_size=8, cmap = 'Spectral')

    plt.show()
    
    

def rank_of_stages(stages):
    '''
    stages: a ordered list 
    return:
        ranks: a list of integers that encode the rank od stages, with the same
        length of `stages`
    '''
    ranks = [0]
    for i, stg in enumerate(stages):
        if i == 0: continue
        if stg == stages[i - 1]:
            ranks += [ranks[i - 1]]
        else:
            ranks += [ranks[i - 1] + 1]
            
    ### Better make a mapping instead !!!
    return ranks

def pos_tree_by_stage(G, ranks=None, mode='ud c'):
    '''
    G: class networkx.Graph 
    ranks:
    
    mode: str; combiniation of "ud"
        ud: (short of "up-to-down")
        rl: (short of "right-to-left")
        c: (short of "center")whether to centher the X-coordinates.
    NOTE:
        There should be at least one node with rank 0 !
    
    '''
    ranks = dict(G.nodes.data('rank')) if ranks is None else ranks
    pos = pd.DataFrame({'x': 0, 'y': ranks, }, index=list(G))
    uni_ranks = np.unique(pos['y'], ) # autometically sorted (ascending!)
    pos.head(10)
    ## setting the root(s) (whose rank are 0)
    roots = pos.index[pos['y'] == 0]
    for i, nd in enumerate(roots): 
        pos['x'][nd] = i
    
    for rk in uni_ranks:
        if rk == 0: continue
        ids = pos['y'] == rk
        nodes = pos.index[ids]
        for i, nd in enumerate(nodes):
            pa = G.nodes[nd]['parent']
            pos['x'][nd] = pos['x'][pa]
        
        for i, v in enumerate(pos['x'][ids].argsort()):# 妈的太烧脑了！
            pos['x'][nodes[v]] = i
    if 'ud' in mode:
        pos['y'] = pos['y'] * (-2)
    if 'rl' in mode: 
        pos['x'] = pos['x'] * (-1)
    if 'c' in mode:
        pos['x'] = pos.groupby('y')['x'].apply(lambda x: x - np.mean(x))
    pos = df2dict(pos, axis=0)
#    pos = {nd: [pos.x[nd], pos.y[nd]] for nd in pos.index}
    return pos

def df2dict(df, axis=0):
    '''
    transform a DataFrame into a dict, with the index as keys, and 
    the lists of values from `df` as the values of the dict
    df: pd.DataFrame
    axis: 0 or 1
        axis=0: take the df.index as the keys of the dict
        axis=1: take the df.columns as the keys of the dict
    '''
    if axis == 0: df = df.T
    
    dct = {c: df[c].values.tolist() for c in df.columns}
    
    return dct



# In[]
    
'''
==========================================
        Single Lineage analysis
==========================================
Step 0:
    Separate groups
    
Step 1:
    Re-embedding: neighbors -> UMAP
    Re-clustering
    Purify (optional)
    DE (marker finding)
    
Step 2:
    
    pseudotime inference


'''



def PH_analysis(ph, n_pcs = 50, nneigh=10, metric='cosine', 
                min_dist=0.25, res=0.6, save=True):

    
    ph.Neighbors(nneigh=nneigh, n_pcs=n_pcs, metric=metric)
    ph.Umap(min_dist=min_dist, plot_color='leiden')
    ph.Clustering(res=res, keep_results=True)
    ph.DE(plot=True, save=save)#, method='wilcoxon')

def PH_paga(ph, groups='leiden'):
    sc.tl.paga(ph.adata, groups=groups)
    sc.pl.paga_compare(ph.adata, save=f'_{ph.name}.png')
    sc.pl.paga(ph.adata, save = f'_{ph.name}.png')
    #ph.Umap(init_pos='paga')
    
def PH_pseudotime(ph, root_group, groupby='leiden', diff_comps=15, n_dcs=10, ):

    sc.tl.diffmap(ph.adata, n_comps=diff_comps)
    sc.pl.diffmap(ph.adata, color=groupby, components=['1,2', '3,4'])
    
    ph.adata.uns['iroot'] = group_center(ph.adata, root_group, groupby)
    sc.tl.dpt(ph.adata, n_dcs=n_dcs)
    PH_vis_group_dpt(ph, groupby=groupby)
    
    
    
def PH_vis_group_dpt(ph, groupby='leiden'):
    ph.vis(color = [groupby, 'dpt_pseudotime'],
           legend_loc='on data',
           save=f'_{groupby}_dpt_{ph.name}.pdf')
    

def PH_vis_paga_compare(ph, ftype='.png'):
    sc.pl.paga_compare(ph.adata, save=f'_{ph.name}{ftype}')
    
def PH_vis_dpt(ph, **kw):
    ph.vis(color ='dpt_pseudotime',
           legend_loc='on data',
           save=f'_dpt_{ph.name}.pdf', **kw)


def PH_dynamic(ph, path_ord, top_n_gene = 5,
               figsize=(8, 16), save_genes=True, 
               saveph=False, **kw):
    
    if isinstance(path_ord, str):
        path_name = path_ord
        path_ord = list(path_ord)
    elif isinstance(path_ord, list):
        path_name = '-'.join(path_ord)
    else:
        print('`path_ord` should be either string or a list!')
        
    genes = ph.markers_top(top_n_gene, groups=path_ord)
    if save_genes:
        pd.Series(genes).to_csv(ph.resdir / ('dynamic_genes_path%s.csv'%path_name), 
                  header=False, index=False)
#    len(genes)
    fig, axs = plt.subplots(figsize=figsize,)
    sc.pl.paga_path(ph.adata, 
                    nodes = list(path_ord),
                    keys = genes,
                    ax = axs,
                    normalize_to_zero_one =True,
#                    groups_key=
#                    palette_groups=fx.get_colors('Spectral', len(path_ord)),
#                    vmax=2.5,
                    n_avg = 50, 
                    show_node_names=False,
                    save = '%s_top%d.pdf'%(path_name, top_n_gene),
#                    save = '%s_top%d_avg50.png'%(path_ord, n),
                    ytick_fontsize=11,
                    legend_fontsize=12,
                    **kw
                    )
    if saveph: ph.save()
    
    
    
def group_center(adata, group, groupby='leiden'):
    ''' 
    find the center of a group, return the index
    temperally take the node with max degree.
    '''
#    indices_ = ph.adata.obs[groupby]  == group
    indices = np.flatnonzero(adata.obs[groupby]  == group)
    A = adata.uns['neighbors']['connectivities'][indices, :]

    # compute the degrees
    d = A.sum(axis=1).A1
    center_id = indices[np.argmax(d)]
    print(f'Auto-selected node {center_id}, with max degree {d.max()}')
    return center_id
    



# In[]
'''             Plotting functions
====================================================
'''  

#def vis_groups(adata, color='leiden'):
#    
#    tt = f'{stage} (adt.shape[0] cells)'
#    sc.pl.umap(adt, color=color, palette = CMAP_CLASS, title = tt,\
#           legend_loc = 'on data', save=figdir / f'{stage}_{new_key}')

def plotMyUmap(adt, color='leiden', stage=None, 
               palette = CMAP_CLASS, # 'tab20'
               tag='', tt=None,
               legend_loc = 'on data', # 'right margin'
               save=True,
               figtype='.pdf'):
    if stage is None:
        stage = adt.obs['stage'][0]
    tt = f'{stage} ({adt.shape[0]} cells)' if tt is None else tt
    fn = f'_{stage}_{color}{tag}{figtype}' if save else None
    sc.pl.umap(adt, color=color, 
               palette = palette, 
               title = tt,
               legend_loc = legend_loc, 
               save=fn)  


def plot_matrix(mat, saveto=None, 
            figsize=(7.5, 5), annot=True, fmt=".2f",
            ax=None, **kwds):
    
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(mat, #font_size=9,
                annot=annot, fmt=fmt, ax=ax, **kwds)
    if saveto:
        fig.savefig(saveto, bbox_inches='tight')
    return ax


def plot_umap_genes(adata, gene_dict, figdir=None,
                    tag='',
                    cmap_val=None, cmap_val_name='diy_RdPu', # 'diy_Blues'
                    save_pdf=True,
                    save_png=True,
                    min_cells=5, 
                    foo_plt = sc.pl.umap,
                    qt = 0.99):
    '''
    gene_dict: a dict mapping `gene-id` to `gene-name`(used as title)
    '''
    if figdir is not None:
        sc.settings.figdir = figdir 
    fx.check_dirs(sc.settings.figdir)
    
    if cmap_val_name.startswith('diy_'):
        _cmap = cmap_val_name.strip('diy_')
        cmap_val = fx.diy_cmap_grey_bg(_cmap)
    else:
        cmap_val = fx.get_colors(cmap_val_name, n=100, to_hex=True)
    
    # plot UMAP highlighting genes (markers)
    for gid, gn in gene_dict.items():
        expressed_genes = adata.var_names if adata.raw is None else adata.raw.var_names
        if gid not in expressed_genes:
            print(f'Not any cells expressed gene {gn}({gid}), skipped')
            continue
        if adata.raw is None:
            gexpr = adata[:, gid].X.flatten()
        else:
            gexpr = adata.raw[:, gid].X.flatten()
        nnz = gexpr[gexpr > 0]
        if len(nnz) <= min_cells:
            print(f'No more than {min_cells} cells expressed gene {gn}({gid}), skipped')
            continue
        _vmax = np.quantile(nnz, qt)
        
        if save_png:
            foo_plt(adata, color=gid, title=gn, cmap=cmap_val, 
               vmax=_vmax,
               save=f'_{tag}_{cmap_val_name}_{gn}_{gid}.png' 
               )
        if save_pdf:
            foo_plt(adata, color=gid, title=gn, cmap=cmap_val, 
               vmax=_vmax,
               save=f'_{tag}_{cmap_val_name}_{gn}_{gid}.pdf' #+ '-sub.png'
               )




















def plotSimiMatrix(sims, figsize = (6, 6), ax=None,
                   cmap=None, save=None, figdir = None, **kw):
    cmap = ['RdBu_r', 'vlag', 'magma_r'][-1]
    
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    sns.clustermap(sims, cmap=cmap, ax=ax,
                   xticklabels=False,
                   linewidths=.5,
#                   row_colors = cl_color_match,
                   row_linkage=None, col_linkage=None,
                   **kw)
#    plt.savefig(figdir / f'corr_{cmap}_top{n_top}mk0.pdf', bbox_inches='tight')
    if save:
        ax.figure.savefig(figdir / f'corr_{cmap}_temp.pdf', bbox_inches='tight')  

    


# In[]
''' helper functions     
'''

def order_contin_df(df, axis=1):
    '''
    re-order contingency matrix (for dataframe)
    axis = 0: re-order the columns
    axis = 1: re-order the rows
    '''
    order = np.argsort(np.argmax(df.values, axis=axis))
    if axis == 1:
        print('Re-order the rows')
        return df.iloc[order, :]
    else:
        print('Re-order the columns')
        return df.iloc[:, order]

def order_contin_matrix(mat, axis = 1):
    '''
    re-order contingency matrix
    axis = 0: re-order the columns
    axis = 1: re-order the rows
    '''
    order = np.argsort(np.argmax(mat, axis=axis))
    if axis == 1:
        print('Re-order the rows')
        return mat[order, :]
    else:
        print('Re-order the columns')
        return mat[:, order]
















