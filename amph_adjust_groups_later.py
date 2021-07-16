# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 18:38:36 2020

@author: xyliu
"""

import os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
import scanpy as sc
import matplotlib.pyplot as plt

WORKDIR = Path(r'D:\Users\xyliu\003')
os.chdir(WORKDIR)

from PipeHandler import PipeHandler 
import Build as B
import funx as fx
import funx_amph as fxa
#from StagewiseNN import stagewise_knn

# In[]
'''     setting result directory
'''

resdir = fxa.DATADIR / 'res-adjust' / '0702-later'
figdir = resdir / 'figs'
fx.check_dirs(figdir)

#genedir = fxa.DATADIR / 'genes'
sc.set_figure_params(fontsize=14, dpi_save=180)


# In[]
'''     load merged `Adata`
'''

dirmg = fxa.DATADIR / 'res-scanpy' / 'stg_nn-hvgs-later'
Adata = sc.read_h5ad(dirmg / 'final-sepPCA50.h5ad')
Adata
sc.pl.umap(Adata, color='stage', palette = 'Spectral')

ph = PipeHandler(Adata, name = 'sepPCA50', resdir = resdir)

# In[]
'''     load stage adatas into a dict
'''
dirsep = fxa.DATADIR / 'res-scanpy' / 'stgwise-0602'

Stages = fxa.Stages
Snames = fxa.sampleNames(Stages)

adatas = {stg: sc.read_h5ad(dirsep / f'{sn}.h5ad') for \
          stg, sn in zip(Stages, Snames)}

#adatas['E6'] = sc.read_h5ad(dirsep / f'E6_polyT.h5ad') 

# In[]
'''     merge geoups (manually deal with each cluster)
'''
def plotMyUmap(adt, color='leiden', stage=None, 
               palette = fxa.CMAP_CLASS, 
               tag='',
               legend_loc = 'on data', # 'right margin'
               save=True,
               figtype='.pdf'):
    if stage is None:
        stage = adt.obs['stage'][0]
    tt = f'{stage} ({adt.shape[0]} cells)'
    fn = f'_{stage}_{color}{tag}{figtype}' if save else None
    sc.pl.umap(adt, color=color, 
               palette = palette, 
               title = tt,
               legend_loc = legend_loc, 
               save=fn)    

def mergeGroupsAdatas(adatas, decisions, key = 'leiden', 
                      new_key = None, rename=True, 
                      plot=True, figtype='.pdf'):
    '''
    '''
    for stage in decisions.keys():
        group_lists = decisions[stage]
        adt = adatas[stage]
        B.MergeGroups(adt, key=key, group_lists=group_lists, 
                      new_key=new_key, rename=rename, copy = False) 
        if plot:
            #sc.pl.umap(adt, color=new_key, palette = fxa.CMAP_CLASS,)
            plotMyUmap(adt, color=new_key, stage=stage, figtype=figtype)



decisions = {
#        'E6': [['1', '3'],],
        'E7': [['0', '1', '3'],],
        'E8': [['0', '1', '2'], ['3', '4', '5'], ['7', '8', '9']],
        'E10': [['2', '3', '7'], ],
        'E11': [['1', '10', '11']],
        'E12': [['1', '3', '11']],
        'E13': [['0', '11']],
#        'E14': [['9', '10']]
        }

#stage = 'E6'
#group_lists = [['1', '3'],]
key = 'leiden'
new_key = f'{key}_new'

mergeGroupsAdatas(adatas, decisions, key, new_key, rename=True)

# In[]
''' Merge by KNN classifier
'''
adt = adatas['E7']
is_query = adt.obs[new_key] == '5'
is_query.value_counts()

fxa.labelByKNN(adt, is_query, new_key, inplace=True, rename_cat=True,)
adt.obs[new_key].value_counts()
plotMyUmap(adt, color=new_key, )

# In[]
'''     save results (stagewise)
'''


# In[]
''' plot stages (uaing newly adjusted labels)
'''

key_colr = 'leiden_new'





# In[]
'''     label `Adata` by newly adjusted labels
'''
def _take_labels(adt, key=new_key, key0 = 'leiden'):
    return adt.obs[key] if key in adt.obs.columns else adt.obs[key0]

new_labels = pd.concat([_take_labels(adatas[stg]) for stg in adatas.keys()])

ph.adata.obs[new_key] = new_labels
#ph.adata.obs[new_key].value_counts()

new_key1 = 'stage_leiden_new'
B.cross_labeling(ph.adata.obs, ['stage', new_key], new_key=new_key1, inplace=True)




# In[]
'''     reconstruction of the tree
'''
import seaborn as sns
#import networkx as nx
import amph_graph2tree as g2t

dirmg = resdir
connect0 = sparse.load_npz(dirmg / 'formal_connect0.npz')
#sparse.save_npz(resdir / 'connect0.npz', connect0)

df = ph.adata.obs
group_lbs = df[new_key1].values
stage_lbs = df['stage'].values
stage_ord = list(fxa.Stages)
key_refined = 'stg_groups_new'

conn_upper = sparse.triu(connect0).tocsc()
adj_max = g2t.max_connection(conn_upper)

ignore_pa = ['E5_0', 'E5_1', 'E5_2', 'E13_4'] #['E6_2', 'E6_0']
###================ main ====================
edgedf, new_group_lbs = g2t.adaptiveTree(adj_max, group_lbs, 
                                         stage_ord=stage_ord,
                                         ignore_pa=ignore_pa)
df[key_refined] = new_group_lbs
edgedf.prop.hist()
edgedf.head(20)
edgedf.tail(20)
#edgedf.iloc[86, 1] = 'E13_4' # E14_9
#edgedf.iloc[82, 1] = 'E13_6' # E14_16
##edgedf.iloc[16, 1] = 'E8_3' # 'E9_7'
#####################[ organize results ]####################
stg_grp_dict = g2t.make_stage_group_dict(new_group_lbs)
cells_per_group = pd.value_counts(new_group_lbs)
all_groups = pd.unique(new_group_lbs)



# In[]
'''
#################[ voting matrices ]###############
'''
voting_props = g2t.agg_group_edge_props( \
        adj_max, new_group_lbs, groups0=all_groups)
voting_props = voting_props.T # tranposed to  son-by-parent


sc.set_figure_params(fontsize=9)
key_foo = lambda x: int(x.split('_')[1])
stage_ord_use = stage_ord
for i, pa in enumerate(stage_ord_use[: -1]):
    son = stage_ord_use[i + 1]
    stg_pair = (son, pa) # son, pa
    print(stg_pair)
    son_grps = sorted(stg_grp_dict[stg_pair[0]], key=key_foo)
    pa_grps = sorted(stg_grp_dict[stg_pair[1]], key=key_foo)
    
    voting_props_sub = fxa.order_contin_df(voting_props.loc[son_grps, pa_grps])
    
    saveto = figdir / 'voting_props_{}-{}.png'.format(*stg_pair)
    ax = fxa.plot_matrix(voting_props_sub, saveto=saveto)


# In[]

'''
####################[ make tree ]#####################
'''
gs1 = edgedf.pivot(index='parent', columns='node', values='prop')
gs1 = gs1.reindex(all_groups, all_groups)
gs1.fillna(0, inplace=True)

G = fxa.TreeByStages(gs1, stages = stage_ord, cross_stg=False)
fxa.plot_tree_by_stages(G, )


#########################################
'''     save init (voting proportions)
'''
treedir = resdir / 'treedata'
fx.check_dirs(treedir)

edgedf.head()
edgedf = edgedf.set_index('node', drop=False)
edgedf['n_cells'] = cells_per_group
edgedf.to_csv(treedir / f'tree_struct_init.csv', index=False, header=True)
voting_props.to_csv(treedir / 'voting_props.csv', index=True, header=True)


# In[]

''' prepare for ggTree
'''


df_tree = fxa.dfFromG(G, default_root = 'root')
df_tree.head()

df_tree.to_csv(treedir/ f'tree_struct.csv', index=False, header=True)

#
### compute group average for all genes
groupby = key_refined
adata = ph.adata
props_all = B.GroupMean(adata, groupby, use_raw=True, binary=True,)
means_all = B.GroupMean(adata, groupby, use_raw=True,)
props_all.to_csv(treedir / f'expr_prop_all.csv', index=True, header=True)
means_all.to_csv(treedir / f'avg_expr_all.csv', index=True, header=True)


# In[]
'''     label stagewise adatas with tree-refined groups
'''
key_refined #= 'stg_groups_new'
key_refined_int = 'refined_group'

sc.set_figure_params(fontsize=14)
for stg in adatas.keys():
    adt = adatas[stg]
    adt.obs[key_refined] = ph.adata.obs[key_refined]
    adt.obs[key_refined_int] = adt.obs[key_refined].apply(\
            lambda x: x.split('_')[1])
    print(adt.obs[[key_refined_int, key_refined]].head())
    _new_key = new_key if new_key in adt.obs.columns else 'leiden'
    plotMyUmap(adt, color=_new_key, )
    plotMyUmap(adt, color=key_refined_int)
    plotMyUmap(adt, color=_new_key, figtype='.png')
    plotMyUmap(adt, color=key_refined_int, figtype='.png')    
    plotMyUmap(adt, color=key_refined_int, figtype='.png', tag='_tab20', palette='tab20')    
    print()
    

#tmp = adt.obs

# In[]
'''     make a column conatining parent index (barcode)
'''
##### these codes have been run when constructing tree
#connect0 = sparse.load_npz(dirmg / 'formal_connect0.npz')
#df = ph.adata.obs
#conn_upper = sparse.triu(connect0).tocsc()
#adj_max = g2t.max_connection(conn_upper)

sons = df.index.take(adj_max.col)   #len(sons.unique())
parents = df.index.take(adj_max.row)

pa_bcd = pd.Series(parents, index=sons)
pa_bcd.head()

key_pa = 'parent_bcd'
df[key_pa] = pa_bcd

#ph.save()
#ph.name

# In[]
'''     label stagewise adatas with parent barcodes and groups
'''
key_pa #= 'parent_bcd'
key_pa_lb = 'parent_group'
key_pa_take = key_refined

for stg in adatas.keys():
    adt = adatas[stg]
    adt.obs[key_pa] = ph.adata.obs[key_pa]
    adt.obs[key_pa_lb] = ph.adata.obs[key_pa_take][adt.obs[key_pa]].values.astype(str)
    
    print(adt.obs[['stage', key_pa_take]].head())
    _new_key = new_key if new_key in adt.obs.columns else 'leiden'
#    plotMyUmap(adt, color=_new_key, )
#    plotMyUmap(adt, color=_new_key, figtype='.png')
    plotMyUmap(adt, color=key_pa_lb, legend_loc='right margin', tag='tab20',
               palette='tab20')
    plotMyUmap(adt, color=key_pa_lb, legend_loc='right margin', tag='tab20',
               figtype='.png', palette='tab20')    
    print()
#sc.pl.scatter
#resdir
#tmp = ph.adata.obs[key_pa_take][adt.obs[key_pa]].values.astype(str)

#adt = adatas['E8']    
#ids_modify = adt.obs[['refined_group', 'parent_group']].apply(\
#                    lambda x: x[0] == '2' and x[1] == 'E7_2', axis=1)
#adt.obs.loc[ids_modify, 'refined_group'] = '3'
#adt.obs.loc[ids_modify, 'leiden_new'] = '3'
##ids_modify.sum()
#sc.pl.umap(adt, color='refined_group')
#sc.pl.umap(adt, color='leiden_new')

# In[]
''' save Stagewise adatas (and metadata / labels)
'''
#fxa.sampleNames(['E8'])

for stg, adt in adatas.items():
    sn = fxa.sampleNames([stg])[0]
    print(stg, ':', sn)
    adt.write(resdir / f'{sn}.h5ad')
    adt.obs.to_csv(resdir / f'{sn}_metadata.csv', index=True, header=True)

#ph.save()
#ph.save_dims('umap', ftype='csv')

#_stg = 'E6'
#
#B.save_embeddings(adatas['E6'], resdir, key='X_umap', tail=_stg)  
#adatas['E6']

# In[]





tmp = fx.get_colors('tab20', n=8)
fx.view_color_map(tmp)










# In[]
''' highlight the clusters in some specific stage
'''

def label_mask_others(df, key, keeps, key_new='tmp', lb_masked='others', copy=False):
    
    if copy:
        df = df.copy()
    df[key_new] = df[key].apply(lambda x: x if x in keeps else lb_masked)
    
    return df if copy else None

stg = 'E6'
key = 'stage_stg_leiden'
key_new = f'{stg}_clusters'
keeps = ph.adata.obs[key][ph.adata.obs[key].str.startswith(stg)].unique()
#keeps = ['E6_' + str(j) for j in [0, 1, 2, 3]]
label_mask_others(ph.adata.obs, key, keeps, key_new)
ph.vis(color = key_new, groups = keeps, save=f'_{key_new}')


# In[]































