# -*- coding: utf-8 -*-
"""
Created on Wed May 20 20:45:15 2020

@author: xyliu
"""


import os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
import scanpy as sc

import matplotlib as mpl
#mpl.use('Agg')
#%matplotlib inline
print(mpl.get_backend()) # module://ipykernel.pylab.backend_inline
print(mpl.matplotlib_fname())
import matplotlib.pyplot as plt


WORKDIR = Path(r'D:\Users\xyliu\003')
os.chdir(WORKDIR)

from PipeHandler import PipeHandler 
import funx as fx
import funx_amph as fxa
import Build as B
import ParamSetting as ps


# In[]
'''      Loading stagewise data into a list
=====================================================
'''

datadir = fxa.DATADIR / 'afterQC_formal'
resdir = datadir / 'Merged'
figdir = resdir / 'figs'
fx.check_dirs(figdir)

# loading adatas into a list
tail = '_afterQC'
Stages = fxa.Stages0
Snames = fxa.sampleNames()
adatas = [sc.read_h5ad(datadir / f'{sn}{tail}.h5ad') for sn in Snames]


 # In[] 
''' HVGs
'''
# HVGs from Seurat
hvgdir = datadir
hvgdir = fxa.DATADIR / 'res-srt'
hvgList = [pd.read_csv(hvgdir / f'{sn}_hvgs.csv', header=None).iloc[:, 0] for sn in Snames]
hvgFreq = pd.concat(hvgList, axis = 0).value_counts()
hvgFreq.value_counts()#hvgFreq.hist(bins=15)
hvgFreq.index = [x.replace('-', '_') for x in hvgFreq.index]
hvgFreq.head()
hvgFreq.to_csv(datadir / 'hvg_freq_srt.csv', index=True, header=False)

hvgs = hvgFreq[hvgFreq >= 5].index.tolist()
hvgs = [g.replace('-', '_') for g in hvgs]


# HVGs from scanpy
hvgs = fx.load_namelist(datadir / 'merged_hvgs.csv')




# In[]
''' add UMAP PCA
'''
drdir = fxa.DATADIR / 'res-srt'

for i, adt in enumerate(adatas):
    sn = Snames[i]
    X_umap = pd.read_csv(drdir / f'umap_{sn}.csv')
    X_pca = pd.read_csv(drdir / f'pca_{sn}.csv')
    adt.obsm['X_umap'] = X_umap.values
    adt.obsm['X_pca'] = X_pca.values




# In[]
''' Normalization
''' 

for adt in adatas:
    
    sc.pp.log1p(adt)
    sc.pp.normalize_total(adt, target_sum=300)


# In[]

''' merge labels (from srt results)
    adatas are NOT involved 
'''
key_lb = 'tree_ident'
lb_list = []
fp_main = fxa.DATADIR / 'res-srt'
for i, sn in enumerate(Snames):

    _fpath = fp_main / f'{sn}_metadata_srt.csv'
    meta_srt = pd.read_csv(_fpath, index_col=0)
#    srt_lbs.head()
    if len(meta_srt['seurat_clusters'].unique()) == 1:
        meta_srt['tree.ident'] = '1'
        
    srt_lbs = meta_srt['tree.ident'] #if 'tree.ident' in meta_srt.columns \
#        else meta_srt['seurat_clusters']
    print(srt_lbs.name)
    lb_list.append(srt_lbs.astype(str))

cated_lbs = pd.concat(lb_list)
cated_lbs.to_csv(fp_main / f'_cated_{key_lb}.csv', index=True, header=True)

# In[]
'''
# reordered labels from Seurat function `BuildClusterTree()`
'''
key_reorder = 'tree_ident'
sn = 'E15'
for i, sn in enumerate(Snames):
    adt = adatas[i]

    _fpath = fxa.DATADIR / 'res-srt' / f'{sn}_metadata_srt.csv'
    meta_srt = pd.read_csv(_fpath, index_col=0)
#    srt_lbs.head()
    if len(meta_srt['seurat_clusters'].unique()) == 1:
        adt.obs[key_reorder] = '1'
#        adt.write(datadir / f'{sn}.h5ad')
        continue
    srt_lbs = meta_srt['tree.ident'] if 'tree.ident' in meta_srt.columns \
        else meta_srt['seurat_clusters']
    adt.obs[key_reorder] = srt_lbs.astype(str)
#    adt.write(datadir / f'{sn}.h5ad')
    
#    adt.obs.head()
    
    
# In[]
    '''     Tree Plot
    '''
import seaborn as sns
## group similarity
keys = key_reorder
metric = ['spearman', 'correlation'][1]
gs, means = B.GroupSimilarityMean(adatas, keys, use_genes=hvgs,
                        metric=metric, output_means=True,
                        use_raw=True, tags=Stages,)

sns.heatmap(gs, )

gs.to_csv(resdir.parent / 'merged' / 'group_sim_cross_stages.csv',)
# construct tree
G = fxa.TreeByStages(gs, )
## plot graph
figsize = (12, 15)
fn_fig = f'lineage_tree_{metric}.pdf'
fxa.plot_tree_by_stages(G, node_size=500, figsize=figsize,
                        save=figdir / fn_fig, show=False)

    
# In[]

'''
   Annotations

''' 
REFDF = pd.read_csv(fxa.DATADIR / 'Summ_markerIds.csv')
REFDF.columns
CellTypes = REFDF['cell type'].unique()
geneIDs = REFDF['ID'].tolist()

cell_types = CellTypes[:5]
gene_ids0 = REFDF.loc[REFDF['cell type'].apply(lambda x: x in cell_types), 'ID'].tolist()


adt = adatas[-1]

gene_ids = [g for g in gene_ids0 if g in adt.var_names]
sc.pl.dotplot(adt, gene_ids, groupby='tree_ident', var_group_labels='e')
    
    


# In[]
'''     stagewise downsampling
 =================================================
'''
adatas

adatas_sub = []
#i = 0
#adt = adatas[i]
for i, adt in enumerate(adatas):
    n = 500 if i < 5 else i * 200
    adt1 = sc.pp.subsample(adt, n_obs=n, copy=True)
    adatas_sub.append(adt1)

adatas_sub

obs_keys = 'stage primer n_genes n_counts tree_ident'.split()
Adata_sub = fx.merge_adatas(adatas_sub, union=True, obs_keys=obs_keys)

Adata_sub.write(fxa.DATADIR / 'afterQC_formal/merged_subsampled.h5ad')
#fx.saveNamedMtx(Adata_sub, fxa.DATADIR / 'merged_subsampled')






















