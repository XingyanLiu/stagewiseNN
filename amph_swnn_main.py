# -*- coding: utf-8 -*-
"""
Created on Sun Jan 10 13:48:19 2021

@author: xyliu
"""


import os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
import scanpy as sc

WORKDIR = Path(r'D:\Users\xyliu\003')
os.chdir(WORKDIR)

from PipeHandler import PipeHandler 
import Build as B
import funx as fx
import funx_amph as fxa
from StagewiseNN import stagewise_knn

datadir = fxa.DATADIR / 'afterQC_formal' 
genedir = fxa.DATADIR / 'genes'
resdir = fxa.DATADIR / 'res-scanpy' / 'stg_nn-20210110'
figdir = resdir / 'figs'
fx.check_dirs(figdir)

sc.settings.figdir = figdir
sc.set_figure_params(fontsize=14)

# In[]
''' Whole data '''
Adata = sc.read_h5ad(datadir / 'merged_afterQC.h5ad')
Adata.obs['stage_name'].value_counts(sort=False)
Adata

# In[]
''' change stage names; merge B0 with B1
'''
stage_id_name = fxa.StageNameDict
def foo_change0(x: str):
    for old, new in stage_id_name.items():
        if old in x:
            return x.replace(old, new)
#            break
    return x
        
change_dict = {
        'B_1': 'B_0',
        'B_2': 'B_1',
        'B_3': 'B_2',
        }
foo_change1 = lambda x: change_dict.get(x, x)

#tmp = 'E6_1'
#tmp1 = foo_change0(tmp)
#foo_change1(tmp1)

Adata.obs['refined_group'].value_counts()
cluster_labels = Adata.obs['refined_group']
cluster_labels = fx.change_names(cluster_labels, foo_change0)
cluster_labels = fx.change_names(cluster_labels, foo_change1)
pd.value_counts(cluster_labels)
Adata.obs['refined_group'] = cluster_labels
Adata.obs['lineage'].value_counts()
Adata.obs['lineage'] = fx.change_names(Adata.obs['lineage'], foo_change1)

### average expression by counts
avg_cnts = B.GroupMean(Adata, 'refined_group')
avg_cnts = avg_cnts.apply(np.log1p)
avg_cnts.to_csv(resdir / f'avg_expr_all-logAvg.csv', index=True, header=True)

### normalize by total sums (NOT necessary)
med = avg_cnts.sum().median()
avg_cnts_scaled = avg_cnts.apply(lambda x: med * x / np.sum(x), axis=0)
avg_cnts_scaled.to_csv(resdir / f'avg_expr_all.csv', index=True, header=True)

### expression proportions
avg_prop = B.GroupMean(Adata, 'refined_group', binary=True)
avg_prop.to_csv(resdir / f'expr_prop_all.csv', index=True, header=True)




# In[]

hvg_freq = pd.read_csv(fxa.DATADIR / 'genes' / 'hvg_freq_scanpy-later.csv', index_col=0, header=None).iloc[:, 0]
used_genes1 = hvg_freq[hvg_freq >= 3].index.tolist()
hvg_freq.value_counts()

candi_markers = fx.load_namelist(genedir / 'candidate_marker_ids.csv', header='infer')

fx.VennPlot([candi_markers, used_genes1])
used_genes = list(set(candi_markers + used_genes1))


# In[]


adata = Adata.copy()

B.normalize_reverse(adata, counts_level=None)
adata = B.set_hvgs(adata, used_genes, slim=True)
adata
adata.raw.X.data[:5]
B.GroupZscoreAdt(adata, groupby='stage_primer')
#sc.tl.pca(adata, n_comps=50, )


# In[]

X = adata.X





