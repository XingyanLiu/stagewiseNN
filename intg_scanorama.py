# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 07:51:51 2019

@author: xyliu
"""

'''

-------------------------------------

# List of data sets (matrices of cells-by-genes):
datasets = [ list of scipy.sparse.csr_matrix or numpy.ndarray ]
# List of gene lists:
genes_list = [ list of list of string ]

>>> import scanorama

# Integration.
>>> integrated, genes = scanorama.integrate(datasets, genes_list)

# Batch correction.
>>> corrected, genes = scanorama.correct(datasets, genes_list)

# Integration and batch correction.
>>> integrated, corrected, genes = scanorama.correct(datasets, genes_list, return_dimred=True)


-------------------------------------
There are also wrappers that make it easy to use Scanorama with scanpy's AnnData object:

# List of data sets:
adatas = [ list of scanpy.AnnData ]

# Integration.
>>> integrated = scanorama.integrate_scanpy(adatas)

# Batch correction.
>>> corrected = scanorama.correct_scanpy(adatas)

# Integration and batch correction.
>>> integrated, corrected = scanorama.correct_scanpy(adatas, return_dimred=True)

'''
import os
from pathlib import Path
import scanpy as sc
import scanorama
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.decomposition import PCA
import umap
import matplotlib.pyplot as plt

WORKDIR = Path(r'D:\Users\xyliu\003')
os.chdir(WORKDIR)

from PipeHandler import PipeHandler 
import funx as fx

DATADIR_main = Path(r'E:\Users\xyliu\data')

Proj = ['amphioxus', 'hippocampus', 'RNF200'][-1]
Tp = ['adults', 'embryos', 'AE', 'NS', 'BF', 'KO'][-1]  
DATADIR = DATADIR_main.joinpath(Proj, Tp)

#ph = PipeHandler(sc.datasets.paul15())
###############[ for series of samples ]##############

Snames = {'adults': ['T%d'%i for i in range(1, 15)], 
          'embryos': ['E%d_02'%i for i in range(1, 16)], 
          'AE': ['AE'], 
          'NS': ['D%d'%i for i in range(1, 7)], 
          'BF': 'B2 B4 B5 B7 B8 F8'.split(),
          'KO': ['KO_1', 'ref']}
#Snames['adults'].insert(11, 'T11_mtx')
Snames['adults'].remove('T6')
snames = Snames[Tp]

resdir0 = DATADIR / ('%s_%s'%('afterQC', 'manu_auto1'))
resdir0 = DATADIR / ('%s_%s'%('afterQC', 'euclidean'))
resdir = resdir0 / 'integrate'
fx.check_dir(resdir)

n_top_de = 50
used_genes = []
for sn in snames:
   
    degs = pd.read_csv('%s/%s_marker_names.csv' % (resdir0, sn))
    top_degs = pd.unique(degs.head(n_top_de).values.flatten()).tolist()
    print('%s: using top %d DE genes, where %d unique genes'%(sn, n_top_de, len(top_degs)))
    used_genes += top_degs
used_genes = pd.unique(used_genes)
print('using %d unique genes'%(len(used_genes)))

#label_period = []
label_batch = []
adatas = []
for i, sn in enumerate(snames):
    adata = sc.read_h5ad('%s/%s_afterQC.h5ad'%(resdir0, sn))
    print('%s\n\tshape of the raw data: '% sn, adata.X.shape)
    gene_indicator = [g in used_genes for g in adata.var_names]
    adata = adata[:, gene_indicator]
    print('\tshape of the filterd data: ', adata.X.shape)
    adatas.append(adata)
    label_batch += [i] *  adata.X.shape[0]    
pd.Series(label_batch).to_csv(resdir / 'meta_batch.csv', header=False, index=False)


print('using %d unique genes for integration'%(len(used_genes)))


alpha = 0.3
dimred = 80
union = False
integrated = scanorama.integrate_scanpy(adatas, union=union, alpha=alpha, 
                                        dimred=dimred)

'''
# ----------- Concatenate reduced matrices -------------
'''
catted = np.vstack(integrated)

Tag = 'scanorama%s_unionTop%s_dim%d'%(alpha, n_top_de, dimred)
print(Tag)
np.save(resdir / ('intg_%s.npy'%(Tag)), catted)
#pd.Series(label_period).to_csv('intg_period.csv', header=False, index=False)
# ---------Reload---------
catted = np.load(resdir / ('intg_%s.npy'%(Tag)))
label_batch = pd.read_csv(resdir / 'meta_batch.csv', header=None).iloc[:, 0]

#ump = umap.UMAP(n_neighbors=30).fit(catted)
#embedding = ump.transform(catted)
nneigh = 10
print('Computing UMAP embedding with %d neighbors'%nneigh)
embedding = umap.UMAP(n_neighbors=nneigh).fit_transform(catted)

#DR = PCA(2)
#embedding = DR.fit_transform(catted)
df_embed = pd.DataFrame({'umap1': embedding[:, 0],
                         'umap2': embedding[:, 1],
                         'batch': label_batch})

df_embed['batch'] = pd.Categorical(df_embed['batch'])
df_embed.to_csv(resdir / 'embed_umap{}.csv'.format(nneigh), header=True, index=None)

'''
#---VIsualization---
'''
name = 'intg_umap{}_{}'.format(nneigh, Tag)
fx.ButtomView(embedding, label_batch, key='umap', sz=0.2, disc=True, 
              save=resdir, name=name)

#from matplotlib import cm
#cmp = cm.get_cmap('RdYlBu', 15)
##cmp = cm.tab20
#fig, ax = plt.subplots()
#for i in df_embed['batch'].unique():
#    ax.scatter('umap1', 'umap2', label=snames[i], color=cmp(i, alpha=0.6), 
#               s=1.5,
#               data = df_embed[df_embed['batch'] == i])
#ax.legend()
#fig.savefig(resdir / 'intg_umap{}_{}_.png'.format(nneigh, Tag))


#---VIsualization---
#plt.scatter(embedding[:, 0], embedding[:, 1], c = label_period)


    




