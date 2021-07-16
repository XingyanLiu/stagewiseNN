# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 19:51:11 2020

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

DATADIR = Path(r'E:\Users\xyliu\data003\amph')


# In[]
'''      Loading stagewise data into a list
=====================================================
'''
datadir0 = Path(r'E:\Users\xyliu\data003\amph\res-raw_ct100')
datadir = datadir0 / 'stagewise_noreg_log1first'
resdir = datadir / 'merged'
figdir = datadir / 'figs_merged'
fx.check_dirs(figdir)

# loading adatas into a list
tail = '_afterQC'
Stages = fxa.Stages0
adatas = [sc.read_h5ad(datadir / f'{sn}{tail}.h5ad') for sn in Stages]

# In[]
'''
# merge adatas into a single AnnData
'''
union= True
obs_keys = 'stage primer n_genes n_counts tree_ident'.split()[:-1]
Adata = fx.merge_adatas(adatas, union=union, obs_keys=obs_keys)
Adata.write(resdir / f'merged{tail}.h5ad')

fx.saveNamedMtx(Adata, resdir / f'merged{tail}_mtx', field='integer')

# In[]
#Adata.obs['n_genes'] = (Adata.X > 0).sum(axis=1).A1
summ_n_genes = Adata.obs.groupby('stage')['n_genes'].describe()
summ_n_genes.to_csv(resdir / 'summ_n_genes.csv')
summ_n_genes
summ_n_counts = Adata.obs.groupby('stage')['n_counts'].describe()
summ_n_counts.to_csv(resdir / 'summ_n_counts.csv')
summ_n_counts

sc.set_figure_params(fontsize=14)

_tt = 'number of genes per cell'
summ_n_genes.plot.bar(y=['min', '50%'], figsize=(6, 3), title=_tt, fontsize=13)
plt.savefig(figdir / f'{_tt}.pdf', bbox_inches='tight')

_tt1 = 'number of counts per cell'
summ_n_counts.plot.bar(y=['min', '50%'], figsize=(6, 3), title=_tt1, fontsize=13)
plt.savefig(figdir / f'{_tt1}.pdf', bbox_inches='tight')




# In[]
'''     Re-loading the marged data
==================================================================
'''


#datadir =  DATADIR / 'res-raw_ct100/stagewise_noreg_log1first/merged'
datadir =  fxa.DATADIR / 'afterQC_formal'

resdir = datadir
figdir = datadir.parent / 'figs_merged'
fx.check_dirs(figdir)

Adata = sc.read_h5ad(datadir / 'merged_afterQC.h5ad')
print(Adata)
hvg_freq = pd.read_csv(datadir/ 'merged_hvg_freq_scanpy.csv', 
                       index_col=0, header=None).iloc[:, 0]
hvg_freq.head()
hvgs = hvg_freq[hvg_freq >= 6].index.tolist()

#hvgs = Adata.var['hvg_freq'][Adata.var['hvg_freq'] >= 6].index.tolist()
#sc.pl.scatter(Adata, x='n_counts', y='n_genes',# palette='RdYlBu',color_map='Spectral',
#              color='stage', save='_gvc')

# In[]
#name = 'lognorm'
name = 'merged0'

# re-analyze
ph = PipeHandler(Adata, name=name, resdir=resdir)
batch_key = 'primer'
metric = 'cosine'
n_pcs = 50
nneigh = 20


ph.NormLog_rev(counts_level=None)


#sc.pp.normalize_total(ph.adata, )
#ph.adata.raw = ph.adata

# In[]
# dowstream analysis
sc.set_figure_params(fontsize=13)
ph.SetHVGs(hvgs, slim=True)
ph.save_hvgs()
ph.scale(do_regress=False, groupby=batch_key)
ph.PCA(n_comps=80, plot=True)
ph.Neighbors(nneigh = nneigh, do_pca = True, use_rep = 'X_pca',
           n_pcs = n_pcs, metric=metric)#'cosine') # 'euclidean' by default

ph.Umap(min_dist = 0.5)
ph.save_dims('umap')
tt = f'total of {ph.adata.shape[0]} cells'
sc.pl.umap(ph.adata, color='stage', palette='Spectral', 
           legend_loc='on data', 
           title=tt, save=f'_all_stage_{name}.pdf')
sc.pl.umap(ph.adata, color='stage', palette='Spectral', 
#           legend_loc='on data', 
           title=tt, save=f'_all_stage_{name}_legend.pdf')
ph.vis(color='primer', key='umap', title=tt, legend_loc='on data', 
       save=f'_all_primer.pdf')


ph.DE('stage', save=True, plot=True,)

# just overview the expression values to decide the max-cutoff
tmp = ph.adata.raw.X.data
plt.hist(tmp, bins=20) 
ph.dot_de(groupby='stage', vmax=2, dendrogram=False)

# formally plot and save
ph.dot_de(save=f'_by_stage_{name}', groupby='stage', vmax=2)#,standard_scale='var')
ph.dot_de(save=f'_by_stage_std_{name}', groupby='stage',standard_scale='var')
ph.save_markers()

# In[]
''' ============= TSNE ============
'''
ph.TSNE()
ph.save_dims('tsne')
tt = f'total of {ph.adata.shape[0]} cells'
ph.vis(color='stage', key='tsne', title=tt, legend_loc='on data', 
       save=f'_all_stage_{name}.pdf')
ph.vis(color='stage', key='tsne', title=tt, #legend_loc='on data', 
       save=f'_all_stage_{name}_legend.pdf')
ph.vis(color='primer', key='tsne', title=tt, legend_loc='on data', 
       save=f'_all_primer_{name}.pdf')

# In[]
# clustering
ph.Clustering(res = 1, do_pca = True, n_pcs = n_pcs, 
              redo=False, pl_umap = True, keep_results=True)

ph.save()

# In[]
''' inspect some certain markers
'''
# PGCs
#Tbx2 (bf_00016583)
#Nanos (bf_00004829)

genes_pgc = {
#        'bf_00016583': 'Tbx2',
#        'bf_00004829': 'Nanos',
        'bf_00008183': 'LRIG3',
        'bf_00020985': 'CHRO',
        }
cmap_g = fx.diy_cmap_grey_bg()

for gid, gn in genes_pgc.items():
    gexpr = ph.adata.raw[:, gid].X.flatten()
    _vmax = np.quantile(gexpr[gexpr > 0], 0.99)
    
    ph.vis(gid, title=gn, cmap=cmap_g, 
           vmax=_vmax,
           save=f'_{gn}')







# In[]
kadd = 'stage_tree_ident'
sc.tl.dendrogram(ph.adata, groupby=kadd, use_rep='X', var_names=hvgs)
ax = sc.pl.dendrogram(ph.adata, groupby=kadd, 
                      figsize=(12, 4), save='_all.png')


sc.pl.correlation_matrix(ph.adata, groupby=kadd, cmap='RdYlBu_r', 
                         save='_all_stages.pdf')



# In[]
#'''     stagewise downsampling
# =================================================
#'''
#adatas
#
#adatas_sub = []
##i = 0
##adt = adatas[i]
#for i, adt in enumerate(adatas):
#    n = 500 if i < 5 else i * 200
#    adt1 = sc.pp.subsample(adt, n_obs=n, copy=True)
#    adatas_sub.append(adt1)
#
#adatas_sub
#
#obs_keys = 'stage primer n_genes n_counts tree_ident'.split()
#Adata_sub = fx.merge_adatas(adatas_sub, union=True, obs_keys=obs_keys)
#
#Adata_sub.write(fxa.DATADIR / 'afterQC_formal/merged_subsampled.h5ad')
#fx.saveNamedMtx(Adata_sub, fxa.DATADIR / 'merged_subsampled')


