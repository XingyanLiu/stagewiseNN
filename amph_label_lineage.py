# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 19:16:53 2020

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

import amph_graph2tree as g2t




datadir = fxa.DATADIR / 'res-adjust' / '0702-later'


resdir = datadir
figdir = resdir / 'figs'
fx.check_dirs(figdir)


sc.set_figure_params(fontsize=14)

# In[]
''' loading merged Adata
'''
Adata = sc.read_h5ad(datadir / 'merged.h5ad')
Adata
sc.pl.umap(Adata, color='stage', palette = 'Spectral')

ph = PipeHandler(Adata, name = 'merged', resdir = resdir)

# In[]
'''     load stage adatas into a dict
'''

Stages = fxa.Stages
Snames = fxa.sampleNames(Stages)

adatas = {stg: sc.read_h5ad(datadir / f'{sn}.h5ad') for \
          stg, sn in zip(Stages, Snames)}

#stage_group_list = g2t.make_stage_group_dict()

# In[]
''' loading tree-structure
'''
_datadir = fxa.DATADIR / 'res-adjust' / '0702-later'
treedir = _datadir / 'treedata'
df_tree = pd.read_csv(treedir / 'tree_struct.csv')
df_tree.head()

# In[]
''' making  lineage dict
    key: lineage_root_node
    value: a list of descendant nodes
'''
## making a dict for looking up descendants 
son_dict = g2t.make_son_dict(df_tree)
fx.save_json_dict(son_dict, treedir / 'son_dict.json')

stage_begin = 'E7'
#lineage_roots = adatas[stage_begin].obs['stg_groups_new'].cat.categories
lineage_roots = [stage_begin + f'_{i}' for i in range(7)]

lineag_dict = {nd: g2t.find_sons(nd, son_dict, with_self=True) \
               for nd in lineage_roots}

rev_lineag_dict = g2t.reverse_dict(lineag_dict)

####################
rev_lineag_dict.update({rt: rt for rt in lineage_roots})
rev_lineag_df = pd.Series(rev_lineag_dict).to_frame('root_group')
rev_lineag_df['lineage'] = rev_lineag_df['root_group'].apply(lambda x: fxa.LineageAnnoAll[x])
rev_lineag_df.to_csv(treedir / 'node_idents_lineage.csv', index=True, index_label='node')

# In[]
'''    
        adding lineage labels 
'''

foo = lambda x: fxa.LineageAnnoAll[rev_lineag_dict.get(x, x)]

foo_s = lambda x: fxa.StageNameDict[x]
# (for stagewise adatas)

key_nd = 'stg_groups_new'
for stg, adt in adatas.items():
    if stg == 'E6': lin_cats = None
    else: lin_cats = fxa.LineageOrd
    
    adt.obs['lineage'] = pd.Categorical(
            adt.obs[key_nd].apply(foo),
            categories=lin_cats)
    print(adt.obs['lineage'].value_counts())
    
    adt.obs['stage_name'] = adt.obs['stage'].apply(foo_s)
    print(adt.obs['stage_name'].value_counts())
    
    
# (for merged Adata)
ph.adata.obs['lineage'] = pd.Categorical(
        ph.adata.obs[key_nd].apply(foo),
        categories=fxa.LineageOrdAll)
ph.adata.obs['stage_name'] = pd.Categorical(
        ph.adata.obs['stage'].apply(foo_s),
        categories = fxa.StageNames)

# In[]
ph.save() # 


# In[]
'''### raw counts'''
fn_merged_raw = fxa.FORMALDIR / 'merged_afterQC.h5ad'
Adata_raw = sc.read_h5ad(fn_merged_raw)

columns = ['stage', 'primer', 'n_genes', 'n_counts', 
           'stage_primer', 'leiden_new', 'refined_group', 
           'parent_bcd', 'lineage', 'stage_name']
for column in columns:
    Adata_raw.obs[column] = Adata.obs[column]

print(Adata_raw)
###Adata_raw.write(fn_merged_raw)
fx.saveNamedMtx(Adata_raw, fxa.FORMALDIR / 'Merged_mtx')


# In[]
''' Test Block (for stage)
''' 
cmap_stg = ['plasma_r', 'Spectral'][1]
_test_colors = fx.get_colors(cmap=cmap_stg, n=9)
fx.view_color_map(_test_colors)

tt = f'Total of {ph.adata.shape[0]} cells'
ph.vis(color='stage_name', palette = cmap_stg, title = tt,
       save = f'_stage_{cmap_stg}.pdf')
#ph.vis(color='stage_name', palette = cmap_stg, title = tt,
#       save = f'_stage_{cmap_stg}.pdf')

 # In[]
''' Test Block (for Lineages)
''' 

cats_lin = ph.adata.obs['lineage'].cat.categories

_tmp = fxa.LineageColorCandi.loc[cats_lin, 'candi1'].tolist()
fx.view_color_map(_tmp)



ph.adata.uns['lineage_colors'] = _tmp
ph.vis(color='lineage', #palette = _test_colors0,
#       title = tt,
       groups = fxa.LineageOrd,
       save='_lineage.pdf')

#ph.vis(color='lineage', groups = ['B_' + str(i) for i in range(3)])
#ph.vis(color='n_counts_disc')


# In[]
''' Test Block (for lineage / sub-plots)
''' 

#cmap_sub = 'Set2_r'
#_test_colors1 = fx.get_colors(cmap=cmap_sub, n=8)
#fx.view_color_map(_test_colors1)
#fx.view_color_map(cmap_sub, n=20)

for stg, adt in adatas.items():
    stg_name = fxa.StageNameDict[stg]
    _cats = adt.obs['lineage'].cat.categories
    adt.uns['lineage_colors'] = fxa.LineageColorCandi.loc[_cats, 'candi0'].tolist()
    tt1 = f'{stg_name} ({adt.shape[0]} cells)'
    sc.pl.umap(adt, color = 'lineage', 
               title = tt1,
               legend_loc = 'right margin',
               save = f'_{stg}_lin.pdf')

    
# In[]
''' Test Block (for expression levels on UMAP)
''' 

cmap_val = ['YlGnBu',
            fx.diy_cmap_grey_bg('Blues'),
             fx.diy_cmap_grey_bg()][-1]

_test_colors = fx.get_colors(cmap=cmap_val, n=100)
fx.view_color_map(_test_colors)

cmap_val_name = cmap_val if isinstance(cmap_val, str) else 'diy'

for gid, gn in fxa.MarkerDictShow.items():
    gexpr = ph.adata.raw[:, gid].X.flatten()
    _vmax = np.quantile(gexpr[gexpr > 0], 0.99)
    
#    ph.vis(gid, title=gn, cmap=cmap_val, 
#           vmax=_vmax,
#           save=f'_show_{cmap_val_name}_{gn}_{gid}.png' #+ '-sub.png'
#           )
    ph.vis(gid, title=gn, cmap=cmap_val, 
           vmax=_vmax,
           save=f'_show_{cmap_val_name}_{gn}_{gid}.pdf' #+ '-sub.png'
           )

gid, gn = "bf_00012984", "CDX"

# In[]
'''     DE between lineages
'''
method_de = ['wilcoxon', 't-test_overestim_var'][1]
ph.DE(groupby = 'lineage', method = method_de)
ph.dot_de(5, groupby='lineage')

''' taking out cell from 7 to 12 (un-run)
'''
#groups_de = ['E_' + str(i) for i in range(7, 13)]
#adt = B.TakeGroups(ph.adata, groups_de, key = 'stage')


# In[]
''' export lineages
'''

_lin_name = fxa.LineageOrd[3]
for _lin_name in fxa.LineageOrd:

    ph.vis(color='lineage', #palette = _test_colors0,
           title = _lin_name,
           groups = [_lin_name],
           save=f'_{_lin_name}.pdf')
    _adt_lin = B.TakeGroups(ph.adata, [_lin_name], key='lineage')
    _adt_lin.write(resdir / f'lin_from_mg_{_lin_name}.h5ad')


# In[]
''' export lineages --- Raw counts
'''

#_lin_name = fxa.LineageOrd[3]
resdir_lin0 = fxa.LINEDIR / '0710'
for _lin_name in fxa.LineageOrd:

    _adt_lin_raw = B.TakeGroups(Adata_raw, [_lin_name], key='lineage')
    
    _adt_lin_raw.write(resdir_lin0 / f'lin_rawcounts_{_lin_name}.h5ad')
    fx.saveNamedMtx(_adt_lin_raw, resdir_lin0 / f'{_lin_name}_mtx', 
                    field='integer')

    _adt_lin_raw.obs[columns].to_csv(resdir_lin0 / f'{_lin_name}_metadata.csv',
         index=True, header=True)

    
# In[]
'''     testing methods for lineage analysis
==========================================================
    
1. directly from the original PC space computed by merging stages
    1.1. neighbor -> UMAP
    1.2. stage NN
2. re-do the pipe-line (normalization, HVG-selection, Z-score, PCA, NN, UMAP)
    2.1. general NN
    2.1 stage NN

-->  embedding performace varies between datasets (lineage):
    for Notochord, better use `2.1` + PAGE + uniform-DPT
    
'''













''' lineage self-embedding
'''

lines = np.take(fxa.LineageOrd, [3])
for lines in fxa.LineageOrd:
    lin_name = lines if isinstance(lines, str) else ' and '.join(lines)
    
    resdir_lin = fxa.LINEDIR / '0710' / f'scanpy/{lin_name}'
    adt_lin = B.TakeGroups(ph.adata, lines, key='lineage')
    adt_lin_raw = B.TakeGroups(Adata_raw, lines, key='lineage')
    
    
    # In[]
    '''     base-line (UMAP on the original PC space)
    '''
    
    from_counts = True
    batch_key = 'primer'
    
    if from_counts:
#        _adt_lin_raw = B.RemoveGroups(adt_lin_raw, ['E14'], key='stage')
        phl = PipeHandler(adt_lin_raw, lin_name, resdir=resdir_lin)
        phl.NormLog_rev(counts_level=300)
        phl.HVGs(min_mean=0.04, min_disp=0.25, # 0.25 default
                batch_key=batch_key,n_bins=30, redo=False,
                plot=True, keep_results=True)
    #    B.ResetX(phl.adata, copy=False)
        phl.scale(groupby=batch_key)
        phl.adata.obs['log_n_counts'] = phl.adata.obs['n_counts'].apply(np.log1p)
        
        phl.PCA(n_comps=50, )
    else:
        phl = PipeHandler(adt_lin, name = lin_name, resdir = resdir_lin)
        phl.vis(color = 'stage', save='000') # original position
    
    metric = ['euclidean', 'cosine'][1]
    nneigh = 15
    n_pcs = 20
    
    if False: 
        B.ResetX(phl.adata)
        phl.scale(groupby = 'primer')
        phl.PCA()
    phl.Neighbors(nneigh = nneigh, do_pca = True, use_rep = 'X_pca',
               n_pcs = n_pcs, metric=metric)#'cosine') # 'euclidean' by default
    phl.Umap(min_dist = 0.2)
    
    ### visualization
    phl.vis(['stage', 'primer'], palette=fxa.CMAP_STAGE0, 
            save=f'_tmp_redo{from_counts}1.png')
    phl.vis(['stage', 'primer'], key='pca', 
            save=f'_tmp_redo{from_counts}1.png',#components='2,3', 
           palette=fxa.CMAP_STAGE)



# In[] 
'''highlight sub-lineage'''
sub_root = 'E9_5'
sub_lin_groups = g2t.find_sons(sub_root, son_dict, with_self=True)

B.label_mask_others(phl.adata.obs, key_nd, keeps=sub_lin_groups)
phl.vis(color = 'tmp', groups = sub_lin_groups)


_gid, _gn = "bf_00012984", "CDX"
phl.vis(color = _gid, cmap = fxa.CMAP_VALUE, vmax = 1.5,
        title = _gn, save=f'_{_gn}_{_gid}')

# In[]
from StagewiseNN import stagewise_knn

X = phl.adata.obsm['X_pca']
#X = phl.adata.X

stage_lbs = phl.adata.obs['stage'][:]
stage_ord = fxa.Stages[1: -1]
    
# ===========================================
leaf_size = 1
#n_pcs = [30] * 5 + [50] * 4
n_pcs = 10
k = [20] * 5 + [10] * 3
#k = [20] + [15] * 4 + [10] * 3
sigma = [_k / 2 for _k in k]

do_pca = False if X.shape[1] < 100 else True
distmat, connect0 = stagewise_knn(X, stage_lbs, stage_ord, 
                                 k=k, leaf_size=leaf_size, do_pca=do_pca,
                                 n_pcs = n_pcs,
                                 norm_dists=True,
                                 sigma=sigma, norm_sample=True)
#sparse.save_npz(resdir_lin / 'formal_distmat.npz', distmat)
#sparse.save_npz(resdir_lin / 'formal_connect0.npz', connect0)

connect = connect0.copy()
connect[connect > 0] = 1
# ===========================================
phl.adata.uns['neighbors'] = {
        'params': {'n_neighbors': k, 'method': 'umap', 'metric': 'cosine'},
        'connectivities': connect, #.to_csr(),
        'distances': distmat,}

sc.tl.umap(phl.adata, min_dist=0.1, spread=1)
phl.vis(color=['stage_name', 'primer'], 
        save='_tmp_stgnn1',
        palette = 'Spectral', )#palette=fxa.CMAP_STAGE,)


# In[]




# In[]
denoise_graph = False
fxa.PH_pseudotime(phl, root_group='E7', groupby='stage', 
                  diff_comps=15, n_dcs=5, )
if denoise_graph:
    phl.Neighbors(nneigh = nneigh * 2, do_pca = True, use_rep = 'X_diffmap',
               n_pcs = 4, metric='euclidean',#metric,
               random_state=2020)#'cosine') # 'euclidean' by default
    phl.Umap(min_dist = 0.3)


phl.vis('refined_group', )





















