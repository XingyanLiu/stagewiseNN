# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 22:20:43 2019

@author: xyliu
"""
import os
from pathlib import Path
from typing import Sequence, Union, Mapping, List

import json
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import scanpy as sc

from scipy import sparse, io

def Snames(old = True, ):
    Snames_dict = {'adults': ['T%d'%i for i in range(1, 15)], 
              'embryos_02': ['E%d_02'%i for i in range(1, 16)], 
              'embryos': ['E%d'%i for i in range(1, 16)], 
              'AE': ['A', 'E'], 
              'NS': ['D%d'%i for i in range(1, 7)], 
              'BF': 'B2 B4 B5 B7 B8'.split(),
              'BF_10G': ['B_10G', 'F_10G'],
              'All': list('AEBF'),
              'KO': ['KO_1'],
              'test': ['T11', 'T11mtx']}
    Snames_dict['adults'].remove('T6')

    if old:
        Snames_dict['adults'][10] = 'T12old'
        Snames_dict['embryos_02'][9] = 'E10old'

    return Snames_dict


###-----------------------------------------------------------------###  
# In[]
''' Manuplating datasets (0)
'''
def change_names(gnames, foo_change, **kw):
    '''
    foo_change: function to map a name-string to a new one
    **kw: other kwargs for foo_change
    
    '''
    return list(map(foo_change, gnames, **kw))
  
def shorter_name(name:str, sep=(',', ' ')):
    '''
    test code
    ----------
    name = 'HNRNPUL2,HNRNPU,HNRNPUL1'
    name = 'Cyclic nucleotide-gated channel cone photoreceptor subunit alpha'
    shorter_name(name)


    '''
    if isinstance(sep, str):
        parts = name.split(sep)
        if len(parts) >= 3:
            name = sep.join(parts[:2] + ['..']) 
    else:
        for s in sep:
            name = shorter_name(name, s)
    
    return name

def reverse_dict_of_list(d: Mapping, ) -> Mapping:
    '''
    the values of the dict must be list-like type
    '''
    d_rev = {}
    for k in d.keys():
        vals = d[k]
        _d = dict.fromkeys(vals, k)
        d_rev.update(_d)
    return d_rev


# In[]
def gene_id_name_mapping(species='amph', todict=True):
    '''
    species: either be "mouse" or "amph"
    return a dictionary with 'gene_id' as keys, and 'gene_name' as values
    '''
    anno_file = Path(f'resources/gene_id_name_{species}.csv')
    mapping = pd.read_csv(anno_file, index_col=0)['gene_name']#.to_dict()
    print('ID-name mapping of length:', len(mapping))
    return mapping.to_dict() if todict else mapping

def rename_gene_mouse(gene_id, mapping = None):
    '''TEST
    rename_gene_mouse('ENSMUSG00000051951')
    rename_gene_mouse('ENSMUSG00000102851')
    '''
    mapping = gene_id_name_mapping('mouse') if mapping is None else mapping
    if gene_id in mapping.keys():
        gname_new = mapping[gene_id]
    else:
#        print('No coresponding name for gene id [%s]!'%gene_id)
        print(gene_id)
        gname_new = gene_id
    return gname_new

def rename_gene_amph(gene_id, mapping = None):
    '''TEST
    rename_gene_amph('BFwtdbg0000020')
    rename_gene_amph('BFwtdbg0000030')
    '''
    mapping = gene_id_name_mapping('amph') if mapping is None else mapping
    gname_new = gene_id
    gname = mapping[gene_id]
    if gname is not np.nan:
        gname_new = gene_id.strip('BFwtdbg') + '_%s' % gname
    return gname_new

##    head = {'amph': 'BFwtdbg', 'mouse': 'ENSMUSG'}


def RenameGenes(gene_ids, species='amph', mapping = None, tolist=False):
    '''
    `gene_ids` is a list of gene-IDs
    '''
    print('the old gene names:', gene_ids[:5])
    mapping = gene_id_name_mapping(species) if mapping is None else mapping
    rename_gene = {
            'amph': rename_gene_amph,
            'mouse': rename_gene_mouse
            }[species]
#    gnames = [rename_gene(gid, mapping) for gid in gene_ids]
    gnames = pd.Series(gene_ids).apply(lambda x: rename_gene(x, mapping))
    print('the new gene names:', gnames[:5])
    return gnames.to_list() if tolist else gnames
###-----------------------------------------------------------------###
# In[]

def _check_path_class(p):
    print('path name: `{}`'.format(p))
    if not isinstance(p, Path):
        return Path(p)
    return p

def check_dirs(*paths):
    for path in paths:
        if os.path.exists(path):
            print('already exsists:\n\t%s'%path)
        else:
            os.makedirs(path)
            print('a new directory made:\n\t%s'%path)

#def csv2list(fname, header=None):
#    col = pd.read_csv(fname, header=header).iloc[:, 0]
#    return col.to_list()
        
def saveNamedMtx(adata, dirname, field=None, raw=True, **kwds):
    check_dirs(dirname)
    ''' better set field='integer' if possible '''
    if adata.raw is not None:
        adata = adata.raw
    
    mtx_file = '%s/matrix.mtx'% (dirname)
    bcd_file = '%s/barcodes.tsv'% (dirname)
    gene_file = '%s/genes.tsv'% (dirname)
    
    genes = adata.var_names.to_frame(index = False, name='genes')
    barcodes = adata.obs_names.to_frame(index = False, name='barcodes')
    print(adata.X.data[:5])
    genes.to_csv(gene_file, index = False, header = False)
    barcodes.to_csv(bcd_file, index = False, header = False)
    sparse.save_npz(f'{dirname}/matrix.npz', adata.X.T)
    io.mmwrite(mtx_file, adata.X.T, field = field, **kwds)
    
    print('Matrix saved to directory `%s`'% dirname)


def saveMtx2df(adata, fname, index=True, header=True, **kwds):
    print("NOTE: saving the dense matrix might take some time. \
          if needed, you can consider the function: `funx.saveNamedMtx`\
          to handle a sparse matrix in a efficient way.")
    df = mtx2df(adata)
    df.to_csv(fname, index=index, header=header)



def mtx2df(adata):
    df = pd.DataFrame(data = adata.X.toarray(),
                      index=adata.obs_names,
                      columns=adata.var_names)
    return df



def load_namelist(fpath, header=None, tolist=True, **kw):
    '''
    only return the first column by default
    '''
    names = pd.read_csv(fpath, header=header, **kw)
    names = names.iloc[:, 0]
    return names.tolist() if tolist else names 

def save_namelist(lst, fname, header=False, index=False, **kw):
    pd.Series(lst).to_csv(fname, header=header, index=index, **kw)
    print('name list seved into:', fname)


def load_dense(fpath, ):
    print(f'loading dense matrix from {fpath}')
    mat = pd.read_csv(fpath, sep='\t', index_col=0)
    return mat


def load_sparse(fpath, backup_npz=True):
    '''
    `fpath` should be ended with '.mtx' or 'npz'
    '''
    fpath = _check_path_class(fpath)
    if fpath.suffix == '.mtx':
        mat = io.mmread(str(fpath))
        if backup_npz:
            print('backup the `.npz` file for speedup loading next time...')
            mat = sparse.csc_matrix(mat)
            sparse.save_npz(fpath.parent / 'matrix.npz', mat)
    elif fpath.suffix == '.npz':
        mat = sparse.load_npz(str(fpath))
    else: 
        raise ValueError("the file path should be ended with '.mtx' or 'npz'")
    return mat

def add_obs_annos(adata: sc.AnnData, 
                  df: Union[Mapping, pd.DataFrame],
                  ignore_index = True,
                  copy=False):
    adata = adata.copy() if copy else adata
    print('adding columns to `adata.obs`:')
    for col in df.keys():
        print(col, end=', ')
        adata.obs[col] = list(df[col]) if ignore_index else df[col]
    print('done!')
    return adata if copy else None

def add_var_annos(adata: sc.AnnData, 
                  df: Union[Mapping, pd.DataFrame],
                  copy=False):
    adata = adata.copy() if copy else adata
    print('adding columns to `adata.var`:')
    for col in df.keys():
        print(col, end=', ')
        adata.var[col] = list(df[col])
    print('done!')
    return adata if copy else None

def make_adata(mat, 
               obs: Union[Mapping, pd.DataFrame, None]=None, 
               var: Union[Mapping, pd.DataFrame, None]=None, 
               fn=None):
    ''' An alternative way to construct AnnData (BUG fixed).
    Something might go wrong when saving the AnnData object constructed 
    by the build-in function `sc.AnnData`.
    '''
    if not sparse.issparse(mat):
        mat = sparse.csr_matrix(mat)
    adata = sc.AnnData(mat)
    if isinstance(obs, pd.DataFrame):
        adata.obs_names = obs.index.values
    if isinstance(var, pd.DataFrame):
        adata.var_names = var.index.values
    if obs is not None:
        add_obs_annos(adata, obs, copy=False)
    if var is not None:
        add_var_annos(adata, var, copy=False)
    
    if fn is not None:
        adata.write(fn)
    return adata

def adataFromRaw(dirname, backup_npz=True, name_mtx='matrix', **kw):
    '''
    read matrix and annotated names from directory `dirname`
    --- Example ---
    |-dirname
        |-matrix.npz
        |-genes.tsv
        |-barcodes.tsv (or meta_cell.tsv)
    ==========================================
    Alternative (for load mtx file):
        `adata = sc.read_10x_mtx(dirname, )`
        
    '''
    fn_mat_npz = f'{name_mtx}.npz'
    files = os.listdir(dirname)
    if fn_mat_npz in files:
        mat = load_sparse(dirname / fn_mat_npz)
    elif f'{name_mtx}.mtx' in files:
        mat = load_sparse(dirname / f'{name_mtx}.mtx', backup_npz=backup_npz)
    else:
        raise FileNotFoundError(f"None of file named `{name_mtx}.npz` or `{name_mtx}.mtx` exists!")
#    if sparse.isspmatrix_coo(mat):
#        mat = sparse.csc_matrix(mat) 
    
    params = dict(header=None, sep='\t', index_col=0)  
    if 'meta_cell.csv' in files:
        barcodes = pd.read_csv(dirname / 'meta_cell.tsv', sep='\t', index_col=0)
    else:
        barcodes = pd.read_csv(dirname / 'barcodes.tsv', **params)
    
    genes = pd.read_csv(dirname / 'genes.tsv', **params)
    adata = sc.AnnData(X=mat.T, obs=barcodes, var=genes)
    adata.X = sparse.csr_matrix(adata.X)
    print(adata)    
    return adata



def AdataSubsetGenes(adt, gene_set, copy = True):
    # gene_set should be a list-like object
    if type(gene_set) is pd.core.series.Series:
        gene_set = gene_set.to_list()
    adt = adt.copy() if copy else adt
    use_genes = [g in gene_set for g in adt.var_names]
    
    adt = adt[:, use_genes]# adt[:, list(gene_set)]
    print('there are %d genes present in the given `gene_set`(%d in total)'%
          (np.sum(use_genes), len(gene_set)))
    return adt


def GetGeneSet(dirname, snames, mode = 'de', ntop=50, ntotal=2000):
    print('dirname: {}\nMode: \t{}'.format(dirname, mode))
    
    used_genes = []
    
    if mode == 'de':
        for i in range(len(snames)):
            sn = snames[i]
            degs = pd.read_csv('%s/%s_marker_names.csv' % (dirname, sn))
            top_degs = pd.unique(degs.head(ntop).values.flatten()).tolist()
            print('%s:\tusing top %d DE genes, where %d unique genes'%(
                    sn, ntop, len(top_degs)))
            used_genes += top_degs
        used_genes = pd.unique(used_genes)
        print('there are %d unique DE genes'%(len(pd.unique(used_genes))))
#        if len(used_genes) > ntotal: 
#            used_genes = used_genes[: ntotal]
            
    elif mode == 'hv':
        for i in range(len(snames)):
            sn = snames[i]
            hvgs = pd.read_csv('%s/%s_hvgs.csv' % (dirname, sn), header=None)
            print('%s:\tusing %d HV genes'%(sn, len(hvgs)))
            hvgs = hvgs.iloc[:, 0].tolist()
            used_genes += hvgs
#        used_genes = pd.unique(used_genes)
        freq = pd.Series(used_genes).value_counts(sort=True, ascending=False)
        if len(freq) > ntotal: 
            print('the top %d frequency gene appears in %d datasets'%(
                    ntotal, freq[ntotal]))
            used_genes = freq.head(ntotal).index.tolist()
        else:
            used_genes = freq.index.tolist()
    print('using [ %d ] genes finaly' % len(used_genes))
    return used_genes

# In[]
'''     merge functions
'''
def merge_adatas(adatas, union=True, obs_keys=None, ):
    ''' (a new adata will be created, use the raw matrix)
    adatas: a list of AnnData objects, assume that adata.raw if not None
    
    '''

    print("merging raw matrix...")
#    merged_mat, genes, bcds = merge_adata_mats(adatas, raw=True, union=union)
    mats = []
    genes = []
    bcds = []
#    print("merging matrix...")
    for adt in adatas:
        data = adt.raw if adt.raw is not None else adt
        mats.append(data.X)
        genes.append(data.var_names)
        bcds.append(data.obs_names)
            
    mats, genes = merge_datasets(mats, genes, union=union)
    merged_mat = sparse.vstack(mats, dtype = mats[0].dtype)
    if obs_keys is None:
        obs = pd.Index(np.concatenate(bcds, axis=0))
    else:
        print("merging metadata...")
        obs = merge_metas(adatas, obs_keys)
    genes = pd.DataFrame(index=genes)
    print(genes.head())
    print(obs.head())
    return sc.AnnData(merged_mat, obs=obs, var=genes)
    
def merge_metas(adatas, obs_keys):
    
    obs_keys = [obs_keys] if isinstance(obs_keys, str) else obs_keys
    obss = [adt.obs[obs_keys] for adt in adatas]
    return pd.concat(obss, axis=0)
    

def merge_hvgs(adatas, freq=True):
    
    HVGs = []
    for adt in adatas:
        hvgs = adt.var_names.to_list()
        HVGs += hvgs
    HVGs = pd.value_counts(HVGs)
    
    return HVGs if freq else HVGs.index.to_list()

def merge_adata_mats(adatas, raw=True, union=True):
        
    mats = []
    genes = []
    bcds = []
#    print("merging matrix...")
    for adt in adatas:
        data = adt.raw if raw and adt.raw is not None else adt
        mats.append(data.X)
        genes.append(data.var_names)
        bcds.append(data.obs_names)
    mats, genes = merge_datasets(mats, genes, union=union)
    merged_mat = sparse.vstack(mats, dtype = mats[0].dtype)
    bcds = np.concatenate(bcds, axis=0)
    return merged_mat, genes, bcds

    

def merge_datasets(datasets, genes, ds_names=None, verbose=True,
                   union=False):
    '''
    This function code is copied from `scanorama.merge_datasets`.

    parameters
    ----------
    datasets:
        a list of cell-by-gene matrices (np.ndarray or sparse.csr_matrix) 
    genes:
        a list of gene-list corresponding to the columns of each matrix in 
        `datasets`
    
    returns
    -------
    datasets:
        a list of cell-by-gene matrices (sparse.csr_matrix) with aligned genes 
        columns.
    ret_genes:
        a gene list shared by all the datasets
    
    '''
    import sys
    if union:
        sys.stderr.write(
            'WARNING: Integrating based on the union of genes is '
            'highly discouraged, consider taking the intersection '
            'or requantifying gene expression.\n'
        )

    # Find genes in common.
    keep_genes = set()
    for idx, gene_list in enumerate(genes):
        if len(keep_genes) == 0:
            keep_genes = set(gene_list)
        elif union:
            keep_genes |= set(gene_list)
        else:
            keep_genes &= set(gene_list)
        if not union and not ds_names is None and verbose:
            print('After {}: {} genes'.format(ds_names[idx], len(keep_genes)))
        if len(keep_genes) == 0:
            print('Error: No genes found in all datasets, exiting...')
            exit(1)
    if verbose:
        print('Found {} genes among all datasets'
              .format(len(keep_genes)))

    if union:
        union_genes = sorted(keep_genes)
        for i in range(len(datasets)):
            if verbose:
                print('Processing data set {}'.format(i))
            X_new = np.zeros((datasets[i].shape[0], len(union_genes)))
            X_old = sparse.csc_matrix(datasets[i])
            gene_to_idx = { gene: idx for idx, gene in enumerate(genes[i]) }
            for j, gene in enumerate(union_genes):
                if gene in gene_to_idx:
                    X_new[:, j] = X_old[:, gene_to_idx[gene]].toarray().flatten()
            datasets[i] = sparse.csr_matrix(X_new)
        ret_genes = np.array(union_genes)
    else:
        # Only keep genes in common.
        ret_genes = np.array(sorted(keep_genes))
        for i in range(len(datasets)):
            # Remove duplicate genes.
            uniq_genes, uniq_idx = np.unique(genes[i], return_index=True)
            datasets[i] = datasets[i][:, uniq_idx]

            # Do gene filtering.
            gene_sort_idx = np.argsort(uniq_genes)
            gene_idx = [ idx for idx in gene_sort_idx
                         if uniq_genes[idx] in keep_genes ]
            datasets[i] = datasets[i][:, gene_idx]
            assert(np.array_equal(uniq_genes[gene_idx], ret_genes))

    return datasets, ret_genes


def merge_metafiles(datadir, fnames, cols=None, 
                    index_col=0, **kwds):
    datadir = Path(datadir)
    print('dir-name:', datadir)
    print('file names:')
    print('\n'.join(fnames))
    if cols is None:
        dfs = [pd.read_csv(datadir / fn, index_col=index_col, **kwds)\
                   for fn in fnames]
    else:
        dfs = [pd.read_csv(datadir / fn, index_col=index_col, **kwds)[cols]\
                   for fn in fnames]
    
    return pd.concat(dfs, axis=0)
#from numba import njit
#from scipy.stats import spearmanr

#@njit()
#def dist_Spearman(u, v, axis=1):
#    '''
#    axis : int or None, optional
#        If axis=0 (default), then each column represents a variable, with
#        observations in the rows. 
#        If axis=1, the relationship is transposed: each row represents a 
#            variable, while the columns contain observations.
#        If axis=None, then both arrays will be raveled.
#    '''
#    rho, pval = spearmanr(u, v, axis=axis)
#    return 1 - rho
# In[]  

def save_json_dict(dct, fname='test_json.json', encoding='utf-8'):
    
    with open(fname, 'w', encoding=encoding) as jsfile:
        json.dump(dct, jsfile, ensure_ascii=False)
    print(fname)

def get_json_dict(fname, encoding = 'utf-8'):
    with open(fname, encoding = encoding) as f:
        dct = json.load(f)
    return dct





# In[]
'''
#########################[ Plotting fumctions ]############################
'''

def view_color_map(cmap='viridis', n=None, figsize=(6, 2), s=150, k=20,
                   grid=False, **kwds):
    '''
    n: total number of colors
    k: number of colors to be plotted on each line.
    
    Examples
    --------
    import funx as fx
    colors = ['Set1', 'viridis', 'Spectral']
    fx.view_color_map(colors[-1], n=20)
    
    cmap = sc.pl.palettes.zeileis_26
    cmap = sc.pl.palettes.default_64
    fx.view_color_map(cmap, k=16)    
    '''
    if not isinstance(cmap, (np.ndarray, list)):
        from matplotlib import cm
        cmp = cm.get_cmap(cmap, n)
        n = 20 if n is None else n
        cmp = [cmp(i) for i in range(n)]
    else:
        n = len(cmap) if n is None else n
        cmp = cmap
        
    plt.figure(figsize=figsize)
    for i in range(n):
        plt.scatter(i % k, i // k, color=cmp[i], s=s)
    plt.grid(b=grid,)#axis='y')
    plt.show()


def get_colors(cmap='Spectral', n=5, to_hex=True):
    '''
    fx.get_colors('Reds', 4)
    fx.view_color_map('Reds', 4)
    '''
#    import matplotlib.colors as mcolors
    cmp = plt.cm.get_cmap(cmap, n)
    colors = [cmp(i) for i in range(n)]
    if to_hex:
        return [mcolors.to_hex(c) for c in colors]
    else:
        return colors

def diy_cmap_grey_bg(name_fg = 'RdPu', low = 0.15, rm_high=0.01, n=100):
    
    s = int(n * low)
    t = max((1, int(n * rm_high)))
    print((s, t))
    candi_colors = ['#d9d9d9']*s + get_colors(name_fg, n)[s: -t]
    cmap = mcolors.ListedColormap(candi_colors)
    return cmap



def buttomViewAdt(adata, key='umap', groupby='batch', size=2, 
                  save=None):
    '''
    adata
    key = {'pca', 'umap', 'tsne', ..}
    groupby = {'leiden', 'batch', ...} (any key in adata.obs.keys())
    save: boll, None, str
    '''
    if groupby not in adata.obs.keys():
        raise KeyError('The groupby-key `{}` is not in `adata.obs.keys()`, which includes {}'.format(
                groupby, adata.obs.keys()))
    group_names = pd.unique(adata.obs[groupby])
    for gn in group_names:
        figname = '{}_{}_{}.png'.format(key, groupby, gn)
        sc.pl.scatter(adata, basis=key, color=groupby, 
                      groups=[gn], save=figname, 
                      size=size,)
     ###[ Note ]--> palette='tab20' is ignored, bug on scanpy.pl.scatter!
#        sc.pl.umap(adata, color=groupby, groups=[gn], 
#                   save=figname, size=size,
#                   palette='tab20')   
    
    

def ButtomViewAdt(adata, key='umap', groupby='batch', sz=0.1, 
                  disc=True, overall=True, cmap=None, tail='',
                  save=None, **kwds):
    '''
    adata
    key = {'pca', 'umap', 'tsne'}
    groupby = {'leiden', 'batch', ...} (any key in adata.obs.keys())
    save: boll, None, str(figure dir)
    '''
    if groupby not in adata.obs.keys():
        raise KeyError('The groupby-key `{}` is not in `adata.obs.keys()`, which includes {}'.format(
                groupby, adata.obs.keys()))
    X = adata.obsm['X_%s'%key]
    labels = adata.obs[groupby]
    ButtomView(X, labels, key = key, sz=sz, disc=disc, 
               cmap=cmap, overall=overall,
               save=save, name='%s_%s_%s'%(key, groupby, tail), **kwds)


def ButtomView(X, labels, unique_labels=None, key='X', 
               sz=0.5, disc = True, overall=True,
               cmap=None, save=None, name='tmp',
               figsize=(8, 8), 
               **kwds):
    '''
    X: np.array
    labels: 1-d np.array or pd.Series
    disc: bool, whether to use discrete color map
    save: boll, None, str(figure dir)
    '''
    if unique_labels is None:
        unique_labels = np.unique(labels) # `pd.unique()` does not sort labels
    
    ## set colormap
    from matplotlib import cm
    if cmap is None:
        if disc:
            cmap = 'tab20'
        else:
            cmap = 'Spectral'#'RdYlBu'
    lut = None if disc else len(unique_labels)
    cmp = cm.get_cmap(cmap, lut)
    
    ## plot for each given labels
    for i, lb in enumerate(unique_labels):
        _X = X[labels==lb, :]
        __X = X[labels != lb, :]
        fig, ax = plt.subplots(figsize=figsize)
        ax.scatter(__X[:, 0], __X[:, 1], label = 'others', s=sz,
                    color=(0.9, 0.9, 0.9), **kwds)
        ax.scatter(_X[:, 0], _X[:, 1], label = lb, s = sz * 2,
                    color=cmp(i, alpha=0.6), **kwds)
        ax.legend(loc='upper right')
        
        plt.grid(False)
        plt.xlabel(key + '1')
        plt.ylabel(key + '2')

        if isinstance(save, (Path, str)):
            figname = '%s/%s_%s.png'%(save, name, lb)
            print('saving figure to %s...'%figname)
            fig.savefig(figname)
            plt.close()
     
    if overall:
        mask = [lb not in unique_labels for lb in labels]
        __X = X[mask, :]
        fig, ax = plt.subplots(figsize=figsize)
        if np.sum(mask) >= 1:
            __X = X[mask, :]
            ax.scatter(__X[:, 0], __X[:, 1], label = 'others', s=sz,
                color=(0.9, 0.9, 0.9), **kwds)

        for i, lb in enumerate(unique_labels):
            _X = X[labels==lb, :]
            ax.scatter(_X[:, 0], _X[:, 1], label = lb, s=sz,
                        color=cmp(i, alpha=0.6),)
        ax.legend(loc='best')
        plt.grid(False)
        plt.xlabel(key + '1')
        plt.ylabel(key + '2')
        if isinstance(save, (Path, str)):
            figname = '%s/%s_overview.png'%(save, name)
            print('saving figure to %s...'%figname)
            fig.savefig(figname)
            plt.close()
    

def ButtomCompareAdt(adata, pair_labels, key='umap',  
                     groupby='batch', sz=0.1, 
                      cmap=None, name=None,
                      save=None, **kwds):
    '''
    adata
    key = {'pca', 'umap', 'tsne'}
    groupby = {'leiden', 'batch', ...} (any key in adata.obs.keys())
    save: boll, None, str(figure dir)
    '''
    if groupby not in adata.obs.keys():
        raise KeyError('The groupby-key `{}` is not in `adata.obs.keys()`, which includes {}'.format(
                groupby, adata.obs.keys()))
    X = adata.obsm['X_%s'%key]
    labels = adata.obs[groupby]
    pname = '(%s)'% ('_'.join(pair_labels))
    name = '%s_%s%s'%(key, groupby, pname) if name is None else name
    ButtomCompare(X, labels, pair_labels, key = key, 
               cmap=cmap, sz=sz, 
               save=save, name=name, **kwds)


def ButtomCompare(X, labels, pair_labels, key='X', 
                   sz=0.5, 
                   cmap=None, save=None, name='tmp',
                   figsize=(8, 8), 
                   **kwds):
    '''
    X: np.array
    labels: 1-d np.array or pd.Series
    pair_labels: a list of 2 label-names to be compared
    disc: bool, whether to use discrete color map
    save: bool, None, str(figure dir)
    '''
#    name = '_'.join(map(str, pair_labels))
    ## set colormap
    from matplotlib import cm
    cmap = cmap if cmap is not None else 'Set1' # 'tab10'#
    cmp = cm.get_cmap(cmap, len(pair_labels))
    
    mask = np.array([lb in pair_labels for lb in labels])
#    mask = pd.Series(labels).apply(lambda x: x in pair_labels).values
#    print('type of X:', type(X), 'shape:', X.shape)
#    print('type of `mask`:', type(mask),)
#    print(mask[:5]) #( for debudding !!! )
    fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(X[~ mask, 0], X[~ mask, 1], label = 'others', s=sz,
                color=(0.9, 0.9, 0.9), **kwds)
#    X_emph = X[mask, :]
    ## plot for each given labels
    for i, lb in enumerate(pair_labels):
        _X = X[labels==lb, :]
        ax.scatter(_X[:, 0], _X[:, 1], label = lb, s = sz * 2,
                    color=cmp(i, alpha=0.6), **kwds)
    
    ax.legend(loc='upper right')
    plt.grid(False)
    plt.xlabel(key + '1')
    plt.ylabel(key + '2')

    if isinstance(save, (Path, str)):
        figname = '%s/%s.png'%(save, name)
        print('saving figure to %s...'%figname)
        fig.savefig(figname)
        plt.close()
        

def hist_log10(Adata, ft_name = 'n_genes', bins = 50):
    plt.figure()
    plt.hist(np.log10(Adata.obs[ft_name].values + 1), bins = bins)
    plt.xlabel('log10_%s'%ft_name)
    plt.ylabel('number of cells')
    plt.title('distribution of log10_%s'%ft_name)
    plt.show()
    
def plot_ranked_log10(Adata, ft_name = 'n_genes'):
    plt.figure()
    plt.plot(np.log10(Adata.obs[ft_name].sort_values().values + 1))
    plt.xlabel('ranked cells')
    plt.ylabel('log10_%s'%ft_name)
    plt.title('ranked values of log10_%s'%ft_name)
    plt.show()
    
def printmm(Adata, ft_name = 'n_genes'):
    # print min \ mean \ median \ max of the given observation feature
    print('The min, median and mean of `%s` are:'%ft_name)
    print('\tmin: %.1f\n\tmedian: %.1f\n\tmean: %.1f\n\tmax: %.1f'%(
            Adata.obs[ft_name].min(), Adata.obs[ft_name].median(), 
            Adata.obs[ft_name].mean(), Adata.obs[ft_name].max()))
    


def plot_combined_log10(Adata, ft_name = 'n_genes', bins = 50, save = None):
    
    fig, axs = plt.subplots(1, 2, figsize = (12, 4), sharey=False)
    axs[0].hist(np.log10(Adata.obs[ft_name].values + 1), bins = bins)
    axs[0].set_xlabel('log10_%s'%ft_name)
    axs[0].set_ylabel('number of cells')
    #ph.adata.obs['n_genes'].sort_values(ascending=False).plot(use_index=False, logx=True, logy=True)
    axs[1].plot(np.log10(Adata.obs[ft_name].sort_values().values + 1))
    axs[1].set_xlabel('ranked cells')
    axs[1].set_ylabel('log10_%s'%ft_name)
    
    fig.suptitle('distribution of log10_%s'%ft_name)
    if isinstance(save, (Path, str)): fig.savefig(save)
    plt.show()

    
def plot_pc_vr(Adata, n = 30, save = None):
    tmp = Adata.uns['pca']['variance_ratio'].copy()
    tmp.sort()
    tmp = tmp[::-1]

    fig, axs = plt.subplots(1, 2, figsize = (12, 4), sharey=False)
    sc.pl.pca(Adata, color = 'n_counts', ax = axs[0], show = False) # should be placed afterwise!!!
    axs[1].plot(np.log(tmp[:n]), 'o')
    axs[1].set_ylabel('log- variance ratio')
    axs[1].set_xlabel('ranked PCs')
    if isinstance(save, str): fig.savefig(save)
    plt.show()
    
#ph.adata.obs['n_genes'].sort_values(ascending=False).plot(use_index=False, logx=True, logy=True)

def DfBar(df, y, title='Bar plot', figsize=(8, 4), save=None, 
          grid = True, **kwds):
    
    df.plot.bar(y = y, title = title, figsize = figsize, grid=grid,  **kwds)
    if save:
        figdir =  save if isinstance(save, (str, Path)) else 'tmp_figs'
#        check_dirs(figdir) 
        plt.savefig('%s/%s.png' % (figdir, title)) 
        
def DfBarTable(df, y, title='Bar plot', figsize=(8, 4), save=None,
          grid = True, **kwds):
    
    fig, ax = plt.subplots(1, 1)
    ax.get_xaxis().set_visible(False)
    df.plot.bar(y = y, title = title, table = True, ax=ax,
                figsize = figsize, grid=grid, **kwds)
    if save:
        figdir =  save if isinstance(save, (str, Path)) else 'tmp_figs'
#        check_dirs(figdir) 
        plt.savefig('%s/%s.png' % (figdir, title)) 
#    plt.close()


########################################################################

def VennPlot(sets, set_labels=None, regular=False,
             tt='Venn plot', saveto=None):
    '''
    sets: iterable
    set_labels: list[str]
    regular: bool, only for sets of strings!
        wether to regularize the strings in sets 
    '''
    from matplotlib_venn import venn2, venn3
    if regular:
        print('Regularizing strings in sets')
        sets = list(map(set, map(lowerStr, sets)))
    if not isinstance(sets[0], set):
        # assure that each item in `sets` is a set object
        sets = list(map(set, sets)) 
    n = len(sets)
    if set_labels is None:
        set_labels = np.arange(n)
    venn = {2: venn2, 3: venn3}[n]
    fig, ax = plt.subplots()
    venn(sets, 
          set_labels = set_labels)
    plt.title(tt)
    plt.show()
    if isinstance(saveto, (str, Path)):
        plt.savefig(saveto / f'{tt}.png')
    
    return ax

def lowerStr(strings):
    return [c.lower() for c in list(strings)]

def upperStr(strings):
    return [c.upper() for c in list(strings)]

########################################################################
def _vk(namespace):
    '''
    namespace: dict, e.g. `global()` or `local()`
    '''
    vk = {} 
    place = namespace.copy()
    for k,v in place.items(): 
        vk[ id(v) ] = k 
    return vk

def var_name(p, namespace,): 
    vk = _vk(namespace) 
    return vk[ id(p) ] 

def namestr(obj, namespace):
    return [name for name in namespace if namespace[name] is obj]
########################################################################




def close_figs(n = 200):
    for i in range(n):
        plt.close()

from sklearn import metrics
def plot_confusion_matrix(y_true, y_pred, symetric=True, #classes,
                          norm = True,
                          axis = 1, figsize=(12, 8),
                          xlabel=None, ylabel=None,
                          title=None, text_on = False, show_matrix= False,
                          cmap=plt.cm.Blues,
                          save=None):
    """
    This function prints and plots the confusion matrix.
    Normalization can be applied by setting `normalize=True`.
    axis: 0 or 1, define the axis for reorder and normalization.
        if axis=0, perform column-wise normalization to unique-sum;
        if axis=1, row-wise normalization
    """
    if not title:
        title = 'Confusion matrix'

    # Compute confusion matrix
    if symetric:
        mat = metrics.confusion_matrix(y_true, y_pred)
        # Only use the labels that appear in the data
        classes0 = np.unique(list(y_true) + list(y_pred))
        classes1 = classes0
        classes0 = None
    else:
        mat = metrics.cluster.contingency_matrix(y_true, y_pred)
        mat = order_contingency_mat(mat, axis=axis)
        classes0 = np.unique(y_true)
        classes1 = np.unique(y_pred)
    if show_matrix: print(mat)
    if norm:
#    cm_norm = mat.astype('float') / mat.sum(axis=axis)[:, np.newaxis]
        mat = pd.DataFrame(mat).apply(lambda x: x / np.sum(x), axis = axis)
        by = {0: 'Column', 1: 'Row'}[axis]
        print(f"Colored by {by}-Normalized confusion matrix")
#        print(mat.sum(axis=axis)[:5])
    if ylabel is None:
        ylabel = y_true.name if hasattr(y_true, 'name') else ''
    if xlabel is None:
        xlabel = y_pred.name if hasattr(y_pred, 'name') else ''
            
    
    kw_plt = dict(
        figsize=figsize,
            text_on=text_on,
            xnames=classes0, ynames=classes1,
            xlabel=xlabel, ylabel=ylabel,
            title=title, 
            save=save,
            )

    return plot_matrix(mat, **kw_plt)

def plot_matrix(mat,  
                xnames=None, ynames=None,
                xlabel='', ylabel='',
                figsize=(10, 8), title='', 
                text_on = False, show_matrix= False,
                cmap=plt.cm.Blues, 
                save=None,
                show=True):
    """
    This function plots the input matrix.
    """
#    cm_norm = mat.astype('float') / mat.sum(axis=1)[:, np.newaxis]
#    print("Colored by Row-Normalized confusion matrix")        
    if show_matrix: print(mat)
    axsettings = dict(
            xticks=np.arange(mat.shape[1]),
            yticks=np.arange(mat.shape[0]),
            ylabel=ylabel, xlabel=xlabel,
            title=title,
            )
    fig, ax = plt.subplots(figsize=figsize)
    if max(mat.shape) >= min(mat.shape) * 2:
        symetric=False
        im = ax.pcolor(mat, cmap=cmap,)
        axsettings['xticks'] = axsettings['xticks'] + 0.5
        axsettings['yticks'] = axsettings['yticks'] + 0.5
    else:
        symetric=True
        im = ax.imshow(mat, interpolation='nearest', cmap=cmap)
    # We want to show all ticks...
    
    ax.figure.colorbar(im, ax=ax, )
    plt.grid(False)
    # ... and label them with the respective list entries
    if xnames is not None:
        axsettings['xticklabels'] = xnames
    elif hasattr(mat, 'columns'):
        axsettings['xticklabels'] = mat.columns
        
    if ynames is not None:
        axsettings['yticklabels'] = ynames
    elif hasattr(mat, 'index'):
        axsettings['yticklabels'] = mat.index
    ax.set(**axsettings)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=90, ha="right",
             rotation_mode="anchor")

    if text_on:
        mat_ = mat.values if hasattr(mat, 'values') else mat
    # Loop over data dimensions and create text annotations.
        fmt = text_on if isinstance(text_on, str) else '{:03.3f}'
        thresh = mat_.max() / 2.
        for i in range(mat_.shape[0]):
            for j in range(mat_.shape[1]):
                i_, j_ = i, j
                if not symetric: i_, j_ = i + 0.5, j + 0.5

                ax.text(j_, i_, fmt.format(mat_[i, j]), #format(mat_[i, j], fmt),
                        ha="center", va="center",
                        color="white" if mat_[i, j] > thresh else "black")

    fig.tight_layout()
#    plt.subplots_adjust(right=0.1)
    
    if isinstance(save, (Path, str)):
        plt.savefig(save)
    
    return ax

def order_contingency_mat(mat, axis = 1):
    '''
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
    
#    order0 = np.argsort(np.argmax(mat, axis=1))
#    cm_ord = mat[order0, :]
#    
#    order1 = np.argsort(np.argmax(mat, axis=0))
#    cm_ord = mat[:, order1]
#    plot_matrix(cm_ord)
        
def order_contingency_df(df: pd.DataFrame, axis = 1):
    '''
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

##############################################################
from sklearn.model_selection import train_test_split
import joblib
import time
class ModelTrainTester(object):
    '''
    Just for convenience, keep the training and testing data in a class
    '''
    def __init__(self, X, y, test_size= 0.1, resdir = 'temp_res',
                 random_state=42):
        
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=test_size, random_state=random_state)
        
        self.X_train = X_train
        self.X_test = X_test
        self.y_train = y_train
        self.y_test = y_test
        self.classes = np.unique(y)
        self.scores = dict()
        self.resdir = Path(resdir)
        check_dirs(resdir)
        print('Shape of the data `X`: ', y.shape)
        print('Shape of the target `y`: ', y.shape)
        print('Classes: ', self.classes)
        print('Testing size: %.2f'% test_size)
        
    def score_init(self, name):
        if name not in self.scores.keys():
            self.scores[name] = {}
    
    
    def train(self, model, name = 'temp_model', save=True,):
        print('=' * 80)
        print('\tTaining Model [ %s ]:'% (name))
        print('-' * 80)
        print(model)
        
        t0 = time.time()
        model.fit(self.X_train, self.y_train)
        tm = time.time() - t0
        print('%s fitted, time used:\t%.4f s'% (name, tm))
        score_train = metrics.accuracy_score(self.y_train, model.predict(self.X_train))
        print("Training accuracy:\t%0.4f" % score_train)
        
        self.score_init(name)
        self.scores[name]['acc_train'] = score_train
        self.scores[name]['time_train'] = tm
        print(model)
        if save:
            resfile = self.resdir / (f'{name}.pkl')
            joblib.dump(model, resfile)
            print(f'Model saved to {resfile}')
        return model
    
    
    def test(self, model, name = 'temp_model', 
             plot_mat = True, return_pred = False):
        t0 = time.time()
        y_pred = model.predict(self.X_test)
        tm = time.time() - t0
        
        print('Prediction on testing data completed, time used:\t%.4f s'% (tm))
        score = metrics.accuracy_score(self.y_test, y_pred)
        print("Testing accuracy:\t%0.4f" % score)
    #    p, r, f1, _ = metrics.precision_recall_fscore_support(y_test, y_pred)
#        print(metrics.classification_report(self.y_test, y_pred))
        if plot_mat:
            plot_confusion_matrix(self.y_test, y_pred, figsize=(10, 10),#classes=class_names,
                                  text_on=False,
                                  title='Normalized confusion matrix (%s)'% name)
            plt.show()
        self.scores[name]['acc_test'] = score
        self.scores[name]['time_test'] = tm
        if return_pred:
            return model.predict_log_proba(self.X_test)
    
    
    def train_test(self, model, name = 'temp_model', save=True,
                   return_model = True, return_pred = False):
        
        model = self.train(model, name = name, save=save)        
        y_pred = self.test(model, return_pred = True, name = name)
    
        if return_pred:
            return model, y_pred
        else:
            return model
    
    def score_df(self, ):
        return pd.DataFrame(self.scores)
    
    
    def save_scores(self, fname='clf_scores', resdir = None, **kwds):
        
        resdir = self.resdir if resdir is None else resdir
        self.score_df().to_csv(resdir / f'{fname}.csv', index=True, header=True)
    
    
    
    