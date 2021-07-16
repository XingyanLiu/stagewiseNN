# -*- coding: utf-8 -*-
"""
Created on Sat Sep 14 23:01:05 2019

@author: xyliu
"""

import os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from scipy.spatial.distance import squareform, pdist, cdist
from scipy import sparse

from matplotlib import pyplot as plt

import scanpy as sc
from sklearn.preprocessing import StandardScaler

# In[]
''' preprocessing
'''


def set_hvgs(adata, gene_list, slim = False):
    '''
    this step does not mean that the given set of genes are real HVGs
    '''
    print('Setting the given set of %d genes as highly variable'%len(gene_list))
    indicator = [g in gene_list for g in adata.var_names]
    adata.var['highly_variable'] = indicator
    if slim:
#            if not adata.raw:
        print('slimming adata to contain only HVGs')
        return adata[:, indicator]
    return adata


def normalize_reverse(adata, counts_level=None, frac=0.05, setraw=True, copy=False):
    adata = adata.copy() if copy else adata
    sc.pp.log1p(adata,)
    sc.pp.normalize_total(adata, target_sum = counts_level, max_fraction=frac)
    adata.raw = adata
    return adata if copy else None




# In[]
'''
    Dimension reduction (mainly for integration)
    - Partial PCA
    - DCCA (PLSSVD)
    - Partial UMAP
'''


def myPCA(adata, key_add = 'X_pca', n_comps=50, 
           copy=False, **kwds):
    ''' this function was not used, just for debugging
    '''
    X = adata.X 

    print('Computing My PCA')
    from sklearn.decomposition import PCA
    pca = PCA(n_components=n_comps, **kwds)
    X_pca = pca.fit_transform(X)
    
    if copy:
        adata = adata.copy()
        
    adata.obsm[key_add] = X_pca
    adata.uns['pca'] = {}
    adata.uns['pca']['variance'] = pca.noise_variance_ 
    adata.uns['pca']['variance_ratio'] = pca.explained_variance_ratio_ 
    adata.varm['PCs'] = pca.components_.T
    
    return adata if copy else None



def partialPCAdt(adata, based_on, key_add = 'X_pca', n_comps=50, 
           copy=False, **kwds):
    '''
    adata.X.shape: n_cells x n_genes
    based_on: a tuple, or a list of length 2, e.g. (key, class_name)
        based_on[0]: key for the class labels
        based_on[1]: class name of the subset of X to fit the PCA
    '''
    labels = adata.obs[based_on[0]]
    X = adata.X 
    labels = pd.Series(labels)
    ind1 = (labels == based_on[1])
    print('Computing PCA on class `{}` grouped by `{}`'.format(*based_on))
    
    X1 = X[ind1, :]
    X2 = X[~ ind1, :]
    X1_pca, X2_pca, pca = partialPCA(X1, X2, n_comps=n_comps, **kwds)
    if copy:
        adata = adata.copy()

    X_pca = np.zeros((X.shape[0], n_comps))
    X_pca[ind1, :] = X1_pca
    X_pca[~ ind1, :] = X2_pca
    adata.obsm[key_add] = X_pca
    adata.uns['pca'] = {}
    adata.uns['pca']['variance'] = pca.noise_variance_ 
    adata.uns['pca']['variance_ratio'] = pca.explained_variance_ratio_ 
    adata.varm['PCs'] = pca.components_.T
    
    return adata if copy else None


def partialPCA(X1, X2, n_comps=2, **kwds):
    '''
    fit PCA on X, and transform both X and Y based on the fitted components
    X1: array-like, shape (n_samples1, n_features)
    X2: array-like, shape (n_samples2, n_features)
    
    return
    X_pca: array-like, shape (n_samples1 + n_samples2, n_features)
    estimator: `PCA` object
    
    '''
    from sklearn.decomposition import PCA
    pca = PCA(n_components=n_comps, **kwds)
    X1_pca = pca.fit_transform(X1)
    X2_pca = pca.transform(X2)
        
    return X1_pca, X2_pca, pca


def partialUMAPAdt(adata, based_on, n_neighbors=20, 
                   metric = 'cosine',
                   key_add = 'X_umap', 
                   copy=False, **kwds):
    '''
    adata.X.shape: n_cells x n_genes
    based_on: a tuple, or a list of length 2, e.g. (key, class_name)
        based_on[0]: key for the class labels
        based_on[1]: class name of the subset of X to fit the UMAP
    '''
    X = adata.X
    labels = adata.obs[based_on[0]]
    train_on = based_on[1]
    print('Computing partial-UMAP on class `{}` grouped by `{}`'.format(*based_on))
    X_umap = partialUMAP(X, labels, train_on, metric=metric,
                         n_neighbors=n_neighbors, **kwds)
    adata.obsm[key_add] = X_umap
    if copy:
        adata = adata.copy()
    return adata if copy else None



def partialUMAP(X, labels, train_on,
                n_neighbors=20, metric = 'euclidean', 
                n_comps = 2, min_dist = 0.5):
    '''
    fit UMAP on X, and transform both X and Y based on the fitted components
    X1: array-like, shape (n_samples1, n_features)
    X2: array-like, shape (n_samples2, n_features)
    
    return
    X_umap: array-like, shape (n_samples1 + n_samples2, n_features)   
    '''
    labels = pd.Series(labels)
    ind_tr = (labels == train_on)
    X0 = X[ind_tr, :]
    X1 = X[~ ind_tr, :]
    X0_umap, X1_umap = partial_umap(X0, X1, n_comps=n_comps,
                                    n_neighbors=n_neighbors, 
                                    min_dist = min_dist, metric = metric)
    
    X_umap = np.zeros((X.shape[0], 2))
    X_umap[ind_tr, :] = X0_umap
    X_umap[~ ind_tr, :] = X1_umap
    print('Combined !')
    return X_umap


def partial_umap(X_train, X_test, 
                 n_comps = 2,
                 n_neighbors=20, 
                 min_dist = 0.5, metric = 'euclidean'):
    from umap import UMAP
    mapper = UMAP(n_neighbors=n_neighbors, n_components=n_comps,
                  min_dist = min_dist, metric = metric)
    print('Mapper constructed:\n', mapper)
    X_train_umap = mapper.fit_transform(X_train)
    print('UMAP Mapper fitted !')
    X_test_umap = mapper.transform(X_test)
    print('Teast data transformed into the low dimensional space !')

    return X_train_umap, X_test_umap



def ccaAdt(adata, groupby, key_add = 'X_cca', n_comps = 50, 
           method = 'plssvd', copy = False, **kwds):
    '''
    
    testing code:
    B.ccaAdt(ph.adata, groupby='RNAcap', n_comps=2, copy=True)
    '''
    labels = adata.obs[groupby]
    X = adata.X 
    print('Computing cross-decomposition on 2 groups by `%s`'% groupby)
    X_cca, CCs = cca(X, labels, n_comps=n_comps, method=method, return_info=False, **kwds)
    if copy:
        adata = adata.copy()
    adata.obsm[key_add] = X_cca
    for lb in CCs.keys():
        adata.varm['CCs_' + lb] = CCs[lb]
    
    return adata if copy else None

def cca(X, labels, n_comps = 2, method = 'plssvd', return_model=False, scale=False, **kwds):
    '''
    inputs
    ======
    Suppose that n_cells = p + q.
    
    X : array-like, shape = [n_cells, n_genes] (will be transposed during computation)
    labels : array-like, shape = [n_cells, 1]
    scale : bool, False by default !
    
    return
    =======
    X_cca : array, [n_cells, n_components] (i.e. [p + q, n_components])
    x_weights_ : array, [p, n_components]
    y_weights_ : array, [q, n_components]
    x_scores_ : array, [n_genes, n_components]
    y_scores_ : array, [n_genes, n_components]
    
    if method is not 'plssvd':
        x_loadings_ : array, [p, n_components]
        y_loadings_ : array, [q, n_components]
    
    NOTE:
        * `n_features` here means `n_cells`
        * `n_samples here means `n_genes`
        * len(pd.unique(labels)) should be 2

    '''
    ## data
    n_cells = X.shape[0]
    X = X.T # should be transposed
    labels = pd.Series(labels)
    unique_lbs = labels.unique()
    ind1 = (labels == unique_lbs[0])
    print(f'unique labels: {unique_lbs}')
    ## model construction
    from sklearn.cross_decomposition import PLSSVD, CCA, PLSCanonical
    print('Using method [ %s ], with [ %d ] components' % (method, n_comps))
    if method == 'plssvd':
        md = PLSSVD(n_components=n_comps, scale=scale, **kwds)
        md.fit(X[:, ind1], X[:, ~ind1])
        X1_, X2_ = md.x_weights_,  md.y_weights_
    elif method == 'cca': # NOT applicatable !!!
        md = CCA(n_components=n_comps, scale=scale, **kwds)
        md.fit(X[:, ind1], X[:, ~ind1])
        X1_, X2_ = md.x_loadings_,  md.y_loadings_
#        X1_, X2_ = md.x_rotations_,  md.y_rotations_
    elif method == 'plsc':
        md = PLSCanonical(n_components=n_comps, scale=scale, **kwds)
        md.fit(X[:, ind1], X[:, ~ind1])
        X1_, X2_ = md.x_weights_,  md.y_weights_
    else:
        raise ValueError('the argument `method` should be either `plssvd` or `cca`')

    X_cca = np.zeros((n_cells, n_comps))
    X_cca[ind1, :] = X1_
    X_cca[~ ind1, :] = X2_
    
    if return_model:
        return X_cca, md
    else:
        CCs = {str(unique_lbs[0]): md.x_scores_,
               str(unique_lbs[1]): md.y_scores_}
        return X_cca, CCs
    
    
# In[]
'''     z-score grouped by ...
'''


def Zscore(X, with_mean = True, scale=True, asdf=False):
    ''' For each column of X, do within-group centering (z-scoring)
    ====
    code borrowed from `scanpy.pp._simple`
    '''
    scaler = StandardScaler(with_mean=with_mean, copy=True).partial_fit(X)
    # user R convention (unbiased estimator)
    if scale:
        scaler.scale_ *= np.sqrt(X.shape[0]/(X.shape[0]-1))
    else:
        scaler.scale_ = 1
    X_new = scaler.transform(X)
    
    if isinstance(X, pd.DataFrame):
        X_new = pd.DataFrame(X_new, index=X.index, columns=X.columns)

    return X_new


def GroupZscore(X, labels, with_mean = True, scale=True, max_value = None):
    '''
    For each column of X, do within-group centering (z-scoring)
    ======
    X: np.array, shape (n_samples, n_features)
        i.e. each row of X is an observation, wile each column is a feature
        
    with_mean: boolean, True by default
        If True, center the data before scaling, and X shoud be a dense matrix.
    '''
    
    X = X.astype(np.float).copy()
    labels = np.asarray(labels)
    unique_labels = np.unique(labels)
    for lb in unique_labels:
        ind = labels == lb
        if sum(ind) == 1:
            print(f'ignoring class {lb} with only one sample.')
            continue
        print(lb)
        X[ind, :] = Zscore(X[ind, :], with_mean=with_mean, scale=scale)
        
    if max_value is not None: 
        X[X > max_value] = max_value
        print('... clipping at max_value', max_value)
    return X


def GroupZscoreAdt(adt, key='counts', groupby='batch', key_new = None, 
                   max_value = None, 
                   with_mean = True, 
                   cover = True, **kwds):
    '''
    adt: AnnData
    key: str in ['X_pca', 'count']
        can be a key from adt.obsm, e.g. `key='X_pca'`
        If key == 'counts', then do scaling on `adt.X` 
        and cover the old count matrix, ignoring the `cover` parameter
    groupby: str; 
        A key from adt.obs, from which the labels are take
    cover: bool
        whether to cover the old X with the scored X
    '''
    labels = adt.obs[groupby]
    if key == 'counts':
        print('doing z-score scaling on count matrix, transformed into a dense array')
        if sparse.issparse(adt.X): adt.X = adt.X.toarray()
        if not cover: adt = adt.copy() 
        adt.X = GroupZscore(adt.X, labels, with_mean = with_mean, 
                            max_value = max_value, **kwds)
#        print('TEST:, adt.X.shape)
    else:
        if cover: 
            key_new = key
        else: 
            key_new = key + '_new' if key_new is None else key_new
        adt.obsm[key_new] = GroupZscore(adt.obsm[key], labels, \
                with_mean=with_mean, max_value = None, **kwds)
#        print('TEST- 1: ', adt.X.shape)
    return adt if not cover else None   


def WrapperScale(adata, zero_center=True, max_value=None,
                 groupby=None, copy=False, **kwds):
    
    '''
    Wrapper function for centering and scaling data matrix `X` in sc.AnnData,
    extended for within-batch cprocessing.
    
    Example
    =======
        WrapperScale(adata, groupby='batch')
    '''
    if groupby is not None:
        print('doing within-group scaling, group by [ %s ]'%groupby)
        GroupZscoreAdt(adata, key = 'counts', 
                       max_value=max_value,
                       groupby = groupby, 
                       with_mean = zero_center,
                       cover = not copy, **kwds) 
    else:
        print('using the build-in function `sc.pp.scale(..)`')
        sc.pp.scale(adata, zero_center=zero_center, 
                    max_value=max_value, copy=copy)
        
        
def normalize_default(adata, target_sum=1e4, copy=False, log_only=False,
                      force_return=False):
    if copy:
        adata = adata.copy()
    print('normalizing datasets with default settings.')
    if not log_only:
        print(f'performing total-sum normalization, target_sum={target_sum}...')
        sc.pp.normalize_total(adata, target_sum=target_sum)
    else:
        print('skipping total-sum normalization')
    sc.pp.log1p(adata)
    return adata if copy or force_return else None
        
# In[]
'''     recurrent clustering 
==========================================================
    to break the whole dataset into many small clusters
'''


def recurrent_clust(adata, 
                    key_new='recurrent_leiden', 
                    key_orig=None,
                    max_cnt = 50,
                    reso = 1,
                    rename = True,
                    figdir='.'):
    if isinstance(key_orig, str) and key_orig in adata.obs.keys():
        adata.obs[key_new] = adata.obs[key_orig]
    else:
        sc.tl.leiden(adata, resolution=reso, key_added = key_new)
    
    group_cnts = adata.obs[key_new].value_counts()
    while True:
        candi_grp = group_cnts.idxmax()
        _n_cells = group_cnts[candi_grp]
        print(_n_cells)
        if _n_cells > max_cnt:
            print(f'breaking cluster `{candi_grp}`...')
            sc.tl.leiden(adata, resolution=reso,
                         restrict_to=(key_new, [candi_grp]),
                         key_added=key_new)
            group_cnts = adata.obs[key_new].value_counts()
        else:
            print(f'no cluster contain more than {max_cnt} cells, stop.')
            break
    n_groups = len(group_cnts)
    _hist_group_counts(group_cnts, figdir = figdir)
    print(f'resulting {n_groups} clusters in total.')
    
    if rename:
        new_cats = list(map(str, range(n_groups)))
        adata.obs[key_new].cat.categories = new_cats
    
    return adata



def _hist_group_counts(group_counts: pd.Series, figdir=None):
    tt = 'number of cells in each aggregation'
    group_counts.hist()
    plt.title(tt)
    plt.xlabel('number of cells')
    plt.ylabel('frequency')
    if figdir is not None:
        plt.savefig(Path(figdir) / f'{tt}.pdf', bbox_inches='tight')
    
    


    
# In[]


def RemoveGroups(adata, group_names, key = None):
    key = _auto_decide_key(adata, key)
#    if isinstance(group_names, (str, int)):
#        group_names = [group_names]
#    
#    labels = adata.obs[key]
#    group_names = _check_type(group_names, ref = labels[0])
#    indicators = [lb not in group_names for lb in labels]
    indicators = take_groups(adata.obs, group_names, key, indicate=True, remove=True)
    return adata[indicators, :].copy()

        


def TakeGroups(adata, group_names, key = None, onlyX=False, copy=False):
    '''key = 'leiden' by default
    '''
    key = _auto_decide_key(adata, key)
    indicators = take_groups(adata.obs, group_names, key, indicate=True)
    if copy:
        return adata[indicators, :].X.copy() if onlyX else adata[indicators, :].copy() 
    else:
        return adata[indicators, :].X if onlyX else adata[indicators, :]


def take_groups(df, group_names, col=None, indicate=False, remove=False):
    '''
    df: pd.Series or pd.DataFrame;
        if Series: ignore `col`
    group_names: names of groups that you want to take out
    col: str, column name in df.columns
    indicate: bool; 
        if True, return a Series of bool indicators of the groups
        else, return a DataFrame 
    '''
#    if isinstance(group_names, tuple):
#        group_names = list(group_names)
    if isinstance(group_names, (str, int)):
        group_names = [group_names]
    if isinstance(df, pd.DataFrame):
        labels = df[col]
    else:
        labels = df
    group_names = _check_type(group_names, ref = labels[0])
#    indicators = [lb in group_names for lb in labels]
    if remove:
        indicators = labels.apply(lambda x: x not in group_names).to_list()
    else:
        indicators = labels.apply(lambda x: x in group_names).to_list()
    if indicate:
        return indicators
    else:
        return df.loc[indicators, :]
    
    
def label_mask_others(df, key, keeps, key_new='tmp', lb_masked='others', copy=False):
    ''' create a new column for highlighting given groups'''
    if copy:
        df = df.copy()
    df[key_new] = df[key].apply(lambda x: x if x in keeps else lb_masked)
    
    return df if copy else None


def MergeGroups(adata, key, group_lists, new_key = None, 
                rename=False, copy=False):
    '''
    merge the given groups into one single group
    which is named as '_'.join(groups[i]) by default
    `group_lists`: a list of lists of group names
    
    === TEST ===
    >>> MergeGroups(adata, 'batch', [list('AB'), list('EF')], copy=True)
    >>> MergeGroups(adata, 'batch', [list('AB'), list('EF')],)
    >>> adata
    '''
    adata = adata.copy() if copy else adata
    labels0 = adata.obs[key].copy()
#    group_lists = [group_lists] if not isinstance(group_lists[0], list) else group_lists
    labels = merge_groups(labels0, group_lists)
    
    new_key = key + '_new' if new_key is None else new_key
    print(pd.value_counts(labels))
    adata.obs[new_key] = pd.Categorical(labels)
    print(f'A new key `{new_key}` added.')
    if rename:
        old_cats = adata.obs[new_key].cat.categories
        n_groups = len(old_cats)
        new_cat = list(map(str, range(n_groups)))
        adata.obs[new_key].cat.categories = new_cat
        print(f'original names:\n{old_cats}')
        print(f'categories are renamed as:\n{new_cat}')
    if copy: return adata


def merge_groups(labels, group_lists):
    '''
    labels: better be a pd.Series object.
    group_lists: a list of lists of group names to be merged
    
    === TEST ===
    >>> merge_groups(adata.obs['batch'], [list('AB'), list('EF')]).unique()
    '''
    if not isinstance(labels, pd.Series):
        labels = pd.Series(list(map(str, labels)))
    
    if not isinstance(group_lists[0], list): 
        ## merge only one set of groups
        groups = list(map(str, group_lists))
        new_name = '_'.join(groups)
        labels = labels.apply(lambda x: new_name if x in groups else x)
        print('groups are merged into a new single group: ', new_name)
#        print(pd.value_counts(labels))
        return labels
    else:
        ## merge multiple sets of groups
        for groups in group_lists:
            new_name = '_'.join(groups)
            labels = merge_groups(labels, groups)
        
        return labels


def MergeByMarkers(adata, key, ):
    
    
    pass

def FindShortcuts(adata, key = None, n=1):
    key = _auto_decide_key(adata, key)

    similarities = GroupSimilarityPaga(adata, key=key, copy=True)
    scores = similarities.sum(axis=0)
    ids = np.argsort(scores)[- n: ]
    
    return ids

# In[]
    
'''     Functions for matching groups (clusters)

    - PAGA
    - Marker
    - correlation (mean expression)
'''
def GroupSimilarityPaga(adata, key=None, copy=True):
    '''
    if copy == True,  run paga on a copy 
    
    return an array
    '''
    key = _auto_decide_key(adata, key)
    adata = sc.tl.paga(adata, groups=key, copy=copy)
    similarities = adata.uns['paga']['connectivities'].toarray()
    
    return similarities


def MatchByMarkersAdt(adata, markers=None, ntop=10):
    '''
    Jaccard Index between corresponding marker sets
    return an array of dimension M*M for M clusters (or groups)
    '''
    if markers is None:
        markers = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(ntop)
    
    score_mat = MatchByMarkers(markers)
    return score_mat


def MatchByMarkers(markers, weighted = False):
    '''
    `markers`:
        a list of M lists of markers, each list corresponds a group;
        or a n-by-M DataFrame, with n markers for each of the M groups;
        
    (Weighted) Jaccard Index between corresponding marker sets
    return an array of dimension M*M for M groups
    '''
    ## format the inputs
    if isinstance(markers, np.recarray):
        markers = pd.DataFrame(markers)
    if isinstance(markers, pd.DataFrame):
        markers = markers.values.T.tolist()
    markers_coded = Encode(markers, weighted = weighted, sparse=False)
    
    ## pair-wise weighted Jaccard similarities
    match_scores = pdist(markers_coded, metric=JaccardIndex, )
    match_scores = squareform(match_scores)
    return match_scores

#    M = len(markers)
#    match_scores = np.zeros(shape = (M, M))
#    for i1 in range(M):
#        for i2 in range(i1 + 1, M):
#            match_scores[i1, i2] = JaccardIndex(markers[i1], markers[i2], weighted)
#    match_scores += match_scores.T

def Encode(X, weighted=False, sparse = True):
    '''
    Transforms lists of feature-value mappings to vectors.
    X:
        a list of sets or dict
    '''
    if not isinstance(X[0], (dict)):
        X = [{x: 1 for x in xc} for xc in X]
#        X = [dict([[x, 1] for x in xc]) for xc in X]
    from sklearn.feature_extraction import DictVectorizer
    code = DictVectorizer(sparse=sparse).fit_transform(X)
    return code    
#tmp = Encode(markers)

def JaccardIndex(v1, v2):
    '''
    Weighted Jaccard Index
    v1, v2:
        2 numeric vectors of the same length
    return:
        Jaccard Index (similarity) between v1 and v2, which ranges
        from 0 (totally different) to 1 (exactly the same)
    '''
    if sparse.issparse(v1):
        stacked = sparse.vstack([v1, v2])
        intersect = stacked.min(axis = 0)
        union = stacked.max(axis = 0)
    else:
        intersect = np.sum(np.min([v1, v2], axis = 0))
        union = np.sum(np.max([v1, v2], axis = 0))
    return intersect / union

def jaccard_index(set1, set2, denominator='union'):
    '''
    Compute Jaccard Index for two sets, which are NOT required 
    to be of the same lengths
    '''
    set1, set2 = tuple(map(set, (set1, set2)))
    intersect = set1.intersection(set2)
    union = set1.union(set2)
    if denominator == 'union':
        d = len(union)
    elif denominator == 'min':
        d = min([len(set1), len(set2)])
    elif denominator == 'max':
        d = max([len(set1), len(set2)])
    elif isinstance(denominator, int):
        d = denominator
    jcd = len(intersect) / d
    return jcd
    
    


def GroupSimilarityMean(adatas, keys, use_genes=None,
                        metric='spearman', output_means=False,
                        use_raw=True, tags=None, ):
    '''
    Computing similarities between clusters in a list of adatas, including
    * clusters within each adata, and
    * clusters from different adatas
    ----------
    Parameters
    ----------
    adatas: sc.AnnData object or a list of AnnDatas
    keys: Union[str, list[str]], i.e. str or a list of strings
    use_genes: list[str], list of gene names
    metric: Union[str, map]
    
    
    '''
    if isinstance(adatas, sc.AnnData): # single case
        key = keys if isinstance(keys, str) else keys[0]
        gmeans = GroupMean(adatas, groupby=key, features=use_genes, use_raw=use_raw)
        clnames = gmeans.columns
    elif hasattr(adatas, '__iter__'):
        n = len(adatas)
        if isinstance(keys, str):
            keys = [keys] * n
        tags = list(map(str, np.arange(n))) if tags is None else tags
        # within-group averages 
        gmeans = [GroupMean(adt, ky, use_genes, use_raw=use_raw) for adt, ky in zip(adatas, keys)]
        gmeans = pd.concat(gmeans, axis=1, keys=tags, sort=False)
        gmeans.fillna(0, inplace=True)
        # Make names
        clnames = ['_'.join(c) for c in gmeans.columns.to_list()]
    
    gmeans.set_axis(clnames, axis=1, inplace=True)
    # Calculates the similarities
#    sim = gmeans.corr(method=method)
    sim = similarity(gmeans, metric=metric, axis=1, name=False)
    sim = pd.DataFrame(sim, index=clnames, columns=clnames)
#    return sim
    return sim, gmeans if output_means else sim



def similarity(X, Y=None, metric='spearman', axis=0, name=True):
    '''
    X: array like, shape (n_samples, n_features) for axis=0
    Y: array like, shape (n_samples1, n_features) for axis=0
        - if Y is None, a symetric matrix is returned, whose shape is 
        (n_samples, n_samples)
        - else, a similarity matix of shape (n_samples, n_samples1) is returned.
    metric: 'spearman', 'correlation', 'cosine' are recommended
    '''
    if axis == 1:
        X = X.T
        if Y is not None: Y = Y.T
        
    print(f'Caculating similarities using matric {metric}...')
    print(f'shape of X: {X.shape}')
#    print(f'shape of Y: {Y.shape}')
    if metric == 'spearman':
        sim, pvals = spearmanr(X, Y, axis=1)
        if Y is not None:
            nx, ny = X.shape[0], Y.shape[0]
            sim = sim[:nx, -ny:]
    else:
        if Y is None:
            dist = squareform(pdist(X, metric=metric, ))
            sim = 1 - dist
        else:
            dist = cdist(X, Y, metric=metric)
            sim = 1 - dist
    print('shape: ', sim.shape)
    print('Done!')
    
    if name and hasattr(X, 'index'): # just keep the names of the dataframe
        index = X.index
        if Y is not None:
            columns = Y.index if hasattr(Y, 'index') else None
        else: columns = X.index
        sim = pd.DataFrame(sim, index=index, columns = columns)
        
    return sim

def similarity_corr(df, method='pearson', center=False):
    ''' correlations among columns
    '''
    if center:
        df = Zscore(df.T, asdf=True).T
    sims = df.corr(method)
    return sims


def concat_columns(dfs, tags=None, sort=False, **kw):
    '''
    a wwraper of `pd.concat()`, but set new column names
    dfs: a list of `pd.DataFrame`
    '''
    n = len(dfs)
    tags = list(map(str, np.arange(n))) if tags is None else tags
    df = pd.concat(dfs, axis=1, keys=tags, sort=sort, **kw)
    df.fillna(0, inplace=True)
    # Make names
    clnames = ['_'.join(c) for c in df.columns.to_list()]
    df.set_axis(clnames, axis=1, inplace=True)
    return df


    
def GroupMeanMultiData(adatas, keys, use_genes=None, binary=False,
                        use_raw=True, tags=None, ):
    n = len(adatas)
    if isinstance(keys, str):
        keys = [keys] * n
    # within-group averages 
    gmeans = [GroupMean(adt, ky, use_genes, binary=binary, use_raw=use_raw)\
              for adt, ky in zip(adatas, keys)]
    tags = list(map(str, np.arange(n))) if tags is None else tags
    gmeans = pd.concat(gmeans, axis=1, keys=tags, sort=False)
    gmeans.fillna(0, inplace=True)
    # Make names
    clnames = ['_'.join(c) for c in gmeans.columns.to_list()]
    gmeans.set_axis(clnames, axis=1, inplace=True)
    return gmeans
    
    

def GroupMean(adata, groupby, features=None, binary=False, 
              use_raw=False, **kwds):
    '''
    groupby: a column name in adata.obs
    features: a subset of names in adata.var_names (or adata.raw.var_names)
    
    return
    ------
    a pd.DataFrame with features as index and groups as columns
    '''
    labels = adata.obs[groupby].values
    print(f'Computing averages grouped by {groupby}')
    if use_raw and adata.raw is not None:
        if features is not None:
            features = [f for f in features if f in adata.raw.var_names]
            X = adata.raw[:, features].X
        else:
            features = adata.raw.var_names
            X = adata.raw.X
    else:
        if features is not None:
            features = [f for f in features if f in adata.var_names]
            X = adata[:, features].X
        else:
            features = adata.var_names
            X = adata.X
    return group_mean(X, labels, binary=binary, features=features, **kwds)


def group_mean(X, labels, binary=False, classes=None, features=None,
               print_groups=True):
    '''
    This function may work better than `df.groupby().mean()` when dealing
    with sparse matrix. (obviously~)
    ---
    X: shape (n_samples, n_features)
    labels: shape (n_samples, )
    classes: optional, names of each group
    features: optional, names of features
        
    '''
    classes = np.unique(labels, ) if classes is None else classes
    if binary:
        X = (X > 0)#.astype('float')
        print('Binarized...the results will be the expression proportions.')
        
    if len(classes) == 1:
        grp_mean = X.mean(axis = 0).T
    else:
        from sklearn.preprocessing import label_binarize
        lb1hot = label_binarize(labels, classes=classes, sparse_output=True)
        if len(classes) == 2:
            lb1hot = lb1hot.toarray()
            lb1hot = np.c_[1 - lb1hot, lb1hot]
        print(f'Calculating feature averages {len(classes)} groups:')
        if print_groups:
            print(classes)
#        if binary:
#            X = (X > 0)#.astype('float')
#            print('Binarized...the results will be the expression proportions.')
        grp_mean = X.T.dot(lb1hot) / lb1hot.sum(axis=0)
    grp_mean = pd.DataFrame(grp_mean, columns = classes, index=features)
    return grp_mean


def moving_average(a, n, pad=True):
    """Moving average over one-dimensional array.

    Parameters
    ----------
    a : np.ndarray
        One-dimensional array.
    n : int
        Number of entries to average over. n=2 means averaging over the currrent
        the previous entry.

    Returns
    -------
    An array view storing the moving average.
    """
    N = len(a)
    if pad:
        if n % 2:
            pad_width = (n - 1) // 2 
        else:
            pad_width = (n - 1) // 2 + 1
        a = np.pad(a, pad_width, mode='edge')
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1: N + n - 1] / n


# In[]
    
def flattened_frequency(df, n=None):
    freq = pd.value_counts(df.values.flatten())
    return freq


def TopMarkers(adata, n = 5, groups=None, unique = True,):
    df = marker_table(adata)
    return top_markers(df, n=n, groups=groups, unique=unique)


def top_markers(marker_df, n = 5, groups=None, unique = True, ):
    '''
    marker_df: a data-frame with cluster names as columns, and genes as values
    groups: a list of cluster names (column names)
    return a flattened marker list
    '''
    groups = marker_df.columns if groups is None else groups
    top = marker_df[groups].iloc[: n].values.T.flatten()
    print('shape: ', top.shape)
    if unique:
        top = pd.unique(top)
        print('shape (unique): ', top.shape)
    return top

def marker_table(adata):
    return pd.DataFrame(adata.uns['rank_genes_groups']['names'])

def marker_info(adata):
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    df = pd.DataFrame(
            {group + '_' + key: result[key][group]
            for group in groups for key in ['names', 'logfoldchanges', 'pvals', 'pvals_adj', 'scores']})
        
    return df

def top_markers_from_srt(mktb, n=20, 
                        groupby = 'cluster',
                        col_take = 'gene',
                        col_sort = 'p_val_adj',
                        ):
    '''handling marker table exported from Seurat.
    mktb: marker table, pd.DataFrame
    '''
    ascending = True
    if col_sort in ('avg_logFC', 'pct.1'):
        ascending = False
        
    _top = mktb.sort_values(col_sort, ascending=ascending)\
        .groupby('cluster')['gene'].apply(lambda x: x.head(n))
    _top = pd.unique(_top) # an ndarray returned
    print('shape (unique): ', _top.shape)
    return _top


def save_embeddings(adata, path, key=None, ndim=2, tail='', ):
    path = Path(path) if not isinstance(path, Path) else path
    if key is None:
        key = 'X_umap' if 'X_umap' in adata.obsm.keys() else 'X_pca'
    if not path.is_file(): 
        path = path / ('%s_%s.csv'% (key, tail))
    
    pd.DataFrame(adata.obsm[key][:, : ndim]).to_csv(path, header=False, index=False)
    print('Embeddings (%s) are saved into: %s'%(key, path))


def cross_labeling(df, keys, fmt=None, new_key=None, inplace=False):
    '''
    Adding a column as combined labels from other columns
    keys: a list of column names of the df
    fmt: format for the new class-names
    inplace:
        if true: a new column within name `new_key` will be added
        else: a Serise will be returned
    '''
#    df = df.copy() if copy else df

    if fmt:
        mapping = lambda x: fmt.format(*x) 
    else:
        mapping = lambda x: '_'.join(map(str, x))
    new_col = df[keys].apply(lambda x: mapping(x), axis=1)
    
    if inplace:
        new_key = '_'.join(keys) if new_key is None else new_key
        df[new_key] = pd.Categorical(new_col)
        print(f'new key `{new_key}` is added.')
    else:
        return new_col
    

# In[]

        
def ResetRaw(adata, ):
    ''' 
    Re-set from the (high dimensional) raw matrix, 
    while keeping the annotations.
    '''    
#    len(ph.adata.raw.var_names)
    _adata = sc.AnnData(X = adata.raw.X.copy(), 
                        obs = adata.obs.copy(), 
                        var = adata.raw.var.copy(), 
                        )
    print('Constructing a new AnnData object...')
    print(_adata)
    
    return _adata


def ResetX(adata, copy = True):
    ''' Generally used for replacing the scaled matrix by the count data.
    Dimension of columns (genes) are not changed.
    Assume that adata.raw is not None
    '''
    _adata = adata.copy() if copy else adata
    var_ids = np.nonzero([g in adata.var_names for g in adata.raw.var_names])[0]
    _adata.X = _adata.raw.X[:, var_ids].copy()
    print('Replacing the scaled matrix by the count data. Dimension of columns (genes) are not changed.')
    print(_adata)
    return _adata if copy else None

# In[]




def subsample(data, prop, groupby=None, random_state=0):
    '''
    data : :class:`~anndata.AnnData`, `np.ndarray`, `sp.sparse`
    '''
    old_n_obs = data.n_obs if isinstance(data, AnnData) else data.shape[0]
    
    from sklearn.utils.random import sample_without_replacement
    sample_without_replacement(100, 10)
    pass

# In[]
    
def rmv_mito(gene_list, tag='mt-'):
    
    remained = [g for g in gene_list if not g.lower().startswith('mt-')]
    return remained

def normalize_col(X, scale_factor = 1, ):
    '''
    make the column elements of X to unit sum
    scale_factor: numeric, None
        if None: use the median of sum level as the scaling factor.
    '''
    sums = X.sum(axis=0)
    is_zero = sums == 0
    if scale_factor is None:
        scale_factor = np.median(sums[~ is_zero])
    sums /= scale_factor
    print(sums.shape)
    # for those rows or columns that summed to 0, just do nothing
    sums[sums == 0] = 1
    norm_ = 1 / sums
    if hasattr(norm_, 'A1'): norm_ = norm_.A1
    
    if sparse.isspmatrix(X):
        print('sparse normalization')
        X_new = X.dot(sparse.diags(norm_))
    else:
        print('dense normalization')
        X_new = X.dot(np.diag(norm_))
        
    if isinstance(X, pd.DataFrame):
        X_new.columns = X.columns
    return X_new 

def col_normalize_and_log1p(X, scale_factor = 1):
    ''' column normalization and log1p
    '''
    X_new = normalize_col(X, scale_factor=scale_factor)
    if isinstance(X, pd.DataFrame):
        X_new = X_new.apply(np.log1p)
    elif sparse.issparse(X_new):
        np.log1p(X_new.data, out=X_new.data)
    else:
        X_new = np.log1p(X_new)
    return X_new

def dist_Spearman(u, v, axis=0):
    '''
    axis : int or None, optional
        If axis=0 (default), then each column represents a variable, with
        observations in the rows. 
        If axis=1, the relationship is transposed: each row represents a 
            variable, while the columns contain observations.
        If axis=None, then both arrays will be raveled.
    '''
    rho, pval = spearmanr(u, v, axis=axis)
    return 1 - np.abs(rho)


def _auto_decide_key(adata, key=None, ):
    if key is None:
        key = 'leiden' if 'leiden' in adata.obs.keys() else 'louvain'
    return key

def _check_type(lst, ref):
    ''' make sure that each element in the list is string
    '''
    if type(lst[0]) != type(ref) and isinstance(ref, str):
        print('type transforming..')
        lst = list(map(str, lst))# [str(x) for x in lst]
    return lst


if __name__ == '__main__':
    

    WORKDIR = Path(r'D:\Users\xyliu\003')
    os.chdir(WORKDIR)
    from PipeHandler import PipeHandler 
    import funx as fx
    DATADIR_main = Path(r'E:\Users\xyliu\data')
    
    

    
#    vv = np.arange(18).reshape((3, 6))
#    
#    vvs = sparse.csr_matrix(vv)
#    ids1 = (vvs.sum(axis=1).A1 < 25)
#    ids2 = (vvs.sum(axis=0).A1 > 25)
#    break_links(vvs, [1], [0])
#    break_links(vvs, ids1, ids2).toarray()
    
    
    



















