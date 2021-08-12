# -*- coding: UTF-8 -*-
"""
@CreateDate: 2020/07/18
@Author: Xingyan Liu
@File: process.py
@Project: stagewiseNN
"""
import os
import sys
from pathlib import Path
from typing import Sequence, Mapping, Optional, Union, Callable
import logging
import pandas as pd
import numpy as np
import scanpy as sc
from scipy import sparse
from sklearn.preprocessing import label_binarize
from ._scale import wrapper_scale


def check_dirs(path):
    if os.path.exists(path):
        print('already exists:\n\t%s' % path)
    else:
        os.makedirs(path)
        print('a new directory made:\n\t%s' % path)


def reverse_dict(d: dict, ):
    """
    the values of the dict must be list-like type
    """
    d_rev = {}
    for k in d.keys():
        vals = d[k]
        _d = dict.fromkeys(vals, k)
        d_rev.update(_d)
    return d_rev


def describe_dataframe(df: pd.DataFrame, **kwargs):
    for c in df.columns:
        print(c.center(40, '-'), **kwargs)
        print(describe_series(df[c], asstr=True), **kwargs)


def describe_series(
        srs: Sequence, max_cats: int = 100,
        asstr: bool = False,
):
    """ inspect data-structure """
    srs = pd.Series(srs)
    if len(srs.unique()) <= max_cats:
        desc_type = 'Value counts'
        result = srs.value_counts(dropna=False)
    elif isinstance(srs[0], (int, float)):
        desc_type = 'Numerical summary'
        result = srs.describe()
    else:
        desc_type = 'Header lines'
        result = srs.head()
    if asstr:
        return f'{desc_type}:\n{result}'
    else:
        return result, desc_type


def make_binary(mat):
    mat_bin = mat.copy()
    mat_bin[mat_bin > 0] = 1.
    return mat_bin


def set_adata_hvgs(
        adata: sc.AnnData,
        gene_list: Optional[Sequence] = None,
        indicator: Optional[Sequence[bool]] = None,
        slim: bool = True,
        copy: bool = False,
):
    """
    Setting the given (may be pre-computed) set of genes as highly variable,
    if `copy` is False, changes will be made to the input adata.
    if slim is True and adata.raw is None, raw data will be backup.
    """
    if copy:
        adata = adata.copy()
    logging.info(
        'Setting the given set of %d genes as highly variable' % len(gene_list))
    if (indicator is None) and (gene_list is not None):
        indicator = [g in gene_list for g in adata.var_names]
    adata.var['highly_variable'] = indicator
    if slim:
        if adata.raw is None:
            adata.raw = adata
        logging.info('slimming adata to contain only HVGs')
        adata = adata[:, adata.var['highly_variable']].copy()
    return adata


def change_names(
        seq: Sequence,
        mapping: Optional[Mapping] = None,
        **kwmaps
) -> list:
    mapping = {} if mapping is None else mapping
    mapping.update(kwmaps)
    func = lambda x: mapping.get(x, x)
    return list(map(func, seq))


def normalize_default(
        adata, target_sum=None,
        copy=False, log_only=False,
):
    """Normalizing datasets with default settings (total-counts normalization
    followed by log(x+1) transform).

    Parameters
    ----------
    adata
        ``AnnData`` object
    target_sum
        scale factor of total-count normalization
    copy
        whether to copy the dataset
    log_only
        whether to skip the "total-counts normalization" and only perform
        log(x+1) transform

    Returns
    -------
    ``AnnData`` or None
    """
    if copy:
        adata = adata.copy()
        logging.info('A copy of AnnData made!')
    else:
        logging.info('No copy was made, the input AnnData will be changed!')
    logging.info('normalizing datasets with default settings.')
    if not log_only:
        logging.info(
            f'performing total-sum normalization, target_sum={target_sum}...')
        sc.pp.normalize_total(adata, target_sum=target_sum)
    else:
        logging.info('skipping total-sum normalization')
    sc.pp.log1p(adata)
    return adata


def normalize_log_then_total(
        adata, target_sum=None,
        copy=False,
):
    """ For SplitSeq data, performing log(x+1) BEFORE total-sum normalization
    will results a better UMAP visualization (e.g. clusters would be less
    confounded by different total-counts ).
    """
    if copy:
        adata = adata.copy()
        logging.info('A copy of AnnData made!')
    sc.pp.log1p(adata, )
    sc.pp.normalize_total(adata, target_sum=target_sum, )
    return adata


def groupwise_hvgs_freq(
        adata,
        groupby='batch',
        return_hvgs: bool = True,
        **hvg_kwds,
):
    """ Separately compute highly variable genes (HVGs) for each group, and
    count the frequencies of genes being selected as HVGs among those groups.

    Parameters
    ----------
    adata
        the ``AnnData`` object
    groupby
        a column name in ``adata.obs`` specifying batches or groups that you
        would like to independently compute HVGs.
    return_hvgs
        whether to return the computed dict of HVG-lists for each group
    hvg_kwds
        Other Parameters for ``sc.pp.highly_variable_genes``
    Returns
    -------
    hvg_freq: dict
        the HVG frequencies
    hvg_dict: dict
        returned only if ``return_hvgs`` is True
    """
    from collections import Counter
    hvg_dict = {}
    hvg_freq = Counter()
    group_labels = adata.obs[groupby]
    for g in group_labels.unique():
        _adt = adata[group_labels == g].copy()
        sc.pp.highly_variable_genes(_adt, **hvg_kwds)
        _hvgs = _adt.var[_adt.var['highly_variable']].index
        hvg_freq += Counter(_hvgs)
        if return_hvgs:
            hvg_dict[g] = list(_hvgs)

    hvg_freq = dict(hvg_freq)
    if return_hvgs:
        return hvg_freq, hvg_dict
    return hvg_freq


def take_high_freq_elements(
        freq: Mapping,
        min_freq: int = 3):
    return list(filter(lambda x: freq[x] >= min_freq, freq.keys()))


def set_precomputed_neighbors(
        adata,
        distances,
        connectivities,
        n_neighbors=15,
        metric='cosine',  # pretended parameter
        method='umap',  # pretended parameter
        metric_kwds=None,  # pretended parameter
        use_rep=None,  # pretended parameter
        n_pcs=None,  # pretended parameter
        key_added=None,  #
):
    if key_added is None:
        key_added = 'neighbors'
        conns_key = 'connectivities'
        dists_key = 'distances'
    else:
        conns_key = key_added + '_connectivities'
        dists_key = key_added + '_distances'

    if connectivities is None:
        connectivities = distances.copy().tocsr()
        connectivities[connectivities > 0] = 1

    adata.obsp[dists_key] = distances
    adata.obsp[conns_key] = connectivities

    adata.uns[key_added] = {}
    neighbors_dict = adata.uns[key_added]
    neighbors_dict['connectivities_key'] = conns_key
    neighbors_dict['distances_key'] = dists_key

    neighbors_dict['params'] = {'n_neighbors': n_neighbors, 'method': method}
    neighbors_dict['params']['metric'] = metric
    if metric_kwds is not None:
        neighbors_dict['params']['metric_kwds'] = metric_kwds
    if use_rep is not None:
        neighbors_dict['params']['use_rep'] = use_rep
    if n_pcs is not None:
        neighbors_dict['params']['n_pcs'] = n_pcs

    return adata


def quick_preprocess_raw(
        adata: sc.AnnData,
        target_sum: Optional[int] = None,
        hvgs: Optional[Sequence] = None,
        batch_key=None,
        copy=True,
        log_first: bool = False,
        **hvg_kwds
) -> sc.AnnData:
    """
    Go through the data-analysis pipeline, including normalization, HVG
    selection, and z-scoring (centering and scaling)

    Parameters
    ----------
    adata
        the ``Anndata`` object
    target_sum
        the target total counts after normalization.
        If `None`, after normalization, each observation (cell) has a total
        count equal to the median of total counts for observations (cells)
        before normalization.
    hvgs
        highly variable genes to be used for dimensionality reduction
        (centering and PCA)
    batch_key
        a column name in ``adata.obs`` specifying the batch labels
    copy
        whether to make a co[y of the input data. if `False`, the data object
        will be change inplace.
    log_first
        for some data distributions, perform log(x+1) before total-count
        normalization might give a better result (e.g. clustering results
        may be less affected by the sequencing depths)
    hvg_kwds
        other key-word parameters for ``sc.pp.highly_variable_genes``

    Returns
    -------

    """
    if copy:
        _adata = adata.copy()
        logging.info('A copy of AnnData made!')
    else:
        _adata = adata
        logging.info('No copy was made, the input AnnData will be changed!')
    # 1: normalization
    if log_first:
        normalize_log_then_total(_adata, target_sum=target_sum)
    else:
        normalize_default(_adata, target_sum=target_sum)
    # 2: HVG selection (skipped if `hvgs` is given)
    if hvgs is None:
        sc.pp.highly_variable_genes(
            _adata, batch_key=batch_key, **hvg_kwds)
        indicator = _adata.var['highly_variable']
        # _adata = _adata[:, _adata.var['highly_variable']]
    else:
        indicator = None
    _adata = set_adata_hvgs(_adata, gene_list=hvgs, indicator=indicator, )
    # 3: z-score
    wrapper_scale(_adata, groupby=batch_key)
    return _adata


def label_binarize_each(labels, classes, sparse_out=True):
    lb1hot = label_binarize(labels, classes=classes, sparse_output=sparse_out)
    if len(classes) == 2:
        lb1hot = lb1hot.toarray()
        lb1hot = np.c_[1 - lb1hot, lb1hot]
        if sparse_out:
            lb1hot = sparse.csc_matrix(lb1hot)
    return lb1hot


def group_mean(X, labels,
               binary=False, classes=None, features=None,
               print_groups=True):
    """
    This function may work with more efficiency than `df.groupby().mean()`
    when handling sparse matrix.

    Parameters
    ----------
    X: shape (n_samples, n_features)
    labels: shape (n_samples, )
    classes: optional
        names of groups
    features: optional
        names of features
    print_groups: bool
        whether to inspect the groups
    """
    classes = np.unique(labels, ) if classes is None else classes
    if binary:
        X = (X > 0)  # .astype('float')
        print('Binarized...the results will be the expression proportions.')

    if len(classes) == 1:
        grp_mean = X.mean(axis=0).T
    else:
        lb1hot = label_binarize_each(labels, classes=classes, sparse_out=True)
        print(f'Calculating feature averages for {len(classes)} groups')
        if print_groups:
            print(classes)
        grp_mean = X.T.dot(lb1hot) / lb1hot.sum(axis=0)
    grp_mean = pd.DataFrame(grp_mean, columns=classes, index=features)
    return grp_mean


def group_mean_dense(
        X, labels, binary=False,
        index_name='group',
        classes=None,
):
    classes = np.unique(labels, ) if classes is None else classes
    if binary:
        X = (X > 0)  # .astype('float')
        logging.info('Binarized...the results will be the expression '
                     'proportions.')
    tmp = pd.DataFrame(X)
    tmp[index_name] = list(labels)
    avgs = tmp.groupby(index_name).mean()
    # print(avgs.shape)
    return avgs.T  # each column as a group


def group_median_dense(
        X, labels, binary=False,
        index_name='group',
        classes=None,
):
    classes = np.unique(labels, ) if classes is None else classes
    if binary:
        X = (X > 0)  # .astype('float')
        print('Binarized...the results will be the expression proportions.')
    tmp = pd.DataFrame(X)
    tmp[index_name] = list(labels)
    avgs = tmp.groupby(index_name).median()
    print(avgs.shape)
    return avgs.T  # each column as a group


def group_mean_adata(adata: sc.AnnData,
                     groupby: str,
                     features=None,
                     binary=False,
                     use_raw=False):
    """
    Compute averaged feature-values for each group

    Parameters
    ----------
    adata: AnnData
    groupby: str
        a column name in adata.obs
    features:
        a subset of names in adata.var_names (or adata.raw.var_names)
    binary: bool
        if True, the results will turn to be the non-zeor proportions for
        all (or the given) features
    use_raw: bool
        whether to access adata.raw to compute the averages.

    Returns
    -------
    a pd.DataFrame with features as index and groups as columns
    """
    labels = adata.obs[groupby]
    logging.debug(f'Computing averages grouped by {groupby}')
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
    if sparse.issparse(X):
        return group_mean(X, labels, binary=binary, features=features)
    else:
        return group_mean_dense(X, labels, binary=binary, )


def __test__():
    pass


if __name__ == '__main__':
    import time

    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s %(filename)s-%(lineno)d-%(funcName)s(): '
               '%(levelname)s\n%(message)s'
    )
    t = time.time()
    __test__()
    print('Done running file: {}\nTime: {}'.format(
        os.path.abspath(__file__), time.time() - t,
    ))
