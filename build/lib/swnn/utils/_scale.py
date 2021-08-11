# -*- coding: UTF-8 -*-
"""
@CreateDate: 2021/07/18
@Author: Xingyan Liu
@File: _scale.py
@Project: stagewiseNN
"""
import os
import sys
from pathlib import Path
from typing import Sequence, Mapping, Optional, Union, Callable
import logging
import pandas as pd
import numpy as np
from scipy import sparse
import scanpy as sc
from sklearn.preprocessing import StandardScaler, label_binarize


def zscore(X, with_mean=True, scale=True, ):
    """ For each column of X, do centering (z-scoring)
    ====
    code borrowed from `scanpy.pp._simple`
    """
    scaler = StandardScaler(with_mean=with_mean, copy=True).partial_fit(X)
    if scale:
        # user R convention (unbiased estimator)
        e_adjust = np.sqrt(X.shape[0] / (X.shape[0] - 1))
        scaler.scale_ *= e_adjust
    else:
        scaler.scale_ = np.array([1] * X.shape[1])
    X_new = scaler.transform(X)
    if isinstance(X, pd.DataFrame):
        X_new = pd.DataFrame(X_new, index=X.index, columns=X.columns)
    return X_new


def group_zscore(X, labels, with_mean=True, scale=True, max_value=None):
    """
    For each column of X, do within-group centering (z-scoring)
    ======
    X: np.array, shape (n_samples, n_features)
        i.e. each row of X is an observation, wile each column is a feature

    with_mean: boolean, True by default
        If True, center the data before scaling, and X shoud be a dense matrix.
    """
    isdf = False
    if isinstance(X, pd.DataFrame):
        isdf = True
        index, columns, X = X.index, X.columns, X.values
    X = X.astype(np.float).copy()
    labels = np.asarray(labels)
    unique_labels = np.unique(labels)
    for lb in unique_labels:
        ind = labels == lb
        if sum(ind) == 1:
            logging.warning(f'ignoring class {lb} with only one sample.')
            continue
        X[ind, :] = zscore(X[ind, :], with_mean=with_mean, scale=scale)

    if max_value is not None:
        X[X > max_value] = max_value
        logging.info('... clipping at max_value', max_value)

    if isdf:
        X = pd.DataFrame(X, index=index, columns=columns)
    return X


def group_zscore_adata(adt, key='counts', groupby='batch', key_new=None,
                       max_value=None,
                       with_mean=True,
                       cover=True, **kwds):
    """
    adt: AnnData
    key: str in ['X_pca', 'count']
        can be a key from adt.obsm, e.g. `key='X_pca'`
        If key == 'counts', then do scaling on `adt.X`
        and cover the old count matrix, ignoring the `cover` parameter
    groupby: str;
        A key from adt.obs, from which the labels are take
    cover: bool
        whether to cover the old X with the scored X
    """
    labels = adt.obs[groupby]
    if key == 'counts':
        logging.info('doing z-score scaling on count matrix, '
                     'transformed into a dense array')
        if sparse.issparse(adt.X):
            X = adt.X.toarray()
        else:
            X = adt.X
        if not cover:
            adt = adt.copy()
        X_scaled = group_zscore(
            X, labels, with_mean=with_mean,
            max_value=max_value, **kwds)
        # logging.debug('X_scaled = %s', X_scaled)
        adt.X = X_scaled
    else:
        if cover:
            key_new = key
        else:
            key_new = key + '_new' if key_new is None else key_new
        adt.obsm[key_new] = group_zscore(
            adt.obsm[key], labels, with_mean=with_mean, max_value=None, **kwds)

    return adt


def wrapper_scale(adata, zero_center=True, max_value=None,
                  groupby=None, copy=False, **kwds):
    """
    Wrapper function for centering and scaling data matrix `X` in sc.AnnData,
    extended for within-batch cprocessing.

    Example
    =======
        wrapper_scale(adata, groupby='batch')
    """
    if groupby is not None:
        logging.info(f'doing within-group scaling, group by [ {groupby} ]')
        return group_zscore_adata(adata, key='counts',
                                  max_value=max_value,
                                  groupby=groupby,
                                  with_mean=zero_center,
                                  cover=not copy,
                                  **kwds)
    else:
        logging.info('using the build-in function `sc.pp.scale(..)`')
        return sc.pp.scale(adata, zero_center=zero_center,
                           max_value=max_value, copy=copy)


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
