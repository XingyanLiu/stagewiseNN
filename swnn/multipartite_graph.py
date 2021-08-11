# -*- coding: UTF-8 -*-
# @CreateDate: 2020/07/18
# @Author: Xingyan Liu
# @File: multipartite_graph.py
# @Project: stagewiseNN

import os
import sys
from pathlib import Path
from typing import Sequence, Mapping, Optional, Union, Callable
import logging
import pandas as pd
import numpy as np
from scipy import sparse
from sklearn.preprocessing import normalize
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors, KDTree


def stagewise_knn(X: np.ndarray,
                  stage_lbs: Sequence,
                  stage_order: Sequence,
                  leaf_size: Optional[int] = 5,
                  n_pcs: Union[Sequence[int], int] = 30,
                  k: Union[Sequence[int], int] = 30,
                  pca_base_on: Optional[str] = 'stack',
                  binary_edge: bool = True,
                  norm_dists=False,
                  **kwargs
                  ):
    """
        Build multipartite KNN-graph stage-by-stage.

        Parameters
        ----------
        X: np.ndarray or sparse matrix
            data matrix, of shape (n_samples, n_features)
        stage_lbs: Sequence
            stage labels for each sample (nodes in `build_graph`)
        stage_order: Sequence
            stage order
        binary_edge: bool (default=True)
            whether to use the binarized edges. Set as True may cause some
            information loss but a more robust result.
        k:
            the number of nearest neighbors to be calculated.
        n_pcs:
            The number of principal components after PCA reduction.
            If `pca_base_on` is None, this will be ignored.
        pca_base_on: str {'x1', 'x2', 'stacked', None} (default='stacked')
            if None, perform KNN on the original data space.
        leaf_size: int (default=5)
            Leaf size passed to BallTree or KDTree, for adjusting the
            approximation level. The higher the faster, while of
            less promises to find the exact nearest neighbors.
            Setting as 1 for brute-force (exact) KNN.
        norm_dists: bool
            whether to normalize the distance for each pair of adjacent-stages.

        Returns
        -------
        distmat: sparse.csr_matrix
            the distance matrix, of shape (n_samples, n_samples)
        connect: sparse.csr_matrix
            the connectivities matrix, of shape (n_samples, n_samples)

    """
    # setting parameters
    n_pcs = [n_pcs] * len(stage_order) if isinstance(n_pcs, int) else n_pcs
    k = [k] * len(stage_order) if isinstance(k, int) else k

    stage_lbs = np.asarray(stage_lbs)
    N = X.shape[0]
    Inds = np.arange(N)

    iis = []
    dds = []
    jjs = []
    conns = []
    for i, stg1 in enumerate(stage_order):
        if i == 0:
            # stg0 = stg1
            inds_earlier = inds_later = stage_lbs == stg1
            X_earlier = X[inds_earlier, :]
            X_later = None
            print(f'perform KNN searching: {stg1} in {stg1}')
        else:
            stg0 = stage_order[i - 1]
            inds_earlier = stage_lbs == stg0
            inds_later = stage_lbs == stg1
            X_earlier = X[inds_earlier, :]
            X_later = X[inds_later, :]
            print(f'perform KNN searching: {stg1} in {stg0}')

        if pca_base_on is not None:
            X_earlier, X_later = pca_transform(
                X_earlier, X_later,
                n_pcs=n_pcs[i], base_on=pca_base_on
            )

        dist, inds0 = approx_knn(X_earlier, X_later, k=k[i],
                                 leaf_size=leaf_size, **kwargs)
        inds = np.take(Inds[inds_earlier], inds0)
        conn = _dist_to_connection(dist, norm=norm_dists, binarize=binary_edge)

        iis.append(np.repeat(Inds[inds_later], k[i]))
        jjs.append(inds.flatten())
        dds.append(dist.flatten())
        conns.append(conn.flatten())

    logging.info('Done on KNN searching, constructing the results distances '
                 'and connections')
    iis, jjs, dds, conns = tuple(map(
        lambda x: np.concatenate(x), (iis, jjs, dds, conns)
    ))

    distmat = sparse.coo_matrix((dds, (iis, jjs)), shape=(N, N))
    distmat += distmat.T
    connect = sparse.coo_matrix((conns, (iis, jjs)), shape=(N, N))
    connect += connect.T
    logging.info(f'distance matrix of shape {distmat.shape} and '
                 f'{distmat.nnz} non-zeros')
    logging.info(f'connection matrix of shape {connect.shape} and '
                 f'{connect.nnz} non-zeros')
    return distmat, connect


def _dist_to_connection(
        dist: np.ndarray,
        norm: bool = True,
        binarize: bool = False,
        sigma: float = 1.,
):
    """
    dist: shape=(n_samples, n_neighbors)
    """
    if binarize:
        return np.ones_like(dist, dtype=float)
    else:
        if norm:
            conn = normalize(
                _connect_heat_kernel(_normalize_dists(dist), sigma),
                'l1', axis=1,
            ) * np.sqrt(dist.shape[1])
        else:
            conn = normalize(
                _connect_heat_kernel(dist, sigma), 'l1', axis=1
            ) * np.sqrt(dist.shape[1])
        return conn


def _normalize_dists(d):
    """
    d: 2-D np.array, normalization will perform on each row
    """
    return np.array(list(map(_normalize_dists_single, d)))


def _normalize_dists_single(d):
    """
    d: 1-D np.array
    """
    d1 = d[d.nonzero()]
    vmin = d1.min() * 0.99
    #    vmax = d1.max() * 0.99
    med = np.median(d1)
    d2 = (d - vmin) / (med - vmin)
    d2[d2 < 0] = 0
    return d2


def _connect_heat_kernel(d, sigma):
    return np.exp(-np.square(d / sigma))


def pca_transform(
        X1: np.ndarray,
        X2: Optional[np.ndarray] = None,
        n_pcs: int = 50,
        base_on: str = 'stacked',
        **kwargs
):
    """
    base_on: {'x1', 'x2', 'stacked'}
    """
    pca = PCA(n_components=n_pcs, **kwargs)
    if X2 is None:
        logging.debug(f'base_on=X1, and X2=None')
        return pca.fit_transform(X1), None
    if base_on.lower() == 'x1':
        logging.debug(f'PCA base_on: {base_on}')
        X1 = pca.fit_transform(X1)
        X2 = pca.transform(X2)
    elif base_on.lower() == 'x2':
        logging.debug(f'PCA base_on: {base_on}')
        X2 = pca.fit_transform(X2)
        X1 = pca.transform(X1)
    else:
        logging.debug(f'PCA on the stacked data matrix (base_on={base_on})')
        n1 = X1.shape[0]
        X_pca = pca.fit_transform(np.vstack([X1, X2]))
        X1 = X_pca[: n1, :]
        X2 = X_pca[n1:, :]
    return X1, X2


def approx_knn(
        X_ref: np.ndarray,
        X_que: Optional[np.ndarray] = None,
        metric='cosine',
        k: int = 5,
        precis: float = 0.1,
        leaf_size: Optional[int] = 5,
        leaf_size_max: int = 30,
        algorithm='kd_tree',
        metric_params=None,
        **kwargs
):
    """
    Parameters
    ----------
    algorithm : {'auto', 'ball_tree', 'kd_tree', 'brute'}
        default='auto'
    leaf_size : int, default=30
        Leaf size passed to BallTree or KDTree.
    """
    if leaf_size is None:
        n = X_ref.shape[0]
        leaf_size = min([int(np.sqrt(n * n) * precis), leaf_size_max])
        leaf_size = max([leaf_size, 1])
    elif leaf_size <= 1:
        algorithm = 'brute'

    if metric == 'cosine':  # pretend cosine matric
        X_ref = normalize(X_ref, axis=1)
        if X_que is not None:
            X_que = normalize(X_que, axis=1)
        metric = 'minkowski'
    indexer = NearestNeighbors(
        n_neighbors=k,
        algorithm=algorithm,
        leaf_size=leaf_size,
        metric=metric, metric_params=metric_params,
        **kwargs
    ).fit(X_ref)
    dist, inds = indexer.kneighbors(X_que, return_distance=True)
    # indexer = KDTree(X_ref, leaf_size=leaf_size, metric=metric)
    # dist, inds = indexer.query(X_que, k=k)
    return dist, inds


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
