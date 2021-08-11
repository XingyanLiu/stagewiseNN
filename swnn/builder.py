# -*- coding: UTF-8 -*-
"""
@CreateDate: 2021/07/25
@Author: Xingyan Liu
@File: builder.py
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

from .utils import quick_preprocess_raw, make_binary
from .multipartite_graph import stagewise_knn
from .graph2tree import max_connection, adaptive_tree


class BuilderParams:

    def __init__(self, **kwargs):

        self._dict = {}
        self.update(**kwargs)

    def update(self, **kwargs):
        self._dict.update(**kwargs)
        return self

    @property
    def keys(self):
        return self._dict.keys()

    def __getattr__(self, key):
        return self._dict[key]


class Builder(object):

    def __init__(
            self,
            stage_order: Sequence,
            **build_params
    ):
        """
        Parameters
        ----------
        stage_order: Sequence
            the order of stages
        """
        self.stage_order = stage_order
        self._params = BuilderParams(**build_params)
        self._distmat = None
        self._connect = None
        self._stage_lbs = None
        self._group_lbs = None
        self._edgedf = None
        self._refined_group_lbs = None

    @property
    def stage_lbs(self):
        return self._stage_lbs

    @property
    def group_lbs(self):
        return self._group_lbs

    # @group_lbs.setter
    # def group_lbs(self, group_lbs):
    #     pass

    @property
    def distmat(self):
        return self._distmat

    @property
    def connect(self):
        return self._connect

    @property
    def connect_bin(self):
        """ binarized edges """
        if self._connect is not None:
            return make_binary(self._connect)
        return None

    @property
    def edgedf(self):
        return self._edgedf

    @property
    def refined_group_lbs(self):
        return self._refined_group_lbs

    def build_graph(
            self,
            X, stage_lbs,
            binary_edge: bool = True,
            ks: Union[Sequence[int], int] = 10,
            n_pcs: Union[Sequence[int], int] = 50,
            pca_base_on: Optional[str] = 'stacked',
            leaf_size: int = 5,
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
        binary_edge: bool (default=True)
            whether to use the binarized edges. Set as True may cause some
            information loss but a more robust result.
        ks:
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
        kwargs:
            other parameters for `stagewise_knn`

        Returns
        -------
        distmat: sparse.csr_matrix
            the distance matrix, of shape (n_samples, n_samples)
        connect: sparse.csr_matrix
            the connectivities matrix, of shape (n_samples, n_samples)
        """
        self._stage_lbs = stage_lbs
        distmat, connect = stagewise_knn(
            X, self.stage_lbs,
            stage_order=self.stage_order,
            k=ks,
            leaf_size=leaf_size,  # 1 for brute-force KNN
            pca_base_on=pca_base_on,
            n_pcs=n_pcs,
            binary_edge=False,
            **kwargs
        )
        self._distmat = distmat
        self._connect = connect
        if binary_edge:
            connect = self.connect_bin

        # record parameters
        self._params.update(
            binary_edge=binary_edge,
            ks=ks,
            n_pcs=n_pcs,
            pca_base_on=pca_base_on,
            leaf_size=leaf_size,
        )

        return distmat, connect

    def build_tree(
            self,
            group_lbs: Sequence,
            stage_lbs: Optional[Sequence] = None,
            ignore_pa=(),
            ext_sep: str = '_',
    ):
        """
        Adaptatively build the developmental tree from the stagewise-KNN graph.

        Parameters
        ----------
        group_lbs: Sequence
            group labels for each sample (nodes in `build_graph`)
        stage_lbs: Sequence
            stage labels for each sample (nodes in `build_graph`)
        ignore_pa: list or set
            parent nodes to be ignored; empty tuple by default.
        ext_sep: str
            parse string for automatically extract the stage-labels from
            `group_lbs`

        Returns
        -------
        edgedf: pd.DataFrame
            pd.DataFrame of columns {'node', 'parent', 'prop'},
            and of the same number of rows as number of total stage-clusters.
            the column 'prop' is the proportion of nodes that have votes for
            the current parent.
        refined_group_lbs:
            refined group labels for each sample (e.g. single-cell)
        """
        # connect-matrix NOT calculated by StagewiseNN may cause un-expected
        # result by using `sparse.triu()`.
        # TODO: define `take_cross_stage_edges(spmatrix)`
        conn_upper = sparse.triu(self.connect)
        adj_max = max_connection(conn_upper)

        self._group_lbs = group_lbs
        if self.stage_lbs is None:
            self._stage_lbs = stage_lbs

        edgedf, refined_group_lbs = adaptive_tree(
            adj_max, self.group_lbs,
            stage_lbs=self.stage_lbs,
            stage_ord=self.stage_order,
            ignore_pa=ignore_pa,
            ext_sep=ext_sep,
        )

        self._edgedf = edgedf
        self._refined_group_lbs = refined_group_lbs

        # record parameters
        self._params.update(
            ignore_pa=ignore_pa,
            ext_sep=ext_sep,
        )

        return edgedf, refined_group_lbs


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
