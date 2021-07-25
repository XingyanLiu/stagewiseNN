# -*- coding: UTF-8 -*-
# @CreateDate: 2020/07/18
# @Author: Xingyan Liu
# @File: graph2tree.py
# @Project: stagewiseNN

import os
import sys
from pathlib import Path
from typing import Sequence, List, Mapping, Optional, Union, Callable, Any
import logging
import pandas as pd
import numpy as np
from scipy import sparse
from sklearn.preprocessing import label_binarize
from .utils import label_binarize_each


def adaptive_tree(adj_max: sparse.spmatrix,
                  group_lbs: Sequence,
                  stage_lbs: Optional[Sequence] = None,
                  stage_ord: Optional[Sequence] = None,
                  ignore_pa: Union[List, set] = (),
                  ext_sep: str = '_',
                  ):
    """
    Adaptatively build the developmental tree from the stagewise-KNN graph.

    Parameters
    ----------
    adj_max:
        sparse.csc_matrix, shape = (n_points, n_points)
        adjacent matrix of single points (cells)
    group_lbs:
        np.array, shape = (n_points,)
        group labels specifying each cluster in each stage.
        e.g. 'stage1_1' specifies the cluster 1 in stage1.
    stage_lbs:
        np.array, shape = (n_points,)
        stage labels, better be explicitly assigned. If None, this will be
        extracted from `group_lbs`, and may cause unexpected results.
    stage_ord:
        np.array, shape = (n_stages,)
        order of stages, better provided by user; if None, it will be decided
        automatically.
    ignore_pa: list or set
        parent nodes to be ignored; empty tuple by default.
    ext_sep: str
        parse string for automatically extract the stage-labels from `group_lbs`

    Returns
    -------
    edgedf: pd.DataFrame
        a DataFrame with each row representing an edge, columns are
        ['node', 'parent', 'prop'], where 'prop' is the proportion of nodes
        that have voted for the current parent.
    group_lbs:
        refined group-labels for each sample (e.g. single-cell)

    Examples
    --------
    >>> edgedf, group_lbs = adaptive_tree(adj_max, group_lbs, stage_ord=stages)

    """
    group_lbs = np.asarray(group_lbs).copy()
    if stage_lbs is None:
        stage_lbs = _extract_field(group_lbs, sep=ext_sep, i=0)
    if stage_ord is None:
        print('stage orders are decided automatically:')
        stage_ord = pd.unique(stage_lbs)
        print(stage_ord)

    stg_grp_dict = make_stage_group_dict(group_lbs, stage_lbs=stage_lbs)

    adj_max = sparse.csc_matrix(adj_max)
    edgedfs = []
    for i, stg0 in enumerate(stage_ord[: -1]):

        groups0 = stg_grp_dict[stg0]
        if len(groups0) < 2:
            continue

        stg1 = stage_ord[i + 1]
        inds0 = np.flatnonzero(stage_lbs == stg0)
        inds1 = np.flatnonzero(stage_lbs == stg1)
        adj_sub = adj_max[:, inds1][inds0, :]
        group_lbs0 = group_lbs[inds0]
        group_lbs1 = group_lbs[inds1]
        #        groups1 = stg_grp_dict[stg1]
        print(f'connecting stage {stg0} and {stg1}')

        edgedf_sub, new_group_lbs1, new_groups1 = connect_near_stages(
            adj_sub,
            group_lbs0, group_lbs1,
            groups0=groups0, groups1=stg_grp_dict[stg1],
            # stg0=stg0, stg1=stg1,
            df2edges=True,
            ignore_pa=ignore_pa,
        )
        group_lbs[inds1] = new_group_lbs1
        stg_grp_dict[stg1] = new_groups1
        edgedfs.append(edgedf_sub)
        print()

    edgedf = pd.concat(edgedfs, axis=0, ignore_index=True)

    return edgedf, group_lbs  # , stg_grp_dict


def connect_near_stages(adj_max,
                        group_lbs0, group_lbs1,
                        groups0=None, groups1=None,
                        # stg0=None, stg1=None,
                        ignore_pa: Sequence = (),
                        df2edges: bool = True,
                        sep: str = '_',
                        ):
    """
    Parameters
    ----------
    adj_max:
        csc_sparse matrix, adjacent matrix of samples (cells),
        shape (n_samples, n_samples)
    group_lbs0:
        group labels of the parent stage
    group_lbs1:
        group labels of the descendent stage

    Returns
    -------
    edgedf:
        a DataFrame with each row representing an edge, columns are
        ['node', 'parent', 'prop'], where 'prop' is the proportion of nodes
        that vote for the current parent.
    new_group_lbs1:
        np.array, refined group labels
    new_groups1:
        unique group names from `new_group_lbs1`

    """
    adj_max = sparse.csc_matrix(adj_max)
    groups0 = pd.unique(group_lbs0) if groups0 is None else groups0
    groups1 = pd.unique(group_lbs1) if groups1 is None else groups1
    # each column normalized to unit sum
    #    group_conn = agg_group_edges(adj_max, group_lbs0, group_lbs1,
    #                                 groups0=groups0, groups1=groups1,)
    #    print(group_conn)
    voting_props = agg_group_edge_props(adj_max, group_lbs0, group_lbs1,
                                        groups0=groups0, groups1=groups1,
                                        axis=0)

    winner_pas = voting_props.idxmax(axis=0)  # a series
    single_pas = [x for x in groups0 if x not in winner_pas.values]
    print('parent nodes that had no descendent:', single_pas)
    for p in ignore_pa:
        if p in single_pas:
            print(f'ignore single parent node {p}')
            single_pas.remove(p)

    # modify the groups labels in group1
    is_strlbs = isinstance(group_lbs1[0], str)
    if is_strlbs:
        stg1 = groups1[0].split(sep)[0]
        new_group_lbs1 = _extract_field(group_lbs1, sep=sep, i=-1).astype(int)
    else:  # integer labels
        new_group_lbs1 = group_lbs1.copy()
    max_num = new_group_lbs1.max()

    print('Taking descendant-points from other nodes (groups)')
    #    adj_max = max_connection(adj_max, axis=0)
    for i, pa in enumerate(single_pas):
        parent_ids = group_lbs0 == pa
        sub_adj = adj_max[parent_ids, :]  # sns.heatmap(sub_adj.toarray())
        rows, cols = sub_adj.nonzero()  # `cols` is the indices of child-nodes
        new_name = max_num + 1 + i
        new_group_lbs1[cols] = new_name

    new_groups1 = np.unique(new_group_lbs1)
    if is_strlbs:
        print('pasting stage labels')
        new_groups1 = [f'{stg1}{sep}{x}' for x in new_groups1]
        new_group_lbs1 = np.array(
            list(map(lambda x: f'{stg1}{sep}{x}', new_group_lbs1)))

    # ========= new voting proportions ============
    voting_props = agg_group_edge_props(adj_max, group_lbs0, new_group_lbs1,
                                        groups0=groups0, groups1=new_groups1,
                                        axis=0, verbose=True)

    edgedf = max_connection(voting_props, axis=0, df2edges=df2edges)
    return edgedf, new_group_lbs1, new_groups1


def agg_group_edge_props(adj, group_lbs0, group_lbs1=None,
                         groups0=None, groups1=None,
                         #                    asdf = True,
                         axis=0, verbose=True):
    group_conn = agg_group_edges(adj, group_lbs0, group_lbs1=group_lbs1,
                                 groups0=groups0, groups1=groups1, asdf=True,
                                 verbose=verbose)
    # normalize each column to unit sum
    voting_props = group_conn.apply(lambda x: x / x.sum(), axis=axis)
    return voting_props


def agg_group_edges(adj, group_lbs0, group_lbs1=None,
                    groups0=None, groups1=None, asdf=True, verbose=True):
    """
    Parameters
    ----------
    adj:
        adjacent matrix of shape (N0, N1), if `group_lbs1` is None, then set N0=N1.
    group_lbs0:
        a list or a np.array of shape (N0,)
    group_lbs1:
        a list or a np.array of shape (N1,)

    Returns
    -------
    group_conn: summation of connected edges between given groups
    """
    #    if sparse.issparse(adj):
    adj = sparse.csc_matrix(adj)

    groups0 = pd.unique(group_lbs0) if groups0 is None else groups0
    lb1hot0 = label_binarize_each(group_lbs0, classes=groups0, sparse_out=True)
    if group_lbs1 is None:
        lb1hot1 = lb1hot0
        groups1 = groups0
    else:
        groups1 = pd.unique(group_lbs1) if groups1 is not None else groups1
        lb1hot1 = label_binarize_each(group_lbs1, classes=groups1,
                                      sparse_out=True)
    if verbose:
        print('---> aggregating edges...')
        print('unique labels of rows:', groups0)
        print('unique labels of columns:', groups1)
        print('grouping elements (edges)')
        print('shape of the one-hot-labels:', lb1hot0.shape, lb1hot1.shape)
    group_conn = lb1hot0.T.dot(adj).dot(lb1hot1)
    #    print(group_conn.shape)
    if asdf:
        group_conn = pd.DataFrame(group_conn.toarray(),
                                  index=groups0, columns=groups1)

    return group_conn


def max_connection(adj, axis=0, df2edges=False):
    """
    keep only max element (connection) for each column (axis=0), and remove
    the other elements (connections)
    """
    if isinstance(adj, pd.DataFrame):
        def keep_max(x):
            cut = x.max()
            x[x < cut] = 0
            return x

        adj_max = adj.apply(keep_max, axis=axis)
        if df2edges:
            adj_max = _make_edgedf(adj_max, col_key='node', row_key='parent',
                                   data_key='prop', )
    else:
        adj = sparse.csc_matrix(adj)
        shape = adj.shape
        vmaxs = adj.max(axis=0).A[0]  # .A[0] for csc matrix
        idxmax = adj.argmax(axis=axis).A[0]  # .A[0] for csc matrix
        adj_max = sparse.coo_matrix(
            (vmaxs,
             (idxmax, np.arange(shape[1]))),
            shape=shape)
    return adj_max


def _make_edgedf(df, col_key='node', row_key='parent',
                 data_key='prop', ):
    coo_data = sparse.coo_matrix(df.values)
    edgedf = pd.DataFrame({
        col_key: df.columns.take(coo_data.col),
        row_key: df.index.take(coo_data.row),
        data_key: coo_data.data,
    })

    return edgedf


def make_stage_group_dict(group_lbs, stage_lbs=None):
    if stage_lbs is None:
        stage_lbs = _extract_field(group_lbs, sep='_', i=0)

    stages = pd.unique(stage_lbs)
    group_lbs = np.asarray(group_lbs)
    dct = {}
    for stg in stages:
        dct[stg] = pd.unique(group_lbs[stage_lbs == stg])
    return dct


def find_children(
        nodes: Sequence,
        children_dict: Mapping[Any, Sequence],
        n: int = 100):
    """
    Parameters
    ----------
    nodes:
        better a list of node(s) to be looked up
    children_dict: dict
        parent (key) -> children (value)
    n:
        max number of iterations
    """
    if isinstance(nodes, (int, str)):
        nodes = [nodes]
    has_children = [nd in children_dict.keys() for nd in nodes]
    has_any_children = any(has_children)
    if not has_any_children:
        return []
    if n < 1:
        return []

    children = []
    for i, nd in enumerate(nodes):
        if has_children[i]:
            children += children_dict[nd]

    n -= 1
    children = children + find_children(children, children_dict, n=n)
    #    if None in children:
    #        children.remove(None)
    return children


def make_children_dict(df_tree, column_pa='parent', column_ch='node'):
    """ making a dict for looking up descendants
    """
    children_dict = df_tree.groupby(column_pa)[column_ch].apply(lambda x: list(x))
    return children_dict.to_dict()


def _extract_field(labels, sep='_', i=0):
    return np.array(list(map(lambda x: x.split(sep)[i], labels)))


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
