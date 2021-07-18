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
from ._scale import wrapper_scale


def check_dirs(path):
    if os.path.exists(path):
        print('already exists:\n\t%s' % path)
    else:
        os.makedirs(path)
        print('a new directory made:\n\t%s' % path)


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


def set_precomputed_neighbors(
        adata,
        distances,
        connectivities=None,
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
