# -*- coding: UTF-8 -*-
"""
@CreateDate: 2020/07/18
@Author: Xingyan Liu
@File: test_main.py
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

ROOT = Path('.') # os.path.dirname(__file__)
sys.path.append(str(ROOT))
import swnn
from swnn import describe_dataframe, set_adata_hvgs, change_names

DATADIR = ROOT / 'sample_data'


def get_adata(datadir=DATADIR, ):
    path = datadir / 'merged_B-L0-0.2.h5ad'
    adata = sc.read_h5ad(path)
    return adata


def get_high_freq_hvgs(min_freq=3, datadir=DATADIR):
    hvg_freq = pd.read_csv(
        datadir / 'hvg_frequencies.csv', index_col=0, header=None)
    hvg_freq = hvg_freq.iloc[:, 0]
    return hvg_freq[hvg_freq >= min_freq].index.tolist()


def formulate_adata(adata, save_path=None):
    adata.obs.columns = change_names(
        adata.obs.columns, stage='stage_id', )
    adata.obs['lineage'] = change_names(
        adata.obs['lineage'],
        {'Tail bud stem cells': 'Unassigned'}
    )
    adata.obs['stage_primer'] = adata.obs[['stage_name', 'primer']].apply(
        lambda x: '_'.join(x), axis=1
    )
    if save_path is not None:
        adata.write(save_path)
    return adata


def main(resdir: Union[str, Path] = '_temp',
         log_file=None):
    resdir = Path(resdir)
    swnn.check_dirs(resdir)

    adata0 = get_adata()
    adata0 = formulate_adata(adata0)
    print(adata0, file=log_file)
    # describe_dataframe(adata.obs, file=log_file)
    # TODO: additional markers of prior knowledge
    hvgs = get_high_freq_hvgs()
    print('Total of %d HVGs are used.', len(hvgs), file=log_file)
    adata = swnn.quick_preprocess_raw(
        adata0, hvgs=hvgs, copy=True, batch_key='stage_primer')

    X = adata.X
    stage_lbs = adata.obs['stage_name']
    stage_order = ("B", "G3", "G4", "G5", "G6", "N0", "N1", "N3", "L0")
    ks = [10] * 7 + [5] + [3]
    n_pcs = [30] * 5 + [50] * 4

    distmat, connect = swnn.stagewise_knn(
        X, stage_lbs, stage_order,
        k=ks,
        leaf_size=1, # 1 for brute-force KNN
        pca_base_on='stacked',
        n_pcs=n_pcs,
        binary_edge=False,
        )
    connect_bin = swnn.make_binary(connect)
    swnn.set_precomputed_neighbors(adata, distmat, connect_bin, )

    ######################################
    # sc.tl.umap(adata, min_dist=0.1)
    # sc.settings.figdir = resdir
    # sc.set_figure_params(fontsize=14)
    #
    # lin_colors = pd.read_csv(
    #     'sample_data/lineage_colors.csv', index_col=0).iloc[:, 0]
    # adata.obs['lineage'] = pd.Categorical(
    #     adata.obs['lineage'], categories=lin_colors.index)
    # adata.uns['lineage_colors'] = lin_colors.tolist()
    # sc.pl.umap(adata, color='lineage', ncols=1, save='_lineage.pdf')
    # sc.pl.umap(adata, color='stage_name', palette='plasma_r', save='_stage.pdf')

    from scipy import sparse
    obs = adata.obs
    group_lbs = obs['stage_stg_leiden'].values
    stage_lbs = obs['stage_name'].values
    KEY_TREE_NODE = 'stg_groups_new'

    conn_upper = sparse.triu(connect).tocsc()
    adj_max = swnn.max_connection(conn_upper)
    edgedf, new_group_lbs = swnn.adaptive_tree(
        adj_max, group_lbs, stage_lbs=stage_lbs, stage_ord=stage_order)

    obs[KEY_TREE_NODE] = new_group_lbs
    edgedf.prop.hist() # voting proportions
    logging.info("edgedf = %s", edgedf)

    df_tree = edgedf[['node', 'parent']].copy()
    df_tree['label'] = df_tree['node'].copy()
    df_tree['stage'] = df_tree['node'].apply(lambda x: x.split('_')[0])
    groupby = KEY_TREE_NODE
    props_all = swnn.group_mean_adata(adata, groupby, use_raw=True, binary=True)
    means_all = swnn.group_mean_adata(adata, groupby, use_raw=True, )
    props_all.to_csv(resdir / f'expr_prop_all.csv', index=True, header=True)
    means_all.to_csv(resdir / f'avg_expr_all.csv', index=True, header=True)


def __test__():
    # TODO: additional markers of prior knowledge
    filename_log = None
    filename_log = 'data_description-formed.txt'
    if isinstance(filename_log, str):
        with open(filename_log, 'w') as f:
            main(log_file=f)
    else:
        main()


if __name__ == '__main__':
    import time

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(filename)s-%(lineno)d-%(funcName)s(): '
               '%(levelname)s\n%(message)s'
    )
    t = time.time()
    __test__()
    print('Done running file: {}\nTime: {}'.format(
        os.path.abspath(__file__), time.time() - t,
    ))
