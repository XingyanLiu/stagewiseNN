# -*- coding: UTF-8 -*-
"""
@CreateDate: 2021/07/25
@Author: Xingyan Liu
@File: test_builder.py
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

ROOT = Path('.')
sys.path.append(str(ROOT))
import swnn

DATADIR = ROOT / 'sample_data'


def get_adata(datadir=DATADIR, ):
    path = datadir / 'subsampled_B-L0-0.2.h5ad'
    adata = sc.read_h5ad(path)
    return adata


def get_high_freq_hvgs(min_freq=3, datadir=DATADIR):
    hvg_freq = pd.read_csv(
        datadir / 'hvg_frequencies.csv', index_col=0, header=None)
    hvg_freq = hvg_freq.iloc[:, 0]
    return hvg_freq[hvg_freq >= min_freq].index.tolist()


def main(resdir: Union[str, Path] = '_temp/builder'):
    resdir = Path(resdir)
    swnn.check_dirs(resdir)

    adata_raw = get_adata()
    hvgs = get_high_freq_hvgs()
    adata = swnn.quick_preprocess_raw(
        adata_raw, hvgs=hvgs, copy=True, batch_key='stage_primer')

    stage_order = ("B", "G3", "G4", "G5", "G6", "N0", "N1", "N3", "L0")
    ks = [20] + [10] * 6 + [5] + [3]
    n_pcs = [30] * 5 + [50] * 4

    builder = swnn.Builder(stage_order=stage_order)
    distmat, connect = builder.build_graph(
            X=adata.X,
            stage_lbs=adata.obs['stage_name'],
            ks=ks, n_pcs=n_pcs, binary_edge=True
        )

    # =========[optional]=========
    swnn.set_precomputed_neighbors(adata, distmat, connect)
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
    # sc.pl.umap(
    #     adata, color='stage_name', palette='plasma_r', save='_stage.pdf')

    obs = adata.obs
    group_lbs = obs['stagewise_cluster'].values
    stage_lbs = obs['stage_name'].values
    KEY_TREE_NODE = 'tree_node'

    builder.build_tree(group_lbs, stage_lbs,)
    obs[KEY_TREE_NODE] = builder.refined_group_lbs

    # =========[optional]=========
    groupby = KEY_TREE_NODE
    props_all = swnn.group_mean_adata(adata, groupby, use_raw=True, binary=True)
    means_all = swnn.group_mean_adata(adata, groupby, use_raw=True, )
    props_all.to_csv(resdir / f'expr_prop_all.csv', index=True, header=True)
    means_all.to_csv(resdir / f'avg_expr_all.csv', index=True, header=True)


def __test__():
    main()


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
