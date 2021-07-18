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

ROOT = os.path.dirname(__file__)
sys.path.append(ROOT)


def get_adata():
    path = f'{ROOT}/sample_data/merged_B-L0-0.2.h5ad'
    adata = sc.read_h5ad(path)
    return adata


def __test__():
    adata = get_adata()
    logging.info(adata)


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
