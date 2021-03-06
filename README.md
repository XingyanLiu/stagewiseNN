StagewiseNN
===========

[//]: # (https://zenodo.org/badge/386473402.svg)
[![DOI](https://zenodo.org/badge/386473402.svg)](https://zenodo.org/badge/latestdoi/386473402)

**StagewiseNN** is a computational tool for constructing
developmental tree from multi-staged single-cell RNA-seq data.

It starts from building a single-cell graph by connecting each cell to its
k-nearest neighbors in the parent stage, followed by the voting-based tree-construction
and an adaptive cluster refinement.

> see [StagewiseNN Documentation ](https://xingyanliu.github.io/stagewiseNN/index.html) for detailed guides

[//]: # (![StagewiseNN]&#40;docs/source/_figs/swnn_overview.png&#41;)
<img src="docs/source/_figs/swnn_overview.png" height="310"/>

The single-cell graph can be further visualized using graph embedding methods, e.g. UMAP, SPRING.

We have used it to build the developmental tree from the **snRNA-seq** of amphioxus embryonic cells, 
across nine developmental stages ("B", "G3", "G4", "G5", "G6", "N0", "N1", "N3", "L0"),
where seven major lineages were recognized.

<img src="docs/source/_figs/umap_rna1.png" height="240"/>

StagewiseNN can also be applied on **scATAC-seq** data sampled at multiple timepoints,
once the peak-by-cell matrix is transformed into a gene-by-cell matrix (i.e., the gene activities).

<img src="docs/source/_figs/umap_atac.png" height="236"/>

It is easy to use:

```python
import swnn

# ====== Inputs ======
# data_matrix = ..
# stage_labels = ..
# group_labels = ..
# stage_order = [f'stage_{i}' for i in range(5)]

builder = swnn.Builder(stage_order=stage_order)

# step1: building a (stage-preserved) single-cell graph
distmat, connect = builder.build_graph(
       X=data_matrix, stage_lbs=stage_labels,
   )
# step2: build a developmental tree from the single-cell graph
edgedf, refined_group_lbs = builder.build_tree(group_labels, stage_labels,)

```

> Here is [the detailed tutorial](https://xingyanliu.github.io/stagewiseNN/tutorial/tutorial_builder_based.html)


Installation
------------

Requirements:

- python >= 3.6
- [scanpy](https://pypi.org/project/scanpy/)
- [scikit-learn](https://pypi.org/project/scikit-learn/)


Install stagewiseNN by running (in the command line):

```shell
pip install swnn
```

or install from source code:

```shell
git clone https://github.com/zhanglabtools/stagewiseNN.git
cd stagewiseNN
python setup.py install
```

Contribute
----------

- Issue Tracker: https://github.com/XingyanLiu/stagewiseNN/issues
- Source Code: 
  - https://github.com/zhanglabtools/stagewiseNN
  - https://github.com/XingyanLiu/stagewiseNN (the developmental version)

Support
-------

If you are having issues, please let us know.
We have a mailing list located at: 

* xingyan@amss.ac.cn
* 544568643@qq.com

Citation
--------
If you find StagewiseNN helps, please cite:

> Pengcheng Ma, Xingyan Liu, Zaoxu Xu et al.,
> **Joint profiling of gene expression and chromatin accessibility of amphioxus 
> development at single-cell resolution**,
> _Cell Reports_ (2022), https://doi.org/10.1016/j.celrep.2022.110979

<img src="docs/source/_figs/graphic_summ.png" height="400"/>