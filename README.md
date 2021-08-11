stagewiseNN
===========

To be completed!

**stagewiseNN** is a computational tool for constructing
developmental tree from Multi-staged single-cell RNA-seq data.

(see [documentation](https://xingyanliu.github.io/stagewiseNN/index.html) for detailed guides)

It is easy to use:

```python
import swnn

# ====== Inputs ======
# data_matrix = ..
# stage_labels = ..
# group_labels = ..
# stage_order = [f'stage_{i}' for i in range(5)]

builder = swnn.Builder(stage_order=stage_order)
# step 1:
# building (stage-wise) single-cell graph
distmat, connect = builder.build_graph(
        X=data_matrix, stage_lbs=stage_labels,
    )
# step 2:
# build developmental tree from single-cell graph
builder.build_tree(group_labels, stage_labels,)
```


Installation
------------

Install stagewiseNN by running (in the command line):

```shell
pip install swnn
```

or install from source code:

```shell
git clone https://github.com/XingyanLiu/stagewiseNN.git
cd stagewiseNN
python setup.py install
```

Contribute
----------

- Issue Tracker: https://github.com/XingyanLiu/stagewiseNN/issues
- Source Code: https://github.com/XingyanLiu/stagewiseNN

Support
-------

If you are having issues, please let us know.
We have a mailing list located at: 

* xingyan@amss.ac.cn
* 544568643@qq.com
