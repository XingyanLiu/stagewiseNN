.. stagewiseNN documentation master file, created by
   sphinx-quickstart on Mon Jul 19 12:54:30 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

stagewiseNN - Building developmental tree from scRNA-seq
========================================================

**stagewiseNN** is a computational tool for constructing
developmental tree from Multi-staged single-cell RNA-seq data.

Install from source code:

.. code:: shell

   git clone https://github.com/XingyanLiu/stagewiseNN.git
   cd stagewiseNN
   python setup.py install


It is easy to use:

.. code:: python3

   import swnn

   # ====== Inputs ======
   # data_matrix = ..
   # stage_labels = ..
   # group_labels = ..
   # stage_order = [f'stage_{i}' for i in range(5)]

   builder = swnn.Builder(stage_order=stage_order)
   # step1:
   # building (stage-wise) single-cell graph
   distmat, connect = builder.build_graph(
           X=data_matrix, stage_lbs=stage_labels,
       )
   # step2:
   # build developmental tree from single-cell graph
   builder.build_tree(group_labels, stage_labels,)


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


.. toctree::
   :caption: Contents
   :maxdepth: 1

   installation
   tutorials
   api
   citation


Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
