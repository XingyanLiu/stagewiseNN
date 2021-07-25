.. module:: swnn
.. automodule:: swnn
   :noindex:

API
===


Import stagewiseNN as:

.. code:: ipython3

   import swnn


.. py:currentmodule:: swnn

Object for Management
---------------------

.. autosummary::
   :toctree: generated/
   :recursive:

   Builder


Data Processing
---------------

.. autosummary::
   :toctree: generated/

   quick_preprocess_raw
   normalize_default
   set_adata_hvgs
   change_names
   group_mean
   group_mean_adata
   wrapper_scale


Make Graph
----------

.. autosummary::
   :toctree: generated/

   stagewise_knn
   adaptive_tree

Others
------

.. autosummary::
   :toctree: generated/

   check_dirs
   set_precomputed_neighbors
   describe_dataframe


