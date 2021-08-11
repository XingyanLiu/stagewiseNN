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

See :doc:`apidoc/swnn.builder` for detailed information.

.. autosummary::
   :recursive:
   :toctree: generated/

   Builder


Data Processing
---------------

.. autosummary::
   :toctree: generated/

   quick_preprocess_raw
   normalize_default
   groupwise_hvgs_freq
   set_adata_hvgs
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


