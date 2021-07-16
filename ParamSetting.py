# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 22:36:18 2019

@author: xyliu
"""

###  ParamSetting



Parameters0 = dict(
        embryos_whole = dict(qc = True,  min_genes = 100, max_genes = 3000,
                  rmv_mito=False, #mito_perc=0.02, 
                  counts_level = None,
                  plot_hvgs=True, min_mean = 0.015, min_disp = 0.25, 
                  do_regress=False, batch_key='primer',
                  n_comps=80, 
                  metric = 'cosine',
                  nneigh = 20, n_pcs=60, 
                  min_dist = 0.5,#0.3,
                  de = True, plot_de=True,
                  cluster=True, 
                  clust_reso = 1,
                  save_middle = True, 
                  save_by_default = True
                  ),
                             
        embryos = dict(qc = True,  min_genes = 100, max_genes = 3000,
                  rmv_mito=False, #mito_perc=0.02, 
                  counts_level = None,
                  plot_hvgs=True, min_mean = 0.015, min_disp = 0.25, 
                  do_regress=False, batch_key=None,
                  n_comps=50, 
                  metric = 'cosine',
                  nneigh = 20, n_pcs=30, 
                  min_dist = 0.3,
                  de = True, plot_de=True,
                  cluster=True, 
                  save_middle = True, 
                  save_by_default = True
                  ),
        adults = dict(qc = True,  min_genes = 100, max_genes = 3000,
                  rmv_mito=False, #mito_perc=0.02, 
                  counts_level = None,
                  plot_hvgs=True, min_mean = 0.02, min_disp = 0.25, 
#                  hvg_batch = 'RNAcap',
                  do_regress=False, batch_key='primer',
                  n_comps=80, 
                  metric = 'cosine',
                  nneigh = 20, n_pcs=50, 
                  min_dist = 0.3,
                  de = True, plot_de=True,
                  cluster=True, 
                  save_middle = True, 
                  save_by_default = True
                  ),
        All = dict(qc = True,  min_genes = 350, max_genes = 5000,
                  rmv_mito = True, mito_perc=0.02, 
                  counts_level = None,
                  plot_hvgs=True, min_mean = 0.035, min_disp = 0.25, 
                  do_regress=False, batch_key=None,
                  n_comps=100, 
                  metric = 'cosine',
                  nneigh = 30, n_pcs = 60, 
                  min_dist = 0.5,
                  de = True, plot_de = True,
                  cluster=True, 
                  save_middle = True, 
                  save_by_default = True
                  ),
        NS = dict(qc = True,  min_genes = 150, max_genes = 2500,
                  rmv_mito = True, mito_perc=0.03, 
                  counts_level = None,
                  plot_hvgs=True, 
                  min_mean = 0.012, 
                  min_disp = 0.25, 
                  do_regress=False, batch_key='primer',
                  n_comps=50, 
                  metric = 'cosine',
                  nneigh = 15, n_pcs = 30, 
                  min_dist = 0.25,
                  de = True, plot_de = True,
                  cluster=True, 
                  save_middle = True, 
                  save_by_default = True
                  ),
                  
        Hipp_NS = dict( # for `ph.AutoDownstream()`
                  counts_level = 1000,
                 # plot_hvgs=True, min_mean = 0.012, min_disp = 0.25, 
                  do_regress=False, batch_key='sample_RNAcap',
                  n_comps=80, 
                  metric = 'cosine',
                  nneigh = 30, n_pcs = 60, 
                  min_dist = 0.25,
                  de = True, plot_de = True,
                  cluster=True, 
                  save_middle = True, 
                  save_by_default = True
                  ),
        Brainstem = dict(qc = True,  min_genes = 250, max_genes = 3500,
                  rmv_mito = True, mito_perc=0.05, 
                  counts_level = None,
                  plot_hvgs=True, min_mean = 0.015, min_disp = 0.25, 
                  do_regress=False, batch_key='primer',
                  n_comps=80, 
                  metric = 'cosine',
                  nneigh = 15, n_pcs = 50, 
                  min_dist = 0.5,
                  de = True, plot_de = True,
                  cluster=True, 
                  save_middle = True, 
                  save_by_default = True
                ),
        ref_KO = dict( # for `ph.AutoDownstream()`
                  counts_level = 1000,
                  do_regress=False, batch_key='batch',
                  n_comps=100, 
                  metric = 'cosine',
                  nneigh = 30, n_pcs = 80, 
                  min_dist = 0.5,
                  de = True, plot_de = True,
                  cluster=True, 
                  save_middle = True, 
                  save_by_default = True
                ),
                    )
        
        
MinGenesQC = dict(
        E1 = 350,
        E2 = 350,
        E3 = 300,
        E4 = 300,
        E5 = 250,
        E6 = 400,
        E7 = 350,
        E8 = 250,
        E9 = 150,
        E10 = 350,
        E11 = 150,
        E12 = 150, 
        E14 = 100,
        E15 = 100,
        )

MinGenes0 = dict(
        embryos = [200, 200, 200, 200, 150, 
                 250, 200, 250, 150, 200,
                 150, 150, 150, 100, 100],
        NS = [20, 20, 20, 30, 30, 20],
        
        )        

NPCs0 = dict(
        embryos = [5, 8, 10, 10, 10,
                   15, 15, 20, 25, 25, 
                   25, 25, 30, 30, 30], 
        NS = [20, 20, 20, 30, 30, 20],
        
        )        

NS = [10, 15, 15, 15, 20, 20]

def Parameters():
    return Parameters0.copy()

def Parameters1(key):
    return Parameters0[key].copy()

def npc_dict():

    return NPCs0.copy()

def qc_n_genes():
    return MinGenes0.copy()

'''
    min_genes = 200
    max_genes = 7000
    n_pcs = 50
    min_mean = 0.012
    min_disp = 0.25
    counts_level = None
    if Tp == 'embryos':
        n_pcs = [20, 20, 20, 20, 20,
                 20, 25, 25, 25, 30, 
                 30, 30, 30, 30, 35][i]
        min_genes = [200, 200, 200, 250, 100, 
                     250, 200, 200, 120, 250,
                     120, 120, 150, 100, 75][i]
        min_mean = 0.025 if i < 10 else 0.012
        min_disp = 0.25
    elif Tp == 'BF':
        n_pcs = 50
        min_genes = 350
        max_genes = 7000
    
    elif Tp == 'adults':
        min_genes = 70
        n_pcs = 40
        min_mean = 0.012
        min_disp = 0.25
    #    elif Tp == '

    elif Tp == 'KO':
        pass
'''