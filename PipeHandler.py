# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 18:12:37 2019

@author: xyliu
"""


import os
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import sparse
from scipy import io
from scipy.stats import spearmanr
from scipy.spatial.distance import pdist
import scanpy as sc
import matplotlib.pyplot as plt

import funx as fx

#from MatHandler import MatHandler

sc.settings.verbosity = 3



def dist_Spearman(u, v):
    rho, _ = spearmanr(u, v)
    return 1 - rho



def discretize(arr, n_bins = 4, categ = False):
    ''' discretize a vector's continuous elements into several bins
    arr: 1-d array
    n_bins: number of bins
    return discretized 1-d array
    '''
    arr_new = np.array(arr).copy()
    m = len(arr)
    ids = np.argsort(arr)#.tolist()
    s = 0
    for i in range(1, n_bins):
        t = int(i * m / n_bins)
        arr_new[ids[s: t]] = i - 1
        s = t
    arr_new[ids[s: ]] = n_bins - 1
    return arr_new
        
def perc_more_than(X, n = 1):
    
    if isinstance(X, list):
        X = np.array(X)
    nnz = len(np.nonzero(X)[0])
    
    return (X > n).sum() / nnz   


class PipeHandler(object):
    
    '''A data procession handler, processing data while record the parameters.
    Based on `scanpy.Anndata` class
    
    Attributes
    ----------
    params:  
        dict, whose keys are formed like 'OPERATION__xxx'
        For example, 'QC__min_genes', 'HVGs__min_mean'
    
    raw_states:
        dict, with keys including:
            n_genes_mean_raw 
            n_genes_median_raw 
            n_counts_mean_raw 
            n_counts_median_raw
            n_cells_raw 
            n_genes_raw
        
    
    Methods
    -------
    summary(self): 
        print summary

            
    summary_d(self): 
        a dictionary will be returned, whose keys includes:
            min_genes
            n_genes_mean_after
            n_genes_median_after
            min_counts
            n_counts_mean_after
            n_counts_median_after
            n_cells_after 
    
    '''
    
#    PARAMS = dict(
#            DE__groupby = 'louvain'
#            )
    
    def __init__(self, adata, name = None, resdir = None, 
                 summ=True, copy = True):
        
        self.adata = adata.copy() if copy else adata
        self.adata.var_names_make_unique()
        self.adata.obs_names_make_unique()
        self.name = 'temp' if name is None else name
        self.resdir = Path('temp_results' if resdir is None else resdir)
        self.params = dict()
        self.raw_states = dict()
        
        fx.check_dirs(self.resdir)
#        if not os.path.exists(self.resdir):
#            os.makedirs(self.resdir)
        self.set_fig_dir()
        
        # to get the atrributes `n_genes` and `n_counts`
        if 'n_genes' not in self.adata.obs_keys():
            sc.pp.filter_cells(self.adata, min_genes = 0, inplace = True)
        if 'n_counts' not in self.adata.obs_keys():
            sc.pp.filter_cells(self.adata, min_counts = 0, inplace = True)
        
        self.raw_states['n_genes_mean_raw'] = self.adata.obs['n_genes'].mean()
        self.raw_states['n_genes_median_raw'] = self.adata.obs['n_genes'].median()
        self.raw_states['n_counts_mean_raw'] = self.adata.obs['n_counts'].mean()
        self.raw_states['n_counts_median_raw'] = self.adata.obs['n_counts'].median()
        self.raw_states['n_cells_raw'] = self.adata.shape[0]
        self.raw_states['n_genes_raw'] = self.adata.shape[1]
        if summ:
            self.summary()
        self.discretize_n(n_bins = 3)

        
    def summary_d(self):
        
        summ = {
                'min_genes': self.adata.obs['n_genes'].min(), 
                'max_genes': self.adata.obs['n_genes'].max(), 
                'n_genes_mean_after': self.adata.obs['n_genes'].mean(), 
                'n_genes_median_after': self.adata.obs['n_genes'].median(), 
                'min_counts': self.adata.obs['n_counts'].min(), 
                'max_counts': self.adata.obs['n_counts'].max(), 
                'n_counts_mean_after': self.adata.obs['n_counts'].mean(), 
                'n_counts_median_after': self.adata.obs['n_counts'].median(), 
                'n_cells_after': self.adata.shape[0],
                'n_genes_after': self.adata.shape[1]
                }
        if 'highly_variable' in self.adata.var_keys():
            summ['n_hvgs'] = np.sum(self.adata.var['highly_variable'])
        return summ

    def df_param(self):
        return pd.DataFrame([self.params], columns = self.params.keys())
    
    def df_summ(self):
        summ = self.summary_d()
        return pd.DataFrame([summ], columns = summ.keys())
    
    def set_res_path(self, resdir):
        self.resdir = resdir
        fx.check_dirs(self.resdir)
        self.set_fig_dir()

    def summary(self):
        print('--- summary of [ %s ] data ---'%self.name)
        print('result directory: %s'%self.resdir)
        print(self.adata)
        print('percentage of [ more-than-one elements ]: %.3f'% (perc_more_than(self.adata.X)))
        fx.printmm(self.adata, 'n_genes')
        fx.printmm(self.adata, 'n_counts')
        
        
    def plotsumm(self, key = 'n_genes', save = None):
        #  plot summary
        
        fx.plot_combined_log10(self.adata, ft_name=key, save = save)
#        plot_combined_log10(self.adata, 'n_counts')
        
    def set_fig_dir(self, dirname = None):
        self.figdir = os.path.join(self.resdir, 'figs') if dirname == None else dirname
        sc.settings.figdir = self.figdir
        fx.check_dirs(self.figdir)
#        if not os.path.exists(self.figdir):
#            os.makedirs(self.figdir)
    
    def RenameGenes(self, species='amph', mapping=None):
        
#        head = self.adata.var.index[0][:7]
#        species= {'BFwtdbg': 'amph', 'ENSMUSG': 'mouse'}[head]
        head = {'amph': 'BFwtdbg', 'mouse': 'ENSMUSG'}[species]
        def check_and_rename(names):
            if names.str.startswith(head).any():
                print('Rename genes by `id_name` format, stripping `%s`'%head)
                return fx.RenameGenes(names, species, mapping)
            else:
                print('Already renamed, doing nothing')
        
        self.adata.var.index = check_and_rename(self.adata.var.index)
        # check raw data
        if hasattr(self.adata, 'raw') and self.adata.raw is not None:
            self.adata.raw.var.index = check_and_rename(self.adata.raw.var.index)
        

    
    #--------------------------------------------------------------------------
    def discretize_n(self, n_bins = 3):
        self.adata.obs['n_genes_disc'] = discretize(self.adata.obs['n_genes'], n_bins = n_bins)
        self.adata.obs['n_counts_disc'] = discretize(self.adata.obs['n_counts'], n_bins = n_bins)

    
    def QC(self, filter_by = 'n_genes',
           min_genes = None, max_genes = None, 
           plot = True, save_plot = True,
           keep_results = None):
    
        min_cells = 3
        
        print('doing [ Quality Control ]...')
        if plot: 
            plt_file = os.path.join(self.figdir, self.name + '_beforeQC.png') if save_plot else None
            self.plotsumm(save = plt_file)
        _adata = self.adata.copy() if not keep_results else self.adata
        
#        sc.pl.scatter(_adata, x='n_counts', y='n_genes', color = 'n_genes')

        # setting thresholds
        min_genes_ref = _adata.obs['n_genes'].quantile(q = 0.33)
        max_genes_ref = _adata.obs['n_genes'].quantile(q = 0.999)
        min_genes = int(input('cell filter --- min n_genes (ref = %d) = '%min_genes_ref)) if min_genes is None else min_genes
        max_genes = int(input('cell filter --- max n_genes (ref = %d) = '%max_genes_ref)) if max_genes is None else max_genes
        
        # do filtering 
        sc.pp.filter_cells(_adata, min_genes = min_genes)
        sc.pp.filter_genes(_adata, min_cells = min_cells)
        _adata = _adata[_adata.obs['n_genes'] < max_genes, :]
        
        # simply calculate the 'n_genes', 'n_counts'
        sc.pp.filter_cells(_adata, min_counts = 0)
        sc.pp.filter_cells(_adata, min_genes = 0)
        
        # show the results
        if plot: 
            _adata.obs.hist(column = ['n_genes', 'n_counts'], bins = 50, 
                            alpha = 0.5, figsize = (12, 4))
            sc.pl.violin(_adata, ['n_genes', 'n_counts'], 
                         log=True)

            
        plt.show()
        
        print('filter thresholds: genes = [%d, %.f], min_cells = %d'%(min_genes, max_genes, min_cells))
        print('shape after filtering: ', _adata.shape)

        fx.printmm(_adata, 'n_genes')
        fx.printmm(_adata, 'n_counts')
        
        if keep_results is None:
            keep_results = input('Keeping the results, cover the old one?[y/n]').lower() == 'y'
        
        if keep_results: 
            self.adata = _adata
            self.discretize_n(n_bins = 3)
            self.params['QC__min_genes'] = min_genes
            self.params['QC__max_genes'] = max_genes
    
    
    def plot_vln(self, keys=['n_genes', 'n_counts'], 
                 groupby=None, log=True, save=None, **kwds):
        if save and not isinstance(save, str):
            save = f'_{keys}_{groupby}.pdf'
        sc.pl.violin(self.adata, keys, groupby=groupby, 
                     log=log, save=save, **kwds)

        
    def ComputeMito(self, tag='mt-'):
        
        is_mito = self.adata.var_names.str.lower().str.startswith(tag)
        self.adata.var['mito_genes'] = is_mito
        self.adata.obs['percent_mito'] = self.adata[:, is_mito].X.sum(axis=1).A1 / self.adata.obs['n_counts']
        
        
            
    def RmvMito(self, do_filter=None, mito_perc=None, plot=False, keep_results=None):
        # Mito genes
        print('Checking [ Mito-genes ]...')
        if 'percent_mito' not in self.adata.obs_names:
            self.ComputeMito()
        _adata = self.adata.copy() if not keep_results else self.adata
        
        
        sc.pl.violin(_adata, ['n_genes', 'n_counts', 'percent_mito'],
                     jitter=0.4, multi_panel=True, save='_vqc_%s.png'%self.name)
        
        if do_filter is None:
            do_filter = input('filter cells with high `percent_mito`? [y/n]').lower() == 'y'
        if do_filter :
            max_perc_mito = float(input('max percentage of mito genes (ref=0.05) = ')) if mito_perc is None else mito_perc
            _adata = _adata[_adata.obs['percent_mito'] < max_perc_mito, :]
        if plot:
            sc.pl.scatter(_adata, x='n_counts', y='percent_mito')
            sc.pl.scatter(_adata, x='n_counts', y='n_genes')
        # simply calculate the 'n_genes', 'n_counts'
        sc.pp.filter_cells(_adata, min_counts = 0)
        sc.pp.filter_cells(_adata, min_genes = 0)
        
        if keep_results is None:
            keep_results = input('Keeping the results, cover the old one?[y/n]').lower() == 'y'        
        if keep_results: 
            self.adata = _adata
            self.discretize_n(n_bins = 3)
            self.params['Mito__removed'] = True
#        else:
#            return _adata
    
    #--------------------------------------------------------------------------
    def NormLog(self, counts_level = None, frac = 0.05):
        '''
        this will make a copy of the raw count data
        '''
        print('doing [ Normalization and Log-transformation ]...')
        
#        counts_level = self.adata.obs['n_counts'].median() * 2 # multiply 2 to rectify expression levels
        sc.pp.normalize_total(self.adata, target_sum = counts_level, max_fraction=frac)
        sc.pp.log1p(self.adata)
        print('normalization and log-transformation completed !')
        self.adata.raw = self.adata # perserve the row data
        self.params['Norm__counts_level'] = counts_level
        self.params['Norm__frac'] = frac
        
    def NormLog_rev(self, counts_level = None, frac = 0.05):
        '''
        this will make a copy of the raw count data
        '''
        print('doing [Log-transformation and Normalization (reversed)]...')
        
#        counts_level = self.adata.obs['n_counts'].median() * 2 # multiply 2 to rectify expression levels
        sc.pp.log1p(self.adata,)
        sc.pp.normalize_total(self.adata, target_sum = counts_level, max_fraction=frac)
        print('normalization and log-transformation completed !')
        self.adata.raw = self.adata # perserve the row data
        self.params['NormRev__counts_level'] = counts_level
        self.params['NormRev__frac'] = frac

    def BC(self, key, covariates = None):
        '''Work on Dense matric! And might produce negative and Inf elements!
        '''
        print('doing [ Batch Correction ] with key `%s`...'%key)
        print('[Warning] Work on Dense matric! And might produce negative and Inf elements!')
        sc.pp.combat(self.adata, key=key, covariates=covariates)
        self.params['BC__key'] = key
        self.params['BC__covar'] = covariates

            
    #--------------------------------------------------------------------------
    def HVGs(self, min_mean = None, min_disp = None, n_bins=20, 
             max_mean = 4,
             batch_key=None,
             plot = True, redo=None, keep_results = None, **kwds):

        print('doing [ HVGs selection ]...')
        
        _adata = self.adata.copy() if not keep_results else self.adata
            
        min_mean = 0.0125 if min_mean is None else min_mean
        min_disp = 0.25 if min_disp is None else min_disp
        #    sc.pp.highly_variable_genes(_adata, n_top_genes = 1000, n_bins = 30, flavor = 'cell_ranger')
        sc.pp.highly_variable_genes(_adata, # flavor = 'seurat',  # by default
                                    min_mean=min_mean, max_mean = max_mean,
                                    min_disp=min_disp, 
                                    n_bins=n_bins,
                                    batch_key=batch_key,
                                    **kwds) 
        _n_hvgs = _adata.var['highly_variable'].sum()
        if plot:
            sc.pl.highly_variable_genes(_adata, save = '_%s_%d.png'% (self.name, _n_hvgs))
        print('Thresholds for gene selection:\n \
              \t`min_mean`: [ {} ]\t`min_disp`: [ {} ]'.format(min_mean, min_disp))
        print('Number of highly vatiable genes: [ %d ]'%_n_hvgs)
        
        if redo is None:
            redo = input('Reset the thresholds?[y/n]').lower() == 'y' 
        if redo:
            min_mean = float(input('setting `min_mean` = '))
            min_disp = float(input('setting `min_disp` = '))
            _adata = self.HVGs(min_mean, min_disp, plot = plot, 
                               keep_results = False, **kwds)
        
        if keep_results is None:
            keep_results = input('Keeping the results, cover the old one?[y/n]').lower() == 'y'        
        if keep_results: 
            self.adata = _adata[:, _adata.var['highly_variable']]
            self.adata.obs['hvg_counts'] = self.adata.X.sum(axis=1)
            
            self.params['HVGs__min_mean'] = min_mean
            self.params['HVGs__min_disp'] = min_disp
        else:
            return _adata
        

        
        
    def SetHVGs(self, gene_list, slim = True):
        '''
        this step does not mean that the given set of genes are real HVGs
        '''
        print('Setting the given set of %d genes as highly variable'%len(gene_list))
        indicator = [g in gene_list for g in self.adata.var_names]
        self.adata.var['highly_variable'] = indicator
        if slim:
#            if not adata.raw:
            print('slimming adata to contain only HVGs')
            self.adata = self.adata[:, indicator]
    
    def hvgs(self, ):
        genes = self.adata.var_names[self.adata.var['highly_variable']]
        return genes.to_series(index = None, name=self.name)
        
    def save_hvgs(self, name = None):
        name = self.name if name is None else name
        hvg_file = '{}/{}_hvgs.csv'.format(self.resdir, name)
        self.hvgs().to_csv(hvg_file, index = False, header = False)
        
        
    def ZeroFilter(self):
        
        print('doing [ Full-Zero Filtering ]...')
        sc.pp.filter_cells(self.adata, min_genes = 0, inplace = True)
        sc.pp.filter_cells(self.adata, min_counts = 0, inplace = True)
        
        self.summary()
        self.discretize_n(n_bins = 3)    
    
    #--------------------------------------------------------------------------
    
    def RegressOut(self, factors=['n_counts'], keep_results = None):
        '''
        ['percent_mito']
        Problems to be handle!!!
        '''
        _adata = self.adata.copy() if not keep_results else self.adata
        
        if factors is None:
            factors_ = ['n_counts', 'n_genes', 'percent_mito']
            print('select the factors that you want to regress out (e.g. "1 2"):')
            ids = input('\t{}[0]; {}[1]; {}[2]'.format(*factors_)).split()
            factors = [factors_[int(i)] for i in ids]
            
        print('doing [ Regressing out unwanted factors ({}) ]...'.format(factors))
        sc.pp.regress_out(_adata, factors, )#n_jobs=8)
        
        
        if keep_results is None:
            keep_results = input('Keeping the results, cover the old one?[y/n]').lower() == 'y'
        
        if keep_results: 
            self.adata = _adata
            self.params['Regressed__factors'] = factors
            
#    def CheckMito(self, ):
#        if _adata.var['mito_genes'].sum() >= 1:
#            mito_regress = input('do you want to regress out the `percent_mito`? [y / n]')
#            if mito_regress == 'y':
#                print('doing regression')
#                sc.pp.regress_out(_adata, ['percent_mito'])
    
            
    def scale(self, do_regress=False, max_value=10, groupby=None,
              zero_center=True, copy=False,  **kwds):
        '''
        `do_regress` can be
            {bool: default `False`
             string or list: indicating the factors to be regressed out
             None: ask wether doing regress
             }
        if groupby is not None: ignore `do_regress`
        '''
        
#        if groupby is None:
        if do_regress is None:
            do_regress = input('regress out unwanted factors before scaling?[y/n]').lower() == 'y'
        if do_regress:
            factors=None
            if isinstance(do_regress, (str, list)):
                factors = do_regress
            self.RegressOut(factors=factors, keep_results=True)            
    
        print('doing [ Scaling ]...')
        
        if groupby is not None:
            from Build import GroupZscoreAdt
            print('doing within-group scaling, group by [ %s ]'%groupby)
            GroupZscoreAdt(self.adata, key = 'counts', 
                           max_value=max_value,
                           groupby = groupby, 
                           with_mean = zero_center,
                           cover = not copy, **kwds) 
#            print('TEST\n', self.adata)
        else:
            print('using the build-in function `sc.pp.scale(..)`')
            sc.pp.scale(self.adata, zero_center=zero_center, 
                        max_value=max_value, copy=copy)
        
        self.params['scale__groupby'] = groupby
        self.params['scale__regress'] = do_regress
    
    def GscalePCA(self, groupby, max_value=10, cover=True):
        '''group z-score on PC space
        '''
        from Build import GroupZscoreAdt
        self.adata = GroupZscoreAdt(self.adata,
                                    key = 'X_pca', max_value=max_value,
                                    groupby = groupby,
                                    cover = cover)
        self.params['GscalePCA__groupby'] = groupby
        
        
    #--------------------------------------------------------------------------
    def PCA(self, n_comps = 50, plot = True, 
            save_plot = True, **kwds):
        print('doing [ PCA ]...')
        
#        from Build import myPCA
#        myPCA(self.adata, n_comps=n_comps, copy=False, **kwds)
        sc.tl.pca(self.adata, n_comps = n_comps, svd_solver='arpack', **kwds) # n_comps = 50 by default
        if plot:
            plt_file = os.path.join(self.figdir, self.name + '_pca.png') if save_plot else None
            fx.plot_pc_vr(self.adata, n = n_comps, save = plt_file)
            
        self.params['PCA__n_comps'] = n_comps
        
        
    def CrossDimReduct(self, groupby, method = 'plssvd', 
                       n_comps=50, key_add = 'X_cca', details = False):
        '''
        cross decomposition betweeen 2 `batches` using CCA-like methods
        method can be one of the following:
            'plssvd', 'cca', 'plsc'
        '''
        import Build as B
        B.ccaAdt(self.adata, groupby=groupby, method = method, 
                 n_comps=n_comps, key_add = key_add, copy=False)
        
        self.params['CCA__n_comps'] = n_comps
        self.params['CCA__method'] = method
        
    def PPCA(self, based_on, key_add = 'X_pca', n_comps=50, plot=True,
           copy=False, **kwds):
        
        print('doing [ Partial PCA ]...')
        from Build import partialPCAdt
        
        partialPCAdt(self.adata, based_on, key_add, n_comps, copy, **kwds)
        if plot:
            plt_file = os.path.join(self.figdir, self.name + '_pca.png') #if save_plot else None
            fx.plot_pc_vr(self.adata, n = n_comps, save = plt_file)
        
        self.params['PCA__based_on'] = based_on
        self.params['PCA__n_comps'] = n_comps

    
    def EmbedNorm(self, key='X_pca', key_new=None, norm='l2', axis=1, **kwds):
        from sklearn.preprocessing import normalize
        X = self.adata.obsm[key]
        Xn = normalize(X, norm=norm, axis=axis, **kwds)
        key_new = f'{key}_n' if key_new is None else key_new
        self.adata.obsm[key_new] = Xn
        
        
    
    
    
    def Neighbors(self, nneigh = None, do_pca = True, n_pcs = None, 
                  metric = 'cosine', use_rep = 'X_pca', **kwds):
        ''' 
        `metric` can be 'euclidean' 'cosine' 'correlation'... 
        or a Numba!!! function computing distances
        
        `use_rep`: str, any key in `self.adata.obsm`
        '''
        
        print('doing [ K Nearest Neighbor searching ]...')
        print('using `%s` distance metric'% metric)
        if nneigh is None:
            nneigh_ref = 20 if self.adata.shape[0] <= 2000 else 30
            nneigh = int(input('choose the number of nearest neighboes (ref=%d): '% nneigh_ref))
            
        self.params['Neighbor__metric'] = metric
        if metric == 'rankcorr':
            metric = dist_Spearman
        if do_pca and use_rep == 'X_pca':
            if 'X_pca' not in self.adata.obsm.keys():
                self.PCA(n_comps = n_pcs, plot = True)
            if n_pcs is None:
                if 'Neighbor__n_pcs' not in self.params.keys():
                    n_pcs = int(input('enter the number of PCs to use: '))
                else: 
                    n_pcs = self.params['Neighbor__n_pcs'] 
            print('Using the top-[ %s ] components from [ %s ]'%(n_pcs, use_rep))
            sc.pp.neighbors(self.adata, n_neighbors=nneigh, use_rep = use_rep,
                            n_pcs=n_pcs, metric=metric, **kwds)
        else: 
            sc.pp.neighbors(self.adata, n_neighbors=nneigh, use_rep = use_rep,
                            metric=metric, **kwds)

        self.params['Neighbor__use_rep'] = use_rep
        self.params['Neighbor__n_pcs'] = n_pcs
        self.params['Neighbor__nneigh'] = nneigh
        
    
    
    def Umap(self, nneigh = None, plot = True, plot_color=None, **kwds):
        print('doing [ UMAP computation ]...')
        if 'neighbors' not in self.adata.uns.keys():
            self.Neighbors(nneigh = nneigh, do_pca = True)
            
        sc.tl.umap(self.adata, **kwds)
        if plot:
            plot_color = 'n_genes_disc' if plot_color is None else plot_color
            self.vis(color = plot_color, save='_%s_%s.png'% (plot_color, self.name))

    def TSNE(self, plot = True, plot_color=None, use_rep='X_pca', **kwds):
        print('doing [ tSNE computation ]...')
            
        sc.tl.tsne(self.adata, use_rep=use_rep, **kwds)
        if plot:
            plot_color = 'n_genes_disc' if plot_color is None else plot_color
            self.vis(color = plot_color, key = 'tsne', save='_%s_%s.png'% (plot_color, self.name))


    def Diffmap(self, **kwds):
        print('doing [ Diffsion Map computation ]...')
        sc.tl.diffmap(self.adata, **kwds)

    def Graph(self, **kwds):
        print('doing [ Force-directed Graph computation ]...')
        sc.tl.draw_graph(self.adata, **kwds)
    
    def Clustering(self, method = 'leiden',
                   res = None, do_pca = True, n_pcs = None,
                    redo=None, pl_umap = True, #startover = False,
                    save_plot = True, keep_results = None,
                    key_added = None,
                    **kw_neigh):
        
        # checking dependencies
        if 'X_pca' not in self.adata.obsm.keys() and do_pca:
            self.PCA(n_comps = 50, scale = True, plot = True)
        
        if 'neighbors' not in self.adata.uns.keys():
            self.Neighbors(do_pca = do_pca, n_pcs = n_pcs, **kw_neigh)
        
        if ('X_umap' not in self.adata.obsm.keys() and pl_umap):
            self.Umap(plot = False, min_dist=0.05)        
        # clustering
        print('doing [ Clustering ]...')
        _adata = self.adata.copy() if not keep_results else self.adata
        
        res = 1 if res is None else res
        print('setting clustering resolution `res` = %.3g'%res)
        key_added = method if key_added is None else key_added
        if method == 'louvain':
            sc.tl.louvain(_adata, resolution = res) # resolution = 1 by default
        else:
            
            sc.tl.leiden(_adata, resolution = res, key_added=key_added) 
        
        if pl_umap:
            n_pcs = self.adata.uns['neighbors']['params'].get('n_pcs', '')
            metric = self.adata.uns['neighbors']['params'].get('metric', '')
            figname = '{}_pc{}_{}_res{}.png'.format(self.name, 
                       n_pcs, metric, res)
#            plt_file = os.path.join(self.figdir, figname) if save_plot else None
            sc.pl.umap(_adata, color = [key_added, 'n_genes_disc'], 
                       legend_loc='on data', save = figname)
        
        if redo is None:
            redo = input('re-do clustering?[y/n]').lower() == 'y' 
        if redo:
            res = float(input('re-setting clustering `res` = '))
            _adata = self.Clustering(res = res, keep_results = False)

        # save
        if keep_results is None:
            keep_results = input('Keeping the results, cover the old one?[y/n]').lower() == 'y'        
        if keep_results: 
            self.adata = _adata
            self.params['Cluster__method'] = method
            self.params['Cluster__%s__res'%method] = res
        else:
            return _adata
        
        
    def RemoveCluster(self, ):
        pass
    
    def cluster_names(self, key = None):
        if key is None:
            key = 'leiden' if 'leiden' in self.adata.obs.keys() else 'louvain'
        labels = self.adata.obs[key]
        if hasattr(labels, 'cat'):
            return labels.cat.categories.copy()
        return labels.unique()
    
    def n_cluster(self, key = None):
        if key is None:
            key = 'leiden' if 'leiden' in self.adata.obs.keys() else 'louvain'
        return len(self.cluster_names(key = key))
    
    #--------------------------------------------------------------------------
    def DE(self, groupby = None, method = 't-test_overestim_var', #'wilcoxon', #
           use_raw = True, save = None, plot = True, save_plot = True, **kwds):
        '''
        method : `{'logreg', 't-test', 'wilcoxon', 't-test_overestim_var'}`, optional (default: 't-test_overestim_var')
        If 't-test', uses t-test, if 'wilcoxon', uses Wilcoxon-Rank-Sum. If
        't-test_overestim_var', overestimates variance of each group. If
        'logreg' uses logistic regression, see [Ntranos18]_, 
        '''
        print('doing [ Marker Testing ]...')
        if groupby is None:
            groupby = 'leiden' if 'leiden' in self.adata.obs.keys() else 'louvain'
        print('using testing method: [ %s ]' % method)
        sc.tl.rank_genes_groups(self.adata, groupby = groupby, method=method, 
                                use_raw = use_raw, **kwds)
        
        if plot is None:
            plot = (input('plot the markers and their scores?[y/n]').lower() == 'y') 
        if plot:
            
            topn = 5 if self.n_cluster(key = groupby) < 10 else 3
#            sc.pl.rank_genes_groups(self.adata, n_genes=15, sharey=False)
#            plt_file = os.path.join(self.figdir, self.name + '_dot_DE.png') if save_plot else None
            self.dot_de(topn=topn, groupby=groupby, save = '_%s_dot_DE.png' % self.name)
            
        self.params['DE__groupby'] = groupby
        self.params['DE__method'] = method
        self.params['DE__use_raw'] = use_raw
        
        if save is None:
            save = input('save markers to a csv file?[y/n]').lower() == 'y'
        if save:
            self.save_markers()
            self.save_markers(True)

    def markers(self):
        # self.adata.uns['rank_genes_groups'] is a dict
        return pd.DataFrame(self.adata.uns['rank_genes_groups']['names'])
    
    
    def markers_top(self, n = 5, groups=None, unique = True):
        '''
        return a flattened marker list
        groups: a list of cluster names (column names)
        '''
        markers = self.markers() 
        groups = markers.columns if groups is None else groups
        top = markers[groups].iloc[: n].values.T.flatten()
        if unique:
            top = pd.unique(top)
        
        return top
        

    def markers_detail(self):
        '''
        keys = ['params', 'scores', 'names', 'logfoldchanges', 'pvals', 'pvals_adj']
        column names are formated as '%s_%s' % (group, key)
            e.g. '0_names', '0_pvals', ..., '1_names', ...
        '''
        result = self.adata.uns['rank_genes_groups']
        groups = result['names'].dtype.names
        df = pd.DataFrame(
                {group + '_' + key: result[key][group]
                for group in groups for key in ['names', 'logfoldchanges', 'pvals', 'pvals_adj']})
            
        return df
    
    def save_markers(self, detailed = False, name = None):
        
        name = self.name if name is None else name
        if detailed:
            marker_file = '{}/{}_markers.csv'.format(self.resdir, name)
            self.markers_detail().to_csv(marker_file, index = False)
            print('markers with detailed info are saved to file `%s`'% marker_file)

        else:
            marker_file = '{}/{}_marker_names.csv'.format(self.resdir, name)
            self.markers().to_csv(marker_file, index = False)
            print('marker names are saved to file `%s`'% marker_file)
            
    def dot_de(self, topn=5, unique = True, groupby = None, 
               standard_scale = None,
               color_map='Reds',
               figsize=None, **kwds):
        gene_list = self.markers_top(topn, unique = unique)
        if groupby is None:
            groupby = 'leiden' if 'leiden' in self.adata.obs.keys() else 'louvain'
        
        sc.pl.dotplot(self.adata, gene_list, 
                      color_map=color_map,
                      standard_scale=standard_scale,
                      groupby=groupby, figsize=figsize, **kwds)

    #--------------------------------------------------------------------------   
    def GoThroughPipe(self, qc = True,
                      rmv_mito=False, mito_perc=None,
                      counts_level = 1000,
                      batch_key=None, 
                      n_top_genes=None, 
                      plot_hvgs=True,  
                      do_regress=None, #batch_key=None,
                      n_comps=50,
                      metric = 'correlation',
                      nneigh = None, 
                      de = None, plot_de=False,
                      cluster=True, 
                      save_middle = None, 
                      save_by_default = True, 
                      **kw_qc):
        
        print('Going Through the Pipeline...')
        
        if qc:
            self.QC(#min_genes = None, max_genes = None, 
                    plot = True, keep_results=True, **kw_qc)
            if rmv_mito:
                self.RmvMito(do_filter=None,  mito_perc = mito_perc, 
                             plot=False, keep_results=True)
                
            self.save(name = '%s_afterQC'%self.name, data_only=True)
            self.save_mtx(name = '%s_afterQC'%self.name)
        
        self.NormLog(counts_level = counts_level, frac = 0.1)
        self.HVGs(min_mean = None, min_disp = None, 
                  batch_key=batch_key, n_top_genes=n_top_genes,
                  plot = plot_hvgs, redo = None, keep_results=True)
        
        if save_middle is None:
            save_middle = (input('Save unscaled middle results?[y/n]').lower() == 'y') 
        if save_middle:
            self.save(name = '%s_beforeZscore'%self.name, data_only=True)
            self.save_hvgs()

        self.scale(do_regress=do_regress, groupby=batch_key)
            
        self.PCA(n_comps = n_comps, plot = True)#, do_regress='n_counts')
        
        if isinstance(metric, list):
            for mtc in metric:
                self.Neighbors(nneigh = nneigh, do_pca = True, n_pcs = None, metric=mtc)
                self.Umap(plot=True)
                self.Clustering(res = None, do_pca = True, n_pcs = None,
                                redo=False, pl_umap = True, keep_results=True)
#                                startover = True,
#                                nneigh = nneigh, metric=mtc)
        else: 
            self.Neighbors(nneigh = nneigh, do_pca = True, n_pcs = None, metric=metric)#'cosine') # 'euclidean' by default
            if cluster:
                self.Clustering(res = None, do_pca = True, n_pcs = None,
                                redo=None, pl_umap = True, keep_results=True)
            
        de = input('find markers for each cluster?[y/n]').lower() == 'y' if de is None else de
        if de:
            self.DE(save = True, plot = plot_de) 
#            self.save_markers(detailed=True)# already saved in DE phase
        
        
        self.summary()
        
        if save_by_default:
            name = self.name
        else:
            if input('save ALL results with default name[0]? or rename[1]?') == '1':
                name = input('enter the name to be saved as:')
            else:
                name = self.name
        self.save(name)
        
    
        
    #--------------------------------------------------------------------------   
    def AutoPipe(self, qc = True, min_genes = 100, max_genes = 5000,
                 save_raw_mtx = True,
                      rmv_mito=False, mito_perc = 0.05,
                      counts_level = 1000, log1first=False,
                      plot_hvgs=True, min_mean = 0.0125, min_disp = 0.25, 
                      batch_key = None, n_top_genes=None,
                      do_regress=False, 
                      cca=False, pca_on = None,
                      n_comps=50, scale_pca_by=None,
                      bknn = False,
                      metric = 'cosine', nneigh = 30, n_pcs=20, 
                      min_dist = 0.3,
                      cluster=True, clust_reso=1,
                      de = True, plot_de=True,
                      save_middle = True, 
                      save_by_default = True, 
                      **kwds):
        
                
        print('Going Through the Pipeline...')
        
        if qc:
            self.QC(min_genes = min_genes, max_genes = max_genes, 
                    plot = True, keep_results=True, )
            if rmv_mito:
                self.RmvMito(do_filter=True, mito_perc = mito_perc, 
                             plot=False, keep_results=True)
            
            self.save(name = '%s_afterQC'%self.name, data_only=True)
            if save_raw_mtx:
                self.save_mtx(name = '%s_afterQC'%self.name)
        else:
            print('Skipping QC phase ...')
        
        # Normalization and HVGs
        if log1first:
            self.NormLog_rev(counts_level = counts_level, frac = 0.1)
        else:
            self.NormLog(counts_level = counts_level, frac = 0.1)
        self.HVGs(min_mean = min_mean, min_disp = min_disp, 
                  batch_key=batch_key, n_top_genes=n_top_genes,
                  plot = plot_hvgs, redo = False, keep_results=True)
        self.save_hvgs()
        if save_middle:
#            self.save_hvgs()
            self.save(name = '%s_beforeZscore'%self.name, data_only=True)
            
        # analysis
        self.scale(do_regress=do_regress, groupby=batch_key)
        
        #--------dimensionality reduction------------
        if isinstance(cca, str): # `cca` canbe 'plssvd' or 'cca'
            self.CrossDimReduct(groupby=batch_key, method=cca, 
                                n_comps=n_pcs,)
            key_neigh = 'X_cca'
        elif pca_on is not None:
            self.PPCA(based_on=pca_on, key_add='X_pca', n_comps=n_comps)
            key_neigh = 'X_pca'
        else:
            self.PCA(n_comps = n_comps, plot = True, )
            if isinstance(scale_pca_by, str):
                self.GscalePCA(groupby=scale_pca_by, cover=True)
            key_neigh = 'X_pca'
        
        # ------------building graph--------------
        if bknn:
            _nb = len(self.adata.obs[batch_key].unique())
            sc.external.pp.bbknn(self.adata, batch_key, n_pcs=n_pcs,
                                 neighbors_within_batch=np.round(nneigh / _nb))
        else:
            self.Neighbors(nneigh = nneigh, do_pca = True, use_rep = key_neigh,
                       n_pcs = n_pcs, metric=metric)#'cosine') # 'euclidean' by default
        
        self.Umap(min_dist = min_dist)
        if cluster:
            self.Clustering(res = clust_reso, do_pca = True, n_pcs = n_pcs,
                            redo=False, pl_umap = True, keep_results=True)
            
            if de:
                self.DE(save = True, plot = plot_de)# .save_markers() # already saved in DE phase
#            self.save_markers(detailed=True)
        
        
        self.summary()
        self.save(self.name)        
        
    def AutoDownstream(self, used_genes,
                       log1first = False,
                       counts_level = 1000,
                       do_regress=False, 
                       batch_key=None, 
                       cca=False, pca_on=None,
                       n_comps=50, scale_pca_by=None,
                       metric = 'cosine', 
                       nneigh = 30, n_pcs=30, 
                       min_dist = 0.3,
                       cluster=True, 
                       de = True, plot_de=True,
                       save_middle = True, # igored
                       save_by_default = True,  # igored
                       **kwds):
        '''downstream analysis using given genes
        i.e. No Quality Control phase 
            and No HVG selection (directly use given gene list)
        '''
        # Normalization and set HVGs
        if log1first:
            self.NormLog_rev(counts_level = counts_level, frac = 0.1)
        else:
            self.NormLog(counts_level = counts_level, frac = 0.1)
        
        self.SetHVGs(used_genes, slim = True)
        self.save_hvgs()        
            
        # scaling (z-score)
        self.scale(do_regress=do_regress, groupby=batch_key)
        
        # dimensionality reduction
        if isinstance(cca, str): # `cca` canbe 'plssvd' or 'cca'
            self.CrossDimReduct(groupby=batch_key, method=cca, 
                                n_comps=n_pcs,)
            key_neigh = 'X_cca'
        elif pca_on is not None:
            self.PPCA(based_on=pca_on, key_add='X_pca', n_comps=n_comps)
            key_neigh = 'X_pca'
        else:
            self.PCA(n_comps = n_comps, plot = True, )
            if isinstance(scale_pca_by, str):
                self.GscalePCA(groupby=scale_pca_by, cover=True)
            key_neigh = 'X_pca'
        
        # building graph
        self.Neighbors(nneigh = nneigh, do_pca = True, 
                       use_rep = key_neigh,
                       n_pcs = n_pcs, metric=metric)#'cosine') # 'euclidean' by default
        
        self.Umap(min_dist = min_dist)
        if cluster:
            self.Clustering(res = 1, do_pca = True, n_pcs = n_pcs,
                            redo=False, pl_umap = True, keep_results=True)
            
        if de:
            self.DE(save = True, plot = plot_de)# .save_markers() # already saved in DE phase
#            self.save_markers(detailed=True)
        
        
        self.summary()
        self.save(self.name) 

    def CustomPipe(self, mtd_list, params = dict(), save = True):
#        print('Setting parameters...')
#################################################        
#        mtd_name = 'QC'
#        mtd = getattr(self, mtd_name)
#        mtd_params = params[mtd_name]
#        mtd(self, **mtd_params)
#################################################
        
        pass
        
    #---------------------------[ Visualization ]--------------------------

    def vis(self, color = None, key='umap', **kwds):
        if color is None:
            by1 = 'leiden' if 'leiden' in self.adata.obs.keys() else 'louvain'
            color = [by1, 'n_genes_disc'] 
            
        if key == 'umap':
            sc.pl.umap(self.adata, color = color, **kwds)
        elif key == 'tsne':
            sc.pl.tsne(self.adata, color = color, **kwds)
        elif key == 'pca':
            sc.pl.pca(self.adata, color = color, **kwds)
        

    def ButtomView(self, key='umap', groupby='batch', sz=0.1, tail=None,
                   disc=True, save=True, **kwds):
        
        tail = self.name if tail is None else tail
        plt_dir = os.path.join(self.figdir) if save else None
        fx.ButtomViewAdt(self.adata, key=key, groupby=groupby, sz=sz, tail=tail,
                         disc=disc, save=plt_dir, **kwds)
        
    def ButtomCompare(self, pair_labels, groupby, 
                      key='umap', sz=0.1, #cmap=cmap,
                      save=True, **kwds):
        plt_dir = os.path.join(self.figdir) if save else None
        fx.ButtomCompareAdt(self.adata, pair_labels, key=key,  
                         groupby=groupby, sz=sz, 
                         save=plt_dir, **kwds)
        
    #------------------------------[ SAVE ]-----------------------------
    
    def save_dims(self, key = 'umap', name=None, ftype='npy', **kwds):
        name = self.name if name is None else name
        dims_file = self.resdir / '{}_{}.{}'.format(key, name, ftype)
        X = self.adata.obsm['X_' + key]
        print('shape of the array:', X.shape)
        if ftype == 'npy':
            np.save(dims_file, X, **kwds)
        elif ftype in ['csv', 'txt', 'tsv']:
            sep = ',' if ftype == 'csv' else ' '
            np.savetxt(dims_file, X, delimiter=sep, **kwds)#, fmt='%.18e', header='',)
        print('Low dimensional data (%s) saved to file `%s`'% (key, dims_file))
        
    
    def save_info(self, name = None, save_hvgs = True):
        
        name = self.name if name is None else name
        param_file = '{}/{}_params.csv'.format(self.resdir, name)
        summ_file = '{}/{}_summary.csv'.format(self.resdir, name)
        
        self.df_param().to_csv(param_file, index = False)
        self.df_summ().to_csv(summ_file, index = False)
        print('parameters and data summary are saved to file `%s` and `%s`'% 
              (param_file, summ_file))
        
    def save_meta(self, name=None, columns=None, **kwds):
        name = self.name if name is None else name
        meta_file = self.resdir / '{}_metadata.csv'.format(name)
        self.adata.obs.to_csv(meta_file, index=True, columns=columns, **kwds)
        print('meta data saved to file `%s`'%meta_file)
        colnames = self.adata.obs_keys() if columns is None else columns
        print('\tcolumns: ', colnames)
        
    def save(self, name = None, data_only = False, save_info=True):
        '''
        if data_only:
            save only the '.h5ad' file
        else:
            besides, save `params`, `summary` and `meta_data` as well.
        '''
        name = self.name if name is None else name
        res_file = '{}/{}.h5ad'.format(self.resdir, name)
        self.adata.write(res_file) # automatically create the directory     
        print('results saved to file `%s`'%res_file)
        
        if not data_only: # save all
            if save_info:
                self.save_info(name = name)
            self.save_meta(name = name)
            
            
    def save_mtx(self, resdir = None, name = None, field = 'integer', **kwds):
        
        resdir = self.resdir if resdir is None else resdir
        name = self.name if name is None else name
        dirname = '{}/{}_mtx'.format(resdir, name)
        fx.saveNamedMtx(self.adata, dirname, field=field, **kwds)
        
#        if not os.path.exists(dirname):
#            os.makedirs(dirname)
#            
#        mtx_file = '%s/matrix.mtx'% (dirname)
#        bcd_file = '%s/barcodes.tsv'% (dirname)
#        gene_file = '%s/genes.tsv'% (dirname)
#        
#        genes = self.adata.var_names.to_frame(index = False, name=self.name)
#        barcodes = self.adata.obs_names.to_frame(index = False, name=self.name)
#        
#        sparse.save_npz(f'{dirname}/matrix.npz', self.adata.X.T)
#        io.mmwrite(mtx_file, self.adata.X.T, field = field, **kwds)
#        genes.to_csv(gene_file, index = False, header = False)
#        barcodes.to_csv(bcd_file, index = False, header = False)
#        
#        print('Matrix saved to directory `%s`'% dirname)
        
        
        
        



if __name__ == '__main__':
    
    pass
##----------------Datadir for amph------------------
#    tt = '_embryos'
#    #tt = '_adults'
##    DATADIR = 'E:/Users/xyliu/data/amphioxus%s'%tt + '/afterQC'
#    DATADIR = 'single_analysis_keep1'
#    snames = ['E%d'%i for i in range(1, 16)]
#    resdir = 'single_analysis_try2'
##----------------Datadir for mouse------------------
    
#    DATADIR = 'E:/Users/xyliu/data/RNF200'   
#    snames = 'HE KO_1 WT'.split()
##    datadir = [sn + '_gene_exon_tagged' for sn in snames]
#    resdir = 'D:/Users/xyliu/003/' + 'afterQC'
#    sname = snames[1]
#    adata = adataFromRaw('%s/%s_star_gene_exon_tagged'%(DATADIR, sname), ftype = '.tsv')
##
###    sname = snames[3]
###    adata = sc.read_h5ad('{}/{}_afterQC.h5ad'.format(DATADIR, sname))
###    adata = sc.read_h5ad('{}/{}.h5ad'.format(DATADIR, sname))
#    ph1 = PipeHandler(ph.adata, sname, resdir)
#    ph = PipeHandler(Adata, resdir='RQ-combined', name='ref_KO_1')

##    for i in range(3):
##        
##        sname = snames[i]
##        adata = sc.read_h5ad('{}/{}_afterQC.h5ad'.format(DATADIR, sname))
##
##        ph = PipeHandler(adata, sname, resdir)
##        ph.save_mtx()
#        
######################[ go through pipeline ]########################### 
#        
#        print('\n', '=' * 80, sep = '')
#        print('\tCheck for `%s`'% sname)
#        print('=' * 80)       
#        
#        ph = PipeHandler(adata, sname, resdir)
#        
#        ph.GoThroughPipe(qc = True, rmv_mito=True, save_middle = True)
##        ph.save()

        
        
#####################[ Visualization ]###########################     
#        if i == 3: 
#            groupby = 'louvain_R'
#        else:
#            groupby='louvain'
##        sc.pl.stacked_violin(ph.adata, gene_list, groupby='louvain', rotation=90)
#        ph.vis([groupby, 'n_genes_disc'])
#        
#######################[ DE and markers ]###########################     
##     #    sc.tl.rank_genes_groups(ph.adata, 'louvain', method = 't-test_overestim_var')#
##        # method : `{'logreg', 't-test', 'wilcoxon', 't-test_overestim_var'}
#        gene_list = ph.markers_top(3)
#        sc.pl.dotplot(ph.adata, gene_list, groupby='louvain')#, figsize = (8, 4))
#        
#        sc.pl.rank_genes_groups_heatmap(adata, vmin = 0, n_genes=3)
#        sc.pl.rank_genes_groups(ph.adata)
    #    sc.pl.rank_genes_groups_violin(ph.adata, n_genes=8)#, scale = 'count')
        
    #    sc.pl.violin(ph.adata, gene_list[::3], groupby='louvain', scale = 'count', figsize = (3, 4))
#####################[  ]###########################


#    ph.adata.obs['zero_hvg_counts'] = pd.Categorical(ph.adata.obs['hvg_counts'] == 0 )
#    sc.pl.heatmap(adata, adata.var_names, groupby = 'louvain', dendrogram=True)
    
    
#    ph.QC()
#    ph.NormLog(frac = 0.02)

#    ph.save('combined_afterBC')
    
#    ph.HVGs(plot = True, redo = None, keep_results=None)
#    ph.ZeroFilter()
#    ph.PCA()
#    sc.external.pp.bbknn(ph.adata, batch_key = 'batch', n_pcs = 40, )
#    
#    
#    sc.tl.tsne(ph.adata)
#    sc.tl.diffmap(ph.adata)
#    sc.pl.diffmap(ph.adata, color = ['periodc'], components=['1,2','1,3'])
#    sc.pl.pca(ph.adata, color = ['periodc'], components=['1,2','1,3'])
#    
#    ph.Clustering(res = 2)#keep_results=True)
#    ph.adata.obs_names
#    
#    
#    ph.vis(['louvain', 'period'])    
#    sc.pl.umap(ph.adata, color = ['louvain', 'period'])#, palette = 'tab20b'
#    sc.pl.tsne(ph.adata, color = ['louvain', 'periodc'])#, palette = 'tab20b'
#    
#    ph.save('combined')
#    sc.tl.draw_graph(ph.adata)


#sc.tl.louvain(ph.adata, resolution = 0.6, restrict_to = ('louvain', ['1']), key_added = 'louvain')
#ph.vis(['louvain_R', 'n_genes_disc']) 
#ph.DE(groupby='louvain_R') #sc.tl.rank_genes_groups(ph.adata, 'louvain_R', method = 't-test_overestim_var')
#    ph.summary_d()
#    ph.df_summ()
#    print(ph.params)
#    ph.plotsumm()

    
    
#    ph1 = PipeHandler(ph.adata, name = 're_' + sname)
#    ph1.ZeroFilter()
#    ph1.summary()
#    ph1.GoThroughPipe()
#    ph1.adata.X[:5, :5].toarray()

#sc.pl.pca(ph.adata, color = ['periodc'], components=['1,2','1,3'])
#sc.pl.umap(ph.adata, color = ['louvain', 'periodc'])
#sc.pl.diffmap(ph.adata, color = ['periodc'], components=['1,2','1,3'])
#sc.pl.draw_graph(ph.adata)







