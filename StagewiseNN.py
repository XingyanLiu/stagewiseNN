# -*- coding: utf-8 -*-
"""
Created on Thu May 28 22:38:56 2020

@author: Administrator
"""
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.preprocessing import normalize
from sklearn.decomposition import PCA
from sklearn.neighbors import BallTree, KDTree


def stagewise_knn(X, stage_lbs, stage_ord, k=2, 
                  leaf_size=None, sigma=1.414, 
                  norm_sample=True, 
                  norm_sample_len = 1,# 3
                  do_pca=False,
                  n_pcs = 30,
                  leaf_size_max = 30,
                  norm_dists=False,
                  precis = 0.1,
                  ):
    '''
    X: the whole data matrix
        np.array; shape = (n_samples, n_vars)
    stage_lbs: 
        list-like, np.array; shape = (n_samples,)
    stage_ord: 
        list-like; shape = (n_stages)
        indicates the order of stages
    
    ================
    other parameters
    ================
    n_pcs: 
        int or a list of integers specifying the number of PCs used in each stage.
    precise: 
        for adjusting the accuracy of the tree-leaves
    '''
    # setting parameters
    leaf_size_input = leaf_size
    n_pcs = [n_pcs] * len(stage_ord) if isinstance(n_pcs, int) else n_pcs
    k = [k] * len(stage_ord) if isinstance(k, int) else k
    sigma = [sigma] * len(stage_ord) if isinstance(sigma, int) else sigma
    
    # checking inputs
    stage_lbs = np.asarray(stage_lbs)
    cell_per_stg = pd.value_counts(stage_lbs)
    print('Number of samples in each stage\n', cell_per_stg)
    print('Stage order:', stage_ord)
#    if norm_sample:# pretend cosine matric
#        X = normalize(X, axis=1) * norm_sample_len#  * 3: avoid Precision overflow
    
    n_points = X.shape[0]
    Inds = np.arange(n_points)
    
    iis = []
    dds = []
    jjs = []
    conns = []
    
    for i, stg1 in enumerate(stage_ord):
        if i == 0:
            stg0 = stg1
        else:
            stg0 = stage_ord[i - 1]
    #    stg0, stg1 = 'E14', 'E15'
        inds_ref = stage_lbs == stg0 # earlier stage
        inds_que = stage_lbs == stg1 # later stage

        X_que = X[inds_que, :]
        X_ref = X[inds_ref, :]   
        if do_pca:
#            print(f'performing PCA on the later stage {stg1}, and transforming.')
#            pca = PCA(n_components=n_pcs[i], )
#            X_que = pca.fit_transform(X_que)
#            X_ref = pca.transform(X_ref)
            print(f'performing PCA on stage {stg0} + {stg1}, n_pcs={n_pcs[i]}.')
            pca = PCA(n_components=n_pcs[i], )
            X_pca = pca.fit_transform(np.vstack([X_ref, X_que]))
            X_ref = X_pca[: cell_per_stg[stg0], :]
            X_que = X_pca[cell_per_stg[stg0]: , :]
            
        
        if leaf_size_input is None:
            # n: Geometric average
            n = np.sqrt(cell_per_stg[stg0] * cell_per_stg[stg0])
            print('n =', n)
            leaf_size = min([int(np.sqrt(n) * precis), leaf_size_max])
            leaf_size = max([leaf_size, 1])
#        k = 
        print(f'searching neighbors of stage {stg1} in {stg0} (k={k[i]})')
        print('leaf-size:', leaf_size)
        
        if norm_sample:# pretend cosine matric
            X_ref = normalize(X_ref, axis=1) * norm_sample_len#  * 3: avoid Precision overflow
            X_que = normalize(X_que, axis=1) * norm_sample_len#  * 3: avoid Precision overflow
#        print(X_ref.shape)
#        print(X_que.shape)
        tree = KDTree(X_ref, leaf_size=leaf_size, )#metric=metric)
        dist, inds0 = tree.query(X_que, k=k[i])
        inds = np.take(Inds[inds_ref], inds0)
        
        # normalized distances
        if norm_dists:
            dist_n = normalize_dists(dist)
#            conn = normalize_dists(connect_heat_kernel(dist_n, sigma[i]))
            conn = normalize(connect_heat_kernel(dist_n, sigma[i]),\
                             'l1', axis=1) * np.sqrt(k[i])
        else:
            conn = normalize(connect_heat_kernel(dist, sigma[i]), \
                             'l1', axis=1)* np.sqrt(k[i])
        #lbs = np.take(class_lbs[inds_ref], inds0, )
        
        iis.append(np.repeat(Inds[inds_que], k[i]))
        jjs.append(inds.flatten())
        dds.append(dist.flatten())
        conns.append(conn.flatten())
        
    print('Done in KNN searching, constructing the results distances and connections')
    iis = np.concatenate(iis)
    jjs = np.concatenate(jjs)
    dds = np.concatenate(dds)
    conns = np.concatenate(conns)
    distmat = sparse.coo_matrix((dds, (iis, jjs)), shape=(n_points, n_points))
    distmat += distmat.T
    print(f'distance matrix of shape {distmat.shape} and {distmat.nnz} non-zeros')
    connect = sparse.coo_matrix((conns, (iis, jjs)), 
                                shape=(n_points, n_points))
    connect += connect.T
    print(f'connection matrix of shape {connect.shape} and {connect.nnz} non-zeros')
    return distmat, connect 



def normalize_dists(d):
    '''
    d: 2-D np.array, normalization will perform on each row
    '''    
    return np.array(list(map(_normalize_dists_single, d)))
    

def _normalize_dists_single(d):
    '''
    d: 1-D np.array
    '''
    d1 = d[d.nonzero()]
    vmin = d1.min() * 0.99
#    vmax = d1.max() * 0.99
    med = np.median(d1)
    d2 = (d - vmin) / (med - vmin)
    d2[d2 < 0] = 0
    return d2
    
    
    


def connect_heat_kernel(d, sigma):
    return np.exp(-np.square(d/sigma))



if __name__== '__main__':
    
    import seaborn as sns
    from sklearn.decomposition import PCA
    data = pd.read_csv('test_data/cluster_avgs_717.csv', index_col=0).iloc[:, 2:]
    
    X = data.values.T
    X = PCA(n_components=30).fit_transform(X)
    
    stage_lbs = [c.split('_')[0] for c in data.columns.values]
    stage_ord = ['E' + str(i+1) for i in range(2, 15)]
    
    k = 2
    sigma = 1.414
    distmat, connect = stagewise_knn(X, stage_lbs, stage_ord, 
                                     k=k, leaf_size=None, sigma=sigma, norm=True)
    
#    sns.heatmap(distmat.toarray())
    sns.heatmap(connect.toarray())
    
    # =============== UMAP embedding the connectivities =====================
    import scanpy as sc
    adata = sc.AnnData(data.T)
    adata.obs['stage'] = stage_lbs
    
    adata.uns['neighbors'] = {
            'params': {'n_neighbors': k, 'method': 'umap', 'metric': 'euclidean'},
            'connectivities': connect, #.to_csr(),
            'distances': distmat,}
    
    sc.tl.umap(adata, min_dist=0.1)
    sc.pl.umap(adata, color='stage', size=100, palette='Spectral')


