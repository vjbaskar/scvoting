#!/usr/bin/env python3

import scanpy as sc
import pandas as pd
import numpy as np
import yaml
import time
import anndata
import sklearn
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt


from sklearn.metrics.pairwise import euclidean_distances
import matplotlib.pyplot as plt
from collections import Counter
from collections import defaultdict

from sklearn.metrics.pairwise import euclidean_distances
import matplotlib.pyplot as plt
from collections import Counter
from collections import defaultdict

class Projector:
    """
    Projector class:
    ## Paired Voting System
    * For every cell in the dataset its NN in the reference dataset gets 1 vote.
    * Smooth out the votes for each cell in the reference dataset by sharing them equally among 100 NN in it. So if its vote is 40 you share 40/100 among each of them.
    * Compute this for WT and mutant datasets. 
    * Plot their difference or logFC etc on the reference landscape (UMAP)
    """
    def __init__(self, adata, ref_adata, ref_hvg, npcs = 25):
        self.adata = adata
        self.ref_adata = ref_adata
        self.ref_hvg = ref_hvg
        self.npcs = npcs
        self.dist_mtx = dict()
    
    def project_on_ref(self):
        """
        Input
        -----
        adata: data to be projected. Needs to have adata.raw instantiated. Ideally it should have all the genes not just hvg.
        ref_adata: reference dataset such as Niki or Nesterowa
        ref_hvg: HVG in reference dataset
        npcs: number of PC
        """
        adata_niki = self.ref_adata
        niki_hvg = self.ref_hvg
        adata = self.adata
        adata = anndata.AnnData(X=np.exp(adata.raw.X.todense())-1,
                            obs=adata.obs, 
                            var=adata.raw.var, 
                            obsm=adata.obsm,
                            uns = adata.uns)

        print("Finding overlapping gene set")
        OLG = np.intersect1d(niki_hvg, adata.var_names)
        print("Genes common to both:")
        print(len(OLG))
        adata = adata[:,OLG].copy()
        adata_niki = adata_niki[:,OLG].copy()

        print("Normalising & log'ing data")
        sc.pp.normalize_per_cell(adata, counts_per_cell_after=10000)
        sc.pp.normalize_per_cell(adata_niki, counts_per_cell_after=10000)
        sc.pp.log1p(adata)
        sc.pp.log1p(adata_niki)

        # scale them together
        print("Combining & Scaling them together")
        data_comb = adata.concatenate(adata_niki)
        sc.pp.scale(data_comb)

        adata = anndata.AnnData(X=data_comb[data_comb.obs['batch'] == '0',:].X, 
                                obs=adata.obs, 
                                var=adata.var, 
                                obsm=adata.obsm, 
                                uns=adata.uns)

        adata_niki = anndata.AnnData(X=data_comb[data_comb.obs['batch'] == '1',:].X, 
                                     obs=adata_niki.obs, 
                                     var=adata_niki.var, 
                                     obsm=adata_niki.obsm, 
                                     uns=adata_niki.uns)
        # PCs for Niki's data
        print("Computing PCA for reference dataset")
        pca_ = sklearn.decomposition.PCA(n_components=50, svd_solver='auto', random_state=0)
        pca_.fit(adata_niki.X)

        # Project data and niki's data onto Niki's PCs
        print("Projecting data on the reference PCA space")
        adata_proj = pca_.transform(adata.X)
        adata_niki_proj = pca_.transform(adata_niki.X)
        self.adata_proj = adata_proj
        self.ref_proj = adata_niki_proj
        self.pca = pca_
    
    def plotPCA(self):
        """
        Some PCA plots. First plot is used to alter pcs
        Second plot is used to determine if projection worked well
        """
        pca_ = self.pca
        adata_niki_proj = self.ref_proj
        adata_proj =  self.adata_proj
        print("Plotting explained variance. Useful to alter npcs parameter")
        plt.plot(pca_.explained_variance_)

        fig = plt.figure()
        print("Plotting the projected cells. They should fall within the reference phases")
        ax1 = fig.add_subplot(111)
        ax1.scatter(adata_niki_proj[:,0], adata_niki_proj[:,1], c='black', alpha=0.5)
        plt.xlabel('PCA1')
        plt.ylabel('PCA2')
        ax1.scatter(adata_proj[:,0], adata_proj[:,1], c='red', alpha=0.5)
        plt.show()
    
    def pwdist(self, rtype='across'):
        """
        Computes pairwise distance between adata and reference data.
        type:   across - between adata (row) and reference data (col)
                adata - between cells of adata
                ref - between cells of reference data

        returns: self.dist_mtx - a dict of {type: numpy mtx}        
        """
        adata_proj = self.adata_proj
        ref_proj = self.ref_proj
        npcs = self.npcs
        if rtype == "across":
            D_sub= euclidean_distances(adata_proj[:,0:npcs], ref_proj[:,0:npcs])
        if rtype == "adata":
            D_sub= euclidean_distances(adata_proj[:,0:npcs])
        if rtype == "ref":
            D_sub= euclidean_distances(ref_proj[:,0:npcs])
        self.dist_mtx[rtype] = D_sub
        self.type_to_work = rtype
    
    def voting(self, tot_votes = 500000):
        """
        voting:
        for each cell in adata, find NN in ref data. 
        add a vote to the each of the cell in NN.
        Needs pwdist(rtype='across') already performed.
        
        """
        adata_niki = self.ref_adata
        D_sub = self.dist_mtx['across']
        adata_proj = self.adata
        tot_cells_adata = adata_proj.shape[0]
        tot_cells_ref = adata_niki.shape[0]
        votes_per_cell = np.int(np.floor(tot_votes/tot_cells_adata))
        print(f"votes per cell = {votes_per_cell}")
        
        ref_names = adata_niki.obs_names
        votes = np.zeros(tot_cells_ref)
        for i in range(D_sub.shape[0]):
            CellDis = D_sub[i,:]
            CellDis_sorted = np.argsort(CellDis)[:votes_per_cell]
            votes[CellDis_sorted] += 1
        print(f"Run for {i} cells in the data") 

        self.tot_votes = tot_votes
        self.votes_per_cell = votes_per_cell
        self.votes_raw = votes
    
    def vote_smoothing(self, NN=100, hscNN = 20, celltype_col="CellSubType"):
        
        adata_niki = self.ref_adata
        ref_dist_matx = self.dist_mtx['ref']
        votes = self.votes_raw.copy()
        
        ref_nn = []
        nn_no = NN
        hsc_nn_no = hscNN
        cell_types = list(adata_niki.obs[celltype_col])
        for i in range(adata_niki.shape[0]):
            CellDis = ref_dist_matx[i,:]
            celltype = cell_types[i]
            if celltype == "HSCs":
                lim = hsc_nn_no
            else:
                lim = nn_no
            CellDis_sorted = list(np.argsort(CellDis)[:lim])
            ref_nn.append(CellDis_sorted)
        
        for i in range(len(ref_nn)):
            nn = ref_nn[i]
            K = len(nn)
            votes_for_i = votes[i]
            to_share = votes_for_i/K
            votes[nn] += to_share
        
        self.NN = 100
        self.hsc_nn = 20
        self.votes_smooth = votes


        
if __name__ == "__main__":
    print("Reading in niki data")
    adata_niki = sc.read_h5ad("/home/vs401/rds1/refdata/10x/mm/niki_landscape/niki_passQC_norm_10K.h5ad")
    #minmax(adata_niki)
    adata_niki.var_names_make_unique()
    adata_niki_all = adata_niki.copy()
    print("Reading hvg data")
    niki_hvg = np.genfromtxt('/home/vs401/rds1/refdata/10x/mm/niki_landscape/gene_names.txt', dtype=str)
    len(niki_hvg)
    adata = sc.read('write/cluster.h5ad')
    temp = Projector(adata = adata, ref_adata=adata_niki, ref_hvg=niki_hvg, npcs=25)
    temp.project_on_ref()
    temp.pwdist()