# scvoting
This is a voting method to obtain the density (smoothed votes) of a 10 sc-RNASeq dataset on a common reference sc-RNASeq data's landscape which can be viewed using umap, fdg etc

The method works with scanpy annData objects only and contains two parts.

## Projecting onto reference landscape
This package follows the Xiaonan Wang's method of projecting a 10x sc-RNASeq dataset onto a reference sc-RNASeq data. 
After projection

## Voting and smoothing
* For every cell in the dataset its NN in the reference dataset gets 1 vote.
* Smooth out the votes for each cell in the reference dataset by sharing them equally among 100 NN in it. So if its vote is 40 you share 40/100 among each of them.
* Compute this for WT and mutant datasets. 
* Plot their difference or logFC etc on the reference landscape (UMAP)

## Example run

```python
from scvoting import Projector
# Initantiation
votp = Projector(adata = adata, ref_adata=adata_niki, ref_hvg=niki_hvg, npcs=25)
# Proj onto ref
votp.project_on_ref()
# Plot PCAs (works only with X11, use jupyter on remote servers)
votp.plotPCA()
# Compute distances between cells in two datasets in PCA space
votp.pwdist()
# Voting starts
votp.voting(tot_votes = 500000)
# Compute distances between cells in reference datasets
votp.pwdist(rtype = 'ref')
# Smooth the voting for better visualisation
votp.vote_smoothing()
# Smoothed votes
results = votp.votes_smooth
```
