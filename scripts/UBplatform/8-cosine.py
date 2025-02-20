import sys
import os
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import squidpy as sq
import networkx as nx
from tqdm import tqdm
import matplotlib.pyplot as plt
import zonist.method as zm


SAMPLE = "UBcalibr2"
#"UBcalibr2"
print(SAMPLE)
mydata_dir = "/home/johanna/spatial_thesis/spatial_metabolomics/data"
#"/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"
subfolder = "UBplatform"
#"/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"



adata = ad.read_h5ad(os.path.join(mydata_dir, subfolder,
                                         f"{SAMPLE}-KEGG_adata.h5ad"))


scores = ad.read_h5ad(os.path.join(mydata_dir, subfolder,
                                         f"{SAMPLE}-PathwayScores.h5ad"))


# TODO rerun all from step 1 to avoid "Pathways in Cancer"


print(scores)
# print(scores.var['pathway_name'].unique())
# for pw in scores.var['pathway_id']:
#     pw_name = scores.var.loc[scores.var['pathway_id'] == pw, "pathway_name"].item()
#     sq.pl.spatial_scatter(scores, color=[pw],
#                                   img=False,
#                           title=pw_name)
#     plt.show()
#
"""
sq.gr.spatial_neighbors(scores, n_neighbors=10,
                        use_rep=None,
                        method='umap', metric='cosine')
"""
sc.pp.neighbors(scores, n_neighbors=600, n_pcs=None, use_rep=None,
                method='umap', metric='cosine', random_state=0)
# creates .obsp["connectivities'] with the cosine dists

sc.tl.leiden(scores,
             adjacency=scores.obsp["connectivities"], # is a squared n_obs x n_obs matrix: can be huge!
             key_added="leiden_clusters", resolution=0.1, flavor='igraph',
             n_iterations=2)
print(scores)
sq.pl.spatial_scatter(scores, color=['leiden_clusters'], img=False)
plt.show()


