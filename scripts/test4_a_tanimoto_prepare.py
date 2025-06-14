import os
import pandas as pd
import numpy as np
import zonist.data as zd
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import squidpy as sq
from local_helpers.local_helpers import provide_annotation_dict
import zonist.processing.measures_methods as zmm




mydata_dir = "../../data/"
SAMPLE = "TD-Lipido-635"
subfolder = SAMPLE

with open("../../data/TD-Lipido-635/TD-Lipido-635-clusters-initial.txt", "r") as f:
    initial_clusters = f.readlines()

initial_clusters.pop(0)
initial_clusters = [float(t.strip().split(",")[1]) for t in initial_clusters]


adata_mz = ad.read_h5ad(filename=os.path.join(
            mydata_dir, subfolder, f'{SAMPLE}.h5ad'))

print(adata_mz.uns["history"])
print("MANDATORY: take away bad quality m/z") # todo: how to make an alert ?
adata_mz = adata_mz[:,
           adata_mz.var.loc[adata_mz.var['final_quality'] == "good",
           "mz"].astype(str).tolist()]

print(adata_mz)
adata_mz.obs['kmeans_5'] = np.array(initial_clusters, dtype=str)

sq.pl.spatial_scatter(adata_mz, color=[f'kmeans_5'],
                             img=False, title="hehe")
plt.show()
print(ababado)

sc.pp.log1p(adata_mz)
print("computing the differential abundances of all the m/z features, cluster to cluster, for all clusters")
sc.tl.rank_genes_groups(adata_mz, groupby='kmeans_5', groups='all',
                            reference='rest',
                            method='wilcoxon' ,  # use logreg if ugly
                            key_added="top_feats_by_cluster",
                            corr_method='benjamini-hochberg')


differential_results_by_cluster = sc.get.rank_genes_groups_df(adata_mz,
                                         key="top_feats_by_cluster",
                                         group=None )  # returns all groups)

differential_results_by_cluster.to_csv("differential_results.csv", index=False)

