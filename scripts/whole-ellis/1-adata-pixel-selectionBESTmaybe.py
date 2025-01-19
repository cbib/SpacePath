import os
import numpy as np
import anndata as ad
import scanpy as sc
import squidpy as sq

from sklearn.utils import sparsefuncs

import matplotlib.pyplot as plt

"""
Custom script to discard the pixels that do not fall
in the tissue, but in the background.
"""
# TODO : test this

# for the whole-ellis dataset

# 0) all variables, strings and dictionaries, to be used
mydata_dir = "/home/johanna/spatial_thesis/spatial_metabolomics/data"
subdir = "whole-ellis"
#"/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"
quantif_file = "whole-ellis_mx.h5ad"

out_file =  "whole-ellis_mx_f1.h5ad"


# 1) given data


# 1.2) matrix of quantifications (==intensities)
print("Opening the quantifications")
adata = ad.read_h5ad(filename=os.path.join(
    mydata_dir, subdir, quantif_file))


##  add attributes to be spatially plottable

# coordinates_arr = generate_coordinates(adata.obs["obs"].tolist())
# adata = ad.AnnData(adata.X,
#                     obs = adata.obs, var= adata.var,
#                     obsm = {"spatial": coordinates_arr},
#                     varm = {"metabolites": adata.var["mz"].to_numpy()}
#                     )


print("\nGet rid of pixels that are not part of the tissue:")
#  use Spatial domains as in https://www.sc-best-practices.org/spatial/domains.html
print("  using most variable features for dimensional reduction+clustering")
means_vars, variance_computed = sparsefuncs.mean_variance_axis(adata.X, axis=0)
adata.var['variance'] = variance_computed

tmp = adata.var.sort_values(by='variance', ascending=False)
sampled_features = tmp.head(500)["mz"].astype(str).tolist()
del tmp
miniadata = adata[:, sampled_features].copy()

print("  computing spatial domains to identify out-of-tissue pixels...")

leiden_resol = 0.4
sc.pp.pca(miniadata)
sc.pp.neighbors(miniadata)  # nearest neighbor graph
sc.tl.leiden(miniadata, resolution=leiden_resol, key_added=f"leiden{leiden_resol}")
nn_graph_genes = miniadata.obsp["connectivities"]
# spatial proximity graph
sq.gr.spatial_neighbors(miniadata)
nn_graph_space = miniadata.obsp["spatial_connectivities"]
print("computing leiden clustering on joint graph...")
alpha = 0.6  # weight for the spatial (coords) nngraph
joint_graph = (1 - alpha) * nn_graph_genes + alpha * nn_graph_space
sc.tl.leiden(miniadata,
             resolution=leiden_resol,
             adjacency=joint_graph,
             key_added=f"squidpy_domains(leiden{leiden_resol})")
print("visualizing the spatial plot of domains before decisions")
sq.pl.spatial_scatter(
    miniadata, color=[f"squidpy_domains(leiden{leiden_resol})"],
    img=False, title=f"before filtering out-of-tissue pixels. (leiden{leiden_resol})")
plt.show()

not_tissue_spatial_domains = [2,7,8]   # DECISION (requires viz first)
not_tissue_spatial_domains = [str(x) for x in not_tissue_spatial_domains ]
pixels_to_keep = miniadata.obs.loc[
    ~miniadata.obs[f"squidpy_domains(leiden{leiden_resol})"].isin(
        not_tissue_spatial_domains), "obs"].tolist()
miniadata = miniadata[pixels_to_keep, :]
domains_list = miniadata.obs[f"squidpy_domains(leiden{leiden_resol})"].tolist()
del miniadata

print("  applying decision (this is a second run)")

adata = adata[pixels_to_keep, :]
adata.obs["spatial_domain"] = domains_list
print("visualizing the spatial plot of domains after decision")
sq.pl.spatial_scatter(
    adata, color=["spatial_domain"],
    img=False, title=f"AFTER filtering out-of-tissue pixels (leiden{leiden_resol})")
plt.show()

outpath = os.path.join(mydata_dir, subdir, out_file)
adata.write(outpath, compression=None, compression_opts=None)
