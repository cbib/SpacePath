import sys
import os
import numpy as np
import pandas as pd
import anndata as ad
import squidpy as sq
import networkx as nx
from tqdm import tqdm
import matplotlib.pyplot as plt
import zonist.method as zm
from scipy.sparse import csr_matrix

from scripts.bank4zonist.quantifs_methods import spatial_standard_computations
from scripts.bank4zonist.visualization import plot_standard_computs

sys.path.append(os.path.join(os.path.dirname(__file__),
                             '../bank4zonist'))
from visualization import *

SAMPLE = "UBcalibr2"
#"UBcalibr2"
print(SAMPLE)
mydata_dir = "/home/johanna/spatial_thesis/spatial_metabolomics/data"
#"/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"
subfolder = "UBplatform"


plot_folder = "../../../figures/platform"

assert os.path.exists(plot_folder), "Cannot continue, the plot_folder does not exist"

k_dict= {"UBcalibr1": [7, 4],
         "UBcalibr2": [5, 3],
         "UBcalibr3": [6, 4],
         "UBcalibr4": [8, 5],
         "UBcalibr5": [6, 4],
         }


##  read objects

# mz
adata_mz = ad.read_h5ad(filename=os.path.join(
    mydata_dir, subfolder, f'{SAMPLE}.h5ad'))

# metaboli
adata = ad.read_h5ad(os.path.join(mydata_dir, subfolder,
                                         f"{SAMPLE}-KEGG_adata.h5ad"))
# pathways
scores_adata = ad.read_h5ad(os.path.join(mydata_dir, subfolder,
                                         f"{SAMPLE}-PathwayScores.h5ad"))
print("scores   :::  ", scores_adata)

## plots

# # mz
# adata_mz = adata_mz[: , scores_adata.uns[
#     "current_pathways_df"]["mz"].astype(str).unique() ]  # first subset good mz
# print(adata_mz)
#
# adata_mz = spatial_standard_computations(
#     adata_mz,
#     k=k_dict[SAMPLE][0], top_variable_features_fraction=0.5)
#
# plot_standard_computs(adata_mz,
#                       plot_folder=plot_folder, prefix=f"{SAMPLE}-MZ")
#
# if SAMPLE == "UBcalibr2":
#     print("doing plots of mz that are highly spatially variable")
#     sq.gr.spatial_neighbors(adata_mz, n_neighs=20)
#     sq.gr.spatial_autocorr(adata_mz, mode="moran", n_perms=100, n_jobs=1)
#     adata_mz.uns['moranI']['pathway_id'] = adata_mz.uns[
#         'moranI'].index.tolist()
#     top_spatial_features = adata_mz.uns['moranI'].head(
#         10
#     )["pathway_id"].tolist()
#     for i_top_mz in top_spatial_features:
#         moran_i = adata_mz.uns['moranI'].loc[i_top_mz, "I"].round(3).item()
#         tmpstr = os.path.join(plot_folder,
#                               f'{SAMPLE}-{np.around(float(i_top_mz), 2)}.png')
#         sq.pl.spatial_scatter(
#             adata_mz, color=[str(i_top_mz)], img=False,
#             title=f"m/z= {np.around(float(i_top_mz), 4)} (Moran's I={moran_i})",
#                               save=os.path.abspath(tmpstr))
#
#
#
# # metaboli
#
# adata = spatial_standard_computations(
#     adata,
#     k=k_dict[SAMPLE][0],
#     top_variable_features_fraction=0.8)
#
# plot_standard_computs(adata,
#      plot_folder=plot_folder, prefix=f"{SAMPLE}-KEGGmetabolites")


# pathways
# note : remember this is bad idea to replace .X by the shifted version, only for test, do not save corrupted object!
# decision: never again do shifting (if so, create shifted_adata).

print(scores_adata)
top_n_by_I = 7 # also did 10
sq.gr.spatial_neighbors(scores_adata, n_neighs=20)  # needed for Moran's I
sq.gr.spatial_autocorr(scores_adata, mode="moran", n_perms=100, n_jobs=1)
scores_adata.uns['moranI']['pathway_id'] = scores_adata.uns['moranI'].index.tolist()
top_spatial_features = scores_adata.uns['moranI'].head(
      top_n_by_I
)["pathway_id"].tolist()
print("subsetting with only spatially relevant features (using Moran'I metrics)")
scores_adata = scores_adata[:, top_spatial_features].copy()

scores_adata = spatial_standard_computations(
    scores_adata,
    pca_ncomps=3,
    k=k_dict[SAMPLE][1],
    as_log=False,   # important !:  in this way,no need to "shift"
    top_variable_features_fraction=1) # this will select all features anyways

plot_standard_computs(scores_adata,
    plot_folder=plot_folder,
    prefix=f"{SAMPLE}-allPathwayScores-{top_n_by_I}topI")
