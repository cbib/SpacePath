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

SAMPLE = "test635_80-400"

print(SAMPLE)
mydata_dir = "/home/johanna/spatial_thesis/spatial_metabolomics/data"
#"/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"
subfolder = "gb-brain"

k_dict= {"test635_80-400": [5, 3],
         "???": [4, 3]}

plot_folder = "../../../figures/gb-brain"

assert os.path.exists(plot_folder), "Cannot continue, the plot_folder does not exist"


##  read objects
#
# # mz
adata_mz = ad.read_h5ad(filename=os.path.join(
    mydata_dir, subfolder, f'{SAMPLE}.h5ad'))

# problematics = ['259.022442561', '346.055808353']
# for malo in problematics:
#     sq.pl.spatial_scatter(adata_mz, color=[malo], img=False)
#     plt.show()

#
# metaboli
adata = ad.read_h5ad(os.path.join(mydata_dir, subfolder,
                                         f"{SAMPLE}-KEGG_adata.h5ad"))
# pathways
scores_adata = ad.read_h5ad(os.path.join(mydata_dir, subfolder,
                                         f"{SAMPLE}-PathwayScores.h5ad"))
#sq.pl.spatial_scatter(scores_adata, color=['hsa00010'], img=False)

print("adata (KEGG)   :::  ", adata)

## plots

# mz
adata_mz = adata_mz[: , scores_adata.uns[
    "current_pathways_df"]["mz"].astype(str).unique() ].copy()  # first subset good mz
print("adata mz : ", adata_mz)

adata_mz = spatial_standard_computations(
    adata_mz, pca_ncomps=6,
    k=k_dict[SAMPLE][0], top_variable_features_fraction=0.5)

plot_standard_computs(adata_mz,
                      plot_folder=plot_folder, prefix=f"{SAMPLE}-MZ")


adata = spatial_standard_computations(
    adata, pca_ncomps=6,
    k=k_dict[SAMPLE][0],
    top_variable_features_fraction=0.5)

plot_standard_computs(adata,
     plot_folder=plot_folder, prefix=f"{SAMPLE}-KEGGmetabolites")

#
# tmp2 = scores_adata.uns['current_pathways_df'][["KEGG", "mz"]]
# tmp2.columns = ["ID", "mz"]
# for pw_name, pw_id in zip(scores_adata.var['pathway_name'].tolist(),
#                           scores_adata.var['pathway_id'].tolist()):
#     print(pw_name)
#     kegg_mtbs_ids = scores_adata.uns['current_pathways_df'].loc[
#         scores_adata.uns['current_pathways_df']['pathway_id'] == pw_id, "KEGG"
#     ].tolist()
#     mtbs_df = adata.var.loc[adata.var['ID'].isin(kegg_mtbs_ids), ["ID", "moleculeNames"]]
#     mtbs_df.reset_index(drop=True, inplace=True)
#     mtbs_df.index.name = "index"
#     haha = tmp2.loc[tmp2["ID"].isin(kegg_mtbs_ids), :]
#     haha = haha.drop_duplicates()
#     haha.index = haha["ID"].tolist()
#     mtbs_df = mtbs_df.join(haha, how="left", on="ID", lsuffix="_x", rsuffix="_y")
#     mtbs_df = mtbs_df.drop(columns=["ID_x", "ID_y"])
#     print(mtbs_df)
#     tmpstr=os.path.join(plot_folder, f'{SAMPLE}-{pw_name}--PwScore.png')
#     sq.pl.spatial_scatter(scores_adata, color=[pw_id], img=False,
#                           title=f"{pw_id} {pw_name}",
#                           save=os.path.abspath(tmpstr)
#                           )
#     #plt.show()
#     #print()


# # pathways

#
# scores_adata = spatial_standard_computations(
#     scores_adata,
#     k=k_dict[SAMPLE][1],
#     as_log=False,   # important !
#     top_variable_features_fraction=0.2)
#
# plot_standard_computs(scores_adata,
#     plot_folder=plot_folder,
#     prefix=f"{SAMPLE}-PathwayScores")



