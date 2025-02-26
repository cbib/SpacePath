import os
import numpy as np
import pandas as pd
import anndata as ad
import zonist.data as zd
import zonist.method as zm

from scipy.sparse import csr_matrix
import sys
import scanpy as sc
import squidpy as sq


sys.path.append(os.path.join(os.path.dirname(__file__),
                             '../bank4zonist'))  # TODO: fix in the same way ailleurs !!
from quantifs_methods import *
from visualization import *

"""
"""
print("oooo")
SAMPLE = "whole-ellis"  # f1.h5ad saved by -pixel-selection.py
mydata_dir = "/home/johanna/spatial_thesis/spatial_metabolomics/data"
#"/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"

subfolder = "whole-ellis"
anno_dict = {'HMDB': "whole-ellis_hmdb.csv",
             'LipidMaps': "whole-ellis_lipidmaps.csv"}
plot_folder = "../../../figures/ellis"

assert os.path.exists(plot_folder), "Cannot continue, the plot_folder does not exist"

pathwaysObj = zd.MetabolicPathways(
    data_dir=mydata_dir,
    subfolder="metabolites-mappings",
    official_file="KEGG_hsa_pathways_compounds_R110.gmt")

pathwaysObj.load()  # just want to have the pathways names: TODO: better save the computed pathways anndata with the var metadata having names of pathways !

kpca_score_adata = ad.read_h5ad(filename=os.path.join(mydata_dir, subfolder,
                                         f"{SAMPLE}-PathwayScores.h5ad"))
kpca_scores = kpca_score_adata.to_df()
print("computed scores:")
print("min:", kpca_scores.min().min())
print("max:", kpca_scores.max().max())


pws_id_names = pathwaysObj.pathways_ids.loc[
    pathwaysObj.pathways_ids['pathway_id'].isin(kpca_scores.columns.tolist()),:]
pws_id_names.index = pws_id_names['pathway_id']
pws_id_names_dict = pws_id_names["pathway_name"].to_dict()
print(len(pws_id_names_dict))

plot_selected_pathways(kpca_score_adata, pws_id_names_dict,
                           f"{SAMPLE}--", plot_folder)

# kpca_scores = (kpca_scores + abs(kpca_scores.min().min())) + 1 # shift to right
# kpca_score_adata.X = csr_matrix(kpca_scores)
kpca_score_adata = spatial_standard_computations(kpca_score_adata,
                                                 k=12, as_log=False)
plot_standard_computs(kpca_score_adata, plot_folder=plot_folder,
                      prefix=f'{SAMPLE}PathwaysScoresV3')



