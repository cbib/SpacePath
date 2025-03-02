
import string
import sys
import os
import numpy as np
import time
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
from typing import Dict, List

import zonist.data as zd

from scipy.sparse import csr_matrix

import sspa as sspa

sys.path.append(os.path.join(os.path.dirname(__file__),
                             '../bank4zonist'))
from quantifs_methods import *
from visualization import *
from datasetcomplement import (yield_specific_typical_pathways)



SAMPLE = "whole-ellis"  # .h5ad saved by       -pixel-selection.py
mydata_dir =  "/home/johanna/spatial_thesis/spatial_metabolomics/data" #"/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"
subfolder = "whole-ellis"
anno_dict = {'HMDB': "whole-ellis_hmdb.csv",
             'LipidMaps': "whole-ellis_lipidmaps.csv"}
plot_folder = "../../../figures/ellis"




## pathways
pathwaysObj = zd.MetabolicPathways(
    data_dir=mydata_dir,
    subfolder="metabolites-mappings",
    official_file="KEGG_hsa_pathways_compounds_R110.gmt")

pathwaysObj.load()



kpca_score_adata = ad.read_h5ad(os.path.join(mydata_dir, subfolder,
                                         f"{SAMPLE}-PathwayScores.h5ad"))

kpca_new_obj = kpca_score_adata.copy()

# todo: delete these two lines:
#tmpX = (kpca_new_obj.X.todense() + abs(kpca_new_obj.X.min().min())) * 1e3
#kpca_new_obj.X = csr_matrix(tmpX)


#kpca_new_obj = spatial_standard_computations(kpca_new_obj, k=12, as_log=False)

hoho = kpca_new_obj.var.copy().astype(str)
hyhy = pathwaysObj.pathways_ids.copy().astype(str)

yes = pd.merge(hoho, hyhy, how='left',
              left_on="pathway_id", right_on="pathway_id",
             left_index=False, right_index=False)
print(yes.head())


kpca_score_adata = compute_moran_i(kpca_score_adata, 12)

a = kpca_score_adata.uns['moranI'].iloc[:60, :].index.tolist()
a = list()
b = kpca_score_adata.uns['moranI'].iloc[107:, :].index.tolist()
spatially_top_pathways = a + b
for m_id in spatially_top_pathways:
    pw_name = pathwaysObj.pathways_ids.loc[
        pathwaysObj.pathways_ids['pathway_id'] == m_id, "pathway_name"].item()

    tmpstr = os.path.join(plot_folder, f'{SAMPLE}-{pw_name}.png')  # TODO: improve destination folder
    try:
        sq.pl.spatial_scatter(kpca_score_adata, color=[m_id],
                              img=False, title=f"{m_id} {pw_name}",
                              save=os.path.abspath(tmpstr))
    except Exception as e:
        print(e)
        print(f"impossible to plot {m_id} {pw_name}")


print("end!")
print(seresere)
lipid_related_pathways = dict() # do one with terpenoid and sphyngo, then add fatty or others :)
plot_selected_pathways(kpca_score_adata, biologically_chosen_pathways,
                  f"{SAMPLE}--Picked", plot_folder)





