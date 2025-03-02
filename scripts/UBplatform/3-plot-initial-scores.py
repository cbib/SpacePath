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
         "UBcalibr2": [4, 3],
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

# mz
adata_mz = adata_mz[: , scores_adata.uns[
    "current_pathways_df"]["mz"].astype(str).unique() ]  # first subset good mz
print(adata_mz)

adata_mz = spatial_standard_computations(
    adata_mz,
    k=k_dict[SAMPLE][0], top_variable_features_fraction=0.2)

plot_standard_computs(adata_mz,
                      plot_folder=plot_folder, prefix=f"{SAMPLE}-MZ")


# metaboli

adata = spatial_standard_computations(
    adata,
    k=k_dict[SAMPLE][0],
    top_variable_features_fraction=0.2)

plot_standard_computs(adata,
     plot_folder=plot_folder, prefix=f"{SAMPLE}-KEGGmetabolites")


# pathways
# note : remember this is bad idea to replace .X by the shifted version, only for test, do not save corrupted object!
# decision: never again do shifting (if so, create shifted_adata).
#  todo: delete :
# scores_adata.X = csr_matrix(scores_adata.X.todense() +
#                            scores_adata.X.todense().min().min() + 1)

scores_adata = spatial_standard_computations(
    scores_adata,
    k=k_dict[SAMPLE][1],
    as_log=False,   # important !:  in this way,no need to "shift"
    top_variable_features_fraction=0.2)

plot_standard_computs(scores_adata,
    plot_folder=plot_folder,
    prefix=f"{SAMPLE}-PathwayScores")
