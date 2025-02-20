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

sys.path.append(os.path.join(os.path.dirname(__file__),
                             '../bank4zonist'))
from visualization import *

SAMPLE = "UBcalibr2"
#"UBcalibr2"
print(SAMPLE)
mydata_dir = "/home/johanna/spatial_thesis/spatial_metabolomics/data"
#"/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"
subfolder = "UBplatform"
#"/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"


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
    mydata_dir, subfolder, f'{SAMPLE}_mx.h5ad'))

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
    k_dict[SAMPLE][0],
    plot_folder, prefix=f"{SAMPLE}-MZ",
    top_variable_features_fraction=0.2)

# metaboli

adata = spatial_standard_computations(
    adata,
    k_dict[SAMPLE][0],
    plot_folder, prefix=f"{SAMPLE}-KEGGmetabolites",
    top_variable_features_fraction=0.2)


# pathways
scores_adata.X = csr_matrix(scores_adata.X.todense() +
                            scores_adata.X.todense().min().min() + 1)

scores_adata = spatial_standard_computations(
    scores_adata,
    k_dict[SAMPLE][1],
    plot_folder, as_log=False,  # new param
    prefix=f"{SAMPLE}-PathwayScores",
    top_variable_features_fraction=0.2)