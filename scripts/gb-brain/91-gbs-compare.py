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
import zonist.data as zd
import zonist.method as zm
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from typing import List, Dict

from zonist.aux.matrix_from_metaspace import build_coord_names

sys.path.append(os.path.join(os.path.dirname(__file__),
                             '../bank4zonist'))
from visualization import custom_scatter_one_feat, plot_with_table
from quantifs_methods import mold_data_for_plot__task1, mold_data_for_plot__task2, build_coords_df





SAMPLE1 = "test635_80-400"  #""  ""
SAMPLE2 = "test635_200-1K"
size_file = "sz2"
mydata_dir = "/home/johanna/spatial_thesis/spatial_metabolomics/data"
#"/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"
subfolder = "gb-brain"
# anno_dict = {'HMDB': f"{SAMPLE}_hmdb.csv",
#              "KEGG": f"{SAMPLE}_kegg.csv",
#              "CoreMetabolome": f"{SAMPLE}_CoreMetabolome.csv"}
# plot_folder = f"../../../figures/gb-brain/{SAMPLE}_{size_file}"

# assert os.path.exists(plot_folder), "Cannot continue, the plot_folder does not exist"

scores_adata1 = ad.read_h5ad(os.path.join(mydata_dir, subfolder,
                                         f"{SAMPLE1}-PathwayScores-{size_file}.h5ad"))
scores_adata2 = ad.read_h5ad(os.path.join(mydata_dir, subfolder,
                                         f"{SAMPLE2}-PathwayScores-{size_file}.h5ad"))

unique_pathways = set(scores_adata1.var['pathway_id'].tolist()).union(set(
scores_adata2.var['pathway_id'].tolist()
))
print(len(unique_pathways))

commonpws = set(scores_adata1.var['pathway_id'].tolist()).intersection(set(
scores_adata2.var['pathway_id'].tolist()
))
print(len(commonpws))
