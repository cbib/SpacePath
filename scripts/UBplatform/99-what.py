import sys
import os
import numpy as np
import pandas as pd
import anndata as ad
import squidpy as sq
import networkx as nx
from tqdm import tqdm
import matplotlib.pyplot as plt
import zonist.data as zd
import zonist.method as zm
from sklearn.metrics.pairwise import cosine_similarity


sys.path.append(os.path.join(os.path.dirname(__file__),
                             '../bank4zonist'))

sys.path.append(os.path.join(os.path.dirname(__file__),
                             '../'))

from experiments_pws_scores import compute_pathways_pseudoscores
from user_routines import prepare_harmonized_objects

from quantifs_methods import *
from visualization import *


SAMPLE =  "UBcalibr5"

print(SAMPLE)
mydata_dir = "/home/johanna/spatial_thesis/spatial_metabolomics/data"
#"/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"
subfolder = "UBplatform"
#"/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"
anno_dict = {'HMDB': f"{SAMPLE}_hmdb.csv",
             "LipidMaps": f"{SAMPLE}_lipidmaps.csv"}
plot_folder = "../../../figures/platform"


assert os.path.exists(plot_folder), "Cannot continue, the plot_folder does not exist"

me_annots = zd.MetabolomeAnnotations(name=SAMPLE,
                            data_dir=mydata_dir,
                            subfolder=subfolder,
                            annotation_files=anno_dict)
me_annots.load()


adata_mz = ad.read_h5ad(filename=os.path.join(
    mydata_dir, subfolder, f'{SAMPLE}_mx.h5ad'))
print(adata_mz)

print(np.array(adata_mz.X.todense()).shape)


pathwaysObj, me_annots, adata_mz, orig_pws_dict = prepare_harmonized_objects(
    mydata_dir, me_annots, adata_mz)
print("mz post sync", adata_mz)

# the KEGG id adata
X_metabolites_probs_df = compute_metabolites_probs_df(adata_mz, me_annots.df)

adata = build_kegg_resolved_adata(adata_mz, X_metabolites_probs_df,
                                  annots_df=me_annots.df)
print("KEGG resolved adata: ", adata)


## compute PSEUDOSCORES
scores_df = compute_pathways_pseudoscores(adata, pathwaysObj, min_size=2, method="average")
scores_adata = score_as_adata(scores_df,
                              adata,
                              pathwaysObj.pathways_ids,
                              pathwaysObj.piled,
                              orig_pws_dict)
print(scores_adata)

## run pseudoscores clustering "efficiently" hahahahha
RADIUS = 1

sq.gr.spatial_neighbors(scores_adata, radius=RADIUS, coord_type="generic")
connectivities_sparse_mx = scores_adata.obsp["spatial_connectivities"]

coordinates_names_str = scores_adata.obs["obs"].tolist()
coords_matcher_dict = {k: v for k, v in zip(
    range(0, len(coordinates_names_str)), coordinates_names_str)}
coords_matcher_rev = {v: k for k, v in coords_matcher_dict.items()}

print(f"parsing spots : for each spot, neighbors at d = {RADIUS} are evaluated")
COSINE_SIMILARITY_THRESHOLD = 0.5
Cgraph= nx.from_scipy_sparse_array(connectivities_sparse_mx)
nodes_list = Cgraph.nodes()

G = nx.Graph()

for v1 in tqdm(nodes_list):
    neighs_of_v1 = Cgraph.neighbors(v1)
    tmp_mx1 = scores_adata[coords_matcher_dict[v1], :].X.todense()
    for v2 in neighs_of_v1:
        tmp_mx2 = scores_adata[coords_matcher_dict[v2], :].X.todense()
        local_cosine_sim = cosine_similarity(np.array(tmp_mx1),
                                             np.array(tmp_mx2))  # 1 x 1 mx
        local_cosine_sim = float(local_cosine_sim.reshape(
            local_cosine_sim.shape[1])[0])
        if local_cosine_sim >= COSINE_SIMILARITY_THRESHOLD:
            G.add_edge(v1, v2, weight=local_cosine_sim)

labels = nx.community.louvain_communities(
    G, weight='weight', resolution=0.2) #
print(len(labels))  # 500 !!!
###https://python.igraph.org/en/stable/api/igraph.GraphBase.html#community_leiden
