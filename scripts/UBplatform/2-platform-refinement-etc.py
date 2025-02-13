import sys
import os
import numpy as np
import pandas as pd
import anndata as ad
import squidpy as sq
import networkx as nx
from tqdm import tqdm
from sklearn.metrics.pairwise import cosine_similarity

from scripts.bank4zonist.quantifs_methods import probabilities_updater

sys.path.append(os.path.join(os.path.dirname(__file__),
                             '../bank4zonist'))
from quantifs_methods import *



SAMPLE = "UBcalibr2"
#"UBcalibr2"
print(SAMPLE)
mydata_dir = "/home/johanna/spatial_thesis/spatial_metabolomics/data"
#"/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"
subfolder = "UBplatform"
#"/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"


plot_folder = "../../../figures/platform"



# kegg resolved adata
adata = ad.read_h5ad(os.path.join(mydata_dir, subfolder,
                                         f"{SAMPLE}-KEGG_adata.h5ad"))

# pathways
scores_adata = ad.read_h5ad(os.path.join(mydata_dir, subfolder,
                                         f"{SAMPLE}-PathwayScores.h5ad"))

refinement_mz_pws_dict = ugly_long_function(scores_adata.uns['current_pathways_df'])

X_metabolites_probs_df = pd.DataFrame(data=adata.layers['X_probabilities'],
                                      index=adata.obs['obs'].tolist(),
                                      columns=adata.var['ID'].tolist())


# spatial parsing
RADIUS = 1
sq.gr.spatial_neighbors(adata, radius=RADIUS, coord_type="generic")
connectivities_sparse_mx = adata.obsp["spatial_connectivities"]

coordinates_names_str = adata.obs["obs"].tolist()
coords_matcher_dict = {k: v for k, v in zip(
    range(0, len(coordinates_names_str)), coordinates_names_str)}
coords_matcher_rev = {v: k for k, v in coords_matcher_dict.items()}

print(f"parsing spots : for each spots, neighbors at d = {RADIUS} are evaluated")
COSINE_SIMILARITY_THRESHOLD = 0.6
G = nx.from_scipy_sparse_array(connectivities_sparse_mx)
nodes_list = G.nodes()

new_probabs = dict()

for v1 in tqdm(nodes_list):
    neighs_of_v1 = G.neighbors(v1)
    tmp_mx1 = adata[coords_matcher_dict[v1], :].X.todense()
    spot1 = coords_matcher_dict[v1]
    for v2 in neighs_of_v1:
        tmp_mx2 = adata[coords_matcher_dict[v2], :].X.todense()
        local_cosine_sim = cosine_similarity(np.array(tmp_mx1),
                                             np.array(tmp_mx2)) # 1 x 1 mx
        local_cosine_sim = float(local_cosine_sim.reshape(
            local_cosine_sim.shape[1])[0])
        if local_cosine_sim >= COSINE_SIMILARITY_THRESHOLD:
            spot2 = coords_matcher_dict[v2]
            tmp_probs_dict = probabilities_updater(refinement_mz_pws_dict,
                                       scores_adata.to_df(),
                    X_metabolites_probs_df, spot1, spot2)
            if len(tmp_probs_dict) != 0:
                new_probabs[spot1] = tmp_probs_dict
