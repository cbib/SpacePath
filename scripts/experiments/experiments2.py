import sys
import os
import numpy as np
import pandas as pd
import anndata as ad

from scipy.sparse import csr_matrix
from sklearn.utils import sparsefuncs
import squidpy as sq
from scipy.sparse import csr_matrix
from sklearn.utils import sparsefuncs
import squidpy as sq
import networkx as nx
from tqdm import tqdm
from sklearn.metrics.pairwise import cosine_similarity

sys.path.append(os.path.join(os.path.dirname(__file__),
                             '../bank4zonist'))
from quantifs_methods import *


X_mz_binary_df = compute_mz_binary_df(adata_mz)
X_metabolites_probs_df = compute_metabolites_probs_df(X_mz_binary_df,
                                                      me_annots_OKdf)
# the KEGG id adata


# uf = compounds_counts_by_pathway(pathwaysObj)

#is_it_well_done(X_metabolites_probs_df, adata_mz, adata, me_annots_OKdf)

refinement_mz_pws_dict = ugly_long_function(pathwaysObj.piled.copy(), pathwaysObj.to_dict())

for t in refinement_mz_pws_dict.keys():
    print(t)
    print("   ", refinement_mz_pws_dict[t])
# Note: the refinement_mz_pws_dict has all sense and anti-sense pathways comparisons (pw1 vs pw2) and (pw2 vs pw1)
the_mz_for_refinement_list = [float(i) for i in list(refinement_mz_pws_dict.keys())]

########
# spatial info
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
            lala = complicated_stuff_2(refinement_mz_pws_dict,
                                       pathways_pseudoscores,
                    X_metabolites_probs_df, spot1, spot2)
            if len(lala) != 0 :
                new_probabs[spot1] = lala

#print(new_probabs)



