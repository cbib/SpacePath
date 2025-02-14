import sys
import os
import numpy as np
import pandas as pd
import anndata as ad
import squidpy as sq
import networkx as nx
from tqdm import tqdm
import zonist.method as zm
import zonist.helpers as zh
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
SF = 0.5   #   stealing_fraction = 0.5  # or 1

sq.gr.spatial_neighbors(adata, radius=RADIUS, coord_type="generic")
connectivities_sparse_mx = adata.obsp["spatial_connectivities"]

coordinates_names_str = adata.obs["obs"].tolist()
coords_matcher_dict = {k: v for k, v in zip(
    range(0, len(coordinates_names_str)), coordinates_names_str)}
coords_matcher_rev = {v: k for k, v in coords_matcher_dict.items()}

print(f"parsing spots : for each spots, neighbors at d = {RADIUS} are evaluated")
COSINE_SIMILARITY_THRESHOLD = 0.9
G = nx.from_scipy_sparse_array(connectivities_sparse_mx)
nodes_list = G.nodes()

new_probabs = dict()

for v1 in tqdm(nodes_list):
    neighs_of_v1 = G.neighbors(v1)

    # to go faster, this is to avoid visiting all neighs
    neighs_of_v1 = [[v for v in neighs_of_v1][0]]   # todo: deactivate
    ## end to go faster
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
                    X_metabolites_probs_df, spot1, spot2,
                    stealing_fraction=SF)
            if len(tmp_probs_dict) != 0:
                new_probabs[spot1] = tmp_probs_dict

# add the new_probabilities to the probabs:
X_probs_new = X_metabolites_probs_df.copy()
for spot in new_probabs.keys():
    for metabolite in new_probabs[spot].keys():
        X_probs_new.loc[spot, metabolite] = new_probabs[spot][metabolite]
print("ok new probs")

adata_mz = ad.read_h5ad(filename=os.path.join(
    mydata_dir, subfolder, f'{SAMPLE}_mx.h5ad'))

adata.var['KEGG'] = adata.var['ID']
adata_new = build_kegg_resolved_adata(adata_mz,
                                      X_probs_new ,
                                  annots_df=adata.var.copy())
print("ok new adata")
adata_new.write_h5ad(os.path.join(mydata_dir, subfolder,
                                         f"{SAMPLE}-crazy_adata.h5ad"))

pathways_dict = zh.pws_piled_to_dict(scores_adata.uns['current_pathways_df'], "pathway_id", mz_or_kegg="KEGG")
new_kpca_scores = zm.ComputePathwaysScores.run(pathways_dict, adata_new)
adata_new.write_h5ad(os.path.join(mydata_dir, subfolder,
                                         f"{SAMPLE}-crazyScores.h5ad"))

