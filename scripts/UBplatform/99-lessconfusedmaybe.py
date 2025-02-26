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
from sklearn.metrics.pairwise import cosine_similarity


sys.path.append(os.path.join(os.path.dirname(__file__),
                             '../bank4zonist'))

sys.path.append(os.path.join(os.path.dirname(__file__),
                             '../'))

from experiments_pws_scores import compute_pathways_pseudoscores
from user_routines import prepare_harmonized_objects


from quantifs_methods import *
from visualization import *


SAMPLE =  "UBcalibr2"

print(SAMPLE)
mydata_dir = "/home/johanna/spatial_thesis/spatial_metabolomics/data"
#"/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"
subfolder = "UBplatform"
anno_dict = {'HMDB': f"{SAMPLE}_hmdb.csv",
             "LipidMaps": f"{SAMPLE}_lipidmaps.csv"}
plot_folder = "../../../figures/platform"


assert os.path.exists(plot_folder), "Cannot continue, the plot_folder does not exist"

#
# me_annots = zd.MetabolomeAnnotations(name=SAMPLE,
#                             data_dir=mydata_dir,
#                             subfolder=subfolder,
#                             annotation_files=anno_dict)
# me_annots.load()
#
# adata_mz = ad.read_h5ad(filename=os.path.join(
#     mydata_dir, subfolder, f'{SAMPLE}_mx.h5ad'))
# print(adata_mz)
#
# print(np.array(adata_mz.X.todense()).shape)
#
#
#
# pathwaysObj, me_annots, adata_mz, orig_pws_dict = prepare_harmonized_objects(
#     mydata_dir, me_annots, adata_mz)
# print("mz post sync", adata_mz)
#
# # the KEGG id adata
# X_metabolites_probs_df = compute_metabolites_probs_df(adata_mz, me_annots.df)
#
# adata = build_kegg_resolved_adata(adata_mz, X_metabolites_probs_df,
#                                   annots_df=me_annots.df)
# print("KEGG resolved adata: ", adata)

## compute PSEUDOSCORES
#scores_df = compute_pathways_pseudoscores(adata, pathwaysObj, min_size=2, method="average")
# scores_df  = zm.ComputePathwaysScores.run(pathwaysObj.to_dict(), adata)
# scores_adata = score_as_adata(scores_df,
#                               adata,
#                               pathwaysObj.pathways_ids,
#                               pathwaysObj.piled,
#                               orig_pws_dict)

scores_adata = ad.read_h5ad(os.path.join(mydata_dir, subfolder,
                                         f"{SAMPLE}-PathwayScores.h5ad"))


def build_coords_df(adata):
    df = pd.DataFrame(data=adata.obsm['spatial'],
                 columns=['x', 'y'],
                 index=adata.obs.index.tolist())
    return df


coords_df = build_coords_df(scores_adata)


print("regular kmeans of pathways scores")
scores_adata = spatial_standard_computations(scores_adata, k=4,
                              as_log=False,
                                  top_variable_features_fraction=0.95)
plot_standard_computs(scores_adata, k=4,
                      plot_folder=None, prefix=f"{SAMPLE}-Pwscore-kmeans"
                      )
plt.show()
print(scores_adata)


##
leiden_resol = 0.4
print("computing neighbors by profile (scanpy standard)")
sc.pp.neighbors(scores_adata, n_neighbors=20, metric='cosine', knn=True, method='umap')  # nearest neighbor graph
nn_graph_profiles = scores_adata.obsp["connectivities"]

# spatial proximity graph
doit=False
if doit:
    print("spatial edges with tutorial method....slowwwww.....")
    sq.gr.spatial_neighbors(scores_adata,
                            n_neighs=20,
                            coord_type="generic")
    nn_graph_space = scores_adata.obsp["spatial_connectivities"]
    print("computing leiden clustering with tutorial method...")
    joint_graph = 0.7 * nn_graph_profiles + 0.3 * nn_graph_space
    sc.tl.leiden(scores_adata,
                 resolution=leiden_resol,
                 adjacency=joint_graph,
                 key_added=f"squidpy_domains(leiden{leiden_resol})")
    print("visualizing the spatial plot of domains")
    sq.pl.spatial_scatter(
        scores_adata, color=[f"squidpy_domains(leiden{leiden_resol})"],
        img=False, title=f"(leiden{leiden_resol})")
    plt.show()
##


print("trying something new to fail as usual")
coordinates_names_str = scores_adata.obs["obs"].tolist()
coords_matcher_dict = {k: v for k, v in zip(
    range(0, len(coordinates_names_str)), coordinates_names_str)}
coords_matcher_rev = {v: k for k, v in coords_matcher_dict.items()}


from scipy.spatial import distance

def is_spatially_distant(scorelocal, cutoff_faraway):
    if scorelocal <= cutoff_faraway:
        return True
    return False

print("checking spatial info to modulate the profiles edges")
G = nx.from_scipy_sparse_array(nn_graph_profiles)
Gderived = nx.Graph()
edges_list = G.edges()

for edge in edges_list:
    u = edge[0]
    v = edge[1]
    spot_name_1 = coords_matcher_dict[u]
    spot_name_2 = coords_matcher_dict[v]
    coords_spot1 = coords_df.loc[spot_name_1, ['x', 'y']].to_numpy().flatten()
    coords_spot2 = coords_df.loc[spot_name_2, ['x', 'y']].to_numpy().flatten()
    print(coords_spot1.reshape(1, -1))
    scorelocal = cosine_similarity(coords_spot1.reshape(1, -1), coords_spot2.reshape(1, -1))
    print(scorelocal)
    uou = G.get_edge_data(u, v)
    edge_weight = uou['weight']
    # if is_spatially_distant(scorelocal, cutoff_faraway=0.4):  # even if faraway they can still be a functional zone, e.g. symetrical and unconnected
    #     local_neighs_spot1 = [coords_matcher_dict[x] for x in G.neighbors(u)]
    #     local_neighs_spot2 = [coords_matcher_dict[w] for w in G.neighbors(v)]
    #
    #     expressions_neighborhood1 = np.array(scores_adata[local_neighs_spot1,:].X.todense()).sum(axis=0)
    #     expressions_neighborhood2 = np.array(scores_adata[local_neighs_spot2,
    #                                     :].X.todense()).sum(axis=0)
    #     neighborhoods_similarity = cosine_similarity(expressions_neighborhood1.reshape(1, -1),
    #                                                  expressions_neighborhood2.reshape(1, -1))
    #     if neighborhoods_similarity <= 0.5:
    #         edge_weight = edge_weight  # decrease its weight if distant,  and if local neigborhoods are dissimilar
    # else:
    #     edge_weight = edge_weight * 2
    local_neighs_spot1 = [coords_matcher_dict[x] for x in G.neighbors(u)]
    local_neighs_spot2 = [coords_matcher_dict[w] for w in G.neighbors(v)]
    expressions_neighborhood1 = np.array(
        scores_adata[local_neighs_spot1, :].X.todense()).sum(axis=0)
    expressions_neighborhood2 = np.array(scores_adata[local_neighs_spot2,
                                      :].X.todense()).sum(axis=0)
    neighborhoods_similarity = cosine_similarity(expressions_neighborhood1.reshape(1, -1),
                                                    expressions_neighborhood2.reshape(1, -1))
    edge_weight = edge_weight + scorelocal + neighborhoods_similarity
    Gderived.add_edge(u, v, weight=edge_weight)

leiden_resol = 0.1
sc.tl.leiden(scores_adata,
             resolution=leiden_resol,
             adjacency= nn_graph_profiles, #nx.to_scipy_sparse_array(Gderived),
             key_added=f"squidpy_domains(leiden{leiden_resol})")

sq.pl.spatial_scatter(
    scores_adata, color=[f"squidpy_domains(leiden{leiden_resol})"],
    img=False, title="meh")
plt.show()





###https://python.igraph.org/en/stable/api/igraph.GraphBase.html#community_leiden
