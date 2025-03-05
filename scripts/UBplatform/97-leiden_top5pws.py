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
from sklearn import decomposition
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from typing import List, Dict
from scipy.spatial.distance import cdist
import igraph as ig




SAMPLE =  "UBcalibr2"

print(SAMPLE)
mydata_dir = "/home/johanna/spatial_thesis/spatial_metabolomics/data"
#"/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"
subfolder = "UBplatform"
anno_dict = {'HMDB': f"{SAMPLE}_hmdb.csv",
             "LipidMaps": f"{SAMPLE}_lipidmaps.csv"}
plot_folder = "../../../figures/platform"


assert os.path.exists(plot_folder), "Cannot continue, the plot_folder does not exist"

scores_adata = ad.read_h5ad(os.path.join(mydata_dir, subfolder,
                                         f"{SAMPLE}-PathwayScores.h5ad"))
#

def build_coords_df(adata):
    df = pd.DataFrame(data=adata.obsm['spatial'],
                 columns=['x', 'y'],
                 index=adata.obs.index.tolist())
    return df

def rescale_arr(my_arr, new_min=1, new_max=20):
    relat_arr = my_arr - my_arr.min() / (my_arr.max() - my_arr.min())
    arr_scaled = relat_arr * (new_max - new_min) + new_min

    return arr_scaled

print(scores_adata)
sq.gr.spatial_neighbors(scores_adata, n_neighs=20)
sq.gr.spatial_autocorr(scores_adata, mode="moran", n_perms=100, n_jobs=1, seed=0)
scores_adata.uns['moranI']['pathway_id'] = scores_adata.uns['moranI'].index.tolist()
top_spatial_features = scores_adata.uns['moranI'].head(
    5  # when using 10 was horrible
)["pathway_id"].tolist()
print("subsetting with only spatially relevant features (using Moran'I metrics)")
scores_adata = scores_adata[:, top_spatial_features]
print(top_spatial_features)


# for ipw in top_ok:
#     pw_name = scores_adata.var.loc[scores_adata.var["pathway_id"] == ipw, 'pathway_name'].item()
#     sq.pl.spatial_scatter(scores_adata, img=False, color=[ipw], title=pw_name)
#     plt.show()


x_y_df = build_coords_df(scores_adata)

eucli_faraway = 50
eucli_close = 15


edges = []
weights = []
# slst = adata_mz.obs['obs'].tolist()
# evenly_sampled_spots = [slst[t] for t in range(0, adata_mz.obs.shape[0], 50)]
even_seeds_indexes = [t for t in range(0, scores_adata.obs.shape[0], 400)]

print("\nCreating graph of Profile's similarities (with cosine similiarities)\n")
#euclidean_distances_all = list()
for i_spot_index in range(len(even_seeds_indexes)):
    i_spot_name = scores_adata.obs.index.tolist()[i_spot_index]
    spot1_arr = scores_adata[i_spot_name, :].X.todense()
    coords1 = x_y_df.loc[i_spot_name, ['x', 'y']].tolist()
    for other_spot_index in tqdm(range(len(scores_adata.obs.index.tolist()))):
        other_spot_name = scores_adata.obs.index.tolist()[other_spot_index]
        if other_spot_name != i_spot_name:
            spot2_arr = scores_adata[other_spot_name, :].X.todense()
            local_cosine_sim = cosine_similarity(np.array(spot1_arr),
                                                 np.array(
                                                     spot2_arr))  # 1 x 1 mx.  E.g. : cosine_similarity(np.array([[8,5]]),np.array([[8,6]])) :  Out: array([[0.99639801]])
            local_cosine_sim = float(local_cosine_sim.reshape(
                local_cosine_sim.shape[1])[0])

            coords2 = x_y_df.loc[other_spot_name, ['x', 'y']].tolist()

            local_spatial_dist = cdist(np.array([coords1]), np.array([coords2]), metric='euclidean')
            local_spatial_dist = float(local_spatial_dist.reshape(local_spatial_dist.shape[1])[0])
            # todo: think: use local_spatial_dist later ?
            edges.append((i_spot_index, other_spot_index))

            weights.append(local_cosine_sim)


G_bad_initial_clusters = ig.Graph()  # of igraph
G_bad_initial_clusters.add_vertices(len(scores_adata.obs.index))
G_bad_initial_clusters.add_edges(edges)
G_bad_initial_clusters.es['weight'] = rescale_arr(np.array(weights), new_min=1, new_max=50)

clustering = G_bad_initial_clusters.community_leiden(
    objective_function='modularity',
    weights="weight",
    resolution_parameter=1  # finest
)

print("clusters  :", np.array(clustering.membership).max())
print([str(j) for j in sorted(set(clustering.membership))])
scores_adata.obs['leiden_igraph'] = pd.Categorical(
    [str(x) for x in clustering.membership],
    categories=[str(j) for j in sorted(set(clustering.membership))])


sq.pl.spatial_scatter(scores_adata, color=['leiden_igraph'], img=False)
plt.show()



##########
#   old chunks

# if local_cosine_sim >= 0.9:
#     weights.append(20)
# if (local_cosine_sim >= 0.7) and (local_cosine_sim < 0.8):
#     weights.append(10)
# elif (local_cosine_sim < 0.7) and (local_cosine_sim > 0.1):
#     weights.append(5)
# elif local_cosine_sim <= 0.1:
#     weights.append(1)