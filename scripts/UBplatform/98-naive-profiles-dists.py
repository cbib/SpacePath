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
from sklearn.metrics import pairwise_distances_argmin_min
from typing import List, Dict

#from quantifs_methods import *
#from visualization import *


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


# for ipw in scores_adata.var["pathway_id"]:
#     pw_name = scores_adata.var.loc[scores_adata.var["pathway_id"] == ipw, 'pathway_name'].item()
#     sq.pl.spatial_scatter(scores_adata, img=False, color=[ipw], title=pw_name)
#     plt.show()




print(scores_adata)
sq.gr.spatial_neighbors(scores_adata, n_neighs=20)
sq.gr.spatial_autocorr(scores_adata, mode="moran", n_perms=100, n_jobs=1, seed=0)
scores_adata.uns['moranI']['pathway_id'] = scores_adata.uns['moranI'].index.tolist()
top_ok = scores_adata.uns['moranI'].head(10)["pathway_id"].tolist()

for ipw in top_ok:
    pw_name = scores_adata.var.loc[scores_adata.var["pathway_id"] == ipw, 'pathway_name'].item()
    sq.pl.spatial_scatter(scores_adata, img=False, color=[ipw], title=pw_name)
    plt.show()



from sklearn.cluster import KMeans

k = 3
pca = decomposition.PCA(n_components=k)
pc = pca.fit_transform(
    scores_adata.to_df()  # pca only accepts dense representations (never sparse csr_matrix)
    #np.asarray(Xnew)
)

pc_df = pd.DataFrame(data = pc ,
        columns = [f'PC{x+1}' for x in range(0,k)])


kmeans = KMeans(n_clusters=k, n_init=10, random_state=0
                ).fit(
    #scores_adata.X
    scores_adata[:, top_ok].copy().X
    #np.asarray(Xnew)
)
scores_adata.obs[f'kmeans'] = kmeans.labels_.astype(str)
closest_to_centroids, _ = pairwise_distances_argmin_min(kmeans.cluster_centers_,
                                                        #scores_adata.X
                                                        scores_adata[:, top_ok].copy().X
                                                        #np.asarray(Xnew)
                                                        )

pc_df['Cluster'] = kmeans.labels_.astype(str)
for inum_close in closest_to_centroids:
    pc_df.iloc[inum_close, -1] = "CENTERS"
#print(pc_df.head(2))
df_var = pd.DataFrame({'var':pca.explained_variance_ratio_,
             'PC':[f'PC{x+1}' for x in range(0,k)]})
sns.lmplot( x="PC1", y="PC2",
  data=pc_df,
  fit_reg=False,
  hue='Cluster', # color by cluster
  legend=True,
  scatter_kws={"s": 30, "alpha": 0.5}) # specify the point size
plt.xlabel(df_var.loc[df_var['PC']== 'PC1', 'var'].item())
plt.ylabel(df_var.loc[df_var['PC']== 'PC2', 'var'].item())
plt.show()

sq.pl.spatial_scatter(scores_adata, img=False, color="kmeans")
plt.show()

x_y_df = build_coords_df(scores_adata)

eucli_faraway = 50
eucli_close = 15



from scipy.spatial.distance import cdist

# k means cluster refinement, will be slow
G_bad_initial_clusters = nx.Graph()
clusters_labels = scores_adata.obs[f"kmeans"].unique()
wait_list = []
new_clus_info: Dict[str, List[str]] = {}
euclidean_distances_all = list()
for i_cluster in range(len(clusters_labels)):
    print("parsing clus: ", i_cluster)
    i_centroid_index = closest_to_centroids[i_cluster]
    cluster_id_str = clusters_labels[i_cluster]
    centroid_spot_name = scores_adata.obs.index.tolist()[i_centroid_index]
    new_clus_info[cluster_id_str] = [centroid_spot_name]   # initialize the list
    spot1_arr = scores_adata[centroid_spot_name, :].X.todense()
    coords1 = x_y_df.loc[centroid_spot_name, ['x', 'y']].tolist()
    spots_curr_cluster = scores_adata.obs.loc[scores_adata.obs['kmeans'] == cluster_id_str, "obs"].tolist()
    for other_spot_name in tqdm(spots_curr_cluster):
        if other_spot_name != centroid_spot_name:
            spot2_arr = scores_adata[other_spot_name, :].X.todense()
            local_cosine_sim = cosine_similarity(np.array(spot1_arr),
                                                 np.array(
                                                     spot2_arr))  # 1 x 1 mx.  E.g. : cosine_similarity(np.array([[8,5]]),np.array([[8,6]])) :  Out: array([[0.99639801]])

            local_cosine_sim = float(local_cosine_sim.reshape(
                local_cosine_sim.shape[1])[0])

            coords2 = x_y_df.loc[other_spot_name, ['x', 'y']].tolist()

            local_spatial_dist = cdist(np.array([coords1]), np.array([coords2]), metric='euclidean')
            local_spatial_dist = float(local_spatial_dist.reshape(local_spatial_dist.shape[1])[0])
            #print("ùùùù", local_spatial_dist)
            euclidean_distances_all.append(local_spatial_dist)

            if local_cosine_sim >= 0.6:
                new_clus_info[cluster_id_str].append(other_spot_name)     # same cluster stay

            elif (local_cosine_sim < 0.6) and (local_cosine_sim >= 0): # moderate similarity
                if local_spatial_dist <= eucli_close:
                    new_clus_info[cluster_id_str].append(other_spot_name)     # same cluster stay
                elif (local_spatial_dist > eucli_close) and local_spatial_dist < eucli_faraway:
                    new_clus_info[cluster_id_str].append(other_spot_name)  # same cluster stay
                elif local_spatial_dist > eucli_faraway:
                    wait_list.append(other_spot_name)

            elif local_cosine_sim < 0:
                wait_list.append(other_spot_name)

sns.displot(x=euclidean_distances_all, bins=30)
plt.show()

for u in tqdm(wait_list):
    u_arr = scores_adata[u, :].X.todense()
    cosines_dists = []
    for i_cluster in range(len(clusters_labels)):
        i_centroid_index = closest_to_centroids[i_cluster]
        cluster_id_str = clusters_labels[i_cluster]
        centroid_spot_name = scores_adata.obs.index.tolist()[i_centroid_index]
        spot1_arr = scores_adata[centroid_spot_name, :].X.todense()
        local_cosine_sim = cosine_similarity(np.array(spot1_arr),
                                             np.array(
                                                 u_arr))

        local_cosine_sim = float(local_cosine_sim.reshape(
            local_cosine_sim.shape[1])[0])
        cosines_dists.append(local_cosine_sim)
    the_closest = np.argmax(np.array(cosines_dists))
    #oo = cosines_dists[the_closest]
    new_clus_info[str(the_closest)].append(u)


scores_adata.obs["newnew"] = ""
for k in new_clus_info.keys():
    obs_list = new_clus_info[k]
    scores_adata.obs.loc[scores_adata.obs['obs'].isin(obs_list), 'newnew'] = k

sq.pl.spatial_scatter(scores_adata, img=False, color="newnew")
plt.show()


    # pick the centroid (a), obtain the profile of it
    # compare each pair of spots (a, b)  within a cluster
    # if profile cosine sim is at least WW, even if spatially distant, make STRONG edge
    # if profile cosine between WW and ZZ (a bit similar but not too much) and spatially close, make weak edge
    # if profile cosine between WW and ZZ (a bit similar but not too much) and spatially distant, NO edge,
    #   and put the 'b' spot in a wait list
    # if profile cosine below ZZ, NO edge
    #    and put the 'b' spot in a wait list

# use the wait_list and the centroid to assign the wait_list to the most likely cluster
# plot !


def stupido():
    print("trying something new to fail as usual")
    coordinates_names_str = scores_adata.obs["obs"].tolist()
    coords_matcher_dict = {k: v for k, v in zip(
        range(0, len(coordinates_names_str)), coordinates_names_str)}
    coords_matcher_rev = {v: k for k, v in coords_matcher_dict.items()}
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
        # end def stupido


