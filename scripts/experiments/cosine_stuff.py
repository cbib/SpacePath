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
from sklearn.metrics.pairwise import cosine_similarity
from sklearn import decomposition
from scipy.sparse import csr_matrix
from sklearn import preprocessing
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from typing import List, Dict
from sklearn.cluster import KMeans

import igraph as ig


#   n obs by group, 10 features each
n = 200
pt1 = np.random.poisson(15, (n, 10))
pt2 = np.random.gamma(3, 2.4, (n, 10))
pt3 = np.random.poisson(0.1,  (n, 10) )

#pt1 = (pt1 + (pt1.max() * 10 )) * 1.7
#pt2 = (pt2 + (pt2.max() * 100))

mat = np.concatenate((pt1, pt2, pt3))
mat = mat
mat = mat * 10
mat = preprocessing.scale(mat, axis=1, with_mean=False, with_std=False, copy=True)

spots_names = [f"spot{u}" for u in range(n*3)]
groups_names = ["g1" for k in range(n)]
groups_names += ["g2" for k in range(n)]
groups_names += ["g3" for k in range(n)]

feats_names = [f"mz{k}" for k in range(10)]

adata = ad.AnnData(X=csr_matrix(mat),
                   var=pd.DataFrame({'var':feats_names},
                                    index=feats_names),
                   obs=pd.DataFrame({'obs': spots_names,
                                     'groups': groups_names},
                                    index=spots_names))
print(adata)


#

sc.pp.pca(adata, n_comps=6)
print(adata)
sc.pl.pca(adata, color="groups")
plt.show()

print("computing kmeans")
kmeans = KMeans(n_clusters=3, n_init=10, random_state=0).fit(
    adata.obsm['X_pca'])
adata.obs[f'kmeans'] = kmeans.labels_.astype(str)
print("plotting ordinary kmeans")
sc.pl.pca(adata, color="kmeans")
plt.show()


# do the graph
print("computing neighbors ")
sc.pp.neighbors(adata,
                n_neighbors=5,
                metric= 'euclidean',#'cosine',
                knn=True,
                method='umap'
                )  # nearest neighbor graph
print("ok computed neighbors")
#sc.tl.umap(adata)
#sc.pl.umap(adata, color="kmeans")
#plt.show()


print("scanpy's leiden")

nn_graph_profiles = adata.obsp["connectivities"]
print("computing Leiden alg...")
leiden_resol = 0.23
sc.tl.leiden(adata,
             resolution=leiden_resol,
             flavor="igraph", n_iterations=10,
             directed=False,
             adjacency= nn_graph_profiles, #nx.to_scipy_sparse_array(Gderived),
             key_added=f"leiden_by_scanpy")
print("done")
print("plot, scanpy strategy")
# sq.pl.spatial_scatter(
#     adata, color=[f"squidpy_domains(leiden{leiden_resol})"],
#     img=False, title="meh")
# plt.show()
sc.pl.pca(adata, color="leiden_by_scanpy")
plt.show()

######
# do the naive method:


def rescale_arr(my_arr, new_min=1, new_max=20):
    relat_arr = my_arr - my_arr.min() / (my_arr.max() - my_arr.min())
    arr_scaled = relat_arr * (new_max - new_min) + new_min

    return arr_scaled


edges = []
weights = []
# slst = adata_mz.obs['obs'].tolist()
even_seeds_indexes = [t for t in range(0, adata.obs.shape[0], 20)]

print("\nCreating graph of Profile's similarities (with cosine similiarities)\n")
#euclidean_distances_all = list()
for i_spot_index in range(len(even_seeds_indexes)):
    i_spot_name = adata.obs.index.tolist()[i_spot_index]
    spot1_arr = adata[i_spot_name, :].X.todense()
    for other_spot_index in tqdm(range(len(adata.obs.index.tolist()))):
        other_spot_name = adata.obs.index.tolist()[other_spot_index]
        if other_spot_name != i_spot_name:
            spot2_arr = adata[other_spot_name, :].X.todense()
            local_cosine_sim = cosine_similarity(np.array(spot1_arr),
                                                 np.array(
                                                     spot2_arr))  # 1 x 1 mx.  E.g. : cosine_similarity(np.array([[8,5]]),np.array([[8,6]])) :  Out: array([[0.99639801]])
            local_cosine_sim = float(local_cosine_sim.reshape(
                local_cosine_sim.shape[1])[0])

            edges.append((i_spot_index, other_spot_index))

            weights.append(local_cosine_sim)


G_bad_initial_clusters = ig.Graph()  # of igraph
G_bad_initial_clusters.add_vertices(len(adata.obs.index))
G_bad_initial_clusters.add_edges(edges)
G_bad_initial_clusters.es['weight'] = rescale_arr(np.array(weights), new_min=1, new_max=50)

clustering = G_bad_initial_clusters.community_leiden(
    objective_function='modularity',
    weights="weight",
    resolution_parameter=0.9
)

print("clusters  :", np.array(clustering.membership).max())
print([str(j) for j in sorted(set(clustering.membership))])
adata.obs['leiden_igraph'] = pd.Categorical(
    [str(x) for x in clustering.membership],
    categories=[str(j) for j in sorted(set(clustering.membership))])


sc.pl.pca(adata, color="leiden_igraph")
plt.show()
