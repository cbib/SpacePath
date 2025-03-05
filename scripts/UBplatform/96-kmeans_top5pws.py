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
from sklearn.cluster import KMeans
from sklearn.metrics import pairwise_distances_argmin_min
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

print(scores_adata)
sq.gr.spatial_neighbors(scores_adata, n_neighs=20)  # needed for Moran's I
sq.gr.spatial_autocorr(scores_adata, mode="moran", n_perms=100, n_jobs=1, seed=0)
scores_adata.uns['moranI']['pathway_id'] = scores_adata.uns['moranI'].index.tolist()
top_spatial_features = scores_adata.uns['moranI'].head(
    5  # when using 10 was horrible
)["pathway_id"].tolist()

print("subsetting with only spatially relevant features (using Moran'I metrics)")
scores_adata = scores_adata[:, top_spatial_features].copy()
print(top_spatial_features)


# for ipw in top_ok:
#     pw_name = scores_adata.var.loc[scores_adata.var["pathway_id"] == ipw, 'pathway_name'].item()
#     sq.pl.spatial_scatter(scores_adata, img=False, color=[ipw], title=pw_name)
#     plt.show()
#
# from sklearn.cluster import KMeans
#
k = 3
# pca = decomposition.PCA(n_components=k)
# pc = pca.fit_transform(
#     scores_adata.to_df()  # pca only accepts dense representations (never sparse csr_matrix)
# )
#
# pc_df = pd.DataFrame(data = pc ,
#         columns = [f'PC{x+1}' for x in range(0,k)])
#
#
kmeans = KMeans(n_clusters=k, n_init=10, random_state=0
                ).fit(
    scores_adata.X  # subset already done
)

scores_adata.obs[f'kmeans'] = kmeans.labels_.astype(str)
print("ok kmeans")
# ## plot dims
# pc_df['Cluster'] = kmeans.labels_.astype(str)
# #print(pc_df.head(2))
# df_var = pd.DataFrame({'var':pca.explained_variance_ratio_,
#              'PC':[f'PC{x+1}' for x in range(0,k)]})
# sns.lmplot( x="PC1", y="PC2",
#   data=pc_df,
#   fit_reg=False,
#   hue='Cluster', # color by cluster
#   legend=True,
#   scatter_kws={"s": 30, "alpha": 0.5}) # specify the point size
# plt.xlabel(df_var.loc[df_var['PC']== 'PC1', 'var'].item())
# plt.ylabel(df_var.loc[df_var['PC']== 'PC2', 'var'].item())
# plt.show()
#
## plot tissue kmeans
sq.pl.spatial_scatter(scores_adata, img=False, color="kmeans")
plt.show()
import sys
sys.exit()
## end kmeans

## leiden as offered by scanpy
# neighbors by profiles:
print("\n ** Checking SCANPY methods **")
print(scores_adata)

print("computing Profile's graph with scanpy method")
sc.pp.neighbors(scores_adata, n_neighbors=1000, # authors say shoud be in range 2 to 100, larger values more global communities
                #knn=True, method="umap",
                knn=True, method="umap", #slower but XXXXXX
                metric="jaccard" # used euclidean default, also minkowski, cosine: all not good
                # jaccard gives relatively good clusters
)
# sc.tl.umap(scores_adata)
# sc.pl.umap(scores_adata, color="kmeans", size=5)
my_flavor="igraph"
print(f"computing Leiden on the Profile's graph (flavor='{my_flavor}')")
sc.tl.leiden(scores_adata,
             resolution=1.2,
             flavor=my_flavor,
             key_added="leiden-scanpy",
             random_state=0,
             use_weights=True,
             n_iterations=2, directed=False)
# sc.pl.umap(scores_adata, color=["leiden-scanpy"], )
plt.show()

sq.pl.spatial_scatter(scores_adata, img=False, color="leiden-scanpy")
plt.show()