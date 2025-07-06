"""very important simulation of small portion of pixels with rather biological
spatial distribution (as a small filiform structure) that demonstrates that
such structures have bad MoranI, unfortunately ( I ~ 0.1) so if we use moran > 0.2
 we risk to exclude them, and for Thomas they can be relevant because
  invasive cells travelling the corpus callosum can be of that shape (small filiform structure)"""
import numpy as np
import pandas as pd
import anndata as ad
import squidpy as sq
from scipy.stats import gamma, uniform
import matplotlib.pyplot as plt
import zonist.processing.measures_methods as zmm

np.random.seed(42)

thefeat = 13
puits = 2500
fake_oho = np.zeros(shape=(puits, thefeat))
print(fake_oho)
coordinates = []
coords_str = []
for u in range(int(np.sqrt(puits))):
    for v in range(int(np.sqrt(puits))):
        coordinates.append([u,v])
        coords_str.append(str(u) + "_" + str(v))

print(len(coordinates))

def simulate_gamma_data(mean, size):
    shape = mean**2 / 1.0  # Assuming scale parameter is 1.0
    return gamma.rvs(a=shape, scale=1.0, size=size)

means = np.linspace(start=1, stop=10, num=puits)
for i in range(0,puits):
    foo = simulate_gamma_data(mean=means[i], size=thefeat)
    fake_oho[i, :] = foo
# means = np.linspace(start=1, stop=10, num=thefeat)
# for i in range(thefeat):
#     foo = simulate_gamma_data(mean=means[i], size=puits)
#     fake_oho[:, i] = foo

fake_oho[:, 0] = 0

df1 = pd.DataFrame(data=fake_oho, index=coords_str, columns=[str(p) for p in range(thefeat)])

chosen_coords = [[30,21], [31,22], [32, 23], [33, 23], [34,24], [34,25], [34,26] ]

# chosen_coords = ['5_6', '7,8']
for j in chosen_coords:
    k = str(j[0]) + '_' + str(j[1])
    df1.loc[k , '0' ] = 2550

chosen_random_features = ['8', '8', '9', '10', '11', '12']
for t in chosen_random_features:
    df1.loc[:, t] = uniform.rvs(loc=100, scale=1, size=puits)

print(df1.shape)
print(len(coords_str))
print(len([str(p) for p in range(thefeat)]))

adata = ad.AnnData(X=df1.values,
                   obs= pd.DataFrame({'obs':coords_str }),
                   var= pd.DataFrame({'var': [str(p) for p in range(thefeat)]}),
                   obsm={'spatial': np.array(coordinates)})
adata.uns['spatial'] = {'my_slice_name': {
        'scalefactors': {'fiducial_diameter_fullres': 100,
                         'spot_diameter_fullres': 1,
                         'tissue_hires_scalef': 0.5,
                         'tissue_lowres_scalef': 0.5}}}
import scanpy as sc
sc.pp.log1p(adata)
sq.pl.spatial_scatter(adata, color=['0'], img=False)
plt.show()

adata = zmm.compute_moran_i(adata)

df2 = adata.uns['moranI']
df2['feature'] = df2.index.tolist()
df2 = df2.sort_values(by='I', ascending=False)
print(df2)  # the isolated pixels that of feature 0 are not top by moran values, unfortunately
tops = df2.head(6)
for to in tops['feature'].tolist():
    sq.pl.spatial_scatter(adata, color=[str(to)], img=False)
    plt.show()



################### other tests lower relevance

import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from scipy.stats import gamma

# Set random seed for reproducibility
np.random.seed(42)

# Parameters
n_observations = 400
n_variables = 600
n_clusters = 4

# Simulate the dataset
def simulate_gamma_cluster(mean, size):
    shape = mean**2 / 1.0  # Assuming scale parameter is 1.0
    return gamma.rvs(a=shape, scale=1.0, size=size)

# Means for the gamma distributions for each cluster
cluster_means = [10, 5, 1, 20]

# Simulate the dataset
data = np.zeros((n_observations, n_variables))
labels = np.zeros(n_observations, dtype=int)

for i in range(n_clusters):
    start_idx = i * (n_observations // n_clusters)
    end_idx = (i + 1) * (n_observations // n_clusters)
    labels[start_idx:end_idx] = i
    for j in range(n_variables):
        data[start_idx:end_idx, j] = simulate_gamma_cluster(cluster_means[i], n_observations // n_clusters)

# K-means clustering directly on the data
kmeans = KMeans(n_clusters=n_clusters, random_state=42)
kmeans_labels = kmeans.fit_predict(data)
silhouette_kmeans = silhouette_score(data, kmeans_labels)

# PCA on the dataset and then K-means on the first 3 principal components
pca = PCA(n_components=3)
principal_components = pca.fit_transform(data)

plt.figure(figsize=(10, 8))
scatter = plt.scatter(principal_components[:, 0], principal_components[:, 1], c=labels, cmap='viridis')
plt.title('Principal Components Analysis')
plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.colorbar(scatter, label='Cluster')
plt.show()


kmeans_pca = KMeans(n_clusters=n_clusters, random_state=42)
kmeans_pca_labels = kmeans_pca.fit_predict(principal_components)
silhouette_kmeans_pca = silhouette_score(principal_components, kmeans_pca_labels)

# Plot the results
fig, axs = plt.subplots(1, 2, figsize=(14, 6))

# Plot for K-means directly on the data
axs[0].scatter(data[:, 0], data[:, 1], c=kmeans_labels, cmap='viridis')
axs[0].set_title(f'K-means Clustering\nSilhouette Score: {silhouette_kmeans:.2f}')
axs[0].set_xlabel('Variable 1')
axs[0].set_ylabel('Variable 2')

# Plot for K-means on PCA components
axs[1].scatter(data[:, 0], data[:, 1], c=kmeans_pca_labels, cmap='viridis')
axs[1].set_title(f'K-means on PCA Components\nSilhouette Score: {silhouette_kmeans_pca:.2f}')
axs[1].set_xlabel('Principal Component 1')
axs[1].set_ylabel('Principal Component 2')

plt.show()

# Print silhouette scores
print(f'Silhouette Score for K-means directly on the data: {silhouette_kmeans:.2f}')
print(f'Silhouette Score for K-means on PCA components: {silhouette_kmeans_pca:.2f}')
