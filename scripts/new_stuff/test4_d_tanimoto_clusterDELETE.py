import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from scipy.cluster.hierarchy import linkage, dendrogram, fcluster

import pickle
import json


with open("tanimoto_similarities.pkl", "rb") as f:
    dsave = pickle.load(f)
##
arr = dsave['arr']
the_identifiers_that_enter = dsave["the_ids"]

myannots = pd.read_csv("identifs_for_tanimoto.csv", index_col=None)

myannots = myannots.loc[myannots['moleculeIds'].isin(the_identifiers_that_enter),
           ['moleculeIds', 'moleculeNames']]
myannots = myannots.assign(name_short= myannots['moleculeNames'].str.slice(0,10))
myannots.index = myannots['moleculeIds'].tolist()
myannots = myannots.T
myannots = myannots[the_identifiers_that_enter]
myannots = myannots.T

arr_df = pd.DataFrame(data=arr,
                      columns=myannots['name_short'].tolist(),
                      index=myannots['name_short'].tolist())


arr_df

sns.clustermap(arr_df, cmap="Blues")
plt.title(f"the clustering of molecules")
plt.tight_layout()
plt.show()

distance_matrix = 1 - arr_df
# make the hierarchical clustering reproducible
np.random.seed(42)
Z = linkage(distance_matrix, method='ward')

# # Plot the dendrogram
# plt.figure(figsize=(10, 5))
# dendrogram(Z)
# plt.title('Hierarchical Clustering Dendrogram')
# plt.xlabel('Sample index or (cluster size)')
# plt.ylabel('Distance')
# plt.show()

# Determine clusters based on a tree cutoff
cutoff = 15
clusters = fcluster(Z, t=cutoff, criterion='distance')
myannots['molecules_cluster'] = clusters

print("the clusters are:", sorted(myannots['molecules_cluster'].unique().tolist()))


o_dict = {} # reinitialize and save as mmt
for clu in sorted(myannots['molecules_cluster'].unique().tolist()):
    mols_in_clus = myannots.loc[myannots['molecules_cluster'] == clu, "moleculeNames"].unique().tolist()
    print(clu, "=> ",len(mols_in_clus))
    first_molecule_group = myannots.loc[myannots['molecules_cluster'] == clu, "moleculeNames"].head(1).item()
    o_dict[f'tg_{clu}'] = [first_molecule_group] + mols_in_clus


max_length = max(len(v) for v in o_dict.values())

for key in o_dict:
    if len(o_dict[key]) < max_length:
        o_dict[key] += [np.nan] * (max_length - len(o_dict[key]))

o_df = pd.DataFrame(o_dict)
o_df.T.to_csv("tanimoto_groupsNAMESDELETE.mmt", sep='\t', header=False)
#


