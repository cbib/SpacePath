import os
import numpy as np
import pandas as pd

from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform

import matplotlib.pyplot as plt
import pickle
import seaborn as sns


mydata_dir = "../../data"
my_csv_table_path = "../small-test-parameters.tsv"
plot_parent_folder = "../../figures/"

cutoff_tree = 1.5

SAMPLE = "LP-example0"
subfolder = SAMPLE
pwtype="WIKI"
o_fi = os.path.join(mydata_dir,
                        subfolder, f'{SAMPLE}-{pwtype}-CleanAnnots.tsv')
df = pd.read_csv(o_fi, sep='\t', header=0)
with open(os.path.join(mydata_dir, subfolder,
                       f"{SAMPLE}-mol-similarities-matrix.pkl") ,"rb") as f:
    mol_similarities_result = pickle.load(f)

np.random.seed(42)

distance_matrix = 1 - mol_similarities_result['array']
# squareform saves memory by representing the matrix with the upper right triangle
computed_dendrogram = hierarchy.linkage(squareform(distance_matrix),
                              method='ward',
                              # it is 'ward' by default
                              metric='cosine'
                              # it is cosine by default
                              )
df = df[['moleculeIds', 'moleculeNames']].copy()
df = df.drop_duplicates()
df.index = df['moleculeIds'].tolist()
df = df.T
df = df[mol_similarities_result['the_ids']]
df = df.T

trunc_names = [i[:30] for i in df['moleculeNames'].tolist()]
arr_df = pd.DataFrame(data=mol_similarities_result['array'],
                          columns=trunc_names,
                          index=trunc_names)

# todo: remember that in true code the clusters are computed,  compulsory
#clusters = hierarchy.fcluster(computed_dendrogram,
#                    t=cutoff_tree, criterion='distance')

## dendro = hierarchy.dendrogram(computed_dendrogram)

# ordered_indices = hierarchy.leaves_list(computed_dendrogram)
#
# # To get the ordered labels:
# ordered_labels = [trunc_names[i] for i in ordered_indices]
# arr_df = arr_df.loc[ordered_labels, :]
# arr_df = arr_df.loc[:, ordered_labels]

print(arr_df.values.shape)



cg = sns.clustermap(arr_df,
                    figsize=(12.5, 11),
               row_linkage=computed_dendrogram,
               col_linkage=computed_dendrogram,
               cmap="Blues",
               dendrogram_ratio=0.15,
               cbar_pos=(0.01, 0.82, .015, .15),
xticklabels=False
              )
cg.ax_row_dendrogram.set_visible(False)
cg.ax_heatmap.set_ylabel("")
plt.savefig("tutut.png", dpi=400)


# did not work:

# # https://stackoverflow.com/questions/67879908/lower-triangle-mask-with-seaborn-clustermap
# values = g.ax_heatmap.collections[0].get_array().reshape(arr_df.shape)
# new_values = np.ma.array(values, mask=mask)
# g.ax_heatmap.collections[0].set_array(new_values)
# g.ax_col_dendrogram.set_visible(False)