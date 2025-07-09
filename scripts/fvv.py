import os
import numpy as np
import pandas as pd

from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform

import matplotlib.pyplot as plt
import pickle
import seaborn as sns

from local_helpers.local_viz import customized_clustermap

mydata_dir = "../../data"
my_csv_table_path = "../small-test-parameters.tsv"
plot_parent_folder = "../../figures/"


cutoff_tree = 1.5

SAMPLE = "LP-example0"
subfolder = SAMPLE
pwtype="WIKI"

plot_folder = os.path.join(plot_parent_folder, SAMPLE)


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


#
# cg = sns.clustermap(arr_df,figsize=(12.5, 11),
#                     row_linkage=computed_dendrogram,
#                     col_linkage=computed_dendrogram,
#                     cmap="Blues",
#                     dendrogram_ratio=0.15,
#                     cbar_pos=(0.07, 0.84, .015, .14),
#                     xticklabels=False)
# plt.text(-1.7,0.23, "similarity index", rotation=90)
#
# cg.ax_row_dendrogram.set_visible(False)
# cg.ax_heatmap.set_ylabel("")

dg = customized_clustermap(arr_df, computed_dendrogram)
#plt.savefig("tutut.png", dpi=400)
plt.savefig(
    os.path.join(plot_folder, f"{SAMPLE}-hierarchy-molecules.png"),
    dpi=400)
plt.close()

# did not work:

# # https://stackoverflow.com/questions/67879908/lower-triangle-mask-with-seaborn-clustermap
# values = g.ax_heatmap.collections[0].get_array().reshape(arr_df.shape)
# new_values = np.ma.array(values, mask=mask)
# g.ax_heatmap.collections[0].set_array(new_values)
# g.ax_col_dendrogram.set_visible(False)