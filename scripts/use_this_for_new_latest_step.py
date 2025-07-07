import sys
import os
import numpy as np
import pandas as pd
import anndata as ad
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import squidpy as sq

import zonist.processing.measures_methods as zmm

from use_full_steps_minimal_data import wrapper_clustering

mydata_dir = "../../data"
my_csv_table_path = "../small-test-parameters.tsv"
plot_parent_folder = "../../figures/"

SAMPLE = "LP-example0"
subfolder = SAMPLE

k = 3
innerk = 5

adata_mz = ad.read_h5ad(filename=os.path.join(
            mydata_dir, subfolder, f'{SAMPLE}.h5ad'))

print(
            "MANDATORY: take away bad quality m/z")  # todo: how to make an alert ?
adata_mz = adata_mz[:,
                   adata_mz.var.loc[adata_mz.var['final_quality'] == "good",
                   "mz"].astype(str).tolist()]

# ##  repeated chunk of stuff below
# chosen_type_pathway = "KEGG" # or WIKI if desired
#
# o_fi = os.path.join(mydata_dir,
#                     subfolder, f'{SAMPLE}-{chosen_type_pathway}-CleanAnnots.tsv')
# anns = pd.read_csv(o_fi, sep='\t', header=0)
#
# scores_pathways = ad.read_h5ad(filename=os.path.join(mydata_dir, subfolder,
#                                                      f'{SAMPLE}-scores-{chosen_type_pathway}.h5ad'))
# chever1 = scores_pathways.var['pathway_name'].unique()
# print(scores_pathways.var['pathway_name'].unique())

##  stuff below

chosen_type_pathway = "WIKI" # or WIKI if desired

o_fi = os.path.join(mydata_dir,
                    subfolder, f'{SAMPLE}-{chosen_type_pathway}-CleanAnnots.tsv')
anns = pd.read_csv(o_fi, sep='\t', header=0)

scores_pathways = ad.read_h5ad(filename=os.path.join(mydata_dir, subfolder,
                                                     f'{SAMPLE}-scores-{chosen_type_pathway}.h5ad'))
chever2 = scores_pathways.var['pathway_name'].unique()
print(scores_pathways.var['pathway_name'].unique())
# end stuff below


scores_families = ad.read_h5ad(filename=os.path.join(mydata_dir, subfolder,
                                                     f'{SAMPLE}-scores-FAMILIES.h5ad'))



""" todo
do sub clustering with mz == : those having moran more than 0.1

on the subclusters done with the mz , project the pathways scores
find the pathways that are different between the subclusters obtained with mz ?
plot these pathways scores :) 
"""



scores_families.var['set_type'] = "families"
scores_pathways.var['set_type'] = "pathways"

scores_full  = ad.concat((scores_pathways, scores_families),
                         join='outer', axis=1)
scores_full.obsm['spatial'] = adata_mz.obsm['spatial']
scores_full.uns['spatial'] = adata_mz.uns['spatial']
# todo : open kegg and families and run:
## for selected clusters, do something with the scores: "nested clster" within the selected cluster
my_selected_clusters =  [1]
i_clus = 1
clusts_computed = pd.read_csv(os.path.join(mydata_dir, subfolder,
                             f'{SAMPLE}-clusters-k-{k}.txt'), header=0)
the_cluster_spots = clusts_computed.loc[clusts_computed[f'kmeans{k}'] == i_clus, "obs"].tolist()
i_data_mz = adata_mz[the_cluster_spots ,:].copy()




i_data_mz = wrapper_clustering(i_data_mz, k=innerk, moran_cutoff_features=0)
sq.pl.spatial_scatter(
    i_data_mz, color=[f'kmeans{innerk}'], img=False,
    title='the mz , subcluster',
    frameon=False,
    # save=os.path.abspath(plot_out_path), dpi=600
)
# plt.close()
plt.show()

##   use scores:
i_scores = scores_full[the_cluster_spots ,:].copy()
#  sc.pp.log1p()    sc.tl.pca()   # no log  nor pca to be done as scores are already PCs
i_scores = zmm.compute_moran_i(i_scores, n_neighs=5, n_perms=100)
i_scores.uns['moranI']['id'] = i_scores.uns['moranI'].index.tolist()
print("max morani : " , i_scores.uns['moranI']['I'].max())
top_morani_vars = i_scores.uns['moranI'].loc[
    i_scores.uns['moranI']['I'] >= 0.1, 'id'
].tolist()
print(top_morani_vars)
i_scores = i_scores[:, top_morani_vars].copy()

# kmeans = KMeans(n_clusters=innerk, random_state=0).fit(
#     i_scores.X)
# the_labels_str_l = kmeans.labels_.astype(str)
# i_scores.obs[f'kmeans{innerk}'] = the_labels_str_l
#
# sq.pl.spatial_scatter(
#     i_scores, color=[f'kmeans{innerk}'], img=False,
#     title='the scores',
#     # save=os.path.abspath(plot_out_path), dpi=600
# )
# plt.close()
plt.show()

from local_helpers.local_viz import custom_scatter_one_feat
#pw = 'hsa00230'
toto = pd.concat([scores_pathways.uns['current_pathways_df'], scores_families.uns['current_pathways_df']])
print(toto.head())
toto = toto[['pathway_id', 'mz']]
toto = toto.drop_duplicates()
tata = toto['pathway_id'].value_counts().rename_axis(
    'pathway_id').reset_index(name='counts')
booh0 = tata.loc[tata['counts'] >= 3 , 'pathway_id'].tolist()

kop = list(set(booh0).intersection(set(top_morani_vars)))
inorder = i_scores.uns['moranI'].loc[kop,:].copy()
inorder = inorder.sort_values(by='I', ascending=False)

iexcell = i_scores[:, inorder['id'].tolist() ].copy()
superpathways = iexcell.var['pathway_id'].tolist()
print("my super cool:")
print(superpathways[:20])



for pw in superpathways[:2]:
    print("")
    if pw.startswith("WP") or pw.startswith("hsa"):
        ufo = scores_pathways.uns['current_pathways_df'].copy()
        themz_in_pw = ufo.loc[ufo['pathway_id'] == pw, 'mz'].unique().tolist()
        pw_name = scores_pathways.var.loc[scores_pathways.var['pathway_id'] == pw, "pathway_name" ].item()
    else:
        ufo = scores_families.uns['current_pathways_df'].copy()
        themz_in_pw = ufo.loc[ufo['pathway_id'] == pw, 'mz'].unique().tolist()
        pw_name = scores_families.var.loc[
            scores_families.var['pathway_id'] == pw, "pathway_name"].item()
    ##
    print(pw_name)
    #tmp, i_mz,
    countermz=0
    for i_mz in themz_in_pw:

        title = f"\n {pw_name[:35]}  \n  (set ID: {pw})  \n\n m/z = {i_mz} \n"

        i_anno_df = anns.loc[anns['mz'] == i_mz, ['moleculeIds', 'moleculeNames' ]]
        i_anno_df['moleculeNames'] = i_anno_df['moleculeNames'].str.slice(0,25)
        sub_df = i_anno_df.loc[i_anno_df['moleculeIds'].isin(
            ufo.loc[(ufo['mz'] == i_mz) & (ufo['pathway_id'] == pw), 'moleculeIds'].tolist())  ,:]
        print(sub_df['moleculeNames'].tolist())
        if i_anno_df.shape[0] > 15:
            i_anno_df = i_anno_df.head(15)  # limit the lenght of the total annots of the current mz
            i_anno_df.loc[15] = ["...", "...(truncated)..."]
        tmp = adata_mz[:,str(i_mz)].to_df()

        tmp['x'] = np.array(adata_mz.obsm['spatial'])[:,0]
        tmp['y'] = np.array(adata_mz.obsm['spatial'])[:,1]
        spot_size=10
        marker="s"
        fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3)
        fig.set_figheight(4.5)
        fig.set_figwidth(13)

        ax1 = custom_scatter_one_feat(tmp, feature=str(i_mz),
                                      title=title, the_ax=ax1,
                                      spot_size=spot_size,
                                      marker=marker, rasterized=True)
        ax2.set_axis_off()
        start_y = 0.09

        auto_height = sub_df.shape[0] * 0.12
        if sub_df.shape[0] == 1:
            auto_height = 0.17
        if auto_height >= 1.2:
            auto_height = 1.2
        auto_width = 0.9

        ax2.text(auto_width, 0.43, "Annotations mapped to set",
                 rotation=270, va='center')

        table = ax2.table(  # https://www.statology.org/seaborn-table/
            cellText=sub_df.values, rowLabels=None,
            colLabels=sub_df.columns,
            bbox=(-0.08, start_y, auto_width, auto_height),
            # [xmin, ymin, width, height]
            colWidths=(0.43, 0.73), cellLoc='center'
        )
        table.auto_set_font_size(False)
        table.set_fontsize(10.5)
        #table.scale(2, 2)  # may help
        ##
        ax3.set_axis_off()
        ax3.text(auto_width, 0.5, "Total annotations of current m/z",
                 rotation=270, va='center')
        auto_height = i_anno_df.shape[0] * 0.12
        if i_anno_df.shape[0] == 1:
            auto_height = 0.17
        if auto_height >= 1.2:
            auto_height = 1.2

        table = ax3.table(  # https://www.statology.org/seaborn-table/
            cellText=i_anno_df.values, rowLabels=None,
            colLabels=i_anno_df.columns,
            bbox=(-0.08, start_y, auto_width, auto_height),
            colWidths=(0.41, 0.75), cellLoc='center'
        )
        table.auto_set_font_size(False)
        table.set_fontsize(9.5)
        ##
        plt.tight_layout()
        #
        plt.savefig(f"{pw}--{countermz}.pdf", dpi=400)
        plt.close()
        countermz +=1

# ## scores families alone zoom
# famizoom = scores_families[the_cluster_spots, :].copy()
# famizoom = zmm.compute_moran_i(famizoom)
# famizoom.var['I'] = famizoom.uns['moranI']['I'].tolist()
# famizoom = wrapper_clustering(famizoom, k=innerk,
#                                      moran_cutoff_features=0, col_features="pathway_id" )
#
# sq.pl.spatial_scatter(famizoom, color=[f'kmeans{innerk}'], img=False, title="famizoom")
# plt.show()
# ##
