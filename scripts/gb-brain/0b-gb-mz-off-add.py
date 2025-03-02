

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
from scipy.sparse import csr_matrix
from sklearn.cluster import KMeans

import zonist.data as zd

sys.path.append(os.path.join(os.path.dirname(__file__),
                             '../bank4zonist'))

from quantifs_methods import *
from visualization import *

"""
Re-writes the rawOffSampleProb of the annotation files
"""
def compare_out_vs_in_groups_of_spots(adata_mz):
    ratios_res = []
    feats_res = []
    obs_group1 = adata_mz.obs.loc[
        adata_mz.obs['location'] == "out", "obs"].tolist()
    obs_group2 = adata_mz.obs.loc[
        adata_mz.obs['location'] == "in", "obs"].tolist()

    for feature in adata_mz.var['mz'].tolist():
        group1_values = adata_mz[obs_group1, str(feature)].X.copy()
        group2_values = adata_mz[obs_group2, str(feature)].X.copy()
        try:
            ratios_res.append(
                (np.array(group1_values.todense()).mean() + 1
                 ) / (np.array(group2_values.todense()).mean()) + 1)
            feats_res.append(feature)
        except:
            ratios_res.append(
                (np.array(group1_values).mean() + 1
                 ) / (np.array(group2_values).mean()) + 1)
            feats_res.append(feature)

    res_df = pd.DataFrame({"ratio_out_in": ratios_res,
                           "log_ratio_out_in": np.log2(np.array(ratios_res)),
                           "mz": feats_res })
    return res_df


def plot_out_vs_in_features(res_df, adata_mz) -> None:
    res_df = res_df.sort_values(by="log_ratio_out_in", ascending=False)
    up_out = res_df.head(5)
    tail_out = res_df.tail(5)
    for mz in up_out['mz'].tolist():
        sq.pl.spatial_scatter(adata_mz, color=[str(mz)], img=False,
                              title="out: yes")
        plt.show()
    for mz in tail_out['mz'].tolist():
        sq.pl.spatial_scatter(adata_mz, color=[str(mz)], img=False,
                              title="out: no")
        plt.show()


def compute_fast_probs(res_df):
    res_df = res_df.sort_values(by="log_ratio_out_in", ascending=False)

    x = res_df['log_ratio_out_in'].to_numpy()
    computed_arr = (x - min(x)) / (max(x) - min(x))
    res_df['rawOffSampleProb'] = computed_arr


    return res_df

## start

SAMPLE = sys.argv[1] # "test635_200-1K"#"test635_80-400" # todo

print(SAMPLE)
mydata_dir = "/home/johanna/spatial_thesis/spatial_metabolomics/data"
#"/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"
subfolder = "gb-brain"
anno_dict = {'HMDB': f"{SAMPLE}_hmdb.csv",
             "KEGG": f"{SAMPLE}_kegg.csv"}
plot_folder = "../../../figures/gb-brain"
kmeans_k = {
    "test635_80-400" : 5,
    "test635_200-1K": 4
}
params_clus_choices = {
    "test635_80-400" : {"in": ['0'], "out":['1', '2', '3', '4']},
    "test635_200-1K" : {"in": ['0', '2'], "out":['1','3']}
}

assert os.path.exists(plot_folder), "Cannot continue, the plot_folder does not exist"

off_sample_is_present = list()
for k_type_annot in anno_dict.keys():
    annot_df = pd.read_csv(os.path.join(
        mydata_dir, subfolder, anno_dict[k_type_annot]
    ), comment='#',  # override any comment (metaspace does)
        sep=",", header=0)
    if annot_df["rawOffSampleProb"].isnull().all():
        off_sample_is_present.append(False)
    else:
        off_sample_is_present.append(True)
    del annot_df

assert np.all(
    np.array(
        off_sample_is_present) == np.repeat(False, len(anno_dict))), (
        "Abort: at least one annotation file has valid OffSample information")

adata_mz = ad.read_h5ad(filename=os.path.join(
    mydata_dir, subfolder, f'{SAMPLE}.h5ad'))
print(adata_mz)


k = kmeans_k[SAMPLE]
adata_mz = spatial_standard_computations(   # ok, random_state=0 in fun def
    adata_mz, k=k,  pca_ncomps=8,
    as_log=True, top_variable_features_fraction=0.7)

adata_mz.X = csr_matrix(adata_mz.X.todense() + 1) #  csr to avoid: ImplicitModificationWarning: X should not be a np.matrix, use np.ndarray instead.

# manually set after viz:
clusters_out = params_clus_choices[SAMPLE]["out"]
clusters_in = params_clus_choices[SAMPLE]["in"]

print('clusters out are:', clusters_out)
print('clusters in are:', clusters_in)
print("check the plot and make sure the choice is correct")
plot_standard_computs(adata_mz,
                      plot_folder=None, prefix='nope')
plt.show()

# using a subset to go faster:
slst = adata_mz.obs['obs'].tolist()
evenly_sampled_spots = [slst[t] for t in range(0, adata_mz.obs.shape[0], 50)]
adata_mz = adata_mz[evenly_sampled_spots,:]


answer = input("Continue?")
if not answer.lower() in ["y", "yes"]:
    print("Aborted")
else:
    adata_mz.obs['location'] = ''
    adata_mz.obs.loc[adata_mz.obs[f'kmeans{k}'].isin(clusters_in), 'location'] = 'in'
    adata_mz.obs.loc[adata_mz.obs[f'kmeans{k}'].isin(clusters_out), 'location'] = 'out'

    res_df = compare_out_vs_in_groups_of_spots(adata_mz)
    res_df = compute_fast_probs(res_df)
    res_df.index = res_df["mz"].tolist()
    # plot_out_vs_in_features(res_df, adata_mz)

    for k_type_annot in anno_dict.keys():
        annot_df = pd.read_csv(os.path.join(
            mydata_dir, subfolder, anno_dict[k_type_annot]
        ), comment='#',  # override any comment (metaspace does)
        sep=",", header=0)

        annot_df = annot_df.drop(columns=["rawOffSampleProb"]) # because has all NaN at the origin
        annot_df['mz'] = np.around(annot_df['mz'].to_numpy(), 9)
        annot_df.index = annot_df["mz"].tolist()
        annot_df = annot_df.join(res_df, on="mz", how="left",
                                 lsuffix='_x', rsuffix='_y')

        # if any NaN in OffSample col, delete the entire row(s)
        annot_df = annot_df.dropna(subset=["rawOffSampleProb"])

        annot_df['mz'] = annot_df['mz_x']
        annot_df = annot_df.drop(columns=["mz_x", "mz_y"])
        annot_df.to_csv(os.path.join(
            mydata_dir, subfolder,
            anno_dict[k_type_annot]),  sep=",", header=True, index=False)
        annot_df = None
        print("saved new annotation :", k_type_annot)
