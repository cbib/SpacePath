import sys
import os
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import squidpy as sq
import networkx as nx
from docutils.nodes import footer
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
Metaspace, unfortunately did not compute rawOffSampleProb for this dataset.
This script:
  Re-writes the rawOffSampleProb of the annotation files (custom approach) 
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

SAMPLE = "test635_200-1K" #sys.argv[1] # "test635_200-1K"#

print(SAMPLE)
mydata_dir = "/home/johanna/spatial_thesis/spatial_metabolomics/data"
#"/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"
subfolder = "gb-brain"
anno_dict = {'HMDB': f"{SAMPLE}_hmdb.csv",
             "KEGG": f"{SAMPLE}_kegg.csv",
             "CoreMetabolome": f"{SAMPLE}_CoreMetabolome.csv"}
plot_folder = "../../../figures/gb-brain"



me_annots = zd.MetabolomeAnnotations(name=SAMPLE,
                                data_dir=mydata_dir,
                                subfolder=subfolder,
                                annotation_files=anno_dict)
me_annots.load()
me_annots.do_mz_quality_df()

foo = me_annots.preloaded.copy()
#fee = me_annots.mz_quality[['mz', 'fdr']]
#foo = foo.join(fee, how='left', on='mz', lsuffix='_x', rsuffix='_y')
foo = foo.sort_values(by=['mz','fdr'])
foo = foo.drop_duplicates(subset=['mz', 'database','fdr'], keep='first' )


bestof = me_annots.mz_quality.loc[me_annots.mz_quality["fdr"] <= 0.2,"mz"].tolist()
print('the_best:', len(bestof))
#foo = me_annots.unfiltered.copy()
#foo = foo.loc[foo["mz"].isin(bestof),:]

print(UIUIUIIUIU)
oo = {}
for k_type_annot in anno_dict.keys():
    oo[k_type_annot] = pd.read_csv(os.path.join(
        mydata_dir, subfolder, anno_dict[k_type_annot]
    ), comment='#',  # override any comment (metaspace does)
        sep=",", header=0)
    print(type(oo[k_type_annot]['mz'].to_numpy()[22]))
    oo[k_type_annot]['mzstr']  = oo[k_type_annot]['mz'].astype(str).tolist()
    #if annot_df["rawOffSampleProb"].isnull().all():


print(len(oo['CoreMetabolome']['mzstr'].tolist()))
my_done_set = set(oo['HMDB']['mzstr'].tolist()).union(set(oo['KEGG']['mzstr'].tolist()))
print(len(my_done_set.intersection(set(
    oo['CoreMetabolome']['mzstr'].tolist())
)))
