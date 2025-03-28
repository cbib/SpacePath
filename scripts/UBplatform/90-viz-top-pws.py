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
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from typing import List, Dict

from zonist.aux.matrix_from_metaspace import build_coord_names

sys.path.append(os.path.join(os.path.dirname(__file__),
                             '../bank4zonist'))
from visualization import custom_scatter_one_feat, plot_with_table
from quantifs_methods import mold_data_for_plot__task1, mold_data_for_plot__task2





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


# mz
adata_mz = ad.read_h5ad(filename=os.path.join(
    mydata_dir, subfolder, f'{SAMPLE}.h5ad'))

# metaboli
adata = ad.read_h5ad(os.path.join(mydata_dir, subfolder,
                                         f"{SAMPLE}-KEGG_adata.h5ad"))

coords_mz_df = build_coords_df(adata_mz)

for pw_id in top_spatial_features:
    pw_name = scores_adata.var.loc[scores_adata.var['pathway_id'] == pw_id, "pathway_name"].item()
    tmpstr = os.path.join(plot_folder,
                          f'{SAMPLE}-{pw_id}.png')

    sq.pl.spatial_scatter(scores_adata, color=[pw_id],
                              img=False,
                             title=f"{pw_id} {pw_name}",
                          legend_loc='lower center',
                              save=os.path.abspath(tmpstr))
    df_mz_metabolites__curr_pw = mold_data_for_plot__task1(pw_id,  scores_adata, adata)
    unique_mzs =  df_mz_metabolites__curr_pw["mz"].unique()
    for i in range(len(unique_mzs)):
        i_mz  = unique_mzs[i]
        sub_df, tmp = mold_data_for_plot__task2(df_mz_metabolites__curr_pw, i_mz, coords_mz_df,
                           adata_mz)

        plot_with_table(tmp, i_mz, '', # third positional arg is title_note
                        sub_df, out_file_path=os.path.join(plot_folder,
                             f'{SAMPLE}-{pw_id}-{i}.png'), spot_size=3.8, dpi=300)













