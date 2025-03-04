import sys
import os
import numpy as np
import pandas as pd
import anndata as ad
import squidpy as sq
import networkx as nx
from tqdm import tqdm
import matplotlib.pyplot as plt
import zonist.method as zm

from scripts.bank4zonist.quantifs_methods import spatial_standard_computations

sys.path.append(os.path.join(os.path.dirname(__file__),
                             '../bank4zonist'))
from visualization import *

SAMPLE = "UBcalibr2"
#"UBcalibr2"
print(SAMPLE)
mydata_dir = "/home/johanna/spatial_thesis/spatial_metabolomics/data"
#"/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"
subfolder = "UBplatform"


plot_folder = "../../../figures/platform"
SF = 1  # set 0.5 first, then 1

adata_new = ad.read_h5ad(os.path.join(mydata_dir, subfolder,
                                         f"{SAMPLE}-theft{SF}-adata.h5ad"))


scores_new = ad.read_h5ad(os.path.join(mydata_dir, subfolder,
                                         f"{SAMPLE}-theft{SF}-Scores.h5ad"))

print(len(adata_new.uns['refined_metabolites']))

scores_new = spatial_standard_computations(scores_new, 3, as_log=False)
plot_standard_computs(scores_new, plot_folder=plot_folder, prefix=f'{SAMPLE}-Post{SF}')


# the former adata and scores:
adata = ad.read_h5ad(os.path.join(mydata_dir, subfolder,
                                         f"{SAMPLE}-KEGG_adata.h5ad"))
scores = ad.read_h5ad(os.path.join(mydata_dir, subfolder,
                                         f"{SAMPLE}-PathwayScores.h5ad"))


impacted_pathways = scores.uns['current_pathways_df'].loc[
        scores.uns['current_pathways_df']['KEGG'].isin(
            adata_new.uns['refined_metabolites']), "pathway_id"].unique()


for pw in impacted_pathways:
    df = pd.DataFrame(data=scores.obsm['spatial'],
                      columns=['x', 'y'])
    pw_name = scores.var.loc[scores.var['pathway_id'] == pw, "pathway_name"].item()
    df[f'{pw} (Before)'] = np.array(scores[:, pw].X.todense())
    df[f'{pw} (After)'] = np.array(scores_new[:, pw].X.todense())
    custom_scatter_two_features(df,
                                f'{pw} (Before)',
                                f'{pw} (After)',
                                width=12, height=3.5,
                                spot_size=10,
                                save_path=os.path.join(
                                    plot_folder, f'{SAMPLE}-{pw_name}-Post{SF}.png'))

mz_stuff_df =  scores.uns['current_pathways_df'].loc[
        scores.uns['current_pathways_df']['KEGG'].isin(
            adata_new.uns['refined_metabolites']), :]
mz_stuff_df = mz_stuff_df.drop_duplicates()

mz_stuff_df.to_csv(os.path.join(plot_folder, f"{SAMPLE}-WhoRefined-{SF}.tsv"),
              sep='\t', header=True, index=True)



# plot metabolites: before and after
for mtb in adata_new.uns['refined_metabolites']:

    df = pd.DataFrame(data=adata.obsm['spatial'], # the np array of x, y coords
                      columns=['x', 'y'])

    df[f'{mtb} (Before)'] = np.array(adata[:, mtb].X.todense())
    df[f'{mtb} (After)'] = np.array(adata_new[:, mtb].X.todense())
    custom_scatter_two_features(df,
                                f'{mtb} (Before)',
                                f'{mtb} (After)',
                                width=12, height=3.5,
                                spot_size=10,
                                save_path=os.path.join(
                                    plot_folder, f'{SAMPLE}-{mtb}-Post{SF}.png'))


# fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2 )
# fig.set_figheight(10)
# fig.set_figwidth(19)
# sq.pl.spatial_scatter(adata, color=mtb,
#                       img=False, title=f"{mtb} before", ax=ax1)
#
# sq.pl.spatial_scatter(adata_new, color=mtb,
#                       img=False, title=f"{mtb} AFTER", ax=ax2)
# plt.show()