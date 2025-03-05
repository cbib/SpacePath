import sys
import os
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import squidpy.pl
import zonist.data as zd
import zonist.method as zm


sys.path.append(os.path.join(os.path.dirname(__file__),
                             '../bank4zonist'))

sys.path.append(os.path.join(os.path.dirname(__file__),
                             '../'))

from user_routines import prepare_harmonized_objects

from quantifs_methods import *
from visualization import *


"""
pathways scores 
"""


SAMPLE = "test635_80-400"   # "test635_200-1K"

print(SAMPLE)
mydata_dir = "/home/johanna/spatial_thesis/spatial_metabolomics/data"
#"/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"
subfolder = "gb-brain"
anno_dict = {'HMDB': f"{SAMPLE}_hmdb.csv",
             "KEGG": f"{SAMPLE}_kegg.csv"}
plot_folder = "../../../figures/gb-brain"
moran_I_threshold = 0.08
moran_padj_cutoff = 0.05
PROB_OUT=0.95

wrong_mz_manually_defined = [259.022442561, 274.033341593, 267.072155657,
                             377.085612821, 292.076136075, 377.085612821,
                             337.012053957, 346.055808353]

assert os.path.exists(plot_folder), "Cannot continue, the plot_folder does not exist"

me_annots = zd.MetabolomeAnnotations(name=SAMPLE,
                            data_dir=mydata_dir,
                            subfolder=subfolder,
                            annotation_files=anno_dict)
me_annots.load()
me_annots.build_feat_classif_df()

adata_mz = ad.read_h5ad(filename=os.path.join(
    mydata_dir, subfolder, f'{SAMPLE}.h5ad'))
print(adata_mz)


print("using offsample info to define in tissue mzs")
assert me_annots.feat_classif[me_annots.offsample_column].max() <= 1, "bad units"
in_tissue_mz =  me_annots.feat_classif.loc[
    me_annots.feat_classif[me_annots.offsample_column] < PROB_OUT,
    "mz"].astype(str).unique()

adata_mz = adata_mz[:, in_tissue_mz].copy()
print("done")

print("using stats to define spatially relevant mzs")

adata_mz = compute_moran_i(adata_mz, n_neighs=20, n_perms=50) # TODO : make sure of doing at least 50 perms to have padj ok
moran_mz_df = adata_mz.uns['moranI'].copy()
moran_mz_df['mz'] = moran_mz_df.index.astype(str).tolist()

# # moran_mz_df.to_csv("test1-80400.csv", sep='\t', header=True, index=False)
#moran_mz_df = pd.read_csv("test1-80400.csv", sep='\t', index_col=None)

good_mz_spatial = moran_mz_df.loc[
    (moran_mz_df['I'] >= moran_I_threshold) & (
     moran_mz_df['pval_sim_fdr_bh'] <= moran_padj_cutoff),
    "mz"].astype(str).tolist()

good_mz_spatial = set(good_mz_spatial).union(set(wrong_mz_manually_defined))

adata_mz = adata_mz[:, list(good_mz_spatial)].copy()
#  adata_mz.var.to_csv("survived.csv")  # todo del
print("done")

print("final object:", adata_mz)

# # save
answer = input("Filtered AnnData object ready to save. Continue?")
if not answer.lower() in ["y", "yes"]:
    print("Aborted")
else:
    adata_mz.write_h5ad(os.path.join(mydata_dir, subfolder,
                                         f"{SAMPLE}.h5ad"))

    print("Saved new data into:", os.path.join(mydata_dir, subfolder))




##   drafts...
#
# import seaborn as sns
# moran_mz_df['-logpadj'] = -1 * (np.log10(moran_mz_df['pval_sim'].to_numpy() + 0.001))
# sns.scatterplot(data=moran_mz_df, x='I', y='-logpadj',
#                 s=3, alpha=0.5)
# plt.show()


#  for imz in wrong_mz_manually_defined:
#     if str(imz) in adata_mz.var['mz'].astype(str).tolist():
#         print("yes!")
#     else:
#         print(imz, "not found")

#
# import squidpy as sq
# j=0
# for imz in wrong_mz_manually_defined:
#     tmpstr = os.path.join(plot_folder, f"Bad-{j}.png")
#     try:
#         try:
#             imorani = moran_mz_df.loc[moran_mz_df['mz'] == str(imz), "I"].item()
#             padj_i = moran_mz_df.loc[
#                     moran_mz_df['mz'] == str(imz), "pval_sim_fdr_bh"].item()
#             imorani = round(imorani, 2)
#             padj_i = round(imorani, 2)
#         except:
#             imorani = 0
#             padj_i = 0
#         sq.pl.spatial_scatter(adata_mz,
#                               color =str(imz), img=False,
#                               title=f'{np.around(float(imz), 3)} (I={imorani}), padj={padj_i}',
#                               save=os.path.abspath(tmpstr), dpi=200)
#     except Exception as e:
#         print(str(imz), "not found", e)
#
#     j +=1
#
# import sys
# sys.exit()