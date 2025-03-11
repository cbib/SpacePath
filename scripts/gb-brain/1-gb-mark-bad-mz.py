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
clear off not-spatially relevant, or misleading, mz
"""

SAMPLE = sys.argv[1]#"test635_80-400"   # "test635_200-1K"  # todo: pending to sweep the larger mz range file


print(SAMPLE)
mydata_dir = "/home/johanna/spatial_thesis/spatial_metabolomics/data"
#"/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"
subfolder = "gb-brain"
anno_dict = {'HMDB': f"{SAMPLE}_hmdb.csv",
             "KEGG": f"{SAMPLE}_kegg.csv",
             "CoreMetabolome": f"{SAMPLE}_CoreMetabolome.csv"}
moran_I_threshold = 0.02 # this is loose, so define by hand those that are visually not clustered at all, or noisy
moran_padj_cutoff = 0.05

do_step_viz_pick_wrong = False # change to false when completed
# those that have exhibit intense artifacts (stripes, etc) or that have random intensities across pixels
# todo: todo maybe???:    transform this into a txt files: one fore each sample name.txt in data/ ?????
wrong_mz_manually_defined_dict = {
    'test635_80-400': [135.029896921, 138.032729825, 140.048379889, 153.043628857, 163.027978793,
                             167.059278921, 166.038877825, 180.054527889, 202.050966713, 204.054527889,
                             219.978325273, 259.022442561,
                             274.033341593, 267.072155657, 277.059672845, 218.992655905,  292.076136075,
                             299.116793797, 330.125943161 , 341.108935089,
                             377.085612821, # bad border
                              351.119774545, 266.070381999,
                              368.092667239, 374.965451185, 274.033341593,
                             349.037828277, 333.042977395,   # small ball
                             295.110645797, 297.126295861, # too much noise
                              267.081538569
                             ],   ##  note: this is nice clustered but with artifact :  346.055808353, 300.048991657,  kept for now
    'test635_200-1K': [] # nothing yet for this other slice
                             }

print(f"Will mark quality of mz features with parameters:\n"
      f"Moran's I value >= {moran_I_threshold}\n"
      f"Adjusted p-value (BH) of Moran's I <= {moran_padj_cutoff}")
if len(wrong_mz_manually_defined_dict[SAMPLE]) > 0:
    print(f"and a list of {len(wrong_mz_manually_defined_dict[SAMPLE])} "
          f"manually defined misleading mz")

adata_mz = ad.read_h5ad(filename=os.path.join(
    mydata_dir, subfolder, f'{SAMPLE}.h5ad'))
print(adata_mz)

# todo: this is temporal
try:
    adata_mz.var = adata_mz.var.drop(columns="mz_quality") # todo: del
    adata_mz.var = adata_mz.var.drop(columns="mz_initial_quality")  # todo: del
except:
    pass

## ok
print("using stats to define spatially relevant mzs")
if not "moranI" in list(adata_mz.uns.keys()):
    adata_mz = compute_moran_i(adata_mz, n_neighs=20, n_perms=50) # TODO : make sure of doing at least 50 perms to have padj ok
    #save with moran I results (avoid to recompute, time expensive)
    adata_mz.write_h5ad(filename=os.path.join(
        mydata_dir, subfolder, f"{SAMPLE}.h5ad"))

moran_mz_df = adata_mz.uns['moranI'].copy()
moran_mz_df['mz'] = moran_mz_df.index.astype(str).tolist()

moran_mz_df.to_csv(f"{SAMPLE}.csv", sep='\t', header=True, index=False)


# ## now check the fdr computed by metaspace  # todo ??????
# todo: rescue the mz that have good fdr ????
#  annotations to be used for several purposes
me_annots = zd.MetabolomeAnnotations(name=SAMPLE,
                                data_dir=mydata_dir,
                                subfolder=subfolder,
                                annotation_files=anno_dict)
me_annots.load()
#

## 1:  "good" by moran's metrics
good_mz_spatial = moran_mz_df.loc[
    (moran_mz_df['I'] >= moran_I_threshold) & (
     moran_mz_df['pval_sim_fdr_bh'] <= moran_padj_cutoff),
    "mz"].astype(str).unique()
print(f" {len(good_mz_spatial)} msz, are within Moran's related defined thresholds")


## aux: visual picking todo: think if doing appart
if do_step_viz_pick_wrong:  # step 1: manually complete wrong list:
    def plotem(adata_mz, me_annots, good_mz_spatial) -> None:
        import scanpy as sc
        import squidpy as sq
        adata_tmp = adata_mz.copy()

        adata_tmp = adata_tmp[:, good_mz_spatial]
        pathwaysObj, me_annots, adata_tmp, initial_pws_dic = prepare_harmonized_objects(
            mydata_dir, me_annots, adata_tmp, prob_out=0.9,
            min_pw_size=2)
        foo = pathwaysObj.piled["mz"].unique()
        foo = np.array(foo)
        print(
            f"visualizing {len(foo)} the features  (only those mapping to any pathway)"
            "after visualizing, please manually complete the list 'wrong_mz_manually_defined'")
        sc.pp.log1p(adata_tmp)
        for i in foo:
            sq.pl.spatial_scatter(adata_tmp, color=str(i), img=False,
                                  title=str(i))
            plt.show()
        # returns none


    plotem(adata_mz, me_annots, good_mz_spatial)

    answer = input("finished step1, continue?")
    if not answer.lower() in ["y", "yes"]:
        print("Aborted")
    # end if
# endif

## 2 :  good by fdr ??? todo ????
me_annots.do_mz_quality_df()
#
# metrics_interest = moran_mz_df.join(me_annots.mz_quality, .....
# todo: make clear we transfer fdr from annotations quality infos to adataMZ
tmpdf = me_annots.mz_quality.copy()
tmpdf.index = [str(x) for x in tmpdf['mz'].tolist()]
tmpdf = tmpdf.loc[adata_mz.var.index.tolist(),:]
adata_mz.var["fdr"] = tmpdf[me_annots.fdr_column].tolist()
del tmpdf
print(adata_mz.var.head())
## 3: "bad" by manual definition
wrong_mz_manually_defined = wrong_mz_manually_defined_dict[SAMPLE]
wrong_mz_manually_defined = [str(z) for z in wrong_mz_manually_defined]
good_mz_spatial = list(set(good_mz_spatial) - set(wrong_mz_manually_defined))

# adata_mz = adata_mz[:, good_mz_spatial].copy()
# print(f"after custom rejection, remain:   {len(good_mz_spatial)} mzs")
# print("final object:", adata_mz)
adata_mz.var["mz_initial_quality"] = "bad"
print(adata_mz.var.columns)

adata_mz.var.loc[adata_mz.var['mz'].isin(
    [float(i) for i in good_mz_spatial]), "mz_initial_quality"] = "good"
print("quality marked as follows:")
print(adata_mz.var.groupby('mz_initial_quality').size().reset_index(name='counts'))
# # save
answer = input("Continue?")#input("Filtered AnnData object ready to save. Continue?")
if not answer.lower() in ["y", "yes"]:
    print("Aborted")
else:
    adata_mz.write_h5ad(filename=os.path.join(
        mydata_dir, subfolder, f"{SAMPLE}.h5ad"))
    print("Saved new data into:", os.path.join(mydata_dir, subfolder))

######
#  outputs
# =
# test635_80-400
# using stats to define spatially relevant mzs
#  407 msz, are within Moran's related defined thresholds
# after custom rejection, remain:   386 mzs
# final object: AnnData object with n_obs × n_vars = 59444 × 386
#     obs: 'obs', 'sum_out_tissue_quants', 'total_quants', 'ratio_out_tissue__total', 'spot_quality'
#     var: 'mz'
#     uns: 'moranI', 'spatial', 'spatial_neighbors', 'spot_quality_colors'
#     obsm: 'spatial'
#     obsp: 'spatial_connectivities', 'spatial_distances'

# =
#
#####


############
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