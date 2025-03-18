import sys
import os
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
import seaborn as sns
import squidpy as sq
import zonist.data as zd
import zonist.method as zm

from scripts.bank4zonist.quantifs_methods import tic_normalization

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

SAMPLE = "MALDI_275_T"


print(SAMPLE)
mydata_dir = "/home/johanna/spatial_thesis/spatial_metabolomics/data"
#"/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"
subfolder = "heiland-msi"
anno_dict = {'HMDB': "Metaspace_annotations_943peaks.csv"}
plot_folder = "../../../figures/heiland"


prob_out=0.2  # 0.8 is too bad, Melanie said
msm_cutoff=0.001
moran_I_threshold = 0.02 # this is loose, so define by hand those that are visually not clustered at all, or noisy
moran_padj_cutoff = 0.05
do_moran_plots=False  # True

print(f"Initial quality of mz features with parameters:\n"
      f"MSM >= {msm_cutoff}\n"
      f"offSample <= {prob_out}")

print(f"And spatial_quality of mz features with parameters:\n"
      f"Moran's I value >= {moran_I_threshold}\n"
      f"Adjusted p-value (BH) of Moran's I <= {moran_padj_cutoff}")

do_step_viz_pick_wrong = False # change to false when completed
# those that have exhibit intense artifacts (stripes, etc) or that have random intensities across pixels
# todo: is this list pertinent??? If yes,  transform this into a txt files: one fore each sample name.txt in data/ ?????


#wx = [str(np.around(j, 6)) for j in wrong_mz_manually_defined_dict[SAMPLE]]
#print(", ".join(wx))


### todo : check! C10H11ClN4 - H]⁻   221.0599    0.001   KEGG    Acetamiprid  ==> is well tumor only

#  annotations to be used for several purposes
me_annots = zd.MetabolomeAnnotations(name=SAMPLE,
                                data_dir=mydata_dir,
                                subfolder=subfolder,
                                annotation_files=anno_dict)
me_annots.load()

me_annots.pull_mz_quality_df()

df = me_annots.mz_quality.copy()

fileoutplotstr = os.path.join(plot_folder, f"{SAMPLE}-msm_&_out.pdf")
sns.set_style("whitegrid")
mypalette = ["#580F41", "#662847", "#907359", "khaki"]
sns.scatterplot(df, x="msm", y="rawOffSampleProb",
                hue="fdr",
                palette=sns.color_palette(mypalette, len(mypalette)),

                ).legend(loc="center right", title="FDR")
plt.savefig(fileoutplotstr)
plt.close()

adata_mz = ad.read_h5ad(filename=os.path.join(
    mydata_dir, subfolder, f'{SAMPLE}.h5ad'))
print(adata_mz)

me_annots.mz_quality['quality_msm_outprob'] = "bad"
me_annots.mz_quality.loc[
    (me_annots.mz_quality['msm'] >= msm_cutoff) & (
        me_annots.mz_quality['rawOffSampleProb'] <= prob_out
    ), "quality_msm_outprob"
] = "good"

## moran's metrics
print("using stats to define spatially relevant mzs")

if not "moranI" in list(adata_mz.uns.keys()):
    print("had to ! ")
    adata_mz = compute_moran_i(adata_mz,
                               n_neighs=10, # todo: exclusive tuned for heiland data !
                               n_perms=50) # TODO : make sure of doing at least 50 perms to have padj ok
    #save with moran I results (avoid to recompute, time expensive)
    adata_mz.write_h5ad(filename=os.path.join(
        mydata_dir, subfolder, f"{SAMPLE}.h5ad"))

moran_mz_df = adata_mz.uns['moranI'].copy()
moran_mz_df['mz'] = moran_mz_df.index.astype(str).tolist()

moran_mz_df.to_csv(f"{SAMPLE}-moranAlone.csv", sep='\t', header=True, index=False)

tmpmo = moran_mz_df.copy()
insiders_mz = me_annots.mz_quality.loc[
    me_annots.mz_quality['quality_msm_outprob'] == "good", 'mz'].astype(str).tolist()
insiders_mz = list(set(insiders_mz).intersection(set(tmpmo.index.tolist())))
print(insiders_mz)
tmpmo = tmpmo.loc[insiders_mz,: ]

tmpmo["log_var_sim"] = np.log(tmpmo["var_sim"].to_numpy())
tmpmo['padj class'] = f" > {moran_padj_cutoff}"
tmpmo.loc[tmpmo["pval_sim_fdr_bh"]<=  moran_padj_cutoff, 'padj class'] = f" <= {moran_padj_cutoff}"
sns.set_style("whitegrid")
fig = sns.scatterplot(tmpmo, x="log_var_sim", y="I",
                hue="padj class",
                style="padj class",
                palette=["teal" ,"slategray"],
                markers=["o", "s"],
                s=20
                )
fig.axhline(moran_I_threshold, ls="--", alpha=0.5)
fig.legend(loc="upper left", title="$p_{adj}$ class"
                         )
fileoutplotstr = os.path.join(plot_folder, f"{SAMPLE}-moran_scores.pdf")
plt.savefig(fileoutplotstr)
plt.close()
print(tmpmo.shape)
del tmpmo

if do_moran_plots:
    only_for_subset_plot_moran(SAMPLE,
        adata_mz, me_annots.mz_quality.loc[
            me_annots.mz_quality['quality_msm_outprob'] == "good", "mz"].tolist(),
                               plot_folder)

#
all_quality_df = adata_mz.var.copy()
all_quality_df["mz"] = all_quality_df["mz"].astype(str).tolist()
me_annots.mz_quality['mz'] = me_annots.mz_quality['mz'].astype(str).tolist()


all_quality_df = custom_join_left(
    all_quality_df, me_annots.mz_quality[['mz', 'msm',
    'rawOffSampleProb', 'fdr', 'quality_msm_outprob']], keycol='mz')

all_quality_df = custom_join_left(all_quality_df, moran_mz_df[['mz', 'I', 'pval_sim_fdr_bh']],
                         keycol='mz')

print(all_quality_df.head())
all_quality_df['spatial_quality'] = "bad"
all_quality_df.loc[
    (all_quality_df['I'] >= moran_I_threshold) & (
     all_quality_df['pval_sim_fdr_bh'] <= moran_padj_cutoff),
    "spatial_quality"] = "good"


# if len(wrong_mz_manually_defined_dict[SAMPLE]) > 0:
#     print(f"and a list of {len(wrong_mz_manually_defined_dict[SAMPLE])} "
#           f"manually defined misleading mz")

all_quality_df['final_quality'] = "bad"
all_quality_df.loc[
    (all_quality_df['quality_msm_outprob'] == "good") & (
            all_quality_df['spatial_quality'] == "good"), 'final_quality'
] = "good"

tmp = all_quality_df.loc[all_quality_df['quality_msm_outprob'] == "good",:]
print("good msm and offsample: ", tmp.shape[0])
tmp = tmp.loc[tmp['spatial_quality'] == "good",:]
print("of these, passed spatial quality thresholds: ", tmp.shape[0])
del tmp

adata_mz.var = all_quality_df

print("FINAL quality marked as follows:")
print(adata_mz.var.groupby('final_quality').size().reset_index(name='counts'))
# todo: some typings conflictual !
adata_mz.var['mz'] = adata_mz.var['mz'].astype(float).to_numpy()
adata_mz.X = csr_matrix(adata_mz.X)
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
# good msm and offsample:  54
# of these, passed spatial quality thresholds:  54
# FINAL quality marked as follows:
#   final_quality  counts
# 0          good      54

# =


