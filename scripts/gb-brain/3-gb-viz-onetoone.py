import sys
import os
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import zonist.data as zd
import zonist.method as zm

sys.path.append(os.path.join(os.path.dirname(__file__),
                             '../bank4zonist'))

sys.path.append(os.path.join(os.path.dirname(__file__),
                             '../'))


from visualization import custom_scatter_one_feat, plot_with_table
from quantifs_methods import task1, task2
from user_routines import prepare_harmonized_objects

from quantifs_methods import *
from visualization import *

"""
pathways scores 
"""

SAMPLE = "test635_200-1K"   # ""  "test635_80-400"
size_file="sz1"
print(SAMPLE)
mydata_dir = "/home/johanna/spatial_thesis/spatial_metabolomics/data"
#"/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"
subfolder = "gb-brain"
anno_dict = {'HMDB': f"{SAMPLE}_hmdb.csv",
             "KEGG": f"{SAMPLE}_kegg.csv",
             "CoreMetabolome": f"{SAMPLE}_CoreMetabolome.csv"}
plot_folder = f"../../../figures/gb-brain/{SAMPLE}_{size_file}"
min_pw_size=1
max_pw_size=1 # 100
assert os.path.exists(plot_folder), "Cannot continue, the plot_folder does not exist"

me_annots = zd.MetabolomeAnnotations(name=SAMPLE,
                            data_dir=mydata_dir,
                            subfolder=subfolder,
                            annotation_files=anno_dict)
me_annots.load()


adata_mz = ad.read_h5ad(filename=os.path.join(
    mydata_dir, subfolder, f'{SAMPLE}.h5ad'))

# todo: exclude from this analysis the mz that are not good
good_mzs = adata_mz.var.loc[adata_mz.var['mz_initial_quality'] == "good", "mz"].tolist()
adata_mz = adata_mz[:, [str(i) for i in good_mzs]].copy()
print("adata (good mzs only) ",adata_mz)


# print("DANGER!!! get rid of mz that have fdr 0.2 or better ")  # this results in poor list of pws
# ok_fdr_mzs =adata_mz.var.loc[adata_mz.var['fdr'] <= 0.2, "mz"]
# adata_mz = adata_mz[:, [str(x) for x in ok_fdr_mzs]].copy()
# print("adata (further filter!!!! by fdr) ",adata_mz)


pathwaysObj, me_annots, adata_mz, orig_pws_dict = prepare_harmonized_objects(
    mydata_dir, me_annots, adata_mz, prob_out=0.95, min_pw_size=min_pw_size,
max_pw_size=max_pw_size)  #
print("mz post sync", adata_mz)
print("pathways remaining:", len(pathwaysObj.piled['pathway_id'].unique()))
# answer = input("ok?")
# if not answer in ['yes', 'y']:
#     raise AssertionError("killed")
# the KEGG id adata
X_metabolites_probs_df = compute_metabolites_probs_df(adata_mz, me_annots.df)

adata = build_kegg_resolved_adata(adata_mz, X_metabolites_probs_df,
                                  annots_df=me_annots.df)

print("KEGG resolved adata: ", adata)

adata.write_h5ad(os.path.join(mydata_dir, subfolder,
                                         f"{SAMPLE}-KEGG_adata-sz{min_pw_size}.h5ad"))

coords_mz_df = build_coords_df(adata_mz)

##
for pw_id in pathwaysObj.piled['pathway_id'].unique():
    pw_name = pathwaysObj.pathways_ids.loc[pathwaysObj.pathways_ids["pathway_id"]==pw_id, "pathway_name"].item()
    foo = adata.var.copy()
    foo = foo[["mz", "KEGG", "moleculeNames"]]

    df_mz_metabolites__curr_pw = pathwaysObj.piled.loc[
        pathwaysObj.piled['pathway_id'] == pw_id, ["KEGG"]] # df with KEGG id, metabolite name and mz for the current pw_id
    df_mz_metabolites__curr_pw = df_mz_metabolites__curr_pw.drop_duplicates()
    df_mz_metabolites__curr_pw = df_mz_metabolites__curr_pw.join(
        foo, on='KEGG', how='left', lsuffix="_x", rsuffix="_y"
    )
    df_mz_metabolites__curr_pw = df_mz_metabolites__curr_pw.drop(
        columns=["KEGG_x", "KEGG_y"]
    )
    df_mz_metabolites__curr_pw.columns = ["KEGG", "mz", "metabolite name"]
    unique_mzs =  df_mz_metabolites__curr_pw["mz"].unique()
    for i in range(len(unique_mzs)):
        i_mz  = unique_mzs[i]
        sub_df, tmp = task2(df_mz_metabolites__curr_pw, i_mz, coords_mz_df,
                           adata_mz)
        fdr = adata_mz.var.loc[str(i_mz), "fdr"].item()
        title_note = f"FDR={np.around(fdr,4)}"
        plot_with_table(tmp, i_mz, title_note, sub_df,
                        out_file_path=os.path.join(plot_folder,
                             f'{SAMPLE}-{pw_id}-{pw_name}-{i}.png'),
                        spot_size=0.5, dpi=400, marker="s")
        plt.close()
