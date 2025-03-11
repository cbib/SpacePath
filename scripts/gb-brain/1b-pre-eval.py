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

from user_routines import prepare_harmonized_objects

from quantifs_methods import *
from visualization import *
import scanpy as sc

"""
pathways scores 
"""

SAMPLE = "test635_200-1K"## ""  "test635_80-400"  #

print(SAMPLE)
mydata_dir = "/home/johanna/spatial_thesis/spatial_metabolomics/data"
#"/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"
subfolder = "gb-brain"
anno_dict = {'HMDB': f"{SAMPLE}_hmdb.csv",
             "KEGG": f"{SAMPLE}_kegg.csv",
             "CoreMetabolome": f"{SAMPLE}_CoreMetabolome.csv"}
plot_folder = "../../../figures/gb-brain"


assert os.path.exists(plot_folder), "Cannot continue, the plot_folder does not exist"

me_annots = zd.MetabolomeAnnotations(name=SAMPLE,
                            data_dir=mydata_dir,
                            subfolder=subfolder,
                            annotation_files=anno_dict)
me_annots.load()

adata_mz = ad.read_h5ad(filename=os.path.join(
    mydata_dir, subfolder, f'{SAMPLE}.h5ad'))


# todo: exclude from this analysis the mz that are not good
#print(adata_mz.uns['moranI'].loc['259.022442561', :])
good_mzs = adata_mz.var.loc[adata_mz.var['mz_initial_quality'] == "good", "mz"].tolist()
adata_mz = adata_mz[:, [str(i) for i in good_mzs]].copy()
print("adata (good mzs only) ",adata_mz)

pathwaysObj, me_annots, adata_mz, orig_pws_dict = prepare_harmonized_objects(
    mydata_dir, me_annots, adata_mz, prob_out=0.95)
print("mz post sync", adata_mz)
print("number of pathways:", len(pathwaysObj.piled['pathway_id'].unique()))

#print("plotting log(values)")
#sc.pp.log1p(adata_mz)
refinement_mz_pws_dict = ugly_long_function(pathwaysObj.piled)

import squidpy as sq
j = 0
for imz in refinement_mz_pws_dict.keys():
    tmpstr = os.path.join(plot_folder, f"Candidate-{j}.png")
    sq.pl.spatial_scatter(adata_mz,
                         color =str(imz), img=False,
                          title=f'{float(imz)}', # (I={imorani}), padj={padj_i}',

                          save=os.path.abspath(tmpstr), dpi=200)
    plt.close()
    j += 1

# TODO: viz the featurs that are out of tissue cause I dunno if my computation is correct for these GB slices!!!!


#endfor
# j=0
# for imz in adata_mz.var['mz'].unique():
#     tmpstr = os.path.join(plot_folder, f"SEEALL-{j}.png")
#     sq.pl.spatial_scatter(adata_mz,
#                          color =str(imz), img=False,
#                           title=f'{float(imz)}', # (I={imorani}), padj={padj_i}',
#
#                           save=os.path.abspath(tmpstr), dpi=200)
#     plt.close()
#     j += 1

"""
test635_80-400

adata (good mzs only)  AnnData object with n_obs × n_vars = 59444 × 386
    obs: 'obs', 'sum_out_tissue_quants', 'total_quants', 'ratio_out_tissue__total', 'spot_quality'
    var: 'mz', 'mz_initial_quality'

Excludes m/z with out-of-tissue probability >= 0.95
INFO:zonist.data:[discard_out_of_tissue_mz]: Excluded 8 m/z features that are "out-of-the-tissue"
mz post sync AnnData object with n_obs × n_vars = 59444 × 26

number of pathways: 23
Total mz:  26
kegg_oneormore__one_pathway:  11
one_kegg_multiple_pathways:  10
multiple_kegg_shared:  4
susceptible_of_refinement:  1

====
test635_200-1K  
adata (good mzs only)  AnnData object with n_obs × n_vars = 58615 × 380
Excludes m/z with out-of-tissue probability >= 0.95
INFO:zonist.data:[discard_out_of_tissue_mz]: Excluded 15 m/z features that are "out-of-the-tissue"

mz post sync AnnData object with n_obs × n_vars = 58615 × 33
number of pathways: 26
Total mz:  22
kegg_oneormore__one_pathway:  9
one_kegg_multiple_pathways:  9
multiple_kegg_shared:  4
susceptible_of_refinement:  0

"""
