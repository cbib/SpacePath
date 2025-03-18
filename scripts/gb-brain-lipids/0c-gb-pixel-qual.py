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

sys.path.append(os.path.join(os.path.dirname(__file__),
                             '../bank4zonist'))

from quantifs_methods import *


SAMPLE =  "test635_200-1K"


choices_dict = {
    "test635_200-1K" : 1,  #  d
}

print(SAMPLE)
mydata_dir = "/home/johanna/spatial_thesis/spatial_metabolomics/data"
#"/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"
subfolder = "gb-brain"
anno_dict = {'HMDB': f"{SAMPLE}_hmdb.csv",
             "KEGG": f"{SAMPLE}_kegg.csv",
             "CoreMetabolome": f"{SAMPLE}_CoreMetabolome.csv",
             "LipidMaps": f"{SAMPLE}_lipidmaps.csv"}
plot_folder = "../../../figures/gb-brain"


assert os.path.exists(plot_folder), "Cannot continue, the plot_folder does not exist"

me_annots = zd.MetabolomeAnnotations(name=SAMPLE,
                            data_dir=mydata_dir,
                            subfolder=subfolder,
                            annotation_files=anno_dict)
me_annots.load()


adata_mz = ad.read_h5ad(filename=os.path.join(
    mydata_dir, subfolder, f'{SAMPLE}.h5ad'))
print(adata_mz)


if "spot_quality" not in adata_mz.obs.columns.tolist():
    print(f"computing spot quality across {adata_mz.obs.shape[0]} pixels... ")
    if adata_mz.obs.shape[0] > 30000:
        print("be patient, large dataset in process...")
    me_annots.pull_mz_quality_df()   ## NEW
    print(me_annots.mz_quality.head(2))
    if not me_annots.mz_quality[me_annots.offsample_column].isnull().all():

        print("marking spot quality with params: ", choices_dict[SAMPLE])
        adata_mz = adata_mz_mark_spot_quality(
            adata_mz, me_annots.mz_quality
            )

        adata_mz.obs['spot_quality'] = "good"
        adata_mz.obs.loc[adata_mz.obs['ratio_in_out'] < choices_dict[SAMPLE],
        'spot_quality'] = "bad"

        tmpstr = os.path.join(plot_folder, f'{SAMPLE}-pixBefore.png')
        sq.pl.spatial_scatter(
            adata_mz, color=["spot_quality"], img=False,
            title="Before pixel selection", save=os.path.abspath(tmpstr))
        plt.show()

        adata_mz = adata_mz[adata_mz.obs["spot_quality"] == "good", :]
        tmpstr = os.path.join(plot_folder,
                              f'{SAMPLE}-pixAfter.png')
        print("visualizing good pixels: will over-write new filtered adata")

        answer = input("Continue?")  # small dialog
        if not answer.lower() in ["y", "yes"]:
            print("Aborted")
            raise AssertionError("killed")

        sq.pl.spatial_scatter(
            adata_mz, color=["spot_quality"], img=False,
            title="After pixel selection", save=os.path.abspath(tmpstr))
        plt.show()

        adata_mz.write_h5ad(os.path.join(
                            mydata_dir, subfolder, f'{SAMPLE}.h5ad'))
        print("saved adata")
    else:
        print("UP:Skipped, the out-of-tissue probabilities are not available")
else:
    print("already computed quality, and selected pixels for adata_mz, nothing to do")
