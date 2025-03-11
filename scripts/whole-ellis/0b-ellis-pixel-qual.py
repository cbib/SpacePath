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


SAMPLE =  "whole-ellis"  # done, success

print(SAMPLE)
mydata_dir = "/home/johanna/spatial_thesis/spatial_metabolomics/data"
#"/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"
subfolder = "whole-ellis"
anno_dict = {'HMDB': f"{SAMPLE}_hmdb.csv",
             "LipidMaps": f"{SAMPLE}_lipidmaps.csv"}
plot_folder = "../../../figures/ellis"


assert os.path.exists(plot_folder), "Cannot continue, the plot_folder does not exist"

me_annots = zd.MetabolomeAnnotations(name=SAMPLE,
                            data_dir=mydata_dir,
                            subfolder=subfolder,
                            annotation_files=anno_dict)
me_annots.load()


adata_mz = ad.read_h5ad(filename=os.path.join(
    mydata_dir, subfolder, f'{SAMPLE}.h5ad'))
print(adata_mz)

me_annots.do_mz_quality_df()

if "spot_quality" not in adata_mz.obs.columns.tolist():
    print(f"computing spot quality across {adata_mz.obs.shape[0]} pixels... ")
    if adata_mz.obs.shape[0] > 30000:
        print("be patient, large dataset in process...")
    me_annots.do_mz_quality_df()   ## NEW
    if not me_annots.mz_quality[me_annots.offsample_column].isnull().all():

        #random_spots = np.random.choice(adata_mz.obs['obs'].to_numpy(), 5000, replace=False)# keep deactivated
        #adata_mz = adata_mz[random_spots,:]# keep deactivated

        adata_mz = adata_mz_mark_spot_quality(
            adata_mz, me_annots.mz_quality,
            cutoff_out_prob=0.5,
            tolerated_fraction=0.05)

        tmpstr = os.path.join(plot_folder,
                              f'{SAMPLE}-pixBefore.png')
        sq.pl.spatial_scatter(
            adata_mz, color=["spot_quality"], img=False,
            title="Before pixel selection", save=os.path.abspath(tmpstr))
        plt.show()

        print("keeping only good pixels...")
        adata_mz = adata_mz[adata_mz.obs["spot_quality"] == "good", :]
        tmpstr = os.path.join(plot_folder,
                              f'{SAMPLE}-pixAfter.png')

        print("visualizing good pixels: will over-write new filtered adata")
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
    print("already computed and selected pixels for adata_mz, nothing to do")
