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
from mini_tic import *



SAMPLE = "UBcalibr2"  #  sys.argv[1]

print(SAMPLE)
mydata_dir = "/home/johanna/spatial_thesis/spatial_metabolomics/data"
#"/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"
subfolder = "UBplatform"
anno_dict = {'HMDB': f"{SAMPLE}_hmdb.csv",
             "LipidMaps": f"{SAMPLE}_lipidmaps.csv"}
plot_folder = "../../../figures/platform"


assert os.path.exists(plot_folder), "Cannot continue, the plot_folder does not exist"

me_annots = zd.MetabolomeAnnotations(name=SAMPLE,
                            data_dir=mydata_dir,
                            subfolder=subfolder,
                            annotation_files=anno_dict)
me_annots.load()



adata_mz = ad.read_h5ad(filename=os.path.join(
    mydata_dir, subfolder, f'{SAMPLE}.h5ad'))
print(adata_mz)

foo = adata_mz.var.copy()
print(foo.loc[(foo['mz'] > 830.569) & (foo['mz'] < 830.57), :])

popo = "830.569432"
import squidpy as sq
sq.pl.spatial_scatter(
            adata_mz, color=popo, img=False, shape="square",
            title="Before tic")
plt.show()

print("normalising")
normalized_mx = tic_normalization(adata_mz.X.copy(),
                                  adata_mz.obs['obs'].tolist(),
                                  adata_mz.var['mz'].tolist())
print("ok")
newada = adata_mz.copy()
newada.X = csr_matrix(normalized_mx)

sq.pl.spatial_scatter(
            newada, color=popo, img=False, shape="square",
            title="after ")
plt.show()