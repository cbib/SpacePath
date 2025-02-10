import sys
import os
import numpy as np
import pandas as pd
import anndata as ad
import zonist.data as zd
from scipy.sparse import csr_matrix
from sklearn.utils import sparsefuncs

sys.path.append(os.path.join(os.path.dirname(__file__),
                             '../bank4zonist'))
from quantifs_methods import *
from datasetcomplement import (build_kegg_resolved_adata,
                               yield_specific_typical_pathways)

"""
verifying *bordeaux proteome example data* 
also suppressing bad pixels of UBcalibr4 data

"""

# variables
SAMPLE = sys.argv[1]  #"UBcalibr4"

mydata_dir = "/home/johanna/spatial_thesis/spatial_metabolomics/data"
subfolder = "UBplatform"
#"/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"
anno_dict = {'HMDB': f"{SAMPLE}_hmdb.csv",
             "LipidMaps": f"{SAMPLE}_lipidmaps.csv"}
plot_folder = "../../../figures/platform"

# start
assert os.path.exists(plot_folder), "Cannot continue, the plot_folder does not exist"

me_annots = zd.MetabolomeAnnotations(name=SAMPLE,
                            data_dir=mydata_dir,
                            subfolder=subfolder,
                            annotation_files=anno_dict)
me_annots.load()

adata_mz = ad.read_h5ad(filename=os.path.join(
    mydata_dir, subfolder, f'{SAMPLE}_mx.h5ad'))

if (SAMPLE == "UBcalibr4") and not os.path.exists(
        os.path.join(mydata_dir, subfolder, f'{SAMPLE}-Old.h5ad')):
    # if not done yet:
    # replace adata saved file "_mx.h5ad" by a new one with selected pixels
    # save the Old version in "-Old.h5ad" file
    adata_mz.write(os.path.join(mydata_dir, subfolder, f'{SAMPLE}-Old.h5ad'),
                   compression=None, compression_opts=None)
    obtained_domains_before_purge_pixels_df = visualizing_spatial_domains_before_decision(adata_mz)
    adata_mz = take_decision(adata_mz, obtained_domains_before_purge_pixels_df,
                             not_tissue_spatial_domains=[4])
    adata_mz.write(os.path.join(mydata_dir, subfolder, f'{SAMPLE}_mx.h5ad'),
                   compression=None, compression_opts=None)

verification(adata_mz, me_annots)

print("finished ! ")



