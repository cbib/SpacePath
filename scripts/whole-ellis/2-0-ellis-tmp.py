import os
import numpy as np
import pandas as pd
import anndata as ad
import zonist.data as zd
from scipy.sparse import csr_matrix
from sklearn.utils import sparsefuncs
import squidpy as sq
import sys

import sspa as sspa

sys.path.append(os.path.join(os.path.dirname(__file__),
                             '../bank4zonist'))  # TODO: fix in the same way ailleurs !!
from statsembryo1 import *
from datasetcomplement import (build_kegg_resolved_adata,
                               yield_specific_typical_pathways)

"""
pathways scores on whole brain data (Dr Ellis) !

"""

SAMPLE = "whole-ellis"  # .h5ad saved by  -pixel-selection.py
mydata_dir = "/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"
subfolder = "whole-ellis"
anno_dict = {'HMDB': "whole-ellis_hmdb.csv",
             'LipidMaps': "whole-ellis_lipidmaps.csv"}
plot_folder = "../../../figures/ellis"

assert os.path.exists(plot_folder), "Cannot continue, the plot_folder does not exist"


me_annots = zd.MetabolomeAnnotations(name="demo",
                                     data_dir=mydata_dir,
                                     subfolder=subfolder,
                                     annotation_files=anno_dict)
me_annots.load()

adata_mz = ad.read_h5ad(filename=os.path.join(
    mydata_dir, subfolder, f'{SAMPLE}_mx_f1.h5ad'))

print(adata_mz.var.head(2))


verification(adata_mz, me_annots)


## reference mappings to get KEGG ids from HMDB and LipidMaps identifiers
me_refmap = zd.MetaboliteIdsMappings(
    data_dir=mydata_dir,
    subfolder="metabolites-mappings",
    mapping_files={'HMDB': "HMDB-to-KEGG_05-07-24.tsv",
                   'LipidMaps': "LipidMaps-to-KEGG_10-18-24.tsv" }
)
me_refmap.load()

me_annots.convert_ids_to_kegg(me_refmap.table)
me_annots.treat_mz_multiplets(adata_mz)

# print(len(me_annots.df["KEGG"].to_list())); print(len(me_annots.df["KEGG"].unique()))

## pathways
pathwaysObj = zd.MetabolicPathways(
    data_dir=mydata_dir,
    subfolder="metabolites-mappings",
    official_file="KEGG_hsa_pathways_compounds_R110.gmt")

pathwaysObj.load()
pathwaysObj.intersect_with_annotations(me_annots)

print("subset adata_mz for only mz mapping to a KEGG metabolite identifiers that in turn map to pathways")
me_annots_OKdf = me_annots.df.loc[
              me_annots.df["KEGG"].isin(pathwaysObj.piled["KEGG"].unique()), :]

adata_mz = adata_mz[:, me_annots_OKdf["mz"].astype(str).unique()].copy()

# the KEGG id adata
adata = build_kegg_resolved_adata(adata_mz, annots_df=me_annots_OKdf)
print("KEGG resolved adata: ", adata)

adata = spatial_standard_computations(adata, pca_ncomps=12,
                                  plot_folder=plot_folder,
                              prefix=f'{SAMPLE}')



print("tmp ended")
