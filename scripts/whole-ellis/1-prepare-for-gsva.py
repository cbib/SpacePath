import os
import numpy as np
import pandas as pd
import csv as csv
import zonist.data as zd
import anndata as ad
import scanpy as sc
import squidpy as sq

from sklearn.utils import sparsefuncs

import matplotlib.pyplot as plt
from sphinx.addnodes import not_smartquotable


# for the whole-ellis dataset

# 0) all variables, strings and dictionaries, to be used
mydata_dir = "/home/johanna/spatial_thesis/spatial_metabolomics/data"
subdir = "whole-ellis"
#"/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"

anno_dict = {'HMDB': "whole-ellis_hmdb.csv",
             'LipidMaps': "whole-ellis_lipidmaps.csv"}
quantif_file = "whole-ellis_mx_f1.h5ad"
my_custom_gmt: str = "metabolite_sets_revised_v1.gmt"

matrix_out_file = "whole-ellis_mx_f2.tsv"
# 1) given data
# 1.1) annotations

me_annots = zd.MetabolomeAnnotations(name="demo",
                            data_dir=mydata_dir,
                            subfolder=subdir,
                            annotation_files=anno_dict)
me_annots.load_annotation()

# 1.2) matrix of quantifications (==intensities)
print("Opening the quantifications")
adata = ad.read_h5ad(filename=os.path.join(
    mydata_dir, subdir, quantif_file))
bool_equal = np.all(np.array(
    sorted(me_annots.annotationtbl["mz"].unique())) == np.array(
    sorted(adata.var["var"].tolist())))
print(f"perfect match between mz in me_a.annotationtbl "
      f"and adata.X: {bool_equal}")

# 2) loading mappings data
# the files are provided and location is the one recommended (in data/)
me_refmap = zd.MetaboliteReferenceMappings(
    data_dir=mydata_dir,
    subfolder="metabolites-mappings",
    mapping_files={'HMDB': "HMDB-to-KEGG_05-07-24.tsv",
                   'LipidMaps': "LipidMaps-to-KEGG_10-18-24.tsv" }
)
me_refmap.load_mappings()
me_annots.convert_ids_to_kegg(me_refmap.table)

# using mapppings, keep only the good features in adata
adata = adata[:, me_annots.annotationtbl["mz"].unique()]

print("saving the dense X matrix for GSVA")
out_df = adata.to_df()
out_df.index.name = "mz"
out_df.to_csv(os.path.join(mydata_dir, subdir, matrix_out_file),
              sep='\t', header=True, index=True)
# 3) GMT
# the files are provided and location is the one recommended (in data/)
me_gmt = zd.MetabolicGmt(
    data_dir=mydata_dir,
    subfolder="metabolites-mappings",
    official_file="KEGG_hsa_pathways_compounds_R110.gmt",
    custom_file=my_custom_gmt)

me_gmt.load()

me_gmt.join_mz(me_annots.annotationtbl[["mz", "KEGG"]])

me_gmt.write_to_json(
    dataset_subfolder=subdir, mz_or_kegg_choice="mz")

