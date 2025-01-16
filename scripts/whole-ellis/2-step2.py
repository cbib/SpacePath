import os
import numpy as np
import anndata as ad

import zonist.data as zd

import scanpy as sc
import squidpy as sq

import pandas as pd


quantif_file = "whole-ellis_mx_f1.h5ad"  # saved by 1-adata-pixel-selection.py
mydata_dir = "/home/johanna/spatial_thesis/spatial_metabolomics/data"
subdir = "whole-ellis"
anno_dict = {'HMDB': "whole-ellis_hmdb.csv",
             'LipidMaps': "whole-ellis_lipidmaps.csv"}

# objects
## annotations as downloaded from Metaspace
me_annots = zd.MetabolomeAnnotations(name="demo",
                            data_dir=mydata_dir,
                            subfolder=subdir,
                            annotation_files=anno_dict)
me_annots.load()

## reference mappings to get KEGG ids from HMDB and LipidMaps identifiers
me_refmap = zd.MetaboliteIdsMappings(
    data_dir=mydata_dir,
    subfolder="metabolites-mappings",
    mapping_files={'HMDB': "HMDB-to-KEGG_05-07-24.tsv",
                   'LipidMaps': "LipidMaps-to-KEGG_10-18-24.tsv" }
)
me_refmap.load()

##  matrix of quantifications (==intensities)
print("Opening the quantifications and verifying mz match")
adata_mz = ad.read_h5ad(filename=os.path.join(
    mydata_dir, subdir, quantif_file))
bool_equal = np.all(np.array(
    sorted(me_annots.unfiltered["mz"].unique())) == np.array(
    sorted(adata_mz.var["mz"].tolist())))
print(f"perfect match between mz in me_a.df "
      f"and adata.X: => {bool_equal} <=")


# combine objects information
## conversion of ids to KEGG
if bool_equal:
    me_annots.convert_ids_to_kegg(me_refmap.table)
else:
    raise ValueError("FATAL: mz do not match between adata_mz and annotation table")

me_annots.treat_mz_multiplets(adata_mz)

print(len(me_annots.df["KEGG"].to_list()))
print(len(me_annots.df["KEGG"].unique()))

me_gmt = zd.MetabolicPathways(
    data_dir=mydata_dir,
    subfolder="metabolites-mappings",
    official_file="KEGG_hsa_pathways_compounds_R110.gmt",
)

me_gmt.load()
me_gmt.intersect_with_annotations(me_annots)

print("ok")




