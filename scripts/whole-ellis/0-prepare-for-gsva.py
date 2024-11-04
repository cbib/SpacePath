import os
import numpy as np
import pandas as pd
import csv as csv
import zonist.data as zd
import anndata as ad

from metaspace_converter import metaspace_to_anndata
import scanpy as sc

# TODO: do all the steps to obtain the files for running GSVA
# for the whole-ellis dataset

# 0) all variables, strings and dictionaries, to be used
mydata_dir = "/home/johanna/spatial_thesis/spatial_metabolomics/data"
subdir = "whole-ellis"
#"/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"

anno_dict = {'HMDB': "whole-ellis_hmdb.csv",
             'LipidMaps': "whole-ellis_lipidmaps.csv"}

my_custom_gmt: str = "user-defined-pathways-or-sets.gmt"

metaspace_registry = "2022-08-16_10h14m23s"
api_versions = {'HMDB': ("HMDB", "v4"),
               'LipidMaps': ("LipidMaps", "2017-12-12")}
api_fdr = 0.5  # to go fast but change

# 1) annotations

me_annots = zd.MetabolomeAnnotations(name="demo",
                            data_dir=mydata_dir,
                            subfolder=subdir,
                            annotation_files=anno_dict)
me_annots.load_annotation()
me_annots.discard_out_of_sample_mz()

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


# 4) import the quantification matrix from the metaspace-converter API
print("import the quantification matrix from the metaspace-converter API")

list_of_mx = []
for key in api_versions.keys():
    adata = metaspace_to_anndata(
        dataset_id=metaspace_registry,
        fdr=api_fdr,
        database=api_versions[key] #("HMDB", "v4")
    )
    mz = adata.var['mz'].tolist()  #
    p1 = "x" + adata.obs['ion_image_pixel_x'].astype(str).to_numpy()
    p2 = "_y" + adata.obs['ion_image_pixel_y'].astype(str).to_numpy()
    p3 = p1 + p2
    Xdf = pd.DataFrame(adata.X)
    Xdf.index = list(p3)
    Xdf.columns = mz
    list_of_mx.append(Xdf)
    # end for

retained_mz = {}


theunion = set(list_of_mx[0].columns.tolist()).union(set(
list_of_mx[1].columns.tolist()
))
commonmz = set(
    me_annots.annotationtbl["mz"].astype(float).tolist()).intersection(set(
    theunion
))
print(len(commonmz))
for k in commonmz:#me_annots.annotationtbl["mz"].astype(float).tolist():
    if str(k) not in retained_mz.keys():
        try:
            retained_mz[str(k)] = list_of_mx[0][[k]].to_numpy()
        except KeyError:
            retained_mz[str(k)] = list_of_mx[1][[k]].to_numpy()


out_mx = pd.DataFrame(retained_mz)
out_mx.index = list_of_mx[1].index  # indistinct if 1 or 0 of the list
out_mx.to_csv(os.path.join(mydata_dir, subdir, "finito.tsv"), sep='\t')

# END
# this is continued with the .sh that runs the Rscript 
# TODO: try to do a rpy2 call of GSVA (.py script in zonist)
