import sys
import os
import numpy as np
import pandas as pd
import anndata as ad
import zonist.data as zd
from scipy.sparse import csr_matrix
from sklearn.utils import sparsefuncs


import sspa as sspa

sys.path.append(os.path.join(os.path.dirname(__file__),
                             '../bank4zonist'))  # TODO: fix in the same way ailleurs !!
from statsembryo1 import *
from datasetcomplement import (build_kegg_resolved_adata,
                               yield_specific_typical_pathways)

"""
verifying RAvi & Heiland data
use
   python <...>.py $SAMPLE

where sample one  of: 
samples = [   "248_T"   , "259_T", "260_T", "262_T", "269_T", "275_T"]
"""



SAMPLE = "248_T"  #sys.argv[1] #  "248_T"

mydata_dir = "/home/johanna/spatial_thesis/spatial_metabolomics/data"
#"/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"
anno_dict = {'HMDB': "Metaspace_annotations_943peaks.csv"}
plot_folder = "../../../figures/heiland"

assert os.path.exists(plot_folder), "Cannot continue, the plot_folder does not exist"

gathered_annots = zd.MetabolomeAnnotations(name="heiland",
                            data_dir=mydata_dir,
                            subfolder="heiland-msi",
                            annotation_files=anno_dict)
gathered_annots.load()


## reference mappings to get KEGG ids from HMDB and LipidMaps identifiers
me_refmap = zd.MetaboliteIdsMappings(
    data_dir=mydata_dir,
    subfolder="metabolites-mappings",
    mapping_files={'HMDB': "HMDB-to-KEGG_05-07-24.tsv",
                   'LipidMaps': "LipidMaps-to-KEGG_10-18-24.tsv" }
)
me_refmap.load()
#

adata_mz = ad.read_h5ad(filename=os.path.join(
    mydata_dir, "heiland-msi", f'MALDI_{SAMPLE}.h5ad'))

# TODO: define if this is more useful later when KEGG adata is created
adata_mz = spatial_standard_computations(adata_mz, 12,
                                         plot_folder=plot_folder, prefix=SAMPLE)
#adata_mz = compute_moran_i(adata_mz, 5)

me_annots = update_annotation_as_adata_mz(gathered_annots, adata_mz)
verification(adata_mz, me_annots)

me_annots.convert_ids_to_kegg(me_refmap.table)
me_annots.treat_mz_multiplets(adata_mz)

print(len(me_annots.df["KEGG"].to_list()))
print(len(me_annots.df["KEGG"].unique()))
print("ok up to here")

print("subset adata_mz for only mz mapping to a KEGG metabolite identifier")
adata_mz = adata_mz[:, me_annots.df["mz"].astype(str).unique()].copy()
print(adata_mz)

## pathways
pathwaysObj = zd.MetabolicPathways(
    data_dir=mydata_dir,
    subfolder="metabolites-mappings",
    official_file="KEGG_hsa_pathways_compounds_R110.gmt",
)

pathwaysObj.load()
pathwaysObj.intersect_with_annotations(me_annots)

bobo = me_annots.df.loc[me_annots.df["mz"].isin(
                           pathwaysObj.piled["mz"].unique()), :].copy()
print(len(bobo["KEGG"].unique()))  # 282
print(len(bobo["KEGG"].to_list()))   # 282
del bobo
# the KEGG id adata
adata = build_kegg_resolved_adata(
    adata_mz,
    annots_df=me_annots.df.loc[
              me_annots.df["mz"].isin(pathwaysObj.piled["mz"].unique()), :])


tdf = pathwaysObj.piled.copy()


# sub_moran = adata_mz.uns['moranI'].loc[
#      adata_mz.uns['moranI'].index.isin([str(x) for x in susceptible_of_refinement_list])
# ]
#
# sub_moran.to_csv(
#     os.path.join(plot_folder, f'{SAMPLE}--moran-subset.tsv'),
#     header=True, index=True, sep='\t'
# )







