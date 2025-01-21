import os
import numpy as np
import pandas as pd
import anndata as ad
import zonist.data as zd
from scipy.sparse import csr_matrix
import sys

sys.path.append(os.path.join(os.path.dirname(__file__),
                             '../bank4zonist'))  # TODO: fix in the same way ailleurs !!
from statsembryo1 import *
from datasetcomplement import (build_kegg_resolved_adata,
                               yield_specific_typical_pathways)

"""
"""

SAMPLE = "whole-ellis"  # f1.h5ad saved by 1-adata-pixel-selection.py
mydata_dir = "/home/johanna/spatial_thesis/spatial_metabolomics/data"
#"/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"

subfolder = "whole-ellis"
anno_dict = {'HMDB': "whole-ellis_hmdb.csv",
             'LipidMaps': "whole-ellis_lipidmaps.csv"}
plot_folder = "../../../figures/ellis"

assert os.path.exists(plot_folder), "Cannot continue, the plot_folder does not exist"


me_annots = zd.MetabolomeAnnotations(name="whole-ellis",
                            data_dir=mydata_dir,
                            subfolder=subfolder,
                            annotation_files=anno_dict)
me_annots.load()

adata_mz = ad.read_h5ad(filename=os.path.join(
    mydata_dir, subfolder, f'{SAMPLE}_mx_f1.h5ad'))  # TODO  : how to avoid forgetting _f1

verification(adata_mz, me_annots)

## reference mappings to get KEGG ids from HMDB and LipidMaps identifiers
me_refmap = zd.MetaboliteIdsMappings(
    data_dir=mydata_dir,
    subfolder="metabolites-mappings",
    mapping_files={'HMDB': "HMDB-to-KEGG_05-07-24.tsv",
                   'LipidMaps': "LipidMaps-to-KEGG_10-18-24.tsv" })

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
              me_annots.df["mz"].isin(pathwaysObj.piled["mz"].unique()), :]

adata_mz = adata_mz[:, me_annots_OKdf["mz"].astype(str).unique()].copy()
adata_mz = compute_moran_i(adata_mz, 12)


# the KEGG id adata
adata = build_kegg_resolved_adata(adata_mz, annots_df=me_annots_OKdf)
print("KEGG resolved adata: ", adata)
adata = spatial_standard_computations(
    adata, 12, plot_folder, prefix=f"{SAMPLE}-KEGGmetabolites")

tdf = pathwaysObj.piled.copy()
susceptible_of_refinement_list = ugly_long_function(tdf)

sub_moran = adata_mz.uns['moranI'].loc[
     adata_mz.uns['moranI'].index.isin(
         [str(x) for x in susceptible_of_refinement_list])]

sub_moran.to_csv(
    os.path.join(plot_folder, f'{SAMPLE}--moran-subset.tsv'),
    header=True, index=True, sep='\t')

tmpstr = os.path.join(plot_folder, f'{SAMPLE}--susceptible_of_refinement.png')
sq.pl.spatial_scatter(adata_mz,
                      color=[str(x) for x in susceptible_of_refinement_list],
                      img=False,
                      save=os.path.abspath(tmpstr))

print("finished ! ")
