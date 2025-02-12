import os
import numpy as np
import pandas as pd
import anndata as ad
import zonist.data as zd
import zonist.method as zm

from scipy.sparse import csr_matrix
import sys
import squidpy as sq
sys.path.append(os.path.join(os.path.dirname(__file__),
                             '../bank4zonist'))  # TODO: fix in the same way ailleurs !!
from quantifs_methods import *
from visualization import *

"""
"""

SAMPLE = "whole-ellis"  # f1.h5ad saved by -pixel-selection.py
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
print(adata_mz)

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

## pathways
pathwaysObj = zd.MetabolicPathways(
    data_dir=mydata_dir,
    subfolder="metabolites-mappings",
    official_file="KEGG_hsa_pathways_compounds_R110.gmt")

pathwaysObj.load()
initial_pws_dict = pathwaysObj.to_dict()   # reserver pour plustard ? simplement avoir les tailles des sets du départ.

pathwaysObj, me_annots = zm.SyncPathwaysAndAnnots.run(pathwaysObj,
                                                      me_annots,
                                                      min_pathway_size=2)  # default max is 500

adata_mz = adata_mz[:, me_annots.df["mz"].astype(str).unique()].copy() # todo this is NOT sync, but like so
print("mz post sync", adata_mz)

# the KEGG id adata
adata = build_kegg_resolved_adata(adata_mz, annots_df=me_annots.df)
print("KEGG resolved adata: ", adata)

# scores
kpca_scores = zm.ComputePathwaysScores.run(pathwaysObj, adata)
# save as anndata to have compression
kpca_score_adata =
kpca_score_adata.write_h5ad(os.path.join(mydata_dir, subfolder,
                                         f"{SAMPLE}-PathwayScores.h5ad"))


