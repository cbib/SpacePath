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

SAMPLE = "whole-ellis"  # .h5ad saved by -pixel-selection.py
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
    mydata_dir, subfolder, f'{SAMPLE}.h5ad'))
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

#####

prob_out = 0.87
pathwaysObj, me_annots, adata_mz, orig_pws_dict = prepare_harmonized_objects(
        mydata_dir, me_annots, adata_mz, prob_out=prob_out, min_pw_size=2)

print("mz post sync", adata_mz)


# the KEGG id adata
X_metabolites_probs_df = compute_metabolites_probs_df(adata_mz, me_annots.df)

adata = build_kegg_resolved_adata(adata_mz, X_metabolites_probs_df,
                                  annots_df=me_annots.df)

print("KEGG resolved adata: ", adata)

adata.write_h5ad(os.path.join(mydata_dir, subfolder,
                                         f"{SAMPLE}-KEGG_adata.h5ad"))

## compute scores
kpca_scores = zm.ComputePathwaysScores.run(pathwaysObj.to_dict(), adata)

# save as anndata with associated metadata
scores_adata = score_as_adata(kpca_scores, adata,
                              pathwaysObj.pathways_ids,
                              pathwaysObj.piled,
                              orig_pws_dict)

#####
scores_adata.write_h5ad(os.path.join(mydata_dir, subfolder,
                                         f"{SAMPLE}-PathwayScores.h5ad"))


