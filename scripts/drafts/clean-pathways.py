import sys
import os
import numpy as np
import pandas as pd
import anndata as ad
import zonist.data as zd
import zonist.method as zm
from scipy.sparse import csr_matrix
from sklearn.utils import sparsefuncs
import squidpy as sq

import sspa as sspa


sys.path.append(os.path.join(os.path.dirname(__file__),
                             '../bank4zonist'))
from quantifs_methods import *
from datasetcomplement import (# build_kegg_resolved_adata,  !!was moved to quantifs_methods
                               yield_specific_typical_pathways)


"""
pathways scores on *bordeaux proteome example data* 
"""

SAMPLE = "UBcalibr2"#sys.argv[1]
#"UBcalibr2"
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
verification(adata_mz, me_annots)  # todo: method me_annots.check_match(adata_mz)

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


pws01, anns01 = zm.SyncPathwaysAndAnnots.run(pathwaysObj.model_copy(),
                                             me_annots.model_copy(),
                                             min_pathway_size=1,
                                             max_pathway_size=1)
print(len(pws01.piled['pathway_id'].tolist()))
print(len(pws01.piled['pathway_id'].unique()))
print("-----------------")
pathwaysObj, me_annots = zm.SyncPathwaysAndAnnots.run(pathwaysObj,
                                                      me_annots,
                                                      min_pathway_size=2)  # default max is 500
print(len(pathwaysObj.piled['pathway_id'].tolist()))
print(len(pathwaysObj.piled['pathway_id'].unique()))

print("subset adata_mz for only mz mapping to a KEGG metabolite identifiers that in turn map to pathways")
print("-----------------")
adata_mz = adata_mz[:, me_annots.df["mz"].astype(str).unique()].copy()

# the KEGG id adata
adata = build_kegg_resolved_adata(adata_mz, annots_df=me_annots.df)
print("KEGG resolved adata: ", adata)

###
# compute kpca
##
kpca_scores = zm.ComputePathwaysScores.run(pathwaysObj.to_dict(), adata)
print(kpca_scores.shape)
print(kpca_scores.iloc[:3,:5])

# # TODO: is it always needed to shift the distr to the right ???

biologically_chosen_pathways = yield_specific_typical_pathways()

#
# kpca_score_adata = ad.AnnData(
#     X=csr_matrix(kpca_scores),
#     var=pd.DataFrame(
#         {'pathway_id' : kpca_scores.columns.tolist()},
#         index=kpca_scores.columns.tolist()),
#     obs=adata.obs,
#     obsm=adata.obsm,
#     uns=adata.uns
# )
#












