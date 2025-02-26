import sys
import os
import numpy as np
import pandas as pd
import anndata as ad
import zonist.data as zd
import zonist.method as zm


sys.path.append(os.path.join(os.path.dirname(__file__),
                             '../bank4zonist'))
from quantifs_methods import *



"""
pathways scores on *bordeaux proteome example data* 
"""

SAMPLE = "UBcalibr4" #sys.argv[1]
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
    mydata_dir, subfolder, f'{SAMPLE}_mx.h5ad'))
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
me_annots.treat_mz_multiplets(adata_mz)  # use the one not yet sync, but run verif  before!

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
#adata_mz = spatial_standard_computations(adata_mz, k=12)
# xxxx(adata_mz, k=12, plot_folder=plot_folder, f"{SAMPLE}KeptMz")


refinement_mz_pws_dict = ugly_long_function(pathwaysObj.piled.copy(),
                                            pathwaysObj.to_dict())

adata = build_kegg_resolved_adata(adata_mz, annots_df=me_annots.df)
#adata = spatial_standard_computations(adata, 12, plot_folder, f"{SAMPLE}-KEGGclean")


