import sys
import os
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import zonist.data as zd
import zonist.method as zm


#from scripts.bank4zonist.quantifs_methods import compute_mz_binary_df

sys.path.append(os.path.join(os.path.dirname(__file__),
                             '../bank4zonist'))
#
# from scripts.bank4zonist.quantifs_methods import (verification,
#                                               compounds_counts_by_pathway,
#                                               pathways_ok__unwanted_metabolites,
#                                               update_really, compute_mz_binary_df,
#                                               compute_metabolites_probs_df,
#                                               compute_pathways_pseudoscores,
#                                               ugly_long_function)
#
# from scripts.bank4zonist.datasetcomplement import (build_kegg_resolved_adata,
#                                                    yield_specific_typical_pathways)

from quantifs_methods import *
from visualization import *


"""
pathways scores on *bordeaux proteome example data* 

Bash call:
mylist=(UBcalibr1 UBcalibr2 UBcalibr3 UBcalibr4 UBcalibr5);for k in ${mylist[*]}; do echo $k; python xxxxx.py $k; done
"""

SAMPLE = sys.argv[1] # "UBcalibr2"
#"UBcalibr2"
print(SAMPLE)
mydata_dir = "/home/johanna/spatial_thesis/spatial_metabolomics/data"
#"/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"
subfolder = "UBplatform"
#"/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"
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

adata_mz = spatial_standard_computations(adata_mz, 12, plot_folder,
                                        prefix=f"{SAMPLE}-MZ-V2")  # ====> !!

refinement_mz_pws_dict = ugly_long_function(pathwaysObj.piled.copy(),
                                            pathwaysObj.to_dict())

# the KEGG id adata
adata = build_kegg_resolved_adata(adata_mz, annots_df=me_annots.df)
print("KEGG resolved adata: ", adata)
adata = spatial_standard_computations(adata, 12, plot_folder, f"{SAMPLE}-KEGGmetabolitesV2")  # ====> !!


