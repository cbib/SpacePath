import sys
import os
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import zonist.data as zd
import zonist.method as zm
import squidpy as sq


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
verifying Ravi & Heiland data
use
   python 1-rh-describe.py $SAMPLE

where sample one  of: 
samples = [   "248_T"   , "259_T", "260_T", "262_T", "269_T", "275_T"]

Or use Bash oneliner:
mylist=(248_T 259_T 260_T 262_T 269_T 275_T);for k in ${mylist[*]}; do echo $k; python 1-rh-describe.py $k; done
"""

SAMPLE = sys.argv[1] #  "248_T"
print("============>", SAMPLE)
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

adata_mz = ad.read_h5ad(filename=os.path.join(
    mydata_dir, "heiland-msi", f'MALDI_{SAMPLE}.h5ad'))

print(adata_mz)
adata_mz = get_rid_of_poor_var_mz(adata_mz)  # TODO : highlight this because is proper to this Ravi&Heiland data
print("now:", adata_mz)

me_annots = update_annotation_as_adata_mz(gathered_annots, adata_mz)  # TODO highlight: proper to Ravi&Heiland data as annotation is 1 file for all 6 slices
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

#adata_mz = spatial_standard_computations(adata_mz, 12, plot_folder, prefix=f"{SAMPLE}-MZ-V2")  # ====> !!

refinement_mz_pws_dict = ugly_long_function(pathwaysObj.piled.copy(),
                                            pathwaysObj.to_dict())

# the KEGG id adata
adata = build_kegg_resolved_adata(adata_mz, annots_df=me_annots.df)
print("KEGG resolved adata: ", adata)
#adata = spatial_standard_computations(adata, 12, plot_folder, f"{SAMPLE}-KEGGmetabolitesV2")  # ====> !!

with_moran = False
if with_moran:
    the_mz_for_refinement_list = list(refinement_mz_pws_dict.keys())
    adata_mz = compute_moran_i(adata_mz, 12)
    sub_moran = adata_mz.uns['moranI'].loc[
        adata_mz.uns['moranI'].index.isin(
            [str(x) for x in the_mz_for_refinement_list])]

    sub_moran.to_csv(
        os.path.join(plot_folder, f'{SAMPLE}--moran-subset.tsv'),
        header=True, index=True, sep='\t')

    tmpstr = os.path.join(plot_folder,
                          f'{SAMPLE}--susceptible_of_refinement.png')
    sq.pl.spatial_scatter(adata_mz,
                          color=[str(x) for x in the_mz_for_refinement_list],
                          img=False,
                          save=os.path.abspath(tmpstr))

    print("finished ! ")