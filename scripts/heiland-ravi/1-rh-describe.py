import sys
import os
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import zonist.data as zd
import zonist.method as zm
import squidpy as sq


from scripts.user_routines import prepare_harmonized_objects

#from scripts.bank4zonist.quantifs_methods import compute_mz_binary_df

sys.path.append(os.path.join(os.path.dirname(__file__),
                             '../bank4zonist'))


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


Note: No need of pixel selection (by spot_quality) because already done by authors.

"""

SAMPLE =  "248_T"  # sys.argv[1] #
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

me_annots = update_annotation_tomatch_adatamz(gathered_annots, adata_mz)  # TODO highlight: proper to Ravi&Heiland data as annotation is 1 file for all 6 slices
verification(adata_mz, me_annots)

from scripts.user_routines import prepare_harmonized_objects
pathwaysObj, me_annots, adata_mz, orig_pws_dict = prepare_harmonized_objects(
    mydata_dir, me_annots, adata_mz)
print("mz post sync", adata_mz)

#adata_mz = spatial_standard_computations(adata_mz, k=12)
#plot_standard_computs(adata_mz, k=12, plot_folder=plot_folder, prefix=f"{SAMPLE}-MZ-V2")

# the KEGG id adata
adata = build_kegg_resolved_adata(adata_mz, annots_df=me_annots.df)
print("KEGG resolved adata: ", adata)

#adata = spatial_standard_computations(adata, k=12)
# plot_standard_computs(adata, k=12, plot_folder=plot_folder, f"{SAMPLE}-KEGGmetabolitesV2")  # ====> !!




