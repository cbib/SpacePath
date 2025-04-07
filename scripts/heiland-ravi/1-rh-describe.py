import sys
import os
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import zonist.data as zd
import zonist.method as zm
import squidpy as sq


#from scripts.user_routines import prepare_harmonized_objects

#from scripts.bank4zonist.quantifs_methods import compute_mz_binary_df

sys.path.append(os.path.join(os.path.dirname(__file__),
                             '../bank4zonist'))

sys.path.append(os.path.join(os.path.dirname(__file__),
                             '../'))

from user_routines import prepare_harmonized_objects
from quantifs_methods import *
from visualization import *


"""
verifying Ravi & Heiland data
use
   python 1-rh-describe.py $SAMPLE

where sample one  of: 
samples = [   "248_T"   , "259_T", "260_T", "262_T", "269_T", "275_T"]

Or use Bash oneliner:
mylist=(248_T 259_T 260_T 262_T 269_T 275_T);for k in ${mylist[*]}; do echo $k; python 1-rh-describe.py MALDI_$k; done


Note: No need of pixel selection (by spot_quality) because already done by authors.

"""

SAMPLE =  "MALDI_275_T"  # sys.argv[1] #
print("============>", SAMPLE)
mydata_dir = "/home/johanna/spatial_thesis/spatial_metabolomics/data"
#"/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"
subfolder="heiland-msi"
anno_dict = {'HMDB': "Metaspace_annotations_943peaks.csv"}
plot_folder = "../../../figures/heiland"

assert os.path.exists(plot_folder), "Cannot continue, the plot_folder does not exist"

me_annots = zd.MetabolomeAnnotations(name="heiland",
                            data_dir=mydata_dir,
                            subfolder=subfolder,
                            annotation_files=anno_dict)
me_annots.load()

adata_mz = ad.read_h5ad(filename=os.path.join(
    mydata_dir, subfolder, f'{SAMPLE}.h5ad'))

print(adata_mz)
# adata_mz = get_rid_of_poor_var_mz(adata_mz)  # TODO: del
print("now:", adata_mz)

verification(adata_mz, me_annots)

me_annots.pull_mz_quality_df()
print("filtering by msm and fdr")
ok_mzs_msm_fdr = me_annots.mz_quality.loc[
    (me_annots.mz_quality['msm'] >= 0.01) & (
        me_annots.mz_quality['fdr'] <= 0.2), "mz"].astype(str).tolist()
ok_mzs_msm_fdr = list(set(adata_mz.var.index.tolist()).intersection(
    set(ok_mzs_msm_fdr)
))
adata_mz = adata_mz[:, ok_mzs_msm_fdr].copy()
print("%%%%%%%%%%%", adata_mz)

pathwaysObj, me_annots, adata_mz, orig_pws_dict = prepare_harmonized_objects(
    mydata_dir, me_annots, adata_mz, prob_out=0.5)
print("mz post sync", adata_mz)


## see Morans autocorrel
print("using stats to define spatially relevant mzs")
if not "moranI" in list(adata_mz.uns.keys()):
    adata_mz = compute_moran_i(adata_mz, n_neighs=10, n_perms=50) # TODO : make sure of doing at least 50 perms to have padj ok
    #save with moran I results (avoid to recompute, time expensive)
    adata_mz.write_h5ad(filename=os.path.join(
        mydata_dir, subfolder, f"{SAMPLE}.h5ad"))
print("plotting Morans I results for current sample")
plot_moran_results(adata_mz, plot_folder, SAMPLE)

truth_mz = pathwaysObj.to_dict(option="mz")
foo = []
fii = []
fee = []
fuu = []
for k in truth_mz.keys():
    foo.append(k)
    fee.append(pathwaysObj.get_name(k))
    fii.append("|".join(truth_mz[k]))
    fuu.append(len(truth_mz[k]))

kaka = pd.DataFrame({'pw_id': foo, 'names': fee, 'n': fuu,  'mzs': fii})
kaka.to_csv("kaka.csv")

#adata_mz = spatial_standard_computations(adata_mz, k=12)
#plot_standard_computs(adata_mz,  plot_folder=plot_folder, prefix=f"{SAMPLE}-MZ-V2")

# the KEGG id adata
X_metabolites_probs_df = compute_metabolites_probs_df(adata_mz, me_annots.df)
adata = build_kegg_resolved_adata(adata_mz, X_metabolites_probs_df, annots_df=me_annots.df)

#print("KEGG resolved adata: ", adata)
print("Total pathways :", len(pathwaysObj.piled['pathway_id'].unique()))

refinement_mz_pws_dict = find_features_for_refinement(pathwaysObj.piled.copy())
#adata = spatial_standard_computations(adata, k=12)
# plot_standard_computs(adata,  plot_folder=plot_folder,
# prefix=f"{SAMPLE}-KEGGmetabolitesV2")  # ====> !!




