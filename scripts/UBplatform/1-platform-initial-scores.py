import sys
import os
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import zonist.data as zd
import zonist.method as zm

sys.path.append(os.path.join(os.path.dirname(__file__),
                             '../bank4zonist'))

sys.path.append(os.path.join(os.path.dirname(__file__),
                             '../'))

from user_routines import prepare_harmonized_objects

from quantifs_methods import *
from visualization import *


"""
pathways scores on *bordeaux proteome example data* 

Bash call:
mylist=(UBcalibr1 UBcalibr2 UBcalibr3 UBcalibr4 UBcalibr5);for k in ${mylist[*]}; do echo $k; python xxxxx.py $k; done
"""

SAMPLE = "UBcalibr2"  #  sys.argv[1]

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
print(adata_mz)

print(np.array(adata_mz.X.todense()).shape)

if SAMPLE == "UBcalibr2":
    prob_out = 1  # §§ ok !!for this small dataset ==> when I excluded even with 0.95, Spatially Relevant mz went away !!!
else:
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

# save all relevant objects
scores_adata.write_h5ad(os.path.join(mydata_dir, subfolder,
                                         f"{SAMPLE}-PathwayScores.h5ad"))


print("Saved new data into:", os.path.join(mydata_dir, subfolder))
