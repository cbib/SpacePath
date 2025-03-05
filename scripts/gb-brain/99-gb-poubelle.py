import sys
import os
import numpy as np
import pandas as pd
import anndata as ad
import squidpy as sq
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
pathways scores 
"""

SAMPLE = "test635_80-400"

print(SAMPLE)
mydata_dir = "/home/johanna/spatial_thesis/spatial_metabolomics/data"
#"/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"
subfolder = "gb-brain"
anno_dict = {'HMDB': f"{SAMPLE}_hmdb.csv",
             "KEGG": f"{SAMPLE}_kegg.csv"}
plot_folder = "../../../figures/gb-brain"


assert os.path.exists(plot_folder), "Cannot continue, the plot_folder does not exist"

me_annots = zd.MetabolomeAnnotations(name=SAMPLE,
                            data_dir=mydata_dir,
                            subfolder=subfolder,
                            annotation_files=anno_dict)
me_annots.load()


adata_mz = ad.read_h5ad(filename=os.path.join(
    mydata_dir, subfolder, f'{SAMPLE}.h5ad'))
print(adata_mz)

pathwaysObj, me_annots, adata_mz, orig_pws_dict = prepare_harmonized_objects(
    mydata_dir, me_annots, adata_mz, prob_out=0.87) ## todo : bad! this was corrected
print("mz post sync", adata_mz)


adata_mz = compute_moran_i(adata_mz, n_neighs=5, n_perms=50)
moran_df = adata_mz.uns['moranI'].copy()
moran_df['mz'] = moran_df.index.tolist()
moran_df = moran_df.sort_values(by=['I'], ascending=False)

positives = moran_df.head(5)
middle = moran_df.loc[(moran_df["I"] < 0.104) & (moran_df["I"] > 0.049),:]
negatives = moran_df.tail(5)

def plot_that(df, title):
    for i_mz in df["mz"].tolist():
        morani = df.loc[i_mz,"I"].item()
        sq.pl.spatial_scatter(adata_mz,
                              color=[i_mz],
                              img=False,
                              title=f'{title} - {i_mz} - {morani}' )
        plt.show()

plot_that(positives, "positives")
plot_that(middle, "middle")
plot_that(negatives, "neg")

"""
Thoughts about Moran's I values:
- the positive ones: 
   true positives ok
   false posittives: if there are any artifactual clusters:
        -important to avoid "bubbles" of fake tissue outside true tissue  (right lower corner of 80-400 data)
        - many noisy 'stripes': i think impossible to eject or attenuate unfortunately
- the middle (around zero):
   - many many noisy 'stripes'
   - random
   -  I found at least 2 false negatives!!!, because clustered but slightly dispersed
- the negative ones:
    - supposed to be dispersed
    - I found at least 2 false negatives!!!, because clustered but slightly dispersed

Check all this in the smaller datasets!
"""


k=7
adata_mz = spatial_standard_computations(
    adata_mz,
    k=k,
    as_log=True,
    top_variable_features_fraction=0.5)

plot_standard_computs(adata_mz,
    plot_folder=None,
    prefix="oo")#f"{SAMPLE}-MZkmeans")
plt.show()

