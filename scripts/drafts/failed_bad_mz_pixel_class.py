import os
import numpy as np
import pandas as pd
import csv as csv
import zonist.data as zd
import anndata as ad
import scanpy as sc
import squidpy as sq

import matplotlib.pyplot as plt

from sklearn.cluster import KMeans



# for the whole-ellis dataset

# 0) all variables, strings and dictionaries, to be used
mydata_dir = "/home/johanna/spatial_thesis/spatial_metabolomics/data"
subdir = "whole-ellis"
#"/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"

anno_dict = {'HMDB': "whole-ellis_hmdb.csv",
             'LipidMaps': "whole-ellis_lipidmaps.csv"}
quantif_file = "whole-ellis_mx.h5ad"
my_custom_gmt: str = ".gmt"


def generate_coordinates(pixels_names):
    coord_l = list()
    for k in pixels_names:   # eg 'x0_y0'
        items = k.split("_")
        coord_l.append(np.array([int(items[0].replace("x", "")),
                       int(items[1].replace("y", ""))]))
    return np.array(coord_l)

    
def yield_submatrix_offsample_mz(me_annots, adata):
    bad_mz_df = me_annots.feat_classif.loc[
                me_annots.feat_classif["offSample"] == True, :]
    adata_subset = adata[:, bad_mz_df["mz"].astype(str).tolist()]
    return adata_subset


# 1) given data
# 1.1) annotations

me_annots = zd.MetabolomeAnnotations(name="demo",
                            data_dir=mydata_dir,
                            subfolder=subdir,
                            annotation_files=anno_dict)
me_annots.load_annotation()



# 1.2) matrix of quantifications (==intensities)
print("opening the quantifications")
adata = ad.read_h5ad(filename=os.path.join(
    mydata_dir, subdir, quantif_file))
print(adata)
bool_equal = np.all(np.array(
    sorted(me_annots.annotationtbl["mz"].unique())) == np.array(
    sorted(adata.var["var"].tolist())))
print(f"perfect match between mz in me_a.annotationtbl "
      f"and adata.X: {bool_equal}")


##  add attributes to be spatially plottable


coordinates_arr = generate_coordinates(adata.obs["obs"].tolist())
adata = ad.AnnData(adata.X,
                    obs = adata.obs, var= adata.var,
                    obsm = {"spatial": coordinates_arr},
                    varm = {"metabolites": adata.var["var"].to_numpy()}
                    )

adata.uns['spatial'] = { 'MY_SLICE7_NAME' : { 'scalefactors': {'fiducial_diameter_fullres': 100,#144.48769000651953,
   'spot_diameter_fullres': 1, #89.44476048022638,
   'tissue_hires_scalef': 0.5, #0.17011142,
   'tissue_lowres_scalef': 0.5  }} }



# try to give a classification based only on sum of intensities of "bad mz"
badmz_adata = yield_submatrix_offsample_mz(me_annots, adata)

badmz_adata.obs["sum"] = badmz_adata.X.sum(axis=1)

badmz_adata.obs["myclass"] = "sum_bad_mz:low"
cutoff_percent_sum_offsample_mz = 70  # %
cut_value = np.percentile(badmz_adata.obs["sum"].to_numpy(),
                       cutoff_percent_sum_offsample_mz)
badmz_adata.obs.loc[badmz_adata.obs["sum"] > cut_value,
    "myclass"] = "sum_bad_mz:high"

sq.pl.spatial_scatter(badmz_adata,
                      color="myclass", img=False)
plt.show()


"""
print(badmz_adata)
sc.pp.pca(badmz_adata)
print(badmz_adata.obsm["X_pca"][:, 0])
t = 33
kmeans = KMeans(n_clusters=4,
                random_state=0).fit(
    badmz_adata.obsm["X_pca"])
badmz_adata.obs["kmean"] = kmeans.labels_.astype(str)

sq.pl.spatial_scatter(badmz_adata,
                      color="kmean", img=False)
"""
