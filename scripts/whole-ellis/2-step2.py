import os
import numpy as np
import anndata as ad

import zonist.data as zd

import scanpy as sc
import squidpy as sq

import pandas as pd


quantif_file = "whole-ellis_mx_f1.h5ad"  # saved by 1-adata-pixel-selection.py
mydata_dir = "/home/johanna/spatial_thesis/spatial_metabolomics/data"
subdir = "whole-ellis"
anno_dict = {'HMDB': "whole-ellis_hmdb.csv",
             'LipidMaps': "whole-ellis_lipidmaps.csv"}

# objects
## annotations as downloaded from Metaspace
me_annots = zd.MetabolomeAnnotations(name="demo",
                            data_dir=mydata_dir,
                            subfolder=subdir,
                            annotation_files=anno_dict)
me_annots.load_annotation()

## reference mappings to get KEGG ids from HMDB and LipidMaps identifiers
me_refmap = zd.MetaboliteReferenceMappings(
    data_dir=mydata_dir,
    subfolder="metabolites-mappings",
    mapping_files={'HMDB': "HMDB-to-KEGG_05-07-24.tsv",
                   'LipidMaps': "LipidMaps-to-KEGG_10-18-24.tsv" }
)
me_refmap.load_mappings()

##  matrix of quantifications (==intensities)
print("Opening the quantifications and verifying mz match")
adata_mz = ad.read_h5ad(filename=os.path.join(
    mydata_dir, subdir, quantif_file))
bool_equal = np.all(np.array(
    sorted(me_annots.annotationtbl["mz"].unique())) == np.array(
    sorted(adata_mz.var["mz"].tolist())))
print(f"perfect match between mz in me_a.annotationtbl "
      f"and adata.X: => {bool_equal} <=")


# combine objects information
## conversion of ids to KEGG
if bool_equal:
    me_annots.convert_ids_to_kegg(me_refmap.table)
else:
    raise ValueError("FATAL: mz do not match between adata_mz and annotation table")

from sklearn.utils import sparsefuncs

def resolve_mz_multi_kegg(adata: ad.AnnData,
                          me_annots: zd.MetabolomeAnnotations) ->:
    """
    HMDB curators do efforts to map one-to-one, to KEGG ids.
    However, still many HMDB ids match to one single KEGG id, yielding
    multiple m/z <> repeated KEGG id in our piled table of annotations.
    To resolve such multi-m/z > 1 KEGG group, this function keeps the
    m/z (of such group) that has the highest variance of the quantifications.
    output: List[tuple] a list of tuples: (mz (float), kegg (str))
    """
    mean_arr, v_arr = sparsefuncs.mean_variance_axis(adata.X, axis=0)
    adata.var['variance'] = v_arr

    df2 = me_annots.annotationtbl.groupby("KEGG", as_index=False).size()
    unique_kegg_ids = df2.loc[df2['size'] == 1, "KEGG"].tolist()
    ok_mz_kegg_tuples = [tuple(t) for t in me_annots.annotationtbl.loc[
        me_annots.annotationtbl['KEGG'].isin(unique_kegg_ids),
        ["mz", "KEGG"]].to_numpy()]
    repeated_kegg_ids = df2.loc[df2['size'] >= 2, "KEGG"].tolist()
    for a_kegg in repeated_kegg_ids:
        multi_annot_df = me_annots.annotationtbl.loc[
                me_annots.annotationtbl["KEGG"] == a_kegg, :]
        df_with_v_metric = adata.var.loc[
               adata.var['mz'].isin(multi_annot_df["mz"].tolist()),:]
        df_with_v_metric = df_with_v_metric.sort_values(
            by='variance', ascending=False)
        mz_here = df_with_v_metric["mz"].tolist()[0] # first item mz
        ok_mz_kegg_tuples.append([mz_here, a_kegg])
        # end for
    return ok_mz_kegg_tuples


