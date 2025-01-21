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
me_annots.load()

## reference mappings to get KEGG ids from HMDB and LipidMaps identifiers
me_refmap = zd.MetaboliteIdsMappings(
    data_dir=mydata_dir,
    subfolder="metabolites-mappings",
    mapping_files={'HMDB': "HMDB-to-KEGG_05-07-24.tsv",
                   'LipidMaps': "LipidMaps-to-KEGG_10-18-24.tsv" }
)
me_refmap.load()

##  matrix of quantifications (==intensities)
print("Opening the quantifications and verifying mz match")
adata_mz = ad.read_h5ad(filename=os.path.join(
    mydata_dir, subdir, quantif_file))
bool_equal = np.all(np.array(
    sorted(me_annots.unfiltered["mz"].unique())) == np.array(
    sorted(adata_mz.var["mz"].tolist())))
print(f"perfect match between mz in me_a.df "
      f"and adata.X: => {bool_equal} <=")


if bool_equal:
    me_annots.convert_ids_to_kegg(me_refmap.table)
else:
    raise ValueError("FATAL: mz do not match between adata_mz and annotation table")


me_annots.treat_mz_multiplets(adata_mz)

print(len(me_annots.df["KEGG"].to_list()))
print(len(me_annots.df["KEGG"].unique()))
print("ok up to here")

print("subset adata_mz for only mz mapping to a KEGG metabolite identifier")
adata_mz = adata_mz[:, me_annots.df["mz"].astype(str).unique()].copy()
print(adata_mz)


## pathways
me_gmt = zd.MetabolicPathways(
    data_dir=mydata_dir,
    subfolder="metabolites-mappings",
    official_file="KEGG_hsa_pathways_compounds_R110.gmt",
)

me_gmt.load()
me_gmt.intersect_with_annotations(me_annots)


## exclude pathways that only have 1 metabolite KEGG:
tdf = me_gmt.piled.copy()

# ## subset adata for useful mz
# no_pathway_score_possible = list()
# for pw in  tdf["pathway_id"].unique():
#     subdf = tdf.loc[tdf["pathway_id"] == pw, :]
#     if len(subdf["KEGG"].unique()) == 1:
#         no_pathway_score_possible.append(pw)
#
# tdf = tdf.loc[~tdf["pathway_id"].isin(no_pathway_score_possible), :]
# unused_v = len(tdf["pathway_id"].unique())


# def get_mz_with_good_moran(adata_mz):
#     """
#     Given the adata_mz object, returns the mz values that have a bad Moran's I score.
#     """
#     sc.pp.highly_variable_genes(adata_mz,
#                                 flavor="seurat_v3",
#                                 check_values=False)
#     sq.gr.spatial_neighbors(adata_mz, n_neighs=3)
#     print("Computing Moran's I")
#     sq.gr.spatial_autocorr(adata_mz, mode="moran", n_perms=100, n_jobs=1)
#     good_mz_str = list(adata_mz.uns['moranI'].loc[
#         adata_mz.uns['moranI']['pval_norm_fdr_bh'] <= 0.05,:].index)
#
#     return good_mz_str
#
#
# good_moran_i_mz_str = get_mz_with_good_moran(adata_mz)
# tdf["mz_str"] = tdf["mz"].astype(str)
# tdf = tdf.loc[tdf["mz_str"].isin(good_moran_i_mz_str), :]
# tdf = tdf.drop(columns=["mz_str"])


### complex cases study
def get_target_mz_pathway(subdf):
    """
    The  subdf represents only *one mz* dataframe (cols: pathway_id, KEGG).
    KEGG is the metabolite identifier.
    For the given distinct KEGG ids and pathway ids in the dataframe,
    first constructs a dict of sets of KEGG ids for each pathway.
    Then, compares each possible pair (k, h) of pathways from the dict,
    to check if there are KEGG ids that are exclusive to k (same for h).
    Returns a dict where the keys are the compared pairs of pathways, e.g.
    {('hsa00230', 'hsa04976'): {'hsa00230': ['C00212', 'C00330'], 'hsa04976': ['C07210']}
    Note: if every (k, h) pathways had same metabolites, output is empty set.
    """
    adict = dict()
    outdict = dict()
    for i, row in subdf.iterrows():
        if row["pathway_id"] not in adict.keys():
            adict[row["pathway_id"]] = {row["KEGG"]}
        adict[row["pathway_id"]].update([row["KEGG"]])
    for k in adict.keys():
        for h in adict.keys():
            if k != h:
                a = adict[k] - adict[h]
                b = adict[h] - adict[k]
                if len(a) > 0 and len(b) > 0:
                    outdict[(k,h)] = {k: list(a), h: list(b)}
    return outdict


one_kegg_one_pathway = 0
one_kegg_multiple_pathways = 0
multiple_kegg_one_pathway = 0
susceptible_of_refinement = 0

for mz in tdf["mz"].unique():
    subdf = tdf.loc[tdf["mz"] == mz, :]
    subdf = subdf.drop(columns=["mz"])
    subdf = subdf.drop_duplicates()
    if subdf.shape[0] == 1:
        one_kegg_one_pathway += 1
    elif (len(subdf["KEGG"].unique()) == 1) and (
            len(subdf["pathway_id"].unique()) > 1):
       one_kegg_multiple_pathways += 1
    elif (len(subdf["KEGG"].unique()) > 1) and (
            len(subdf["pathway_id"].unique()) == 1):
        multiple_kegg_one_pathway += 1
    elif (len(subdf["KEGG"].unique()) > 1) and (
            len(subdf["pathway_id"].unique()) > 1):
        subdf = subdf.drop_duplicates("KEGG")
        out = get_target_mz_pathway(subdf)
        if len(out) > 0:
            susceptible_of_refinement += 1
        else:  # not real multi-kegg multi-pathway, see target_mz_pathway
            multiple_kegg_one_pathway += 1

print(len(tdf["mz"].unique()))
## ok
print("one_kegg_one_pathway: ", one_kegg_one_pathway)
print("one_kegg_multiple_pathways: ", one_kegg_multiple_pathways)
print("multiple_kegg_one_pathway: ", multiple_kegg_one_pathway)
print("susceptible_of_refinement: ", susceptible_of_refinement)

print("ok")



