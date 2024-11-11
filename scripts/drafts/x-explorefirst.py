import os
import numpy as np
import pandas as pd
import csv as csv
import zonist.data as zd
import anndata as ad

from metaspace_converter import metaspace_to_anndata
import scanpy as sc

# comparing the anots from Download vs. from metaspaceConverter

# 0) all variables, strings and dictionaries, to be used
mydata_dir = "/home/johanna/spatial_thesis/spatial_metabolomics/data"
subdir = "tmp"
#"/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"

anno_dict = {'HMDB': "whole-ellis_hmdb.csv",
             'LipidMaps': "whole-ellis_lipidmaps.csv"}

my_custom_gmt: str = "user-defined-pathways-or-sets.gmt"

metaspace_registry = "2022-08-16_10h14m23s"
api_versions = {'HMDB': ("HMDB", "v4"),
               'LipidMaps': ("LipidMaps", "2017-12-12")}
api_fdr = 0.1  # to go fast but change


# **) blah
if not os.path.exists(os.path.join(mydata_dir, subdir, f"HMDB-CONVANNOTS.tsv")):
#     adata.var.to_csv(os.path.join(mydata_dir, subdir, f"{key}-CONVANNOTS.tsv"), sep='\t')
    for key in api_versions.keys():
         adata = metaspace_to_anndata(
             dataset_id=metaspace_registry,
             fdr=api_fdr,
             database=api_versions[key] #("HMDB", "v4")
         )
         adata.var.to_csv(os.path.join(mydata_dir, subdir, f"{key}-CONVANNOTS.tsv"), sep='\t')
#
# print("yeah")
df1 = pd.read_csv(os.path.join(mydata_dir, subdir, f"HMDB-CONVANNOTS.tsv"), sep='\t')
df2 = pd.read_csv(os.path.join(mydata_dir, subdir, f"LipidMaps-CONVANNOTS.tsv"), sep='\t')

df1['mz_round'] = np.around(df1['mz'].to_numpy(), 9)
df2['mz_round'] = np.around(df2['mz'].to_numpy(), 9)
print(df1.shape)
print(df2.shape)
print(df2.columns)
# 1) annotations

me_annots = zd.MetabolomeAnnotations(name="demo",
                            data_dir=mydata_dir,
                            subfolder=subdir,
                            annotation_files=anno_dict)
me_annots.load_annotation()

me_annots.cutoff_prob_out = 0.9
me_annots.discard_out_of_sample_mz()  # deactivate just to verify equality

me_annots.annotationtbl["mz_round"] = np.around(
    me_annots.annotationtbl["mz"].to_numpy(),
    9)


# 99) xxxxxx
print('this is bad:')
theunion = set(df1["mz_round"].tolist()).union(set(
df2["mz_round"].tolist()
))
commonmz = set(
    me_annots.annotationtbl["mz"].astype(float).tolist()).intersection(set(
    theunion
))
print(len(me_annots.annotationtbl["mz"].unique().astype(float).tolist()))
print(len(commonmz))

print("this is good (rounded mz)")
theunion = set(df1["mz_round"].tolist()).union(set(
df2["mz_round"].tolist()
))
commonmz = set(
    me_annots.annotationtbl["mz_round"].astype(float).tolist()).intersection(set(
    theunion
))
print(len(me_annots.annotationtbl["mz_round"].unique().astype(float).tolist()))
print(len(commonmz))
# ***)
print("this is also good (formula)")

theunionF = set(df1["formula"].tolist()).union(set(
df2["formula"].tolist()
))
commonF = set(
    me_annots.annotationtbl["formula"].tolist()).intersection(set(
    theunionF
))
print(len(me_annots.annotationtbl["formula"].unique().tolist()))
print(len(commonF))
print(f"I used off sample cutoff ====> {me_annots.cutoff_prob_out}")
print("?????")
# ouf = list(commonF)[:5]
# ko = me_annots.annotationtbl[["mz", "formula", "rawOffSampleProb"]]
# ko = ko.drop_duplicates()
# la = ko.join(df1[["mz","formula", ""]])
# s= 1
# END
print("whatever")
def print_this(dflist):
    beforeformulas = set()
    formulasuniques = set()
    for df in dflist:
        beforeformulas = beforeformulas.union(
            set(df['formula'].unique().tolist())
        )
        df = df.loc[df['offSample'] != True, :]  # keep that are not off sample
        formulasuniques = formulasuniques.union(set(
            df['formula'].tolist()
        ))
    print("before", len(beforeformulas))

    print("after", len(formulasuniques))

print_this([df1, df2])
print_this([df1.loc[df1['fdr'] <= 0.2,:],
            df2.loc[df2['fdr'] <= 0.2,:]])


