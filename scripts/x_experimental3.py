import os
import pandas as pd
import numpy as np
import zonist.data as zd
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import squidpy as sq
from local_helpers.local_helpers import provide_annotation_dict
import zonist.processing.measures_methods as zmm

print("")
print("DOING WITH HMDB, but using the table taht I created in CompoundsIdentifiers repo, passing to CHEBI with my conversion thing")
print("")
# udf = pd.read_csv(
#     "/home/johanna/spatial_thesis/moreIdsAndOntologiesMetabolites/lipids_from_swiss_database.tsv",
#     encoding="latin1", sep="\t", header=0)
# print(udf.columns)

udf = pd.read_csv(
    "/home/johanna/spatial_thesis/moreIdsAndOntologiesMetabolites/tables/HMDB-to-CHEBI_05-07-24.tsv",
     sep="\t", header=0)
print(udf.columns)

jjj = udf[['HMDB', 'CHEBI']]
jjj = jjj.dropna()
print('shape of the HMDB »« chebi table')
print(jjj.shape)

#################
with open("../../data/TD-Lipido-635/TD-Lipido-635-WIKI-NO.txt", "r") as f:
    WIKI_NO = f.readlines()
WIKI_NO.pop(0)
WIKI_NO = [float(t.strip()) for t in WIKI_NO]

with open("../../data/TD-Lipido-635/TD-Lipido-635-clusters-initial.txt", "r") as f:
    initial_clusters = f.readlines()

initial_clusters.pop(0)
initial_clusters = [float(t.strip().split(",")[1]) for t in initial_clusters]


mydata_dir = "../../data/"
SAMPLE = "TD-Lipido-635"
subfolder = SAMPLE

anno_dict = provide_annotation_dict(SAMPLE, ["HMDB"])
me_annots = zd.MetabolomeAnnotations(
    name=SAMPLE, data_dir=mydata_dir,
    subfolder=subfolder, annotation_files=anno_dict)
me_annots.load()
me_annots.usable_id = "moleculeIds"
me_annots.verify()
#me_annots.keep_unique_mz_id_tuples(adata_mz.var)
#me_annots.remove_redundant_aliases()

oldyAnnot = me_annots.model_copy()
print(oldyAnnot.df.shape)

me_annots, mz_to_keep_part1 = zmm.part_a__mz_quality_mark(
    me_annots, msm_cutoff=0.1, offSample_cutoff=0.35, fdr_cutoff=0.5)

good_mz = me_annots.mz_quality.loc[
    me_annots.mz_quality['quality_msm_off_fdr'] == "good",
    "mz"].astype(float).unique()

me_annots.df = me_annots.df.loc[me_annots.df['mz'].astype(float).isin(
  good_mz),:]

print(me_annots.df.shape)
print(len(me_annots.df["mz"].unique()))

adata_mz = ad.read_h5ad(filename=os.path.join(
            mydata_dir, subfolder, f'{SAMPLE}.h5ad'))

print(adata_mz.uns["history"])
print("MANDATORY: take away bad quality m/z") # todo: how to make an alert ?
adata_mz = adata_mz[:,
           adata_mz.var.loc[adata_mz.var['final_quality'] == "good",
           "mz"].astype(str).tolist()]

print(adata_mz)
me_annots, adata_mz = zmm.match_adata_mz_and_annots(me_annots, adata_mz)

xperimentCSV = pd.read_csv("xperimentCSV.csv", index_col=None)

print("markers unfiltered", xperimentCSV.shape)
xperimentCSV = xperimentCSV.loc[
    (abs(xperimentCSV['logfoldchanges']) >= 0.8) & (xperimentCSV['pvals_adj'] <=0.05), :
].copy()
print("markers after filter $$", xperimentCSV.shape)

print("intersect with the spatially relevant, top moranI, and save ")
spatially_interesting = adata_mz.var.copy()
spatially_interesting = spatially_interesting.loc[spatially_interesting['I'] >= 0.4,:].copy()

xperimentCSV = xperimentCSV.loc[xperimentCSV['names'].astype(float).isin(
    spatially_interesting['mz'].astype(float)
),:]

xperimentCSV['names'] = xperimentCSV['names'].astype(float)

print(len(xperimentCSV['names'].unique()))

# take those annots that are NOT in pathwyas:
mydf = me_annots.df.loc[me_annots.df['mz'].isin(
 list(set(WIKI_NO))), :]

# and that are top markers:
mydf = mydf.loc[mydf['mz'].isin(
xperimentCSV['names'].unique()
),:]

print(len(mydf['mz'].unique()))
print(mydf.shape)
print(mydf.columns)

print("the thing !!!!!!!!!!!!!!!!")
print("of a total of these ids that are *top* AND are *not in pathways*:")
print(len(mydf["moleculeIds"].unique()))
print("that represent a total of 'n' unique m/z:")
print(len(mydf["mz"].unique()))

print(" i could obtain this nb of ids: ")
the_intersection_ids = set(mydf["moleculeIds"].unique().tolist()).intersection(
    set(jjj['HMDB'].unique().tolist()))
print(len(the_intersection_ids))
print("and that represents a total 'n' of unique m/z:")
print(len(mydf.loc[mydf['moleculeIds'].isin(the_intersection_ids), "mz"].unique()))


uu = 9
t = 5
import sys
sys.exit()
