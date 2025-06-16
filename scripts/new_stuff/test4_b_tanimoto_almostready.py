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



smiles_lipidmaps = pd.read_csv(
    "../resources/LipidMaps-SMILES_06-08-25.tsv",
     sep='\t',
    header=0)

smiles_hmdb =  pd.read_csv(
    "../resources/HMDB-SMILES_05-07-24.tsv",
     sep='\t',
    header=0)

smiles_df = pd.concat([smiles_lipidmaps, smiles_hmdb])


#################
with open("../../data/TD-Lipido-635/TD-Lipido-635-WIKI-NO.txt", "r") as f:
    WIKI_NO = f.readlines()
WIKI_NO.pop(0)
WIKI_NO = [float(t.strip()) for t in WIKI_NO]

with open("../../data/TD-Lipido-635/TD-Lipido-635-KEGG-NO.txt", "r") as f:
    KEGG_NO = f.readlines()
KEGG_NO.pop(0)
KEGG_NO = [float(t.strip()) for t in KEGG_NO]

ALL_MZ_NO = list(set(WIKI_NO).union(set(KEGG_NO)))


mydata_dir = "../../data/"
SAMPLE = "TD-Lipido-635"
subfolder = SAMPLE

print("loading adata, and it already has quality in .var attr")
adata_mz = ad.read_h5ad(filename=os.path.join(
            mydata_dir, subfolder, f'{SAMPLE}.h5ad'))

print(adata_mz.uns["history"])
print("MANDATORY: take away bad quality m/z") # todo: how to make an alert ?
adata_mz = adata_mz[:,
           adata_mz.var.loc[adata_mz.var['final_quality'] == "good",
           "mz"].astype(str).tolist()]

anno_dict = provide_annotation_dict(SAMPLE,
                                    ["LipidMaps", "HMDB", "CoreMetabolome"])
me_annots = zd.MetabolomeAnnotations(
    name=SAMPLE, data_dir=mydata_dir,
    subfolder=subfolder, annotation_files=anno_dict)
me_annots.load()
me_annots.usable_id = "moleculeIds"
me_annots.verify()

print("updating me_annots to be as in adata_mz")
me_annots, adata_mz = zmm.match_adata_mz_and_annots(me_annots, adata_mz)

me_annots.keep_unique_mz_id_tuples(adata_mz.var)
me_annots.remove_redundant_aliases()  # todo highlight this as superimportant

print("re-updating me_annots to be as in adata_mz")  # todo check how to avoid to re-do this
me_annots, adata_mz = zmm.match_adata_mz_and_annots(me_annots, adata_mz)

print("intersect with the spatially relevant, top moranI, and save ")
spatially_interesting = adata_mz.var.copy()
spatially_interesting = spatially_interesting.loc[spatially_interesting['I'] >= 0.4,:].copy()

differential_results = pd.read_csv("differential_results.csv", index_col=None)

differential_results = differential_results.loc[differential_results['names'].astype(float).isin(
    spatially_interesting['mz'].astype(float)
),:]

differential_results = differential_results.loc[
    (abs(differential_results['logfoldchanges']) >= 0.8) & (differential_results['pvals_adj'] <=0.05), :
].copy()

for clus in differential_results['group'].unique():
    print("cluster", str(clus),"top featurs NOT WIKI / top features this cluster ")
    # TODO: how to recall that "names" in the result of rank_genes_groups is the same as mz
    the_elems_this_cluster = differential_results.loc[differential_results['group']==clus,"names"].unique().tolist()
    not_wiki_here = set(the_elems_this_cluster).intersection(set(ALL_MZ_NO))
    n_this_cluster = len(the_elems_this_cluster)
    print(len(not_wiki_here), "/", n_this_cluster)
    print("")

##

differential_results['names'] = differential_results['names'].astype(float)

# take those annots that are NOT in pathwyas:
mydf = me_annots.df.loc[me_annots.df['mz'].astype(float).isin(ALL_MZ_NO), :]

# and that are top markers:
if max(abs(differential_results['logfoldchanges'].to_numpy())) >= 0.8:
    mydf = mydf.loc[mydf['mz'].astype(float).isin(
        differential_results['names'].unique()),:]

print("of a total of these ids that are *top* AND are *not in pathways*:")
print(len(mydf["moleculeIds"].unique()))
print("that represent a total of 'n' unique m/z:")
print(len(mydf["mz"].unique()))


print(" i could obtain this nb of SMILES: ")
tmp_df = smiles_df.loc[~smiles_df['SMILES'].isnull(),:].copy()
the_intersection_ids = set(mydf["moleculeIds"].unique().tolist()).intersection(
    set(tmp_df['moleculeIds'].unique().tolist()))
print(len(the_intersection_ids))
print("and that represents a total 'n' of unique m/z:")
print(len(mydf.loc[mydf['moleculeIds'].isin(the_intersection_ids), "mz"].unique()))

print("on the other hand for CHEBI, I could obtain this number of non empty CHEBI ids")
tmp_df2 = smiles_df.loc[~smiles_df['CHEBI'].isnull(),:].copy()
the_intersection_ids2 = set(mydf["moleculeIds"].unique().tolist()).intersection(
    set(tmp_df2['moleculeIds'].unique().tolist()))
print(len(the_intersection_ids2))
print("and that represents a total 'n' of unique m/z:")
print(len(mydf.loc[mydf['moleculeIds'].isin(the_intersection_ids2), "mz"].unique()))

print("\n I have decided to use tanimoto with smiles, will export the annotations to compute in another script")
out_df = mydf.loc[mydf['moleculeIds'].isin(tmp_df['moleculeIds']),:]
out_df.to_csv('identifs_for_tanimoto.csv', index=False)
uu = 9
t = 5

