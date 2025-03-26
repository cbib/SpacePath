import sys
import os
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import zonist.data as zd
import zonist.method as zm
import zonist.helpers as zh
sys.path.append(os.path.join(os.path.dirname(__file__),
                             '../bank4zonist'))

sys.path.append(os.path.join(os.path.dirname(__file__),
                             '../'))

from user_routines import prepare_harmonized_objects

from quantifs_methods import *
from visualization import *
import scanpy as sc

"""
evaluate 
"""

SAMPLE =  "test635_200-1K"

print(SAMPLE)
mydata_dir = "/home/johanna/spatial_thesis/spatial_metabolomics/data"
#"/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"
subfolder = "gb-brain"
anno_dict = {'HMDB': f"{SAMPLE}_hmdb.csv",
             "KEGG": f"{SAMPLE}_kegg.csv",
             "CoreMetabolome": f"{SAMPLE}_CoreMetabolome.csv",
             "LipidMaps": f"{SAMPLE}_lipidmaps.csv"}
plot_folder = "../../../figures/gb-brain"
prob_out = 0.4
min_pathway_size = 1
assert os.path.exists(plot_folder), "Cannot continue, the plot_folder does not exist"

me_annots = zd.MetabolomeAnnotations(name=SAMPLE,
                            data_dir=mydata_dir,
                            subfolder=subfolder,
                            annotation_files=anno_dict)
me_annots.load()

adata_mz = ad.read_h5ad(filename=os.path.join(
    mydata_dir, subfolder, f'{SAMPLE}.h5ad'))

okmz = adata_mz.var.loc[ adata_mz.var['final_quality'] == "good",  "mz"].astype(str).tolist()

##
adata_mz = adata_mz[:, okmz]
print(adata_mz)

unmaps = adata_mz.var.copy()  # always declare before further subset

me_refmap = zd.MetaboliteIdsMappings(
        data_dir=mydata_dir,
        subfolder="metabolites-mappings",
        mapping_files={'HMDB': "HMDB-to-KEGG_05-07-24.tsv",
                       'LipidMaps': "LipidMaps-to-KEGG_10-18-24.tsv",
                       'CoreMetabolome': "CoreMetabolome-to-KEGG_05-07-24.tsv"})
me_refmap.load()


me_annots = update_annotation_tomatch_adatamz(me_annots, adata_mz)  # todo: mark as new
print(me_annots, "!!!!!!!!!!!!!!")
verification(adata_mz, me_annots)  # todo del
print("Excludes m/z with out-of-tissue probability >=", prob_out)
me_annots.discard_out_of_tissue_mz(prob_out=prob_out)
me_annots.convert_ids_to_kegg(me_refmap.table)
me_annots.treat_mz_multiplets(adata_mz)  # cleans .df and also creates property .mz_redundant

## pathways are also the same for all my projects
pathwaysObj = zd.MetabolicPathways(
    data_dir=mydata_dir,
    subfolder="metabolites-mappings",
    official_file="KEGG_hsa_pathways_compounds_R110.gmt")

pathwaysObj.load()
initial_pws_dict = pathwaysObj.to_dict()  # reserver pour plustard ? simplement avoir les tailles des sets du départ.
pathwaysObj, me_annots = zm.SyncPathwaysAndAnnots.run(pathwaysObj.model_copy(),
                                                          me_annots.model_copy(),
                                                          min_pathway_size=min_pathway_size,
                                                          max_pathway_size=200)

ll = [str(z) for z in sorted(me_annots.df["mz"].astype(float).unique())]
adata_mz = adata_mz[:, ll].copy()  # todo this is NOT sync, but like so

print("Saving things not mapping to any pathways")

print(f"when mapped to KEGG, with min pathwaysize {min_pathway_size}, "
      f"remained {me_annots.df.shape[0]} different KEGG id")

unmaps = unmaps.loc[~unmaps['mz'].astype(str).isin(ll),:]
unmaps = unmaps.loc[~unmaps["mz"].astype(str).isin(me_annots.mz_redundant), :]
unmaps = unmaps.drop_duplicates(subset='mz')
tmp = me_annots.preloaded[["mz", "moleculeIds", "moleculeNames","database"]].copy()
tmp['mzETdb'] = tmp['mz'].astype(str).str.cat(tmp['database'], sep="&")
tmp2 = me_annots.mz_quality.copy()
tmp2['mzETdb'] = tmp2['mz'].astype(str).str.cat(tmp2['database'], sep="&")
tmp = tmp.loc[tmp['mzETdb'].isin(tmp2['mzETdb'].unique()), :]
tmp = tmp.drop_duplicates(subset="mz")
tmp['mz'] = tmp["mz"].astype(str).tolist()
intersectionmzs = list(set(unmaps["mz"].astype(str).tolist()).intersection(
    set(tmp["mz"].astype(str).tolist())))
tmp = tmp.drop(columns='mzETdb')
tmp = tmp.set_index("mz")
unmaps = unmaps.loc[unmaps["mz"].astype(str).isin(intersectionmzs)]
tmp = tmp.loc[intersectionmzs,:]
for col in tmp.columns:
    unmaps[col] = tmp[col].tolist()
del tmp

unmaps = unmaps.sort_values(by=[ "I", "rawOffSampleProb"], ascending=[ False, True])
print(f"Total of unique mz not mapped at all: {unmaps.shape[0]}")
unmaps.to_csv(f"{SAMPLE}-notMappableAll.csv", index=False)

unmaps.loc[(unmaps["rawOffSampleProb"] <= 0.1
            ) & (unmaps["msm"] >= 0.1), :].to_csv(f"{SAMPLE}-notMappableTOP.csv", index=False)

# Pathways mapping
##   Total of unique mz not mapped at all: 144



"""
"""
