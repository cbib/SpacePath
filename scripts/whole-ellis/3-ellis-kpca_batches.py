import string
import sys
import os
import numpy as np
import time
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
from typing import Dict, List


import zonist.data as zd
import zonist.helpers as zh

# why do I have dask here ? TODO check if really used :
from dask.config import update

from scipy.sparse import csr_matrix
from sklearn.utils import sparsefuncs
import squidpy as sq
import networkx as nx
from tqdm import tqdm

import sspa as sspa


sys.path.append(os.path.join(os.path.dirname(__file__),
                             '../bank4zonist'))
from statsembryo1 import *
from vizembryo1 import *
from pathways_methods import *
from datasetcomplement import (build_kegg_resolved_adata,
                               yield_specific_typical_pathways)



## start
# only use random batch

SAMPLE = "whole-ellis"  # .h5ad saved by       -pixel-selection.py
mydata_dir =  "/home/johanna/spatial_thesis/spatial_metabolomics/data" #"/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"
subfolder = "whole-ellis"
anno_dict = {'HMDB': "whole-ellis_hmdb.csv",
             'LipidMaps': "whole-ellis_lipidmaps.csv"}
plot_folder = "../../../figures/ellis"

my_option = "RandomBatch"
my_n_grous = 30

assert os.path.exists(plot_folder), ("Cannot continue, the plot_folder does "
                                     "not exist")

me_annots = zd.MetabolomeAnnotations(name=SAMPLE,
                            data_dir=mydata_dir,
                            subfolder=subfolder,
                            annotation_files=anno_dict)
me_annots.load()


adata_mz = ad.read_h5ad(filename=os.path.join(
    mydata_dir, subfolder, f'{SAMPLE}_mx.h5ad'))


## reference mappings to get KEGG ids from HMDB and LipidMaps identifiers
me_refmap = zd.MetaboliteIdsMappings(
    data_dir=mydata_dir,
    subfolder="metabolites-mappings",
    mapping_files={'HMDB': "HMDB-to-KEGG_05-07-24.tsv",
                   'LipidMaps': "LipidMaps-to-KEGG_10-18-24.tsv" }
)
me_refmap.load()
me_annots.convert_ids_to_kegg(me_refmap.table)
me_annots.treat_mz_multiplets(adata_mz)


## pathways
pathwaysObj = zd.MetabolicPathways(
    data_dir=mydata_dir,
    subfolder="metabolites-mappings",
    official_file="KEGG_hsa_pathways_compounds_R110.gmt")

pathwaysObj.load()
pathwaysObj.intersect_with_annotations(me_annots)   # updated pathwaysObj

print("subset adata_mz for only mz mapping to a KEGG metabolite identifiers that in turn map to pathways")
me_annots_OKdf = me_annots.df.loc[
              me_annots.df["KEGG"].isin(pathwaysObj.piled["KEGG"].unique()), :]  # now update the annotations

adata_mz = adata_mz[:, me_annots_OKdf["mz"].astype(str).unique()].copy()

# the KEGG id adata
adata = build_kegg_resolved_adata(adata_mz, annots_df=me_annots_OKdf)
print("KEGG resolved adata: ", adata)

del adata_mz
import gc
gc.collect()


# exclude all metabolites that do not comply with the policy "at least 2 metabolites by pathway"
good_pathways, metabolites_to_throw_away = pathways_ok__unwanted_metabolites(
    pathwaysObj , min_pathway_size=2)
me_annots, pathwaysObj = update_really(me_annots, pathwaysObj, good_pathways, metabolites_to_throw_away)


pw_as_df = pd.DataFrame(
    {key: pd.Series(value) for key, value in pathwaysObj.pathways.items()})
pw_as_df = pw_as_df.transpose()
pw_as_df['Pathway_name'] = 0  # sspa expects this column, which is dropped

groups_d = dict()


    groups_d = spots_groups_by_random_choice(
            adata.obs['obs'].to_numpy(), n_groups=my_n_groups)

    t0 = time.time()
    kpca_scores_d = compute_kpca_on_batches(pw_as_df, adata, groups_d)
    hours, minutes, seconds = time_formatting(time.time() - t0)
    print(f"HH:MM:SS : {hours}:{minutes}:{seconds}")
    print("gathering results...")

    kpca_scores = verticalassembly(kpca_scores_d, ordered_spots_names=adata.obs['obs'].tolist())
    del kpca_scores_d
    kpca_score_adata = ad.AnnData(
        X=csr_matrix(kpca_scores),
        var=pd.DataFrame(
            {'pathway_id' : kpca_scores.columns.tolist()},
            index=kpca_scores.columns.tolist()),
        obs=adata.obs,
        obsm=adata.obsm,
        uns=adata.uns
    )

    kpca_scores.write_h5ad(os.path.join(mydata_dir, subfolder,
                                        f"{SAMPLE}-{my_option}{my_n_grous}-ScoresKPCA.h5ad"))

    print("saved ! ")
    biologically_chosen_pathways = yield_specific_typical_pathways()

    plot_selected_pathways(kpca_score_adata, biologically_chosen_pathways,
                      f"{SAMPLE}-{my_option}{my_n_grous}-ScoresKPCA",
                           plot_folder)










