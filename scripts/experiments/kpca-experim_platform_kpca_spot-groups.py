import sys
import os
import numpy as np

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

from sklearn.metrics.pairwise import cosine_similarity


sys.path.append(os.path.join(os.path.dirname(__file__),
                             '../bank4zonist'))

from quantifs_methods import *
from visualization import *
from pathways_methods import *  # kpca ....
from datasetcomplement import (build_kegg_resolved_adata,
                               yield_specific_typical_pathways)

"""
pathways scores on *bordeaux proteome example data* 

Bash call:
mylist=(UBcalibr1 UBcalibr2 UBcalibr3 UBcalibr4 UBcalibr5);for k in ${mylist[*]}; do echo $k; python xxxxx.py $k; done

"""


## start
CHOICE_STYLE_BATCH =  ["RandomGrps"] # ["RandomGrps", "RectGrps"]
SAMPLE = "UBcalibr2" #"UBcalibr2"

print(SAMPLE)
mydata_dir = "/home/johanna/spatial_thesis/spatial_metabolomics/data"
#"/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"
subfolder = "UBplatform"
#"/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"
anno_dict = {'HMDB': f"{SAMPLE}_hmdb.csv",
             "LipidMaps": f"{SAMPLE}_lipidmaps.csv"}
plot_folder = "../../../figures/platform"


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


# exclude all metabolites that do not comply with the policy "at least 2 metabolites by pathway"
good_pathways, metabolites_to_throw_away = pathways_ok__unwanted_metabolites(
    pathwaysObj , min_pathway_size=2)

hyhy = pathwaysObj.pathways_ids.copy().astype(str)  # todo delete

me_annots, pathwaysObj = update_really(me_annots, pathwaysObj, good_pathways, metabolites_to_throw_away)

hoho = pd.DataFrame({'pathway_id': list(pathwaysObj.to_dict().keys())}) # todo delete
yes = pd.merge(hoho, hyhy, how='left',
              left_on="pathway_id", right_on="pathway_id",
             left_index=False, right_index=False)

pw_as_df = reshape_pathways_for_kpca(pathwaysObj)

for my_option in CHOICE_STYLE_BATCH:
    if my_option == "RectGrps":
        kpca_scores =  pathways_scores_by_rectangular_cuts(pathwaysObj,
                                                           adata, n_groups=4)
    elif my_option == "RandomGrps":
        kpca_scores = compute_pathways_scores(pathwaysObj,
                                adata, method="sspa_kpca", min_entity=2,
                                n_groups_of_spots=4)

    kpca_score_adata = ad.AnnData(
        X=csr_matrix(kpca_scores),
        var=pd.DataFrame(
            {'pathway_id' : kpca_scores.columns.tolist()},
            index=kpca_scores.columns.tolist()),
        obs=adata.obs,
        obsm=adata.obsm,
        uns=adata.uns
    )

    hoho = kpca_score_adata.var.copy().astype(str)
    hyhy = pathwaysObj.pathways_ids.copy().astype(str)

    yes = pd.merge(hoho, hyhy, how='left',
                   left_on="pathway_id", right_on="pathway_id",
                   left_index=False, right_index=False)
    yes.to_csv("yeeees.txt")
    biologically_chosen_pathways = yield_specific_typical_pathways()

    plot_selected_pathways(kpca_score_adata, biologically_chosen_pathways,
                      f"{SAMPLE}-{my_option}--Picked", plot_folder)
