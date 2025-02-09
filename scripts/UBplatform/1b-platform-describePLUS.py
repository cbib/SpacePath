import sys
import os
import numpy as np

import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import zonist.data as zd
from dask.config import update
from scipy.sparse import csr_matrix
from sklearn.utils import sparsefuncs
import squidpy as sq
import networkx as nx
from tqdm import tqdm

import sspa as sspa

from sklearn.metrics.pairwise import cosine_similarity

#from scripts.bank4zonist.statsembryo1 import compute_mz_binary_df

sys.path.append(os.path.join(os.path.dirname(__file__),
                             '../bank4zonist'))
#
# from scripts.bank4zonist.statsembryo1 import (verification,
#                                               compounds_counts_by_pathway,
#                                               pathways_ok__unwanted_metabolites,
#                                               update_really, compute_mz_binary_df,
#                                               compute_metabolites_probs_df,
#                                               compute_pathways_pseudoscores,
#                                               ugly_long_function)
#
# from scripts.bank4zonist.datasetcomplement import (build_kegg_resolved_adata,
#                                                    yield_specific_typical_pathways)

from statsembryo1 import *
from vizembryo1 import *
from datasetcomplement import (build_kegg_resolved_adata,
                               yield_specific_typical_pathways,
                               reshape_pathways_for_kpca)



"""
pathways scores on *bordeaux proteome example data* 

Bash call:
mylist=(UBcalibr1 UBcalibr2 UBcalibr3 UBcalibr4 UBcalibr5);for k in ${mylist[*]}; do echo $k; python xxxxx.py $k; done

"""

SAMPLE = "UBcalibr2"
#"UBcalibr2"
print(SAMPLE)
mydata_dir = "/home/johanna/spatial_thesis/spatial_metabolomics/data"
#"/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"
subfolder = "UBplatform"
#"/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"
anno_dict = {'HMDB': f"{SAMPLE}_hmdb.csv",
             "LipidMaps": f"{SAMPLE}_lipidmaps.csv"}
plot_folder = "../../../figures/platform"


assert os.path.exists(plot_folder), "Cannot continue, the plot_folder does not exist"

me_annots = zd.MetabolomeAnnotations(name=SAMPLE,
                            data_dir=mydata_dir,
                            subfolder=subfolder,
                            annotation_files=anno_dict)
me_annots.load()


adata_mz = ad.read_h5ad(filename=os.path.join(
    mydata_dir, subfolder, f'{SAMPLE}_mx.h5ad'))
print(adata_mz)

# ## check with and without tese two lines:
# adata_mz = adata_mz_filter_negligible_abundance(adata_mz, spots_non_zero=5)  # this filter is soft: ~ 2 mz are excluded, no reason to use  # TODO delete

verification(adata_mz, me_annots)

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
initial_kegg_pathways = pathwaysObj.to_dict()  # if I want to keep the initial composition of Metabolic pws as downloaded from KEGG db

pathwaysObj.intersect_with_annotations(me_annots)   # updated pathwaysObj   #

print("subset adata_mz for only mz mapping to a KEGG metabolite identifiers that in turn map to pathways")
me_annots_OKdf = me_annots.df.loc[
              me_annots.df["KEGG"].isin(pathwaysObj.piled["KEGG"].unique()), :]  # now update the annotations

adata_mz = adata_mz[:, me_annots_OKdf["mz"].astype(str).unique()].copy()

# # NEw!
adata_mz = spatial_standard_computations(adata_mz, 12, plot_folder,
                                        prefix=f"{SAMPLE}-MZ")
print("mz that do have KEGG compound ids : ", adata_mz) # to know how many mz are left

X_mz_binary_df = compute_mz_binary_df(adata_mz)
X_metabolites_probs_df = compute_metabolites_probs_df(X_mz_binary_df,
                                                      me_annots_OKdf)
# the KEGG id adata
adata = build_kegg_resolved_adata(adata_mz, annots_df=me_annots_OKdf)
print("KEGG resolved adata: ", adata)
uf = compounds_counts_by_pathway(pathwaysObj)
uf.to_csv("uf1.csv")

# is_it_well_done(X_metabolites_probs_df, adata_mz, adata, me_annots_OKdf)

# exclude all metabolites that do not comply with the policy "at least 2 metabolites by pathway"
good_pathways, metabolites_to_throw_away = pathways_ok__unwanted_metabolites(
    pathwaysObj , min_pathway_size=2)
me_annots, pathwaysObj = update_really(me_annots, pathwaysObj, good_pathways, metabolites_to_throw_away)


refinement_mz_pws_dict = ugly_long_function(pathwaysObj.piled.copy())
for t in refinement_mz_pws_dict.keys():
    print(t)
    print("   ", refinement_mz_pws_dict[t])
# Note: the refinement_mz_pws_dict has all sense and anti-sense pathways comparisons (pw1 vs pw2) and (pw2 vs pw1)
the_mz_for_refinement_list = [float(i) for i in list(refinement_mz_pws_dict.keys())]

adata = adata[:, adata.var.loc[adata.var[
    "ID"].isin(me_annots.df["KEGG"].unique()),"ID"].tolist()].copy()
print(adata)

uf = compounds_counts_by_pathway(pathwaysObj)
uf.to_csv("uf2.csv")

# for the moment use pseudoscores, TODO : in a separate script do kpca by grouping
pathways_pseudoscores = compute_pathways_pseudoscores(adata, pathwaysObj,
                                                      min_size=2,
                                                      method="sum")  # sum is preferred, average is also implemented

X_mz_binary_df = compute_mz_binary_df(adata_mz)
X_metabolites_probs_df = compute_metabolites_probs_df(X_mz_binary_df,
                                                      me_annots_OKdf)

########
##  new stuff

new_probabs = dict()

spot1 = "x72_y45"
spot2 = "x71_y45"
tmp_mx1 = adata_mz[spot1, :].X.todense()  # already 2D array as required by cosine_similarity
tmp_mx2 = adata_mz[spot2, :].X.todense()
local_cosine_sim = cosine_similarity(np.array(tmp_mx1),
                                             np.array(tmp_mx2))

def subset_m_probs_curr_pw(metabolites_probs_here, pathwaysObj, pw):  # TODO is this used ?

    theor_compounds_in_pw = pathwaysObj.to_dict()[pw]
    real_compounds_in_pw = set(theor_compounds_in_pw).intersection(
        set(metabolites_probs_here.columns.tolist()))

    return metabolites_probs_here[list(real_compounds_in_pw)].T

# if pathways_pseudoscores.shape[0] < 5000:
#     import seaborn as sns
#     sns.displot(x=np.array(pathways_pseudoscores).reshape(
#         pathways_pseudoscores.size))
#     plt.show()

# new_probabs[spot1] = dict()
#
# #X_metabolites_probs_df
# for one_mz in the_mz_for_refinement_list:
#     pws_tuples_l = list(refinement_mz_pws_dict[str(one_mz)].keys())
#     for pws_tuple in pws_tuples_l:
#         pw1 = pws_tuple[0] ;      pw2 = pws_tuple[1]   # todo: what if I say that I choose the second element as pw1 ? (inverse sense?)
#         metabolite1 = refinement_mz_pws_dict[str(one_mz)][pws_tuple][pw1][0]
#         metabolite2 = refinement_mz_pws_dict[str(one_mz)][pws_tuple][pw2][0]
#         pw_score = pathways_pseudoscores.loc[[spot1, spot2],  [pw1, pw2]]
#         met_prob = X_metabolites_probs_df.loc[[spot1, spot2],  [metabolite1, metabolite2]]
#         # green_light_pw1 = False
#         # green_light_pw2 = False
#
#         if (met_prob.loc[spot1, metabolite1] + met_prob.loc[
#             spot1, metabolite2] > 0):  # first check if probs are tunable: must be distinct than zero both
#                 happy1, happy2 = complicated_stuff_1(pw_score, met_prob, spot1, spot2, pw1, pw2, metabolite1, metabolite2)
#                 if happy1 + happy2 > 0:
#                     new_probabs[spot1] = { (one_mz,pws_tuple): {metabolite1: (met_prob.loc[spot1, metabolite1], happy1),
#                                       metabolite2: (met_prob.loc[spot1, metabolite1], happy2)}}
#
#                 # we do not touch spot 2

# spatial info
sq.gr.spatial_neighbors(adata, radius=1, coord_type="generic")
connectivities_sparse_mx = adata.obsp["spatial_connectivities"]
#dists = adata.obsp["spatial_distances"]

coordinates_names_str = adata_mz.obs["obs"].tolist()
coords_matcher_dict = {k: v for k, v in zip(
    range(0, len(coordinates_names_str)), coordinates_names_str)}
coords_matcher_rev = {v: k for k, v in coords_matcher_dict.items()}

print("ok here")

COSINE_SIMILARITY_THRESHOLD = 0.5
G = nx.from_scipy_sparse_array(connectivities_sparse_mx)
nodes_list = G.nodes()

new_probabs = dict()



for v1 in tqdm(nodes_list):
    neighs_of_v1 = G.neighbors(v1)
    tmp_mx1 = adata_mz[coords_matcher_dict[v1], :].X.todense()
    spot1 = coords_matcher_dict[v1]
    for v2 in neighs_of_v1:
        tmp_mx2 = adata_mz[coords_matcher_dict[v2], :].X.todense()
        local_cosine_sim = cosine_similarity(np.array(tmp_mx1),
                                             np.array(tmp_mx2))
        local_cosine_sim = local_cosine_sim.reshape(
            local_cosine_sim.shape[1])[0]
        if local_cosine_sim > COSINE_SIMILARITY_THRESHOLD:
            spot2 = coords_matcher_dict[v2]
            lala = complicated_stuff_2(refinement_mz_pws_dict, pathways_pseudoscores,
                    X_metabolites_probs_df, spot1, spot2)
            if len(lala) != 0 :
                new_probabs[spot1] = lala

print(new_probabs)
##
print(pleasestop)
#
# COSINE_SIMILARITY_THRESHOLD = 0.85
# G = nx.from_scipy_sparse_array(connectivities_sparse_mx)
# nodes_list = G.nodes()
# neighs_of10 = G.neighbors(10)
# print(coords_matcher_dict[10])
# for ki in tqdm(neighs_of10):
#     print(coords_matcher_dict[ki])
#
#
# local_cosine_sim = local_cosine_sim.reshape(
#             local_cosine_sim.shape[1])[0]
#
#
# spatial_and_cos_sim_G = nx.Graph()
# for v1 in tqdm(nodes_list):
#     neighs_of_v1 = G.neighbors(v1)
#     tmp_mx1 = adata_mz[coords_matcher_dict[v1], :].X.todense()  # already 2D array as required by cosine_similarity
#     # v1_expr_arr = np.array(tmp_mx1).reshape(tmp_mx1.shape[1])
#     for v2 in neighs_of_v1:
#         tmp_mx2 = adata_mz[coords_matcher_dict[v2], :].X.todense()
#         # v2_expr_arr = np.array(tmp_mx2).reshape(tmp_mx2.shape[1])
#         local_cosine_sim = cosine_similarity(np.array(tmp_mx1),
#                                              np.array(tmp_mx2))
#         local_cosine_sim = local_cosine_sim.reshape(
#             local_cosine_sim.shape[1])[0]
#         if local_cosine_sim > COSINE_SIMILARITY_THRESHOLD:
#             spatial_and_cos_sim_G.add_edge(v1, v2, weight=local_cosine_sim)
#             #print(coords_matcher_dict[v1], coords_matcher_dict[v2])
#
#
# new_nodes_list = spatial_and_cos_sim_G.nodes()
# for v1 in new_nodes_list:
#     neighs_of_v1 =  spatial_and_cos_sim_G.neighbors(v1)
#     spot1_name = coords_matcher_dict[v1]
#     for v2 in neighs_of_v1:
#         spot2_name = coords_matcher_dict[v2]



"""
barre_pw1 = 10 # min(pw_score.loc[spot1, pw1] / 1.2,  pw_score.loc[spot2, pw1] / 1.2)
            dummy_ammount = 100  # bad threshold TODO fix this  because values have wide range ! 0 - 60K
            if max(pw_score.loc[spot1, pw1],
                   pw_score.loc[spot2, pw1]) > dummy_ammount:  # pw1 in at least one spot is strong
                if abs(pw_score.loc[spot1, pw1] - pw_score.loc[spot2, pw1]) <= barre_pw1: # and both have pw1 similarly strong
                    green_light_pw1 = True
            # we want to 'steal' weight from pw2, but first we have to see if pw2 is weaker to pw1 to be able to do so, in spot 1
            if pw_score.loc[spot1, pw2] < pw_score.loc[spot1, pw1]:
                green_light_pw2 = True

            if green_light_pw1 and green_light_pw2:
                gain_or_pain = 0.1
                if met_prob.loc[spot1, metabolite2] > 0 :
                    happy2 = float(0)
                    if met_prob.loc[spot1, metabolite2] >= gain_or_pain:
                        happy2 = met_prob.loc[spot1, metabolite2] - gain_or_pain
                    else:
                        gain_or_pain = met_prob.loc[spot1, metabolite2]
                    happy1 = met_prob.loc[spot1, metabolite1] + gain_or_pain
"""

"""
bandita = max( met_prob.loc[spot1, metabolite1] / 2 ,
           met_prob.loc[spot1, metabolite2] / 2)
if (abs(met_prob.loc[spot1, metabolite1] - pw_score.loc[spot1, metabolite2]
    ) > bandita):
 pass
 uoi = max( pw_score.loc[spot1, pw1] / 2 ,
           pw_score.loc[spot2, pw1] / 2)
 if pw_score.loc[spot1, pw1] > 
"""


