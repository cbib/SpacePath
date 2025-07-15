import os
from typing import Dict, List, Union

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
from scipy.sparse import csr_matrix
from sklearn.cluster import KMeans
from sklearn.utils import sparsefuncs
from zonist.data import MetabolicPathways, MetabolomeAnnotations
import zonist.helpers as zh
import zonist.method as zm



def get_top_n(group, *args):
    ntop_per_cluster = args[0]
    sorted_group = group.sort_values(
        by=['pvals_adj', 'abslfc'], ascending=[True, False])
    return sorted_group.head(ntop_per_cluster)


def differential_results_filter(differential_results, lfc: float=0.8,
                                pvals_adj: float=0.05,
                                ntop_per_cluster: Union[str, int]="auto",
                                total_n_features:int=250):
    if ntop_per_cluster == "auto":
        ntop_per_cluster = int(total_n_features / len(
            differential_results['group'].unique()))

    tmp = differential_results.loc[(abs(differential_results[
             'logfoldchanges']) >= lfc) & (
            differential_results['pvals_adj'] <= pvals_adj), :].copy()

    tmp = tmp.assign(abslfc=abs(tmp['logfoldchanges']),
                     cluster=tmp['group'].astype(str))
    tmp = tmp.sort_values(by=['group', 'pvals_adj'], ascending=[True, False])
    tmp = tmp.groupby('cluster', group_keys=False, as_index=False
                      ).apply(get_top_n, ntop_per_cluster,
                              include_groups=False  # avoid DeprecationWarning
                              ).reset_index(drop=True)
    return tmp




def score_as_adata(kpca_scores, adata,
                   sets_ids_df_df: pd.DataFrame,
                   pathways_piled: pd.DataFrame):
    score_adata = ad.AnnData(
        X=csr_matrix(kpca_scores),
        var=pd.DataFrame(
            {'pathway_id': kpca_scores.columns.tolist()},
            index=kpca_scores.columns.tolist()),
        obs=adata.obs,
        obsm=adata.obsm
    )
    score_adata.uns['spatial'] = {'my_slice_name': {
        'scalefactors': {'fiducial_diameter_fullres': 100,
                         'spot_diameter_fullres': 1,
                         'tissue_hires_scalef': 0.5,
                         'tissue_lowres_scalef': 0.5}}}
    tmp = sets_ids_df_df
    tmp.index = tmp['pathway_id']
    tmp = tmp.loc[score_adata.var['pathway_id'],:]
    score_adata.var = score_adata.var.assign(
        pathway_name=tmp['pathway_name'].tolist())

    score_adata.uns['current_pathways_df'] = pathways_piled.copy() # # has m/z, KEGG/moleculeIds and pathway id columns

    return score_adata


def compute_mz_binary_df(adata_mz: ad.AnnData):
    X_presence_absence_mask = adata_mz.X.todense() > 0
    X_mz_binary_df = pd.DataFrame(
        np.array(X_presence_absence_mask, dtype=int),
        columns= adata_mz.var['mz'].to_numpy(),
        index=adata_mz.obs.index.tolist())

    return X_mz_binary_df


def compute_metabolites_probs_df(adata_mz: ad.AnnData, annots_df,
                                 usable_id, mz_col="mz"):
    X_mz_binary_df = compute_mz_binary_df(adata_mz)
    preX_probs = dict()
    for a_mz in X_mz_binary_df.columns:
        i_arr = X_mz_binary_df[a_mz].to_numpy()
        the_ids_here = annots_df.loc[
            annots_df[mz_col] == a_mz, usable_id].unique().tolist()
        n = len(the_ids_here)
        for kid in the_ids_here:
            registered_prob = 0.0
            if kid in preX_probs.keys(): # if already registered from another mz, as it can happen, checked
                registered_prob = preX_probs[kid][0]
            compound_probability = max(registered_prob, 1 / n)  # attribute the best probability between registered (if any) and current
            res = i_arr * np.repeat(compound_probability, len(i_arr))
            preX_probs[kid] = np.array(res, dtype=np.float32)

    X_metabolites_probs_df = pd.DataFrame(preX_probs)
    X_metabolites_probs_df.index = X_mz_binary_df.index
    return X_metabolites_probs_df


def build_identifiers_adata(adata_mz,  X_metabolites_probs_df,
                            annots_df, usable_id):
    mz_ids = annots_df[['mz', usable_id]].copy()
    mz_ids = mz_ids.drop_duplicates()
    mz_ids = mz_ids.sort_values(by=usable_id)
    unique_ids = mz_ids[usable_id].unique().tolist()
    tmp_X = dict()
    for u_id in unique_ids:
        mz_one_or_more = mz_ids.loc[mz_ids[usable_id]==u_id, 'mz'].astype(str).tolist()
        sub_adata = adata_mz[:,mz_one_or_more ].copy()
        mz_arr = np.array(sub_adata.X.todense().T).flatten()
        if len(mz_one_or_more) > 1:
            tmp = mz_arr.sum(axis=0)
            mz_arr = tmp

        curr_prod_arr = X_metabolites_probs_df[u_id].to_numpy() * mz_arr

        tmp_X[u_id] = np.array(curr_prod_arr, dtype=np.float32)

    identifiers_list = list(tmp_X.keys())
    tmp_X = np.array(pd.DataFrame(tmp_X), dtype=np.float32)
    adataIds = ad.AnnData(X=csr_matrix(tmp_X),
                          obs=adata_mz.obs,
                          var=pd.DataFrame({"ID": identifiers_list},
                                           index=identifiers_list),
                          obsm=adata_mz.obsm)
    return adataIds


def run_pathway_scores_steps(adata_mz, pws_obj, annots_df, usable_id):
    adata_mz = adata_mz[:, annots_df['mz'].astype(str).unique()].copy()
    X_df = compute_metabolites_probs_df(adata_mz, annots_df, usable_id)
    adataIds = build_identifiers_adata(adata_mz, X_df, annots_df, usable_id)
    del adata_mz # free some space for computation
    scores_df = zm.ComputePathwaysScores.run(
        pws_obj.to_dict(), adataIds
    )
    scores_adata = score_as_adata(scores_df, adataIds,
                   pws_obj.sets_ids_df,
                   pws_obj.piled)

    return scores_adata