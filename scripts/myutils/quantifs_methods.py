import os
from typing import Dict, List, Union

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix
from sklearn.cluster import KMeans
from sklearn.utils import sparsefuncs
from smpath.data import MetabolicPathways, MetabolomeAnnotations
import smpath.method as spm
from smpath.method import ComputePathwaysScoresConfig


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





# todo: important to reserve for documnetation:
# note: id do not build anymore anndata with the assignments, but only dataframe
# I can add to doc an example for puting it into adataIds to do plots and so on
#   ewample:
#      adataIds = ad.AnnData(X=csr_matrix(data_df),
#                       obs=adata_mz.obs,
#                       var=pd.DataFrame({"ID": data_df.columns.tolist()},
#                                        index=data_df.columns.tolist()),
#                       obsm=adata_mz.obsm)
#        adataIds.uns['spatial'] = adata_mz.uns['spatial']




def score_as_adata(scores_df,
                   sets_ids_df_df: pd.DataFrame,
                   pathways_piled: pd.DataFrame, coords: np.array):
    score_adata = ad.AnnData(
        X=csr_matrix(scores_df),
        var=pd.DataFrame(
            {'pathway_id': scores_df.columns.tolist()},
            index=scores_df.columns.tolist()),
        obs=pd.DataFrame({'obs': scores_df.index.tolist()},
                         index=scores_df.index.tolist()),
        obsm=coords
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


def run_pathway_scores_steps_produce_anndata(adataIds, coords,
                                             pws_obj,
                                             min_set_size,
                                             type_abbr):

    updated_config = spm.ComputePathwaysScoresConfig(type_abbr=type_abbr,
                                                     min_set_size=min_set_size)
    scores_df = spm.ComputePathwaysScores.run(
        pws_obj.to_dict(), adataIds, config=updated_config
    )
    scores_adata = score_as_adata(scores_df,
                   pws_obj.sets_ids_df,
                   pws_obj.piled,
                   coords)

    return scores_adata
