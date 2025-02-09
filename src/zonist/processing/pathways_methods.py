import os
import logging
import numpy as np

from typing import List, Dict, Literal, Optional, Set, Union

import pandas as pd
import anndata as ad
import sspa as sspa

from pydantic import BaseModel as PydanticBaseModel
from pydantic import Extra
from .quantifications_methods import is_log_transformed


def test_function(arr, scale_val, threshold):
    scaled_array = arr * scale_val
    return np.where(scaled_array > threshold, 1, 0)


def spots_groups_by_random_choice(spots_names: List[str],
                                  n_groups: int=4) -> Dict[int, List[str]]:
    seen_list = []
    n_spots_by_group = int(len(spots_names) / n_groups)
    residual_n = len(spots_names) % n_groups
    groups_dict = dict()
    for i in range(0, (n_groups - 1)):
        remaining_spots_l = list(set(spots_names) - set(seen_list))
        sps = np.random.choice(
            np.array(remaining_spots_l), n_spots_by_group, replace=False)
        groups_dict[i] = list(sps)
        seen_list += list(sps)

    n_spots_last_group = n_spots_by_group + residual_n   # size of last group
    remaining_spots_l = list(set(spots_names) - set(seen_list))
    sps = np.random.choice(
        np.array(remaining_spots_l), n_spots_last_group, replace=False)
    groups_dict[n_groups - 1] = list(sps)

    return groups_dict


def compute_kpca(pw_as_df, adata, min_entity=2):
    print("Computing Pathway Scores with kpca from sspa package")
    kpca_model = sspa.sspa_kpca.sspa_KPCA(
        pw_as_df, min_entity=min_entity, random_state=0)

    if is_log_transformed(adata.X):
        df_log_transformed = adata.to_df()
    else:
        df_log_transformed = np.log(adata.to_df() + 1) # as sc.pp.log1p does

    kpca_scores = kpca_model.fit_transform(df_log_transformed)

    return kpca_scores


def reshape_pathways_for_kpca(pathways_dict):
    pw_as_df = pd.DataFrame({key: pd.Series(
        m_list) for key, m_list in pathways_dict.items()})
    pw_as_df = pw_as_df.transpose()
    pw_as_df['Pathway_name'] = 'none'  # sspa expects this column, which is dropped

    return pw_as_df


def compute_kpca_on_spot_groups(pw_as_df, adata, groups_d: Dict[int, list],
                                min_entity=2) -> Dict[int, dict]:
    print(f"Computing Pathway Scores for {len(groups_d)} groups of spots...")
    kpca_scores_d = dict()
    for k in groups_d.keys():
        print("Group:", str(k + 1))
        adata_here = adata[groups_d[k], :].copy()
        kpca_scores_d[k] = compute_kpca(pw_as_df, adata_here, min_entity)

    return kpca_scores_d


def gather_scores_dfs_dict(kpca_scores_d: Dict[int, pd.DataFrame],
                           ordered_spots_names: List[str]) -> pd.DataFrame:
    """return the dataframe of computed scores"""
    def add_missing_cols(df, expected_cols: set):
        absent_expected_cols = list(
            set(expected_cols) - set(df.columns.tolist()))
        for aec in absent_expected_cols:
            df[aec] = np.nan
        return df

    columns_set = set()
    for k in kpca_scores_d.keys():
        columns_set = columns_set.union(set(
            kpca_scores_d[k].columns.tolist()))
    # end for
    tmp = add_missing_cols(kpca_scores_d[0], columns_set)
    tmp = tmp[sorted(list(columns_set))]
    for k in range(1, len(kpca_scores_d.keys())):
        tmp_k = add_missing_cols(kpca_scores_d[k], columns_set)
        tmp_k = tmp_k[sorted(list(columns_set))]
        tmp  = pd.concat([tmp, tmp_k], axis=0)
    tmp = tmp.fillna(tmp.min().min())
    tmp = tmp.loc[ordered_spots_names, :]

    return tmp


def compute_pathways_scores(pathways_dict,
                            adata: ad.AnnData, method: str, min_entity: int,
                            n_groups_of_spots: int,
                            threshold_n_spots_is_high: int)-> pd.DataFrame:
    if method == "sspa_kpca":
        pw_as_df = reshape_pathways_for_kpca(pathways_dict)
        if adata.obs.shape[0] < threshold_n_spots_is_high:
            kpca_scores = compute_kpca(pw_as_df, adata, min_entity=min_entity)
        else:
            groups_d = spots_groups_by_random_choice(
                adata.obs['obs'].tolist(), n_groups=n_groups_of_spots)
            kpca_scores_d = compute_kpca_on_spot_groups(
                pw_as_df, adata, groups_d, min_entity=min_entity)
            kpca_scores = gather_scores_dfs_dict(
                kpca_scores_d, ordered_spots_names=adata.obs['obs'].tolist())
    else:
        raise NotImplementedError

    return kpca_scores


def find_metabolites_to_exclude(good_pathways: List[str], pathways_dict):
    """TODO: ddetailed doc, very important"""  # todo !
    mts = {k: v for k, v in pathways_dict.items() if k in good_pathways}
    good_compounds_list = [i for m_id_list in mts.values() for i in m_id_list]

    metabolites_to_eject = set()
    unwanted_pathways = list(
        set(list(pathways_dict.keys())) - set(good_pathways))
    for pw in unwanted_pathways:
        tmp_set = set(pathways_dict[pw]) - set(good_compounds_list)
        metabolites_to_eject = metabolites_to_eject.union(tmp_set)

    return list(set(metabolites_to_eject))


def yield_compliant_pathways(pathways_dict,
                             min_pathway_size=2,
                             max_pathway_size=500) -> List[str]:
    """TODO doc"""
    pathways_sizes = {k: len(pathways_dict[k]) for k in pathways_dict.keys()}
    good_pathways_sizes = {k: v for k, v in pathways_sizes.items()
                  if ((v >= min_pathway_size) and (v <= max_pathway_size)) }
    good_pathways: List[str] = list(good_pathways_sizes.keys())

    return good_pathways


