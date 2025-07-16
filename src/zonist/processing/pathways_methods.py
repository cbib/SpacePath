import sys
import os
import numpy as np

from typing import List, Dict, Literal, Optional, Set, Union

import pandas as pd
import anndata as ad
import sspa as sspa
import scanpy as sc

from pydantic import BaseModel as PydanticBaseModel
from pydantic import Extra
import zipfile
from zonist.helpers import is_log_transformed



def save_families_as_mmt_file(dedupl_df, molecules_ids_df,
                               SAMPLE, subfolder, mydata_dir) -> dict:
    molecs_families_dict = {}  # initialize and save as mmt
    for cluster in sorted(dedupl_df['molecules_cluster'].unique().tolist()):
        mols_in_cluster = dedupl_df.loc[
            dedupl_df['molecules_cluster'] == cluster,
            "moleculeIds"].unique().tolist()
        first_molecule_of_group = molecules_ids_df.loc[
            molecules_ids_df['moleculeIds'] == mols_in_cluster[
                0], "moleculeNames"].values[0]  # works better than .item()

        molecs_families_dict[
            f'MolFam_{cluster}'] = [first_molecule_of_group] + mols_in_cluster

    max_length = max(len(v) for v in molecs_families_dict.values())

    for key in molecs_families_dict:
        if len(molecs_families_dict[key]) < max_length:
            molecs_families_dict[key] += [np.nan] * (
                    max_length - len(molecs_families_dict[key]))

    molecs_families_mmt_df = pd.DataFrame(molecs_families_dict)
    molecs_families_mmt_df = molecs_families_mmt_df.T
    out_file = f"{SAMPLE}-molecules_families.mmt"
    molecs_families_mmt_df.to_csv(os.path.join(
        mydata_dir, subfolder, out_file), sep='\t', header=False)

    return molecs_families_dict

def compounds_and_smiles_ordered(smiles_df, molecules_ids_list):
    mycopy = smiles_df.copy()
    mycopy = mycopy.drop_duplicates(subset=['moleculeIds', 'SMILES'],
                                    keep='first')
    the_identifiers_allowed = list(
        set(mycopy['moleculeIds'].unique().tolist()).intersection(
            set(molecules_ids_list))
    )
    mycopy = mycopy[['moleculeIds', 'SMILES']].copy()
    mycopy.index = mycopy['moleculeIds'].tolist()
    mycopy = mycopy.T
    mycopy = mycopy[list(set(the_identifiers_allowed))]
    mycopy = mycopy.T
    mycopy = mycopy.drop_duplicates()

    return mycopy['moleculeIds'].tolist(), mycopy['SMILES'].tolist()


def unzip_if_not_tsv(library_path, smiles_files)-> None:
    smiles_folder = os.path.join(library_path, "resources", "SMILES")
    for type_smile_file in smiles_files.keys():
        smile_file = smiles_files[type_smile_file]
        f_basename, f_ext = os.path.splitext(smile_file)
        assert f_ext == ".tsv", "Error, smiles must be .tsv files"
        tsv_file_path = os.path.join(smiles_folder, f'{f_basename}{f_ext}')
        zip_file_path = os.path.join(smiles_folder, f'{f_basename}.zip')
        if not os.path.exists(tsv_file_path):
            with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
                zip_ref.extractall(smiles_folder)


def spots_groups_by_random_choice(spots_names: List[str],
                                  n_groups: int=4) -> Dict[int, List[str]]:
    np.random.seed(0)
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


def compute_kpca(pw_as_df: pd.DataFrame,
                 adata: ad.AnnData, min_entity: int=2)-> pd.DataFrame:
    """
    Computes kPCA on the compounds quantifications matrix
     to yield scores per pathway (or set), per spot.
    Credits of the sspa package: C Wieder et al., 2022

    :param pw_as_df: pathways or sets df, each row is a set (see sspa)
    :param adata: the data with the quantifications in its .X property
    :parm min_entity: the minimal size of the set to do the computation

    :returns: Computed scores as a dataframe
    """
    kpca_model = sspa.sspa_kpca.sspa_KPCA(
        pw_as_df, min_entity=min_entity, random_state=0)

    kpca_scores = kpca_model.fit_transform(adata.to_df())

    return kpca_scores


def reshape_pathways_for_kpca(pathways_dict):
    pw_as_df = pd.DataFrame({key: pd.Series(
        m_list) for key, m_list in pathways_dict.items()})
    pw_as_df = pw_as_df.transpose()
    pw_as_df['Pathway_name'] = 'none'  # sspa expects this column, which is dropped

    return pw_as_df


def compute_kpca_on_spot_groups(pw_as_df, adata, groups_d: Dict[int, list],
                                min_entity=2) -> Dict[int, dict]:
    kpca_scores_d = dict()
    for k in groups_d.keys():
        sys.stdout.write(".")
        sys.stdout.flush()
        adata_here = adata[groups_d[k], :].copy()
        kpca_scores_d[k] = compute_kpca(pw_as_df, adata_here, min_entity)

    return kpca_scores_d


def gather_scores_dfs_dict(kpca_scores_d: Dict[int, pd.DataFrame],
                           ordered_spots_names: List[str]) -> pd.DataFrame:
    """return the dataframe of gathered computed scores"""
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

        if not is_log_transformed(adata.X):
            sc.pp.log1p(adata)

        if adata.obs.shape[0] < threshold_n_spots_is_high:
            scores_res = compute_kpca(pw_as_df, adata, min_entity=min_entity)
        else:
            groups_dict = spots_groups_by_random_choice(
                adata.obs['obs'].tolist(), n_groups=n_groups_of_spots)
            scores__dict = compute_kpca_on_spot_groups(
                pw_as_df, adata, groups_dict, min_entity=min_entity)
            scores_res = gather_scores_dfs_dict(
                scores__dict, ordered_spots_names=adata.obs['obs'].tolist())
    else:
        raise NotImplementedError

    return scores_res


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
                             max_pathway_size=200) -> List[str]:
    """TODO doc"""
    pathways_sizes = {k: len(pathways_dict[k]) for k in pathways_dict.keys()}
    if min_pathway_size != max_pathway_size:
        good_pathways_sizes = {k: v for k, v in pathways_sizes.items()
                  if ((v >= min_pathway_size) and (v <= max_pathway_size)) }
    else:
        good_pathways_sizes = {k: v for k, v in pathways_sizes.items()
                               if v == min_pathway_size}
    good_pathways: List[str] = list(good_pathways_sizes.keys())

    return good_pathways
