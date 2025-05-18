import sys
import logging
from typing import List, Dict, Union
import pandas as pd
import numpy as np
import json
import anndata as ad
#import scanpy as sc
#import squidpy as sq
from tqdm import tqdm
from sklearn.utils import sparsefuncs


def reset_log_config(logger):
    """parameters for the log file"""
    out_log_file_name = "log.log"
    logging.basicConfig(encoding='utf-8', stream=sys.stdout)
    file_handler = logging.FileHandler(out_log_file_name)
    logger.addHandler(file_handler)
    logger.setLevel(logging.DEBUG)
    return logger


def combine_frames_identical_columns(dfs_dict : Dict[str, pd.DataFrame]
) -> pd.DataFrame:
    list_of_frames = []
    for k in dfs_dict.keys():
        list_of_frames.append(dfs_dict[k])

    df_out = pd.concat(list_of_frames)
    df_out = df_out.drop_duplicates()

    return df_out


def spreadout_annotation_df(df: pd.DataFrame, mz_col: str, formula_col: str, 
    molecule_ids_col: str, molecule_names_col: str):
    """
    Transforms this format:

    mz      formula      moleculeIds                  moleculeNames                 database
    720.55	C39H78NO8P   LMGP02011252, LMGP01010532	  PE(13:0/21:0), PC(15:0/16:0)  LipidMaps
    
    Into:

    moleculeIds   mz            formula     moleculeNames   database
    LMGP02011252	720.553782	C39H78NO8P	PE(13:0/21:0)	LipidMaps
    LMGP01010532	720.553782	C39H78NO8P	PC(15:0/16:0)	LipidMaps
    """
    pd.options.mode.chained_assignment = None  # mute SettingWithCopyWarning

    delimiter_colon_plus_space = ", " # please do not remove the space

    tmp1 = df[[mz_col, formula_col, molecule_ids_col]]
    tmp2 = df[[molecule_names_col, 'database']]
    tmp1[molecule_ids_col] = tmp1[molecule_ids_col
                   ].str.split(delimiter_colon_plus_space)
    tmp2[molecule_names_col] = tmp2[molecule_names_col
                   ].str.split(delimiter_colon_plus_space)

    tmp1 = tmp1.explode(molecule_ids_col)
    tmp2 = tmp2.explode(molecule_names_col)
    df_out = pd.concat([tmp1, tmp2], axis=1)

    return df_out


def pick_mz_quality_infos(df, mz_col: str, offsample_col: str,
                         total_intens_col: str, fdr_column: str,
                         msm_column: str,  formula_col: str) -> pd.DataFrame:
    """ """
    df = df[[mz_col, offsample_col, total_intens_col, fdr_column,
             msm_column,  formula_col]]
    df = df.sort_values(by=[mz_col, offsample_col, fdr_column, msm_column],
                        ascending=[True, True, True, False])
    df = df.drop_duplicates(subset=[mz_col], keep='first')

    return df


def convert_metabo_annots_tokegg(
        annotdf: pd.DataFrame,  mappings_ref_ids: pd.DataFrame,
        molecule_ids_col: str, mz_col: str, usable_id: str="KEGG") -> pd.DataFrame:
    assert usable_id == 'KEGG', 'Error, usable_id must be "KEGG". Aborted.'
    dfs_list = []
    common_notkegg_database = set(annotdf['database'].unique()).intersection(
        set(mappings_ref_ids['database'].unique()))
    for cd in list(common_notkegg_database):
        tmp_df = annotdf.loc[annotdf['database'] == cd, :]
        tmp_ref = mappings_ref_ids.loc[mappings_ref_ids['database'] == cd, :]
        joined = pd.merge(tmp_df, tmp_ref, left_on=molecule_ids_col,
                          right_on='initial_id_unconverted', how='left')
        joined = joined.assign(database=joined.database_y)
        joined = joined.drop(columns=['database_x', 'database_y',
                                      'initial_id_unconverted']).dropna()
        dfs_list.append(joined)
    # end for
    if usable_id in annotdf['database'].unique():
        kegg_tmp_df = annotdf.loc[annotdf['database'] == usable_id, :]
        kegg_tmp_df = kegg_tmp_df.assign(KEGG=kegg_tmp_df[molecule_ids_col])
        dfs_list.append(kegg_tmp_df)

    final = pd.concat(dfs_list)
    final = final.sort_values(by=molecule_ids_col, ascending=True)
    final = final.drop_duplicates(subset=[mz_col, 'KEGG'], keep='first')
    final = final.sort_values(by=mz_col, ascending=True)

    return final


def identify_multi_mz_to_one_compound(annotdf: pd.DataFrame,
                                      mz_col: str, id_col: str) -> bool:
    """True if multiple mz for one Compound id exist in the annotations df"""
    out_bool = False
    df2 = annotdf[[mz_col, id_col]]
    df2 = df2.drop_duplicates()
    df2 = df2.groupby(id_col, as_index=False).size()
    df2 = df2.loc[df2['size'] >= 2, :].head()
    if df2.shape[0] > 0:
        out_bool = True

    return out_bool


def compute_compounds_n_per_mz(var_df, ann_df, mz_col, id_col):
    nn = list()
    for i, row in var_df.iterrows():
        mzfloat = float(row[mz_col])
        ids_l = ann_df.loc[ann_df[mz_col].astype(float) == mzfloat, id_col].unique()
        nn.append(len(ids_l))

    var_df['compound_ids_n__per_mz'] = nn

    return var_df


def get_approved_mz_id_pairs(ann_df: pd.DataFrame, var_df: pd.DataFrame,
                             mz_col: str, id_col: str, msm_col: str, fdr_col:str, offsam: str) -> List[tuple]:
    """
    Provides a final list of approved (m/z, compoundID) tuples.
    Note:
        In our piled table of annotations, an annotation error can yield
        (multiple m/z, one compound) pairs or tuples. E.g:
            (188.4545,  Fucose), (559.869, Fucose), (345.8556, Fucose).
        in such case, retains only the best pairs by quality criteria.
    output: List[tuple] a list of tuples: (mz (float), id (str))
    """
    ann_df = ann_df.loc[ann_df[mz_col].astype(float).isin(var_df[mz_col].astype(float)), :]
    df2 = ann_df.groupby(id_col, as_index=False).size()
    unique_compound_ids = df2.loc[df2['size'] == 1, id_col].tolist()
    ok_mz_id_tuples = [tuple(t) for t in ann_df.loc[
        ann_df[id_col].isin(unique_compound_ids), [mz_col, id_col]].to_numpy()]

    var_df = compute_compounds_n_per_mz(var_df, ann_df, mz_col, id_col)
    compounds_multiassigned = df2.loc[df2['size'] >= 2, id_col].unique()
    var_df = var_df.loc[var_df['compound_ids_n__per_mz'] >= 1, :]
    for one_id in compounds_multiassigned:
        multi_annot_df = ann_df.loc[ann_df[id_col] == one_id, :]
        curr_qual = var_df.loc[var_df[mz_col].astype(float).isin(
            multi_annot_df[mz_col].astype(float)), :].copy()
        if 'I' in curr_qual.columns :
            curr_qual = curr_qual.sort_values(
                by=['compound_ids_n__per_mz', msm_col, fdr_col, offsam, 'I'], # always respect this order !
                ascending=[True, False, True, True, False])
        else:
            curr_qual = curr_qual.sort_values(
                by=['compound_ids_n__per_mz', msm_col, fdr_col, offsam], # in this order !
                ascending=[True, False, True, True])
        its_best_mz = float(curr_qual.iloc[0][mz_col])
        ok_mz_id_tuples.append((its_best_mz, one_id))

    return ok_mz_id_tuples


def DEPRECATED_multimz_oneid(adata: ad.AnnData, df: pd.DataFrame,  # todo: delete this function
                             mz_col: str, id_col: str, criterion: str) -> List[tuple]:
    """
    In our piled table of annotations, an annotation error can yield
    multiple m/z -> one compound
       E.g: 188.4545 -> Fucose; 559.869 -> Fucose; 3345.8556 -> Fucose.
    To resolve such errors, this function keeps the
    m/z (of such compound) that has the highest variance (or Moran's I)
    of the quantifications.
    output: List[tuple] a list of tuples: (mz (float), kegg (str))
    """
    if criterion == "variance":
        mean_arr, v_arr = sparsefuncs.mean_variance_axis(adata.X, axis=0)
        adata.var = adata.var.assign(computed_variance = v_arr)
    elif criterion == "Moran":
        if "moranI" not in list(adata.uns.keys()):
            raise AssertionError("Could not find 'Moran' in adata.uns keys")
    df2 = df.groupby(id_col, as_index=False).size()
    unique_compound_ids = df2.loc[df2['size'] == 1, id_col].tolist()
    ok_mz_id_tuples = [tuple(t) for t in df.loc[
        df[id_col].isin(unique_compound_ids), [mz_col, id_col]].to_numpy()]
    repeated_compound_ids = df2.loc[df2['size'] >= 2, id_col].unique()

    for one_id in repeated_compound_ids:
        multi_annot_df = df.loc[df[id_col] == one_id, :]
        if criterion == "variance":
            df_with_v_metric = adata.var.loc[
                   adata.var[mz_col].isin(multi_annot_df[mz_col].tolist()),:]
            df_with_v_metric = df_with_v_metric.sort_values(
                by='computed_variance', ascending=False)
        elif criterion == "Moran":
            df_with_v_metric = adata.uns['moranI'].copy()
            df_with_v_metric[mz_col] = df_with_v_metric.index.astype(float).tolist()
            df_with_v_metric = df_with_v_metric.loc[
                df_with_v_metric[mz_col].isin(multi_annot_df[mz_col].tolist()), :]
            df_with_v_metric = df_with_v_metric.sort_values(by='I', ascending=False)

        mz_best = df_with_v_metric[mz_col].tolist()[0] # first item mz keep
        ok_mz_id_tuples.append([mz_best, one_id])
        # end for
    return ok_mz_id_tuples


def warrant__one_to_one_mz_id(annotdf: pd.DataFrame,
                             var_df, #
                              mz_col: str, id_col: str, msm_col: str, fdr_col:str, offsam: str
                              ) -> pd.DataFrame:
    ok_mz_kegg_tuples = get_approved_mz_id_pairs(
        annotdf, var_df, mz_col, id_col, msm_col, fdr_col, offsam)

    annotdf["keep"] = 0
    for tup in ok_mz_kegg_tuples:
        mz_i = tup[0]
        keggid = tup[1]
        annotdf.loc[(annotdf[mz_col] == mz_i) & (
                    annotdf[id_col] == keggid), "keep"] = 1

    rejected = annotdf.loc[annotdf["keep"] == 0, :].copy()
    annotdf = annotdf.loc[annotdf["keep"] == 1, :].copy()
    annotdf = annotdf.drop(columns=["keep"])

    return annotdf


def prettify_gmt(gmt_df, undesired_name_substring: Union[str, None]):
    gmt_df[0] = gmt_df[0].str.strip()
    gmt_df[1] = gmt_df[1].str.strip()
    gmt_df["to_delete"] = 0
    # delete all rows that do not have any molecule (3rd col)
    gmt_df.loc[gmt_df.iloc[:, 2].isna(), "to_delete"] = 1
    gmt_df = gmt_df.loc[gmt_df["to_delete"] != 1, :]
    gmt_df = gmt_df.drop(columns=["to_delete"])
    gmt_df.iloc[:, 1] = gmt_df.iloc[:, 1].str.replace(" / ", ", ")  # required
    if undesired_name_substring is not None:
        gmt_df.iloc[:, 1] = gmt_df.iloc[:, 1].str.replace(
            undesired_name_substring, "")
    return gmt_df


def purge_pathways(gmt_df, not_metabolic_reactions_pathways: List[str],
                   inconvenient_id_name_df: Union[None, pd.DataFrame]=None):
    # part 1:
    if inconvenient_id_name_df is not None:
        gmt_df["id"] = gmt_df.iloc[:, 0].tolist()  # column 0 has pathway_id
        gmt_df = gmt_df.loc[~gmt_df["id"].isin(
            inconvenient_id_name_df['pathway_id'].unique()), :].copy()

    gmt_df = gmt_df.drop(columns=["id"])
    gmt_df = gmt_df.dropna(axis=1, how="all")  # drop if entire column is NaN

    # part 2:
    gmt_df["bad"] = 0
    gmt_df["lower"] = gmt_df.iloc[:, 1].str.lower()  # column 1 has the names
    gmt_df.loc[gmt_df["lower"].str.contains(
        "|".join(not_metabolic_reactions_pathways), regex=True), "bad"] = 1
    gmt_df = gmt_df.loc[gmt_df["bad"] != 1, :]
    gmt_df = gmt_df.drop(columns=["bad", "lower"])

    return gmt_df


def gmt_to_piled_form(gmt_df, pathway_id_col: str, molecule_id_col: str):
    """return df with 2 columns:  pathway name, molecule kegg identifiers"""
    l_1 = []
    l_2 = []
    for i, r in gmt_df.iterrows():
        p_id = r.iloc[0]
        molecules = r.iloc[2:].dropna().tolist()
        l_1 += [p_id for i in range(len(molecules))]
        l_2 += molecules

    piled_df = pd.DataFrame({pathway_id_col: l_1, molecule_id_col: l_2})

    return piled_df


def gmts_merge_pile_splitid(dfs_dict, pathway_id_col: str,
                            pathway_name_col: str, molec_ids_pws_col: str):
    pathways_identities_dfs = {}
    for k in dfs_dict.keys():
        tmp = dfs_dict[k]
        pathways_identities_dfs[k] = pd.DataFrame(
            {pathway_id_col: tmp.iloc[:,0].tolist(),
             pathway_name_col: tmp.iloc[:,1].tolist()})
        dfs_dict[k] = gmt_to_piled_form(     # is a table in piled form
            tmp, pathway_id_col, molec_ids_pws_col)
    # end for
    pws_piled = combine_frames_identical_columns(dfs_dict)
    id_name_df = combine_frames_identical_columns(pathways_identities_dfs)

    return id_name_df, pws_piled


def pathways_and_annots_join(pws_piled: pd.DataFrame,
                             anno_mz_ids: pd.DataFrame, mol_id_in_pws: str) -> pd.DataFrame:
    """

    """
    assert len(list(pws_piled.columns)) == 2, \
        "too many columns in gmt piled table, bad format"
    assert len(list(anno_mz_ids.columns)) == 2, \
        "too many annotation columns: only mz and identifiers cols admitted"
    assert mol_id_in_pws in list(anno_mz_ids.columns), \
        f"no column {mol_id_in_pws}  in annotations table"

    anno_mz_ids = anno_mz_ids.drop_duplicates()
    accepted_ids_intersect = set(pws_piled[mol_id_in_pws].tolist()).intersection(
        set(anno_mz_ids[mol_id_in_pws].tolist()))
    pws_piled_ok = pd.merge(pws_piled, anno_mz_ids, on=mol_id_in_pws,
                            how="left")
    # only the intersection is kept:
    pws_piled_ok = pws_piled_ok.loc[pws_piled_ok[mol_id_in_pws].isin(list(
        accepted_ids_intersect)), :]

    return pws_piled_ok


def pws_piled_to_dict(pws_piled_mz, pathway_id: str, mz_or_identity_col: str):
    keptdf = pws_piled_mz[[pathway_id, mz_or_identity_col]]
    pathways_uniqs = keptdf[pathway_id].unique().tolist()

    out = {}
    for w in pathways_uniqs:
        items_w = keptdf.loc[keptdf[pathway_id] == w, 
                             mz_or_identity_col].astype(str).unique().tolist()
        out[w] = items_w

    return out


def generate_coordinates_arr(pixels_names_l: List[str]) -> np.array:
    """expects each name in the form 'xXX_yYY', e.g. 'x88_y99'"""
    coord_l = list()
    for k in pixels_names_l:   # eg k = 'x00_y00'
        items = k.split("_")
        coord_l.append(np.array([int(items[0].replace("x", "")),
                       int(items[1].replace("y", ""))]))

    return np.array(coord_l)


def is_log_transformed(sparse_matrix):

    if sparse_matrix.min().min() < 0:
        return True

    return False


def allocate_aliases_into_groups(aliases_df) -> (list, list):
    accepted_ids = set()
    precluded_ids = set()
    for i, row in aliases_df.iterrows():
        if row['LipidMaps'] != "empty":
            precluded_ids = precluded_ids.union({row['KEGG'], row['HMDB']})
            accepted_ids = accepted_ids.union(
                {row['LipidMaps']})  # preference for lipidmaps
        elif (row['LipidMaps'] == "empty") and (row['HMDB'] != "empty"):
            precluded_ids = precluded_ids.union({row['KEGG']})
            accepted_ids = accepted_ids.union(
                {row['HMDB']})  # second pref hmdb
        elif (row['LipidMaps'] == "empty") and (row['HMDB'] == "empty"):
            accepted_ids = accepted_ids.union(
                {row['KEGG']})  # kegg id : only if the other two are empty

    return list(accepted_ids), list(precluded_ids)


def aliases_groups(aliases_df, anno_df, usable_id):
    m_ids = anno_df[usable_id].unique()
    h = [u for u in m_ids if u.startswith("HMDB")]
    lm = [u for u in m_ids if u.startswith("LM")]
    kg = [u for u in m_ids if u.startswith("C")]
    aliases_df = aliases_df.loc[(aliases_df["HMDB"].isin(h)
                                 ) | (aliases_df["LipidMaps"].isin(lm)) | (
                                    aliases_df["KEGG"].isin(kg)), :]
    extra_kegg = set(kg) - set(aliases_df["KEGG"].unique())
    tail_df = pd.DataFrame({}, columns=aliases_df.columns)
    tail_df["KEGG"] = list(extra_kegg)
    aliases_df = pd.concat([aliases_df, tail_df])
    aliases_df = aliases_df.fillna("empty")
    aliases_df = aliases_df.replace(to_replace='', value="empty")

    accepted_ids, precluded = allocate_aliases_into_groups(aliases_df)

    dout = {'accepted' : accepted_ids, 'precluded': precluded}

    return dout
















