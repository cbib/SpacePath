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

    return df_out


def spreadout_annotation_df(df: pd.DataFrame, col_mz: str, col_formula: str, 
    col_molecule_ids: str, col_molecule_names: str):
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

    comma_plus_space = ", " ## caution, do not remove the space

    tmp1 = df[[col_mz, col_formula, col_molecule_ids]]
    tmp2 = df[[col_molecule_names, 'database']]
    tmp1[col_molecule_ids] = tmp1[col_molecule_ids
                   ].str.split(comma_plus_space)
    tmp2[col_molecule_names] = tmp2[col_molecule_names
                   ].str.split(comma_plus_space)

    tmp1 = tmp1.explode(col_molecule_ids)
    tmp2 = tmp2.explode(col_molecule_names)
    df_out = pd.concat([tmp1, tmp2], axis=1)

    return df_out

def pull_offsample_info(df, col_mz: str, col_offsample: str):
    df = df[[col_mz, col_offsample]]
    df[col_offsample] = df[col_offsample].astype(bool)
    df = df.sort_values(by=[col_mz, col_offsample],
                        ascending=[True, True])
    df = df.drop_duplicates(subset=[col_mz, col_offsample],
                            keep='first')
    return df


def convert_metabo_annots_tokegg(
        annotdf: pd.DataFrame,  mappings_ref_ids: pd.DataFrame,
        col_molecule_ids: str, col_mz: str) -> pd.DataFrame:
    dfs_list = []
    common_notkegg_database = set(annotdf['database'].unique()).intersection(
        set(mappings_ref_ids['database'].unique()))
    for cd in list(common_notkegg_database):
        tmp_df = annotdf.loc[annotdf['database'] == cd, :]
        tmp_ref = mappings_ref_ids.loc[mappings_ref_ids['database'] == cd, :]
        joined = pd.merge(tmp_df, tmp_ref, left_on=col_molecule_ids,
                          right_on='not_kegg_id', how='left')
        joined = joined.assign(database=joined.database_y)
        joined = joined.drop(columns=['database_x', 'database_y',
                                      'not_kegg_id']).dropna()
        dfs_list.append(joined)
    # end for
    if "KEGG" in annotdf['database'].unique():
        kegg_tmp_df = annotdf.loc[annotdf['database'] == "KEGG", :]
        kegg_tmp_df = kegg_tmp_df.assign(KEGG=kegg_tmp_df[col_molecule_ids])
        dfs_list.append(kegg_tmp_df)

    final = pd.concat(dfs_list)
    final = final.sort_values(by=col_molecule_ids, ascending=True)
    final = final.drop_duplicates(subset=[col_mz, 'KEGG'], keep='first')
    final = final.sort_values(by=col_mz, ascending=True)

    return final


def identify_multi_mz_to_one_kegg(annotdf: pd.DataFrame,
                                  col_mz: str) -> bool:
    """True if multiple mz for one KEGG id exist in the annotation table"""
    out_bool = False
    df2 = annotdf[[col_mz, "KEGG"]]
    df2 = df2.drop_duplicates()
    df2 = df2.groupby("KEGG", as_index=False).size()
    df2 = df2.loc[df2['size'] >= 2, :].head()
    if df2.shape[0] > 0:
        out_bool = True

    return out_bool


def resolve_mz_multi_kegg(adata: ad.AnnData, df: pd.DataFrame) -> List[tuple]:
    """
    HMDB curators do efforts to map one-to-one, to KEGG ids.
    However, still many HMDB ids match to one single KEGG id, yielding
    multiple m/z <> repeated KEGG id in our piled table of annotations.
    To resolve such multi-m/z > 1 KEGG group, this function keeps the
    m/z (of such group) that has the highest variance of the quantifications.
    output: List[tuple] a list of tuples: (mz (float), kegg (str))
    """
    mean_arr, v_arr = sparsefuncs.mean_variance_axis(adata.X, axis=0)
    adata.var = adata.var.assign(computed_variance = v_arr)
    df2 = df.groupby("KEGG", as_index=False).size()
    unique_kegg_ids = df2.loc[df2['size'] == 1, "KEGG"].tolist()
    ok_mz_kegg_tuples = [tuple(t) for t in df.loc[
        df['KEGG'].isin(unique_kegg_ids), ["mz", "KEGG"]].to_numpy()]
    repeated_kegg_ids = df2.loc[df2['size'] >= 2, "KEGG"].unique()
    for a_kegg in repeated_kegg_ids:
        multi_annot_df = df.loc[df["KEGG"] == a_kegg, :]
        df_with_v_metric = adata.var.loc[
               adata.var['mz'].isin(multi_annot_df["mz"].tolist()),:]
        df_with_v_metric = df_with_v_metric.sort_values(
            by='computed_variance', ascending=False)
        mz_best = df_with_v_metric["mz"].tolist()[0] # first item mz
        ok_mz_kegg_tuples.append([mz_best, a_kegg])
        # end for
    return ok_mz_kegg_tuples


def compute_one_to_one_mz_keggid_df(annotdf: pd.DataFrame, adata_mz: ad.AnnData,
                                       col_mz: str) -> pd.DataFrame:
    ok_mz_kegg_tuples = resolve_mz_multi_kegg(adata_mz, annotdf)
    annotdf["keep"] = 0
    for tup in ok_mz_kegg_tuples:
        mz_i = tup[0]
        keggid = tup[1]
        annotdf.loc[(annotdf["mz"] == mz_i) & (
                    annotdf["KEGG"] == keggid), "keep"] = 1

    annotdf = annotdf.loc[annotdf["keep"] == 1]
    annotdf = annotdf.drop(columns=["keep"])

    return annotdf


def prettify_gmt(gmt_df, undesired_name_substring: Union[str, None]):
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


def purge_pathways(gmt_df, not_metabolic_reactions_pathways: List[str]):
    gmt_df["bad"] = 0
    gmt_df["lower"] = gmt_df.iloc[:, 1].str.lower()  # column 1 has the names
    gmt_df.loc[gmt_df["lower"].str.contains(
        "|".join(not_metabolic_reactions_pathways), regex=True), "bad"] = 1
    gmt_df = gmt_df.loc[gmt_df["bad"] != 1, :]
    gmt_df = gmt_df.drop(columns=["bad", "lower"])

    return gmt_df


def gmt_to_piled_form(gmt_df, col_pathway_id: str, col_molecule_id: str):
    """return df with 2 columns:  pathway name, molecule kegg identifiers"""
    l_1 = []
    l_2 = []
    for i, r in gmt_df.iterrows():
        p_id = r.iloc[0]
        molecules = r.iloc[2:].dropna().tolist()
        l_1 += [p_id for i in range(len(molecules))]
        l_2 += molecules

    piled_df = pd.DataFrame({col_pathway_id: l_1, col_molecule_id: l_2})

    return piled_df


def gmts_merge_pile_splitid(dfs_dict, col_pathway_id: str,
                            col_pathway_name: str, col_molecule_ids: str):
    pathways_identities_dfs = {}
    for k in dfs_dict.keys():
        tmp = dfs_dict[k]
        pathways_identities_dfs[k] = pd.DataFrame(
            {col_pathway_id: tmp.iloc[:,0].tolist(),
             col_pathway_name: tmp.iloc[:,1].tolist()})
        dfs_dict[k] = gmt_to_piled_form(     # is a table in piled form
            tmp, col_pathway_id, col_molecule_ids)
    # end for
    gmt_piled = combine_frames_identical_columns(dfs_dict)
    id_name_df = combine_frames_identical_columns(pathways_identities_dfs)

    return id_name_df, gmt_piled


def pathways_and_annots_join(gmt_piled: pd.DataFrame,
         anno_mz_kegg: pd.DataFrame, col_molecule_id: str) -> pd.DataFrame:
    assert len(list(gmt_piled.columns)) == 2, \
        "too many columns in gmt piled table, bad format"
    assert len(list(anno_mz_kegg.columns)) == 2, \
        "too many annotation columns: only the mz and KEGG are admitted"
    assert col_molecule_id in list(anno_mz_kegg.columns), \
        f"no column {col_molecule_id}  in annotations table"

    anno_mz_kegg = anno_mz_kegg.drop_duplicates()
    accepted_kegg_intersect = set(gmt_piled["KEGG"].tolist()).intersection(
        set(anno_mz_kegg["KEGG"].tolist()))
    gmt_piled_ok = pd.merge(gmt_piled, anno_mz_kegg, on="KEGG", how="left")
    # only the intersection is kept:
    gmt_piled_ok = gmt_piled_ok.loc[gmt_piled_ok["KEGG"].isin(list(
        accepted_kegg_intersect)), :]

    return gmt_piled_ok


def pws_piled_to_dict(gmt_piled_mz, pathway_id: str, mz_or_kegg: str):
    keptdf = gmt_piled_mz[[pathway_id, mz_or_kegg]]
    pathways_uniqs = keptdf[pathway_id].unique().tolist()

    out = {}
    for w in pathways_uniqs:
        items_w = keptdf.loc[keptdf[pathway_id] == w, 
                             mz_or_kegg].astype(str).unique().tolist()
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











