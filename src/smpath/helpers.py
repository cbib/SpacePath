import logging
import os
import sys
from datetime import datetime
from typing import Dict, List, Union, Tuple

import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem

import scipy.sparse


LOGFILE = datetime.now().strftime("%y-%m-%-d_%H-%M-%S-%f.log")


def setup_logger(logger, prefix=''):
    """parameters for the log file"""
    if not os.path.exists(os.path.join(os.getcwd(), "logs")):
        os.makedirs(os.path.join(os.getcwd(), "logs"))

    logging.basicConfig(encoding='utf-8', stream=sys.stdout)
    file_handler = logging.FileHandler(
        os.path.join("logs", f'{prefix}{LOGFILE}'))
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    return logger


def secured_trios(identifiers_validated: List[str],
                  smiles_validated: List[str]) -> Tuple[
    List[str], List[str], List[DataStructs.cDataStructs.ExplicitBitVect]
]:
    """
    returns the lists of safely validated elements:
        identifiers : e.g. ['HMDB005474523, ...]
        smiles: e.g. ['COOHCOOHCOOH', ...]
        fingerprints: list of RDKit.DataStructs objects (each object is a fingerprint)
    the lists represent matched information (zippable)
    """
    fpgen = AllChem.GetRDKitFPGenerator()
    smiles_final = list()
    identifiers_final = list()
    fingerprints_final = list()
    for id, smi in zip(identifiers_validated, smiles_validated):
        try:
            chem_molecule = Chem.MolFromSmiles(smi)
            fingerprint = fpgen.GetFingerprint(chem_molecule)
            smiles_final.append(smi)
            identifiers_final.append(id)
            fingerprints_final.append(fingerprint)
        except Exception as e:
            print(e, "inappropriate smile for %s %s, skip", id, smi)

    return identifiers_final, smiles_final, fingerprints_final


def combine_frames_identical_columns(dfs_dict: Dict[str, pd.DataFrame]
                                     ) -> pd.DataFrame:
    list_of_frames = []
    for k in dfs_dict.keys():
        list_of_frames.append(dfs_dict[k])

    df_out = pd.concat(list_of_frames)
    df_out = df_out.drop_duplicates()

    return df_out


def spreadout_annotation_df(df: pd.DataFrame, mz_col: str, formula_col: str,
                            molecule_ids_col: str,
                            molecule_names_col: str) -> pd.DataFrame:
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

    delimiter_colon_plus_space = ", "  # please do not remove the space

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


def purge_annotations(df_notranslated: pd.DataFrame,
                      discard_pharmaceuticals: bool,
                      discard_plant_compounds: bool,
                      discard_other_non_animal: bool,
                      library_path: str) -> pd.DataFrame:
    ids_col = "moleculeIds"
    f_picker = {"pharma": {'istrue': discard_pharmaceuticals,
                           'file': "pharmaceuticals.tsv"},
                "plant": {'istrue': discard_plant_compounds,
                          'file': "plant_compounds.tsv"},
                "other": {'istrue': discard_other_non_animal,
                          'file': "non_animal_other.tsv"}
                }
    for k in f_picker:
        if f_picker[k]['istrue']:
            tmp_df = pd.read_csv(
                os.path.join(library_path, "resources", f_picker[k]['file']),
                sep="\t", header=0, index_col=None)
            df_notranslated = df_notranslated.loc[
                ~df_notranslated[ids_col].isin(tmp_df[ids_col]),
                :].copy()

    return df_notranslated


def pick_mz_quality_infos(df, mz_col: str, offsample_col: str,
                          total_intens_col: str, fdr_column: str,
                          msm_column: str, formula_col: str) -> pd.DataFrame:
    """ """
    df = df[[mz_col, offsample_col, total_intens_col, fdr_column,
             msm_column, formula_col]]
    df = df.sort_values(by=[mz_col, offsample_col, fdr_column, msm_column],
                        ascending=[True, True, True, False])
    df = df.drop_duplicates(subset=[mz_col], keep='first')

    return df


def convert_metabo_annots_tokegg(
        annotdf: pd.DataFrame, mappings_ref_ids: pd.DataFrame,
        molecule_ids_col: str, mz_col: str,
        usable_id: str = "KEGG") -> pd.DataFrame:
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
        ids_l = ann_df.loc[
            ann_df[mz_col].astype(float) == mzfloat, id_col].unique()
        nn.append(len(ids_l))

    var_df['compound_ids_n__per_mz'] = nn

    return var_df


def get_unique_pairs_list(ann_df, id_col, mz_col):
    """counts the occurrences of compound identifiers in the annotation,
    and finds the pairs of (mz, compoundID) that are unique"""
    count_df = ann_df.groupby(id_col, as_index=False).size()  # count ocurr
    unique_compound_ids = count_df.loc[count_df['size'] == 1, id_col].tolist()
    approved_mz_id_tuples = [tuple(t) for t in ann_df.loc[
        ann_df[id_col].isin(unique_compound_ids), [mz_col,
                                                   id_col]].to_numpy()]

    return count_df, approved_mz_id_tuples


def get_approved_mz_id_pairs(ann_df: pd.DataFrame, var_df: pd.DataFrame,
                             mz_col: str, id_col: str, msm_col: str,
                             fdr_col: str, offsam: str) -> List[tuple]:
    """
    Provides a final list of approved (m/z, compoundID) tuples.
    Note:
        In our piled table of annotations, there are cases of reverse
        multimapping (multiple m/z, one compound) in pairs or tuples.
           E.g: (188.4545,  Fucose), (559.869, Fucose), (345.8556, Fucose).
        in such case, retains the two top better pairs by quality criteria.
    output: List[tuple] a list of tuples: (mz (float), id (str))
    """
    ann_df = ann_df.loc[ann_df[mz_col].astype(float).isin(
        var_df[mz_col].astype(float)), :]
    count_df, approved_mz_id_tuples = get_unique_pairs_list(ann_df, id_col,
                                                            mz_col)

    var_df = compute_compounds_n_per_mz(var_df, ann_df, mz_col, id_col)
    compounds_multiassigned = count_df.loc[
        count_df['size'] >= 2, id_col].unique()
    var_df = var_df.loc[var_df['compound_ids_n__per_mz'] >= 1, :]
    for one_id in compounds_multiassigned:
        multi_annot_df = ann_df.loc[ann_df[id_col] == one_id, :]
        curr_qual = var_df.loc[var_df[mz_col].astype(float).isin(
            multi_annot_df[mz_col].astype(float)), :].copy()
        if 'I' in curr_qual.columns:
            curr_qual = curr_qual.sort_values(
                by=['compound_ids_n__per_mz', msm_col, fdr_col, offsam, 'I'],
                # always respect this order !
                ascending=[True, False, True, True, False])
        else:
            curr_qual = curr_qual.sort_values(
                by=['compound_ids_n__per_mz', msm_col, fdr_col, offsam],
                # in this order !
                ascending=[True, False, True, True])
        the_top_mzs = curr_qual.head(2)[mz_col].tolist()
        for a_top_mz in the_top_mzs:
            approved_mz_id_tuples.append((a_top_mz, one_id))

    return approved_mz_id_tuples


def best_mz_id_pairs_df(annotdf: pd.DataFrame,
                        var_df: pd.DataFrame,
                        mz_col: str, id_col: str, msm_col: str, fdr_col: str,
                        offsam: str
                        ) -> pd.DataFrame:
    ok_mz_kegg_tuples = get_approved_mz_id_pairs(
        annotdf, var_df, mz_col, id_col, msm_col, fdr_col, offsam)

    annotdf["keep"] = 0
    for tup in ok_mz_kegg_tuples:
        mz_i = tup[0]
        metabolite_id = tup[1]
        annotdf.loc[(annotdf[mz_col] == mz_i) & (
                annotdf[id_col] == metabolite_id), "keep"] = 1

    # rejected = annotdf.loc[annotdf["keep"] == 0, :].copy()  # for debugging
    annotdf = annotdf.loc[annotdf["keep"] == 1, :].copy()
    annotdf = annotdf.drop(columns=["keep"])
    annotdf = annotdf.sort_values(by='database', ascending=True)
    annotdf = annotdf.drop_duplicates(subset=[mz_col, id_col],
                                      keep='first')  # favors CoreMetabolome, ok

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
                   inconvenient_id_name_df: Union[None, pd.DataFrame] = None):
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


def gmt_to_piled_form(gmt_df, set_id_col: str, molecule_id_col: str):
    """return df with 2 columns:  pathway name, molecule kegg identifiers"""
    l_1 = []
    l_2 = []
    for i, r in gmt_df.iterrows():
        p_id = r.iloc[0]
        molecules = r.iloc[2:].dropna().tolist()
        l_1 += [p_id for i in range(len(molecules))]
        l_2 += molecules

    piled_df = pd.DataFrame({set_id_col: l_1, molecule_id_col: l_2})

    return piled_df


def gmts_merge_pile_splitid(dfs_dict: Dict[str, pd.DataFrame],
                            pathway_id_col: str,
                            set_name_col: str, molec_ids_pws_col: str
                            ) -> Tuple[pd.DataFrame, pd.DataFrame]:
    pathways_identities_dfs = {}
    for k in dfs_dict.keys():
        tmp = dfs_dict[k]
        pathways_identities_dfs[k] = pd.DataFrame(
            {pathway_id_col: tmp.iloc[:, 0].tolist(),
             set_name_col: tmp.iloc[:, 1].tolist()})
        dfs_dict[k] = gmt_to_piled_form(  # is a table in piled form
            tmp, pathway_id_col, molec_ids_pws_col)
    # end for
    pws_piled = combine_frames_identical_columns(dfs_dict)
    id_name_df = combine_frames_identical_columns(pathways_identities_dfs)

    return id_name_df, pws_piled


def pathways_and_annots_join(pws_piled: pd.DataFrame,
                             anno_mz_ids: pd.DataFrame,
                             mol_id_in_pws: str) -> pd.DataFrame:
    """
    Filters the pathways piled table, by retaining only the elements
    in the intersection with the annotations dataframe.
    This also adds the 'mz' column to the pathways piled table.
    :pws_piled: pathways df with columns: moleculeIds (or KEGG), pathway_id
    :anno_mz_ids: annotations df with columns: mz, moleculeIds (or KEGG)
    :mol_id_in_pws: "KEGG" or "moleculeIds"
    """
    assert len(list(pws_piled.columns)) == 2, \
        "too many columns in gmt piled table, bad format"
    assert len(list(anno_mz_ids.columns)) == 2, \
        "too many annotation columns: only mz and identifiers cols admitted"
    assert mol_id_in_pws in list(anno_mz_ids.columns), \
        f"no column {mol_id_in_pws}  in annotations table"

    anno_mz_ids = anno_mz_ids.drop_duplicates()
    accepted_ids_intersect = set(
        pws_piled[mol_id_in_pws].tolist()).intersection(
        set(anno_mz_ids[mol_id_in_pws].tolist()))
    pws_piled_ok = pd.merge(pws_piled, anno_mz_ids, on=mol_id_in_pws,
                            how="left")
    # only the intersection is kept:
    pws_piled_ok = pws_piled_ok.loc[pws_piled_ok[mol_id_in_pws].isin(list(
        accepted_ids_intersect)), :]

    return pws_piled_ok


def pws_piled_to_dict(pws_piled_mz, pathway_id: str, mz_or_identity_col: str
                      ) -> Dict[str, List[str]]:
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
    for k in pixels_names_l:  # eg k = 'x00_y00'
        items = k.split("_")
        coord_l.append(np.array([int(items[0].replace("x", "")),
                                 int(items[1].replace("y", ""))]))

    return np.array(coord_l)


def is_log_transformed(
        sparse_matrix: Union[scipy.sparse.csr_matrix, np.array]) -> bool:
    if sparse_matrix.min().min() < 0:
        return True

    return False


def yield_clean_noalias_ids(annots_with_aliases_df) -> List[str]:
    accepted_ids = set()
    for i, row in annots_with_aliases_df.iterrows():
        if row['LipidMaps'] != "empty":
            accepted_ids = accepted_ids.union(
                {row['LipidMaps']})  # preference for lipidmaps
        elif (row['LipidMaps'] == "empty") and (row['HMDB'] != "empty"):
            accepted_ids = accepted_ids.union(
                {row['HMDB']})  # second pref hmdb
        elif (row['LipidMaps'] == "empty") and (row['HMDB'] == "empty"):
            accepted_ids = accepted_ids.union(
                {row['KEGG']})  # kegg id : only if the other two are empty

    return list(accepted_ids)


def ids_not_in_internal_aliases_table(aliases_df, hmdb_ids, lm_idf, kegg_ids
                                      ) -> pd.DataFrame:
    exclusive_h = list(
        set(hmdb_ids) - set(aliases_df["HMDB"].unique().tolist()))
    exclusive_lm = list(
        set(lm_idf) - set(aliases_df["LipidMaps"].unique().tolist()))
    exclusive_kg = list(
        set(kegg_ids) - set(aliases_df["KEGG"].unique().tolist()))

    annots_nomate_ids = pd.DataFrame({
        "HMDB": exclusive_h + ['' for i in range(len(exclusive_lm))] + [
            '' for i in range(len(exclusive_kg))],
        "LipidMaps": ['' for i in range(len(exclusive_h))
                      ] + exclusive_lm + ['' for i in
                                          range(len(exclusive_kg))],
        "KEGG": ['' for i in range(len(exclusive_h))
                 ] + ['' for i in range(len(exclusive_lm))] + exclusive_kg
    })

    return annots_nomate_ids


def aliases_handler(aliases_df, anno_df, usable_id: str) -> List[str]:
    m_ids = anno_df[usable_id].unique()
    hmdb_ids = [u for u in m_ids if u.startswith("HMDB")]
    lm_ids = [u for u in m_ids if u.startswith("LM")]
    kegg_ids = [u for u in m_ids if u.startswith("C")]

    annots_nomate_ids = ids_not_in_internal_aliases_table(aliases_df,
                                                          hmdb_ids, lm_ids,
                                                          kegg_ids)
    # intersection of aliases_df and current annots
    annots_aliases = aliases_df.loc[(aliases_df["HMDB"].isin(hmdb_ids)
                                     ) | (aliases_df["LipidMaps"].isin(
        lm_ids)) | (aliases_df["KEGG"].isin(kegg_ids)), :].copy()
    del aliases_df
    annots_aliases = pd.concat([annots_aliases, annots_nomate_ids])
    annots_aliases = annots_aliases.fillna("empty")
    annots_aliases = annots_aliases.replace(to_replace='', value="empty")
    accepted_ids = yield_clean_noalias_ids(annots_aliases)

    return accepted_ids
