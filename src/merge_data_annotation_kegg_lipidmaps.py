# johaGL 2024
# usage in Python (not bash):
# from XXX import combine_kegg_lipidmaps_metaspace_annotations
# ... # you have loaded all 3 dataframes and defined the FLOATPOSITIONS int
#  res = combine_kegg_lipidmaps_metaspace_annotations(
#         pre_k2l, kegg_df, lipim_df, FLOATPOSITIONS)
# returns object (pd.DataFrame)

import os
import pandas as pd
import numpy as np


def treat_conversion_table_kegg_lipidmaps(pre_k2l):
    """
    Several lipidmaps id may match with one kegg id,
    make rows to be unique combinations kegg <-> lipidmap id
    k2l (kegg to lipidmaps) dataframe
    """
    print('Stacking conversion table.')
    pre_k2l = pre_k2l.dropna()  # delete rows that are NA for KEGG id
    new_kegg_list = list()
    new_lipidmaps_list = list()
    for i, r in pre_k2l.iterrows():
        if len(r['LIPIDMAPS'].split(" ")) >= 2:
            lipidmaps_ids_list = r['LIPIDMAPS'].split(" ")
            for elem in lipidmaps_ids_list:
                new_kegg_list.append(r['KEGG'].strip())
                new_lipidmaps_list.append(elem.strip())
        else:
            new_kegg_list.append(r['KEGG'].strip())
            new_lipidmaps_list.append(r['LIPIDMAPS'].strip())

    k2l = pd.DataFrame(
        {'KEGG_from_conversion': new_kegg_list,
         'LIPIDMAPS_from_conversion': new_lipidmaps_list})
    return k2l


def process_metaspace_annotations_tables(df, table_type: str):
    """ annotations downloaded from metaspace follow the logics :
    mz      fdr     moleculeNames        moleculeIds
    111.33  0.05    mol1, mol2, mol3    id1, id2, id3
    Transform the dataframe to have unique combinations of mz, name and id
    """
    print("Stacking input table of type:", table_type,
          "-> Input 1st record: ", df["moleculeIds"].tolist()[0])
    req_columns = ['mz', 'fdr', 'moleculeNames', 'moleculeIds']
    if len(set(df.columns.tolist()).intersection(
            set(req_columns))) != 4:
        raise ValueError(f"Error, required columns: {req_columns} not found")
    df = df[req_columns]
    names_list = list()
    ids_list = list()
    mz_list = list()
    fdr_list = list()
    for i, r in df.iterrows():
        names = r['moleculeNames'].split(", ")
        ids = r['moleculeIds'].split(", ")
        for n, k in zip(names, ids):
            names_list.append(n.strip())
            ids_list.append(k.strip())
            mz_list.append(r['mz'])
            fdr_list.append(r['fdr'])
    df_out = pd.DataFrame({
        'mz': mz_list, 'fdr': fdr_list,
        'molecule_name': names_list, 'moleculeIds': ids_list})
    return df_out


def combine_tables_as_kegg(k2l, kegg_df, lipim_df, FLOATPOSITIONS):
    print("Combining annotations information.")
    # rename columns
    kegg_df.columns = ['mz', 'fdr', 'molecule_name', 'KEGG_ID']
    lipim_df.columns = ['mz', 'fdr', 'molecule_name', 'LIPIDMAPS_ID']
    # transform lipidmaps into kegg id:
    foo = pd.merge(lipim_df, k2l, left_on='LIPIDMAPS_ID',
                   right_on='LIPIDMAPS_from_conversion', how='left')
    # only retain those having KEGG
    foo = foo.loc[~foo['KEGG_from_conversion'].isna(), :]
    kegg_df['mz'] = np.around(kegg_df['mz'].to_numpy(), FLOATPOSITIONS)
    foo = foo.assign(mz=np.around(foo['mz'].to_numpy(), FLOATPOSITIONS))
    foo = foo[['mz', 'fdr', 'molecule_name', 'KEGG_from_conversion']]
    foo = foo.assign(db_origin='lipidmaps')
    kegg_df = kegg_df.assign(db_origin='kegg')
    foo.columns = ['mz', 'fdr', 'molecule_name', 'KEGG_ID', 'db_origin']
    res = pd.concat([kegg_df, foo], ignore_index=True)  # vertically

    return res


def clean_combined_annotation(res):
    # sort by smallest to biggest fdr; if duplic kegg id, keep smallest fdr
    res = res.sort_values('fdr', ascending=True)
    # do not allow repeated kegg id, retain the one with min fdr (first)
    res = res.drop_duplicates(subset=['KEGG_ID'], keep="first")
    # debug_var = res.loc[res['KEGG_ID'] == 'C00059',:].copy()
    res = res.sort_values('mz', ascending=True)

    return res


def combine_kegg_lipidmaps_metaspace_annotations(
        pre_k2l: pd.DataFrame, kegg_df: pd.DataFrame, lipim_df: pd.DataFrame,
        FLOATPOSITIONS: int) -> pd.DataFrame:
    """
    call steps to run the combination of kegg and lipidmaps annotations
    downloaded from metaspace web results page
    """
    print("Running combination of KEGG and Lipidmaps annotations.")
    k2l = treat_conversion_table_kegg_lipidmaps(pre_k2l)
    kegg_df = process_metaspace_annotations_tables(kegg_df, "KEGG")
    lipim_df = process_metaspace_annotations_tables(
        lipim_df, "Lipidmaps")

    res = combine_tables_as_kegg(k2l, kegg_df,
                                   lipim_df, FLOATPOSITIONS)

    res = clean_combined_annotation(res)
    print("Done")
    return res


if __name__ == "__main__" :
    FLOATPOSITIONS = 8
    os.chdir("~/spatial_thesis/spatial_metabolomics")
    # this part would go into a jupyter when using notebooks
    kegg_f = "perfbrain1-5_04m_nedc_14x_50um_recal--metaspace_kegg.csv"
    lipim_f = "perfbrain1-5_04m_nedc_14x_50um_recal--metaspace_lipidmaps.csv"
    kegg_to_lipidmaps_file = "table2_kegg_lipidmaps_1716991921.9568334.txt"

    # load the correspondances of kegg <-> lipidmaps, file 2 directories above:
    pre_k2l = pd.read_csv(os.path.join("data/db_tables",
                          kegg_to_lipidmaps_file), sep='\t')

    kegg_df = pd.read_csv(os.path.join("data/annotations_tables", kegg_f))
    lipim_df = pd.read_csv(os.path.join("data/annotations_tables", lipim_f))

    res = combine_kegg_lipidmaps_metaspace_annotations(
        pre_k2l, kegg_df, lipim_df, FLOATPOSITIONS)








