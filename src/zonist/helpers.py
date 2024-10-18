import sys
import logging
from typing import List, Dict
import pandas as pd
import numpy as np
import anndata as ad
import gc
import scanpy as sc
import harmonypy as hm
import squidpy as sq
from tqdm import tqdm


def reset_log_config(logger):
    """parameters for the log file"""
    out_log_file_name = "log.log"
    logging.basicConfig(encoding='utf-8', stream=sys.stdout)
    file_handler = logging.FileHandler(out_log_file_name)
    logger.addHandler(file_handler)
    logger.setLevel(logging.DEBUG)
    return logger
    

def combine_annotations(dfs_dict : Dict[str, pd.DataFrame]
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


def compute_metrics_for_variables(adata):
    adata.var['mz'] = adata.var.index.astype(float).to_numpy()
    #adata.var['rounded_mz'] = np.around(adata.var['mz'].to_numpy(),
    #                                    decimals=9)
    adata.var['variance'] = np.var(adata.X, axis=0)
    adata.var['mean'] = np.mean(adata.X, axis=0)
    adata.var['sum'] = np.sum(adata.X, axis=0)

    return adata


def get_rid_of_low_variance_mz(
        adata, cutoff_percentile_variance_tolerated: float):
    cutoff_mz_value = np.percentile(
        a=adata.var['variance'].to_numpy(),
        q=cutoff_percentile_variance_tolerated * 100)
    mz_keep = adata.var.loc[
        adata.var['variance'] >= cutoff_mz_value, "mz"].tolist()
    new_adata = adata[:, [str(o) for o in mz_keep]].copy()
    return new_adata


def select_n_most_variable_features(adata, n_top: int):
    assert 'variance' in adata.var.columns.tolist(), (
        "Error, 'variance' column not in adata.var")
    sorted_variance_df = adata.var.copy()
    sorted_variance_df = sorted_variance_df.sort_values(
        by='variance', ascending=True)  # .reset_index()
    top_features = sorted_variance_df.head(n_top)["mz"].tolist()  #iloc[0:n_top, "mz"].tolist()

    return adata[:, [str(o) for o in top_features]]


def find_common_mz__two_adata(adata1, adata2):
    mzs1 = adata1.var.index.astype(float).to_numpy()
    mzs2 = adata2.var.index.astype(float).to_numpy()
    common_mzs = set(mzs1).intersection(set(mzs2))
    if len(common_mzs) <= min(len(mzs1), len(mzs2)):
        print("Warning! the common mz are too few, "
              "consider an approximate matching")
    return common_mzs


def approximate_matching(mzs1: np.array, mzs2: np.array, TOLERANCE: float):
    print("running approximate matching")
    all_features = list()
    match_in_1 = list()
    match_in_2 = list()
    for i in range(len(mzs1)):
        all_features.append(mzs1[i])
        match_in_1.append(mzs1[i])
        u = np.argmin(np.abs(mzs2.copy() - mzs1[i]))
        if np.abs(mzs2[u] - mzs1[i]) <= TOLERANCE:
            match_in_2.append(mzs2[u])
        else:
            match_in_2.append(None)

    mz_of_2_not_found = [float(k) for k in mzs2 if float(k) not in match_in_2]
    for p2 in mz_of_2_not_found:
        all_features.append(p2)
        match_in_1.append(None)
        match_in_2.append(p2)

    return pd.DataFrame({"all_features": all_features,
                         "mzs1": match_in_1, "mzs2": match_in_2})


def purge_approximate_matching(df ):
    """purges cases of multiple approximate match, by keeping only the match
    with the minimum distance, e.g  mzs2 has repeated match 63.028146683:
    all_features     mzs1          mzs2
    63.029669847	63.029669847	63.028146683
    63.03229616	   63.03229616   	63.028146683
    63.027043535	63.027043535	63.028146683
    63.034922472	63.034922472	63.028146683

    output:
    63.029669847	63.029669847		NaN
    63.03229616	63.03229616		Nan
    63.027043535	63.027043535	63.028146683
    63.034922472	63.034922472		NaN
    """
    multiple_counts = df['mzs2'].value_counts(dropna=True)
    multiple_counts_mz = multiple_counts[multiple_counts >= 2].index.tolist()
    df['aux'] = 0  # not to change
    for n in multiple_counts_mz:
        arr_values_1 = df.loc[df['mzs2'] == n, 'mzs1'].to_numpy()
        the_item = np.argmin(np.abs(arr_values_1 - n))
        this_array = np.ones(len(arr_values_1))
        this_array[the_item] = 0
        df.loc[df['mzs2'] == n, 'aux'] = this_array

    df.loc[df['aux'] == 1, "mzs2"] = None
    df = df.drop(columns=['aux'])
    return df


def run_approximate_match_two_adata(adata1, adata2, TOLERANCE: float):
    mzs1 = adata1.var['mz'].astype(float).to_numpy()
    mzs2 = adata2.var['mz'].astype(float).to_numpy()
    df = approximate_matching(mzs1, mzs2, TOLERANCE)
    df = purge_approximate_matching(df)
    return df


def plot_gross_initial_clusters(adata, out_substring: str, leiden_resolution) -> None:
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=leiden_resolution)
    adata.uns['spatial'] = {'my_slice_name': {
        'scalefactors': {'fiducial_diameter_fullres': 100,
                         'spot_diameter_fullres': 1,
                         'tissue_hires_scalef': 0.5,
                         'tissue_lowres_scalef': 0.5}}}
    adata.obs = adata.obs.assign(leiden=adata.obs['leiden'].str.cat(
        np.repeat(out_substring, repeats=adata.obs.shape[0]), sep="|"
    ))
    sq.pl.spatial_scatter(adata, img=False, color='leiden', shape="square",
                          save=f"{out_substring}-leiden.pdf")


def plot_post_filtered_clusters(adata, out_substring):
    assert 'leiden' in adata.obs.columns.tolist(), "Error, leiden not in .obs"
    adata.uns['spatial'] = {'my_slice_name': {
        'scalefactors': {'fiducial_diameter_fullres': 100,
                         'spot_diameter_fullres': 1,
                         'tissue_hires_scalef': 0.5,
                         'tissue_lowres_scalef': 0.5}}}
    adata.obs = adata.obs.assign(leiden=adata.obs['leiden'].str.cat(
        np.repeat(out_substring, repeats=adata.obs.shape[0]), sep="|"
    ))
    sq.pl.spatial_scatter(adata, img=False, color='leiden', shape="square",
                          save=f"{out_substring}-leiden.pdf")



def harmonize_mz_identifiers_data2(adata2, df):
    # preference by the adata1 mz identifiers
    df2 = df.loc[~df['mzs2'].isna(), :]
    def order_good(df2, adata2):
        # order df2 by mzs2 identical as in adata2.var['mz']
        df2.index = df2["mzs2"]
        tmp = df2.T
        tmp = tmp[adata2.var["mz"].tolist()]
        return tmp.T

    df2 = df2.assign(harmo_mzs2=df2["mzs2"])
    for i, r in df2.iterrows():
        if not np.isnan(r["mzs1"]):
            df2.loc[i, "harmo_mzs2"] = r["mzs1"]
    df2 = order_good(df2, adata2)
    adata2.var['harmo_mzs2'] = df2["harmo_mzs2"].tolist()
    adata2.var.index = adata2.var['harmo_mzs2'].tolist()
    adata2.varm['transcripts'] = np.array(adata2.var['harmo_mzs2'].tolist())
    return adata2, df2






