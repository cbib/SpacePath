import logging
from typing import List, Tuple, Union

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix
import smpath.helpers as helpers
from smpath.data import MetabolomeAnnotations
import squidpy as sq

logger = logging.getLogger(__name__)
logger = helpers.setup_logger(logger)

""""
Functions related to mz intensities, formulae quality metrics, and to weights of the assignments
"""


# weights of the assignments

def build_identifiers_weights_df(adata_mz: ad.AnnData,
                                 X_metabolites_probs_df: pd.DataFrame,
                                 clean_annotations: pd.DataFrame,
                                 usable_id: str) -> pd.DataFrame:
    mz_ids = clean_annotations[['mz', usable_id]].copy()
    mz_ids = mz_ids.drop_duplicates()
    mz_ids = mz_ids.sort_values(by=usable_id)
    unique_compounds_ids = mz_ids[usable_id].unique().tolist()
    tmp_X = dict()
    for u_id in unique_compounds_ids:
        mz_one_or_more = mz_ids.loc[mz_ids[usable_id] == u_id, 'mz'].astype(
            str).tolist()
        sub_adata = adata_mz[:, mz_one_or_more].copy()
        mz_arr = np.array(sub_adata.X.todense().T).flatten()
        if len(mz_one_or_more) > 1:
            tmp = mz_arr.sum(axis=0)
            mz_arr = tmp

        curr_product_arr = X_metabolites_probs_df[u_id].to_numpy() * mz_arr

        tmp_X[u_id] = np.array(curr_product_arr, dtype=np.float32)

    metabolites_computed_w = pd.DataFrame(data=tmp_X,
                                          columns=list(tmp_X.keys()),
                                          # compounds ids
                                          index=adata_mz.obs.index.tolist())  # spots ids
    return metabolites_computed_w


def compute_mz_binary_df(adata_mz: ad.AnnData) -> pd.DataFrame:
    X_presence_absence_mask = adata_mz.X.todense() > 0
    X_mz_binary_df = pd.DataFrame(
        np.array(X_presence_absence_mask, dtype=int),
        columns=adata_mz.var['mz'].to_numpy(),
        index=adata_mz.obs.index.tolist())

    return X_mz_binary_df


def compute_metabolites_probs_df(adata_mz: ad.AnnData,
                                 annots_df: pd.DataFrame,
                                 usable_id, mz_col="mz") -> pd.DataFrame:
    """
    computes the matrix-like dataframe of assignments' probabilities
    the probability of one assignment (compound id) is computed as 1/n,
     where n is the number of clean assignments (clean compounds identifiers).

    Note:
    for a current id in the nested loop, already registered from another mz,
    as it can happen (checked), retains the highest probability
    """
    X_mz_binary_df = compute_mz_binary_df(adata_mz)
    preX_probs = dict()
    for a_mz in X_mz_binary_df.columns:
        i_arr = X_mz_binary_df[a_mz].to_numpy()
        the_ids_here = annots_df.loc[
            annots_df[mz_col] == a_mz, usable_id].unique().tolist()
        n = len(the_ids_here)
        for kid in the_ids_here:
            registered_prob = 0.0
            if kid in preX_probs.keys():
                registered_prob = preX_probs[kid][0]
            # attribute the best probability between registered (if any) and current
            compound_probability = max(registered_prob, 1 / n)
            res = i_arr * np.repeat(compound_probability, len(i_arr))
            preX_probs[kid] = np.array(res, dtype=np.float32)

    X_metabolites_probs_df = pd.DataFrame(preX_probs)
    X_metabolites_probs_df.index = X_mz_binary_df.index

    return X_metabolites_probs_df


# normalization of m/z's (ions) and autocorrelation

def tic_normalization(
        intensity_matrix: Union[csr_matrix, np.array]) -> np.matrix:
    """
    Perform Total Ion Count (TIC) normalization on a raw 2D intensity matrix
    of m/z features
    :param intensity_matrix: rows are spots, columns are mz features
    :return:  TIC normalized intensity matrix
    :return type: numpy.matrix (verified)
    """
    logger.info(
        "Applying Total Ion Count (TIC) normalization: processes the m/z "
        "features, across spots")

    total_ion_count = intensity_matrix.sum(axis=1)
    scale_factor = np.median(np.array(total_ion_count))
    # Avoid division by zero
    total_ion_count = np.array(
        total_ion_count) + 1  # avoid any zero error div
    # Normalize the intensity matrix
    normalized_matrix = scale_factor * (intensity_matrix / total_ion_count)

    return normalized_matrix


def compute_moran_i(adata: ad.AnnData, n_neighs: int = 20,
                    n_perms: int = 5) -> ad.AnnData:
    """
    computes Moran's I with the squidpy implementation (spatial_autocorr).

     If not `spatial_connectivities` in .obsp attr of the AnnData object,
     the squidpy's 'spatial_neighbors' is also applied to satisfy spatial_autocorr
    """
    if 'spatial_connectivities' not in list(adata.obsp.keys()):
        logger.info('computing spatial neighbors, with n_neighs: %s',
                    n_neighs)
        sq.gr.spatial_neighbors(adata, n_neighs=n_neighs)
    logger.info("Computing Moran's I autocorrelation "
                "on data with shape (obs,vars) %s", adata.shape)
    sq.gr.spatial_autocorr(adata, mode="moran", n_perms=n_perms, n_jobs=1,
                           seed=0)

    return adata


# phases of classification for m/z

def part_a__mz_quality_mark(annotations_obj, msm_cutoff, offSample_cutoff,
                            fdr_cutoff
                            ) -> Tuple[MetabolomeAnnotations, List[str]]:
    annotations_obj.pull_mz_quality_df()
    annotations_obj.mz_quality['quality_msm_off_fdr'] = "bad"
    annotations_obj.mz_quality.loc[
        (annotations_obj.mz_quality[annotations_obj.msm_column] >= msm_cutoff
         ) & (annotations_obj.mz_quality[
                  annotations_obj.offsample_column] <= offSample_cutoff  #
              ) & (annotations_obj.mz_quality[
                       annotations_obj.fdr_column] <= fdr_cutoff),
        'quality_msm_off_fdr'] = "good"

    mz_to_keep_part1 = annotations_obj.mz_quality.loc[
        annotations_obj.mz_quality[
            'quality_msm_off_fdr'] == "good", "mz"].astype(
        str).tolist()

    return annotations_obj, mz_to_keep_part1


def custom_join_left(df1, df2, keycol: str = 'mz') -> Union[
    pd.DataFrame, None]:
    df1 = df1.copy()
    df2 = df2.copy()
    df2[keycol] = df2[keycol].astype(str)
    df2 = df2.drop_duplicates(subset=keycol)
    df2 = df2.set_index(keycol)
    try:
        df2 = df2.loc[
            df1[keycol].astype(str).tolist(), :]  # same row order as 1
        # np.all(df2.index.astype(str).to_numpy(
        # ) == df1[keycol].astype(str).to_numpy())
        for col in df2.columns:
            df1[col] = df2[col]
        return df1
    except KeyError:
        logger.error("fatal (no match), custom_join_left impossible!")
        return None
    except Exception as e:
        logger.error("fatal (no match), custom_join_left impossible! %s", e)
        return None


def reset_anndata_var_columns(adata_mz):
    """resets adata.var columns if user re-runs analysis """
    try:
        columns_to_reset = set(adata_mz.var.columns.tolist()) - set(["mz"])
        adata_mz.var = adata_mz.var.drop(columns=list(columns_to_reset))
    except KeyError:
        pass

    return adata_mz


def part_b__mz_quality_by_autocorr(adata_mz: ad.AnnData,
                                   mz_good_from_part_a: List[str],
                                   annotations_obj: MetabolomeAnnotations,
                                   moran_I_cutoff: float,
                                   moran_padj_cutoff: float,
                                   var_sim_cutoff: float,
                                   n_neighs: int, n_perms: int
                                   ) -> Tuple[
    pd.DataFrame, List[str], pd.DataFrame]:
    # prepare quality merged df
    adata_mz = reset_anndata_var_columns(adata_mz)
    all_quality_df = adata_mz.var.copy()
    all_quality_df["mz"] = all_quality_df["mz"].astype(str).tolist()
    annotations_obj.mz_quality['mz'] = annotations_obj.mz_quality[
        'mz'].astype(
        str).tolist()
    all_quality_df = custom_join_left(
        all_quality_df, annotations_obj.mz_quality[
            ['mz', annotations_obj.msm_column,
             annotations_obj.offsample_column,
             annotations_obj.fdr_column, annotations_obj.formula_col,
             "quality_msm_off_fdr"]],
        keycol='mz')

    copy_mz = adata_mz[
        :, mz_good_from_part_a].copy()  # only compute moran for the ok_part1 features
    try:
        if copy_mz.uns['moranI'].shape[0] != len(mz_good_from_part_a):
            del copy_mz.uns['moranI']
    except KeyError:
        pass
    except Exception as e:
        logger.error(e)
        pass

    if "moranI" not in list(copy_mz.uns.keys()):
        if not helpers.is_log_transformed(copy_mz.X):
            sc.pp.log1p(
                copy_mz)  # log-transform matrix, as data is heavily skewed.
            # this log transformed 'copy_mz' is not returned by the function
        copy_mz = compute_moran_i(copy_mz, n_neighs=n_neighs, n_perms=n_perms)
    else:

        logger.info(
            f"using the available .uns['moranI'] with dimension: {copy_mz.uns['moranI'].shape}")
        assert copy_mz.uns['moranI'].shape[0] == len(mz_good_from_part_a), (
            "Fatal: moranI table dimensions do not match with length of ions list passing the 'Part a'."
            "Killed! (solution: delete .uns['moranI'] dataframe and try again)")
    moran_mz_df = copy_mz.uns['moranI'].copy()
    # discard features with extremely low simulated variance
    moran_mz_df = moran_mz_df.loc[
        moran_mz_df['var_sim'] >= var_sim_cutoff].copy()

    moran_mz_df['mz'] = moran_mz_df.index.astype(str).tolist()
    all_quality_df = pd.merge(all_quality_df,
                              moran_mz_df[
                                  ['mz', 'I', 'pval_sim_fdr_bh', 'var_sim']],
                              how="left", left_on="mz", right_on="mz")
    mz_to_keep_part2 = all_quality_df.loc[
        (all_quality_df['I'] >= moran_I_cutoff) & (
                all_quality_df['pval_sim_fdr_bh'] <= moran_padj_cutoff),
        "mz"].astype(str).tolist()

    all_quality_df['quality_autocorr'] = 'bad'
    all_quality_df.loc[
        all_quality_df['mz'].astype(str).isin(mz_to_keep_part2),
        'quality_autocorr'] = "good"

    return all_quality_df, mz_to_keep_part2, moran_mz_df


def adata_mz_insert_full_quality(adata_mz: ad.AnnData, moran_mz_df,
                                 all_quality_df) -> ad.AnnData:
    adata_mz.uns['moranI'] = moran_mz_df

    all_quality_df['final_quality'] = "bad"
    all_quality_df.loc[(all_quality_df['quality_msm_off_fdr'] == "good") & (
            all_quality_df['quality_autocorr'] == "good"), 'final_quality'
    ] = "good"

    if np.all(np.array(adata_mz.var.index.tolist()) == np.array(
            all_quality_df['mz'].tolist())):
        all_quality_df.index = all_quality_df['mz'].astype(
            str).tolist()
        all_quality_df['mz'] = all_quality_df['mz'].astype(
            float).tolist()
        adata_mz.var = all_quality_df

    else:
        raise AssertionError(
            "Fatal: could not assembly final object. Aborted.")

    return adata_mz


# matchers

def verification(adata_mz: ad.AnnData,
                 annotations_obj: MetabolomeAnnotations) -> None:
    try:
        bool_equal = np.all(np.array(
            sorted(annotations_obj.no_translated["mz"].unique())) == np.array(
            sorted(adata_mz.var["mz"].tolist())))
        logger.info("perfect match between mz in me_a.df "
                    f"and adata.X: => %s <=", bool_equal)
    except ValueError:
        logger.error("mz do not match between adata_mz and annotation table")
        logger.error("unique annotated mz: %s",
                     len(annotations_obj.no_translated["mz"].unique(
                     )))
        logger.error("unique mz in adata_mz: %s",
                     len(adata_mz.var["mz"].unique()))


def match_adata_mz_and_annots(annotations_obj: MetabolomeAnnotations,
                              adata_mz: ad.AnnData
                              ) -> Tuple[MetabolomeAnnotations, ad.AnnData]:
    """
    Updates:
     - adata_mz
     - the following attributes of the annotation object:
        .preloaded
        .no_translated
        .df (if exists)
    making them match perfectly, by taking the m/z in the intersection only.
    """
    if not isinstance(annotations_obj, MetabolomeAnnotations):
        logger.error(
            "error: the first argument must be the annotation object!")
        raise AssertionError
    annotations_obj = annotations_obj.model_copy()  # copy of a pydantic class
    okmz = list(set(adata_mz.var['mz'].astype(float).tolist()).intersection(
        set(annotations_obj.no_translated["mz"].astype(float).unique()))
    )
    annotations_obj.preloaded = annotations_obj.preloaded.loc[
        annotations_obj.preloaded['mz'].isin(okmz), :].copy()
    annotations_obj.no_translated = annotations_obj.no_translated.loc[
        annotations_obj.no_translated['mz'].isin(okmz), :].copy()
    if annotations_obj.df.shape[0] != 0:  # if .df attr is not an empty df
        annotations_obj.df = annotations_obj.df.loc[
            annotations_obj.df['mz'].isin(okmz), :].copy()
    adata_mz = adata_mz[:, np.array(sorted(okmz), dtype=str)].copy()

    return annotations_obj, adata_mz


#  spots quality functions

def adata_mz_mark_spot_quality(adata_mz, features_classification_df,
                               max_n_tails: int = 100,
                               fdr_cutoff: float = 0.1,
                               mz_col: str = 'mz',
                               off_col: str = 'rawOffSampleProb',
                               fdr_column='fdr',
                               intensity_col='totalIntensity'
                               ) -> ad.AnnData:
    """
    :return:
    """

    def give_sum_array(adata_mz, mz_list) -> np.array:
        tmp = adata_mz[:, mz_list].copy()
        sarray = tmp.X.sum(axis=1)

        return np.array(sarray).flatten()

    df = features_classification_df[
        [mz_col, off_col, intensity_col, fdr_column]]
    if not df[off_col].max() <= 1:
        raise AssertionError(f"values of {off_col} must be in [0, 1]")
    df = df.loc[df[fdr_column] <= fdr_cutoff, :].sort_values(by=off_col)
    df_in = df.loc[df[off_col] <= 0.2, :].copy()
    n_sample = max_n_tails
    if df_in.shape[0] < n_sample:
        n_sample = df_in.shape[0]
    df_in = df_in.sort_values(by=[off_col, intensity_col, fdr_column],
                              ascending=[True, False, True])
    df_in = df_in.head(n=n_sample)
    df_out = df.loc[df[off_col] >= 0.8, :].copy()
    df_out = df_out.sort_values(by=[off_col, intensity_col, fdr_column],
                                ascending=[True, True, False])  # not eq to in
    df_out = df_out.tail(n=n_sample)
    in_sum_arr = give_sum_array(adata_mz, df_in['mz'].astype(str).tolist())
    out_sum_arr = give_sum_array(adata_mz, df_out['mz'].astype(str).tolist())
    ratios_arr = (in_sum_arr + 1) / (out_sum_arr + 1)
    adata_mz.obs = adata_mz.obs.assign(
        sum_out_tissue_quants=out_sum_arr,
        sum_in_tissue_quants=in_sum_arr,
        ratio_in_out=ratios_arr)

    return adata_mz


def compute_spot_quality(adata_mz, annotations_obj: MetabolomeAnnotations,
                         allowed_in_out_ratio: float,
                         fdr_cutoff: float) -> ad.AnnData:
    logger.info("marking spot quality with params: %s", allowed_in_out_ratio)
    adata_mz = adata_mz_mark_spot_quality(
        adata_mz, annotations_obj.mz_quality, fdr_cutoff=fdr_cutoff
    )
    adata_mz.obs['spot_quality'] = "good"
    adata_mz.obs.loc[adata_mz.obs['ratio_in_out'] < allowed_in_out_ratio,
    'spot_quality'] = "bad"
    try:
        if len(adata_mz.uns['spot_quality_colors']
               ) != len(adata_mz.obs['spot_quality'].unique()):
            del adata_mz.uns['spot_quality_colors']
    except KeyError:
        pass
    return adata_mz


def mark_weak_signal_spots(adata_mz, cutoff: float = 10) -> ad.AnnData:
    adata_mz.obs['totalIntensitySpot'] = adata_mz.X.sum(axis=1)
    adata_mz.obs.loc[
        adata_mz.obs['totalIntensitySpot'] <= cutoff, 'spot_quality'] = 'bad'

    return adata_mz
