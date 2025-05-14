import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
from scipy.sparse import csr_matrix
from sklearn.cluster import KMeans
from sklearn.utils import sparsefuncs
import zonist.helpers as zh
from zonist.data import MetabolicPathways, MetabolomeAnnotations
from typing import Tuple, List, Union
import logging

logger = logging.getLogger(__name__)
logger = zh.reset_log_config(logger)


# related to mz annotations (Metabolome Annotations) and quantifications (andata) object


# matchers

def verification(adata_mz: ad.AnnData, me_annots: MetabolomeAnnotations) -> None:
    try:
        bool_equal = np.all(np.array(
            sorted(me_annots.no_translated["mz"].unique())) == np.array(
            sorted(adata_mz.var["mz"].tolist())))
        logger.info("perfect match between mz in me_a.df "
              f"and adata.X: => %s <=", bool_equal)
    except ValueError:
        logger.error("mz do not match between adata_mz and annotation table")
        logger.error("unique annotated mz: %s", len(me_annots.no_translated["mz"].unique(
        )))
        logger.error("unique mz in adata_mz: %s", len(adata_mz.var["mz"].unique()))


def match_adata_mz_and_annots(annotationsObj: MetabolomeAnnotations,
                              adata_mz: ad.AnnData
                              )->  Tuple[MetabolomeAnnotations, ad.AnnData]:
    """
    Updates:
     - adata_mz
     - the  following attributes of the annotations object:
        .preloaded
        .no_translated
        .df (if exists)
    making them match perfectly, by taking the m/z in the intersection only.
    """
    if not isinstance(annotationsObj, MetabolomeAnnotations):
        logger.error("error: the first argument must be the annotation object!")
        raise AssertionError
    me_annots = annotationsObj.model_copy()  # copy of a pydantic class
    okmz = list(set(adata_mz.var['mz'].astype(float).tolist()).intersection(
        set(me_annots.no_translated["mz"].astype(float).unique()))
    )
    me_annots.preloaded = me_annots.preloaded.loc[
        me_annots.preloaded['mz'].isin(okmz),:].copy()
    me_annots.no_translated = me_annots.no_translated.loc[
        me_annots.no_translated['mz'].isin(okmz),:].copy()
    try:
        me_annots.df = me_annots.df.loc[
                                  me_annots.df['mz'].isin(okmz), :].copy()
    except:
        pass
    adata_mz = adata_mz[:, np.array(sorted(okmz), dtype=str)].copy()

    return me_annots, adata_mz


#  spots quality functions

def adata_mz_mark_spot_quality(adata_mz, features_classification_df,
        max_n_tails: int = 100, fdr_cutoff: float = 0.1,
        mz_col: str='mz', off_col: str='rawOffSampleProb', fdr_column= 'fdr',
                               intensity_col='totalIntensity'
)-> ad.AnnData:
    """
    :return:
    """
    def give_sum_array(adata_mz, mz_list):
        tmp = adata_mz[:, mz_list].copy()
        sarray = tmp.X.sum(axis=1)

        return np.array(sarray).flatten()

    df = features_classification_df[[mz_col, off_col, intensity_col, fdr_column]]
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
                                ascending=[True, True, False]) # not eq to in
    df_out = df_out.tail(n=n_sample)
    in_sum_arr = give_sum_array(adata_mz, df_in['mz'].astype(str).tolist())
    out_sum_arr = give_sum_array(adata_mz, df_out['mz'].astype(str).tolist())
    ratios_arr = (in_sum_arr + 1) / (out_sum_arr + 1)
    adata_mz.obs = adata_mz.obs.assign(
        sum_out_tissue_quants=in_sum_arr,
        sum_in_tissue_quants=out_sum_arr,
        ratio_in_out=ratios_arr)

    return adata_mz


def compute_spot_quality(adata_mz, me_annots: MetabolomeAnnotations,
                         allowed_in_out_ratio: float, fdr_cutoff: float):
    logger.info("marking spot quality with params: %s", allowed_in_out_ratio)
    adata_mz = adata_mz_mark_spot_quality(
        adata_mz, me_annots.mz_quality, fdr_cutoff=fdr_cutoff
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


def mark_weak_signal_spots(adata_mz, cutoff: float = 10):
    adata_mz.obs['totalIntensitySpot'] =  adata_mz.X.sum(axis=1)
    adata_mz.obs.loc[adata_mz.obs['totalIntensitySpot'] <= cutoff, 'spot_quality'] = 'bad'

    return adata_mz


# #  normalization and autocorrelation


def tic_normalization(intensity_matrix: Union[csr_matrix, np.array],
                      obs_names: List[str], mz_names: List[float]) -> np.matrix:
    """
    Perform Total Ion Count (TIC) normalization on a raw 2D intensity matrix
    of m/z features
    :param intensity_matrix: rows are spots, columns are mz features
    :return:  TIC normalized intensity matrix
    """
    logger.info("Applying Total Ion Count (TIC) normalization: processes the m/z "
          "features, across spots")
    assert intensity_matrix.shape == (
        len(obs_names), len(mz_names)), f"wrong intensity array dimensions"
    total_ion_count = intensity_matrix.sum(axis=1)
    scale_factor = np.median(np.array(total_ion_count))
    # Avoid division by zero
    total_ion_count = np.array(total_ion_count) + 1  # avoid any zero error div
    # Normalize the intensity matrix
    normalized_matrix = scale_factor * (intensity_matrix / total_ion_count)

    return normalized_matrix


def compute_moran_i(adata, n_neighs: int=20, n_perms: int=5):
    """note: to run Moran's I, through squidpy, the `spatial_connectivities`
     in `adata.obsp` are required, so spatial neighbors must be run first"""
    logger.info("computing spatial neighbors")
    logger.info('with n_neighs: %s', n_neighs)
    sq.gr.spatial_neighbors(adata, n_neighs=n_neighs)
    logger.info("Computing Moran's I")
    sq.gr.spatial_autocorr(adata, mode="moran", n_perms=n_perms, n_jobs=1, seed=0)

    return adata

# phases of classification for m/z

def part_a_filter(me_annots, msm_cutoff, offSample_cutoff, fdr_cutoff):
    me_annots.pull_mz_quality_df()
    me_annots.mz_quality['quality_msm_off_fdr'] = "bad"
    me_annots.mz_quality.loc[
        (me_annots.mz_quality[me_annots.msm_column] >= msm_cutoff
         ) & (me_annots.mz_quality[me_annots.offsample_column] <= offSample_cutoff  #
         ) & (me_annots.mz_quality[me_annots.fdr_column] <= fdr_cutoff),
        'quality_msm_off_fdr'] = "good"

    mz_to_keep_part1 = me_annots.mz_quality.loc[
        me_annots.mz_quality['quality_msm_off_fdr'] == "good", "mz"].astype(
        str).tolist()

    return me_annots, mz_to_keep_part1


def custom_join_left(df1, df2, keycol: str='mz') -> Union[pd.DataFrame, None]:
    df1 = df1.copy()
    df2 = df2.copy()
    df2[keycol] = df2[keycol].astype(str)
    df2 = df2.drop_duplicates(subset=keycol)
    df2 = df2.set_index(keycol)
    try:
        df2 = df2.loc[df1[keycol].astype(str).tolist(), :]  # same row order as 1
        # np.all(df2.index.astype(str).to_numpy(
        # ) == df1[keycol].astype(str).to_numpy())
        for col in df2.columns:
            df1[col] = df2[col]
        return df1
    except KeyError:
        print("fatal (no match), custom_join_left impossible!")
        return None
    except Exception as e:
        print(e, "fatal (no match), custom_join_left impossible!")
        return None


def part_b_filter_with_autocorr(adata_mz, ok_part1, me_annots,
                                moran_I_cutoff, moran_padj_cutoff, n_neighs, n_perms):  #
    # prepare quality merged df
    all_quality_df = adata_mz.var.copy()
    all_quality_df["mz"] = all_quality_df["mz"].astype(str).tolist()
    me_annots.mz_quality['mz'] = me_annots.mz_quality['mz'].astype(
        str).tolist()
    all_quality_df = custom_join_left(
        all_quality_df, me_annots.mz_quality[
            ['mz', me_annots.msm_column, me_annots.offsample_column,
             me_annots.fdr_column, me_annots.formula_col, "quality_msm_off_fdr"]], keycol='mz')

    copy_mz = adata_mz[:, ok_part1].copy()   # only compute moran for the ok_part1 features
    try:
        if copy_mz.uns['moranI'].shape[0] != len(ok_part1):
            del copy_mz.uns['moranI']
    except Exception as e:
        pass

    if not "moranI" in list(copy_mz.uns.keys()) :
        copy_mz = compute_moran_i(copy_mz, n_neighs=n_neighs, n_perms=n_perms)
    else:
        print(f"will use the existing .uns['moranI'] with dimension: {copy_mz.uns['moranI'].shape}")
        assert copy_mz.uns['moranI'].shape[0] == len(ok_part1), "Fatal: existing moranI dimensions do not match with length of m/z list passing a) filters. Killed! (solution: delete .uns['moranI'] dataframe and try again)"
    moran_mz_df = copy_mz.uns['moranI'].copy()
    moran_mz_df['mz'] = moran_mz_df.index.astype(str).tolist()
    all_quality_df = pd.merge(all_quality_df,
                              moran_mz_df[['mz', 'I', 'pval_sim_fdr_bh']],
                              how="left", left_on="mz", right_on="mz")
    mz_to_keep_part2 = all_quality_df.loc[
        (all_quality_df['I'] >= moran_I_cutoff) & (
                all_quality_df['pval_sim_fdr_bh'] <= moran_padj_cutoff),
        "mz"].astype(str).tolist()

    all_quality_df['quality_autocorr'] = 'bad'
    all_quality_df.loc[all_quality_df['mz'].astype(str).isin(mz_to_keep_part2),
     'quality_autocorr'] = "good"

    return all_quality_df, mz_to_keep_part2, moran_mz_df


def adata_mz_insert_full_quality(adata_mz, moran_mz_df,
                                 all_quality_df):
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






