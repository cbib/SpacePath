import os
import numpy as np
import anndata as ad
import pandas as pd
import scanpy as sc
from typing import List
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from sklearn.utils import sparsefuncs
import squidpy as sq
import zonist.processing.measures_methods as zmm

"""
get the matrix of KEGG given the annotations, and the adata_mz
"""



def provide_annotation_dict(SAMPLE, annotation_types_list):
    out_d = {}
    for anntype in annotation_types_list:
        out_d[anntype] = f"{SAMPLE}-{anntype}.csv"

    return out_d


#
# def part_a_filtersOLD(me_annots, msm_cutoff, offSample_cutoff, fdr_cutoff,
#                    SAMPLE, plot_folder):
#     print(f">> Part a): Select m/z features with parameters: "
#           f"MSM >= {msm_cutoff}, offSample <= {offSample_cutoff}, FDR <= {fdr_cutoff}")
#
#     me_annots.pull_mz_quality_df()
#     me_annots.mz_quality['quality_msm_off_fdr'] = "bad"
#     me_annots.mz_quality.loc[
#         (me_annots.mz_quality['msm'] >= msm_cutoff
#          ) & (me_annots.mz_quality['rawOffSampleProb'] <= offSample_cutoff
#          ) & (me_annots.mz_quality['fdr'] <= fdr_cutoff),
#         'quality_msm_off_fdr'] = "good"
#
#     ok_mz_keep_part1 = me_annots.mz_quality.loc[
#         me_annots.mz_quality['quality_msm_off_fdr'] == "good", "mz"].astype(
#         str).tolist()
#     text_notice = f"valid $m/z$\n n = {len(ok_mz_keep_part1)}"
#     plot_mz_given_quality_metrics(me_annots.mz_quality, SAMPLE,  # todo: take the plot part outside, and call correctly in all usages
#                                   plot_folder,
#                                   msm_cutoff, offSample_cutoff, text_notice)
#     print(f"A total of {len(set(ok_mz_keep_part1))} m/z pass MSM, offSample and FDR cutoffs.")
#
#     return me_annots, ok_mz_keep_part1


def part_a_filters(me_annots, msm_cutoff, offSample_cutoff, fdr_cutoff):
    print(f">> Part a): Select m/z features with parameters: "
          f"MSM >= {msm_cutoff}, offSample <= {offSample_cutoff}, FDR <= {fdr_cutoff}")

    me_annots.pull_mz_quality_df()
    me_annots.mz_quality['quality_msm_off_fdr'] = "bad"
    me_annots.mz_quality.loc[
        (me_annots.mz_quality['msm'] >= msm_cutoff
         ) & (me_annots.mz_quality['rawOffSampleProb'] <= offSample_cutoff
         ) & (me_annots.mz_quality['fdr'] <= fdr_cutoff),
        'quality_msm_off_fdr'] = "good"

    ok_mz_keep_part1 = me_annots.mz_quality.loc[
        me_annots.mz_quality['quality_msm_off_fdr'] == "good", "mz"].astype(
        str).tolist()

    print(f"A total of {len(set(ok_mz_keep_part1))} m/z pass MSM, offSample and FDR cutoffs.")

    return me_annots, ok_mz_keep_part1



def custom_join_left(df1, df2, keycol='mz'):
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


def part_b_filters_with_autocorr(adata_mz, ok_part1, me_annots, moran_I_cutoff, moran_padj_cutoff):
    print(f">> Part b): Classifying m/z features, by :\n"
          f"Moran's I value >= {moran_I_cutoff}, "
          f"(BH) adjusted p-value <= {moran_padj_cutoff}")
    # prepare quality merged df
    all_quality_df = adata_mz.var.copy()
    all_quality_df["mz"] = all_quality_df["mz"].astype(str).tolist()
    me_annots.mz_quality['mz'] = me_annots.mz_quality['mz'].astype(
        str).tolist()
    all_quality_df = custom_join_left(
        all_quality_df, me_annots.mz_quality[
            ['mz', 'msm', 'rawOffSampleProb', 'fdr',
             'quality_msm_off_fdr' # todo: add formula? this will go into .var
             ]], keycol='mz')

    print("(Moran's metrics are only computed for the m/z that passed the a) filters)")
    copy_mz = adata_mz[:, ok_part1].copy()   # only compute moran for the ok_part1 features
    try:
        if copy_mz.uns['moranI'].shape[0] != len(ok_part1):
            del copy_mz.uns['moranI']
    except Exception as e:
        pass

    if not "moranI" in list(copy_mz.uns.keys()):
        n_perms = 50
        if copy_mz.obs.shape[0] > 100000:
            n_perms = 3
        copy_mz = zmm.compute_moran_i(copy_mz, n_neighs=20, n_perms=n_perms)
    else:
        print(f"will use the existing .uns['moranI'] with dimension: {copy_mz.uns['moranI'].shape}")
        assert copy_mz.uns['moranI'].shape[0] == len(ok_part1), "Fatal: existing moranI dimensions do not match with length of m/z list passing a) filters. Killed! (solution: delete .uns['moranI'] dataframe and try again)"
    moran_mz_df = copy_mz.uns['moranI'].copy()
    moran_mz_df['mz'] = moran_mz_df.index.astype(str).tolist()
    all_quality_df = pd.merge(all_quality_df,
                              moran_mz_df[['mz', 'I', 'pval_sim_fdr_bh']],
                              how="left", left_on="mz", right_on="mz")
    keep_part2 = all_quality_df.loc[
        (all_quality_df['I'] >= moran_I_cutoff) & (
                all_quality_df['pval_sim_fdr_bh'] <= moran_padj_cutoff),
        "mz"].astype(str).tolist()

    all_quality_df['quality_autocorr'] = 'bad'
    all_quality_df.loc[all_quality_df['mz'].astype(str).isin(keep_part2),
     'quality_autocorr'] = "good"

    return all_quality_df, copy_mz, keep_part2, moran_mz_df










# def yield_specific_typical_pathways() -> dict:
#     biologically_chosen_pathways = {
#         "hsa00600": "Sphingolipid metabolism",
#         "hsa00030": "Pentose phosphate pathway",
#         "hsa00010": "Glycolysis Gluconeogenesis", # attention, name here has been edited to avoid "/" in the middle
#         "hsa00620" : "Pyruvate metabolism",
#         "hsa00020"  : "TCA cycle",
#         "hsa00900": "Terpenoid backbone biosynthesis",
#         "hsa04727": "GABAergic synapse"}
#     return biologically_chosen_pathways


