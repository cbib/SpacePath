import os
from typing import Dict, List

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
from scipy.sparse import csr_matrix
from sklearn.cluster import KMeans
from sklearn.utils import sparsefuncs
from zonist.data import MetabolicPathways, MetabolomeAnnotations
import zonist.helpers as zh


def closest_match_mz(mz_list1, mz_list2, tolerance=0.002):
    """ mz_list1 are the mzs as they appear in annotations
        mz_list2 are the mzs as they appear in intensities matrix"""
    df = pd.DataFrame(
        {"v1_f": np.array(mz_list1, dtype=float),
         "v2_f": [np.nan for i in mz_list1],
         "v1_str": mz_list1,
         "v2_str": [np.nan for i in mz_list1]}
    )
    v2_floats = np.array(mz_list2, dtype=float)
    for i, row in df.iterrows():
        tmp = abs(v2_floats - row["v1_f"])
        which_min = np.argmin(tmp)
        if tmp[which_min] <= tolerance:
            df.iloc[i, 1] = v2_floats[which_min]
            df.iloc[i, 3] = str(v2_floats[which_min])
    print("pre", df.shape)
    df = df.dropna()
    print("after", df.shape)
    return df


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


def tic_normalization(intensity_matrix, obs_names, mz_names):
    """
    Perform Total Ion Count (TIC) normalization on a raw 2D intensity matrix
    :param intensity_matrix: rows are spots, columns are mz features
    :return:  TIC normalized intensity matrix
    """
    print("Applying Total Ion Count (TIC) normalization")
    assert intensity_matrix.shape == (
        len(obs_names), len(mz_names)), f"wrong intensity array dimensions"
    total_ion_count = intensity_matrix.sum(axis=1)
    scale_factor = np.median(np.array(total_ion_count))
    # Avoid division by zero
    total_ion_count = np.array(total_ion_count) + 1  # avoid any zero error div
    # Normalize the intensity matrix
    normalized_matrix = scale_factor * (intensity_matrix / total_ion_count)

    return normalized_matrix


def mold_data_for_plot__task1(pw_id, scores_adata, adata):
    """
        First phase of prep for plot with table.
    :param pw_id: the current pathway identifier
    :param scores_adata: the pathways identifiers and scores object
    :param adata: the KEGG compound identifiers and quantities object
    :return: df with KEGG id, metabolite name and mz for the current pw_id
    """
    df_mz_metabolites__curr_pw  = scores_adata.uns['current_pathways_df'].loc[
        scores_adata.uns['current_pathways_df']['pathway_id'] == pw_id,
        ["mz", "KEGG"]
    ]
    df_mz_metabolites__curr_pw = df_mz_metabolites__curr_pw.drop_duplicates()
    df2 = adata.var[["ID", "moleculeNames"]].copy()
    df2.columns = ["KEGG", "metabolite name"]
    df_mz_metabolites__curr_pw = df_mz_metabolites__curr_pw.join(
        df2, how="left", on="KEGG", lsuffix="_x", rsuffix="_y"
    )
    df_mz_metabolites__curr_pw = df_mz_metabolites__curr_pw.drop(
        columns=["KEGG_x", "KEGG_y"]
    )
    return df_mz_metabolites__curr_pw



def mold_data_for_plot__task2(df_mz_metabolites__curr_pw, i_mz, coords_mz_df, adata_mz):
    """
        Second phase of prep for plot with table.

    :param df_mz_metabolites__curr_pw:
    :param i_mz:
    :param adata_mz:
    :return:
    """
    sub_df = df_mz_metabolites__curr_pw.loc[
             df_mz_metabolites__curr_pw['mz'] == i_mz, :]
    sub_df = sub_df.drop(columns=["mz"])
    sub_df = sub_df.drop_duplicates()
    sub_df['metabolite name'] = sub_df['metabolite name'].str.slice(0, 25)
    tmp = coords_mz_df.copy()
    tmp[i_mz] = np.array(adata_mz[:, str(i_mz)].X.todense()).flatten()

    return sub_df, tmp



def score_as_adata(kpca_scores, adata,
                   #pathwaysObj: MetabolicPathways,
                   pathways_ids_df,
                   pathways_piled,
                   original_pathways_dict):
    score_adata = ad.AnnData(
        X=csr_matrix(kpca_scores),
        var=pd.DataFrame(
            {'pathway_id': kpca_scores.columns.tolist()},
            index=kpca_scores.columns.tolist()),
        obs=adata.obs,
        obsm=adata.obsm,
        uns=adata.uns
    )
    tmp = pathways_ids_df #pathwaysObj.pathways_ids.copy()
    tmp.index = tmp['pathway_id']
    tmp = tmp.loc[score_adata.var['pathway_id'],:]
    score_adata.var = score_adata.var.assign(
        pathway_name=tmp['pathway_name'].tolist())
    score_adata.uns['former_pathways_dict'] = original_pathways_dict
    score_adata.uns['current_pathways_df'] = pathways_piled #pathwaysObj.piled  # has m/z, KEGG and pathway id columns

    return score_adata


def compute_mz_binary_df(adata_mz: ad.AnnData):
    X_presence_absence_mask = adata_mz.X.todense() > 0
    X_mz_binary_df = pd.DataFrame(
        np.array(X_presence_absence_mask, dtype=int),
        columns= adata_mz.var['mz'].to_numpy(),
        index=adata_mz.obs.index.tolist())

    return X_mz_binary_df


def compute_metabolites_probs_df(adata_mz: ad.AnnData, annots_OKdf):
    X_mz_binary_df = compute_mz_binary_df(adata_mz)
    preX_probs = []
    col_names_kegg = []
    for a_mz in X_mz_binary_df.columns:
        i_arr = X_mz_binary_df[a_mz].to_numpy()
        the_kegg_here = annots_OKdf.loc[
            annots_OKdf['mz'] == a_mz, "KEGG"].tolist()
        n = len(the_kegg_here)
        for kid in the_kegg_here:
            col_names_kegg.append(kid)
            initial_prob = 1 / n
            # res = i_arr * np.array([initial_prob for v in range(len(i_arr))])
            res = i_arr * np.repeat(initial_prob, len(i_arr))
            preX_probs.append(np.array(res, dtype=np.float32))

    X_metabolites_probs_df = pd.DataFrame(np.array(preX_probs).T,
                                          columns=col_names_kegg,
                                          index=X_mz_binary_df.index)
    return X_metabolites_probs_df


def attach_annots_to_adata(adata, annots_df):
    """"""
    annots_df = annots_df.drop_duplicates()
    assert len(set(adata.var['ID'].tolist()).intersection(
        set(annots_df['KEGG'].tolist()))) == len(
        set(annots_df['KEGG'].tolist())), ("Error, adata.var ID column must "
                                           "match KEGG column of annotations")
    annots_df.index = annots_df['KEGG'].tolist()
    df_ok = annots_df.loc[adata.var['ID'], :]  # rows in same order than in adata.var
    expected_index = adata.var.index.tolist()
    adata.var = pd.merge(adata.var, df_ok,  how='left' , left_on='ID', right_on='KEGG')
    #adata.var.index = adata.var['ID'].tolist()
    adata.var.set_index('ID',
                        drop=False, # whether to drop the ID column
                        inplace=True)
    assert np.all(np.array(expected_index) == np.array(
        adata.var.index.tolist())), ("all messed up "
                                     "when attaching annots to adata")

    return adata


def build_kegg_resolved_adata(adata_mz, X_metabolites_probs_df, annots_df):
    """
    Given the adata_mz and the annots_df, returns a new adata object
    """
    tmp_X = dict()
    for a_mz in adata_mz.var['mz'].to_numpy():
        the_kegg_here = annots_df.loc[
            annots_df['mz'] == a_mz, "KEGG"].tolist()
        sub_adata = adata_mz[:, str(a_mz)].copy()
        i_mz_arr = np.array(sub_adata.X.todense().T).flatten()
        for a_kg in the_kegg_here:
            curr_prod_arr = \
                X_metabolites_probs_df[a_kg].to_numpy() * i_mz_arr
            tmp_X[a_kg] = np.array(curr_prod_arr, dtype=np.float32)
    # end for
    kegg_list = list(tmp_X.keys())
    tmp_X = np.array(pd.DataFrame(tmp_X), dtype=np.float32)
    adata = ad.AnnData(X=csr_matrix(tmp_X),
                       obs=adata_mz.obs,
                       var=pd.DataFrame({"ID": kegg_list},
                                        index=kegg_list),
                       obsm=adata_mz.obsm)
    adata.uns['spatial'] = {'my_slice_name': {
        'scalefactors': {'fiducial_diameter_fullres': 100,
                         'spot_diameter_fullres': 1,
                         'tissue_hires_scalef': 0.5,
                         'tissue_lowres_scalef': 0.5}}}

    adata.layers["X_probabilities"] = X_metabolites_probs_df

    adata = attach_annots_to_adata(adata, annots_df) ## new !

    return adata


def adata_mz_mark_spot_quality(adata_mz, features_classification_df,
        max_n_tails: int = 100, fdr_cutoff: float = 0.1,
        mz_col: str='mz', off_col: str='rawOffSampleProb', fdr_column= 'fdr',
                               intensity_col='totalIntensity'
):
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

#

def match_adata_mz_and_annots(annotationsObj: MetabolomeAnnotations, adata_mz: ad.AnnData):  # zd.MetabolomeAnnotations
    """
    Updates:
     - adata_mz
     - the two following attributes of the annotations object:
        .preloaded
        .no_translated
    making them match perfectly, in terms of  m/z features.
    """
    # todo: add a check of adata_mz to verify that only good quality features are present in adata_mz (checks if all values in adata_mz.var['final_quality'] are "good")
    if not isinstance(annotationsObj, MetabolomeAnnotations):
        raise AssertionError("error: the first argument must be the annotation object!")
    me_annots = annotationsObj.model_copy()  # copy is deprecated for pydantic classes
    okmz = list(set(adata_mz.var['mz'].astype(float).tolist()).intersection(
        set(me_annots.no_translated["mz"].astype(float).unique()))
    )
    me_annots.preloaded = me_annots.preloaded.loc[
        me_annots.preloaded['mz'].isin(okmz),:]
    me_annots.no_translated = me_annots.no_translated.loc[
        me_annots.no_translated['mz'].isin(okmz),:]
    try:
        me_annots.df = me_annots.df.loc[
                                  me_annots.df['mz'].isin(okmz), :]
    except:
        pass
    adata_mz = adata_mz[:, np.array(sorted(okmz), dtype=str)]

    return me_annots, adata_mz


def update_annotation_tomatch_adatamz(annotationsObj: MetabolomeAnnotations, adata_mz: ad.AnnData):  # zd.MetabolomeAnnotations
    # update the annotations to match with the adata_mz mz features
    """
    Updates two attributes of the annotations object:
        .preloaded
        .no_translated
      making them match with the mz features of adata_mz object
    """
    me_annots = annotationsObj.model_copy()  # copy is deprecated for pydantic classes
    okmz = adata_mz.var['mz'].astype(float).tolist()
    me_annots.preloaded = me_annots.preloaded.loc[
        me_annots.preloaded['mz'].isin(okmz),:]
    me_annots.no_translated = me_annots.no_translated.loc[
        me_annots.no_translated['mz'].isin(okmz),:]

    return me_annots


def verification(adata_mz, me_annots: MetabolomeAnnotations):
    print("verifying mz match")
    try:
        bool_equal = np.all(np.array(
            sorted(me_annots.no_translated["mz"].unique())) == np.array(
            sorted(adata_mz.var["mz"].tolist())))
        print(f"perfect match between mz in me_a.df "
              f"and adata.X: => {bool_equal} <=")
    except ValueError:
        print("mz do not match between adata_mz and annotation table")
        print("unique annotated mz: ", len(me_annots.no_translated["mz"].unique(
        )))
        print("unique mz in adata_mz: ", len(adata_mz.var["mz"].unique()))


# def get_rid_of_poor_var_mz(adata_mz): # todo: del
#     print("Excluding mz features with extremely low variance")
#     means_vars, variance_computed = sparsefuncs.mean_variance_axis(
#         adata_mz.X, axis=0)
#     #logger.info("variance info: %s", pd.DataFrame(variance_computed).describe())
#     min_accept_tile = np.percentile(variance_computed, 5)
#     #low_variance = variance_computed[variance_computed >= min_accept_centile]
#     adata_mz.var["computed_var"] = variance_computed
#     good_mz = adata_mz.var.loc[adata_mz.var["computed_var"
#                      ] > min_accept_tile, "mz"].to_numpy()
#     adata_mz = adata_mz[:, [str(i) for i in good_mz]]
#
#     return adata_mz


def spatial_standard_computations(
        adata, k: int=8,  pca_ncomps: int=12, as_log: bool=True,
        top_variable_features_fraction: float=0.5):
    if (adata.var.shape[0] > pca_ncomps) and (adata.var.shape[0] > 20):
        if as_log:
            adata.layers['X_raw'] = adata.X.copy()
            sc.pp.log1p(adata)
        # end if as_log
        topn = int(adata.var.shape[0] * top_variable_features_fraction)
        try:
            sc.pp.highly_variable_genes(
                adata, flavor="seurat_v3",
                check_values=False, n_top_genes=topn)
        except ValueError:
            print("Could not compute highly variable feats with scanpy. Computing with sparsefuncs")
            try:
                means_vars, variance_computed = sparsefuncs.mean_variance_axis(
                    adata.X, axis=0)
                adata.var['variance'] = variance_computed
                adata.var['tmp_index'] = range(adata.var.shape[0])
                tmpdf = adata.var.copy()
                tmpdf = tmpdf.sort_values(by='variance', ascending=False)
                tmpdf = tmpdf.head(topn).copy()
                adata.var['highly_variable'] = False
                adata.var.loc[adata.var['tmp_index'].isin(tmpdf['tmp_index'].tolist()), 'highly_variable'] = True
                adata.var = adata.var.drop(columns=['tmp_index'])
            except ValueError:
                print("impossible to compute highly variable features, marked all as True")
                adata.var['highly_variable'] = True
        print("computing pca")
        sc.pp.pca(adata, n_comps=pca_ncomps)  # use_highly_variable is true by default, but .var['highly_variable'] col was True in all values
        print("computing kmeans")
        kmeans = KMeans(n_clusters=k, n_init=10, random_state=0).fit(
            adata.obsm['X_pca'])
    else:
        print("skipping PCA: too few features")
        X = np.array(adata.X.todense())
        if as_log:
            X = np.log1p(X)
        kmeans = KMeans(n_clusters=k, n_init=10, random_state=0
                        ).fit(X)
    adata.obs[f'kmeans{k}'] = kmeans.labels_.astype(str)

    return adata


def compute_moran_i(adata, n_neighs=10, n_perms=5):
    """note: to run Moran's I, through squidpy, the `spatial_connectivities`
     in `adata.obsp` are required, so spatial neighbors must be run first"""
    print("computing spatial neighbors")
    print('with n_neighs: ', n_neighs)
    sq.gr.spatial_neighbors(adata, n_neighs=n_neighs)
    print("Computing Moran's I")
    sq.gr.spatial_autocorr(adata, mode="moran", n_perms=n_perms, n_jobs=1, seed=0)

    return adata


def get_target_pathways_of_one_mz(my_mz: str, subdf, pathways_dict, pws_mz_dict):
    """
    :param: my_mz: current mz
    :param subdf: represents only *one mz* dataframe (cols: pathway_id, KEGG),
       being KEGG the metabolite identifier.
    :param: pathways_dict: keys are pathway ids, values are all the kegg ids
    :param: pws_mz_dict: keys are pathway ids, values are all the mzs
    :return: a dict where the keys are the compared pairs of pathways, e.g.
       {('hsa0230', 'hsa4976'): {'hsa0230': ['C00212'], 'hsa4976': ['C07210']}

    Explanation:

    For the given distinct KEGG ids and pathway ids in the dataframe,
    first constructs a dict (of the current mz) of sets of KEGG ids for each pathway.
    Then, compares each possible pair (k, h) of pathways from the 'mz' dict,
    to check if there are KEGG ids that are exclusive to k (same for h).
    Also takes into account the total KEGG ids of k and h pathways (not ony the current mz).

    If (k, h) pathways had an empty intersection related to the kegg of the current mz,
               !?                 if there are, in each pathway, a total of kegg ids that is superior or equal to the
                 ?               n+1 kegg ids of the current mz,
    and if mzs other than current mz are in each pathway,
    then both (k, h) pathways go into target dict.
    """
    adict_curr_mz: Dict[str, set] = dict()
    pathways_to_target = dict()
    for i, row in subdf.iterrows():
        if row["pathway_id"] not in adict_curr_mz.keys():
            adict_curr_mz[row["pathway_id"]] = {row["KEGG"]}
        adict_curr_mz[row["pathway_id"]].update([row["KEGG"]])
    for k in adict_curr_mz.keys():  # k is a pathway id
        for h in adict_curr_mz.keys():  # h is also a pathway id
            if k != h:
                other_mz_pw1 = set(pws_mz_dict[k]) - set([my_mz])
                other_mz_pw2 = set(pws_mz_dict[h]) - set([my_mz])
                a: set = adict_curr_mz[k] - adict_curr_mz[h]  # a is a set of kegg ids
                b: set = adict_curr_mz[h] - adict_curr_mz[k]  # b is also a set of kegg ids
                intersect_k_h :set = adict_curr_mz[h].intersection(adict_curr_mz[k])
                if (len(a) >= 1) and (len(b) >= 1) and (len(intersect_k_h) == 0):
                    if ((len(other_mz_pw1) >= 1) and (len(other_mz_pw2) >= 2) or (
                            len(other_mz_pw1) >= 2) and (len(other_mz_pw2) >= 1)
                    )and (len(other_mz_pw1.intersection(other_mz_pw2)) == 0):
                        pathways_to_target[(k,h)] = {k: list(a), h: list(b)}

    return pathways_to_target



def message_factory_1(total_nb, kegg_oneormore__one_pathway,
                       one_kegg_multiple_pathways,  multiple_kegg_shared,
                      susceptible_of_refinement) -> None:
    print("Total mz: ", total_nb)
    print("kegg_oneormore__one_pathway: ", kegg_oneormore__one_pathway)
    print("one_kegg_multiple_pathways: ", one_kegg_multiple_pathways)
    print("multiple_kegg_shared: ", multiple_kegg_shared)
    # print("susceptible_of_refinement: ", susceptible_of_refinement)


def find_features_for_refinement(mz_met_pathway_piled: pd.DataFrame):
    """
    :param mz_met_pathway_piled: The df with the correspondences of mz,
        metabolites KEGG ids, and pathways KEGG ids.
    :return: a dictionary with the mz and pairs of pathways that
        are candidates for metabolites assignments refinement, e.g.
        '555.5555999' : {('hsa0230', 'hsa4976'): {
                          'hsa0230': ['C00212'], 'hsa4976': ['C07210']}
    """
    pathways_dict = zh.pws_piled_to_dict(mz_met_pathway_piled,
                                         'pathway_id', 'KEGG')
    pws_mz_dict = zh.pws_piled_to_dict(mz_met_pathway_piled,
                                         'pathway_id', 'mz')
    kegg_oneormore__one_pathway = 0
    one_kegg_multiple_pathways = 0
    multiple_kegg_shared = 0
    susceptible_of_refinement = 0
    final_dict: Dict[str, Dict[tuple, Dict[str, list]]] = dict()
    for mz in mz_met_pathway_piled["mz"].unique():
        subdf = mz_met_pathway_piled.loc[mz_met_pathway_piled["mz"] == mz, :]
        subdf = subdf.drop(columns=["mz"])
        subdf = subdf.drop_duplicates()
        if subdf.shape[0] == 1:
            kegg_oneormore__one_pathway += 1
        elif (len(subdf["KEGG"].unique()) == 1) and (
                len(subdf["pathway_id"].unique()) > 1):
           one_kegg_multiple_pathways += 1
        elif (len(subdf["KEGG"].unique()) > 1) and (
                len(subdf["pathway_id"].unique()) == 1):
            kegg_oneormore__one_pathway += 1
        elif (len(subdf["KEGG"].unique()) > 1) and (
                len(subdf["pathway_id"].unique()) > 1):
            subdf = subdf.drop_duplicates("KEGG")
            pws_ok = get_target_pathways_of_one_mz(str(mz), subdf, pathways_dict, pws_mz_dict)
            if len(pws_ok) > 0:
                susceptible_of_refinement += 1
                final_dict[str(mz)] = pws_ok
            else:  # not real multi-kegg multi-pathway, see target_mz_pathway
                multiple_kegg_shared += 1

    message_factory_1(len(mz_met_pathway_piled["mz"].unique()), kegg_oneormore__one_pathway,
                      one_kegg_multiple_pathways, multiple_kegg_shared,
                      susceptible_of_refinement)

    return final_dict






def inner_fun_probs_updater(pw_score: pd.DataFrame, spot1: str, spot2: str, pw1: str, pw2: str,
                            metabolites1: List[str], metabolites2: List[str],
                            met_prob: pd.DataFrame,
                            stealing_fraction: float,
                            evaluated_keggs_curr_spot_curr_mz: List[str],
                            updated_probs_spot1: Dict[str,float]):
    """
    :param pw_score: e.g. :
             hsa00100  hsa00120
        x44_y0  0.576450 -0.056398
        x43_y0 -0.287701  0.670142

    :param spot1:
    :param spot2:
    :param pw1:
    :param pw2:
    :param metabolites1:
    :param metabolites2:
    :param met_prob: e.g.:
                    C01673  C17339  C17331  C17336  C17332
        x44_y0     0.2     0.2     0.2     0.2     0.2
        x43_y0     0.0     0.0     0.0     0.0     0.0


    :param stealing_fraction:
    :param evaluated_keggs_curr_spot_curr_mz:
    :param updated_probs_spot1:
    :return:

    note: as metabolites1 and metabolits2 are often of different lenghts,
    the stealing is done in the pooled values

    """

    if pw_score.loc[spot2, pw1] >= pw_score.loc[spot1, pw1]:
        if pw_score.loc[spot1, pw1] > pw_score.loc[spot1, pw2]:
            steal_prop_of_mt2 = met_prob.loc[
                                    spot1, metabolites2].to_numpy() * stealing_fraction  # * 0.5 or 1
            steal_pool = steal_prop_of_mt2.sum() / len(met_prob.loc[spot1, metabolites1].to_numpy())
            new_pr_2s = met_prob.loc[spot1, metabolites2].to_numpy() - steal_prop_of_mt2
            new_pr_1s = met_prob.loc[spot1, metabolites1].to_numpy() + steal_pool

            # print("the old was:", met_prob.loc[spot1, metabolites2].to_numpy())
            for a in range(len(metabolites1)):
                if metabolites1[a] not in evaluated_keggs_curr_spot_curr_mz:
                    updated_probs_spot1[metabolites1[a]] = new_pr_1s[a]
                    evaluated_keggs_curr_spot_curr_mz.append(metabolites1[a])
            for b in range(len(metabolites2)):
                if metabolites2[b] not in evaluated_keggs_curr_spot_curr_mz:
                    updated_probs_spot1[metabolites2[b]] = new_pr_2s[b]
                    evaluated_keggs_curr_spot_curr_mz.append(metabolites2[b])
            # print("the new group1 is: ", new_pr_1s, metabolites1)
            # print("the new group2 is: ", new_pr_2s, metabolites2)
            # print("")

    return updated_probs_spot1, evaluated_keggs_curr_spot_curr_mz


def probabilities_updater(refinement_mz_pws_dict, pathways_scores_df,
                          X_metabolites_probs_df, spot1, spot2, stealing_fraction: float):
    the_mz_for_refinement_list = list(refinement_mz_pws_dict.keys())
    updated_probs_spot1 = dict()
    for i_mz in the_mz_for_refinement_list:
        pws_tuples_l = list(refinement_mz_pws_dict[str(i_mz)].keys())
        x = 0
        evaluated_keggs_curr_spot_curr_mz = list()
        while (x < len(pws_tuples_l)):
            pws_tuple = pws_tuples_l[x]
            pw1 = pws_tuple[0]
            pw2 = pws_tuple[1]
            metabolites1 = refinement_mz_pws_dict[str(i_mz)][pws_tuple][pw1]
            metabolites2 = refinement_mz_pws_dict[str(i_mz)][pws_tuple][pw2]
            pw_score = pathways_scores_df.loc[[spot1, spot2], [pw1, pw2]]
            met_prob = X_metabolites_probs_df.loc[
                [spot1, spot2], metabolites1 + metabolites2]

            # in spot2: if pw1 is largely superior to pw2 (at least twice)
            if (pw_score.loc[spot2, pw1] >= (pw_score.loc[spot2, pw2] * 2)) and (
                    (met_prob.loc[spot1, metabolites1].sum() + met_prob.loc[
                spot1, metabolites2].sum()) > 0):

                updated_probs_spot1, evaluated_keggs_curr_spot_curr_mz = inner_fun_probs_updater(
                    pw_score, spot1, spot2, pw1, pw2,
                    metabolites1, metabolites2, met_prob, stealing_fraction,
                    evaluated_keggs_curr_spot_curr_mz, updated_probs_spot1)

            x += 1
        # end while

    return updated_probs_spot1



def build_coords_df(adata):
    df = pd.DataFrame(data=adata.obsm['spatial'],
                 columns=['x', 'y'],
                 index=adata.obs.index.tolist())
    return df


if __name__ == "__main__":
    intensity_matrix = np.array([[100.44, 400.55, 700.33],
                                 [200.2, 500.7, 800.5],
                                 [0., 0., 0.],
                                 [300.4, 600.1,
                                  900.3]])  # 3 features(peaks), 4 spots (each spectrum is a row (spot))

    intensity_matrix = csr_matrix(intensity_matrix)

    # Perform TIC normalization
    normalized_matrix = tic_normalization(intensity_matrix,
                                          ["s1", "s2", "s3", "s4"],
                                          [84.66, 88.31, 150.55])

    print("Original Intensity Matrix:")
    print(intensity_matrix.todense())  # , shape=(3,4)))
    print("\nTIC Normalized Intensity Matrix:")
    print(normalized_matrix)
