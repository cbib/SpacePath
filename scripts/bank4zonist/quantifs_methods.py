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



# todo: a new module named assignments_refine in :  zonist.src.zonist.processing



def task1(pw_id, scores_adata, adata):
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



def task2(df_mz_metabolites__curr_pw, i_mz, coords_mz_df, adata_mz):
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

#
def adata_mz_mark_spot_quality(# for future: put this as aux_adata_methods because is tooo slow for big nb of spots
        adata_mz, features_classification_df, cutoff_out_prob: float=0.87,
        tolerated_fraction: float=0.5,
        mz_col: str='mz', off_col: str='rawOffSampleProb'):
    """
    Marks good and bad pixels in two steps:
    a)  Defines a quality ratio for each spot: it is a ratio computed as
    the sum of the quantities (signal) of the 'out-of-tissue' m/z (mz) features
    over the total of the quantities of all mz, for the current spot.

    Note : The out-of-tissue decision is defined by cutoff_out_prob which is
           the cutoff of the probability of an m/z being outside the tissue,
           computed by Metaspace.

    The smaller the ratio, the better the quality of the spot. e.g.:
    ratio_out_tissue__total = 0.5 means that sum of signal of out-of-tissue
    mzs represent the half of the total sum of signal of all mzs.

    b) The "good" pixel score is inferior or equal to the tolerated_fraction.

    Results stored in adata_mz.obs
    """
    df = features_classification_df[[mz_col, off_col]]
    if not df[off_col].max() <= 1:
        raise AssertionError(f"values of {off_col} must be in [0, 1]")
    df_tr = df.T
    df_tr.index = [mz_col, off_col]
    df_tr.columns = df[mz_col].tolist()
    df_tr = df_tr[adata_mz.var[mz_col].tolist()]  # columns in appropriate order
    allright = np.all(np.array(df_tr.columns.tolist()) == adata_mz.var['mz'].to_numpy())
    sum_out_tissue_quantities_arr = np.zeros(shape=adata_mz.obs.shape[0])
    total_arr = np.zeros(shape=adata_mz.obs.shape[0])
    if allright:
        for i in range(0, adata_mz.obs.shape[0]):
            x = adata_mz.obs.index.tolist()[i]
            tmp = pd.DataFrame({
                "quantifications": np.array(adata_mz[x, :].X.todense()).flatten(),
                off_col : df_tr.loc[off_col, :].to_numpy()
            })
            sum_out_tissue_quantities_arr[i] = tmp.loc[
                tmp[off_col] >= cutoff_out_prob,
                "quantifications"].sum()
            total_arr[i] =  tmp["quantifications"].sum()
        # end for

        ratios_arr = (sum_out_tissue_quantities_arr + 1) / (total_arr + 1)  # avoid zero division error
        adata_mz.obs = adata_mz.obs.assign(
            sum_out_tissue_quants=sum_out_tissue_quantities_arr,
            total_quants=total_arr,
            ratio_out_tissue__total=ratios_arr)
        adata_mz.obs['spot_quality'] = "good"
        adata_mz.obs.loc[adata_mz.obs['ratio_out_tissue__total'] > tolerated_fraction,
            'spot_quality'] = "bad"

    else:
        print("impossible to add the spot quality")

    return adata_mz


def adata_mz_filter_negligible_abundance(adata_mz, spots_non_zero=5):   #TODO change to 10?
    # TODO: evaluate the nb of zeroes in entire matrix, if too few , fix a minval for converting to zero all values <= minval ??
    non_zero_counts = adata_mz.X.getnnz(axis=0)
    mz_df = adata_mz.var.copy()
    mz_df['bool_arr'] = [1 if x >= spots_non_zero
                         else 0 for x in non_zero_counts]
    adata_mz = adata_mz[:,
               mz_df.loc[mz_df['bool_arr'] == 1, "mz"].astype(str).tolist()]
    return adata_mz


def update_annotation_tomatch_adatamz(annotationsObj: MetabolomeAnnotations, adata_mz):  # zd.MetabolomeAnnotations
    # update the annotations to match with the adata_mz mz features
    me_annots = annotationsObj.model_copy()  # copy is deprecated
    okmz = adata_mz.var['mz']
    me_annots.preloaded = me_annots.preloaded.loc[
        me_annots.preloaded['mz'].isin(okmz),:]
    me_annots.unfiltered = me_annots.unfiltered.loc[
        me_annots.unfiltered['mz'].isin(okmz),:]

    return me_annots


def verification(adata_mz, me_annots: MetabolomeAnnotations):
    print("verifying mz match")
    try:
        bool_equal = np.all(np.array(
            sorted(me_annots.unfiltered["mz"].unique())) == np.array(
            sorted(adata_mz.var["mz"].tolist())))
        print(f"perfect match between mz in me_a.df "
              f"and adata.X: => {bool_equal} <=")
    except ValueError:
        print("mz do not match between adata_mz and annotation table")
        print("unique annotated mz: ", len(me_annots.unfiltered["mz"].unique(
        )))
        print("unique mz in adata_mz: ", len(adata_mz.var["mz"].unique()))


def get_rid_of_poor_var_mz(adata_mz):
    print("Excluding mz features with extremely low variance")
    means_vars, variance_computed = sparsefuncs.mean_variance_axis(
        adata_mz.X, axis=0)
    #logger.info("variance info: %s", pd.DataFrame(variance_computed).describe())
    min_accept_tile = np.percentile(variance_computed, 5)
    #low_variance = variance_computed[variance_computed >= min_accept_centile]
    adata_mz.var["computed_var"] = variance_computed
    good_mz = adata_mz.var.loc[adata_mz.var["computed_var"
                     ] > min_accept_tile, "mz"].to_numpy()
    adata_mz = adata_mz[:, [str(i) for i in good_mz]]

    return adata_mz


def closeness_by_fraction(value1: float, value2: float, fraction: float=0.8) -> bool:
    """
    Checks if two values are close to each other, by their max/min quotient

    :param value1: First value
    :param value2: Second value
    :param fraction: The x in [0.5, 1] expected for max_value/min_value ratio
    :return: True if the values are close (ratio >= fraction), False otherwise
    """
    assert 0.5 <= fraction <= 1, "'fraction' must be in the interval [0.5, 1]"
    max_value = max(value1, value2)
    min_value = min(value1, value2)

    return (max_value / min_value) >= fraction



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
    print(adata.uns['moranI'].head())

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
    print("susceptible_of_refinement: ", susceptible_of_refinement)


def ugly_long_function(mz_met_pathway_piled: pd.DataFrame):
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


def my_switch(x: int):
    match x:
        case x if x == 0:
            return "unchanged"
        case x if (x >= 1) and (x <= 5):
            return "1to5Kegg"
        case x if (x >= 6) and (x <= 10):
            return "6to10Kegg"
        case x if x > 10:
            return "11andmoreKegg"
        case _:
            return "FAILED"


def find_category_by_nb_kegg(x: int):
    out = ""
    if x == 0:
        out = "unchanged"
    elif (x >= 1) and (x <= 5):
        out = "1to5Kegg"
    elif (x >= 6) and (x <= 10):
        out = "6to10Kegg"
    elif x > 10:
        out = "11andmoreKegg"
    return out


def evaluate_spot1_change(refinement_mz_pws_dict, Xbinary_pathways,
                             X_metabolites_probs_df, spot1, spot2):

    present_pathway = 1
    absent_pathway = 0
    the_mz_for_refinement_list = list(refinement_mz_pws_dict.keys())
    modifiable_keggs = set()
    for one_mz in the_mz_for_refinement_list:
        pws_tuples_l = list(refinement_mz_pws_dict[str(one_mz)].keys())
        for pws_tuple in pws_tuples_l:
            pw1 = pws_tuple[0];   pw2 = pws_tuple[1]
            metabolites_pw1 = refinement_mz_pws_dict[str(one_mz)][pws_tuple][pw1]
            metabolites_pw2 = refinement_mz_pws_dict[str(one_mz)][pws_tuple][pw2]
            pw1andpw2_metabolites = list(set(
                metabolites_pw1).union(set(metabolites_pw2)))
            pw_score = Xbinary_pathways.loc[[spot1, spot2], [pw1, pw2]]
            metabo_probs_local = X_metabolites_probs_df.loc[
                [spot1, spot2], pw1andpw2_metabolites]
            if (metabo_probs_local.sum().sum()) == 0:  # impossible to proceed, all values are zero
                continue
            else:
                if pw_score.loc[spot1, pw1] == present_pathway:
                    if pw_score.loc[spot2, pw1] == present_pathway:
                        modifiable_keggs = modifiable_keggs.union(set(metabolites_pw1))

    return modifiable_keggs


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


def subset_m_probs_curr_pw(metabolites_probs_here, pathwaysObj, pw):  # TODO is this used ?

    theor_compounds_in_pw = pathwaysObj.to_dict()[pw]
    real_compounds_in_pw = set(theor_compounds_in_pw).intersection(
        set(metabolites_probs_here.columns.tolist()))

    return metabolites_probs_here[list(real_compounds_in_pw)].T


if __name__ == "__main__":
    import time
    t0 = time.time()
    print(my_switch(5))
    elapsed1 = time.time() - t0
    t0 = time.time()
    print(find_category_by_nb_kegg(5))
    elapsed2 = time.time() - t0
    if elapsed1 < elapsed2:
        print(elapsed1)
        print("switch is faster")
    elif elapsed2 < elapsed1:
        print("switch is slooooower")
    print(elapsed1)
    print(elapsed2)