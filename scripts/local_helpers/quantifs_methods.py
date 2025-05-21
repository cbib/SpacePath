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
                   pathways_ids_df: pd.DataFrame,
                   pathways_piled: pd.DataFrame):  # deleted arg: original_pathways_dict
    score_adata = ad.AnnData(
        X=csr_matrix(kpca_scores),
        var=pd.DataFrame(
            {'pathway_id': kpca_scores.columns.tolist()},
            index=kpca_scores.columns.tolist()),
        obs=adata.obs,
        obsm=adata.obsm
    )
    score_adata.uns['spatial'] = {'my_slice_name': {
        'scalefactors': {'fiducial_diameter_fullres': 100,
                         'spot_diameter_fullres': 1,
                         'tissue_hires_scalef': 0.5,
                         'tissue_lowres_scalef': 0.5}}}
    tmp = pathways_ids_df #pathwaysObj.pathways_ids.copy()
    tmp.index = tmp['pathway_id']
    tmp = tmp.loc[score_adata.var['pathway_id'],:]
    score_adata.var = score_adata.var.assign(
        pathway_name=tmp['pathway_name'].tolist())
    # score_adata.uns['former_pathways_dict'] = original_pathways_dict    # deleted
    score_adata.uns['current_pathways_df'] = pathways_piled.copy() #pathwaysObj.piled  # has m/z, KEGG and pathway id columns

    return score_adata


def compute_mz_binary_df(adata_mz: ad.AnnData):
    X_presence_absence_mask = adata_mz.X.todense() > 0
    X_mz_binary_df = pd.DataFrame(
        np.array(X_presence_absence_mask, dtype=int),
        columns= adata_mz.var['mz'].to_numpy(),
        index=adata_mz.obs.index.tolist())

    return X_mz_binary_df


def compute_metabolites_probs_df(adata_mz: ad.AnnData, annots):
    X_mz_binary_df = compute_mz_binary_df(adata_mz)
    preX_probs = []
    col_names_kegg = []
    for a_mz in X_mz_binary_df.columns:
        i_arr = X_mz_binary_df[a_mz].to_numpy()
        the_ids_here = annots.df.loc[
            annots.df['mz'] == a_mz, annots.usable_id].tolist()
        n = len(the_ids_here)
        for kid in the_ids_here:
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


# def build_coords_df(adata):  # todo: if never used, del
#     df = pd.DataFrame(data=adata.obsm['spatial'],
#                  columns=['x', 'y'],
#                  index=adata.obs.index.tolist())
#     return df


#if __name__ == "__main__":

