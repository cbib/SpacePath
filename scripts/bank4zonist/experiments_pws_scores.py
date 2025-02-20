
import numpy as np
import anndata as ad
import pandas as pd
from typing import List, Union, Dict
import zonist.helpers as zh


def compounds_counts_by_pathway(pathwaysObj):
    """
    Generates a table with the number of compounds by pathway
    :return: e.g.

    """
    pw_id_l = []
    pw_name = []
    pw_compounds = []
    pw_count = []
    for pw_id in pathwaysObj.to_dict().keys():
        metabolites = pathwaysObj.to_dict()[pw_id]
        pw_id_l.append(pw_id)
        pw_name.append(pathwaysObj.pathways_ids.loc[
                           pathwaysObj.pathways_ids["pathway_id"] == pw_id,
                           "pathway_name"].item())
        pw_compounds.append("; ".join(metabolites))
        pw_count.append(len(metabolites))
    uf = pd.DataFrame({'pathway_id': pw_id_l, 'name': pw_name,
                       'cmpslist': pw_compounds, 'count': pw_count})
    uf = uf.sort_values(by='count', ascending=False)
    return uf


def assemble_adata_from_scores_result(adatakegg, scores_df):  # TODO: see if it is practical or uselesss (not yet used)
    adata_scores = ad.AnnData(
            X=csr_matrix(scores_df),
            var=pd.DataFrame(
                {'pathway_id': scores_df.columns.tolist()},
                index=scores_df.columns.tolist()),
            obs=adatakegg.obs,
            obsm=adatakegg.obsm,
            uns=adatakegg.uns
        )
    return adata_scores


def time_formatting(timeseconds):
    hours = int(timeseconds / 3600)
    remaining_seconds = timeseconds % 3600
    minutes = int(remaining_seconds / 60)
    seconds = np.around(remaining_seconds % 60, 1)
    return hours, minutes, seconds


def add_groups_to_spots_df(df, groups_dict, col_name="group"):
    df[col_name] = ""
    for igr in groups_dict.keys():
        df.loc[df["name"].isin(groups_dict[igr]), col_name] = str(igr)

    return df


def pathways_scores_by_rectangular_cuts(pathwaysObj, adata, min_entity=2,
                                        n_groups=4):
    """only experimental purpose, it is not better than random sorting"""
    pw_as_df = reshape_pathways_for_kpca(pathwaysObj)
    groups_d, df = spots_groups_by_rectangular_cuts(
        adata.obs['obs'].tolist(), n_groups=n_groups)
    kpca_scores_d = compute_kpca_on_spot_groups(
        pw_as_df, adata, groups_d, min_entity=min_entity)
    kpca_scores = verticalassembly(
        kpca_scores_d, ordered_spots_names=adata.obs['obs'].tolist())

    return kpca_scores


def spots_groups_by_rectangular_cuts(spots_names: List[str], n_groups=4):
    """expects each name in the form xXX_yYY, e.g x88_y99"""
    coord_arr = zh.generate_coordinates_arr(spots_names)
    assert n_groups in [2, 4], "only 2 or 4 groups are supported"
    df = pd.DataFrame(data=coord_arr, columns=["x", "y"])
    df["name"] = spots_names
    groups_dict = dict()
    span_x = max(coord_arr[:, 0]) - min(coord_arr[:, 0])
    span_y = max(coord_arr[:, 1]) - min(coord_arr[:, 1])
    x_cut_point = int(span_x / 2) + min(coord_arr[:, 0])
    y_cut_point = int(span_y / 2) + min(coord_arr[:, 1])
    if n_groups == 2:
        if span_x > span_y:
            g1 = df.loc[df["x"] <= x_cut_point, "name"].tolist()
            g2 = df.loc[df["x"] > x_cut_point, "name"].tolist()
        elif span_x <= span_y:
            g1 = df.loc[df["y"] <= y_cut_point , "name"].tolist()
            g2 = df.loc[df["y"] > y_cut_point, "name"].tolist()
        groups_dict[0] = g1
        groups_dict[1] = g2
    elif n_groups == 4:
        groups_dict[0] = df.loc[(df["x"] <= x_cut_point) & (
                df["y"] <= y_cut_point), "name"].tolist()
        groups_dict[1] = df.loc[(df["x"] <= x_cut_point) & (
                df["y"] > y_cut_point), "name"].tolist()
        groups_dict[2] = df.loc[(df["x"] > x_cut_point) & (
                df["y"] <= y_cut_point), "name"].tolist()
        groups_dict[3] = df.loc[(df["x"] > x_cut_point) & (
                df["y"] > y_cut_point), "name"].tolist()

    df = add_groups_to_spots_df(df, groups_dict)

    return groups_dict, df


def compute_pathways_pseudoscores(adata, pathwaysObj, min_size=2, method="sum"):
    out_arr_l = []
    pw_ok = []
    for pw in pathwaysObj.to_dict().keys():
        kegg_metabolites_l = pathwaysObj.to_dict()[pw]
        sub_df = adata[:, kegg_metabolites_l].X.todense()
        if sub_df.shape[1] >= min_size:
            if method == "sum":
                res_arr = np.array(sub_df).sum(axis=1)
            elif method == "average":
                res_arr = np.array(sub_df).mean(axis=1)
            else:
                raise NotImplementedError

            out_arr_l.append(res_arr)
            pw_ok.append(pw)
        # end if
    out_df = pd.DataFrame(data=np.array(out_arr_l).T,
                          index=adata.obs.index.tolist(), columns=pw_ok)

    return out_df


def compute_ssgsea(adata, pathwaysObj):
    """will take maybe days if spots > 1e+6. Same problem as with kpca"""
    print("using gseapy")
    ssgsea_res = gseapy.ssgsea(data=adata.to_df().T,
                    gene_sets=pathwaysObj.to_dict(),
                    min_size=2,
                    outdir=None,
                    sample_norm_method='rank',
                    no_plot=True)
    ssgsea_scores = ssgsea_res.res2d.pivot(index='Name', columns='Term', values='NES')
    ssgsea_scores.index.name = "obs"
    ssgsea_scores.columns.name = "ID"

    return ssgsea_scores


def modules_scores_by_scanpy(adata, pathways: dict, do_log=False):
    adatacopy = adata.copy()
    if do_log:
        if adatacopy.X.min().min() < 0:
            print("Warning: negative numbers before log")
            print("applying shift of the log values (all to be positives)")
            adatacopy.X = adatacopy.X + abs(adatacopy.X.min().min())
        sc.pp.log1p(adatacopy)

    pathway_names = list()
    scores_lists = list()
    for pw in pathways.keys():
        pathway_names.append(pw)
        mygenelist = pathways[pw]
        miniadata = sc.tl.score_genes(adatacopy, gene_list=mygenelist,
                                      ctrl_size=len(mygenelist), copy=True)
        scores_lists.append(miniadata.obs['score'].tolist())

    df = pd.DataFrame(np.array(scores_lists),
                      index=pathway_names,
                      columns=adatacopy.obs['obs'].tolist(),
                      )
    return df.T


def build_kegg_resolved_adataV0(adata_mz, annots_df):
    """
    Given the adata_mz and the annots_df, returns a new adata object
    OLD version
    """
    mz_list = annots_df["mz"].unique()
    intensities_arrs_list = []
    kegg_list = list()
    for one_mz in mz_list:
        local_kegg = annots_df.loc[
            annots_df["mz"] == one_mz, "KEGG"].tolist()
        local_1D_adata_mz = adata_mz[ :, str(one_mz)].copy()
        for a_kegg in local_kegg:
            intensities_arrs_list.append(
                np.array(local_1D_adata_mz.X.todense().T).flatten() / len(local_kegg))
            kegg_list.append(a_kegg)
    adata = ad.AnnData(X=csr_matrix(np.array(intensities_arrs_list).T),
                       obs=adata_mz.obs,
                       var=pd.DataFrame({"ID": kegg_list},
                                        index=kegg_list),
                       obsm=adata_mz.obsm)
    adata.uns['spatial'] = {'my_slice_name': {
        'scalefactors': {'fiducial_diameter_fullres': 100,
                         'spot_diameter_fullres': 1,
                         'tissue_hires_scalef': 0.5,
                         'tissue_lowres_scalef': 0.5}}}


    adata = attach_annots_to_adata(adata, annots_df) ## new !

    return adata