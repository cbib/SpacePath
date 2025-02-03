import numpy as np
import anndata as ad
import scanpy as sc
import pandas as pd
from typing import Dict, List
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
import sspa as sspa
import gseapy
import zonist.helpers as zh

"""
** Note on pathway scoring problem **:

 tested  ssgsea on the data that has 100K spots : did not work fast
 also used implementation done by the sspa package:  did not work fast neither
OR
pending:
use decoupler ?, see : https://liana-py.readthedocs.io/en/latest/notebooks/misty.html  
"""

def spots_groups_by_random_choice(spots_names: np.array([str]), n_groups=4):
    seen_list = []
    n_spots_by_group = int(len(spots_names) / n_groups)
    residual_n = len(spots_names) % n_groups
    groups_dict = dict()
    for i in range(0, (n_groups - 1)):
        remaining_spots_l = list(set(spots_names) - set(seen_list))
        sps = np.random.choice(np.array(remaining_spots_l), n_spots_by_group, replace=False)
        groups_dict[i] = list(sps)
        seen_list += list(sps)

    n_spots_last_group = n_spots_by_group + residual_n   # size of last group
    remaining_spots_l = list(set(spots_names) - set(seen_list))
    sps = np.random.choice(np.array(remaining_spots_l), n_spots_last_group, replace=False)
    groups_dict[n_groups - 1] = list(sps)

    return groups_dict


def add_groups_to_spots_df(df, groups_dict, col_name="group"):
    df[col_name] = ""
    for igr in groups_dict.keys():
        df.loc[df["name"].isin(groups_dict[igr]), col_name] = str(igr)

    return df


def spots_groups_by_rectangular_cuts(spots_names: List[str], n_groups=4):
    """expects each name in the form xXX_yYY, e.g x88_y99"""
    coord_arr = zh.generate_coordinates_arr(spots_names)  # TODO: replace by zd.helpers.
    assert n_groups in [2, 4], "only 2 or 4 groups supported to date"
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

    return groups_dict, df  # TODO in the future avoid output df


def simple_test_groups_xy():
    better = []
    for k in range(4, 150):
        for h in range(4, 93):
            better.append(f"x{k}_y{h}")

    rsult, df = spots_groups_by_rectangular_cuts(better, 2)  # TODO avoid output df

    import seaborn as sns
    sns.scatterplot(data=df, x="x", y="y", hue="group", s=6)
    plt.show()

    rsu2 = spots_groups_by_random_choice(np.array(better), 4)
    df = add_groups_to_spots_df(df, rsu2, "group2")
    sns.scatterplot(data=df, x="x", y="y", hue="group2", s=6)
    plt.show()


def compute_kpca(pw_as_df, adata):
    print("Computing Pathway Scores with kpca from sspa package")
    kpca_model = sspa.sspa_kpca.sspa_KPCA(
        pw_as_df, min_entity=2, random_state=0)
    df_log_transformed = np.log(adata.to_df() + 1)  # same as sc.pp.log1p does
    kpca_scores = kpca_model.fit_transform(df_log_transformed)

    return kpca_scores


def compute_kpca_on_batches(pw_as_df, adata, groups_d: Dict[int, list]) -> Dict[int, dict]:
    print(f"Computing Pathway Scores using {len(groups_d)} batches...")
    kpca_scores_d = dict()
    for k in groups_d.keys():
        print("Batch:", str(k))
        adata_here = adata[groups_d[k], :].copy()
        kpca_scores_d[k] = compute_kpca(pw_as_df, adata_here)

    return kpca_scores_d


def add_missing_cols(df, expected_cols: set):
    absent_expected_cols = list(set(expected_cols) - set(df.columns.tolist()))
    for aec in absent_expected_cols:
        df[aec] = np.nan
    return df


def verticalassembly(kpca_scores_d, ordered_spots_names) -> pd.DataFrame:
    """"""
    columns_set = set()
    for k in kpca_scores_d.keys():
        columns_set = columns_set.union(set(
            kpca_scores_d[k].columns.tolist()
        ))
    ordered_columns_l = sorted(list(columns_set))
    tmp = add_missing_cols(kpca_scores_d[0], columns_set)
    tmp = tmp[ordered_columns_l]
    print([u for u in kpca_scores_d.keys()])
    for k in range(1, len(kpca_scores_d.keys())):
        tmp_k = add_missing_cols(kpca_scores_d[k], columns_set)
        tmp_k = tmp_k[ordered_columns_l]
        tmp  = pd.concat([tmp, tmp_k], axis=0) #, sort=True)

    tmp = tmp.fillna(tmp.min().min())
    tmp = tmp.loc[ordered_spots_names, :]
    return tmp


def compute_pathways_pseudoscores(adata, pathwaysObj, min_size=2, method="sum"):
    out_arr_l = []
    pw_ok = []
    for pw in pathwaysObj.pathways:
        kegg_metabolites_l = pathwaysObj.pathways[pw]
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



def assemble_adata_from_scores_result(adatakegg, scores_df):
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


def compute_ssgsea(adata, pathwaysObj):
    """will take maybe days if spots > 1e+6. Same problem as with kpca"""
    print("using gseapy")
    ssgsea_res = gseapy.ssgsea(data=adata.to_df().T,
                    gene_sets=pathwaysObj.pathways,
                    min_size=2,
                    outdir=None,
                    sample_norm_method='rank',
                    no_plot=True)
    ssgsea_scores = ssgsea_res.res2d.pivot(index='Name', columns='Term', values='NES')
    ssgsea_scores.index.name = "obs"
    ssgsea_scores.columns.name = "ID"

    return ssgsea_scores

