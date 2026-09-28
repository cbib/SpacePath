import os
import time
from typing import Dict, List, Tuple, Union

import anndata as ad
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import SpacePath.core.data as spd
import SpacePath.core.method as spm
import squidpy as sq
from scipy.sparse import csr_matrix



def join_scores_adatas(
        pathways_types: List[str], SAMPLE: str,
        subfolder: str, data_centralized_folder: str) -> ad.AnnData:
    """joins computed scores (anndata objects) in one single anndata object"""
    if not "FAMILIES" in pathways_types:
        pathways_types.append("FAMILIES")

    if not os.path.exists(os.path.join(data_centralized_folder, subfolder,
                                       f'{SAMPLE}-scores-FAMILIES.h5ad')):
        pathways_types = list(set(pathways_types) - {'FAMILIES'})

    list_of_scores_objects = list()
    list_of_dfs = list()
    for u in list(set(pathways_types)):
        adata_scores = ad.read_h5ad(
            filename=os.path.join(data_centralized_folder, subfolder,
                                  f'{SAMPLE}-scores-{u}.h5ad'))
        adata_scores.var['set_type'] = u
        list_of_scores_objects.append(adata_scores)
        df_uns = adata_scores.uns['current_pathways_df'].copy()
        df_uns.columns = ["pathway_id", "moleculeIds", "mz"]
        list_of_dfs.append(df_uns)

    joined = ad.concat(adatas=list_of_scores_objects, join='outer', axis=1)

    joined.uns['current_pathways_df'] = pd.concat(
        list_of_dfs, join='outer', axis=0, ignore_index=True)

    joined.obsm = adata_scores.obsm
    joined.uns['spatial'] = {'my_slice_name': {
        'scalefactors': {'fiducial_diameter_fullres': 100,
                         'spot_diameter_fullres': 1,
                         'tissue_hires_scalef': 0.5,
                         'tissue_lowres_scalef': 0.5}}}

    return joined


def score_as_adata(scores_df: pd.DataFrame,
                   sets_ids_df: pd.DataFrame,
                   pathways_piled: pd.DataFrame,
                   coords: np.ndarray) -> ad.AnnData:
    """
    wraps the computed scores into an anndata object
    also:
     - stores the pathways' information in .uns['current_pathways_df'],
        which is the table with column names:
           - pathway_id
           - moleculeIds, or KEGG
           - mz
     - stores coordinates dict {'spatial': <np.array>} into .obsm
     - adds dictionary of parameters for plotting into .uns['spatial']
    """
    score_adata = ad.AnnData(
        X=csr_matrix(scores_df),
        var=pd.DataFrame(
            {'pathway_id': scores_df.columns.tolist()},
            index=scores_df.columns.tolist()),
        obs=pd.DataFrame({'obs': scores_df.index.tolist()},
                         index=scores_df.index.tolist()),
        obsm={'spatial': coords}
    )
    score_adata.uns['spatial'] = {'my_slice_name': {
        'scalefactors': {'fiducial_diameter_fullres': 100,
                         'spot_diameter_fullres': 1,
                         'tissue_hires_scalef': 0.5,
                         'tissue_lowres_scalef': 0.5}}}
    tmp = sets_ids_df
    tmp.index = tmp['pathway_id']
    tmp = tmp.loc[score_adata.var['pathway_id'], :]
    score_adata.var = score_adata.var.assign(
        pathway_name=tmp['pathway_name'].tolist())

    score_adata.uns[
        'current_pathways_df'] = pathways_piled.copy()  # has m/z, KEGG/moleculeIds and pathwayid columns

    return score_adata


def run_pathway_scores_steps_produce_anndata(metabolites_computed_df,
                                             coords: np.ndarray,
                                             pws_obj: spd.MetabolicPathways,
                                             min_set_size: int,
                                             type_abbr: str) -> Tuple[
    ad.AnnData, float]:
    """
    Uses core library to compute scores dataframe,
       and wraps the scores into a new anndata object

    returns: - anndata object
             - computation time in seconds
    """
    updated_config = spm.ComputePathwaysScoresConfig(type_abbr=type_abbr,
                                                     min_entity=min_set_size)
    t0 = time.time()
    scores_df = spm.ComputePathwaysScores.run(
        pws_obj.to_dict(), metabolites_computed_df, config=updated_config
    )
    time_elapsed = time.time() - t0
    scores_adata = score_as_adata(scores_df,
                                  pws_obj.sets_ids_df,
                                  pws_obj.piled,
                                  coords)

    return scores_adata, time_elapsed


## differential routine

def get_top_n(group: pd.DataFrame, *args) -> pd.DataFrame:
    n_top_per_cluster = args[0]
    direction_type = args[1]
    sorted_group = pd.DataFrame()
    if direction_type == "up":
        group = group.loc[group['logfoldchanges'] > 0, :].copy()
        sorted_group = group.sort_values(
            by=['pvals_adj', 'logfoldchanges'], ascending=[True, False])
    elif direction_type == "down":
        group = group.loc[group['logfoldchanges'] < 0, :].copy()
        sorted_group = group.sort_values(
            by=['pvals_adj', 'logfoldchanges'], ascending=[True, True])

    return sorted_group.head(n_top_per_cluster)


def differential_results_filter(differential_results: pd.DataFrame,
                                abs_lfc: float = 0.8,
                                pvals_adj: float = 0.05,
                                n_top_per_cluster: Union[str, int] = "auto",
                                total_n_features_out: int = 250,
                                up_down_proportion: tuple = (0.7, 0.3)
                                ) -> pd.DataFrame:
    if n_top_per_cluster == "auto":
        n_top_per_cluster = int(total_n_features_out / len(
            differential_results['group'].unique()))
    up_n = int(n_top_per_cluster * up_down_proportion[0])
    down_n = int(n_top_per_cluster * up_down_proportion[1])
    tmp = differential_results.loc[(abs(differential_results[
                                            'logfoldchanges']) >= abs_lfc) & (
                                           differential_results[
                                               'pvals_adj'] <= pvals_adj), :].copy()

    tmp = tmp.assign(abslfc=abs(tmp['logfoldchanges']),
                     cluster=tmp['group'].astype(str))
    tmp = tmp.sort_values(by=['group', 'pvals_adj'], ascending=[True, False])
    tmp_up = tmp.groupby('cluster', group_keys=False, as_index=False
                         ).apply(get_top_n, up_n, "up",
                                 include_groups=False
                                 # avoid DeprecationWarning
                                 ).reset_index(drop=True).copy()
    tmp_down = tmp.groupby('cluster', group_keys=False, as_index=False
                           ).apply(get_top_n, down_n, "down",
                                   include_groups=False
                                   # avoid DeprecationWarning
                                   ).reset_index(drop=True).copy()
    tmp = pd.concat([tmp_up, tmp_down], axis=0)

    return tmp


def squidpy_differential_analysis(adata: ad.AnnData,
                                  group_str: str) -> pd.DataFrame:
    """
    computes differential analysis using the clusters registered in adata.obs

    group_str (str) is the name of the column in anndata.obs (e.g. 'kmeans-5')
    """
    sc.tl.rank_genes_groups(adata, groupby=group_str,
                            groups='all',
                            reference='rest',
                            method='wilcoxon',
                            key_added="feats_by_cluster",
                            corr_method='benjamini-hochberg')

    differential_results_by_cluster = sc.get.rank_genes_groups_df(
        adata, key="feats_by_cluster",
        group=None)  # returns all groups

    return differential_results_by_cluster


def safe_insert_clusters_column(adata: ad.AnnData, clusters_df: pd.DataFrame,
                                column_interest: str='kmeans-5', spots='obs'):
    """
    inserts clusters column into anndata.obs from an external dataframe
    :param adata: anndata.AnnData
    :param clusters_df: dataframe with two columns, e.g. ['obs', 'kmeans-5']
    :param column_interest: string column name to insert into anndata.obs
    :param spots: string column name that has the spots names as strings
    """
    adata.obs = adata.obs.assign(**{spots: adata.obs.index.tolist()})
    clusters_df = clusters_df.drop_duplicates()
    clusters_df.index = clusters_df[spots].tolist()
    # sort entire dataframe to match with adata.obs
    clusters_df = clusters_df.T[adata.obs[spots].tolist()].T.copy()

    ordered_categories_clusters: List[int] = np.sort(
        clusters_df[column_interest].unique()
    ).tolist()

    adata.obs = adata.obs.assign(**{column_interest : pd.Categorical(
        clusters_df[column_interest].tolist(),
        categories=ordered_categories_clusters)})

    return adata


def plot_clustering_results(adata_mz, clustering_parameters, SAMPLE,
                            plot_folder, colors_spots: List[str]) -> None:
    """Saves figure of clustering results"""
    str_out = f"{clustering_parameters['choice']}-{clustering_parameters['numeric_parameter']}"

    plot_out_path = os.path.join(plot_folder,
                                 f'{SAMPLE}-clusters-{str_out}.pdf')

    colors_spots = colors_spots[:adata_mz.obs[str_out].nunique()]

    fig, ax = plt.subplots(1, 1)
    sq.pl.spatial_scatter(
        adata_mz, color=[str_out], img=False,
        title="Clustering of spots",
        linewidth=0, edgecolors=(0, 0, 0, 0),
        shape="circle", size=1.15,  # 'circle', 'square', 'hex'
        palette=ListedColormap(colors_spots),
        ax=ax
    )
    ax.collections[0].set_rasterized(True)

    plt.savefig(os.path.abspath(plot_out_path),
                format="pdf",
                dpi=600,
                bbox_inches='tight')
    plt.close()

    print("saved figure")


def compute_leiden__weighted_by_spatial_graph(
        datacopy, resolution_param, alpha_gr=0.25
) -> Tuple[ad.AnnData, List[str]]:
    """
    Computes Leiden community detection algorithm for the spots clustering

    :param datacopy: a copy of the anndata object with the norm intensities
    :param resolution_param: resolution parameter for the Leiden algorithm
    :param alpha_gr: the weight of the spatial graph

    follows 'spatial domains' instructions as in www.sc-best-practices.org
    """
    try:
        sc.pp.neighbors(datacopy, n_pcs=20)
    except:
        sc.pp.neighbors(datacopy)
    nn_graph_genes = datacopy.obsp["connectivities"]
    sq.gr.spatial_neighbors(datacopy)
    nn_graph_space = datacopy.obsp["spatial_connectivities"]

    alpha = alpha_gr
    joint_graph = (1 - alpha) * nn_graph_genes + alpha * nn_graph_space
    sc.tl.leiden(datacopy,
                 adjacency=joint_graph,
                 n_iterations=2,
                 flavor="igraph",
                 directed=False, key_added="Leiden_sg",
                 resolution=resolution_param)
    the_labels_str_l = datacopy.obs['Leiden_sg'].astype(str).tolist()

    return datacopy, the_labels_str_l


## more custom routines

def safe_load_processed_data(SAMPLE: str, subfolder: str,
                             data_centralized_folder: str) -> ad.AnnData:
    """safely load mzs intensities (anndata) with only good quality ions"""
    adata_mz = ad.read_h5ad(filename=os.path.join(
        data_centralized_folder, subfolder, f'{SAMPLE}.h5ad'))
    adata_mz = adata_mz[:,
    adata_mz.var.loc[adata_mz.var['final_quality'] == "good",
    "mz"].astype(str).tolist()]

    return adata_mz


def read_clean_annotations(SAMPLE: str, pathways_types: List[str],
                           subfolder: str,
                           data_centralized_folder: str) -> Dict[
    str, pd.DataFrame]:
    """provide the dictionary of annotations dataframes,
     where the keys are the pathways_types
     Note:  ignores FAMILIES for annotations loading, as the annotations
             in such case are the same that served to 'WIKI' type of pathway
     """
    pathways_types = list(set(pathways_types) - {"FAMILIES"})

    annots_dict = {}
    for pw_type in pathways_types:
        tmpfi = os.path.join(data_centralized_folder,
                             subfolder,
                             f'{SAMPLE}-{pw_type}-CleanAnnotations.tsv')
        annots_dict[pw_type] = pd.read_csv(tmpfi,
                                           sep='\t', header=0, index_col=None)

    return annots_dict


def safe_load_data_and_clean_annotations(SAMPLE: str,
                                         pathways_types: List[str],
                                         subfolder: str,
                                         data_centralized_folder: str
) -> Tuple[ad.AnnData, Dict[str,pd.DataFrame]]:
    """
    safely load mzs intensities (anndata) and annotations (dataframe):
     only good quality ions AND only clean annotations, for both objects
    """
    pathways_types = list(set(pathways_types) - {"FAMILIES"})

    adata_mz = safe_load_processed_data(SAMPLE, subfolder,
                                        data_centralized_folder)
    annots_dict = read_clean_annotations(SAMPLE, pathways_types, subfolder,
                                         data_centralized_folder)

    adequate_mz = set()
    for pw_type in pathways_types:
        adequate_mz = adequate_mz.union(
            set(annots_dict[pw_type]['mz'].astype(str).unique())
        )
    adata_mz = adata_mz[:, list(adequate_mz)].copy()

    return adata_mz, annots_dict


def prepare_molecular_ids_df(top_differential: pd.DataFrame,
                             clean_ann_df: pd.DataFrame) -> pd.DataFrame:
    """yields clean df with columns 'moleculeIds', 'mz', 'moleculeNames'"""
    top_differential_list = top_differential['names'
    ].astype(float).unique().tolist()  # names are mzs in the diff result
    molecules_ids_df = clean_ann_df.loc[clean_ann_df['mz'
    ].isin(top_differential_list),
    ['moleculeIds', 'mz', 'moleculeNames']].copy()

    molecules_ids_df.dropna(how='any', inplace=True)
    molecules_ids_df.drop_duplicates(inplace=True)
    molecules_ids_df.reset_index(inplace=True)

    return molecules_ids_df
