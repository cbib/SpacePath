import os
import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
import time

from typing import Dict, List, Tuple, Union

from scipy.sparse import csr_matrix
import smpath.method as spm



def join_scores_adatas(
        pathways_types: List[str], SAMPLE: str,
        subfolder:str, data_centralized_folder: str) -> ad.AnnData:
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
                   pathways_piled: pd.DataFrame, coords: np.array) -> ad.AnnData:
    """
    wraps the computed scores into an anndata object
    also:
     - stores the pathways' information in .uns['current_pathways_df'],
        which is the table with column names:
           - pathway_id
           - moleculeIds, or KEGG
           - mz
     - stores coordinates array into .obsm
     - adds dictionary of parameters for plotting into .uns['spatial']
    """
    score_adata = ad.AnnData(
        X=csr_matrix(scores_df),
        var=pd.DataFrame(
            {'pathway_id': scores_df.columns.tolist()},
            index=scores_df.columns.tolist()),
        obs=pd.DataFrame({'obs': scores_df.index.tolist()},
                         index=scores_df.index.tolist()),
        obsm=coords
    )
    score_adata.uns['spatial'] = {'my_slice_name': {
        'scalefactors': {'fiducial_diameter_fullres': 100,
                         'spot_diameter_fullres': 1,
                         'tissue_hires_scalef': 0.5,
                         'tissue_lowres_scalef': 0.5}}}
    tmp = sets_ids_df
    tmp.index = tmp['pathway_id']
    tmp = tmp.loc[score_adata.var['pathway_id'],:]
    score_adata.var = score_adata.var.assign(
        pathway_name=tmp['pathway_name'].tolist())

    score_adata.uns['current_pathways_df'] = pathways_piled.copy() # # has m/z, KEGG/moleculeIds and pathway id columns

    return score_adata


def run_pathway_scores_steps_produce_anndata(metabolites_computed_df, coords,
                                             pws_obj,
                                             min_set_size,
                                             type_abbr) -> Tuple[ad.AnnData, float]:
    """
    Uses smpath library to compute scores dataframe,
       and wraps the scores into a new anndata object

    returns: - anndata object
             - computation time in seconds
    """
    updated_config = spm.ComputePathwaysScoresConfig(type_abbr=type_abbr,
                                                     min_set_size=min_set_size)
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

def get_top_n_simple(group: pd.DataFrame, *args) -> pd.DataFrame:
    n_top_per_cluster = args[0]
    sorted_group = group.sort_values(
        by=['pvals_adj', 'abslfc'], ascending=[True, False])
    return sorted_group.head(n_top_per_cluster)


def get_top_n(group: pd.DataFrame, *args) -> pd.DataFrame:
    n_top_per_cluster = args[0]
    direction_type = args[1]
    sorted_group= pd.DataFrame()
    if direction_type == "up":
        group = group.loc[group['logfoldchanges'] > 0,:].copy()
        sorted_group = group.sort_values(
            by=['pvals_adj', 'logfoldchanges'], ascending=[True, False])
    elif direction_type == "down":
        group = group.loc[group['logfoldchanges'] < 0, :].copy()
        sorted_group = group.sort_values(
            by=['pvals_adj', 'logfoldchanges'], ascending=[True, True])

    return sorted_group.head(n_top_per_cluster)


def differential_results_filter(differential_results,
                                abs_lfc: float=0.8,
                                pvals_adj: float=0.05,
                                n_top_per_cluster: Union[str, int]="auto",
                                total_n_features:int=250,
                                up_down_proportion:tuple=(0.7, 0.3)
                                )-> pd.DataFrame:
    if n_top_per_cluster == "auto":
        n_top_per_cluster = int(total_n_features / len(
            differential_results['group'].unique()))
    up_n = int(n_top_per_cluster * up_down_proportion[0])
    down_n = int(n_top_per_cluster * up_down_proportion[1])
    tmp = differential_results.loc[(abs(differential_results[
             'logfoldchanges']) >= abs_lfc) & (
            differential_results['pvals_adj'] <= pvals_adj), :].copy()

    tmp = tmp.assign(abslfc=abs(tmp['logfoldchanges']),
                     cluster=tmp['group'].astype(str))
    tmp = tmp.sort_values(by=['group', 'pvals_adj'], ascending=[True, False])
    tmp_up = tmp.groupby('cluster', group_keys=False, as_index=False
                      ).apply(get_top_n, up_n, "up",
                              include_groups=False  # avoid DeprecationWarning
                              ).reset_index(drop=True).copy()
    tmp_down = tmp.groupby('cluster', group_keys=False, as_index=False
                      ).apply(get_top_n, down_n, "down",
                              include_groups=False  # avoid DeprecationWarning
                              ).reset_index(drop=True).copy()
    tmp = pd.concat([tmp_up, tmp_down], axis=0)

    return tmp


def squidpy_differential_analysis(adata, group_str)-> pd.DataFrame:
    sc.tl.rank_genes_groups(adata, groupby=group_str,
                            groups='all',
                            reference='rest',
                            method='wilcoxon',
                            key_added="top_feats_by_cluster",
                            corr_method='benjamini-hochberg')

    differential_results_by_cluster = sc.get.rank_genes_groups_df(
        adata, key="top_feats_by_cluster",
        group=None)  # returns all groups

    return differential_results_by_cluster


## more custom routines

def safe_load_processed_data(SAMPLE, subfolder, data_centralized_folder)->ad.AnnData:
    """safely load mzs intensities (anndata) with only good quality ions"""
    adata_mz = ad.read_h5ad(filename=os.path.join(
        data_centralized_folder, subfolder, f'{SAMPLE}.h5ad'))
    adata_mz = adata_mz[:,
               adata_mz.var.loc[adata_mz.var['final_quality'] == "good",
               "mz"].astype(str).tolist()]

    return adata_mz


def read_clean_annotations(SAMPLE, pathways_types, subfolder,
                           data_centralized_folder)-> Dict[str, pd.DataFrame]:
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
                                         subfolder:str, data_centralized_folder:str) ->(ad.AnnData, pd.DataFrame):
    """
    safely load mzs intensities (anndata) and annotations (dataframe):
     only good quality ions AND only clean annotations, for both objects
    """
    adata_mz = safe_load_processed_data(SAMPLE, subfolder, data_centralized_folder)
    annots_dict = read_clean_annotations(SAMPLE, pathways_types, subfolder, data_centralized_folder)

    adequate_mz = set()
    for pw_type in pathways_types:
        adequate_mz = adequate_mz.union(
            set(annots_dict[pw_type]['mz'].astype(str).unique())
        )
    adata_mz = adata_mz[:, list(adequate_mz)].copy()

    return adata_mz, annots_dict


def prepare_molecular_ids_df(top_differential:pd.DataFrame,
                             clean_ann_df:pd.DataFrame) -> pd.DataFrame:
    """yields clean df with columns 'moleculeIds', 'mz', 'moleculeNames'"""
    top_differential_list = top_differential['names'
        ].astype(float).unique().tolist()  # names are mzs in the diff result
    molecules_ids_df = clean_ann_df.loc[clean_ann_df['mz'
        ].isin(top_differential_list),
        ['moleculeIds', 'mz', 'moleculeNames'] ].copy()

    molecules_ids_df.dropna(how='any', inplace=True)
    molecules_ids_df.drop_duplicates(inplace=True)
    molecules_ids_df.reset_index(inplace=True)

    return molecules_ids_df




