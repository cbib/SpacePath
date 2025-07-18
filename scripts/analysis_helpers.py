import os
import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc

from typing import Dict, List, Union

from scipy.sparse import csr_matrix
import smpath.method as spm


## score outer level
# todo: "send to smpath"
def score_as_adata(scores_df,
                   sets_ids_df_df: pd.DataFrame,
                   pathways_piled: pd.DataFrame, coords: np.array) -> ad.AnnData:
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
    tmp = sets_ids_df_df
    tmp.index = tmp['pathway_id']
    tmp = tmp.loc[score_adata.var['pathway_id'],:]
    score_adata.var = score_adata.var.assign(
        pathway_name=tmp['pathway_name'].tolist())

    score_adata.uns['current_pathways_df'] = pathways_piled.copy() # # has m/z, KEGG/moleculeIds and pathway id columns

    return score_adata




def run_pathway_scores_steps_produce_anndata(metabolites_computed_df, coords,
                                             pws_obj,
                                             min_set_size,
                                             type_abbr) -> ad.AnnData:

    updated_config = spm.ComputePathwaysScoresConfig(type_abbr=type_abbr, min_set_size=min_set_size)
    scores_df = spm.ComputePathwaysScores.run(
        pws_obj.to_dict(), metabolites_computed_df, config=updated_config
    )
    scores_adata = score_as_adata(scores_df,
                   pws_obj.sets_ids_df,
                   pws_obj.piled,
                   coords)

    return scores_adata

## differential stuff


def get_top_n(group, *args):
    ntop_per_cluster = args[0]
    sorted_group = group.sort_values(
        by=['pvals_adj', 'abslfc'], ascending=[True, False])
    return sorted_group.head(ntop_per_cluster)


def differential_results_filter(differential_results, lfc: float=0.8,
                                pvals_adj: float=0.05,
                                ntop_per_cluster: Union[str, int]="auto",
                                total_n_features:int=250)-> pd.DataFrame:
    if ntop_per_cluster == "auto":
        ntop_per_cluster = int(total_n_features / len(
            differential_results['group'].unique()))

    tmp = differential_results.loc[(abs(differential_results[
             'logfoldchanges']) >= lfc) & (
            differential_results['pvals_adj'] <= pvals_adj), :].copy()

    tmp = tmp.assign(abslfc=abs(tmp['logfoldchanges']),
                     cluster=tmp['group'].astype(str))
    tmp = tmp.sort_values(by=['group', 'pvals_adj'], ascending=[True, False])
    tmp = tmp.groupby('cluster', group_keys=False, as_index=False
                      ).apply(get_top_n, ntop_per_cluster,
                              include_groups=False  # avoid DeprecationWarning
                              ).reset_index(drop=True)
    return tmp


def extract_top_differential_features(adata, group_str, ntop_per_cluster
                                      )-> pd.DataFrame:
    sc.tl.rank_genes_groups(adata, groupby=group_str,
                            groups='all',
                            reference='rest',
                            method='wilcoxon',
                            key_added="top_feats_by_cluster",
                            corr_method='benjamini-hochberg')

    differential_results_by_cluster = sc.get.rank_genes_groups_df(
        adata, key="top_feats_by_cluster",
        group=None)  # returns all groups

    top_differential = differential_results_filter(
        differential_results_by_cluster, ntop_per_cluster=ntop_per_cluster)

    return top_differential



## more routines


def safe_load_processed_data(SAMPLE, subfolder, data_centralized_folder):

    adata_mz = ad.read_h5ad(filename=os.path.join(
        data_centralized_folder, subfolder, f'{SAMPLE}.h5ad'))
    adata_mz = adata_mz[:,
               adata_mz.var.loc[adata_mz.var['final_quality'] == "good",
               "mz"].astype(str).tolist()]


    return adata_mz


def read_clean_annotations(SAMPLE, pathways_types, subfolder, data_centralized_folder):
    annots_dict = {}
    for pw_type in pathways_types:
        tmpfi = os.path.join(data_centralized_folder,
                             subfolder,
                             f'{SAMPLE}-{pw_type}-CleanAnnotations.tsv')
        annots_dict[pw_type] = pd.read_csv(tmpfi,
                                          sep='\t', header=0, index_col=None)

    return annots_dict


def safe_load_data_and_clean_annotations(SAMPLE,
                                         pathways_types,
                                         subfolder, data_centralized_folder):

    adata_mz = safe_load_processed_data(SAMPLE, subfolder, data_centralized_folder)
    annots_dict = read_clean_annotations(SAMPLE, pathways_types, subfolder, data_centralized_folder)

    adequate_mz = set()
    for pw_type in pathways_types:
        adequate_mz = adequate_mz.union(
            set(annots_dict[pw_type]['mz'].astype(str).unique())
        )
    adata_mz = adata_mz[:, list(adequate_mz)].copy()

    return adata_mz, annots_dict




def prepare_molecular_ids_df(top_differential, clean_ann_df):
    top_differential_list = top_differential['names'
        ].astype(float).unique().tolist()  # names are mzs in the diff result
    molecules_ids_df = clean_ann_df.loc[clean_ann_df['mz'
        ].isin(top_differential_list),
        ['moleculeIds', 'mz', 'moleculeNames'] ].copy()

    molecules_ids_df.dropna(how='any', inplace=True)
    molecules_ids_df.drop_duplicates(inplace=True)
    molecules_ids_df.reset_index(inplace=True)

    return molecules_ids_df


def files_of_pathways__location_info():
    # stored in the same centralized data folder (same level of the tissues' folders)
    pathways_files_defaults = {
        "KEGG" : {'file': "KEGG_hsa_pathways_compounds_R110.gmt",
                      'subfolder': "metabolites-mappings",
                      'id_type': 'KEGG'},
        "WIKI" : {'file': "WIKIPATHWAYS_human_multi_id_redundant.mmt",
                          'subfolder': "metabolites-mappings",
                          'id_type': 'moleculeIds'}}
    return pathways_files_defaults


