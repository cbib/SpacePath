# script for exclusive KEGG+WIKI pathways alluvial plot
#  pathways that are important to show and have significance (padj<=0.05)
#  but they are not top-listed because the families have much better moran I
#  e.g. the healthy brain.

import os
from scipy.sparse import csr_matrix
import numpy as np
import pandas as pd

from SpacePath.src.helpers import (join_scores_adatas,
                                   safe_insert_clusters_column,
                                   differential_results_filter,
                                   squidpy_differential_analysis)
import SpacePath.core.processing.measures_methods as measures_methods

from viz_results_helpers import run_alluvial_plot


def pathways_extractor(scores_merged, top_differential, relaxed_moran_for_pathways=0.02):
    # usually moran threshold of 0.2,  but for some pathways is too stringent
    top_sets_by_moran_i = scores_merged.uns['moranI'].loc[
        scores_merged.uns['moranI'][
            'I'] >= relaxed_moran_for_pathways, :].copy().index.tolist()
    top_differential = top_differential.loc[
        top_differential['names'].isin(top_sets_by_moran_i), :]
    top_differential = top_differential.loc[~top_differential['names'].str.startswith('MolFam'),:]

    return top_differential['names'].tolist()


def families_extractor(scores_merged, top_differential,
                       moran_cutoff=0.08    #  slightly lower than 0.2 for this dataset
):

    top_sets_by_moran_i = scores_merged.uns['moranI'].loc[
        scores_merged.uns['moranI'][
            'I'] >= moran_cutoff, :].copy().index.tolist()
    top_differential = top_differential.loc[
        top_differential['names'].isin(top_sets_by_moran_i), :]
    top_differential = top_differential.loc[
        top_differential['names'].str.startswith('MolFam'), :]

    return top_differential['names'].tolist()


def modified_top_selection(pathways_types, clustering_parameters, SAMPLE, subfolder,
                                       data_centralized_folder,
                           cluster_not_tissue=None,
                           n_top_per_cluster=200):
    str_out = f"{clustering_parameters['choice']}-{clustering_parameters['numeric_parameter']}"

    scores_merged = join_scores_adatas(pathways_types, SAMPLE, subfolder,
                                       data_centralized_folder)

    clusters_df = pd.read_csv(os.path.join(data_centralized_folder, subfolder,
                                           f'{SAMPLE}-clusters-{str_out}.txt'),
                              header=0)

    scores_merged = safe_insert_clusters_column(scores_merged,
                                                clusters_df,
                                                str_out)
    #  shift to avoid negative values
    scores_merged.X = csr_matrix(
        scores_merged.X.todense() + abs(scores_merged.X.min().min()))

    differential_results_by_cluster = squidpy_differential_analysis(
        scores_merged, group_str=str_out)

    top_differential = differential_results_filter(
        differential_results_by_cluster, abs_lfc=0.01, #  relax abs_lfc to allow KEGG pathways to emerge
        pvals_adj=0.05, n_top_per_cluster=n_top_per_cluster, total_n_features_out = 300)

    if cluster_not_tissue is not None:
        top_differential = top_differential.loc[~top_differential[
            'group'].astype(int).isin(cluster_not_tissue), :].copy()

    scores_merged = measures_methods.compute_moran_i(scores_merged)

    scores_merged.var['I'] = scores_merged.uns['moranI'].loc[
        scores_merged.var.index.tolist(), 'I'].astype(float).tolist()

    # moran filter variable for pathways and families (stringent only for families)
    selected_pathways = pathways_extractor(scores_merged, top_differential)
    selected_families = families_extractor(scores_merged, top_differential)

    selected_sets = list(set(selected_families).union(set(selected_pathways)))

    top_differential = top_differential.loc[
        top_differential['names'].isin(selected_sets), :]


    top_differential = pd.merge(
        top_differential,
        scores_merged.var[['pathway_id', 'pathway_name', 'I']].copy(),
        how='left', left_on="names", right_on="pathway_id")

    return top_differential, scores_merged.uns['current_pathways_df']



if __name__ == '__main__':

    figures_centralized_folder = "../../../extra-figures"
    SAMPLE = "Cerebellum-mouse"  #"Brain-healthy-mouse"

    pathways_types = ['KEGG', 'WIKI']
    clustering_parameters= {"Brain-healthy-mouse" : {
                                'choice': 'kmeans', 'numeric_parameter': 6},
                            "Cerebellum-mouse" : {
                                'choice': 'kmeans', 'numeric_parameter': 4}}

    subfolder = SAMPLE
    data_centralized_folder = "../../../apollo-data"

    cluster_not_tissue = {"Brain-healthy-mouse": [0],
                          "Cerebellum-mouse" : [1]}

    # each of these families were inspected, they were less relevant for publication
    sets_to_hide = { "Brain-healthy-mouse": ['MolFam_5', 'MolFam_8', 'MolFam_10',
                                   'MolFam_12','MolFam_17', 'MolFam_13',
                                     'MolFam_18',
                                    'hsa05231', 'hsa04723', 'WP4313'],
                     "Cerebellum-mouse" : ['MolFam_14','MolFam_13',
                                           'MolFam_19', 'MolFam_19',
                                           'MolFam_24', 'MolFam_23', 'MolFam_21',
                                           'MolFam_18', 'MolFam_16',
                                           'MolFam_6', 'MolFam_7',
                                            'MolFam_2', 'MolFam_1']
    }

    try:
        os.makedirs(os.path.join(figures_centralized_folder, subfolder))
    except OSError:
        pass

    df, mapping_info = modified_top_selection(
        pathways_types, clustering_parameters[SAMPLE],
        SAMPLE, subfolder, data_centralized_folder,
        cluster_not_tissue=cluster_not_tissue[SAMPLE])

    df = df.loc[~df['names'].isin(sets_to_hide[SAMPLE]),:]

    df = df.sort_values(by='names', ascending=False)

    diff_mz = pd.read_csv(os.path.join(data_centralized_folder, subfolder,
                                       f'{SAMPLE}-top-differential-mzs.tsv'),
                          sep='\t', header=0)


    diff_mz = diff_mz.sort_values(['group', 'logfoldchanges'], ascending=[True,False])
    diff_mz = diff_mz.groupby('group').head(30).reset_index(drop=True).copy()
    # note: also did expanded version using .head(80) for supplementary material

    run_alluvial_plot(SAMPLE, df, mapping_info, diff_mz,
    figures_centralized_folder, subfolder)
