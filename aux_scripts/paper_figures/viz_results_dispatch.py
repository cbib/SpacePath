
import sys
import os
import warnings
from pathlib import Path
import anndata as ad

#__path__ = __import__('../../src').extend_path(__path__, __name__)

sys.path.append(str(Path(__file__).resolve().parents[2]))

from src.single_slice_helpers import join_scores_adatas
from src.process_single_slice import (join_clean_annotations_from_dict,
                                      safe_load_data_and_clean_annotations
                                      )
from viz_results_helpers import run_alluvial_plot, run_individual_mz_scatters, run_individual_tabular_plots


from src.dispatch_tissues_slices import extract_params_per_dataset

warnings.simplefilter('ignore',category=FutureWarning)


import numpy as np
import pandas as pd

import squidpy as sq
from matplotlib.axes import Axes
from typing import Union

import matplotlib.pyplot as plt
from matplotlib.figure import Figure


np.random.seed(42)  # important to make the computations be reproducible


def filter_sets_per_cluster_for_plot(diff_copy_df, cluster_to_check,
                                     top_for_plotting=2):
    """
    selects differential pathways being the n top by lfc, for plotting
       forces selection by positive and negative log2FCs, to avoid having
       all down or all up.

    Note that 'down' or 'up' are theoretical, not yet biologically verified
    """
    diff_copy_df = diff_copy_df.loc[
                 diff_copy_df['group'].astype(str) == str(cluster_to_check), :].copy()

    up = diff_copy_df.sort_values(
        by=['logfoldchanges', 'pvals_adj'],
        ascending=[False, True]).head(top_for_plotting).copy()
    down = diff_copy_df.sort_values(
        by=['logfoldchanges', 'pvals_adj'],
        ascending=[True, True]).head(top_for_plotting).copy()

    diff_copy_df = pd.concat([up, down], axis=0)

    return diff_copy_df



def run_visualization(clustering_parameters,
                      pathways_types, min_pathway_size, cutoff_tree,
                      min_family_size, SAMPLE,
                      data_centralized_folder,
                      figures_centralized_folder, sets_relevant_n_limit,
                      mz_spatially_relevant_n_limit,
                      mode):
    subfolder = SAMPLE
    print("\n", SAMPLE)
    str_out = f"{clustering_parameters['choice']}-{clustering_parameters['numeric_parameter']}"

    diff_sets = pd.read_csv(os.path.join(data_centralized_folder, subfolder,
                            f'{SAMPLE}-top-differential-SETS.tsv'),
                            sep='\t', header=0)
    diff_sets = diff_sets.sort_values(
        by='I', ascending=False).head(sets_relevant_n_limit).copy()

    # print(diff_sets.head())
    scores_concatenated = join_scores_adatas(pathways_types, SAMPLE, subfolder,
                                data_centralized_folder)

    mapping_info = scores_concatenated.uns['current_pathways_df']
    # print(mapping_info.head())
    """ 
             pathway_id   moleculeIds          mz
        0     WP1403   HMDB0003540  461.165659
        1     WP1403   HMDB0000122  294.165962
        2      WP167  LMFA03010021  430.306419
        3      WP167  LMFA03020023  432.322069
        """
    #print("clusters : ", diff_sets['group'].unique().tolist())
    ## check how many unique pathways come out in general
    print("number of unique pathways that passed significance filters:")
    print(len(diff_sets['group'].unique().tolist()))

    diff_mz = pd.read_csv(os.path.join(data_centralized_folder, subfolder,
                          f'{SAMPLE}-top-differential-mzs.tsv'),
                          sep='\t', header=0)
    pathways_types__for_annots = list(set(pathways_types) - {"FAMILIES"})

    #concatenated_clean_annotations = join_clean_annotations_from_files(
    #    pathways_types__for_annots, SAMPLE, subfolder, data_centralized_folder
    #)

    adata_mz, annots_dict = safe_load_data_and_clean_annotations(SAMPLE,
                                         pathways_types__for_annots,
                                         subfolder,
                                         data_centralized_folder)
    concatenated_clean_annotations = join_clean_annotations_from_dict(
        annots_dict)
    # print(concatenated_clean_annotations.head())
    """
         moleculeIds                                 moleculeNames
    0      C11284                          3-Indolebutyric acid
    1      C07359                            3-keto-Digoxigenin
    2      C04230                                  LysoPC(16:0)
    3      C05474  3a,11b,21-Trihydroxy-20-oxo-5b-pregnan-18-al
    4      C05471                               Dihydrocortisol
    """

    ## check how many unique mz come out in general
    print("unique mzs that passed significance filters:")
    print(len(diff_mz['mz'].unique().tolist()))

    ## check mzs both mapped and significant and with adequate moran I
    diff_mz_relevant = diff_mz.loc[ diff_mz['mz'].astype(float).isin(
        mapping_info['mz'].astype(float).tolist()
    ), : ].copy()

    # ============ only good moran I mzs  ===========
    # mzs_moran_good = adata_mz.var.loc[adata_mz.var['I'] >= 0.3, 'mz'].astype(float).tolist()
    # print("!!!!", len(mzs_moran_good))

    # adata_mz.var.loc[adata_mz.var['I'] >= 0.3, 'mz'].astype(float).tolist()
    # # Note: instead imposing a hard filter of I value, take the top 25 mzs

    top_moran_mzs = adata_mz.var.copy().sort_values(
        by='I', ascending=False).head(mz_spatially_relevant_n_limit).copy()

    """
             mz       msm  rawOffSampleProb   fdr      formula quality_msm_off_fdr         I  pval_sim_fdr_bh   var_sim quality_autocorr final_quality
760.548697  760.548697  0.528119          0.000296  0.50   C41H78NO9P                good  0.741693         0.019608  0.000172             good          good
512.45744   512.457440  0.941250          0.000016  0.05      C28H46O                good  0.689096         0.019608  0.000136             good          good
834.621861  834.621861  0.545073          0.000657  0.50  C45H88NO10P                good  0.669166         0.019608  0.000190             good          good
    """

    diff_mz_relevant = diff_mz_relevant.loc[diff_mz_relevant['mz'].astype(float).isin(
        top_moran_mzs['mz'].astype(float).tolist() ),:].copy()

    print("unique mzs that passed significance filters AND MAPPED AND good Moran:")
    print(len(diff_mz_relevant['mz'].unique().tolist()))
    ##                  do sankey diagrams of cluster and pathways and mz
    try:
        run_alluvial_plot(SAMPLE, diff_sets, mapping_info, diff_mz_relevant,
                      figures_centralized_folder, subfolder)
    except Exception as e:
        print(e, "failed run alluvial")

    plot_folder = os.path.join(figures_centralized_folder,
                               subfolder, "prominent")
    if not os.path.exists(plot_folder):
        os.makedirs(plot_folder)

    try:
        run_individual_mz_scatters(SAMPLE, diff_sets, mapping_info, diff_mz_relevant,
                               adata_mz,
                      figures_centralized_folder, subfolder)
    except Exception as e:
        print(e, "failed run scatters")

    try:
        run_individual_tabular_plots(
            SAMPLE, diff_sets, mapping_info, diff_mz_relevant,
             concatenated_clean_annotations,
             figures_centralized_folder, subfolder
         )
    except Exception as e:
        print(e, "failed run tabular plots")
    # # DONE: in a distinct funcs    ##  mission: >>>>>> plot the top mz plots alone


# use all above

def visualize_results_for_tissues(input_file, data_dir, figures_dir, mode=None):
    """
    Assumes that SPACEPATH was already computed across the specified datasets

      - parses the input dataframe having the specifics of the input datasets
            for each row representing a dataset.
      - runs the visualizations of the SPACEPATH results
        """
    data_centralized_folder = data_dir
    figures_centralized_folder = figures_dir

    datasets_parameters_df = pd.read_csv(input_file,
                                         header=0, index_col=None, sep='\t')

    for i, row in datasets_parameters_df.iterrows():
        # each row represents a dataset and its parameters
        (SAMPLE, annotations_dict, offSample_cutoff, msm_cutoff, fdr_cutoff,
         moran_I_cutoff, moran_padj_cutoff, clustering_parameters,
         pathways_types, min_pathway_size, cutoff_tree,
         min_family_size
         ) = extract_params_per_dataset(row)

        run_visualization(clustering_parameters,
                          pathways_types, min_pathway_size, cutoff_tree,
                          min_family_size, SAMPLE,
                          data_centralized_folder,
                          figures_centralized_folder, sets_relevant_n_limit=100,
                          mz_spatially_relevant_n_limit=100,  # todo : chekc if augmenting it is ok or not
                          mode=None)


# if __name__ == "__main__":

    # diff_sets = diff_sets.loc[
    #     diff_sets['group'].astype(str) == cluster_to_check, :].copy()
    # diff_sets = diff_sets.sort_values(by='abslfc', ascending=False).head(5)
    # print(diff_sets.shape)   # <<< it was not useful to filter by fold change

    #pid = diff_sets['pathway_id'].tolist()
    #pnm = diff_sets['pathway_name'].tolist()

    # pid = joined.var['pathway_id'].tolist()
    # pnm = [t for t in joined.var['pathway_name'].astype(str).to_numpy()]
