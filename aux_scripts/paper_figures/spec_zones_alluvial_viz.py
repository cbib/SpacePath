# script for selection of specific zones of the tissue for visualization
#  reuses already computed objects for the selected SAMPLE

import os
import pandas as pd
import numpy as np

from SpacePath.src.helpers import (join_scores_adatas,
                                   safe_load_data_and_clean_annotations)
from viz_results_dispatch import join_clean_annotations_from_dict

from viz_results_helpers import (run_alluvial_plot,
                                 run_individual_mz_scatters,
                                 run_individual_tabular_plots
                                 )

SAMPLE = "GB-x-lipids-1"
data_centralized_folder = "../../../apollo-data"
pathways_types = ['KEGG', 'WIKI']
clusters_to_keep = [1, 3, 4, 5]  # spots' clusters for viz for *this* data

# fine-tuned to make fam 10  appear as well:
sets_relevant_n_limit = 105
mz_spatially_relevant_n_limit = 210

figures_centralized_folder = "../../../extra-figures"

if not os.path.exists(os.path.join(figures_centralized_folder, SAMPLE)):
    os.makedirs(os.path.join(figures_centralized_folder, SAMPLE))


subfolder = SAMPLE
print("\n", SAMPLE)

diff_sets = pd.read_csv(os.path.join(data_centralized_folder, subfolder,
                        f'{SAMPLE}-top-differential-SETS.tsv'),
                        sep='\t', header=0)

tum_clus = diff_sets.loc[diff_sets['group'].astype(int) == 5, "pathway_id"].unique()
print(tum_clus)


diff_sets = diff_sets.loc[diff_sets['group'].astype(
    int).isin(clusters_to_keep), :].copy()

print(diff_sets['pathway_id'].unique())

diff_sets = diff_sets.sort_values(
    by='I', ascending=False).head(sets_relevant_n_limit).copy()


scores_concatenated = join_scores_adatas(pathways_types, SAMPLE, subfolder,
                            data_centralized_folder)

mapping_info = scores_concatenated.uns['current_pathways_df']
""" mapping info:
         pathway_id   moleculeIds          mz
    0     WP1403   HMDB0003540  461.165659
    1     WP1403   HMDB0000122  294.165962
    2      WP167  LMFA03010021  430.306419
    3      WP167  LMFA03020023  432.322069
    """
## check how many unique pathways come out in general
print("number of unique pathways that passed significance filters:")
print(len(diff_sets['group'].unique().tolist()))

diff_mz = pd.read_csv(os.path.join(data_centralized_folder, subfolder,
                      f'{SAMPLE}-top-differential-mzs.tsv'),
                      sep='\t', header=0)
pathways_types__for_annots = list(set(pathways_types) - {"FAMILIES"})

adata_mz, annots_dict = safe_load_data_and_clean_annotations(SAMPLE,
                                     pathways_types__for_annots,
                                     subfolder,
                                     data_centralized_folder)
concatenated_clean_annotations = join_clean_annotations_from_dict(
    annots_dict)
""" concatenated_clean_annotations:
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

# instead imposing a hard filter of I value, take the top n mzs:
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
# sankey diagrams of cluster and pathways and mz
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