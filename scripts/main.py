import sys
import os
import pandas as pd
from analysis_workflow import run_analysis_one_tissue, extract_params_per_dataset

"""
Running using a table listing the datasets and their parameters,
 for the analyses presented in the paper.

Note that this code is an example, users can adapt the way the 
workflow is applied to their data.
"""



datasets_params_file = sys.argv[1]
datasets_params_dir = sys.argv[2]

data_centralized_folder = sys.argv[3]
figures_centralized_folder = sys.argv[4]

datasets_parameters_df = pd.read_csv(
    os.path.join(datasets_params_dir, datasets_params_file ),
    header=0,  index_col=None, sep='\t' )

for i, row in datasets_parameters_df.iterrows():
    # each row represents a dataset and its parameters
    (SAMPLE, annotations_dict, offSample_cutoff, msm_cutoff, fdr_cutoff,
      moran_I_cutoff, moran_padj_cutoff, clustering_parameters,
      pathways_types, min_pathway_size, cutoff_tree,
      min_family_size
     ) = extract_params_per_dataset(row)

    run_analysis_one_tissue(SAMPLE, annotations_dict, offSample_cutoff,
                            msm_cutoff, fdr_cutoff,
                            moran_I_cutoff, moran_padj_cutoff,
                            clustering_parameters,
                            pathways_types, min_pathway_size, cutoff_tree,
                            min_family_size,
                            data_centralized_folder, figures_centralized_folder)