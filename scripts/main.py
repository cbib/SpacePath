import sys
import os
import pandas as pd
from typing import Dict, List, Tuple
from analysis_workflow import run_analysis_one_tissue

"""
Running using a table listing the datasets and their parameters,
 for the analyses presented in the paper.

Note that this code is an example, users can adapt
 the way the workflow is applied to their data.
"""

def extract_params_per_dataset(row
)-> Tuple[
    str, Dict[str, str], float, float, float, float, float,
    Dict[str, str|int|float],
    List[str], int, float, int
]:
    SAMPLE = row["dataset_name"]
    annotations_to_use = row["annotations_to_use"].split(" ")
    annotations_dict = annotation_types_to_file_names_dict(SAMPLE, annotations_to_use)
    offSample_cutoff = row["offSample_cutoff"]
    msm_cutoff = row["msm_cutoff"]
    fdr_cutoff = row["fdr_cutoff"]
    moran_I_cutoff = row["Moran_I_cutoff"]
    moran_padj_cutoff = row["Moran_padj_cutoff"]
    clustering_params = row['cluster_params'].split(":")
    clustering_parameters = {'choice': clustering_params[0].lower(),
                             'numeric_parameter': clustering_params[1]}
    pathways_types = row["pathways_types"].split(" ")
    min_pathway_size = row["min_pathway_size"]
    cutoff_tree = row['cut_tree_mol_families']
    min_family_size = row["min_family_size"]

    return (SAMPLE, annotations_dict, offSample_cutoff,
            msm_cutoff, fdr_cutoff, moran_I_cutoff,
            moran_padj_cutoff, clustering_parameters,
            pathways_types, min_pathway_size, cutoff_tree,
            min_family_size)


def annotation_types_to_file_names_dict(SAMPLE, annotation_types_list):
    out_d = {}
    for annotation_type in annotation_types_list:
        out_d[annotation_type] = f"{SAMPLE}-{annotation_type}.csv"

    return out_d

# start

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

    run_analysis_one_tissue(annotations_dict, offSample_cutoff,
                            msm_cutoff, fdr_cutoff,
                            moran_I_cutoff, moran_padj_cutoff,
                            clustering_parameters,
                            pathways_types, min_pathway_size, cutoff_tree,
                            min_family_size, SAMPLE,
                            data_centralized_folder, figures_centralized_folder)
