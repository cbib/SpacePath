from .process_single_slice import *
import warnings
warnings.simplefilter('ignore',category=FutureWarning)



def annotation_types_to_file_names_dict(SAMPLE, annotation_types_list):
    out_d = {}
    for annotation_type in annotation_types_list:
        out_d[annotation_type] = f"{SAMPLE}-{annotation_type}.csv"

    return out_d


def extract_params_per_dataset(row)-> (
        str, Dict[str, str], float, float, float, float, float,
        Dict[str, str|int|float], List[str], int, float, int):
    """
    extracts the parameters from the current row
    (for one dataset) within the input dataframe
    """
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






def process_tissues(input_file:str, data_dir:str, figures_dir:str)-> None:
    """
    parses the input dataframe having the specifics of the datasets
     for each row representing a dataset, runs one analysis
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


        run_analysis_one_tissue(annotations_dict, offSample_cutoff,
                                msm_cutoff, fdr_cutoff,
                                moran_I_cutoff, moran_padj_cutoff,
                                clustering_parameters,
                                pathways_types, min_pathway_size, cutoff_tree,
                                min_family_size, SAMPLE,
                                data_centralized_folder,
                                figures_centralized_folder)
