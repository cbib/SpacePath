from .tissue_workflow import *
import smpath.data as spd
import warnings
warnings.simplefilter('ignore',category=FutureWarning)
import squidpy as sq


"""
All this workflow uses the library smpath
"""

def annotation_types_to_file_names_dict(SAMPLE, annotation_types_list):
    out_d = {}
    for annotation_type in annotation_types_list:
        out_d[annotation_type] = f"{SAMPLE}-{annotation_type}.csv"

    return out_d


def extract_params_per_dataset(row
)-> Tuple[str, Dict[str, str], float, float, float, float, float,
    Dict[str, str|int|float], List[str], int, float, int]:
    """
    extracts the parameters from the current row
    (for one dataset), read from the input dataframe
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


def run_analysis_one_tissue(annotations_dict: Dict[str, str],
                            offSample_cutoff: float ,
                            msm_cutoff: float, fdr_cutoff: float,
                            moran_I_cutoff: float, moran_padj_cutoff: float,
                            clustering_parameters: Dict[str, str|int|float],
                            pathways_types: List[str], min_pathway_size: int,
                            cutoff_tree: float,
                            min_family_size: int,  SAMPLE: str,
                            data_centralized_folder: str,
                            figures_centralized_folder: str)-> None:
    """runs step by step analysis for a given dataset with its parameters"""
    subfolder = SAMPLE
    plot_folder = os.path.join(figures_centralized_folder, SAMPLE)
    if not os.path.exists(plot_folder):
            os.mkdir(plot_folder)

    annots_obj = spd.MetabolomeAnnotations(
        name=SAMPLE, data_dir=data_centralized_folder,
        subfolder=subfolder, annotation_files=annotations_dict)
    annots_obj.load()  # load annotations with smpath library

    normalize_classify_and_save_anndata_mz(
        annots_obj, offSample_cutoff, msm_cutoff, fdr_cutoff, moran_I_cutoff,
        moran_padj_cutoff, SAMPLE, subfolder,
        data_centralized_folder, plot_folder)

    clean_annotations_and_save_tsv(annots_obj, pathways_types, SAMPLE, subfolder,
                                   data_centralized_folder); del annots_obj  # obj no longer used: use clean .tsv files from now

    compute_clustering_of_spots(pathways_types, clustering_parameters,
        SAMPLE, subfolder, data_centralized_folder, plot_folder)

    compute_differential_analysis_of_mzs(pathways_types,
         clustering_parameters, SAMPLE, subfolder, data_centralized_folder)

    compute_pathway_scores( pathways_types, min_pathway_size, SAMPLE,
         subfolder, data_centralized_folder)

    compute_molecules_families(cutoff_tree,
           SAMPLE, subfolder, data_centralized_folder, plot_folder)

    compute_families_scores( min_family_size,
          SAMPLE, subfolder, data_centralized_folder)

    compute_differential_families_and_pathways(pathways_types,
        clustering_parameters, SAMPLE, subfolder, data_centralized_folder)

    # todo the folder metabolites-mappings in the exampleData !!! &todo for documentation: explain the structure of data/ folder,


def process_tissues(input_file:str, data_dir:str, figures_dir:str)-> None:
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
