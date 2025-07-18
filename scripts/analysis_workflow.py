from analysis_quasidispatcher import *
import warnings
warnings.simplefilter('ignore',category=FutureWarning)
import squidpy as sq


"""
All this workflow uses the library smpath

"""
# the function that calls all of the above

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
    """runs every step for a given dataset and its parameters"""
    subfolder = SAMPLE
    plot_folder = os.path.join(figures_centralized_folder, SAMPLE)
    if not os.path.exists(plot_folder):
            os.mkdir(plot_folder)

    annots_obj = spd.MetabolomeAnnotations(
        name=SAMPLE, data_dir=data_centralized_folder,
        subfolder=subfolder, annotation_files=annotations_dict)
    annots_obj.load()  # load annotations

    normalize_classify_and_save_anndata_mz(annots_obj, offSample_cutoff,
                                           msm_cutoff, fdr_cutoff, moran_I_cutoff, moran_padj_cutoff, SAMPLE, subfolder,
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

    # todo for documentation: explain the structure of data/ folder,
    #  give the data/ folder with minimal data and with the folder metabolites-mappings !!!!

    # EL FIN !!! 