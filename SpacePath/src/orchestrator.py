# Orchestration of operations for the CLI interface
# This also demonstrates how to use the SpacePath Workflow as an API
# johaGL 2026
import gc
import warnings

warnings.simplefilter('ignore', category=FutureWarning)

from SpacePath.api import *  # all the required callable operations
import matplotlib as mpl
import pandas as pd
from typing import Dict, List, Tuple, Any

logger = logging.getLogger(__name__)
setup_logger(logger)

np.random.seed(42)  # important to make the computations be reproducible


class SingleTissueWorkflow:
    """
    Class that wraps the SpacePath analysis for one single tissue.
    (tissue and sample are exchangeable terms)
    """

    def __init__(self, sample: str,
                 annotations_dict: Dict[str, str],
                 offSample_cutoff: float, msm_cutoff: float,
                 fdr_cutoff: float,
                 moran_I_cutoff: float, moran_padj_cutoff: float,
                 clustering_parameters: Dict[str, Any],
                 pathways_types: List[str], min_pathway_size: int,
                 cutoff_tree: float, min_family_size: int,

                 data_dir: str, figures_dir: str, mode: str):
        """
        Initialize the  class.
        """
        self.sample = sample
        self.annotations_dict = annotations_dict
        self.offSample_cutoff = offSample_cutoff
        self.msm_cutoff = msm_cutoff
        self.fdr_cutoff = fdr_cutoff
        self.moran_I_cutoff = moran_I_cutoff
        self.moran_padj_cutoff = moran_padj_cutoff
        self.clustering_parameters = clustering_parameters
        self.pathways_types = pathways_types
        self.min_pathway_size = min_pathway_size
        self.min_family_size = min_family_size
        self.cutoff_tree = cutoff_tree
        self.min_family_size = min_family_size

        self.data_dir = data_dir
        self.figures_dir = figures_dir
        self.subfolder = sample
        self.mode = mode

        self.plot_folder = os.path.join(figures_dir, sample)

        # Create plot folder if it doesn't exist
        if not os.path.exists(self.plot_folder):
            os.makedirs(self.plot_folder)

    def run(self) -> None:
        """
        Step by step analysis for a single specified tissue slice dataset.
        """
        logger.info('\nRunning sample: %s', self.sample)

        ##### these three lines are mandatory to comply with Nature journals #####

        mpl.rcParams['font.family'] = "sans-serif"
        mpl.rcParams['font.sans-serif'] = "Arial"

        if self.mode in ["all", "1", "2"]:
            annots_obj = spd.MetabolomeAnnotations(
                name=self.sample, data_dir=self.data_dir,
                subfolder=self.subfolder,
                annotation_files=self.annotations_dict)
            annots_obj.load()  # load annotations with core library

        if self.mode in ["all", "1"]:
            plt.rcParams.update({'font.size': 22})
            normalize_classify_and_save_anndata_mz(
                annots_obj, self.offSample_cutoff,
                self.msm_cutoff, self.fdr_cutoff, self.moran_I_cutoff,
                self.moran_padj_cutoff, self.sample, self.subfolder,
                self.data_dir, self.plot_folder)

        if self.mode in ["all", "2"]:
            clean_annotations_and_save_tsv(annots_obj, self.pathways_types,
                                           self.sample, self.subfolder,
                                           self.data_dir
                                           );
            del annots_obj  # obj no longer used: use clean .tsv files from now

        if self.mode in ["all", "3"]:
            plt.rcParams.update({'font.size': 20})
            compute_clustering_of_spots(self.pathways_types,
                                        self.clustering_parameters,
                                        self.sample, self.subfolder,
                                        self.data_dir, self.plot_folder)

        if self.mode in ["all", "4"]:
            compute_differential_analysis_of_mzs(self.pathways_types,
                                                 self.clustering_parameters,
                                                 self.sample, self.subfolder,
                                                 self.data_dir)

        if self.mode in ["all", "5"]:
            try:
                compute_pathway_scores(self.pathways_types,
                                       self.min_pathway_size,
                                       self.sample, self.subfolder,
                                       self.data_dir)
            except Exception as e:
                logger.warning("failed compute_pathway_scores, %s", e)

        if self.mode in ["all", "6"]:
            try:
                plt.rcParams.update({'font.size': 14})
                compute_molecules_families(self.cutoff_tree,
                                           self.sample, self.subfolder,
                                           self.data_dir, self.plot_folder)
            except Exception as e:
                logger.warning("failed compute_molecules_families , %s", e)

        if self.mode in ["all", "7"]:
            try:
                compute_families_scores(self.min_family_size,
                                        self.sample, self.subfolder,
                                        self.data_dir)
            except Exception as e:
                logger.warning("failed compute_families_scores , %s", e)

        if self.mode in ["all", "8"]:
            try:
                compute_differential_families_and_pathways(
                    self.pathways_types,
                    self.clustering_parameters, self.sample, self.subfolder,
                    self.data_dir)
            except Exception as e:
                logger.warning("failed differential families %s", e)

        logger.info("processed sample %s \n\n", self.sample)

    # end for


# end


class DatasetsProcessor:
    def __init__(self, input_file: str, data_dir: str, figures_dir: str,
                 mode: str):
        """
        Initialize the SpacePathAnalysisMultipleDatasets class,
        that handles the analysis of multiple datasets.

        Args:
            input_file (str): Path to the input file containing datasets specifics.
            data_dir (str): Directory where data files are stored.
            figures_dir (str): Directory where figures will be saved.
            mode (str): "all" or the string corresponding to the step to run.
        """
        self.input_file = input_file
        self.data_dir = data_dir
        self.figures_dir = figures_dir
        self.mode = mode
        self.datasets_parameters_df = pd.read_csv(input_file, header=0,
                                                  index_col=None, sep='\t')

    @staticmethod
    def annotation_types_to_file_names_dict(
            sample: str, annotation_types_list: List[str]) -> Dict[str, str]:
        """
        Create a dictionary of annotation files with standardized file names.

        Args:
            sample (str): The sample name.
            annotation_types_list (List[str]): List of annotation types.

        Returns:
            Dict[str, str]: Dictionary mapping annotation types to file names.
        """
        annot_dict = {annotation_type: f"{sample}-{annotation_type}.csv"
                      for annotation_type in annotation_types_list}
        return annot_dict

    @staticmethod
    def extract_params_per_dataset(row: pd.Series
                                   ) -> Tuple[
        str, Dict[str, str], float, float, float, float, float,
        Dict[str, Any], List[str], int, float, int]:
        """
        Extract parameters from a row of the input dataframe.

        Args:
            row (pd.Series): A row from the input dataframe.

        Returns:
            Tuple: Extracted parameters.
        """
        sample = row["dataset_name"]
        annotations_to_use = row["annotations_to_use"].split(" ")
        annotations_dict = DatasetsProcessor.annotation_types_to_file_names_dict(
            sample, annotations_to_use)
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

        return (sample, annotations_dict, offSample_cutoff, msm_cutoff,
                fdr_cutoff, moran_I_cutoff,
                moran_padj_cutoff, clustering_parameters, pathways_types,
                min_pathway_size, cutoff_tree, min_family_size)

    def process_tissues(self) -> None:
        """
        Process each dataset according to the specified mode.
        """
        for _, row in self.datasets_parameters_df.iterrows():
            (sample, annotations_dict, offSample_cutoff, msm_cutoff,
             fdr_cutoff, moran_I_cutoff,
             moran_padj_cutoff, clustering_parameters, pathways_types,
             min_pathway_size, cutoff_tree,
             min_family_size) = self.extract_params_per_dataset(row)

            analysis_one_tissue = SingleTissueWorkflow(
                sample, annotations_dict,
                offSample_cutoff, msm_cutoff, fdr_cutoff,
                moran_I_cutoff, moran_padj_cutoff,
                clustering_parameters, pathways_types, min_pathway_size,
                cutoff_tree, min_family_size,
                self.data_dir, self.figures_dir, self.mode
            )
            analysis_one_tissue.run()

            gc.collect()
