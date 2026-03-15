# johaGL 2025
# this is method/__init__
import logging
import os
from pathlib import Path
from typing import Dict, List, Tuple, Union

import anndata as ad
import numpy as np
import pandas as pd

from pydantic import BaseModel as PydanticBaseModel
from pydantic import ConfigDict

from rdkit import DataStructs
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.sparse import csr_matrix
from scipy.spatial.distance import squareform
import smpath.helpers as helpers
import smpath.processing.measures_methods as measures_methods
import smpath.processing.pathways_methods as pathways_methods
import smpath.visualization.visualization as sviz
from smpath.data import MetabolicPathways, MetabolomeAnnotations
import yaml

logger = logging.getLogger(__name__)
logger = helpers.setup_logger(logger)

module_path = Path(__file__).resolve()
library_path = module_path.parent.parent.parent.parent


class BaseModel(PydanticBaseModel):
    # class Config:  # deprecated since pydantic v2
    #     arbitrary_types_allowed = True
    model_config = ConfigDict(arbitrary_types_allowed=True)


# 1) classes of the configurations

class MethodConfig(BaseModel):
    class_type: str = "method"


class NormalizeIntensitiesConfig(MethodConfig):
    description: str = "normalization of ions (m/z) intensities"
    method: str = "TotalIonCount"


class ClassifyMZProtocolConfig(MethodConfig):
    description: str = ("protocol to classify the ions (m/z) by the quality "
                        "of molecular formulae annotations, "
                        "and the spatial distributions")
    n_neighs_for_graph: int = 20  # for neighbors graph calc before Moran's I
    n_perms_moran_calc: Union[int, str] = "auto"  # for Moran's I
    msm_cutoff: float = 0.5
    offSample_cutoff: float = 0.2
    fdr_cutoff: float = 0.2,
    moran_I_cutoff: float = 0.1
    moran_padj_cutoff: float = 0.1
    var_sim_cutoff: float = 1e-20


class ComputePathwaysScoresConfig(MethodConfig):
    description: str = "parameters to compute the pathways scores"
    method: str = "sspa_kpca"
    type_abbr: str = ""  # WIKI, KEGG or FAMILIES are to date supported
    min_entity: int = 2
    n_groups_of_spots: int = 30
    threshold_n_spots_is_high: int = 15000


class ComputeMoleculesSimilarityMatrixConfig(MethodConfig):
    description: str = "compute similarity matrix of compounds"
    method: str = "TanimotoSimilarity"
    smiles_files: Dict[str, str] = {
        'HMDB': "HMDB-SMILES_05-07-24.tsv",
        'LipidMaps': "LipidMaps-SMILES_06-08-25.tsv"}  # internal in library
    threshold_similarity: float | None = None


class ComputeMoleculesFamiliesConfig(MethodConfig):
    description: str = ("compute and save families file, "
                        "from a distance-transformed similarity matrix")
    method: str = "HierarchicalClustering"
    cutoff_tree: float = 1
    method_linkage: str = "ward"
    metric_linkage: str = "cosine"


# 2) the methods using the configs above


class NormalizeIntensities:
    config: NormalizeIntensitiesConfig

    @staticmethod
    def run(adata_mz: ad.AnnData, config=NormalizeIntensitiesConfig()
            ) -> ad.AnnData:
        """
        Updates the adata_mz object:
         adata_mz.X with the normalized result
         adata_mz.layers['raw_X'] with the raw data
         adata_mz.uns['norm'] with the flag of the applied method
        """
        if config.method == "TotalIonCount":
            adata_mz.uns['tmp_slot'] = 0  # guarantees `.uns` is accessible
            try:
                if 'norm' not in list(adata_mz.uns.keys()):
                    logger.info("%s", config.description)
                    adata_mz.layers['raw_X'] = adata_mz.X.copy()
                    normalised_X = measures_methods.tic_normalization(
                        adata_mz.X.copy())
                    adata_mz.X = csr_matrix(normalised_X)
                    adata_mz.uns['norm'] = config.method
                else:
                    logger.info("Intensities already normalised, "
                                "check .uns['norm'] and .layers['raw_X']")
            except Exception as e:
                logger.critical("Error %s \ncould not normalise", e)
        else:
            logger.error("%s is not supported", config.method)

        del adata_mz.uns['tmp_slot']

        return adata_mz


class ClassifyMZProtocol:
    config: ClassifyMZProtocolConfig

    @staticmethod
    def run(adata_mz: ad.AnnData, annotations_obj: MetabolomeAnnotations,
            SAMPLE: str, plot_folder: str, config=ClassifyMZProtocolConfig()
            ) -> ad.AnnData:
        logger.info("Running %s", config.description)
        logger.info(f">> Part a: classifying the ions with thresholds: "
                    f"MSM >= %s, offSample <= %s, FDR <= %s",
                    config.msm_cutoff,
                    config.offSample_cutoff, config.fdr_cutoff)
        annotations_obj, mz_good_from_part_a = measures_methods.part_a__mz_quality_mark(
            annotations_obj, config.msm_cutoff,
            config.offSample_cutoff, config.fdr_cutoff)

        logger.info(
            f"  a total of %s ions pass MSM, offSample and FDR cutoffs.",
            len(set(mz_good_from_part_a)))
        text_notice = f"valid $m/z$\n n = {len(mz_good_from_part_a)}"
        sviz.plot_mz_given_quality_metrics(annotations_obj.mz_quality, SAMPLE,
                                           plot_folder, config.msm_cutoff,
                                           config.offSample_cutoff,
                                           text_notice)
        logger.info("  saved plot (MSM, offSample, FDR) to %s", plot_folder)
        adata_mz.uns["history"]['features_n_after_mz_clean_1'] = len(
            mz_good_from_part_a)

        if config.n_perms_moran_calc == "auto":
            config.n_perms_moran_calc = 50
            if (adata_mz.shape[0] > 9e4) and (adata_mz.shape[1] > 1e3):
                config.n_perms_moran_calc = 10

        logger.info(f">> Part b: classifying ions based on Moran's I: "
                    f"thresholds: I >= %s,  (BH) adjusted p-value <= %s",
                    config.moran_I_cutoff, config.moran_padj_cutoff)
        logger.info(
            "  metrics are only computed for the ions that passed the 'Part a'.")
        all_quality_df, mz_good_part_b, moran_mz_df = measures_methods.part_b__mz_quality_by_autocorr(
            adata_mz, mz_good_from_part_a, annotations_obj,
            config.moran_I_cutoff, config.moran_padj_cutoff,
            config.var_sim_cutoff,
            config.n_neighs_for_graph,
            config.n_perms_moran_calc)

        text_notice = f"valid $m/z$: n = {len(mz_good_part_b)}"
        sviz.plot__moran_info(moran_mz_df, plot_folder,
                              SAMPLE,
                              text_notice, config.moran_I_cutoff,
                              config.moran_padj_cutoff)
        logger.info("  saved plot (autocorrelation) to %s", plot_folder)

        adata_mz.uns["history"]['features_n_after_mz_clean_2'] = len(
            mz_good_part_b)

        adata_mz = measures_methods.adata_mz_insert_full_quality(adata_mz,
                                                                 moran_mz_df,
                                                                 all_quality_df)
        n_good_mz = adata_mz.var.loc[
            adata_mz.var['final_quality'] == "good", :].shape[0]
        logger.info("A total of %s ions exhibit good final quality "
                    "(check column 'final_quality' of adata_mz.var)",
                    n_good_mz)

        return adata_mz


class SyncPathwaysAndAnnots:
    @staticmethod
    def run(pathwaysObj: MetabolicPathways,
            annotations_df: pd.DataFrame,
            usable_id: str = "",
            mz_col: str = "mz",
            min_set_size: int = 2,
            max_set_size: int = 200,
            mode: str = "mz_aware"  # or "default", but mz_aware is preferred
            ) -> Tuple[MetabolicPathways, pd.DataFrame]:
        if not usable_id == pathwaysObj.molecule_id_type:
            logger.error(
                "No same types of compound ids: %s and %s, respectively",
                pathwaysObj.molecule_id_type, usable_id)
            raise AssertionError("see the error in .log file. Aborted.")
        if max_set_size == 1:
            mode = "default"
        pathwaysObj.intersect_with_annotations(annotations_df, usable_id,
                                               mz_col)  # this updates pathwaysObj only
        pathwaysObj.constrain(min_set_size, max_set_size, mode)

        annotations_df = annotations_df.loc[
            annotations_df[usable_id].isin(
                pathwaysObj.piled[pathwaysObj.molecule_id_type].unique()), :]

        logger.info("Pathways and Annotations objects synchronized, with"
                    " min set size = %s", min_set_size)

        return pathwaysObj, annotations_df


class ComputeAssignmentsWeights:
    @staticmethod
    def run(adata_mz: ad.AnnData, clean_annotations: pd.DataFrame,
            usable_id: str) -> Tuple[pd.DataFrame, np.ndarray]:
        adata_mz = adata_mz[
            :, clean_annotations['mz'].astype(str).unique()].copy()
        X_df = measures_methods.compute_metabolites_probs_df(
            adata_mz, clean_annotations, usable_id)
        adataIds = measures_methods.build_identifiers_weights_df(
            adata_mz, X_df, clean_annotations, usable_id)
        coords = adata_mz.obsm['spatial'].copy()

        return adataIds, coords


class ComputePathwaysScores:
    config: ComputePathwaysScoresConfig

    @staticmethod
    def run(pathways_dict: Dict[str, List[str]],
            metabolites_computed_w: pd.DataFrame,
            config=ComputePathwaysScoresConfig()
            ) -> pd.DataFrame:
        logger.info("\n Will compute spot-wise Pathways Scores "
                    "with the following config: %s",
                    yaml.dump(config.model_dump()))
        if config.method == "sspa_kpca":
            logger.info(" using sspa library, credits: Wieder C et al., 2022")
        scores_df = pathways_methods.compute_pathways_scores(
            pathways_dict, metabolites_computed_w,
            config.method, config.min_entity,
            config.n_groups_of_spots, config.threshold_n_spots_is_high)

        logger.info("successfully computed scores for the specified sets"
                    " (abbr type): %s", config.type_abbr)

        return scores_df


class ComputeMoleculesSimilarityMatrix:
    config: ComputeMoleculesSimilarityMatrixConfig

    @staticmethod
    def run(molecules_ids_list: List[str],
            config=ComputeMoleculesSimilarityMatrixConfig()
            ) -> Dict[str, np.array | List[str] | List[str]]:
        logger.info("computing similarity matrix with config %s",
                    yaml.dump(config.model_dump()))
        molecules_ids_list = list(set(molecules_ids_list))
        smiles_dfs = []
        pathways_methods.unzip_if_not_tsv(library_path, config.smiles_files)
        for file_name in config.smiles_files.values():
            df = pd.read_csv(
                os.path.join(library_path, "resources", "SMILES", file_name),
                sep='\t', header=0)
            smiles_dfs.append(df)

        smiles_df = pd.concat(smiles_dfs)
        smiles_df = smiles_df.loc[~smiles_df['SMILES'].isnull(), :].copy()

        (identifiers_validated, smiles_validated
         ) = pathways_methods.compounds_and_smiles_ordered(
            smiles_df, molecules_ids_list)

        identifiers_validated, smiles_validated, fps = helpers.secured_trios(
            identifiers_validated, smiles_validated
        )

        arr = np.zeros(shape=(len(identifiers_validated),
                              len(identifiers_validated)))

        logger.info("computing similarities for all pairs of SMILES")
        for u in range(len(identifiers_validated)):
            for v in range(u + 1, len(identifiers_validated)):
                tanimoto_similarity = DataStructs.TanimotoSimilarity(
                    fps[u], fps[v])
                arr[u, v] = tanimoto_similarity
                arr[v, u] = tanimoto_similarity

        np.fill_diagonal(arr, 1)
        if config.threshold_similarity is not None:
            arr[arr < config.threshold_similarity] = 0

        mol_similarities_result = {"array": arr,
                                   "identifiers": identifiers_validated,
                                   "smiles": smiles_validated}

        return mol_similarities_result


class ComputeMoleculesFamilies:
    config: ComputeMoleculesFamiliesConfig

    @staticmethod
    def run(mol_similarities_result: Dict[
        str, np.array | List[str] | List[str]],
            molecules_ids_df: pd.DataFrame,
            SAMPLE: str, subfolder: str, mydata_dir: str, plot_folder: str,
            config=ComputeMoleculesFamiliesConfig()
            ) -> Dict[str, np.ndarray | pd.DataFrame | float]:
        logger.info("computing families of molecules, with config %s",
                    yaml.dump(config.model_dump()))
        # make the hierarchical clustering reproducible
        np.random.seed(0)

        distance_matrix = 1 - mol_similarities_result['array']
        # squareform saves memory by representing the matrix with the upper right triangle
        computed_dendrogram = linkage(squareform(distance_matrix),
                                      method=config.method_linkage,
                                      # it is 'ward' by default
                                      metric=config.metric_linkage
                                      # it is cosine by default
                                      )
        clusters = fcluster(computed_dendrogram,
                            t=config.cutoff_tree, criterion='distance')

        id_to_cluster = pd.DataFrame(
            {'moleculeIds': mol_similarities_result['identifiers'],
             'molecules_cluster': clusters})

        unique_mz_nb = len(molecules_ids_df.loc[
                               molecules_ids_df['moleculeIds'].isin(
                                   mol_similarities_result[
                                       'identifiers']), 'mz'].unique())

        logger.info("the clusters are: %s",
                    sorted(
                        id_to_cluster['molecules_cluster'].unique().tolist()))
        logger.info("number of used compound identifiers: %s",
                    len(mol_similarities_result['identifiers']))
        logger.info("representing distinct mzs: %s", unique_mz_nb)

        pathways_methods.save_families_as_mmt_file(
            id_to_cluster, molecules_ids_df, SAMPLE, subfolder, mydata_dir
        )
        logger.info("saved .mmt file to %s %s", mydata_dir, subfolder)

        fig_obj = sviz.plot_families(molecules_ids_df,
                                     mol_similarities_result,
                                     computed_dendrogram)
        sviz.save_families_plot(fig_obj, SAMPLE, plot_folder)
        logger.info("saved figure to %s", plot_folder)

        return {"computed_dendrogram": computed_dendrogram,
                "id_to_cluster_df": id_to_cluster,
                "cutoff_tree": config.cutoff_tree}
