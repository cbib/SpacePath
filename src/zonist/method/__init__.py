# johaGL 2025
# this is method/__init__
import os
import logging
import numpy as np
from typing import List, Dict, Literal, Optional, Set, Tuple, Union

import pandas as pd
import anndata as ad

from pydantic import BaseModel as PydanticBaseModel
from pydantic import ConfigDict
from pydantic import Extra

from pathlib import Path

from rdkit import DataStructs
from rdkit import Chem
from rdkit.Chem import AllChem

from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform

import zonist.helpers as helpers
from zonist.data import MetabolicPathways, MetabolomeAnnotations
import zonist.processing.measures_methods as zmm
import zonist.processing.pathways_methods as pathways_methods
import zonist.visualization.visualization as zviz


logger = logging.getLogger(__name__)
logger = helpers.reset_log_config(logger)

module_path = Path(__file__).resolve()
library_path = module_path.parent.parent.parent.parent



class BaseModel(PydanticBaseModel):
    # class Config:  # deprecated since pydantic v2
    #     arbitrary_types_allowed = True
    model_config = ConfigDict(arbitrary_types_allowed=True)


# 1 methods configs
class MethodConfig(BaseModel):
    class_type: str = "method"


class ClassifyMZProtocolConfig(MethodConfig):
    description: str = "phased protocol to filter m/z by the quality of molecular formulae annotations, and the spatial distributions"
    n_neighs_for_graph: int = 20   # for neighbors graph computation before Moran's I
    n_perms_moran_calc: Union[int, str] = "auto"   # for Moran's I


class SyncPathwaysAndAnnotsConfig(MethodConfig): # todo: use this config ?
    description: str = "parameters to synchronize pathways -or families of compounds- and annotations"


class ComputePathwaysScoresConfig(MethodConfig):
    description: str = "parameters to compute the pathways scores"
    method: str = "sspa_kpca"
    min_entity: int = 2
    n_groups_of_spots: int = 30
    threshold_n_spots_is_high: int = 15000   # todo : in the stable version use 20K


class ComputeMoleculesSimilarityMatrixConfig(MethodConfig):
    description : str = "...."
    #method: str=""


class ComputeMoleculesFamiliesConfig(MethodConfig):
    description : str = "...."
    #method: str=""


def compounds_and_smiles_ordered(smiles_df, molecules_ids_list):
    mycopy = smiles_df.copy()
    mycopy = mycopy.drop_duplicates(subset=['moleculeIds', 'SMILES'],
                                    keep='first')
    the_identifiers_allowed = list(
        set(mycopy['moleculeIds'].unique().tolist()).intersection(
            set(molecules_ids_list))
    )
    mycopy = mycopy[['moleculeIds', 'SMILES']].copy()
    mycopy.index = mycopy['moleculeIds'].tolist()
    mycopy = mycopy.T
    mycopy = mycopy[list(set(the_identifiers_allowed))]
    mycopy = mycopy.T
    mycopy = mycopy.drop_duplicates()

    return mycopy['moleculeIds'].tolist(), mycopy['SMILES'].tolist()

class ComputeMoleculesSimilarityMatrix:
    config = ComputeMoleculesSimilarityMatrixConfig
    @staticmethod
    def run(molecules_ids_list, threshold_similarity=None):
        smiles_dfs  = []
        for file_name in ["LipidMaps-SMILES_06-08-25.tsv",
                          "HMDB-SMILES_05-07-24.tsv"]:
            df = pd.read_csv(
                os.path.join(library_path, "resources", "SMILES", file_name),
             sep='\t', header=0)
            smiles_dfs.append(df)

        smiles_df = pd.concat(smiles_dfs)
        smiles_df = smiles_df.loc[~smiles_df['SMILES'].isnull(), :].copy()

        the_identifiers_that_entered, these_smiles_entered = compounds_and_smiles_ordered(
            smiles_df, molecules_ids_list)

        fpgen = AllChem.GetRDKitFPGenerator()
        the_ms = [Chem.MolFromSmiles(k) for k in these_smiles_entered]
        fps = [fpgen.GetFingerprint(x) for x in the_ms]

        arr = np.zeros(shape=(len(the_identifiers_that_entered),
                              len(the_identifiers_that_entered)))

        logger.info("filling array with Tanimoto similarities, computed for all pairs of SMILES")
        for u in range(len(the_identifiers_that_entered)):
            for v in range(u + 1, len(the_identifiers_that_entered)):
                tanimoto_similarity = DataStructs.TanimotoSimilarity(
                    fps[u],  fps[v])
                arr[u, v] = tanimoto_similarity
                arr[v, u] = tanimoto_similarity

        np.fill_diagonal(arr, 1)
        if threshold_similarity is not None:
            arr[arr < threshold_similarity] = 0  # todo: verify this with unittest

        mol_similarities_result = {"array": arr,
                                   "the_ids": the_identifiers_that_entered,
                                   "the_smiles": these_smiles_entered}

        return mol_similarities_result


class ComputeMoleculesFamilies:
    config = ComputeMoleculesFamiliesConfig
    @staticmethod
    def run(mol_similarities_result, molecules_ids_df,
            mydata_dir, subfolder, SAMPLE,
            cutoff_tree=15, method_linkage="ward", metric_linkage="cosine") -> dict:
        print("now doing the hierarchical clustering of the computed matrix")
        # make the hierarchical clustering reproducible
        np.random.seed(42)

        distance_matrix = 1 - mol_similarities_result['array']
        # squareform saves memory by representing the matrix with the upper right triangle
        computed_dendrogram = linkage(squareform(distance_matrix),
                                      method=method_linkage,
                                      metric=metric_linkage
                                      )
        clusters = fcluster(computed_dendrogram,
                            t=cutoff_tree, criterion='distance')
        myannots = molecules_ids_df.loc[molecules_ids_df["moleculeIds"].isin(
            mol_similarities_result['the_ids']),:]
        myannots = myannots.assign(molecules_cluster=clusters)

        logger.info("the clusters are: %s",
              sorted(myannots['molecules_cluster'].unique().tolist()))

        molecs_families_dict = {}  # initialize and save as mmt
        for clu in sorted(myannots['molecules_cluster'].unique().tolist()):
            mols_in_clus = myannots.loc[myannots[
                                            'molecules_cluster'] == clu, "moleculeIds"].unique().tolist()
            first_molecule_group = myannots.loc[
                myannots['molecules_cluster'] == clu, "moleculeNames"].head(
                1).item()
            molecs_families_dict[
                f'tg_{clu}'] = [first_molecule_group] + mols_in_clus

        max_length = max(len(v) for v in molecs_families_dict.values())

        for key in molecs_families_dict:
            if len(molecs_families_dict[key]) < max_length:
                molecs_families_dict[key] += [np.nan] * (
                            max_length - len(molecs_families_dict[key]))

        molecs_families_mmt_df = pd.DataFrame(molecs_families_dict)
        molecs_families_mmt_df = molecs_families_mmt_df.T
        molecs_families_mmt_df.to_csv(os.path.join(
            mydata_dir, subfolder, f"{SAMPLE}-molecules_families.mmt"),
            sep='\t', header=False)
        logger.info("saved .mmt file")

        return {"computed_dendrogram" : computed_dendrogram,
                "molecules_families": molecs_families_dict,
                "cutoff_tree": cutoff_tree}



class FakeMethodConfig(MethodConfig):
    description : str = "dummy test for config and method verification"
    dummy_cutoff: float = 0.5
    dummy_scale: float = 1


# 2 the methods using above configs
class FakeMethod:  # no-inheritance avoids instantiation: more practical
    config: FakeMethodConfig
    @staticmethod
    def run(arr: np.ndarray, config=FakeMethodConfig()) -> np.ndarray:
        """
        example of how to declare the method
        """
        logger.info("Will run FakeMethod "
                            "with the following config: %s", config)

        return pathways_methods.test_function(arr, config.dummy_scale,
                                              config.dummy_cutoff)

class ClassifyMZProtocol:
    config: ClassifyMZProtocolConfig
    @staticmethod
    def run(adata_mz, me_annots, msm_cutoff, offSample_cutoff, fdr_cutoff,
            moran_I_cutoff, moran_padj_cutoff,
            SAMPLE, plot_folder, config=ClassifyMZProtocolConfig()):
        logger.info("Running %s", config.description)
        logger.info(f">> Part a): Select m/z features with parameters: "
          f"MSM >= %s, offSample <= %s, FDR <= %s", msm_cutoff, offSample_cutoff,fdr_cutoff)
        me_annots, mz_good_from_part_a = zmm.part_a__mz_quality_mark(
            me_annots, msm_cutoff, offSample_cutoff, fdr_cutoff)

        logger.info(f"A total of %s m/z pass MSM, offSample and FDR cutoffs.", len(set(mz_good_from_part_a)))
        text_notice = f"valid $m/z$\n n = {len(mz_good_from_part_a)}"
        zviz.plot_mz_given_quality_metrics(me_annots.mz_quality, SAMPLE,
                                      plot_folder, msm_cutoff,
                                      offSample_cutoff, text_notice)
        logger.info("saved plot (MSM, offSample, FDR) to %s", plot_folder)
        adata_mz.uns["history"]['features_n_after_mz_clean_1'] = len(
            mz_good_from_part_a)

        if config.n_perms_moran_calc == "auto":
            config.n_perms_moran_calc = 50
            if adata_mz.obs.shape[0] > 100000:
                config.n_perms_moran_calc = 3

        logger.info(f">> Part b): Classifying m/z features, by :\n"
              f"Moran's I value >= %s,  (BH) adjusted p-value <= %s", moran_I_cutoff, moran_padj_cutoff)
        logger.info("Moran's metrics are only computed for the m/z that passed the a) filters)")
        all_quality_df, mz_good_part_b, moran_mz_df = zmm.part_b__mz_quality_by_autocorr(
            adata_mz, mz_good_from_part_a, me_annots,
            moran_I_cutoff, moran_padj_cutoff, config.n_neighs_for_graph, config.n_perms_moran_calc)

        text_notice = f"valid $m/z$: n = {len(mz_good_part_b)}"
        zviz.plot__moran_info(moran_mz_df, plot_folder,
                         SAMPLE,
                         text_notice, moran_I_cutoff, moran_padj_cutoff)
        logger.info("saved plot (autocorrelation) to %s", plot_folder)

        adata_mz.uns["history"]['features_n_after_mz_clean_2'] = len(
            mz_good_part_b)

        adata_mz = zmm.adata_mz_insert_full_quality(adata_mz,
                                                    moran_mz_df,
                                                    all_quality_df)
        n_good_mz = adata_mz.var.loc[
                adata_mz.var['final_quality'] == "good", :].shape[0]
        logger.info("a total of %s are good quality m/z features (check column 'final_quality' of adata_mz.var)", n_good_mz)

        return adata_mz


class SyncPathwaysAndAnnots:
    #config: SyncPathwaysAndAnnotsConfig  TODO : use ?
    @staticmethod
    def run(pathwaysObj: MetabolicPathways,
            annots_df: pd.DataFrame,  # TODO make it to accept also TranscriptomeAnnotations and create such TranscriptomeAnnotations (should be easy)
            usable_id: str="",
            mz_col: str="mz",
            min_set_size: int=2,
            max_set_size: int=500,
            mode: str="mz_aware"  # or "default", but mz_aware is preferred
            )-> Tuple[MetabolicPathways, MetabolomeAnnotations]:
        if not usable_id == pathwaysObj.molecule_id_type:
            logger.error(
                "No same types of compound ids: %s and %s, respectively",
                pathwaysObj.molecule_id_type, usable_id)
            raise AssertionError("see the error in .log file. Aborted.")
        if max_set_size == 1:
            mode = "default"
        pathwaysObj.intersect_with_annotations(annots_df, usable_id, mz_col)  # this updates pathwaysObj only
        pathwaysObj.constrain(min_set_size, max_set_size, mode)

        annots_df = annots_df.loc[
            annots_df[usable_id].isin(
                pathwaysObj.piled[pathwaysObj.molecule_id_type].unique()), :]

        logger.info("Pathways and Annotations objects synchronized, with"
                    " min set size = %s", min_set_size)

        return pathwaysObj, annots_df


class ComputePathwaysScores:
    config: ComputePathwaysScoresConfig

    @staticmethod
    def run(pathways_dict: Dict[str, List[str]], adata: ad.AnnData,
            config=ComputePathwaysScoresConfig()
            ) -> pd.DataFrame:
        logger.info(" Will compute spot-wise Pathways Scores "
                    "with the following config: %s", config)
        scores_df = pathways_methods.compute_pathways_scores(
            pathways_dict, adata, config.method, config.min_entity,
            config.n_groups_of_spots, config.threshold_n_spots_is_high)

        return scores_df


if __name__ == "__main__":
    myarr = np.array([0.1, 0.5, 0.7, 0.9])
    # Call the static method
    updated_config = FakeMethodConfig(dummy_scale=3) # instantiate the config object, if any new config to apply
    result = FakeMethod.run(myarr, config=updated_config)
    print("Original array:", myarr)
    print("Processed array:", result)
    result = FakeMethod.run(myarr)  # using defaults of config
    print("Original array:", myarr)
    print("Processed array:", result)