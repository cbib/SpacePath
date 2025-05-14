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


import zonist.helpers as helpers
from zonist.data import MetabolicPathways, MetabolomeAnnotations
import zonist.processing.measures_methods as zmm
import zonist.processing.pathways_methods as pathways_methods
import zonist.visualization.visualization as zviz


logger = logging.getLogger(__name__)
logger = helpers.reset_log_config(logger)


class BaseModel(PydanticBaseModel):
    # class Config:  # deprecated since pydantic v2
    #     arbitrary_types_allowed = True
    model_config = ConfigDict(arbitrary_types_allowed=True)


# 1 methods configs
class MethodConfig(BaseModel):
    class_type: str = "method"


class ClassifyMZProtocolConfig(MethodConfig):
    description: str = "phased protocol to filter m/z by the quality of molecular formulae annotations, and the spatial distributions"
    n_neighs: int = 20   # for neighbors graph computation before Moran's I
    n_perms: Union[int, str] = "auto"   # for Moran's I


class SyncPathwaysAndAnnotsConfig(MethodConfig): # todo: use this config ?
    description: str = "parameters to synchronize pathways and annotations"


class ComputePathwaysScoresConfig(MethodConfig):
    description: str = "parameters to compute the pathways scores"
    method: str = "sspa_kpca"
    min_entity: int = 2
    n_groups_of_spots: int = 30
    threshold_n_spots_is_high: int = 15000   # todo : in the stable version use 20K


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
          f"MSM >= %s, offSample <= %s, FDR <= %s", msm_cutoff, offSample_cutoff,fdr_cutoff )
        me_annots, mz_to_keep_part1 = zmm.part_a_filter(
            me_annots, msm_cutoff, offSample_cutoff, fdr_cutoff)

        logger.info(f"A total of %s m/z pass MSM, offSample and FDR cutoffs.", len(set(mz_to_keep_part1)))
        text_notice = f"valid $m/z$\n n = {len(mz_to_keep_part1)}"
        zviz.plot_mz_given_quality_metrics(me_annots.mz_quality, SAMPLE,
                                      plot_folder, msm_cutoff,
                                      offSample_cutoff, text_notice)
        logger.info("saved plot (MSM, offSample, FDR) to %s", plot_folder)
        adata_mz.uns["history"]['features_n_after_mz_clean_1'] = len(
            mz_to_keep_part1)

        if config.n_perms == "auto":
            config.n_perms = 50
            if adata_mz.obs.shape[0] > 100000:
                config.n_perms = 3

        logger.info(f">> Part b): Classifying m/z features, by :\n"
              f"Moran's I value >= %s,  (BH) adjusted p-value <= %s", moran_I_cutoff, moran_padj_cutoff)
        logger.info("Moran's metrics are only computed for the m/z that passed the a) filters)")
        all_quality_df, mz_to_keep_part2, moran_mz_df = zmm.part_b_filter_with_autocorr(
            adata_mz, mz_to_keep_part1, me_annots,
            moran_I_cutoff, moran_padj_cutoff, config.n_neighs, config.n_perms)

        text_notice = f"valid $m/z$: n = {len(mz_to_keep_part2)}"
        zviz.plot__moran_info(moran_mz_df, plot_folder,
                         SAMPLE,
                         text_notice, moran_I_cutoff, moran_padj_cutoff)
        logger.info("saved plot (autocorrelation) to %s", plot_folder)

        adata_mz.uns["history"]['features_n_after_mz_clean_2'] = len(
            mz_to_keep_part2)

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
    def run(pathwaysObj: MetabolicPathways, annotsObj: MetabolomeAnnotations,  # TODO make it to accept also TranscriptomeAnnotations and create such TranscriptomeAnnotations (should be easy)
            min_pathway_size: int=2,
            max_pathway_size: int=500,
            mode: str="default"  # or mz_aware  # todo?: make automatic? "mz_aware" ??if an object of type MetabolomeAnnotations is the second argument  ??
            )-> Tuple[MetabolicPathways, MetabolomeAnnotations]:
        if not (annotsObj.usable_id == pathwaysObj.molecule_id_type):
            logger.error(
                "No same types of compound ids: %s and %s, respectively",
                pathwaysObj.molecule_id_type, annotsObj.usable_id)
            raise AssertionError("see the error in .log file. Aborted.")
        if (max_pathway_size == 1) and (mode != "default"):
            logger.info("max pathway size is 1, auto-set default mode")
            mode = "default"
        pathwaysObj.intersect_with_annotations(annotsObj)  # this updates pathwaysObj
        pathwaysObj.constrain(min_pathway_size, max_pathway_size, mode)

        # todo: check if this .df attribute update is enough (other attrs to update?) for annotations:
        annotsObj.df = annotsObj.df.loc[
            annotsObj.df[annotsObj.usable_id].isin(
                        pathwaysObj.piled[pathwaysObj.molecule_id_type].unique()), :]

        logger.info("Pathways and Annotations objects syncronized with success")
        return pathwaysObj, annotsObj


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