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
import zonist.processing.pathways_methods as pathways_methods

logger = logging.getLogger(__name__)
logger = helpers.reset_log_config(logger)


class BaseModel(PydanticBaseModel):
    # class Config:  # deprecated since pydantic v2
    #     arbitrary_types_allowed = True
    model_config = ConfigDict(arbitrary_types_allowed=True)


# 1 methods configs
class MethodConfig(BaseModel):
    class_type: str = "method"

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


class SyncPathwaysAndAnnots:
    @staticmethod
    def run(pathwaysObj: MetabolicPathways, annotsObj: MetabolomeAnnotations,
            min_pathway_size=2,
            max_pathway_size=500
            )-> Tuple[MetabolicPathways, MetabolomeAnnotations]:

        pathwaysObj.intersect_with_annotations(annotsObj)  # this updates pathwaysObj
        pathwaysObj.constrain(min_pathway_size, max_pathway_size)

        # todo: check if this .df attribute update is enough (other attrs to update?) for annotations:
        annotsObj.df = annotsObj.df.loc[
                    annotsObj.df["KEGG"].isin(
                        pathwaysObj.piled["KEGG"].unique()), :]

        logger.info("syncronized sucessfuly")
        return pathwaysObj, annotsObj


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