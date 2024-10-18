
import os
import logging
import time
from pathlib import Path
from typing import List, Dict, Literal, Optional, Set, Union

import pandas as pd


from pydantic import BaseModel as PydanticBaseModel
from pydantic import Extra

import zonist.helpers as helpers

# from zonist.helpers import (approximate_matching,
#                             purge_approximate_matching)

from zonist.constants import (supported_mappings,
                              supported_annotations)


logger = logging.getLogger(__name__)
logger = helpers.reset_log_config(logger)
logger.info("\nTime: {}".format(time.strftime('%Y%m%d-%H.%M.%S')))


class BaseModel(PydanticBaseModel):
    class Config:
        arbitrary_types_allowed = True


class Dataset(BaseModel, extra=Extra.allow):
    name : str
    data_dir: str  # /m/n/o/prj/st/spatial_metabo/data    
    subfolder: str  # where either imzml+ibd files, or spatialmatrix  ! 


class DataAnnotations(Dataset):
    annotation_files : Dict[str, str]
    delimiter : str = ","  # the metaspace current delimiter is comma
    col_molecule_ids: str = "moleculeIds"  
    col_mz : str = "mz"
    col_formula : str = "formula"
    col_molecule_names: str = "moleculeNames"

    def load_annotation(self):
        logger.info("Loading annotations")

        self.check_expectations_1()

        dfs_dict = {}
        admitted_keys = list(set(supported_annotations).intersection(
            set(list(self.annotation_files.keys()))))
        for anno_key in admitted_keys:
            logger.info(f"  loading {anno_key}")
            try:
                df = pd.read_csv(
                    os.path.join(
                        self.data_dir, self.subfolder, 
                        self.annotation_files[anno_key]),
                    comment='#',  # override any comment (metaspace comments)
                    sep=self.delimiter,
                    header=0
                )
            except Exception as e:
                logger.info(f"Error (read_csv): \n {e}")

            self.check_expectations_2(df)

            df = df[[self.col_mz, self.col_formula, 
                     self.col_molecule_ids, self.col_molecule_names]]
            
            df['database'] = anno_key
            dfs_dict[anno_key] = df

            self.annotation_pretable = helpers.combine_annotations(
                dfs_dict
                )
            self.annotation_table = helpers.spreadout_annotation_df(
                self.annotation_pretable, self.col_mz, self.col_formula, 
                self.col_molecule_ids, self.col_molecule_names
                )
     
    def check_expectations_1(self):
        if len(set(list(self.annotation_files.keys())).intersection(
            set(supported_annotations))) == 0:
            logger.error(f"Error, none of the annotations is supported," 
            f" or the dictionary keys are misspelled. Only: " 
            f"{supported_annotations}) are supported")
        for anno_key in list(self.annotation_files.keys()):
            if anno_key not in supported_annotations:
                logger.warning(f"warning: '{anno_key}' annotation is not" 
                " yet supported and will be ignored")

    def check_expectations_2(self, df):
        expected_columns_by_df = [self.col_mz, self.col_formula, 
                                  self.col_molecule_ids, 
                                  self.col_molecule_names]

        if len(list(df.columns)) <= 1:
            logger.error(f"Error: insuficient columns in annotation dataframe,"
                         f" field delimiter by default is '{self.delimiter}'" 
                         f" (you can change it by using .delimiter attribute)")

        if len(set(expected_columns_by_df)) != len(set(
            list(df.columns)).intersection(set(expected_columns_by_df))):
            logger.error(f"one or more columns of in your annotations are missing, "
                         f"expected ones are exactly {expected_columns_by_df}")

           
class DataImzmlIbd(Dataset):
    data_imz_ibd_basename : str # one
    # TODO: complete with real data later


class DataSpatialMatrix(Dataset):
    matrix_filename: str  # one only
    col_mz : str = "mz"
    style_coordinates_names: str = "metaspace"
    style_coordinates_explanation: Dict[str, str] = {
        'metaspace' : "xR_yS where R and S are integers"}

    def preload(self):
        print("loading")
    
    def check_expectations(self):
        print("checking")


class ReferenceMappings(BaseModel, extra=Extra.allow):
    data_dir: str
    subfolder : str  
    delimiter : str = '\t'  # we provide tab delimited .tsv viles
    mapping_files: Dict[str, str] # list with basename.extension

    def load_mappings(self):
        logger.info("Loading reference files for mapping compounds identifiers")
        self.check_expectations_map1()

        dfs_dict = {}
        admitted_keys = list(set(supported_mappings).intersection(
            set(list(self.mapping_files.keys()))))
        for map_key in admitted_keys:
            logger.info(f"  loading {map_key}")
            try:
                df = pd.read_csv(
                    os.path.join(
                        self.data_dir, self.subfolder, 
                        self.mapping_files[map_key]),
                    sep=self.delimiter,
                    header=0
                )
            except Exception as e:
                logger.info(f"Error (read_csv): \n {e}")

            #  self.check_expectations_map2(df)  # TODO: define

            df = df[[self.col_mz, self.col_formula, 
                     self.col_molecule_ids, self.col_molecule_names]]
            
            df['database'] = map_key
            dfs_dict[map_key] = df
        
        self.mapping_df_dict = dfs_dict


    def check_expectations_map1(self):
        if len(set(list(self.mapping_files.keys())).intersection(
            set(supported_mappings))) == 0:
            logger.error(f"Error, none of the mappings is supported," 
            f" or the dictionary keys are misspelled. Only: " 
            f"{supported_mappings}) are supported")
        for k in list(self.mapping_files.keys()):
            if k not in supported_mappings:
                logger.warning(f"warning: '{k}' mapping is not" 
                " yet supported and will be ignored")


      

    



    

    
    

