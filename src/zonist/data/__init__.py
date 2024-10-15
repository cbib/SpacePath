import logging
import os
from pathlib import Path
from typing import List, Dict, Literal, Optional, Set, Union

import pandas as pd


from pydantic import BaseModel as PydanticBaseModel

from zonist.helpers import (approximate_matching,
                            purge_approximate_matching)

class BaseModel(PydanticBaseModel):
    class Config:
        arbitrary_types_allowed = True


class Dataset(BaseModel):
    #config: DatasetConfig
    name : str
    parent_dir: str  # until the level of data dir, e.g. : /m/n/o/prj/st/spatial_metabo/data    
    subfolder_data: str  #  contains either imzml+ibd files, or spatialmatrices tsv files ! 

    # their annotations alongside as well
    annotations_filenames : Dict[str, str]  # 1 or more

class References(BaseModel):
    parent_dir: str
    subfolder_mappings : str  # containts the files of identifiers mappings to kegg
    mappings_filenames: Dict[str, str] # list with full namees if any files there, with their extension !

    
class DataImzmlIbd(Dataset):
    data_imz_ibd_basename : str # one


class DataSpatialMatrices(Dataset):
    matrices_filenames: Dict[str, str]  # one or more

  

    



    

    
    

