# johaGL 2024
import os
import logging
import numpy as np
import time
from pathlib import Path
from typing import List, Dict, Literal, Optional, Set, Union

import pandas as pd
import anndata as ad
import sspa as sspa

from pydantic import BaseModel as PydanticBaseModel
from pydantic import Extra

import zonist.helpers as helpers

from zonist.constants import (supported_mappings,
                              supported_annotations,
                              columns_ref_map_expected)


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



class MetaboliteIdsMappings(BaseModel, extra=Extra.allow):
    data_dir: str
    subfolder: str
    delimiter: str = '\t'  # we provide tab delimited .tsv viles
    mapping_files: Dict[str, str]  # list with basename.extension
    table: pd.DataFrame = pd.DataFrame()

    def load(self):
        logger.info(
            "Loading reference files for mapping compounds identifiers")
        self.check_expectations_map1()

        dfs_dict = {}
        admitted_keys = list(set(supported_mappings).intersection(
            set(list(self.mapping_files.keys()))))
        for map_key in admitted_keys:
            logger.info("  loading %s", map_key)
            try:
                df = pd.read_csv(os.path.join(
                        self.data_dir, self.subfolder,
                        self.mapping_files[map_key]),
                    sep=self.delimiter, header=0)
            except Exception as e:
                logger.info("Error (read_csv): \n %s", e)

            self.check_expectations_map2(df, map_key)
            #logger.info(" your loaded dataframe :\n {}".format(df.head(2)))
            logger.info(" your loaded dataframe :\n %s", df.head(2))
            df = df.dropna()
            df = df.loc[df["KEGG"] != "", :]
            df = df.loc[df["KEGG"] != "NA", :]
            df = df.loc[df["KEGG"] != "NaN", :]
            df['database'] = map_key
            df.columns = ["not_kegg_id", "KEGG", "database"]
            dfs_dict[map_key] = df

        self.table = helpers.combine_frames_identical_columns(dfs_dict)

    def check_expectations_map1(self):
        if len(set(list(self.mapping_files.keys())).intersection(
                set(supported_mappings))) == 0:
            logger.error("Error: none of the mappings is supported, "
                         "or the dictionary keys are misspelled. Only: "
                         "%s are supported", supported_mappings)
        for k in list(self.mapping_files.keys()):
            if k not in supported_mappings:
                logger.warning("warning: '%s' mapping is not"
                               " yet supported and will be ignored", k)

    def check_expectations_map2(self, df, map_key: str):
        cols = list(df.columns)
        if not np.all(np.array(cols) == np.array([map_key, "KEGG"])):
            logger.error("Error: the columns of this table are not "
                         "as expected, they must be exactly : "
                         "'%s' and 'KEGG'", map_key)


class MetabolomeAnnotations(Dataset):
    annotation_files : Dict[str, str]
    delimiter : str = ","  # the metaspace delimiter is comma
    col_molecule_ids: str = "moleculeIds"
    col_mz : str = "mz"
    col_formula : str = "formula"
    col_molecule_names: str = "moleculeNames"
    offsample_column: str = "offSample"  # has bool values in feat_classif df
    unfiltered: pd.DataFrame = pd.DataFrame()
    df: pd.DataFrame = pd.DataFrame()
    feat_classif: pd.DataFrame = pd.DataFrame()

    def load(self):
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
                    comment='#',  # override any comment (metaspace does)
                    sep=self.delimiter, header=0)
            except Exception as e:
                logger.info("Error (read_csv): \n %s", e)

            self.check_expectations_2(df)

            try:
                df = df[[self.col_mz, self.col_formula, self.col_molecule_ids,
                         self.col_molecule_names, self.offsample_column]]
            except KeyError:
                df = df[[self.col_mz, self.col_formula, self.col_molecule_ids,
                         self.col_molecule_names]]
            df[self.col_mz] = np.around(df[self.col_mz].to_numpy(), decimals=9)
            df = df.dropna(how="all")
            df['database'] = anno_key
            if anno_key == "custom_kegg":
                df['database'] = "KEGG"
            dfs_dict[anno_key] = df

        preloaded_annots = \
                helpers.combine_frames_identical_columns(dfs_dict)

        self.preloaded = preloaded_annots

        self.unfiltered = helpers.spreadout_annotation_df(
                preloaded_annots, self.col_mz, self.col_formula,
                self.col_molecule_ids, self.col_molecule_names)

        self.unfiltered[self.col_mz] = np.around(
            self.unfiltered[self.col_mz].to_numpy(), decimals=9)

    def run_feat_classif(self):
        """optional, if offsample column exists"""
        self.feat_classif = helpers.pull_offsample_info(self.preloaded,
                                                        self.col_mz,
                                                        self.offsample_column)

    def discard_out_of_sample_mz(self):
        """optional, to discard the mz that are out-of-tissue"""
         # Tested in 'user-pipeline/drafts': it is NOT possible TO FIND THE PIXELS THAT ARE RICH IN THOSE MZ-TO-BE-DISCARDED unfortunately
        mz_keep = self.feat_classif.loc[self.feat_classif[self.offsample_column] == False,
                    self.col_mz].tolist()
        updated_df = self.unfiltered.loc[
            self.unfiltered[self.col_mz].isin(mz_keep)]
        self.unfiltered = updated_df
        logger.info(' [discard_out_of_sample_mz]> Excluded %s unique '
                    'features that are "out-of-the-tissue". Retained: %s.',
                    len(set(self.feat_classif[
                                self.col_mz].unique()) - set(mz_keep)),
                    len(set(mz_keep)))

    def convert_ids_to_kegg(self, ref_table: pd.DataFrame):
        logger.info('Translating ids in the metabolome annotation'
                    ' table, to KEGG')
        assert np.all(np.array(list(ref_table.columns)) == np.array(
                columns_ref_map_expected)), f"Metabolite reference"
        f" mappings columns must be {columns_ref_map_expected}"
        try:
            self.df = helpers.convert_metabo_annots_tokegg(
                self.unfiltered, ref_table,
                self.col_molecule_ids, self.col_mz)
            excess_mz_bool: bool = helpers.identify_multi_mz_to_one_kegg(
                self.df,  self.col_mz)
            if excess_mz_bool:
                message_multimz = (
                    "Identified 'multi m/z to one KEGGid' cases!: Please run "
                    "'<obj>.treat_mz_multiplets(anndata_mz)'.")
                logger.warning(message_multimz)
                print(message_multimz)
        except Exception as e:
            logger.warning("Error: %s", e)

    def treat_mz_multiplets(self, adata_mz: ad.AnnData):
        ## TODO : add a check to see if the adata_mz has the same mz values as the ones in the annotation
        logger.info('[treat_mz_multiplets]: if "multi m/z to one KEGGid", for'
                    ' such KEGGid, the m/z with the top variance is kept.')
        self.df = helpers.compute_one_to_one_mz_keggid_df(
            self.df, adata_mz, self.col_mz)

    def check_expectations_1(self):
        if len(set(list(self.annotation_files.keys())).intersection(
            set(supported_annotations))) == 0:
            logger.error("Error, none of the annotations is supported,"
                " or the dictionary keys are misspelled. Only: " 
                "%s are supported", supported_annotations)
        for anno_key in list(self.annotation_files.keys()):
            if anno_key not in supported_annotations:
                logger.warning("warning: %s annotation is not" 
                " yet supported and will be ignored", anno_key)

    def check_expectations_2(self, df):
        expected_columns_by_df = [
            self.col_mz, self.col_formula,
           self.col_molecule_ids, self.col_molecule_names]

        if len(list(df.columns)) <= 1:
            logger.error("Error: insuficient columns in annotation dataframe,"
                         " field delimiter by default is '%s' (you can change" 
                         " it by using .delimiter attribute)", self.delimiter)

        if len(set(expected_columns_by_df)) != len(set(
            list(df.columns)).intersection(set(expected_columns_by_df))):
            logger.error("one or more columns of in your annotations are missing, "
                         "expected ones are exactly %s", expected_columns_by_df)


class MetabolicPathways(BaseModel, extra=Extra.allow):
    data_dir: str
    subfolder: str
    official_file: str
    custom_file: Optional[str] = ''
    file_style: str = "gmt-sspa"  #  official and custom must be same style
    undesired_name_substring: str = " - Homo sapiens (human)" # in names
    col_pathway_id: str = "pathway_id"
    col_pathway_name: str = "pathway_name"
    col_molecule_ids: str = "KEGG"   #  must be KEGG to match joined annots
    piled: pd.DataFrame = pd.DataFrame()
    pathways_ids: pd.DataFrame = pd.DataFrame()
    pathways: Dict[str, List[str]] = dict()

    def preload_official(self):
        if not os.path.exists(os.path.join(
                self.data_dir, self.subfolder, self.official_file)):
            if self.file_style == "gmt-sspa":
                logger.info("writing sspa official gmt to data dir, wait ...")
                sspa.process_pathways.process_kegg(
                    organism="hsa", infile=None, download_latest=True,
                    filepath=os.path.join(self.data_dir, self.subfolder))
                logger.info("saved!")
            else:
                logger.info("only 'gmt-sspa' file_style is currently supported")
        try:
            official_gmt = pd.read_csv(
                os.path.join(self.data_dir, self.subfolder,
                             self.official_file),
                header=None, sep='\t', index_col=False)
        except Exception as e:
            logger.error(e)

        self.official = helpers.prettify_gmt(
            official_gmt, self.undesired_name_substring)

    def preload_custom(self):
        try:
            custom_gmt = pd.read_csv(
                os.path.join(self.data_dir, self.subfolder,
                             self.custom_file),
                header=None, sep='\t', index_col=False)
        except Exception as e:
            logger.error(e)
        first_row = custom_gmt.iloc[0, :]
        cus_message = (f" The first line of your custom gmt file:\n"
                      f"  abbreviation: '{first_row.iloc[0]}'; "
                      f"  name : '{first_row.iloc[1]}';  metabolites"
                      f"  (first three) : '{first_row.iloc[2:5].tolist()}'.")
        logger.info("%s", cus_message)
        self.custom = helpers.prettify_gmt(custom_gmt, None)

    def load(self):
        self.preload_official()
        dfs_dict = {'official': self.official}
        if self.custom_file != '':
            self.preload_custom()
            dfs_dict['custom'] = self.custom

        id_name_pathways_df, gmt_piled = helpers.gmts_merge_pile_splitid(
            dfs_dict, self.col_pathway_id, self.col_pathway_name,
            self.col_molecule_ids)

        self.pathways_ids = id_name_pathways_df
        self.piled = gmt_piled

        logger.info("success, attributes : piled, pathways_ids")

    def intersect_with_annotations(self, annots_obj: MetabolomeAnnotations):  # pd.DataFrame
        """adds:
         - the mz column to the same .piled dataframe
         - the attribute .pathways, a dictionary with pathway_id as key"""
        annotations_df_mz_kegg =  annots_obj.df[[annots_obj.col_mz, "KEGG"]]
        self.piled = helpers.pathways_and_annots_join(
            self.piled, annotations_df_mz_kegg, self.col_molecule_ids)

        self.pathways = helpers.gmt_piled_to_dict(
            self.piled, self.col_pathway_id, mz_or_kegg=self.col_molecule_ids)

    # def write_to_json(self, dataset_subfolder: str, mz_or_kegg_choice: str = "mz"):
    #     full_path_file_json = os.path.join(
    #             self.data_dir, dataset_subfolder,
    #             f"gmt_{mz_or_kegg_choice}.json")
    #
    #     gmt_dict = helpers.gmt_dict_write_json(
    #         self.piled, self.col_pathway_id,
    #         mz_or_kegg_choice, full_path_file_json)
    #     self.gmt_dict = gmt_dict
    #     logger.info("saved json file to : %s", full_path_file_json)


class MetabolomeSpatialMatrix(Dataset):
    # TODO : is this class is useful? or use adata object directly instead ?
    adata_filename: str  # one only
    # col_mz : str = "mz"
    #style_coordinates_names: str = "metaspace"
    #style_coordinates_explanation: Dict[str, str] = {
    #    'metaspace' : "xR_yS where R and S are integers"}

    def preload(self):
        print("not yet supported")
        # TODO
        #print("loading")
    
    def check_expectations(self):
        print("not yet supported")
        # TODO
        # print("checking")





      

    



    

    
    

