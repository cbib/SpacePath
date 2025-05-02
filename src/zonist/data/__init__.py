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
from pydantic import ConfigDict
from pandas.api.types import is_numeric_dtype

import zonist.helpers as helpers
import zonist.processing.pathways_methods as pathways_methods

from zonist.constants import (supported_mappings, supported_annotations,
                              columns_ref_map_expected,
                              not_metabolic_reactions_pathways,
                              mz_decimal_places)

logger = logging.getLogger(__name__)
logger = helpers.reset_log_config(logger)
#logger.info("\nTime: {}".format(time.strftime('%Y%m%d-%H.%M.%S')))


class BaseModel(PydanticBaseModel):
    # class Config:  # deprecated since pydantic v2
    #     arbitrary_types_allowed = True
    model_config = ConfigDict(arbitrary_types_allowed=True)


class Dataset(BaseModel, extra='allow'):
    name : str
    data_dir: str  # /main/dir/of/spatial/data
    subfolder: str    # e.g. projectName/


class MetaboliteIdsMappings(BaseModel, extra='allow'):
    data_dir: str
    subfolder: str
    delimiter: str = '\t'  # we provide tab delimited .tsv viles
    mapping_files: Dict[str, str]  # list with basename.extension
    table: pd.DataFrame = pd.DataFrame()
    style_desired: str = "KEGG"
    #

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
            df = df.loc[df[self.style_desired] != "", :]
            df = df.loc[~df[self.style_desired].isin(["NA", "NaN"]), :]
            if map_key == "LIPID_MAPS": # old version of lipidmaps, assimilate
                map_key = "LipidMaps"

            df['database'] = map_key
            df.columns = ["initial_id_unconverted", self.style_desired, "database"]
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
        if not np.all(np.array(cols) == np.array([map_key, self.style_desired])):
            logger.error("Error: the columns of this table are not "
                         "as expected, they must be exactly : "
                         "'%s' and '%s'", map_key, self.style_desired)


class MetabolomeAnnotations(Dataset):
    annotation_files : Dict[str, str]
    delimiter : str = ","  # the metaspace delimiter is comma
    molecule_ids_col: str = "moleculeIds"
    mz_col : str = "mz"
    formula_col : str = "formula"
    molecule_names_col: str = "moleculeNames"
    offsample_column: str = "rawOffSampleProb"
    total_intens_col: str = "totalIntensity"
    fdr_column: str = "fdr"
    msm_column: str = "msm"
    no_translated: pd.DataFrame = pd.DataFrame()
    usable_id: str = "KEGG"  # or moleculeIds if not using conversion to kegg
    df: pd.DataFrame = pd.DataFrame()
    mz_quality: pd.DataFrame = pd.DataFrame()

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
                df = df[[self.mz_col, self.formula_col, self.molecule_ids_col,
                         self.molecule_names_col, self.offsample_column,
                         self.total_intens_col, self.fdr_column, self.msm_column]]
            except KeyError:
                df = df[[self.mz_col, self.formula_col, self.molecule_ids_col,
                         self.molecule_names_col]]

            df[self.mz_col] = np.around(df[self.mz_col].to_numpy(),
                                        decimals=mz_decimal_places)
            df = df.dropna(how="all")
            if anno_key == "custom_kegg":
                anno_key = "KEGG"
            elif anno_key == "LIPID_MAPS":  # assimilate
                anno_key = "LipidMaps"
            df['database'] = anno_key
            dfs_dict[anno_key] = df

        self.preloaded = \
                helpers.combine_frames_identical_columns(dfs_dict)

        self.no_translated = helpers.spreadout_annotation_df(
                self.preloaded, self.mz_col, self.formula_col,
                self.molecule_ids_col, self.molecule_names_col)

        self.no_translated[self.mz_col] = np.around(
            self.no_translated[self.mz_col].to_numpy(), decimals=mz_decimal_places)

    def pull_mz_quality_df(self):
        """Builds df, if quality columns exist"""
        try:
            self.mz_quality = helpers.pick_mz_quality_infos(
                self.preloaded, self.mz_col, self.offsample_column,
                self.total_intens_col, self.fdr_column, self.msm_column,
                self.formula_col
            )
            if self.mz_quality[self.offsample_column].isnull().all():
                logger.warning("Useless: %s is NaN", self.offsample_column)
            else:
                if not self.mz_quality[self.offsample_column].max() <= 1:
                    raise ValueError("Values of %s must be comprised "
                                     "between 0 and 1", self.offsample_column)
        except KeyError:
            logger.error("Failed, one or more of required columns (%s) are missing", [
                self.offsample_column, self.total_intens_col, self.fdr_column, self.msm_column
                ])

    def discard_out_of_tissue_mz(self, prob_out=0.87,
                                 fraction_to_exclude=1):
        # todo: check if not needed anymore this discard_out_of_tissue_mz: I am using a more complex strategy: see repo user-pipeline/scripts/classify_normalize_features

        """Recommended (optional), to discard the mz that are out-of-tissue.
        Among the off-sample features (those having >= prob_out)
        ordered by decreasing total intensity, exclude the defined fraction"""
        self.pull_mz_quality_df()
        # if column is not empty NaN entirely, skip:
        if self.mz_quality[self.offsample_column].isnull().all():
            logger.warning("Useless: %s is NaN", self.offsample_column)
        else:
            if not self.mz_quality[self.offsample_column].max() <= 1:
                raise AssertionError("% must be between 0 and 1", self.offsample_column)

            bad_feat_df = self.mz_quality.loc[self.mz_quality[
                self.offsample_column] >= prob_out, :]

            bad_feat_df = bad_feat_df.sort_values(by=[
                self.mz_col, self.total_intens_col], ascending=[True, False])
            n_bad = int(bad_feat_df.shape[0] * fraction_to_exclude)
            mz_out = bad_feat_df[self.mz_col].tolist()[:n_bad]

            self.no_translated = self.no_translated.loc[~self.no_translated[
                self.mz_col].isin(mz_out)]
            logger.info('[discard_out_of_tissue_mz]: Excluded %s m/z '
                        'features that are "out-of-the-tissue"', len(mz_out))

    def keep_unique_mz_id_tuples(self, var_df: pd.DataFrame):
        logger.info(
            'Resolving multiple m/z to one-compound-id (degenerate pairs):\n'
            '    Keeps the best (m/z, id) unique pairs, by quality criteria.')
        basic_expected_cols = [self.mz_col, self.fdr_column,
                               self.offsample_column, self.msm_column]
        if not len(set(basic_expected_cols).intersection(
                set(var_df.columns))) == len(basic_expected_cols):
            logger.error("Error, %s are expected in var_df",
                         basic_expected_cols)
            raise AssertionError
        if "I" in var_df.columns:
            logger.info("    Moran's 'I' in var_df will be also used.")
        assert is_numeric_dtype(var_df[
                                    self.mz_col]), "the 'mz' column in adata.var must be of type float"
        assert len(set(var_df[self.mz_col].tolist()).intersection(
            set(self.df[self.mz_col].unique()))
        ) == len(self.df[self.mz_col].unique()), "Error, no perfect match with var_df mz"

        self.df = helpers.warrant__one_to_one_mz_id(self.df, var_df,
            self.mz_col, self.usable_id, self.msm_column,
            self.fdr_column, self.offsample_column)

    def convert_ids_to_kegg(self, ref_table: pd.DataFrame):
        logger.info('Translating ids in the metabolome annotation'
                    ' table, to KEGG (into .df attribute)')
        if not self.usable_id == "KEGG":
            logger.error('The .usable_id must be "KEGG". Aborted.')
            raise ValueError
        assert np.all(np.array(list(ref_table.columns)) == np.array(
                columns_ref_map_expected)), f"Metabolite reference"
        f" mappings columns must be {columns_ref_map_expected}"
        try:
            self.df = helpers.convert_metabo_annots_tokegg(
                self.no_translated, ref_table,
                self.molecule_ids_col, self.mz_col, self.usable_id)
            self.verify()
        except Exception as e:
            logger.warning("Error: %s", e)

    def verify(self):
        if self.usable_id != "KEGG":
            self.df = self.no_translated
        try:
            excess_mz_bool: bool = helpers.identify_multi_mz_to_one_compound(
                self.df, self.mz_col, self.usable_id)
            if excess_mz_bool:
                # logger.warning("Please run keep_unique_mz_id_tuples(anndata_mz)")
                # print(message_multimz)
                message_multimz = (
                    "Identified 'multi m/z to one KEGGid' cases!: Please run "
                    "'<metabo_annotations>.keep_unique_mz_id_tuples(anndata_mz)'.")
                logger.warning(message_multimz)
                print(message_multimz)
        except Exception as e:
            logger.warning("Error: %s", e)

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
            self.mz_col, self.formula_col,
           self.molecule_ids_col, self.molecule_names_col]

        if len(list(df.columns)) <= 1:
            logger.error("Error: insuficient columns in annotation dataframe,"
                         " field delimiter by default is '%s' (you can change" 
                         " it by using .delimiter attribute)", self.delimiter)

        if len(set(expected_columns_by_df)) != len(set(
            list(df.columns)).intersection(set(expected_columns_by_df))):
            logger.error("one or more columns of in your annotations are missing, "
                         "expected ones are exactly %s", expected_columns_by_df)


class MetabolicPathways(BaseModel, extra='allow'):
    data_dir: str
    subfolder: str
    official_file: str
    custom_file: Optional[str] = ''
    download: Optional[str | None] = None  # or "sspa-kegg" if download is desired
    undesired_substring: str = " - Homo sapiens (human)" # in names
    pathway_id_col: str = "pathway_id"
    pathway_name_col: str = "pathway_name"
    molecule_id_type: str = "KEGG"  # or moleculeIds (identical to annots.usable_id)
    piled: pd.DataFrame = pd.DataFrame()
    pathways_ids: pd.DataFrame = pd.DataFrame()
    purge: bool = True

    def preload_official(self):
        if not os.path.exists(os.path.join(
                self.data_dir, self.subfolder, self.official_file)):
            if self.download == "sspa-kegg":
                logger.info("writing sspa official gmt to data dir, wait ...")
                sspa.process_pathways.process_kegg(
                    organism="hsa", infile=None, download_latest=True,
                    filepath=os.path.join(self.data_dir, self.subfolder))
                logger.info("saved!")
            else:
                logger.info("only 'sspa-kegg' download is currently supported")
        try:
            official_gmt = pd.read_csv(
                os.path.join(self.data_dir, self.subfolder,
                             self.official_file),
                header=None, sep='\t', index_col=False)
        except Exception as e:
            logger.error(e)

        official_gmt = helpers.prettify_gmt(official_gmt,
                                            self.undesired_substring)
        if self.purge:
            official_gmt = helpers.purge_pathways(
                official_gmt, not_metabolic_reactions_pathways)
        self.official = official_gmt

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

        id_name_pathways_df, pws_piled = helpers.gmts_merge_pile_splitid(
            dfs_dict, self.pathway_id_col, self.pathway_name_col,
            self.molecule_id_type)

        self.pathways_ids = id_name_pathways_df
        self.piled = pws_piled

        logger.info("success, attributes : piled, pathways_ids")

    def get_name(self, a_pathway_id):
        pathway_name = self.pathways_ids.loc[
            self.pathways_ids[self.pathway_id_col] == a_pathway_id,
            self.pathway_name_col].item()
        return pathway_name

    def to_dict(self, option="id") -> Dict[str, List[str]]:
        if option=="id":
            mz_or_id = self.molecule_id_type
        elif option=="mz":
            mz_or_id = "mz"
        else:
            logger.error("unrecognized option, must be 'id' or 'mz'")
        return helpers.pws_piled_to_dict(
            self.piled, self.pathway_id_col, mz_or_identity_col=mz_or_id)

    def intersect_with_annotations(self, annots_obj: MetabolomeAnnotations):
        """
         - adds the mz column to the .piled df in the MetabolicPathways obj
         - retains only the elements that are in the intersection of both objs
        """
        if not self.molecule_id_type == annots_obj.usable_id:
            logger.error("%s and %s are not identical. Aborted.", self.molecule_id_type, annots_obj.usable_id)
            raise AssertionError
        self.piled = helpers.pathways_and_annots_join(
            self.piled[[self.pathway_id_col, self.molecule_id_type]],
            annots_obj.df[[annots_obj.mz_col, annots_obj.usable_id]],
            self.molecule_id_type)

    def constrain(self, min_pathway_size, max_pathway_size, mode):
        assert mode in ["default", "mz_aware"], "Error, mode param not recognized"
        if mode == "default" :
            good_pathways_l = pathways_methods.yield_compliant_pathways(
                self.to_dict(option="id"), min_pathway_size, max_pathway_size)
        elif mode == "mz_aware":
            good_pathways_l = pathways_methods.yield_compliant_pathways(
                self.to_dict(option="mz"), min_pathway_size, max_pathway_size)

        metabolites_to_eject = pathways_methods.find_metabolites_to_exclude(
            good_pathways_l, self.to_dict())
        self.piled = self.piled.loc[
            self.piled[self.pathway_id_col].isin(good_pathways_l)]
        self.piled = self.piled.loc[~self.piled[
            self.molecule_id_type].isin(metabolites_to_eject), :]
        self.pathways_ids = self.pathways_ids.loc[
                                   self.pathways_ids[
                                       self.pathway_id_col].isin(good_pathways_l), :]
