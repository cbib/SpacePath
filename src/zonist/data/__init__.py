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

import zonist.helpers as helpers
import zonist.processing.pathways_methods as pathways_methods

from zonist.constants import (supported_mappings, supported_annotations,
                              columns_ref_map_expected,
                              not_metabolic_reactions_pathways)

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
    fdr_column: str= "fdr"
    unfiltered: pd.DataFrame = pd.DataFrame()
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
                         self.total_intens_col, self.fdr_column]]
            except KeyError:
                df = df[[self.mz_col, self.formula_col, self.molecule_ids_col,
                         self.molecule_names_col, self.fdr_column]]
                auto_remove_offsample = False  # not enough columns to be True

            df[self.mz_col] = np.around(df[self.mz_col].to_numpy(), decimals=9)
            df = df.dropna(how="all")
            df['database'] = anno_key
            if anno_key == "custom_kegg":
                df['database'] = "KEGG"  # ok
            #elif anno_key == "CoreMetabolome":
            #    df['database'] = "HMDB"  # ok
            dfs_dict[anno_key] = df

        preloaded_annots = \
                helpers.combine_frames_identical_columns(dfs_dict)

        self.preloaded = preloaded_annots

        self.unfiltered = helpers.spreadout_annotation_df(
                preloaded_annots, self.mz_col, self.formula_col,
                self.molecule_ids_col, self.molecule_names_col)

        self.unfiltered[self.mz_col] = np.around(
            self.unfiltered[self.mz_col].to_numpy(), decimals=9)

    def do_mz_quality_df(self):
        """Builds df, if probability of being outside-of-tissue exists"""
        if self.offsample_column in self.preloaded.columns.tolist():
            self.mz_quality = helpers.pull_mz_quality_info(
                self.preloaded, self.mz_col, self.offsample_column,
                self.total_intens_col, self.fdr_column)
            if self.mz_quality[self.offsample_column].isnull().all():
                logger.warning("Useless: %s is NaN", self.offsample_column)
            else:
                if not self.mz_quality[self.offsample_column].max() <= 1:
                    raise ValueError("Values of %s must be comprised "
                                     "between 0 and 1", self.offsample_column)
        else:
            logger.error("Failed, no %s column", self.offsample_column)

    def discard_out_of_tissue_mz(self, prob_out=0.87,
                                 fraction_to_exclude=1):
        """Recommended (optional), to discard the mz that are out-of-tissue.
        Among the off-sample features (those having >= prob_out)
        ordered by decreasing total intensity, exclude the defined fraction"""
        self.do_mz_quality_df()
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

            self.unfiltered = self.unfiltered.loc[~self.unfiltered[
                self.mz_col].isin(mz_out)]
            logger.info('[discard_out_of_tissue_mz]: Excluded %s m/z '
                        'features that are "out-of-the-tissue"', len(mz_out))

    def convert_ids_to_kegg(self, ref_table: pd.DataFrame):   # todo: is it possible to make it generalizable and not only KEGG dependent ?
        logger.info('Translating ids in the metabolome annotation'
                    ' table, to KEGG')
        assert np.all(np.array(list(ref_table.columns)) == np.array(
                columns_ref_map_expected)), f"Metabolite reference"
        f" mappings columns must be {columns_ref_map_expected}"
        try:
            self.df = helpers.convert_metabo_annots_tokegg(
                self.unfiltered, ref_table,
                self.molecule_ids_col, self.mz_col)
            excess_mz_bool: bool = helpers.identify_multi_mz_to_one_kegg(
                self.df,  self.mz_col)
            if excess_mz_bool:
                message_multimz = (
                    "Identified 'multi m/z to one KEGGid' cases!: Please run "
                    "'<metabo_annotations>.treat_mz_multiplets(anndata_mz)'.")
                logger.warning(message_multimz)
                print(message_multimz)
        except Exception as e:
            logger.warning("Error: %s", e)

    def treat_mz_multiplets(self, adata_mz: ad.AnnData):
        ## TODO : add a check to see if the adata_mz has the same mz values as the ones in the annotation
        logger.info('[treat_mz_multiplets]: if "multi m/z to one KEGGid", for'
                    ' such KEGGid, the m/z with the top variance is kept.')
        self.df = helpers.compute_one_to_one_mz_keggid_df(
            self.df, adata_mz, self.mz_col)

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
    file_style: str = "sspa-kegg"  # for both official and custom pathways
    undesired_substring: str = " - Homo sapiens (human)" # in names
    pathway_id_col: str = "pathway_id"
    pathway_name_col: str = "pathway_name"
    molec_ids_in_pathways_col: str = "KEGG"   #  must match joined annots
    piled: pd.DataFrame = pd.DataFrame()
    pathways_ids: pd.DataFrame = pd.DataFrame()

    def preload_official(self):
        if not os.path.exists(os.path.join(
                self.data_dir, self.subfolder, self.official_file)):
            if self.file_style == "sspa-kegg":
                logger.info("writing sspa official gmt to data dir, wait ...")
                sspa.process_pathways.process_kegg(
                    organism="hsa", infile=None, download_latest=True,
                    filepath=os.path.join(self.data_dir, self.subfolder))
                logger.info("saved!")
            else:
                logger.info("only 'sspa-kegg' file_style is currently supported")
        try:
            official_gmt = pd.read_csv(
                os.path.join(self.data_dir, self.subfolder,
                             self.official_file),
                header=None, sep='\t', index_col=False)
        except Exception as e:
            logger.error(e)

        official_gmt = helpers.prettify_gmt(official_gmt,
                                            self.undesired_substring)
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
            self.molec_ids_in_pathways_col)

        self.pathways_ids = id_name_pathways_df
        self.piled = pws_piled

        logger.info("success, attributes : piled, pathways_ids")

    def to_dict(self) -> Dict[str, List[str]]:
        return helpers.pws_piled_to_dict(
            self.piled, self.pathway_id_col,
            mz_or_identity_col=self.molec_ids_in_pathways_col)

    def intersect_with_annotations(self, annots_obj: MetabolomeAnnotations):
        """
         - adds the mz column to the .piled df in the MetabolicPathways obj
         - retains only the elements that are in the intersection of both objs
        """
        self.piled = helpers.pathways_and_annots_join(
            self.piled,
            annots_obj.df[[annots_obj.mz_col, self.molec_ids_in_pathways_col]],
            self.molec_ids_in_pathways_col)

    def constrain(self, min_pathway_size, max_pathway_size):
        good_pathways_l = pathways_methods.yield_compliant_pathways(
            self.to_dict(), min_pathway_size, max_pathway_size)
        metabolites_to_eject = pathways_methods.find_metabolites_to_exclude(
            good_pathways_l, self.to_dict())
        self.piled = self.piled.loc[
            self.piled["pathway_id"].isin(good_pathways_l)]
        self.piled = self.piled.loc[~self.piled[
            self.molec_ids_in_pathways_col].isin(metabolites_to_eject), :]
        self.pathways_ids = self.pathways_ids.loc[
                                   self.pathways_ids[
                                       "pathway_id"].isin(good_pathways_l), :]
