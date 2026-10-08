# johaGL 2024
import logging
import os
from pathlib import Path
from typing import Dict, List, Optional, Union
import numpy as np
import pandas as pd
from pandas.api.types import is_numeric_dtype

from pydantic import BaseModel as PydanticBaseModel
from pydantic import ConfigDict

import SpacePath.core.helpers as helpers
import SpacePath.core.processing.pathways_methods as pathways_methods
from SpacePath.core.constants import (columns_kegg_conversion_expected,
                                      message_multimz,
                                      mz_decimal_places,
                                      not_metabolic_reactions_pathways,
                                      supported_annotations,
                                      supported_mappings)

import sspa as sspa

logger = logging.getLogger(__name__)
logger = helpers.setup_logger(logger)

module_path = Path(__file__).resolve()
library_path = module_path.parent.parent.parent  # .parent


class BaseModel(PydanticBaseModel):
    """
    New usage of pydantic base model
    Note that 'class Config:' is deprecated since pydantic v2.
    """
    model_config = ConfigDict(arbitrary_types_allowed=True)


class Dataset(BaseModel, extra='allow'):
    name: str
    data_dir: str  # data/ is recommended
    subfolder: str  # e.g. dataset-01/


class MetaboliteIdsMappings(BaseModel, extra='allow'):
    data_dir: str
    subfolder: str
    delimiter: str = '\t'  # we provide tab delimited .tsv files
    mapping_files: Dict[str, str]  # list with basename.extension
    table: pd.DataFrame = pd.DataFrame()
    style_desired: str = "KEGG"

    def load(self):
        logger.info(
            "Loading reference tables of %s identifiers",
            self.style_desired)
        self.check_expectations_map1()

        dfs_dict = {}
        admitted_keys = list(set(supported_mappings).intersection(
            set(list(self.mapping_files.keys()))))
        for map_key in admitted_keys:
            logger.info(" loading %s", map_key)
            file_path = os.path.join(self.data_dir, self.subfolder,
                                     self.mapping_files[map_key])
            try:
                df = pd.read_csv(file_path, sep=self.delimiter, header=0)
            except FileNotFoundError:
                logger.error("Mapping file not found: %s", file_path)
                raise
            except pd.errors.EmptyDataError:
                logger.error("Mapping file is empty: %s", file_path)
                raise
            except Exception as e:
                logger.error("Failed to read mapping file '%s': %s", map_key,
                             e)
                raise

            self.check_expectations_map2(df, map_key)
            df = df.dropna()
            df = df.loc[df[self.style_desired] != "", :]
            df = df.loc[~df[self.style_desired].isin(["NA", "NaN"]), :]
            if map_key == "LIPID_MAPS":  # old lipidmaps version, assimilate
                map_key = "LipidMaps"

            df = df.assign(database=map_key)
            df.columns = ["initial_id_unconverted", self.style_desired,
                          "database"]
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
        if not np.all(
                np.array(cols) == np.array([map_key, self.style_desired])):
            logger.error("Error: the columns of this table are not "
                         "as expected, they must be exactly : "
                         "'%s' and '%s'", map_key, self.style_desired)


class MetabolomeAnnotations(Dataset):
    annotation_files: Dict[str, str]
    delimiter: str = ","  # the metaspace delimiter is comma
    molecule_ids_col: str = "moleculeIds"
    mz_col: str = "mz"
    formula_col: str = "formula"
    molecule_names_col: str = "moleculeNames"
    offsample_column: str = "rawOffSampleProb"
    total_intens_col: str = "totalIntensity"
    fdr_column: str = "fdr"
    msm_column: str = "msm"
    no_translated: pd.DataFrame = pd.DataFrame()
    usable_id: str = "KEGG"  # or moleculeIds if not using conversion to kegg
    df: pd.DataFrame = pd.DataFrame()
    mz_quality: pd.DataFrame = pd.DataFrame()
    discard_plant_compounds: bool = True
    discard_pharmaceuticals: bool = True
    discard_other_non_animal: bool = True
    discard_odd_n_carb_lipids: bool = True
    cutoff_n_carb: int = 17  # used in combination with odd lipids filtering

    def load(self):
        logger.info("Loading annotations for dataset: '%s'", self.name)

        self.check_expectations_1()

        dfs_dict = {}
        admitted_keys = list(set(supported_annotations).intersection(
            set(list(self.annotation_files.keys()))))

        for anno_key in admitted_keys:
            file_path = os.path.join(self.data_dir, self.subfolder,
                                     self.annotation_files[anno_key])
            try:
                df = pd.read_csv(
                    file_path,
                    comment='#',  # override any comment (metaspace does)
                    sep=self.delimiter, header=0)
            except FileNotFoundError:
                logger.error("Annotation file not found: %s", file_path)
                raise
            except pd.errors.EmptyDataError:
                logger.error("Annotation file is empty: %s", file_path)
                raise
            except Exception as e:
                logger.error("Failed to read annotation file '%s': %s",
                             anno_key, e)
                raise

            self.check_expectations_2(df)
            try:
                df = df[[self.mz_col, self.formula_col, self.molecule_ids_col,
                         self.molecule_names_col, self.offsample_column,
                         self.total_intens_col, self.fdr_column,
                         self.msm_column]]
            except KeyError:
                df = df[[self.mz_col, self.formula_col, self.molecule_ids_col,
                         self.molecule_names_col]]

            df = df.assign(
                **{self.mz_col: np.around(df[self.mz_col].to_numpy(),
                                          decimals=mz_decimal_places)})
            df = df.dropna(how="all")
            if anno_key == "custom_kegg":
                anno_key = "KEGG"
            elif anno_key == "LIPID_MAPS":  # assimilate to LipidMaps
                anno_key = "LipidMaps"
            df = df.assign(database=anno_key)
            dfs_dict[anno_key] = df

        self.preloaded = \
            helpers.combine_frames_identical_columns(dfs_dict)

        self.no_translated = helpers.spreadout_annotation_df(
            self.preloaded, self.mz_col, self.formula_col,
            self.molecule_ids_col, self.molecule_names_col)

        self.no_translated[self.mz_col] = np.around(
            self.no_translated[self.mz_col].to_numpy(),
            decimals=mz_decimal_places)

        self.no_translated = helpers.purge_annotations(
            self.no_translated, self.discard_pharmaceuticals,
            self.discard_plant_compounds, self.discard_other_non_animal,
            self.discard_odd_n_carb_lipids, self.cutoff_n_carb,
            self.molecule_ids_col, self.molecule_names_col,
            library_path=library_path)
        logger.info(" success, annotations loaded (attr `.no_translated`)")

    def pull_mz_quality_df(self):
        """Builds df, if quality columns exist"""
        logger.info("gets quality of the annotations (per ion and formula)")
        try:
            self.mz_quality = helpers.pick_mz_quality_infos(
                self.preloaded, self.mz_col, self.offsample_column,
                self.total_intens_col, self.fdr_column, self.msm_column,
                self.formula_col
            )
            if self.mz_quality[self.offsample_column].isnull().all():
                logger.warning("Useless: %s is NaN",
                               self.offsample_column)
            else:
                if not self.mz_quality[self.offsample_column].max() <= 1:
                    raise ValueError("Values of %s must be comprised "
                                     "between 0 and 1", self.offsample_column)
        except KeyError:
            logger.error(
                "Failed, one or more of required columns (%s) are missing",
                [self.offsample_column, self.total_intens_col,
                 self.fdr_column, self.msm_column])

    def verify(self):
        if self.usable_id != "KEGG":
            self.df = self.no_translated
        try:
            excess_mz_bool: bool = helpers.identify_multi_mz_to_one_compound(
                self.df, self.mz_col, self.usable_id)
            if excess_mz_bool:
                logger.warning(message_multimz)
        except (KeyError, ValueError, TypeError):
            logger.warning("failed to identify multi mz cases")
        except Exception as e:
            logger.warning("Error: %s", e)

    def convert_ids_to_kegg(self, ref_table: pd.DataFrame):
        logger.info('Translating ids in the metabolome annotation'
                    ' table, to KEGG (into .df attribute)')
        if not self.usable_id == "KEGG":
            logger.error('The .usable_id must be "KEGG". Aborted.')
            raise ValueError
        assert np.all(np.array(list(ref_table.columns)) == np.array(
            columns_kegg_conversion_expected)), (
            f"Metabolite reference mappings columns must be {columns_kegg_conversion_expected}"
        )
        try:
            self.df = helpers.convert_metabo_annots_tokegg(
                self.no_translated, ref_table,
                self.molecule_ids_col, self.mz_col, self.usable_id)
            self.verify()
        except (KeyError, ValueError, TypeError):
            logger.warning("failed to convert identifiers to KEGG compounds")
        except Exception as e:
            logger.warning("Error: %s", e)

    def keep_best_mz_id_pairs(self, var_df: pd.DataFrame):
        logger.info(
            "Resolving 'multiple m/z' that map to 'one compound id':\n"
            "  keeps the two top better (m/z, id) pairs by quality criteria.")
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
        ) == len(self.df[
                     self.mz_col].unique()), "Error, no perfect match with var_df mz"

        self.df = helpers.best_mz_id_pairs_df(self.df, var_df,
                                              self.mz_col, self.usable_id,
                                              self.msm_column,
                                              self.fdr_column,
                                              self.offsample_column)

    def remove_redundant_aliases(self):
        logger.info("collapsing metabolites' identifiers aliases")
        if self.usable_id != "KEGG":
            aliases_df = pd.read_csv(os.path.join(library_path,
                                                  "resources",
                                                  "aliases-HMDB-05-07-24-LipidMaps-04-30-25.tsv"),
                                     sep="\t", header=0, index_col=None,
                                     low_memory=False)

            accepted_ids = helpers.aliases_handler(
                aliases_df, self.df, self.usable_id)

            self.df = self.df.loc[
                self.df[self.usable_id].isin(accepted_ids), :].copy()
            self.no_translated = self.no_translated.loc[
                self.no_translated[self.usable_id].isin(
                    accepted_ids), :].copy()
        else:
            logger.info("all is %s ids, no alias", self.usable_id)

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
            logger.error(
                "Error: insufficient columns in annotation dataframe,"
                " field delimiter by default is '%s' (you can change"
                " it by using .delimiter attribute)", self.delimiter)

        if len(set(expected_columns_by_df)) != len(set(
                list(df.columns)).intersection(set(expected_columns_by_df))):
            logger.error(
                "one or more columns of in your annotations are missing, "
                "expected ones are exactly %s", expected_columns_by_df)


class MetabolicPathways(BaseModel, extra='allow'):
    """Class for any type of metabolites set, whether pathways, or families"""
    data_dir: str
    subfolder: str
    file: str
    supplementary_file: Optional[str] = ''
    download: Union[str, None] = None  # "sspa-kegg" if download is desired
    undesired_substring: str = " - Homo sapiens (human)"  # in names
    set_id_col: str = "pathway_id"
    set_name_col: str = "pathway_name"
    molecule_id_type: str = "KEGG"  # or moleculeIds (identical to annots.usable_id)
    piled: pd.DataFrame = pd.DataFrame()
    sets_ids_df: pd.DataFrame = pd.DataFrame()
    purge: bool = True

    def preload_official(self):
        if not os.path.exists(os.path.join(
                self.data_dir, self.subfolder, self.file)):
            if self.download == "sspa-kegg":
                logger.info("writing sspa official gmt to data dir, wait ...")
                sspa.process_pathways.process_kegg(
                    organism="hsa", infile=None, download_latest=True,
                    filepath=os.path.join(self.data_dir, self.subfolder))
                logger.info("done")
            else:
                logger.info(
                    "only 'sspa-kegg' download is currently supported")
        file_path = os.path.join(self.data_dir, self.subfolder, self.file)
        try:
            official_gmt = pd.read_csv(
                file_path, header=None, sep='\t', index_col=False,
                low_memory=False)
            official_gmt = official_gmt.dropna(axis=1, how="all")
        except FileNotFoundError:
            logger.error("Pathway file not found: %s", file_path)
            raise
        except pd.errors.EmptyDataError:
            logger.error("Pathway file is empty: %s", file_path)
            raise
        except Exception as e:
            logger.error("Failed to load pathway file '%s': %s", self.file, e)
            raise

        official_gmt = helpers.prettify_gmt(official_gmt,
                                            self.undesired_substring)
        inconvenient_df = pd.read_csv(os.path.join(
            library_path, "resources", "inconvenient_list_wikipathways.csv"
        ), header=0, index_col=None)

        if self.purge:
            official_gmt = helpers.purge_pathways(
                official_gmt, not_metabolic_reactions_pathways,
                inconvenient_id_name_df=inconvenient_df
            )
        self.official = official_gmt

    def preload_custom(self):
        try:
            custom_gmt = pd.read_csv(
                os.path.join(self.data_dir, self.subfolder,
                             self.supplementary_file),
                header=None, sep='\t', index_col=False)
        except Exception as e:
            logger.error(e)
        first_row = custom_gmt.iloc[0, :]
        cus_message = (f" The first line of your custom gmt or mmt file:\n"
                       f"  abbreviation: '{first_row.iloc[0]}'; "
                       f"  name : '{first_row.iloc[1]}';  metabolites"
                       f"  (first three) : '{first_row.iloc[2:5].tolist()}'.")
        logger.info("%s", cus_message)
        self.custom = helpers.prettify_gmt(custom_gmt, None)

    def load(self):
        self.preload_official()
        dfs_dict = {'official': self.official}
        if self.supplementary_file != '':
            self.preload_custom()
            dfs_dict['custom'] = self.custom

        id_name_pathways_df, pws_piled = helpers.gmts_merge_pile_splitid(
            dfs_dict, self.set_id_col, self.set_name_col,
            self.molecule_id_type)

        self.sets_ids_df = id_name_pathways_df
        self.piled = pws_piled

        logger.info("success, attributes : piled, sets_ids_df")

    def get_name(self, a_pathway_id):
        pathway_name = self.sets_ids_df.loc[
            self.sets_ids_df[self.set_id_col] == a_pathway_id,
            self.set_name_col].item()
        return pathway_name

    def to_dict(self, option="id") -> Dict[str, List[str]]:
        if option == "id":
            mz_or_id = self.molecule_id_type
        elif option == "mz":
            mz_or_id = "mz"
        else:
            logger.error("unrecognized option, must be 'id' or 'mz'")
        return helpers.pws_piled_to_dict(
            self.piled, self.set_id_col, mz_or_identity_col=mz_or_id)

    def intersect_with_annotations(self, annots_df: pd.DataFrame,
                                   usable_id: str, mz_col: str):
        """
        Modifies in the Pathways object: only the elements that are
        in its intersection with the annotations dataframe.

        In other words, the annotated compounds that are mapped
        to the pathways are kept within this Pathways object.

        Also adds the 'mz' column to the .piled df within the Pathways obj
        """
        if not self.molecule_id_type == usable_id:
            logger.error("%s and %s are not identical. Aborted.",
                         self.molecule_id_type, usable_id)
            raise AssertionError
        self.piled = helpers.pathways_and_annots_join(
            self.piled[[self.set_id_col, self.molecule_id_type]],
            annots_df[[mz_col, usable_id]],
            self.molecule_id_type)
        self.intersected: bool = True

    def constrain(self, min_pathway_size, max_pathway_size, mode):
        """filters the Pathways object after intersection with annotations"""
        assert mode in ["default", "mz_aware"], "Error, mode not recognized"
        assert self.intersected, "Must first run: .intersect_with_annotations"
        if mode == "default":
            good_pathways_l = pathways_methods.yield_compliant_pathways(
                self.to_dict(option="id"), min_pathway_size, max_pathway_size)
        elif mode == "mz_aware":
            good_pathways_l = pathways_methods.yield_compliant_pathways(
                self.to_dict(option="mz"), min_pathway_size, max_pathway_size)

        metabolites_to_eject = pathways_methods.find_metabolites_to_exclude(
            good_pathways_l, self.to_dict())
        self.piled = self.piled.loc[
            self.piled[self.set_id_col].isin(good_pathways_l)]
        self.piled = self.piled.loc[~self.piled[
            self.molecule_id_type].isin(metabolites_to_eject), :]
        self.sets_ids_df = self.sets_ids_df.loc[
            self.sets_ids_df[
                self.set_id_col].isin(good_pathways_l), :]

    def eject_molecules(self, ejection_list: List[str | float], column):
        """
        From sets object, delete specified mz, or molecules identifiers
        :input column: str, one of "KEGG", "moleculeIds", "mz"
        :input ejection_list: List of identifiers (list of strings)
         or mz (list of floats) to be removed from the object

        Note: this can be useful if having a list of compounds that were 'mapped-to-pathways',
        and willing to exclude them from the families sets object at any moment,
        but that is facultative
        """
        if column not in ["KEGG", "moleculeIds", "mz"]:
            logger.error("column '%s' not supported", column)
        self.piled = self.piled.loc[~self.piled[
            column].isin(ejection_list), :]
        self.sets_ids_df = self.sets_ids_df.loc[self.sets_ids_df[
            self.set_id_col].isin(self.piled[self.set_id_col].unique()), :]
