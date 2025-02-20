from unittest import TestCase
from typing import List, Dict
import src.zonist.helpers as helpers

import pandas as pd
import numpy as np


class Test(TestCase):

    def test_convert_metabo_annots(self):
        metabolome_annot_df = pd.DataFrame({
            "mz": ["102.0317", "132.0899", "913.821868", "132.0899"],
            "formula": ["C4H6O3", "C5H12N2O2", "C59H108O6", "C5H12N2O2"],
            "moleculeIds": ["HMDB0001259", "C00077", "LMGL03010536",
                            "HMDB0000214"],
            "moleculeNames": ["Succinic semialdehyde", "L-Ornithine",
                              "TG (17..)..", "Ornithine"],
            "database": ["HMDB", "KEGG", "LipidMaps", "HMDB"]
        })
        ref_ids_map_df = pd.DataFrame({
            "initial_id_unconverted": ["HMDB0001259", "HMDB0000214", "LMGL03010536",
                            "LMTG9999", "LMPE8888"],
            "KEGG": ["C00232", "C00077", "CXXXXYY", "C0101010", "C999999"],
            "database": ["HMDB", "HMDB", "LipidMaps", "LipidMaps", "HMDB"]
        })
        result = helpers.convert_metabo_annots_tokegg(
            metabolome_annot_df, ref_ids_map_df,
            molecule_ids_col="moleculeIds", mz_col="mz")
        self.assertListEqual(
            list(result.columns),
            ['mz', 'formula', 'moleculeIds', 'moleculeNames',
             'KEGG', 'database'])
        self.assertListEqual(
            result['moleculeNames'].tolist(),
            ["Succinic semialdehyde", "L-Ornithine", "TG (17..).."]
        )
        self.assertListEqual(result['KEGG'].tolist(),
                             ["C00232", "C00077", "CXXXXYY"])

        self.assertListEqual(result['database'].tolist(),
                             ["HMDB", "KEGG", "LipidMaps"])

    def test_prettify_gmt(self):
        data = [
            ['IDpath1', 'pathway-one-!?', 'm1', 'm4', 'm6', np.nan, np.nan],
            ['IDpath2', 'pathway-two-!?', 'm2', 'm3', np.nan, np.nan, np.nan],
            ['IDpath3', 'pathway-three-!?', 'm3', 'm4', 'm5',
             'm7', 'm8'],
            ['IDpath4', 'pathway--no-elements-!?', np.nan, np.nan, np.nan,
             np.nan, np.nan]]

        gmt = pd.DataFrame(data, columns=range(0, 7))
        result = helpers.prettify_gmt(gmt_df=gmt,
                                      undesired_name_substring='-!?')
        self.assertTrue('IDpath4' not in result.iloc[:, 0].tolist())
        self.assertListEqual(['IDpath1', 'IDpath2', 'IDpath3'],
                             result.iloc[:, 0].tolist())
        self.assertListEqual(
            ['pathway-one', 'pathway-two', 'pathway-three'],
            result.iloc[:, 1].tolist())

    def test_purge_pathways(self):  #
        data = [
            ['IDpath1', 'Fatty Oxidation', 'm1', 'm3', 'm4', 'm6',  np.nan],
            ['IDpath2', 'Meiosis', 'm2', 'm3', np.nan, np.nan, np.nan],
            ['IDpath3', 'TCA cycle', 'm3', 'm4', 'm5', 'm7', np.nan],
            ['IDpath4', 'Yersinia Infection', 'm8', 'm9', 'm15', 'm1', 'm2']]
        not_strictly_metabolical =  ["meiosis", "infection"]
        numerical_column_names = range(0, 7)
        gmt = pd.DataFrame(data, columns=numerical_column_names)
        result = helpers.purge_pathways(gmt, not_strictly_metabolical)
        self.assertListEqual(['IDpath1',  'IDpath3'],
                             result.iloc[:, 0].tolist())
        self.assertListEqual(
            ['Fatty Oxidation', 'TCA cycle'],
            result.iloc[:, 1].tolist())

    def test_gmt_to_piled_form(self):
        data = [
            ['IDpath1', 'pathway-one', 'm1', 'm4', 'm6', np.nan, np.nan],
            ['IDpath2', 'pathway-two', 'm2', 'm3', np.nan, np.nan, np.nan],
            ['IDpath3', 'pathway-three', 'm3', 'm4', 'm5',
             'm7', 'm8']]
        gmt = pd.DataFrame(data, columns=range(0, 7))
        result = helpers.gmt_to_piled_form(gmt, "pathway_id", "KEGG")
        self.assertTrue(result.shape[1] == 2)
        self.assertTrue(result.shape[0] == 10)
        self.assertTrue(np.all(np.array(["IDpath3", "m8"]) == (
            result.iloc[9, :].to_numpy()
        )))

    def test_gmt_join_mz(self):
        pws_piled = pd.DataFrame({
            "pathway_id": ["one", "one", "one", "one", "two", "two"],
            "KEGG": ["C00232", "C09090", "C0", "C00077", "CKKKK", "CXXXXYY"]
        })
        my_annotations = pd.DataFrame({
            "mz": ["102.0317", "132.0899", "588.5665543", "913.821868"],
            "KEGG": ["C00232", "C00077", "C09090", "CXXXXYY"]
        })
        result = helpers.pathways_and_annots_join(pws_piled, my_annotations,
                                                  mol_id_in_pws="KEGG")
        # result : 
        #    pathway_id    KEGG      mz
        #    one        C00232    102.0317
        #    one        C09090    588.5665543
        #    one        C00077    132.0899
        #    two        CXXXXYY   913.821868
        self.assertListEqual(result.iloc[3, :].tolist(),
                             ["two", "CXXXXYY", "913.821868"])

    def test_pws_piled_to_dict(self):
        data = [["one", "C00232", "102.0317"],
                ["one", "C09090", "588.56"],
                ["one", "C00077", "132.0899"],
                ["two", "CXXXXYY", "913.821868"]]
        pws_piled = pd.DataFrame(data, columns=["pathway_id", "KEGG", "mz"])

        result = helpers.pws_piled_to_dict(pws_piled, pathway_id="pathway_id",
                                           mz_or_identity_col="mz")
        self.assertListEqual(list(result.keys()),
                             ["one", "two"])
        self.assertTrue(np.all(np.array(result["one"]) == np.array(
                ["102.0317", "588.56", "132.0899"])))

    def test_generate_coordinates_arr(self):
        coords_l: List[str] = ['x4_y1', 'x4_y2', 'x5_y3', 'x6_y8', 'x6_y9']
        result = helpers.generate_coordinates_arr(coords_l)
        uf = result[:, 0]
        self.assertTrue(np.all(result[:, 0] == np.array([4, 4, 5, 6, 6])))
        self.assertTrue(np.all(result[:, 1] == np.array([1, 2, 3, 8, 9])))


