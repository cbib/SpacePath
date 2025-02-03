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
        ref_ids_map_kegg = pd.DataFrame({
            "not_kegg_id": ["HMDB0001259", "HMDB0000214", "LMGL03010536",
                            "LMTG9999", "LMPE8888"],
            "KEGG": ["C00232", "C00077", "CXXXXYY", "C0101010", "C999999"],
            "database": ["HMDB", "HMDB", "LipidMaps", "LipidMaps", "HMDB"]
        })
        result = helpers.convert_metabo_annots_tokegg(
            metabolome_annot_df, ref_ids_map_kegg,
            col_molecule_ids="moleculeIds", col_mz="mz")
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

    def test_clean_gmt(self):
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

    def test_gmt_to_piled_form(self):
        data = [
            ['IDpath1', 'pathway-one', 'm1', 'm4', 'm6', np.nan, np.nan],
            ['IDpath2', 'pathway-two', 'm2', 'm3', np.nan, np.nan, np.nan],
            ['IDpath3', 'pathway-three', 'm3', 'm4', 'm5',
             'm7', 'm8']]
        gmt = pd.DataFrame(data, columns=range(0, 7))
        result = helpers.gmt_to_piled_form(gmt, "path_id", "KEGG")
        self.assertTrue(result.shape[1] == 2)
        self.assertTrue(result.shape[0] == 10)
        self.assertTrue(np.all(np.array(["IDpath3", "m8"]) == (
            result.iloc[9, :].to_numpy()
        )))

    def test_gmt_join_mz(self):
        gmt_piled = pd.DataFrame({
            "path_id": ["one", "one", "one", "one", "two", "two"],
            "KEGG": ["C00232", "C09090", "C0", "C00077", "CKKKK", "CXXXXYY"]
        })
        anno_mz_kegg = pd.DataFrame({
            "mz": ["102.0317", "132.0899", "588.5665543", "913.821868"],
            "KEGG": ["C00232", "C00077", "C09090", "CXXXXYY"]
        })
        result = helpers.pathways_and_annots_join(gmt_piled, anno_mz_kegg,
                                                  col_molecule_id="KEGG")
        # result : 
        #    path_id    KEGG      mz
        #    one        C00232    102.0317
        #    one        C09090    588.5665543
        #    one        C00077    132.0899
        #    two        CXXXXYY   913.821868
        self.assertListEqual(result.iloc[3, :].tolist(),
                             ["two", "CXXXXYY", "913.821868"])

    def test_gmt_piled_to_dict(self):
        data = [["one", "C00232", "102.0317"],
                ["one", "C09090", "588.56"],
                ["one", "C00077", "132.0899"],
                ["two", "CXXXXYY", "913.821868"]]
        gmt_piled = pd.DataFrame(data, columns=["path_id", "KEGG", "mz"])

        result = helpers.gmt_piled_to_dict(gmt_piled,
                                           pathway_id="path_id",
                                           mz_or_kegg="mz")
        self.assertListEqual(list(result.keys()),
                             ["one", "two"])
        self.assertTrue(np.all(
            np.array(result["one"]) == np.array(
                ["102.0317", "588.56", "132.0899"])))

    def test_generate_coordinates_arr(self):
        coords_l: List[str] = ['x4_y1', 'x4_y2', 'x5_y3', 'x6_y8', 'x6_y9']
        result = helpers.generate_coordinates_arr(coords_l)
        uf = result[:, 0]
        self.assertTrue(np.all(result[:, 0] == np.array([4,4,5,6,6])))
        self.assertTrue(np.all(result[:, 1] == np.array([1,2,3,8,9])))


