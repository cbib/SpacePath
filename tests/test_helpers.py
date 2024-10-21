from unittest import TestCase
import src.zonist.helpers as helpers
import pandas as pd


class Test(TestCase):
    
    def test_translate_metabo_annots(self):
        metabolome_annot_df = pd.DataFrame({
            "mz": [102.0317, 132.0899, 913.821868, 132.0899],
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
        result = helpers.translate_metabo_annots(metabolome_annot_df,
                                ref_ids_map_kegg,
                                col_molecule_ids="moleculeIds", col_mz="mz")
        self.assertListEqual(
            list(self.result.columns),
            ['mz', 'formula', 'moleculeIds', 'moleculeNames',
             'KEGG', 'database'])
        self.assertListEqual(
            self.result['moleculeNames'].tolist(),
            ["Succinic semialdehyde", "L-Ornithine", "TG (17..).."]
        )
        self.assertListEqual(self.result['KEGG'].tolist(),
                             ["C00232", "C00077", "CXXXXYY"])

        self.assertListEqual(self.result['database'].tolist(),
                             ["HMDB", "KEGG", "LipidMaps"])
        
    def test_clean_gmt(self):
        data = [
            ['IDpath1', 'pathway-one', 'm1', 'm4', 'm6', '', ''], 
            ['IDpath2', 'pathway-two', 'm2', 'm3', '', '', ''], 
            ['IDpath3', 'pathway-three', 'm3', 'm4', 'm5',
                'm7', 'm8'],
            ['IDpath4', 'pathway-four', '', '', '', '', '']]

        gmt = pd.DataFrame(data, columns=range(0, 7))
        result = helpers.clean_gmt(gmt_df=gmt, undesired_substring='')
        self.assertTrue('IDpath4' not in gmt.iloc[:, 0].tolist())
        self.assertListEqual(['IDpath1', 'IDpath2', 'IDpath3' ],
                             gmt.iloc[:, 0].tolist())
    
    def test_pull_pathway_ids_name_df(self):
        data = [
            ['IDpath1', 'pathway-one', 'm1', 'm4', 'm6', '', ''], 
            ['IDpath2', 'pathway-two', 'm2', 'm3', '', '', ''], 
            ['IDpath3', 'pathway-three', 'm3', 'm4', 'm5',
                'm7', 'm8']]

        gmt = pd.DataFrame(data, columns=range(0, 7))
        result = helpers.pull_pathway_ids_name_df(gmt_df=gmt)
        self.assertTrue(result.shape[1] == 2)
        self.assertListEqual(["pathwayId", "pathwayName"],
                             list(result.columns))



