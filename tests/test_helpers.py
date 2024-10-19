from unittest import TestCase
import src.zonist.helpers as helpers
import pandas as pd


class Test(TestCase):
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
    result = helpers.newfun(metabolome_annot_df,
                            ref_ids_map_kegg,
                            col_molecule_ids="moleculeIds", col_mz="mz")
    def test_newfun(self):
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

