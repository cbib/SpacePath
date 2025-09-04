from unittest import TestCase
import pandas as pd
from SPACEPATH.src.single_slice_helpers import (prepare_molecular_ids_df,
                                                get_top_n_simple,
                                                get_top_n,
                                                differential_results_filter)
from typing import List, Dict
import numpy as np


class Test(TestCase):
    def test_prepare_molecular_ids_df(self):
        top_differential = pd.DataFrame({
            'names': ['100.0', '200.0', '300.0', '100.0'],  # String mzs
            'logfoldchanges': [2, 1.9, 3, 6], 'padj': [0, 0, 0, 0]
        })

        clean_ann_df = pd.DataFrame({
            'moleculeIds': ['mol1', 'mol2', 'mol3', 'mol4', None],
            'mz': [100.0, 200.0, 400.0, 300.0, 999.9],
            'moleculeNames': ['name1', 'name2', 'name3', 'name4', 'name5']
        })
        result = prepare_molecular_ids_df(top_differential, clean_ann_df)
        self.assertIsInstance(result, pd.DataFrame)

        expected_mzs = [100.0, 200.0, 300.0]
        self.assertTrue(all(mz in expected_mzs for mz in result['mz']))
        self.assertNotIn(None, result['mz'].values)
        self.assertEqual(result.shape[0], len(result['mz'].unique()))

    def test_get_top_n_simple(self):
        data = {
            'pvals_adj': [0.01, 0.9, 0.02, 0.01, 0.0001],
            'abslfc': [10, 20, 15, 25, 5],
            'cluster': ['A', 'A', 'A', 'A', 'A']
        }
        df = pd.DataFrame(data)
        result = get_top_n_simple(df, 3)
        self.assertEqual(len(result), 3)
        self.assertEqual(result.iloc[0]['pvals_adj'],
                         min(df['pvals_adj']))
        sorted_pvals = sorted(df['pvals_adj'])
        self.assertEqual(result.iloc[1]['pvals_adj'], sorted_pvals[1])

    def test_get_top_n(self):
        data = {
            'pvals_adj': [0.01, 0.05, 0.02, 0.01, 0.0001, 0.01, 0.01, 0.01],
            'logfoldchanges': [-10, -1, 15, -2, -15, -3, -5, 7],
            'cluster': ['A', 'A', 'A', 'A', 'A', 'A', 'A', 'A']
        }
        df = pd.DataFrame(data)
        result1 = get_top_n(df, 3, "down")
        result2 = get_top_n(df, 7, "up")
        print(result2)
        self.assertTrue(np.all(np.around(
            result1['logfoldchanges'].to_numpy(), 0) == np.array(
            [-15, -10, -5]
        )))
        self.assertTrue(np.all(np.around(
            result2['logfoldchanges'].to_numpy(), 0) == np.array(
            [7, 15]
        )))

    def test_differential_results_filter(self):
        data = {
            'group': ['A', 'A', 'B', 'B', 'C', 'C'],
            'logfoldchanges': [0.9, 0.5, -1.2, 0.3, 1.0, -0.7],
            'pvals_adj': [0.01, 0.2, 0.04, 0.06, 0.03, 0.5]
        }
        df = pd.DataFrame(data)

        result = differential_results_filter(df, lfc=0.8, pvals_adj=0.05,
                                             total_n_features=3)
        #  group   logfoldchanges  pvals_adj  abslfc
        #    A            0.9       0.01      0.9
        #    B           -1.2       0.04      1.2
        #    C            1.0       0.03      1.0

        self.assertTrue(all(abs(result['logfoldchanges']) >= 0.8))
        self.assertTrue(all(result['pvals_adj'] <= 0.05))



