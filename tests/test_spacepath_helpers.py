from unittest import TestCase
import pandas as pd
from SpacePath.src.helpers import (prepare_molecular_ids_df,
                           get_top_n,
                           differential_results_filter,
                           score_as_adata,
                           squidpy_differential_analysis)
import anndata as ad
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

        result = differential_results_filter(df, abs_lfc=0.8, pvals_adj=0.05,
                                             total_n_features_out=3)
        #  group   logfoldchanges  pvals_adj  abslfc
        #    A            0.9       0.01      0.9
        #    B           -1.2       0.04      1.2
        #    C            1.0       0.03      1.0

        self.assertTrue(all(abs(result['logfoldchanges']) >= 0.8))
        self.assertTrue(all(result['pvals_adj'] <= 0.05))

    def test_score_as_adata(self):
        scores_df = pd.DataFrame(data=np.array([[-2, 0, 1.4],
                                                [0, 2, -4], [-2, 1, 0]]),
                                 index=["x1y0", "x0y0", "x1y1"],
                                 columns=['MolFam1', 'hsa97', 'WP99'])
        sets_ids_df = pd.DataFrame(
            {'pathway_id': ['hsa97', 'MolFam1', 'WP99'],
             'pathway_name': ['fatty ac digest', 'PC(22:2)', 'oxylipins']})

        pathways_piled = pd.DataFrame(
            {'pathway_id': ['hsa97', 'MolFam1', 'WP99'],
            'moleculeIds':  ['LM4545', 'LM1167', 'HMDB00'],
             'mz': [888, 765, 934]})

        coords = np.array(([1,0],[0,0],[1,1]))

        result = score_as_adata(scores_df, sets_ids_df,
            pathways_piled, coords)

        self.assertTrue(
            result.var.columns.tolist() == ['pathway_id', 'pathway_name'])
        self.assertEqual(result.obsm['spatial'].shape[0], coords.shape[0])

    def test_differential_analysis(self):
        X = np.random.rand(20, 20)  # 20 cells, 20 lipids
        obsdf = pd.DataFrame({
            'kmeans3': np.random.choice(['cluster1', 'cluster2'], size=20)
        })
        adata = ad.AnnData(X=X, obs=obsdf)
        result = squidpy_differential_analysis(
            adata, group_str='kmeans3')

        self.assertIsInstance(result, pd.DataFrame)
        self.assertGreater(len(result), 0)



