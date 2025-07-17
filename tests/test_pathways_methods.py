from unittest import TestCase
from typing import Dict, List
import numpy as np
import pandas as pd
# import anndata as ad
import src.smpath.helpers as helpers
import src.smpath.processing.pathways_methods as pathways_methods


class PathwaysForTests:
    """
    class to use in replacement of MetabolicPathways for faster testing
    """

    def __init__(self):
        self.piled = pd.DataFrame({
            "pathway_id": ["p1", "p1", "p1", "p1", "p2", "p2", "p2"],
            "KEGG": ["9", "8", "7", "6", "5", "4", "3"]
        }).astype(str)

    def to_dict(self, mz_or_id="KEGG") -> Dict[str, List[str]]:
        return helpers.pws_piled_to_dict(
            self.piled, "pathway_id", mz_or_identity_col=mz_or_id)


class Test(TestCase):
    def test_spots_groups_by_random_choice(self):
        spots_names = [str(i) for i in np.arange(13)]
        result = pathways_methods.spots_groups_by_random_choice(
            spots_names, n_groups=4)
        self.assertEqual(len(result[0]), 3)
        self.assertEqual(len(result[1]), 3)
        self.assertEqual(len(result[2]), 3)
        self.assertEqual(len(result[3]), 4)

    def test_reshape_pathways_for_kpca(self):
        pathways_dict = {"p1": ["9", "8", "7", "6", ],
                         "p2": ["5", "4", "3"]}
        result = pathways_methods.reshape_pathways_for_kpca(pathways_dict)
        self.assertIsInstance(result, pd.DataFrame)
        self.assertListEqual(result.index.tolist(), ['p1', 'p2'])
        self.assertTrue(np.all(result.loc["p1", :].to_numpy() ==
                               np.array(["9", "8", "7", "6", "none"])))

    def test_verticalassembly_of_dicts(self):
        df1 = pd.DataFrame(data=np.array([[55, 88, 6, -1], [0, 2, 4, 0],
                                          [78, -2, 0, 0], [78, 0, 0, 0]]),
                           columns=["c05", "c04", "c11", "c12"],
                           index=['a', 'u', 'v', 'o'])
        df2 = pd.DataFrame(data=np.array([[1, 0, 0], [0, 2, 0], [78, 0, 0]]),
                           columns=["c10", "c11", "c12"],
                           index=['e', 'i', 'w'])
        dfs_dicts = {0: df1, 1: df2}
        ordered_spots_names = ['a', 'e', 'i', 'o', 'u', 'v', 'w']
        result = pathways_methods.gather_scores_dfs_dict(
            dfs_dicts, ordered_spots_names)
        self.assertIsInstance(result, pd.DataFrame)
        self.assertListEqual(result.columns.tolist(),
                             ['c04', 'c05', 'c10', 'c11', 'c12'])
        self.assertListEqual(result.index.tolist(), ordered_spots_names)
        self.assertTrue(result.min().min() == -2)

    def test_compute_pathways_scores(self):
        data_df = pd.DataFrame(data=np.array([[600, 6477, 4545, 1e4],
                                             [1e4, 23423, 7000, 1e3],
                                             [1e3, 95345, 1e4, 4353]]),
                      columns=["lipid1", "lipid2", "lipid3", "lipid4"],
                      index=["spotA", "spotB", "spotC"])
        sets_dict = {'Set1': ["lipid3", "lipid1"], 'Set2': ["lipid2", "lip6"],
                     'Set3': ["lipid2", "lipid4", "lipid3"], 'Set4': ["lip9"]}
        result = pathways_methods.compute_pathways_scores(
            sets_dict, data_df, method='sspa_kpca', min_entity=2,
            n_groups_of_spots=30, threshold_n_spots_is_high=1000
        )
        result = result.round(4)
        self.assertListEqual(result.columns.tolist(),
                             ['Set1', 'Set3'])
        self.assertTrue(np.all(
            result.loc['spotA',:].to_numpy() == np.array([-0.4661, 0.8059]))
        )
        self.assertTrue(np.all(
            result.loc['spotB', :].to_numpy() == np.array([0.8787, -0.3404])),
        )
        self.assertTrue(np.all(
            result.loc['spotC', :].to_numpy() == np.array([-0.4126, -0.4654])),
        )
