from unittest import TestCase

import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

import src.smpath.processing.measures_methods as pmm


class Test(TestCase):
    def test_tic_normalization(self):
        intensities_matrix = csr_matrix(np.array(
            [[100.44, 400.55, 700.33],
             [200.2, 500.7, 800.5],
             [0., 0., 0.],
             [300.4, 600.1, 900.3]]))

        result = pmm.tic_normalization(intensities_matrix)
        result = result.toarray()  # yes this method is correct, test ok

        self.assertTrue(np.all(
            np.array(np.around(result[0, :], 4)).flatten(
            ) == np.array([112.8906, 450.2023, 787.1431])))
        self.assertTrue(np.all(
            np.array(np.around(result[1, :], 4)).flatten(
            ) == np.array([180.0734, 450.3634, 720.0237])))
        self.assertTrue(np.all(
            np.array(np.around(result[2, :], 4)).flatten(
            ) == np.array([0, 0, 0])))
        self.assertTrue(np.all(
            np.array(np.around(result[3, :], 4)).flatten(
            ) == np.array([225.3017, 450.0783, 675.23])))

    def test_compute_metabolite_probs(self):
        a_mat = np.array([[100, 90, 300], [200, 45, 90], [100, 25, 81],
                          [39, 17, 20], [13, 40, 71]])
        obs = pd.DataFrame(
            index=["spot1", "spot2", "spot3", "spot4", "spot5"])
        my_mzs = ['400.55', '577.5', '700.33']
        var = pd.DataFrame({'mz': my_mzs}, index=my_mzs)
        adata_mz = ad.AnnData(X=csr_matrix(a_mat), obs=obs, var=var)
        annots_df = pd.DataFrame({
            "mz": np.array(['400.55', '400.55', '577.5', '700.33', '700.33']),
            "moleculeIds": ['hmdbA', 'hmdbB', 'hmdbY', 'lmpV', 'lmpW']
        })
        result = pmm.compute_metabolites_probs_df(
            adata_mz, annots_df, "moleculeIds", "mz")
        #          hmdbA  hmdbB  hmdbY  lmpV  lmpW
        # spot1    0.5    0.5    1.0   0.5   0.5
        # spot2    0.5    0.5    1.0   0.5   0.5
        # spot3    0.5    0.5    1.0   0.5   0.5
        # spot4    0.5    0.5    1.0   0.5   0.5
        # spot5    0.5    0.5    1.0   0.5   0.5
        self.assertEqual(result.shape, (5, 5))
        self.assertTrue(
            np.all(result.loc['spot1', :].to_numpy() == np.array(
                [0.5, 0.5, 1, 0.5, 0.5])))
