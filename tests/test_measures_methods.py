from unittest import TestCase
from typing import Dict, List
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix
import src.zonist.helpers as helpers
import src.zonist.processing.measures_methods as pmm


class Test(TestCase):
    def test_tic_normalization(self):
        intensities_matrix = csr_matrix(np.array(
            [[100.44, 400.55, 700.33],
             [200.2, 500.7, 800.5],
             [0., 0., 0.],
             [300.4, 600.1, 900.3]]))
             # 3 features (m/z) x 4 spots (each spectrum is a row (spot))
        spots_names = ["s1", "s2", "s3", "s4"]
        features_names = [84.66, 88.31, 150.55]

        result = pmm.tic_normalization(intensities_matrix,
                                                  spots_names,
                                                  features_names)
        hoho = np.array(np.around(result[0, :], 4)
                     ).flatten()
        print(hoho, "::::")
        print(type(hoho))
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


