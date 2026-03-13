# johaGL 2025
# script for cleaning spots in semimanual fashion for specific tissues

import os
import pandas as pd
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
from typing import List
import scipy.sparse.linalg as linalg
import smpath.method as spm
import json

import warnings

warnings.simplefilter('ignore', category=FutureWarning)
import squidpy as sq

from spots_purge import show_message
from spots_preclassif_TD import run_custom_spot_qual

data_centralized_folder = "../../../data"
my_csv_table_path = "../../datasets_parameters/public-datasets-parameters.tsv"
figures_centralized_folder = "../../../figures/"

# start
datasets_parameters_df = pd.read_csv(my_csv_table_path, header=0,
                                     index_col=None, sep='\t')

datasets_parameters_df = datasets_parameters_df.loc[
    datasets_parameters_df['dataset_name'].isin(
        ["Brain-healthy-mouse", "Brain-lmrange-mouse"]), :].copy()  #

weights_spatial_graph = {
    "Brain-healthy-mouse": 0.3,
    "Brain-lmrange-mouse" : 0.3
}
answer = input("continue? y/n")
if answer in ["y", "yes"]:
    run_custom_spot_qual(datasets_parameters_df,
                data_centralized_folder,
                figures_centralized_folder, weights_spatial_graph)

