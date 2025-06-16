import os
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import squidpy as sq
from tqdm import tqdm
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from typing import List

import json


import zonist.data as zd
import zonist.processing.measures_methods as zmm



if __name__ == "__main__":
    # define here my input
    mydata_dir = "../../../data"
    my_csv_table_path = "../../privateNEW-datasets-parameters.tsv"
    plot_parent_folder = "../../figures/"

    # start
    datasets_parameters_df = pd.read_csv(my_csv_table_path, header=0,
                      index_col=None, sep='\t')

    datasets_parameters_df = datasets_parameters_df.loc[
        datasets_parameters_df['dataset_name'].isin(
            ["TD-Lipido-636"]), :]  #  "TD-Lipido-635", "TD-Lipido-836", "TD-Lipido-827"


    print("\nRunning almost no relevant test")

    for i, row in datasets_parameters_df.iterrows():
        SAMPLE = row["dataset_name"]
        allowed_in_out_ratio = row["in/out intensities for spot class"]
        select_spots_yes_no = row["perform_spots_selection"]
        print("\nSample: ", SAMPLE)
        subfolder = SAMPLE
        plot_folder = os.path.join(plot_parent_folder, SAMPLE)
        if not os.path.exists(plot_folder):
            os.mkdir(plot_folder)

        adata_mz = ad.read_h5ad(filename=os.path.join(
            mydata_dir, subfolder, f'{SAMPLE}.h5ad'))

        adata_mz.X = csr_matrix(zmm.tic_normalization(
            adata_mz.X.copy(), adata_mz.obs['obs'].tolist(),
            adata_mz.var['mz'].tolist()))

        sc.pp.highly_variable_genes(
            adata_mz, flavor="seurat_v3",
            check_values=False)
        sc.tl.pca(adata_mz, svd_solver="arpack")
        sc.pp.neighbors(adata_mz)
        sc.tl.umap(adata_mz)
        sc.tl.leiden(adata_mz, resolution=1, flavor="igraph",
                     n_iterations=2, directed=False)

        plot_out_path = os.path.join(plot_folder,
                                     f'{SAMPLE}-lala.png')

        sq.pl.spatial_scatter(
            adata_mz, color=["leiden"], img=False,
            title="mmmmm", save=os.path.abspath(plot_out_path), dpi=600)
        plt.close()
        print(f"plot saved to: {plot_out_path}")












