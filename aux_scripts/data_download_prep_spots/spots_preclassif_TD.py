# johaGL 2025
# script for cleaning spots in semimanual fashion for specific tissues
# usage example in `data_dps.sh`
import os
import pandas as pd
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
from typing import List
import scipy.sparse.linalg as linalg
import SpacePath.core.method as spm
import json


import squidpy as sq

from spots_purge import show_message

def outside_tissue_spots__to_json(
        spots_domains_df, out_tissue_domains: List[str],
        data_centralized_folder, subfolder, SAMPLE)-> None:
    spots_domains_df['obs'] = spots_domains_df['obs'].astype(str).tolist()
    spots_out_list = spots_domains_df.loc[
        spots_domains_df['squidpy_domains'].isin(out_tissue_domains),
        "obs"].unique()

    output_d =  {'out_tissue_spots': list(spots_out_list)}

    with open(os.path.join(data_centralized_folder, subfolder,
                           f'{SAMPLE}_out_tissue.json'), "w") as f:
        json.dump(output_d, f)



def run_custom_spot_qual(datasets_parameters_df,
                         data_centralized_folder,
                         figures_centralized_folder, weights_spatial_graph):
    show_message("\nRunning custom spot classification as being OUTSIDE the tissue,"
          "this script does not changes the AnnData object, but only saves .json file with the list of bad spots")

    for i, row in datasets_parameters_df.iterrows():
        SAMPLE = row["dataset_name"]
        show_message(f"\nSample: {SAMPLE}")
        subfolder = SAMPLE
        plot_folder = os.path.join(figures_centralized_folder, SAMPLE)
        if not os.path.exists(plot_folder):
            os.mkdir(plot_folder)

        if not os.path.exists(os.path.join(data_centralized_folder, subfolder,
                            f'{SAMPLE}_out_tissue.json')):
            adata_mz = ad.read_h5ad(filename=os.path.join(
                data_centralized_folder, subfolder, f'{SAMPLE}.h5ad'))
            adata_mz = spm.NormalizeIntensities.run(adata_mz)  # just for computations, but not saved

            try:
                sc.tl.pca(adata_mz, svd_solver="arpack")
            except linalg._eigen.arpack.ArpackError:
                # unfrequent error, when the first column vector variance is zero.
                # to rescue, dropping the zero variance features
                Xdf = adata_mz.to_df()
                feats_nonzero_var =  Xdf.loc[:, Xdf.var() != 0].columns.tolist()
                adata_mz = adata_mz[:, feats_nonzero_var].copy()
                del Xdf
                sc.tl.pca(adata_mz, svd_solver="arpack")

            sc.pp.neighbors(adata_mz)
            nn_graph_genes = adata_mz.obsp["connectivities"]
            sq.gr.spatial_neighbors(adata_mz)
            nn_graph_space = adata_mz.obsp["spatial_connectivities"]

            alpha = weights_spatial_graph[SAMPLE]

            joint_graph = (1 - alpha) * nn_graph_genes + alpha * nn_graph_space
            sc.tl.leiden(adata_mz,
                         adjacency=joint_graph,
                         n_iterations=2,
                         flavor="igraph",
                         directed=False,
                         key_added="squidpy_domains", resolution=0.8)

            plot_out_path = os.path.join(plot_folder,
                                         f'{SAMPLE}-spotsCustomViz.png')

            sq.pl.spatial_scatter(
                adata_mz, color=["squidpy_domains"], img=False,
                title="check spots ", save=os.path.abspath(plot_out_path), dpi=600)
            plt.close()
            show_message(f"plot saved to: {plot_out_path}, revise it")

            show_message(adata_mz.obs[['obs', 'squidpy_domains']].head())

            answer = input(
                "enter the BAD (outside tissue) domains numbers (example: 8,6,7)"
                " separated by commas without space : ")
            out_tissue_domains = [str(n) for n in list(set(answer.strip().split(",")))]

            outside_tissue_spots__to_json(adata_mz.obs[["obs", "squidpy_domains"]],
                                          out_tissue_domains,
                                          data_centralized_folder, subfolder, SAMPLE)

            show_message("Spots outside tissue saved to .JSON file")
        else:
            show_message("nothing to do")


if __name__ == "__main__":
    # inputs must be immutable, custom procedure of this script
    data_centralized_folder = "../../../data"
    my_csv_table_path = "../../datasets_parameters/private-datasets-parameters.tsv"
    figures_centralized_folder = "../../../figures/"

    # start
    datasets_parameters_df = pd.read_csv(my_csv_table_path, header=0,
                      index_col=None, sep='\t')

    datasets_parameters_df = datasets_parameters_df.loc[
        datasets_parameters_df['dataset_name'].isin(
            ["TD-Lipido-635", "TD-Lipido-636", "TD-Lipido-836", "TD-Lipido-827"]), :]  #

    weights_spatial_graph = {
        "TD-Lipido-635" : 0.5,
        "TD-Lipido-636" : 0.5,
        "TD-Lipido-836" : 0.7,  # had to use this to force the internal areas not to fall in same cluster of external areas
        "TD-Lipido-827" : 0.6   # same as for line above
    }

    answer = input("continue? y/n")
    if answer in ["y", "yes"]:
        run_custom_spot_qual(datasets_parameters_df,
                             data_centralized_folder,
                             figures_centralized_folder, weights_spatial_graph )
