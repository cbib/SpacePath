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
import zonist.method as zm
import json


def outside_tissue_spots__to_json(
        spots_domains_df, out_tissue_domains: List[str],
        mydata_dir, subfolder, SAMPLE)->None:
    spots_domains_df['obs'] = spots_domains_df['obs'].astype(str).tolist()
    spots_out_list = spots_domains_df.loc[
        spots_domains_df['squidpy_domains'].isin(out_tissue_domains),
        "obs"].unique()

    output_d =  {'out_tissue_spots': list(spots_out_list)}

    with open(os.path.join(mydata_dir, subfolder,
                           f'{SAMPLE}_out_tissue.json'), "w") as f:
        json.dump(output_d, f)



if __name__ == "__main__":
    # define here my input
    mydata_dir = "../../../data"
    my_csv_table_path = "../../private-datasets-parameters.tsv"
    plot_parent_folder = "../../../figures/"

    # start
    datasets_parameters_df = pd.read_csv(my_csv_table_path, header=0,
                      index_col=None, sep='\t')

    datasets_parameters_df = datasets_parameters_df.loc[
        datasets_parameters_df['dataset_name'].isin(
            [ "TD-Lipido-635", "TD-Lipido-636"  , "TD-Lipido-836", "TD-Lipido-827"]), :]  #

    weights_spatial_graph = {
        "TD-Lipido-635" : 0.5,
        "TD-Lipido-636" : 0.5,
        "TD-Lipido-836" : 0.7,  # had to use this to force the internal areas not to fall in same cluster of external areas
        "TD-Lipido-827" : 0.6   # same as for line above
    }

    print("\nRunning custom spot classification as being OUTSIDE the tissue,"
          "this script does not changes the AnnData object")

    for i, row in datasets_parameters_df.iterrows():

        SAMPLE = row["dataset_name"]
        allowed_in_out_ratio = row["in/out intensities for spot class"]
        select_spots_yes_no = row["perform_spots_selection"]
        print("\nSample: ", SAMPLE)
        subfolder = SAMPLE
        plot_folder = os.path.join(plot_parent_folder, SAMPLE)
        if not os.path.exists(plot_folder):
            os.mkdir(plot_folder)

        if not os.path.exists(os.path.join(mydata_dir, subfolder,
                            f'{SAMPLE}_out_tissue.json')):
            adata_mz = ad.read_h5ad(filename=os.path.join(
                mydata_dir, subfolder, f'{SAMPLE}.h5ad'))
            adata_mz = zm.NormalizeIntensities.run(adata_mz)

            # PCA
            print("PCA")
            sc.tl.pca(adata_mz, svd_solver="arpack")
            print("nearest neighbor graph")
            sc.pp.neighbors(adata_mz)
            nn_graph_genes = adata_mz.obsp["connectivities"]
            print("spatial proximity graph")
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
            print(f"plot saved to: {plot_out_path}, revise it")

            print(adata_mz.obs[['obs', 'squidpy_domains']].head())

            answer = input(
                "enter the BAD (outside tissue) domains numbers (example: 8,6,7)"
                " separated by commas without space : ")
            out_tissue_domains = [str(n) for n in list(set(answer.strip().split(",")))]

            outside_tissue_spots__to_json(adata_mz.obs[["obs", "squidpy_domains"]],
                                          out_tissue_domains,
                                          mydata_dir, subfolder, SAMPLE)

            print("Spots outside tissue saved to .JSON file")
        else:
            print("nothing to do")
