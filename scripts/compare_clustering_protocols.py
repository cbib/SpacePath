import os
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import scanpy as sc
import squidpy as sq
from sklearn.cluster import KMeans
from sklearn.utils import sparsefuncs
import time

np.random.seed(42)


# todo: add to the environment igraph !!!  did micromamba install -c conda-forge python-igraph


def give_time(elapsed_seconds):
    minutes = int(elapsed_seconds // 60)
    seconds = int(elapsed_seconds % 60)
    elapsed_remaining =  elapsed_seconds - (seconds + (minutes * 60))
    resm  = np.around(elapsed_remaining * 1000, 2)

    return f"{minutes}:{seconds}:{resm}"


""""
"""
if __name__ == "__main__":

    # define here my files and dirs

    mydata_dir = "../../data"
    plot_parent_folder = "../../figures/"

    my_csv_table_path_1 = "../old-public-datasets-parameters.tsv"  # "../../small_test.csv"  # todo: correct this file ?
    my_csv_table_path_2 = "../private-datasets-parameters.tsv"  # "../../small_test.csv"

    t0 = time.time()

    df1 = pd.read_csv(my_csv_table_path_1, header=0,
                      index_col=None, sep='\t')
    df2 = pd.read_csv(my_csv_table_path_2, header=0,
                      index_col=None, sep='\t')
    datasets_parameters_df = pd.concat([df1, df2])

    datasets_parameters_df = datasets_parameters_df.loc[
        datasets_parameters_df['dataset_name'].isin(
            ["TD-Lipido-635"])
    ]
    t1 = time.time()

    print(t1 - t0)

    time_str = give_time(elapsed_seconds=t1 - t0)

    print(time_str, "initial test time")

    # start
    print("\nDo .....")

    for i, row in datasets_parameters_df.iterrows():
        SAMPLE = row["dataset_name"]
        min_pathway_size = row["min_pathway_size"]
        k = 8 # row["k_for_kmeans"]  # todo: correct this

        moran_cutoff_features = 0.2

        print("\nSample: ", SAMPLE)
        subfolder = SAMPLE
        plot_folder = os.path.join(plot_parent_folder, SAMPLE)

        mydata = ad.read_h5ad(filename=os.path.join(
            mydata_dir, subfolder, f'{SAMPLE}.h5ad'))
        print(mydata)

        print(mydata.uns["history"])
        print(
            "MANDATORY: take away bad quality m/z")  # todo: how to make an alert ?
        mydata = mydata[:,
                 mydata.var.loc[mydata.var['final_quality'] == "good",
                 "mz"].astype(str).tolist()].copy()

        top_moran_i = mydata.var.loc[
                mydata.var['I'] >= moran_cutoff_features, "mz"].astype(
                str).tolist()
        mydata = mydata[:, top_moran_i].copy()

        print("vars of the data for clustering:", mydata.var.shape[0])


        # X = np.array(mydata.X.todense())
        # X = np.log1p(X)  # ok log the initial X

        X = np.array(mydata.X.todense())
        X = np.log1p(X)  # ok log the initial X
        t0 = time.time()
        # here the action
        kmeans = KMeans(n_clusters=k, n_init=10, random_state=0).fit(X)
        t1 = time.time()
        print("direct kmeans", give_time(elapsed_seconds=t1-t0))

        the_labels_str_l = kmeans.labels_.astype(str)
        mydata.obs["kmeans_direct"] = the_labels_str_l
        plot_out_path = os.path.join(plot_folder,
                                     f'{SAMPLE}-directkmeans-{k}.png')
        sq.pl.spatial_scatter(
            mydata, color=["kmeans_direct"], img=False,
            title='All m/z (good quality)',
            save=os.path.abspath(plot_out_path), dpi=600
        )

        ##
        # new log
        sc.pp.log1p(mydata)
        sc.tl.pca(mydata)

        t0 = time.time()
        # here the action
        kmeans = KMeans(n_clusters=k, random_state=0).fit(
            mydata.obsm['X_pca'])  # fix the random state for reproducibility
        t1 = time.time()
        print("pca + kmeans", give_time(elapsed_seconds=t1 - t0))

        the_labels_str_l = kmeans.labels_.astype(str)
        mydata.obs["kmeans_pca"] = the_labels_str_l
        plot_out_path = os.path.join(plot_folder,
                                     f'{SAMPLE}-pca-kmeans-{k}.png')
        sq.pl.spatial_scatter(
            mydata, color=["kmeans_pca"], img=False,
            title='All m/z (good quality): pca + kmeans',
            save=os.path.abspath(plot_out_path), dpi=600
        )

        ##  time computing neighbors and umap, needed for the following steps
        t0 = time.time()
        sc.pp.neighbors(mydata, metric="cosine")
        sc.tl.umap(mydata)
        t1 = time.time()
        print("umap and neighbors computation", give_time(t1- t0))

        ##
        t0 = time.time()
        # here the action

        kmeans = KMeans(n_clusters=k, random_state=0).fit(
            mydata.obsm['X_umap'])
        t1 = time.time()
        print("umap + kmeans", give_time(elapsed_seconds=t1 - t0) )

        the_labels_str_l = kmeans.labels_.astype(str)
        mydata.obs["kmeans_umap"] = the_labels_str_l
        plot_out_path = os.path.join(plot_folder,
                                     f'{SAMPLE}-umap-kmeans-{k}.png')
        sq.pl.spatial_scatter(
            mydata, color=["kmeans_umap"], img=False,
            title='All m/z (good quality): umap + kmeans',
            save=os.path.abspath(plot_out_path), dpi=600
        )
        ##
        my_resolution_leiden = 0.7 # 0.4
        my_resolution_leiden_spa = 0.7  # 0.7

        ##
        t0 = time.time()
        # here the action, recompute umap to take the time whole step
        sc.tl.leiden(mydata,
                     directed=False, key_added="leiden",
                     flavor="igraph",
                     n_iterations=2,
                     resolution=my_resolution_leiden)
        t1 = time.time()
        print("leiden", give_time(elapsed_seconds=t1 - t0))

        plot_out_path = os.path.join(plot_folder,
                                     f'{SAMPLE}-leiden-{my_resolution_leiden}.png')
        sq.pl.spatial_scatter(
            mydata, color=["leiden"], img=False,
            title='All m/z (good quality): leiden',
            save=os.path.abspath(plot_out_path), dpi=600
        )

        ##

        t0 = time.time()
        sc.pp.neighbors(mydata)
        nn_graph_genes = mydata.obsp["connectivities"]
        sq.gr.spatial_neighbors(mydata)
        nn_graph_space = mydata.obsp["spatial_connectivities"]

        alpha = 0.35  # weight spatial graph

        joint_graph = (1 - alpha) * nn_graph_genes + alpha * nn_graph_space
        sc.tl.leiden(mydata,
                     adjacency=joint_graph,
                     n_iterations=2,
                     flavor="igraph",
                     directed=False,
                     key_added="leiden_w_spatial_gr",
                     resolution=my_resolution_leiden_spa)
        t1 = time.time()
        print("leiden with weighted spatial graph joined", give_time(elapsed_seconds=t1 - t0))
        plot_out_path = os.path.join(plot_folder,
                                     f'{SAMPLE}-leiden_spa-{my_resolution_leiden_spa}.png')
        sq.pl.spatial_scatter(
            mydata, color=["leiden_w_spatial_gr"], img=False,
            title='All m/z (good quality): leiden with spatial graph',
            save=os.path.abspath(plot_out_path), dpi=600
        )
















