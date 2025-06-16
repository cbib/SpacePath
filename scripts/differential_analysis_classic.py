import os
import pandas as pd
import numpy as np
import zonist.data as zd
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import squidpy as sq


""""
"""
if __name__ == "__main__":

   # define here my files and dirs

    mydata_dir = "../../data"
    plot_parent_folder = "../../figures/"

    my_csv_table_path_1 = "../public-datasets-parameters.tsv"  # "../../small_test.csv"
    my_csv_table_path_2 = "../private-datasets-parameters.tsv"  # "../../small_test.csv"

    df1 = pd.read_csv(my_csv_table_path_1, header=0,
                      index_col=None, sep='\t')
    df2 = pd.read_csv(my_csv_table_path_2, header=0,
                      index_col=None, sep='\t')
    datasets_parameters_df = pd.concat([df1, df2])

    # start
    print("\nDo .....")

    for i, row in datasets_parameters_df.iterrows():
        SAMPLE = row["dataset_name"]
        print("\nSample: ", SAMPLE)
        subfolder = SAMPLE
        plot_folder = os.path.join(plot_parent_folder, SAMPLE)

        with open(os.path.join(mydata_dir, subfolder,
                f"{SAMPLE}-clusters-initial.txt"), "r") as f:
            initial_clusters = f.readlines()

        initial_clusters.pop(0) # the header, away
        initial_clusters = [float(t.strip().split(",")[1]) for t in initial_clusters]

        adata_mz = ad.read_h5ad(filename=os.path.join(
                    mydata_dir, subfolder, f'{SAMPLE}.h5ad'))

        print(adata_mz.uns["history"])
        print("MANDATORY: take away bad quality m/z") # todo: how to make an alert ?
        adata_mz = adata_mz[:,
                   adata_mz.var.loc[adata_mz.var['final_quality'] == "good",
                   "mz"].astype(str).tolist()]

        print(adata_mz)
        adata_mz.obs['kmeans_5'] = np.array(initial_clusters, dtype=str)

        sq.pl.spatial_scatter(adata_mz, color=[f'kmeans_5'],
                                     img=False, title="pop up (verification)")
        #plt.show()
        #plt.close()

        sc.pp.log1p(adata_mz)
        print("computing the differential abundances of all the m/z features, cluster to cluster, for all clusters")
        sc.tl.rank_genes_groups(adata_mz, groupby='kmeans_5', groups='all',
                                    reference='rest',
                                    method='wilcoxon' ,  # use logreg if ugly
                                    key_added="top_feats_by_cluster",
                                    corr_method='benjamini-hochberg')

        differential_results_by_cluster = sc.get.rank_genes_groups_df(adata_mz,
                                                 key="top_feats_by_cluster",
                                                 group=None )  # returns all groups)

        differential_results_by_cluster.to_csv(
            os.path.join(
                mydata_dir, subfolder, f"{SAMPLE}-differential-results.csv"),
            index=False)