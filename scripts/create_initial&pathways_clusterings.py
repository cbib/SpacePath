import os
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import squidpy as sq
from sklearn.cluster import KMeans

# TODO: solve problem with k definition (file not taken into account ????)

def save_clusters_constant_file_rootname_distinct_k(obs_df, k, mydata_dir, subfolder, SAMPLE):
    df = obs_df[['obs', f'kmeans_{k}']]
    df.to_csv(os.path.join(mydata_dir, subfolder, f'{SAMPLE}-clusters-k_{k}.tsv'),
              sep='\t', header=True, index=False)
    print("saved clusters as tsv file :)")

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
    dict_of_k = { 'NO': 5,   'Scores': 5}
    k_for_original_dataset = 5
    print("\nDo .....")

    for i, row in datasets_parameters_df.iterrows():
        SAMPLE = row["dataset_name"]
        min_pathway_size = row["min_pathway_size"]
        k = row["k_for_kmeans"]
        if k not in [2, 30]:
            k = k_for_original_dataset
        print("\nSample: ", SAMPLE)
        subfolder = SAMPLE
        plot_folder = os.path.join(plot_parent_folder, SAMPLE)

        mydata = ad.read_h5ad(filename=os.path.join(
            mydata_dir, subfolder, f'{SAMPLE}.h5ad'))

        print(mydata.uns["history"])
        print(
            "MANDATORY: take away bad quality m/z")  # todo: how to make an alert ?
        mydata = mydata[:,
                 mydata.var.loc[mydata.var['final_quality'] == "good",
                 "mz"].astype(str).tolist()].copy()

        if not os.path.exists(os.path.join(
                mydata_dir, subfolder,
                f'{SAMPLE}-clusters-initial.txt')):

            data_sp_valid = mydata[:, mydata.var.loc[mydata.var['I'] >= 0.2,
                "mz"].astype(str).tolist()].copy()

            X = np.array(data_sp_valid.X.todense())
            X = np.log1p(X)  # ok log the initial X
            k = k_for_original_dataset

            kmeans = KMeans(n_clusters=k, n_init=10, random_state=0).fit(X)
            the_labels_str_l = kmeans.labels_.astype(str)

            pd.DataFrame(
                {'obs': mydata.obs['obs'].tolist(),
                 f'kmeans_{k}': the_labels_str_l}
            ).to_csv(
                os.path.join(
                    mydata_dir, subfolder,
                    f'{SAMPLE}-clusters-initial.txt'), index=False)
        else:
            print(f" for {SAMPLE}, initial kmeans exists")
            tdf = pd.read_csv(os.path.join(
                    mydata_dir, subfolder,
                    f'{SAMPLE}-clusters-initial.txt'), index_col=None)
            try:
                the_labels_str_l = tdf[f'kmeans_{k}'].astype(str).tolist()
            except:
                the_labels_str_l = tdf[f'kmeans_5'].astype(str).tolist()
            del tdf


        mydata.obs[f'kmeans_{k}'] = the_labels_str_l
        plot_out_path = os.path.join(plot_folder,
                                     f'{SAMPLE}-Initial-k-{k}.png')
        sq.pl.spatial_scatter(
            mydata, color=[f'kmeans_{k}'], img=False,
            title='Initial', save=os.path.abspath(plot_out_path), dpi=600)
        plt.close()

        for pwtype in ['WIKI', 'KEGG']:
            adata = mydata.copy()
            for SUFFIX in ['NO', 'Scores']:  # strict order to respect
                if SUFFIX == 'NO':
                    no_pws_mzs = pd.read_csv(
                        os.path.join(mydata_dir, subfolder,
                                     f'{SAMPLE}-{pwtype}-{SUFFIX}.txt'),
                        index_col=None)
                    adata = adata[:, no_pws_mzs['mz'].astype(str).tolist()].copy()
                elif SUFFIX == 'Scores':
                    adata = ad.read_h5ad(
                os.path.join(mydata_dir, subfolder, f'{SAMPLE}-scores-{pwtype}.h5ad'))

                X = np.array(adata.X.todense())

                if SUFFIX != "Scores":
                    X = np.log1p(X)   # do not log  the pathway scores, ok

                k = dict_of_k[SUFFIX]

                if not os.path.exists(os.path.join(
                            mydata_dir, subfolder,
                            f'{SAMPLE}-{pwtype}-clusters-{SUFFIX}-k-{k}.txt')):

                    kmeans = KMeans(n_clusters=k, n_init=10, random_state=0 ).fit(X)
                    the_labels_str_l =  kmeans.labels_.astype(str)

                    pd.DataFrame(
                        {'obs': adata.obs['obs'].tolist(),
                         f'kmeans_{k}': the_labels_str_l}
                    ).to_csv(
                        os.path.join(
                            mydata_dir, subfolder,
                            f'{SAMPLE}-{pwtype}-clusters-{SUFFIX}-k-{k}.txt'),
                        index=False)
                else:
                    print(f" for {SAMPLE}: {pwtype}: {SUFFIX}: {k}, kmeans exists")
                    tdf = pd.read_csv(os.path.join(
                        mydata_dir, subfolder,
                        f'{SAMPLE}-{pwtype}-clusters-{SUFFIX}-k-{k}.txt'), index_col=None)
                    the_labels_str_l = tdf[f'kmeans_{k}'].astype(str).tolist()
                    del tdf

                adata.obs[f'kmeans_{k}'] = the_labels_str_l
                plot_out_path = os.path.join(plot_folder,
                                             f'{SAMPLE}-{pwtype}-{SUFFIX}-k-{k}.png')
                sq.pl.spatial_scatter(
                    adata, color=[f'kmeans_{k}'], img=False,
                    title=f'{pwtype} {SUFFIX} k={k}',
                    save=os.path.abspath(plot_out_path),
                    dpi=600
                )
                plt.close()

                # save_clusters_constant_file_rootname_distinct_k(adata.obs.copy(), k, mydata_dir, subfolder, SAMPLE)

            # end for SUFFIX
            del adata
        # end for pwtype