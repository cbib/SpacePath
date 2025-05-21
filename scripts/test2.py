import sys
import os
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
import seaborn as sns
import squidpy as sq
import zonist.data as zd
import zonist.method as zm
import zonist.processing.measures_methods as zmm
import zonist.helpers as zh
from sklearn.cluster import KMeans

from local_helpers.local_helpers import provide_annotation_dict


""""
"""
if __name__ == "__main__":

   # define here my files and dirs

    mydata_dir = "../../data"
    plot_parent_folder = "../../figures/"

    my_csv_table_path_1 = "../subset1.tsv"  # "../../small_test.csv"

    df2 = pd.read_csv(my_csv_table_path_1, header=0,
                      index_col=None, sep='\t')
    datasets_parameters_df = pd.read_csv(my_csv_table_path_1, header=0,
                      index_col=None, sep='\t')

    # start
    dict_of_k = { 'NO': 5,
                  'Scores':5}

    k_mz_initial = 5

    print("\nDo .....")

    for i, row in datasets_parameters_df.iterrows():
        SAMPLE = row["dataset_name"]
        min_pathway_size = row["min_pathway_size"]
        print("\nSample: ", SAMPLE)
        subfolder = SAMPLE
        plot_folder = os.path.join(plot_parent_folder, SAMPLE)
        annotations_to_use = row["annotations_to_use"].split(" ")
        anno_dict = provide_annotation_dict(SAMPLE, annotations_to_use)

        mydata = ad.read_h5ad(filename=os.path.join(
            mydata_dir, subfolder, f'{SAMPLE}.h5ad'))

        print(mydata.uns["history"])
        print(
            "MANDATORY: take away bad quality m/z")  # todo: how to make an alert ?
        mydata = mydata[:,
                 mydata.var.loc[mydata.var['final_quality'] == "good",
                 "mz"].astype(str).tolist()].copy()

        X = np.array(mydata.X.todense())
        X = np.log1p(X)  # do not log  the pathway scores, ok
        k = k_mz_initial
        kmeans = KMeans(n_clusters=k, n_init=10, random_state=0).fit(X)
        the_labels_str_l = kmeans.labels_.astype(str)

        pd.DataFrame(
            {'obs': mydata.obs['obs'].tolist(),
             f'kmeans_{k}': the_labels_str_l}
        ).to_csv(
            os.path.join(
                mydata_dir, subfolder,
                f'{SAMPLE}-clusters-initial.txt'), index=False)

        mydata.obs[f'kmeans_{k}'] = the_labels_str_l
        plot_out_path = os.path.join(plot_folder,
                                     f'{SAMPLE}-Initial-k-{k}')
        sq.pl.spatial_scatter(
            mydata, color=[f'kmeans_{k}'], img=False,
            title='Initial', save=os.path.abspath(plot_out_path))
        plt.close()
        del X
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

                adata.obs[f'kmeans_{k}'] = the_labels_str_l
                plot_out_path = os.path.join(plot_folder,
                                             f'{SAMPLE}-{pwtype}-{SUFFIX}-k-{k}')
                sq.pl.spatial_scatter(
                    adata, color=[f'kmeans_{k}'], img=False,
                    title=f'{pwtype} {SUFFIX} k={k}',
                    save=os.path.abspath(plot_out_path))
                plt.close()

                if SUFFIX == 'Scores':
                    k -= 1
                    while k >= 3:
                        kmeans = KMeans(n_clusters=k, n_init=10,
                                        random_state=0).fit(X)
                        the_labels_str_l = kmeans.labels_.astype(str)

                        pd.DataFrame(
                            {'obs': adata.obs['obs'].tolist(),
                             f'kmeans_{k}': the_labels_str_l}
                        ).to_csv(
                            os.path.join(
                                mydata_dir, subfolder,
                                f'{SAMPLE}-{pwtype}-clusters-{SUFFIX}-k-{k}.txt'),
                            index=False)

                        adata.obs[f'kmeans_{k}'] = the_labels_str_l
                        plot_out_path = os.path.join(plot_folder,
                                                     f'{SAMPLE}-{pwtype}-{SUFFIX}-k-{k}')
                        sq.pl.spatial_scatter(
                            adata, color=[f'kmeans_{k}'], img=False,
                            title=f'{pwtype} {SUFFIX} k={k}',
                            save=os.path.abspath(plot_out_path))
                        plt.close()
                        k -= 1
                    # end while
                # end if
            # end for SUFFIX
            del adata
        # end for pwtype