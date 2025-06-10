import sys
import os
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
import squidpy as sq
from sklearn.cluster import KMeans

from local_helpers.local_helpers import provide_annotation_dict




if __name__ == "__main__":
    mydata_dir = "../../data"
    plot_parent_folder = "../../figures/"

    my_csv_table_path_1 = "../public-datasets-parameters.tsv"  # "../../small_test.csv"

    datasets_parameters_df = pd.read_csv(my_csv_table_path_1, header=0,
                      index_col=None, sep='\t')

    # start
    dict_of_k = { 'NO': 5,
                  'Scores':5}
    k_for_original_dataset = 5
    print("\nDo .....")

    #for i, row in datasets_parameters_df.iterrows():
    for no_row in range(1):
        SAMPLE = "TD-Lipido-635"  #row["dataset_name"]
        min_pathway_size = 2   #row["min_pathway_size"]
        print("\nSample: ", SAMPLE)
        subfolder = SAMPLE
        plot_folder = os.path.join(plot_parent_folder, SAMPLE)
        annotations_to_use = ["LipidMaps", "HMDB", "CoreMetabolome"]  #row["annotations_to_use"].split(" ")
        anno_dict = provide_annotation_dict(SAMPLE, annotations_to_use)   # todo: check if never used here (also see test2_kmeans)

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
        k = k_for_original_dataset

        # if not os.path.exists(os.path.join(
        #         mydata_dir, subfolder,
        #         f'{SAMPLE}-clusters-initial.txt')):
        #     kmeans = KMeans(n_clusters=k, n_init=10, random_state=0).fit(X)
        #     the_labels_str_l = kmeans.labels_.astype(str)
        #
        #     pd.DataFrame(
        #         {'obs': mydata.obs['obs'].tolist(),
        #          f'kmeans_{k}': the_labels_str_l}
        #     ).to_csv(
        #         os.path.join(
        #             mydata_dir, subfolder,
        #             f'{SAMPLE}-clusters-initial.txt'), index=False)
        # else:
        #     print(f" for {SAMPLE}, initial kmeans exists")
        #     tdf = pd.read_csv(os.path.join(
        #             mydata_dir, subfolder,
        #             f'{SAMPLE}-clusters-initial.txt'), index_col=None)
        #     the_labels_str_l = tdf[f'kmeans_{k}'].astype(str).tolist()
        #     del tdf
        #
        #
        # mydata.obs[f'kmeans_{k}'] = the_labels_str_l
        # plot_out_path = os.path.join(plot_folder,
        #                              f'{SAMPLE}-Initial-k-{k}')
        # sq.pl.spatial_scatter(
        #     mydata, color=[f'kmeans_{k}'], img=False,
        #     title='Initial', save=os.path.abspath(plot_out_path), dpi=600)
        # plt.close()
        del X
        for pwtype in ['TANIMOTO']:
            adata = mydata.copy()
            for SUFFIX in ['Scores']:  # strict order to respect
                if SUFFIX == 'NO':  #  todo: stop saving the NO part, because it corresponds to  :  yes pathways OR not spatially differential
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
                    X = np.log1p(X)   # only log iif Not pathway scores, ok

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

            del adata
            break  # todo : this is violent

    print("end")