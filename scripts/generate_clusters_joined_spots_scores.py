import os
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import squidpy as sq
from sklearn.cluster import KMeans

# TODO: solve problem with k definition (file not taken into account ????)

def cluster_families_alone(adata, k, mydata_dir, subfolder, SAMPLE):
    X = np.array(adata.X.todense())

    kmeans = KMeans(n_clusters=k, n_init=10, random_state=0).fit(X)
    the_labels_str_l = kmeans.labels_.astype(str)

    adata.obs[f'kmeans_{k}']   =  the_labels_str_l

    plot_out_path = os.path.join(plot_folder,
                                 f'{SAMPLE}-clusters-FamiliesScores-k-{k}.png')

    sq.pl.spatial_scatter(
        adata, color=[f'kmeans_{k}'], img=False,
        title=f'Families alone k={k}',
        # save=os.path.abspath(plot_out_path),
        # dpi=600
    )
    plt.show()

    df = adata.obs[['obs', f'kmeans_{k}']].copy()
    df.to_csv(os.path.join(mydata_dir, subfolder, f'{SAMPLE}-clusters-FamiliesScores-k_{k}.txt'),
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

    # my_csv_table_path_1 = "../small-test-parameters.tsv"  # "../../small_test.csv"
    #
    # datasets_parameters_df = pd.read_csv(my_csv_table_path_1, header=0,
    #                                      index_col=None, sep='\t')

# start
    print("\nDo .....")

    for i, row in datasets_parameters_df.iterrows():
        SAMPLE = row["dataset_name"]
        min_pathway_size = row["min_pathway_size"]
        k = row["k_for_kmeans"]
        if k not in [2, 30]:
            k = 5

        print("\nSample: ", SAMPLE)
        subfolder = SAMPLE
        plot_folder = os.path.join(plot_parent_folder, SAMPLE)

        o_fam = os.path.join(mydata_dir, subfolder,
                             f'{SAMPLE}-scores-Families.h5ad')
        o_kegg = os.path.join(mydata_dir, subfolder,
                             f'{SAMPLE}-scores-KEGG.h5ad')
        dafam = ad.read_h5ad(o_fam)

        cluster_families_alone(dafam,  k, mydata_dir, subfolder, SAMPLE)


        dapw = ad.read_h5ad(o_kegg)
        mydata = ad.concat((dapw, dafam), join='outer', axis=1)
        mydata.obs['obs'] = mydata.obs.index.tolist()

        # for plot purposes
        mydata.obsm = dafam.obsm
        mydata.uns['spatial'] = dafam.uns['spatial']

        print(mydata.var.head())
        print(mydata.obs.head())

        print(mydata)

        yep = True
        if yep :
        #if not os.path.exists(os.path.join( mydata_dir, subfolder,f'{SAMPLE}-clusters-initial.txt')):
            X = np.array(mydata.X.todense())

            kmeans = KMeans(n_clusters=k, n_init=10, random_state=0).fit(X)
            the_labels_str_l = kmeans.labels_.astype(str)

            pd.DataFrame(
                {'obs': mydata.obs['obs'].tolist(),
                 f'kmeans_{k}': the_labels_str_l}
            ).to_csv(
                os.path.join(
                    mydata_dir, subfolder,
                    f'{SAMPLE}-clusters-CombiFamilPathScores.txt'), index=False)
        else:
            tdf = pd.read_csv(os.path.join(
                    mydata_dir, subfolder,
                    f'{SAMPLE}-clusters-CombiFamilPathScores.txt'), index_col=None)
            the_labels_str_l = tdf[f'kmeans_{k}'].astype(str).tolist()
            del tdf

        mydata.obs[f'kmeans_{k}'] = the_labels_str_l
        plot_out_path = os.path.join(plot_folder,
                                     f'{SAMPLE}-clusters-CombiFamilPathScores-k-{k}.png')

        sq.pl.spatial_scatter(
            mydata, color=[f'kmeans_{k}'], img=False,
            title=f'Families & pathways combined k={k}',
            #save=os.path.abspath(plot_out_path),
            #dpi=600
        )
        plt.show()
        plt.close()


