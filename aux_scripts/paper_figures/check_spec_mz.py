import os

import anndata as ad
import pandas as pd
import squidpy as sq
import matplotlib.pyplot as plt
import seaborn as sns


if __name__ == "__main__":

    mysample = "LP-brain-mouse-5"
    mymz = "369.351578"

    # define here my files and dirs
    data_centralized_folder = "../../../data"

    my_csv_table_path_1 = "../../datasets_parameters/public-datasets-parameters.tsv"  #
    my_csv_table_path_2 = "../../datasets_parameters/private-datasets-parameters.tsv"  #

    df1 = pd.read_csv(my_csv_table_path_1, header=0,
                      index_col=None, sep='\t')
    df2 = pd.read_csv(my_csv_table_path_2, header=0,
                      index_col=None, sep='\t')
    datasets_parameters_df = pd.concat([df1, df2])

    datasets_parameters_df = datasets_parameters_df.loc[
        datasets_parameters_df["dataset_name"] == mysample,
    :].copy()

    for i, row in datasets_parameters_df.iterrows():
        SAMPLE = row["dataset_name"]
        min_pathway_size = row["min_pathway_size"]
        #print("\nSample: ", SAMPLE)
        subfolder = SAMPLE
        #try:
        try:
            adata_mz = ad.read_h5ad(filename=os.path.join(
                data_centralized_folder, subfolder, f'{SAMPLE}.h5ad'))
            # adata_mz = adata_mz[:,
            # adata_mz.var.loc[adata_mz.var['final_quality'] == "good",
            # "mz"].astype(str).tolist()]
            print("attention, bad and good quality mzs are present in the object !")

            print(adata_mz.uns['norm'])
            sns.set_context("paper", font_scale=1.5)

            if adata_mz.uns['norm'] == "TotalIonCount":
                sq.pl.spatial_scatter(
                    adata_mz, color=[str(mymz)], img=False,
                    title=f'$m/z$: {mymz}',
                    linewidth=0, edgecolors=(0, 0, 0, 0),
                    shape="circle", size=1.2,
                    palette='viridis')
                plt.show()
            else:
                print("Error: only normalized data accepted ! ")
        except FileNotFoundError:
            print(SAMPLE, ": could not compute top differential SETS, insufficient data quality or annotations")



