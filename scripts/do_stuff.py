import os
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import squidpy as sq
import numpy as np
import zonist.helpers as zh
from do_steps_picked_data import differential_results_filter
from do_steps_picked_data import wrapper_differential

# define here my files and dirs
mydata_dir = "../../data"
plot_parent_folder = "../../figures/"

picked_sample = "TD-Lipido-635"  # brain 5, TD and LP-brain-rat-2 ,

my_csv_table_path_1 = "../public-datasets-parameters.tsv"  # "../../small_test.csv"
my_csv_table_path_2 = "../private-datasets-parameters.tsv"  # "../../small_test.csv"

haha = zh.is_log_transformed(np.array([[4, -6, 8],[8,4, -3]]))
print(haha)

print(inlijnlikj)
df1 = pd.read_csv(my_csv_table_path_1, header=0,
                  index_col=None, sep='\t')
df2 = pd.read_csv(my_csv_table_path_2, header=0,
                  index_col=None, sep='\t')

datasets_parameters_df = pd.concat([df1, df2])
force_for_these = []
# start

datasets_parameters_df = datasets_parameters_df.loc[
    datasets_parameters_df['dataset_name'] == picked_sample
    ].copy()


print("\nfull pipeline my picked data")

for i, row in datasets_parameters_df.iterrows():
    SAMPLE = row["dataset_name"]
    normalize_yes_no = row["TIC_normalization"]
    offSample_cutoff = row["offSample_cutoff"]
    msm_cutoff = row["msm_cutoff"]
    fdr_cutoff = row["fdr_cutoff"]
    moran_I_cutoff = row["Moran_I_cutoff"]
    moran_padj_cutoff = row["Moran_padj_cutoff"]
    cutoff_tree = row['cut_tree_mol_families']

    print("\nSample: ", SAMPLE)

    subfolder = SAMPLE

    user_k = 8

    plot_folder = os.path.join(plot_parent_folder, SAMPLE)

    adata_mz = ad.read_h5ad(filename=os.path.join(
        mydata_dir, subfolder, f'{SAMPLE}.h5ad'))
    print("Loaded object", adata_mz.uns["history"])

    ##  1. classify mzs as good or bad quality: already done

    ## 2. Clean annotations that will be used:  already done by the clean annots to tables

    # for each patwhay type: kegg or wiki
    print(adata_mz.uns["history"])
    print(
        "MANDATORY: take away bad quality m/z")  # todo: how to make an alert ?

    adata_mz = adata_mz[:,
               adata_mz.var.loc[adata_mz.var['final_quality'] == "good",
               "mz"].astype(str).tolist()]

    cluss = pd.read_csv(
        os.path.join(mydata_dir, subfolder,
                     f'{SAMPLE}-REGRESSIVES-K--{user_k}.txt'), header=0)

    adata_mz.obs[f"kmeans{user_k}"] = cluss[f'kmeans{user_k}'].astype(str).tolist()
    #
    # sq.pl.spatial_scatter(
    #     adata_mz, color=[f'kmeans{user_k}'], img=False,
    #     title=f'k={user_k}',
    #     #save=os.path.abspath(plot_out_path),
    #     )
    # plt.show()


    diffres = wrapper_differential(adata_mz, k=user_k)

    diff_filt = differential_results_filter(differential_results=diffres, ntop_per_cluster=10)

    print(diff_filt.head())
    my_intere = diff_filt.loc[diff_filt['group']==str(2), 'names'].tolist()
    print(my_intere)
    for imz in my_intere:
        sq.pl.spatial_scatter(
                adata_mz, color=[str(imz)], img=False,
                 title=str(imz),
                 #save=os.path.abspath(plot_out_path),
                 )
        plt.show()

    leritos = 88
    yuyuy = 00

