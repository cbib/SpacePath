
import os
import pandas as pd
import anndata as ad


mydata_dir = "../../data"

print("Reset all objects so the features classification is removed")
print("  but keeps .X as TIC normalized data")

# define here my input
my_csv_table_path_1 = "../public-datasets-parameters.tsv"   #"../../small_test.csv"
my_csv_table_path_2 = "../private-datasets-parameters.tsv"   #"../../small_test.csv"


df1 = pd.read_csv(my_csv_table_path_1, header=0,
                  index_col=None, sep='\t')
df2 = pd.read_csv(my_csv_table_path_2, header=0,
                  index_col=None, sep='\t')
datasets_parameters_df = pd.concat([df1,df2])


for i, r in datasets_parameters_df.iterrows():
    SAMPLE = r["dataset_name"]

    subfolder = SAMPLE
    adata_mz = ad.read_h5ad(filename=os.path.join(
        mydata_dir, subfolder, f'{SAMPLE}.h5ad'))


    subfolder = SAMPLE

    adata_mz = ad.read_h5ad(filename=os.path.join(
           mydata_dir, subfolder, f'{SAMPLE}.h5ad'))
    try:
        del adata_mz.uns['history']['features_n_after_mz_clean_1']
    except:
        pass
    try:
        del adata_mz.uns['history']['features_n_after_mz_clean_2']
    except:
        pass
    try:
        del adata_mz.uns['moranI']
    except:
        pass
    if adata_mz.var.shape[1] > 2 :
        these_columns = adata_mz.var.columns.tolist()[1:]
        adata_mz.var = adata_mz.var.drop(columns=these_columns)

    adata_mz.write_h5ad(os.path.join(mydata_dir, subfolder, f'{SAMPLE}.h5ad'))

