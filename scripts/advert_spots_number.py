import os
import pandas as pd
import anndata as ad


mydata_dir = "../../data"

# define here my input
my_csv_table_path_1 = "../public-datasets-parameters.tsv"   #"../../small_test.csv"
my_csv_table_path_2 = "../private-datasets-parameters.tsv"   #"../../small_test.csv"


df1 = pd.read_csv(my_csv_table_path_1, header=0,
                  index_col=None, sep='\t')
df2 = pd.read_csv(my_csv_table_path_2, header=0,
                  index_col=None, sep='\t')
datasets_parameters_df = pd.concat([df1,df2])

name_l = []
spot_n = []
features_n = []
for i, r in datasets_parameters_df.iterrows():
    SAMPLE = r["dataset_name"]
    #annotations_to_use = r["annotations_to_use"].split(" ")
    #annot_how_dwld = r["annotations_download"]
    subfolder = SAMPLE
    adata_mz = ad.read_h5ad(filename=os.path.join(
        mydata_dir, subfolder, f'{SAMPLE}.h5ad'))
    try:
        nspot = adata_mz.uns['history']['spot_number_initial']
    except:
        nspot = adata_mz.obs.shape[0]

    try:
        nfeat = adata_mz.uns['history']['features_n_initial']
    except:
        nfeat = adata_mz.var.shape[0]

    print("sample : ", SAMPLE)
    name_l.append(SAMPLE)
    spot_n.append(nspot)
    features_n.append(nfeat)
    print("")

out_df = pd.DataFrame({'dataset_name': name_l,
                       'spot_number': spot_n,
                       'features_n' : features_n })

out_df.to_csv("../figs/initial_data_dimensions.tsv", index=False, sep='\t')