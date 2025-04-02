import os
import pandas as pd
import anndata as ad


mydata_dir = "../../data"

# define here my input
my_csv_table_path = "../datasets-parameters.tsv"   #"../../small_test.csv"


datasets_parameters_df = pd.read_csv(my_csv_table_path, header=0,
                  index_col=None, sep='\t')

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
    print("sample : ", SAMPLE)
    name_l.append(SAMPLE)
    spot_n.append(adata_mz.obs.shape[0])
    features_n.append(adata_mz.var.shape[0])
    print("")

out_df = pd.DataFrame({'dataset_name': name_l,
                       'spot_number': spot_n,
                       'features_n' : features_n })

out_df.to_csv("../initial_data_dimensions.tsv", index=False, sep='\t')