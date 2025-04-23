import os
from operator import index
import pandas as pd

import zonist.aux.matrix_from_metaspace as mfm

""""
download any dataset using info in a .csv table
"""

# define here my files and dirs
my_csv_table_path = "../../public-datasets-parameters.tsv"   #"../../small_test.csv"


datasets_parameters_df = pd.read_csv(my_csv_table_path, header=0,
                  index_col=None, sep='\t')

out_dir = "../../../data/"

# constants
whole_api_versions = {'HMDB': ("HMDB", "v4"),
                'LipidMaps': ("LipidMaps", "2017-12-12"),
                'LIPID_MAPS' : ("LIPID_MAPS", "2016"),
                'KEGG': ("KEGG", "v1"),
                'CoreMetabolome': ('CoreMetabolome', 'v3')
                }
api_fdr = 0.5

# start
for i, r in datasets_parameters_df.iterrows():
    data_name = r["dataset_name"]
    metaspace_registry = r["metaspace_registry"]
    annotations_to_use = r["annotations_to_use"].split(" ")
    annot_how_dwld = r["annotations_download"]
    subfolder = data_name

    print("\n", data_name)

    path_out = os.path.join(out_dir, subfolder)
    expected_out_file = os.path.join(path_out, f'{data_name}.h5ad')

    if not os.path.exists(path_out):
        os.makedirs(path_out)

    if not os.path.exists(expected_out_file):

        if not len(set(annotations_to_use).intersection(
            set(list(whole_api_versions.keys()))
        )        ) == len(set(annotations_to_use)):
            print("will fail: one or more annotations_to_use is not supported")

        api_versions = {k: whole_api_versions[k] for k in annotations_to_use}
        print("    will fetch quantification matrices of:", list(api_versions.keys()))

        print(f"Processing: {metaspace_registry}")

        out_df, annotations_dfs_dict = mfm.steps_get_united_df(
            metaspace_registry, api_versions, api_fdr)

        # new : save annotations on the fly:
        if annot_how_dwld == "automatic":
            for k in api_versions.keys():
                annotations_dfs_dict[k].to_csv(
                    os.path.join(path_out,  f'{data_name}-{k}.csv'),
                    sep=",", header=True, index=False
                )
            print("saved annotations automatically")

        del annotations_dfs_dict

        mfm.save_quantifications_anndata(
            out_df, os.path.join(path_out, data_name))  # adds automatically .h5ad extension

        print(f"Saved: {expected_out_file}")
        del out_df
    else:
        print(f"{expected_out_file} exists, nothing to do")









