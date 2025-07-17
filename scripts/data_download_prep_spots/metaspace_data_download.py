import sys
import os
import pandas as pd
from spots_purge import show_message
import smpath.anndata_aux.matrix_from_metaspace as mfm
import warnings
""""
download any dataset using info in a .csv table
"""

warnings.simplefilter('ignore', category=FutureWarning)

# constants
def get_whole_api_versions():
    whole_api_versions = {'HMDB': ("HMDB", "v4"),
                    'LipidMaps': ("LipidMaps", "2017-12-12"),
                    'LIPID_MAPS' : ("LIPID_MAPS", "2016"),
                    'KEGG': ("KEGG", "v1"),
                    'CoreMetabolome': ('CoreMetabolome', 'v3')
                    }
    return whole_api_versions


# start
def run_data_download(datasets_parameters_df, out_dir)-> None:
    api_fdr = 0.5
    for i, r in datasets_parameters_df.iterrows():
        SAMPLE = r["dataset_name"]
        metaspace_registry = r["metaspace_registry"]
        annotations_to_use = r["annotations_to_use"].split(" ")
        subfolder = SAMPLE

        show_message(f"\n {SAMPLE}")

        path_out = os.path.join(out_dir, subfolder)
        expected_out_file = os.path.join(path_out, f'{SAMPLE}.h5ad')
        whole_api_versions = get_whole_api_versions()

        if not os.path.exists(path_out):
            os.makedirs(path_out)

        if not os.path.exists(expected_out_file):

            assert len(set(annotations_to_use).intersection(
                set(list(whole_api_versions.keys()))
            )) == len(set(annotations_to_use)), \
                "will fail: one or more annotations_to_use is not supported"

            api_versions = {k: whole_api_versions[k] for k in annotations_to_use}

            out_df, annotations_dfs_dict = mfm.get_unified_intensities_df(
                metaspace_registry, api_versions, api_fdr)
            del annotations_dfs_dict
            mfm.save_quantifications_anndata(
                out_df, os.path.join(path_out, SAMPLE))  # adds automatically .h5ad extension

            show_message(f"Saved: {expected_out_file}")
            del out_df
        else:
            show_message(f"{expected_out_file} exists, nothing to do")


if __name__ == "__main__":
    # define here my files and dirs
    #my_csv_table_path = "../../datasets_parameters/example-datasets-parameters.tsv"
    #out_dir = "../../../data/"
    my_csv_table_path = os.path.join(sys.argv[2], sys.argv[1])
    out_dir = sys.argv[3]

    datasets_parameters_df = pd.read_csv(my_csv_table_path, header=0,
                                         index_col=None, sep='\t')

    run_data_download(datasets_parameters_df, out_dir)