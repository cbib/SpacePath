# johaGL 2025
# script for downloading tissue slices of Spatial Metabolomics from METASPACE
# usage example in `data_dps.sh`

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

        if (metaspace_registry != "not provided"
        ) and not os.path.exists(expected_out_file):

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
            show_message(f"Data already present!\n >> {expected_out_file}")


if __name__ == "__main__":
    # define here my files and dirs
    #input_file = "../../datasets_parameters/example-datasets-parameters.tsv"
    #data_dir = "../../../data/"
    input_file = sys.argv[1]
    data_dir = sys.argv[2]

    datasets_parameters_df = pd.read_csv(input_file, header=0,
                                         index_col=None, sep='\t')

    run_data_download(datasets_parameters_df, data_dir)
