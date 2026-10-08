# johaGL 2025
# script for cleaning spots in a tissue slice of Spatial Metabolomics
# usage example in `data_dps.sh`
#
# further usage examples:
# - running without forcing on specific tissues
#     python3 spots_purge.py "$input_directory/$input_file" ../../../data ../../../figures
# or
# - forcing two tissues to run again:
#     python3 spots_purge.py "$input_directory/$input_file" ../../../data ../../../figures "LP-brain-mouse-7;LP-kidney3"
import sys
import logging
import os
import pandas as pd
import anndata as ad
import json
from typing import List
import SpacePath.core.data as spd
import SpacePath.core.processing.measures_methods as measures_methods
import matplotlib.pyplot as plt

import squidpy as sq


logger = logging.getLogger(__name__)

def annotation_types_to_file_names_dict(SAMPLE, annotation_types_list):
    out_d = {}
    for annotation_type in annotation_types_list:
        out_d[annotation_type] = f"{SAMPLE}-{annotation_type}.csv"

    return out_d


def plot_spots(adata_mz, plot_out_path,
               title='my title') -> ad.AnnData:

    if 'spot_quality_colors' in adata_mz.uns.keys():
        # avoid ValueError: Expected palette to be of length `2`, found `1`
        adata_mz.uns['spot_quality_colors'] = None  # compulsory to unlock key
        del adata_mz.uns['spot_quality_colors']

    adata_mz.obs['spot_quality'] = pd.Categorical(
        adata_mz.obs['spot_quality'].tolist(), categories=['bad', 'good'])
    show_message("Saving plot of spots")
    sq.pl.spatial_scatter(
        adata_mz, color=["spot_quality"], img=False,
        title=title, save=os.path.abspath(plot_out_path),
        dpi=600)
    plt.close()

    return adata_mz


def show_message(message)-> None:
    logger.info(message)
    print(message)


def custom_out_tissue(adata_mz, data_centralized_folder, subfolder, SAMPLE):
    file_custom_list_good_spots = os.path.join(
        data_centralized_folder, subfolder, f'{SAMPLE}_out_tissue.json'
    )
    if os.path.exists(file_custom_list_good_spots):
        with open(file_custom_list_good_spots, "r") as f:
            mydict = json.load(f)

        bad_spots = mydict['out_tissue_spots']
        adata_mz.obs.loc[adata_mz.obs["obs"].isin(bad_spots), "spot_quality"] = "bad"
    else:
        raise FileNotFoundError

    return adata_mz

def run_final_spots_preparation(datasets_parameters_df, data_centralized_folder,
                                figures_centralized_folder, force_for_these):

    show_message(
        "\nRunning spots classification (with option to select or not) across specified datasets")

    for i, row in datasets_parameters_df.iterrows():
        SAMPLE = row["dataset_name"]
        allowed_in_out_ratio = row["in/out intensities for spot class"]
        select_spots_yes_no = row["perform_spots_selection"]
        show_message(f"\nSample: {SAMPLE}")
        subfolder = SAMPLE
        plot_folder = os.path.join(figures_centralized_folder, SAMPLE)
        if not os.path.exists(plot_folder):
            os.mkdir(plot_folder)

        adata_mz = ad.read_h5ad(filename=os.path.join(
            data_centralized_folder, subfolder, f'{SAMPLE}.h5ad'))

        if not 'history' in adata_mz.uns.keys():
            adata_mz.uns['history'] = {
                'spot_number_initial': adata_mz.obs.shape[0],
                'features_n_initial': adata_mz.var.shape[0]}

        if (SAMPLE in force_for_these) or (
                "spot_quality" not in adata_mz.obs.columns.tolist()) or (
                select_spots_yes_no.lower() in ["yes", "y"]
        ) :
            annotations_to_use = row["annotations_to_use"].split(" ")
            annotations_dict = annotation_types_to_file_names_dict(SAMPLE,
                                                       annotations_to_use)
            annots_obj = spd.MetabolomeAnnotations(name=SAMPLE,
                                                  data_dir=data_centralized_folder,
                                                  subfolder=subfolder,
                                                  annotation_files=annotations_dict)
            annots_obj.load()
            annots_obj, adata_mz = measures_methods.match_adata_mz_and_annots(annots_obj,
                                                                 adata_mz)

            annots_obj.pull_mz_quality_df()  ## NEW
            if not annots_obj.mz_quality[
                annots_obj.offsample_column].isnull().all():
                fdr_cutoff_features_to_compute_spot_qual = 0.1
                if SAMPLE in ["LP-brain-mouse-5"]:
                    fdr_cutoff_features_to_compute_spot_qual = 0.2
                x = fdr_cutoff_features_to_compute_spot_qual
                adata_mz = measures_methods.compute_spot_quality(adata_mz, annots_obj,
                                                    allowed_in_out_ratio,
                                                    fdr_cutoff=x)
                adata_mz = measures_methods.mark_weak_signal_spots(adata_mz, cutoff=100)

                try:
                    adata_mz = custom_out_tissue(
                        adata_mz, data_centralized_folder, subfolder, SAMPLE)
                except:
                    pass

                plot_out_path = os.path.join(plot_folder,
                                             f'{SAMPLE}-spotsBefore.png')
                if not os.path.exists(plot_out_path):
                    plot_spots(adata_mz, plot_out_path,
                               title="Before spots selection")
                    show_message(f"plot of initial spots saved to: {plot_out_path}")

                if select_spots_yes_no.lower() in ["yes", "y"]:
                    adata_mz = adata_mz[
                               adata_mz.obs["spot_quality"] == "good", :]
                    plot_out_path = os.path.join(plot_folder,
                                                 f'{SAMPLE}-spotsAfter.png')

                    plot_spots(adata_mz, plot_out_path,
                               title="After spots selection")
                    show_message(
                        "Discarded spots of bad quality, will overwrite new .h5ad")
                else:
                    show_message(
                        "Spot quality added to object, but original spots preserved")

                adata_mz.uns['history']['spot_number_retained'] = \
                adata_mz.obs.shape[0]
                adata_mz.write_h5ad(
                    os.path.join(data_centralized_folder, subfolder,
                                 f'{SAMPLE}.h5ad'))
                show_message("saved object to .h5ad file")

            else:
                show_message(
                    "Skipped!, the out-of-tissue probabilities are not available")
        else:
            show_message(
                "already computed quality of spots for adata_mz, nothing to do")


if __name__ == "__main__":
    # input_file = "../../datasets_parameters/example-datasets-parameters.tsv"  #  "private-datasets-parameters.tsv"
    # data_centralized_folder = "../../../data"
    # figures_centralized_folder = "../../../figures/"
    # my_csv_table_path = os.path.join(input_dir, input_file)

    my_csv_table_path = sys.argv[1]
    data_centralized_folder = sys.argv[2]
    figures_centralized_folder = sys.argv[3]
    force_for_these = []
    try:
        force_for_these = sys.argv[4].split(";")
    except IndexError:
        pass

    datasets_parameters_df = pd.read_csv(my_csv_table_path, header=0,
                                         index_col=None, sep='\t')
    run_final_spots_preparation(datasets_parameters_df, data_centralized_folder,
                               figures_centralized_folder,
                                force_for_these=force_for_these)

