import sys
import os
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import squidpy as sq
import networkx as nx
from tqdm import tqdm
import matplotlib.pyplot as plt
import json

import zonist.data as zd
import zonist.processing.measures_methods as zmm


from local_helpers.local_helpers import provide_annotation_dict

from local_helpers.local_viz import plot_spots




def custom_out_tissue(adata_mz, mydata_dir, subfolder, SAMPLE):
    file_custom_list_good_spots = os.path.join(
        mydata_dir, subfolder, f'{SAMPLE}_out_tissue.json'
    )
    if os.path.exists(file_custom_list_good_spots):
        with open(file_custom_list_good_spots, "r") as f:
            mydict = json.load(f)

        bad_spots = mydict['out_tissue_spots']
        adata_mz.obs.loc[adata_mz.obs["obs"].isin(bad_spots), "spot_quality"] = "bad"
    else:
        raise FileNotFoundError

    return adata_mz


if __name__ == "__main__":
    # define here my input
    mydata_dir = "../../data"
    my_csv_table_path = "../privateNEW-datasets-parameters.tsv"
    plot_parent_folder = "../../figures/"

    # start
    datasets_parameters_df = pd.read_csv(my_csv_table_path, header=0,
                      index_col=None, sep='\t')

    #  , "TD-Lipido-836", "TD-Lipido-827"
    force_for_these = [ "TD-Lipido-635", "TD-Lipido-636"] # leave empty if no samples to force here   # e.g. ["TD-hemi-635-1-A", "TD-hemi-635-1-B"]    "TD-635-1-B"

    print("\nRunning spots classification (with option to select or not) across specified datasets")

    for i, row in datasets_parameters_df.iterrows():
        SAMPLE = row["dataset_name"]
        allowed_in_out_ratio = row["in/out intensities for spot class"]
        select_spots_yes_no = row["perform_spots_selection"]
        print("\nSample: ", SAMPLE)
        subfolder = SAMPLE
        plot_folder = os.path.join(plot_parent_folder, SAMPLE)
        if not os.path.exists(plot_folder):
            os.mkdir(plot_folder)

        adata_mz = ad.read_h5ad(filename=os.path.join(
            mydata_dir, subfolder, f'{SAMPLE}.h5ad'))

        if not 'history' in adata_mz.uns.keys():
            adata_mz.uns['history'] = {'spot_number_initial': adata_mz.obs.shape[0],
                                   'features_n_initial': adata_mz.var.shape[0]}


        if (SAMPLE in force_for_these) or "spot_quality" not in adata_mz.obs.columns.tolist():
            annotations_to_use = row["annotations_to_use"].split(" ")
            anno_dict = provide_annotation_dict(SAMPLE, annotations_to_use)
            me_annots = zd.MetabolomeAnnotations(name=SAMPLE,
                                                 data_dir=mydata_dir,
                                                 subfolder=subfolder,
                                                 annotation_files=anno_dict)
            me_annots.load()  # load annotations
            print(adata_mz.uns)
            me_annots, adata_mz = zmm.match_adata_mz_and_annots(me_annots, adata_mz)

            print("pulling quality of m/z features from annotations")
            me_annots.pull_mz_quality_df()  ## NEW
            if not me_annots.mz_quality[me_annots.offsample_column].isnull().all():
                fdr_cutoff_features_to_compute_spot_qual = 0.1
                if SAMPLE in ["LP-brain-mouse-5"]:
                    fdr_cutoff_features_to_compute_spot_qual = 0.2
                x =  fdr_cutoff_features_to_compute_spot_qual
                adata_mz = zmm.compute_spot_quality(adata_mz, me_annots,
                                                   allowed_in_out_ratio, fdr_cutoff=x)
                adata_mz = zmm.mark_weak_signal_spots(adata_mz, cutoff=100)

                try:
                    adata_mz = custom_out_tissue(
                        adata_mz, mydata_dir, subfolder, SAMPLE)
                except:
                    pass

                plot_out_path = os.path.join(plot_folder, f'{SAMPLE}-spotsBefore.png')

                plot_spots(adata_mz, plot_out_path, title="Before spots selection")
                print(f"plot of initial spots saved to: {plot_out_path}")

                if select_spots_yes_no.lower() in ["yes", "y"]:
                    adata_mz = adata_mz[adata_mz.obs["spot_quality"] == "good", :]
                    plot_out_path = os.path.join(plot_folder,
                                                 f'{SAMPLE}-spotsAfter.png')
                    plot_spots(adata_mz, plot_out_path, title="After spots selection")
                else:
                    print("Spot quality added to object, but original spots preserved")

                adata_mz.uns['history']['spot_number_retained'] = adata_mz.obs.shape[0]
                adata_mz.write_h5ad(os.path.join(mydata_dir, subfolder, f'{SAMPLE}.h5ad'))

            else:
                print(
                    "Skipped!, the out-of-tissue probabilities are not available")
        else:
            print(
                "already computed quality of spots for adata_mz, nothing to do")