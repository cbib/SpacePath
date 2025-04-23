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
import zonist.data as zd

from bank4zonist.quantifs_methods import *
from bank4zonist.constants import provide_annotation_dict


def compute_spot_quality(adata_mz, me_annots,
                         allowed_in_out_ratio, fdr_cutoff):
    print("marking spot quality with params: ", allowed_in_out_ratio)
    adata_mz = adata_mz_mark_spot_quality(
        adata_mz, me_annots.mz_quality, fdr_cutoff=fdr_cutoff
    )
    adata_mz.obs['spot_quality'] = "good"
    adata_mz.obs.loc[adata_mz.obs['ratio_in_out'] < allowed_in_out_ratio,
    'spot_quality'] = "bad"
    try:
        if len(adata_mz.uns['spot_quality_colors']
               ) != len(adata_mz.obs['spot_quality'].unique()):
            del adata_mz.uns['spot_quality_colors']
    except KeyError:
        pass
    return adata_mz


def mark_weak_spots(adata_mz, cutoff: float = 1000):  # todo: deprecated
    adata_mz.obs['totalIntensitySpot'] =  adata_mz.X.sum(axis=1)
    adata_mz.obs.loc[adata_mz.obs['totalIntensitySpot'] <= cutoff, 'spot_quality'] = 'bad'
    return adata_mz


def plot_spots(adata_mz, plot_out_path, title='my title'):
    print("Saving plot of spots")
    sq.pl.spatial_scatter(
        adata_mz, color=["spot_quality"], img=False,
        title=title, save=os.path.abspath(plot_out_path))
    plt.close()

    return adata_mz


# define here my input
mydata_dir = "../../data"
my_csv_table_path = "../private-datasets-parameters.tsv"
plot_parent_folder = "../../figures/"

# start
datasets_parameters_df = pd.read_csv(my_csv_table_path, header=0,
                  index_col=None, sep='\t')

force_for_these = [] # leave empty if no samples to force here   # e.g. ["TD-hemi-635-1-A", "TD-hemi-635-1-B"]    "TD-635-1-B"


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
        me_annots, adata_mz = match_adata_mz_and_annots(me_annots, adata_mz)

        print("pulling quality of m/z features from annotations")
        me_annots.pull_mz_quality_df()  ## NEW
        if not me_annots.mz_quality[me_annots.offsample_column].isnull().all():
            #if "spot_quality" not in adata_mz.obs.columns.tolist():
            fdr_cutoff_features_to_compute_spot_qual = 0.1
            if SAMPLE in ["LP-brain-mouse-5"]:
                fdr_cutoff_features_to_compute_spot_qual = 0.2
            x =  fdr_cutoff_features_to_compute_spot_qual
            adata_mz = compute_spot_quality(adata_mz, me_annots,
                                 allowed_in_out_ratio, fdr_cutoff=x)
            # adata_mz = mark_weak_spots(adata_mz, cutoff=1e4)  #  useless, no effect
            plot_out_path = os.path.join(plot_folder, f'{SAMPLE}-spotsBefore.png')

            plot_spots(adata_mz, plot_out_path, title="Before spots selection")
            print(f"plot of initial spots saved to: {plot_out_path}")
            if select_spots_yes_no.lower() in ["yes", "y"]:
                adata_mz = adata_mz[adata_mz.obs["spot_quality"] == "good", :]
                plot_out_path = os.path.join(plot_folder,
                                             f'{SAMPLE}-spotsAfter.png')
                plot_spots(adata_mz, plot_out_path, title="After spots selection")
            else:
                print("Spot quality added to object, but original spots kept")

            adata_mz.uns['history']['spot_number_retained'] = adata_mz.obs.shape[0]
            adata_mz.write_h5ad(os.path.join(mydata_dir, subfolder, f'{SAMPLE}.h5ad'))

        else:
            print(
                "Skipped!, the out-of-tissue probabilities are not available")
    else:
        print(
            "already computed quality of spots for adata_mz, nothing to do")