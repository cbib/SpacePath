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
                         allowed_in_out_ratio, plot_folder):
    print("marking spot quality with params: ", allowed_in_out_ratio)
    adata_mz = adata_mz_mark_spot_quality(
        adata_mz, me_annots.mz_quality
    )
    adata_mz.obs['spot_quality'] = "good"
    adata_mz.obs.loc[adata_mz.obs['ratio_in_out'] < allowed_in_out_ratio,
    'spot_quality'] = "bad"

    tmpstr = os.path.join(plot_folder, f'{SAMPLE}-spotsBefore.png')
    sq.pl.spatial_scatter(
        adata_mz, color=["spot_quality"], img=False,
        title="Before spots selection", save=os.path.abspath(tmpstr))
    plt.close()
    print(f"plot of initial spots in: {tmpstr}")
    return adata_mz


def keep_only_good_spots(adata_mz, plot_folder):
    adata_mz = adata_mz[adata_mz.obs["spot_quality"] == "good", :]
    print("Will save plot and over-write adata object")
    # sq.pl.spatial_scatter(
    #     adata_mz, color=["spot_quality"], img=False,
    #     title="After spots selection")
    # plt.show()

    tmpstr = os.path.join(plot_folder,
                          f'{SAMPLE}-spotsAfter.png')
    sq.pl.spatial_scatter(
        adata_mz, color=["spot_quality"], img=False,
        title="After spots selection", save=os.path.abspath(tmpstr))
    plt.close()

    return adata_mz


# define here my input
mydata_dir = "../../data"
my_csv_table_path = "../datasets-parameters.tsv"
plot_parent_folder = "../../figures/"

# start
datasets_parameters_df = pd.read_csv(my_csv_table_path, header=0,
                  index_col=None, sep='\t')

force_for_these = [] # todo: defin if needed #
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

        me_annots, adata_mz = match_adata_mz_and_annots(me_annots, adata_mz)

        print("pulling quality of m/z features from annotations")
        me_annots.pull_mz_quality_df()  ## NEW
        if not me_annots.mz_quality[me_annots.offsample_column].isnull().all():
            if "spot_quality" not in adata_mz.obs.columns.tolist():
                compute_spot_quality(adata_mz, me_annots,
                                     allowed_in_out_ratio, plot_folder)
            if select_spots_yes_no.lower() in ["yes", "y"]:
                adata_mz = keep_only_good_spots(adata_mz, plot_folder)
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



# any_anno_dict = {'HMDB': f"{SAMPLE}-HMDB.csv",
#              "KEGG": f"{SAMPLE}-KEGG.csv",
#              "LipidMaps": f"{SAMPLE}-LipidMaps.csv",
#              "CoreMetabolome": f"{SAMPLE}_CoreMetabolome.csv",
#              }