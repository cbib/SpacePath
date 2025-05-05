import os
import pandas as pd
import anndata as ad

import zonist.data as zd

from local_helpers.local_helpers import (provide_annotation_dict,
                                   part_a_filters)
from local_helpers.visualization import plot_mz_given_quality_metrics
#from bank4zonist.visualization import *
#from classify_normalize_features import part_a_filters

"""
bad or good quality of m/z features
"""
# define here my files and dirs
mydata_dir = "../../data"
my_csv_table_path = "../small-test-parameters.tsv"
plot_parent_folder = "../../figures/"

# start
datasets_parameters_df = pd.read_csv(my_csv_table_path, header=0,
                  index_col=None, sep='\t')

print("\nCheck quality m/z across specified datasets (no-modify object)")

for i, row in datasets_parameters_df.iterrows():
    SAMPLE = row["dataset_name"]
    normalize_yes_no = row["TIC_normalization"]
    offSample_cutoff = row["offSample_cutoff"]
    msm_cutoff = row["msm_cutoff"]
    fdr_cutoff = row["fdr_cutoff"]
    moran_I_cutoff = row["Moran_I_cutoff"]
    moran_padj_cutoff = row["Moran_padj_cutoff"]
    print("\nSample: ", SAMPLE)

    subfolder = SAMPLE
    plot_folder = os.path.join(plot_parent_folder, SAMPLE)
    if not os.path.exists(plot_folder):
        os.mkdir(plot_folder)

    annotations_to_use = row["annotations_to_use"].split(" ")
    anno_dict = provide_annotation_dict(SAMPLE, annotations_to_use)
    me_annots = zd.MetabolomeAnnotations(
        name=SAMPLE, data_dir=mydata_dir,
        subfolder=subfolder, annotation_files=anno_dict)

    me_annots.load()


    adata_mz = ad.read_h5ad(filename=os.path.join(
        mydata_dir, subfolder, f'{SAMPLE}.h5ad'))

    try:
        columns_to_reset = set(adata_mz.var.columns.tolist()) - set(["mz"])
        adata_mz.var = adata_mz.var.drop(columns=list(columns_to_reset))
    except KeyError:
        pass
    print("Loaded object", adata_mz.uns["history"])

    if 'features_n_after_mz_clean_2' in list(adata_mz.uns['history'].keys()):
        print("Nothing to do for this sample.")
    elif not "is_X_TIC_normalized" in adata_mz.uns.keys():
        print("no normalization yet applied")
    me_annots, ok_mz_keep_part1 = part_a_filters(
            me_annots, msm_cutoff, offSample_cutoff, fdr_cutoff
            )
    text_notice = f"valid $m/z$\n n = {len(ok_mz_keep_part1)}"
    plot_mz_given_quality_metrics(me_annots.mz_quality, SAMPLE,
                                 plot_folder,
                                msm_cutoff, offSample_cutoff, text_notice)

    print("only figures produced, adata object intact")
