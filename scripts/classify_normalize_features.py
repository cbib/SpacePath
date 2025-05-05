import sys
import os
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
import seaborn as sns
import squidpy as sq
import zonist.data as zd

import zonist.processing.measures_methods as zmm

from local_helpers.local_helpers import (provide_annotation_dict,
  part_a_filters,
    part_b_filters_with_autocorr)

from local_helpers.visualization import plot_mz_given_quality_metrics, plot__moran_info

#from local_helpers.visualization import *



if __name__ == "__main__":
    """
    bad or good quality of m/z features
    """
    # define here my files and dirs
    mydata_dir = "../../data"
    my_csv_table_path = "../small-test-parameters.tsv"   # "../public-datasets-parameters.tsv"
    plot_parent_folder = "../../figures/"

    force_for_these = ["LP-example0"]    # keep empty if no sample to force   # e.g.  , "TD-635-1-A", "TD-635-1-B"   "TD-hemi-635-1-A", "TD-hemi-635-1-B"

    # start
    datasets_parameters_df = pd.read_csv(my_csv_table_path, header=0,
                      index_col=None, sep='\t')

    print("\nMarks bad quality m/z across specified datasets (.var property)")

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

        adata_mz = ad.read_h5ad(filename=os.path.join(
            mydata_dir, subfolder, f'{SAMPLE}.h5ad'))

        odfsodf = 99
        # todo qjdfqlslqk
        try:
            adata_mz.X = adata_mz.layers['raw_X']
            del adata_mz.uns['is_X_TIC_normalized']
        except:
            print("will continue")
            pass
        # end of todo
        print("continues here")

        try:
            columns_to_reset = set(adata_mz.var.columns.tolist()) - set(["mz"])
            adata_mz.var = adata_mz.var.drop(columns=list(columns_to_reset))
        except KeyError:
            pass
        print("Loaded object", adata_mz.uns["history"])

        if (SAMPLE not in force_for_these) and ('features_n_after_mz_clean_2' in list(adata_mz.uns['history'].keys())):
            print("Nothing to do for this sample.")
        else:
            me_annots.load()  # load annotations
            me_annots, adata_mz = zmm.match_adata_mz_and_annots(me_annots,
                                                            adata_mz)
            if (normalize_yes_no.lower() in ['y', 'yes']) and (
                    not "is_X_TIC_normalized" in adata_mz.uns.keys()):
                adata_mz.X = csr_matrix(zmm.tic_normalization(
                    adata_mz.X.copy(), adata_mz.obs['obs'].tolist(),
                    adata_mz.var['mz'].tolist()))

            me_annots, ok_mz_keep_part1 = part_a_filters(
                me_annots, msm_cutoff, offSample_cutoff, fdr_cutoff )
            text_notice = f"valid $m/z$\n n = {len(ok_mz_keep_part1)}"
            plot_mz_given_quality_metrics(me_annots.mz_quality, SAMPLE,
                                          plot_folder,
                                          msm_cutoff, offSample_cutoff,
                                          text_notice)

            adata_mz.uns["history"]['features_n_after_mz_clean_1'] = len(ok_mz_keep_part1)

            all_quality_df, copy_mz, keep_part2, moran_mz_df = part_b_filters_with_autocorr(
                adata_mz, ok_mz_keep_part1, me_annots,
                moran_I_cutoff, moran_padj_cutoff)

            text_notice = f"valid $m/z$: n = {len(keep_part2)}"

            plot__moran_info(copy_mz.uns['moranI'].copy(), plot_folder,
                             SAMPLE,
                             text_notice, moran_I_cutoff, moran_padj_cutoff)


            del copy_mz
            adata_mz.uns['moranI'] = moran_mz_df
            adata_mz.uns["history"]['features_n_after_mz_clean_2'] = len(keep_part2)
            all_quality_df['final_quality'] = "bad"
            all_quality_df.loc[(all_quality_df['quality_msm_off_fdr'] == "good") & (
                all_quality_df['quality_autocorr'] == "good"), 'final_quality'
            ] = "good"
            okok = all_quality_df.loc[all_quality_df['final_quality']=="good",:].shape[0]
            if np.all(np.array(adata_mz.var.index.tolist()) == np.array(
                    all_quality_df['mz'].tolist())):
                all_quality_df.index = all_quality_df['mz'].astype(
                    str).tolist()
                all_quality_df['mz'] = all_quality_df['mz'].astype(
                    float).tolist()
                adata_mz.var = all_quality_df
                print(adata_mz)
                print(adata_mz.uns['history'])
            else:
                raise AssertionError("Fatal: could not assembly final object. Aborted.")
            answer = input(
                f"A total of {okok} m/z have good final_quality (pass both a) and b) filters)."
                f"\nDo you agree to continue?")
            if not answer.lower() in ["y", "yes"]:
                print("going to the next, nothing saved")
                continue
            else:
                print("Usable m/z elements (n of features):", okok)
                adata_mz.write_h5ad(filename=os.path.join(
                   mydata_dir, subfolder, f"{SAMPLE}.h5ad"))
                print("saved object")

