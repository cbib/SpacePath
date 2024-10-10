import os
import pickle
from typing import List
import pandas as pd
import numpy as np
import anndata as ad
import pyimzml.ImzMLWriter as imzmlw
import gc
import scanpy as sc
import harmonypy as hm
import squidpy as sq
from tqdm import tqdm
from helpers import (select_n_most_variable_features, get_rid_of_low_variance_mz,
                     find_common_mz__two_adata,
                     run_approximate_match_two_adata, rounded_mz_as_variable_index_adata,
                     plot_gross_initial_clusters, plot_post_filtered_clusters,
                     harmonize_mz_identifiers_data2, concatenate_finally,
                     save_visualizations_after_concat, run_harmony)
# scanpy uses harmonypy (pip install)
# #https://www.sc-best-practices.org/surface_protein/batch_correction.html



def only_use_few_peaks_for_initial_tests(adata1, adata2, df):
    h1 = df.loc[df['mzs1'].isna(), :].sample(n=20)  # iloc[0:20, :]
    h2 = df.loc[df['mzs2'].isna(), :].sample(n=20)  # iloc[0:20, :]
    h3 = df.loc[~df['mzs1'].isna(), :].loc[
         ~df['mzs2'].isna(), :].sample(n=800)
    uf = pd.concat((h1, h2, h3), axis=0)

    def ioio(adata, list_of_mz):
        dfvar = adata.var.copy().T
        U = pd.DataFrame(adata.X)
        U.columns = dfvar.columns.tolist()
        dfvar = dfvar[list_of_mz]
        U = U[list_of_mz]

        out = ad.AnnData(np.array(U),
                   obsm={"spatial": np.array(adata.obsm['spatial'])},  # coords
                   varm={"transcripts" : np.array(list_of_mz)})
        out.obs = adata.obs.copy()
        out.var = dfvar.T
        return out

    tmpl = uf.loc[~uf['mzs1'].isna(), "mzs1"].tolist()
    adata1 = ioio(adata1, tmpl)
    tmpl = uf.loc[~uf['mzs2'].isna(), "mzs2"].tolist()
    adata2 = ioio(adata2, tmpl)
    gc.collect()

    return adata1, adata2, uf


def save_adata_to_imzml(adata, out_file_name):
    out_file_name_full = os.path.join("output/objects", out_file_name)
    coordinates_arr = adata.obsm['spatial']
    with imzmlw.ImzMLWriter(out_file_name_full,
                            mode="processed") as imzml:
        for n in tqdm(range(0, len(coordinates_arr))):
            imzml.addSpectrum(mzs=adata.var['rounded_mz'], #  adata.var.index.astype(float),
                              intensities=adata.X[n],
                              coords=coordinates_arr[n])
    print("saved imzML!\n")


def open_adata__filter__save_imzml(filename, n_top,
                                   leiden_resolution,
                                   bad_leiden_clusters: List[int],
                                   out_substring0, out_substring1, out_file
                                   ) -> None:
    print("opening...")
    with open(filename, 'rb') as fi:
        adata = pickle.load(fi)

    plot_gross_initial_clusters(adata, out_substring0, leiden_resolution)
    #h0 = adata.X[1000, :].sum()
    gc.collect()
    print("pca... neigh... leiden")
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    #h1 = adata.X[1000, :].sum()
    sc.tl.leiden(adata, resolution=leiden_resolution)

    # adata.obs['pixel_name'] = adata.obs_names.tolist()     #  ** too heavy **

    # good_pixels = adata.obs.loc[
    #     ~adata.obs['leiden'].astype(int).isin(bad_leiden_clusters),
    #     "pixel_name"]

    good_pixels = adata.obs.loc[
             ~adata.obs['leiden'].astype(int).isin(bad_leiden_clusters),:].index.tolist()
    adata = adata[good_pixels, :]

    # feat extraction (simple, with variance)
    print("add stuff")
    adata.var['variance'] = np.var(adata.X, axis=0)
    adata.var['mz'] = adata.var.index.astype(float).to_numpy()
    adata.var['rounded_mz'] = np.around(adata.var['mz'].to_numpy(),
                                        decimals=9)
    if len(adata.var['rounded_mz'].tolist()) > len(
            adata.var['rounded_mz'].unique()):
        print("oh no ! ")
        print(len(adata.var['rounded_mz'].unique(),
                  "--> unique rounded_mz has less elements that mz list"))
        raise ValueError
    print("selecting n most var feat")
    adata = select_n_most_variable_features(adata, n_top=n_top)
    adata.var = adata.var.drop(columns=["mz", "variance"])  # to avoid overloading memory

    plot_post_filtered_clusters(adata, out_substring1)
    print("plot ok ! ")
    save_adata_to_imzml(adata, out_file)



def open_adata__filterB__save_imzml(filename,
                                   leiden_resolution,
                                   bad_leiden_clusters: List[int],
                                   out_substring0, out_substring1, out_file_name
                                   ) -> None:
    """retain all features, but exclude pixels of the chemical matrix"""
    with open(filename, 'rb') as fi:
        adata = pickle.load(fi)

    adata.obs['pixel_name'] = adata.obs_names.tolist()
    plot_gross_initial_clusters(adata, out_substring0, leiden_resolution)

    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=leiden_resolution)

    good_pixels = adata.obs.loc[
        ~adata.obs['leiden'].astype(int).isin(bad_leiden_clusters),
        "pixel_name"]
    adata = adata[good_pixels, :]
    plot_post_filtered_clusters(adata, out_substring1)
    # save_adata_to_imzml(adata, out_file)   # bad as expects 'rounded_mz' in .var
    out_file_name_full = os.path.join("output/objects", out_file_name)
    coordinates_arr = adata.obsm['spatial']
    with imzmlw.ImzMLWriter(out_file_name_full,
                            mode="processed") as imzml:
        for n in tqdm(range(0, len(coordinates_arr))):
            imzml.addSpectrum(mzs=adata.var.index.astype(float).to_numpy(),
                              intensities=adata.X[n],
                              coords=coordinates_arr[n])
    print("saved imzML !")


if __name__ == "__main__":
    wdir = '/home/johanna/spatial_thesis/spatial_metabolomics/'
    ctrl_f = "output/objects/Ctrl--anndata.pk"
    ko_f = "output/objects/KO--anndata.pk"
    # cutoff_percentile_variance_tolerated = 0.7
    os.chdir(wdir)

    # with open(ko_f, 'rb') as fi:
    #     ko_ad = pickle.load(fi)
    #
    # with open(ctrl_f, 'rb') as f:
    #     ctrl_ad = pickle.load(f)

    # j = find_common_mz__two_adata(ctrl_ad, ko_ad)

    # ctrl_ad = compute_metrics_for_variables(ctrl_ad)
    # ko_ad = compute_metrics_for_variables(ko_ad)

    # ctrl_ad = get_rid_of_low_variance_mz(
    #     ctrl_ad, cutoff_percentile_variance_tolerated)
    # ko_ad = get_rid_of_low_variance_mz(
    #     ko_ad, cutoff_percentile_variance_tolerated)

    # ctrl_ad = select_n_most_variable_features(ctrl_ad, n_top=12000)
    # ko_ad = select_n_most_variable_features(ko_ad, n_top=12000)

    # run_blah(ctrl_ad, "control_ft1")
    # run_blah(ko_ad, "ko_ft1")

    # # save '.8' means 80% features (the most variable) were retained, bad pixels excluded too.

    open_adata__filter__save_imzml(ko_f, n_top=30640, # total is 38300
                                   leiden_resolution=0.2,
                                   bad_leiden_clusters=[2,4], #  BUT  Only 0 if ctrl
                                   out_substring0="ko.8",
                                           out_substring1="ko.f.8",
                                           out_file="ko_ft")


    open_adata__filter__save_imzml(ctrl_f, n_top=36560,  # total is 45700
                                   leiden_resolution=0.25,
                                   bad_leiden_clusters=[0],
                                       out_substring0="ctrl.8",
                                       out_substring1="ctrl.f.8",
                                       out_file="ctrl_ft")



    # save all features but exclude bad pixels
    # open_adata__filterB__save_imzml(ko_f,
    #                                leiden_resolution=0.2,
    #                                bad_leiden_clusters=[2, 4],
    #                                # 0 6 .. for kko.  Only 0 if ctrl
    #                                out_substring0="ko",
    #                                out_substring1="ko.f",
    #                                out_file_name="ko_ft2")
    #
    # open_adata__filterB__save_imzml(ctrl_f,
    #                                leiden_resolution=0.25,
    #                                bad_leiden_clusters=[0],
    #                                out_substring0="ctrl",
    #                                out_substring1="ctrl.f",
    #                                out_file_name="ctrl_ft2")
    import sys
    #sys.exit()  # TODO: organize:NOT appropriate for src/. Moreover, Combination has not yet been done


    # # TODO: code from here must migrate somewhere else:
    # # ctrl_ad.obs['pixel_name'] = ctrl_ad.obs_names.tolist()
    # # ko_ad.obs['pixel_name'] = ko_ad.obs_names.tolist()
    # df = run_approximate_match_two_adata(ctrl_ad, ko_ad, 0.01)
    # df.to_csv("fantastic.csv")
    # ctrl_ad = rounded_mz_as_variable_index_adata(ctrl_ad)
    #
    # ko_ad = rounded_mz_as_variable_index_adata(ko_ad)
    #
    # ctrl_ad, ko_ad, df = only_use_few_peaks_for_initial_tests(ctrl_ad, ko_ad, df)  # TODO: delete this line when real outputs desired
    #
    # #plot_gross_initial_clusters(ctrl_ad, out_substring="ctrl")
    # #plot_gross_initial_clusters(ko_ad, out_substring="ko")
    #
    # ctrl_ad.obs['pixel_name'] = ctrl_ad.obs_names.tolist()
    # ko_ad.obs['pixel_name'] = ko_ad.obs_names.tolist()
    #
    # ko_ad, df2 = harmonize_mz_identifiers_data2(ko_ad, df)
    #
    # # coords of each dataset for later use
    # coords_ctrl = ctrl_ad.obsm['spatial']
    # coords_ko = ko_ad.obsm['spatial']
    #
    # bigo = concatenate_finally(ctrl_ad, ko_ad)
    # save_visualizations_after_concat(bigo)  # this shows that there is a cluster of the pixels that represent the chemical support.
    #
    # del ctrl_ad
    # del ko_ad
    # gc.collect()
    # p = 33
    # bigo = run_harmony(bigo)
    #
    # print('nerviosssss... harmony  to do ! ')
    # print(bigo)
    # # harmony uses X_pca slot
    #
    # print("ended")
    # r = 55




