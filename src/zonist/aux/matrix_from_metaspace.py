# johaGL
# Get the aligned & recalibrated matrix of all pixels spectra from Metaspace
# (preferred if the given raw imzML+ibd are unaligned & not-recalibrated)
import os
from typing import List

from tqdm import tqdm

import numpy as np
import pandas as pd
import anndata as ad

from metaspace_converter import metaspace_to_anndata
from scipy.sparse import csr_matrix



def build_coord_names(df):
    p1 = "x" + df['ion_image_pixel_x'].astype(str).to_numpy()
    p2 = "_y" + df['ion_image_pixel_y'].astype(str).to_numpy()
    p3 = p1 + p2
    return list(p3)


def adatas_to_matrices_list(metaspace_registry, api_versions, api_fdr):
    mx_dfs = []
    all_mz_with_multiplets = list()
    for key in api_versions.keys():
        adata = metaspace_to_anndata(
            dataset_id=metaspace_registry,
            fdr=api_fdr,
            database=api_versions[key]  # ("HMDB", "v4")
        )
        Xdf = pd.DataFrame(adata.X)
        ok = np.around(adata.var["mz"].astype(float).to_numpy(), 9)
        Xdf.columns = ok
        coordinates_names_list = build_coord_names(adata.obs)
        Xdf.index = coordinates_names_list
        mx_dfs.append(Xdf)
        all_mz_with_multiplets += list(ok)
        del adata

    return set(all_mz_with_multiplets), mx_dfs


def homogenize_rows_order_mx_dfs(mx_dfs):
    result = [mx_dfs[0]]  # initialize with the first df
    ref_rownames = mx_dfs[0].index.tolist()
    for i in range(1, len(mx_dfs)):
        tmp = mx_dfs[i].T  # transpose to select faster by cols
        tmp = tmp[ref_rownames]
        result.append(tmp.T)  # transpose to yield original shape

    return result


def yield_matrix_dataframe_united(ordered_unique_mzs, matrices_dfs_list):
    """yields matrix-dataframe with shape: coordinates x features"""
    united_mx_df = {}
    for mz in tqdm(ordered_unique_mzs):
        if not mz in united_mx_df.keys():
            found = False
            n = 0
            while (n < len(matrices_dfs_list)) and (found != True):
                try:
                    found_mz = matrices_dfs_list[n][mz].to_numpy()
                    found = True
                except KeyError:
                    pass
                n += 1
            # end while
            united_mx_df[mz] = found_mz
        # end if
    # end for
    out_df = pd.DataFrame(united_mx_df)
    out_df.index = matrices_dfs_list[0].index.tolist()
    out_df = out_df.T
    out_df.index.name = "mz"

    return out_df.T


def steps_get_united_df(metaspace_registry: str, api_versions: dict,
         api_fdr: float) -> pd.DataFrame:
    """returns a df with shape: coordinates (rows) x features (columns)"""
    print("fetching data with metaspace-converter (T Alexandrov lab package)")
    unique_mzs, matrices_dfs_list = adatas_to_matrices_list(
        metaspace_registry, api_versions, api_fdr)

    matrices_dfs_list = homogenize_rows_order_mx_dfs(matrices_dfs_list)
    ordered_unique_mzs = sorted(list(unique_mzs))

    outdf = yield_matrix_dataframe_united(ordered_unique_mzs,
                                          matrices_dfs_list)
    return outdf


def save_quantifications_anndata(outdf, outfile) -> None:
    print(f"saving to anndata file, with .X sparse matrix of shape: "
          f"({outdf.shape[0]} rows x {outdf.shape[1]} columns)")
    out_obj = ad.AnnData(X=csr_matrix(outdf),
                         obs=pd.DataFrame({"obs":outdf.index.tolist()},
                                          index=outdf.index.tolist()),
                         var=pd.DataFrame({"var": outdf.columns.tolist()},
                                          index=outdf.columns.tolist()))
    outfile = outfile + ".h5ad"
    out_obj.write(outfile, compression=None, compression_opts=None)
    print(f"End: file saved to '{outfile}'")


def save_quantifications_to_tsv(outdf, outfile) -> None:
    """slow, use instead save anndata above"""
    print(f"saving to tab delimited file of shape: ({outdf.shape[0]} rows x"
          f" {outdf.shape[1]} columns)")
    outfile = outfile + ".tsv"
    outdf.to_csv(outfile, sep='\t', index=True, header=True )
    print(f"End: file saved to '{outfile}'")


def write_united_matrix_to_file(ordered_unique_mzs,
                                matrices_dfs_list, out_file) -> None:
    """
    DEPRECATED
       Deprecated because opening a file with 170K columns is not efficient

    instead saving matrix-df pandas (time consuming), saves each row to file:
    mx  x0_y0  x0_y1  x ....
    33.333  0  0  0 ...
    Output: None
        but writes matrix of shape: features (rows) x coords (cols)
    """
    taken = []
    f = open(out_file, "a")
    firstline = ["mz"] + matrices_dfs_list[0].index.tolist()
    f.write('\t'.join(firstline) + "\n")
    for mz in tqdm(ordered_unique_mzs):
        if not mz in taken: # united_mx_df.keys():
            found = False
            n = 0
            while (n < len(matrices_dfs_list)) and (found != True):
                try:
                    found_mz = matrices_dfs_list[n][mz].to_numpy()
                    found = True
                except KeyError:
                    pass
                n += 1
            # end while
            taken.append(mz)
            iline = [str(mz)] + [str(x) for x in found_mz]
            f.write('\t'.join(iline) + "\n")
        # end if
    # end for
    f.close()


if __name__ == "__main__":
    #  my vars
    metaspace_registry = "2022-08-16_10h14m23s"
    api_versions = {'HMDB': ("HMDB", "v4"),
                   'LipidMaps': ("LipidMaps", "2017-12-12")}
    api_fdr = 0.5
    out_dir = "~/spatial_thesis/spatial_metabolomics/data/whole-ellis/"
    filename = "whole-ellis_mx" # no extension

    out_file = os.path.join(os.path.expanduser(out_dir),
                            filename)
    # start
    out_df = steps_get_united_df(metaspace_registry, api_versions, api_fdr)
    #save_quantifications_to_tsv(out_df, out_file)  # is slow
    save_quantifications_anndata(out_df, out_file)