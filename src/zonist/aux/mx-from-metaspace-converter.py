import os
from typing import List

from tqdm import tqdm

import numpy as np
import pandas as pd

from metaspace_converter import metaspace_to_anndata


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
    """yields matrix-dataframe with shape features x coordinates"""
    foo = {}
    for mz in tqdm(ordered_unique_mzs):
        if not mz in foo.keys():
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
            foo[mz] = found_mz
        # end if
    # end for
    # yield features (rows) x coords (cols) as Metaspace's Intensities mx
    out_df = pd.DataFrame(foo)
    out_df.index = matrices_dfs_list[0].index.tolist()
    out_df = out_df.T
    out_df.index.name = "mz"

    return out_df


def write_united_matrix_to_file(ordered_unique_mzs,
                                matrices_dfs_list, out_file) -> None:
    """
    instead loading a matrix-df in memory, saves each row to file:
    mx  x0_y0  x0_y1  x ....
    33.333  0  0  0 ...
    """
    taken = []
    f = open(out_file, "a")
    firstline = ["mz"] + matrices_dfs_list[0].index.tolist()
    f.write("\t".join(firstline) + "\n")
    for mz in tqdm(ordered_unique_mzs):
        if not mz in taken:#foo.keys():
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
            f.write("\t".join(iline) + "\n")
        # end if
    # end for
    # yield features (rows) x coords (cols) as Metaspace's Intensities mx
    f.close()

def run_mx_from_metaspace_converter(metaspace_registry, api_versions, api_fdr,
                                    out_file) -> None:
    print("fetching data with metaspace-converter (T Alexandrov lab package)")
    unique_mzs, matrices_dfs_list = adatas_to_matrices_list(
        metaspace_registry, api_versions, api_fdr)

    matrices_dfs_list = homogenize_rows_order_mx_dfs(matrices_dfs_list)
    ordered_unique_mzs = sorted(list(unique_mzs))

    print(f"writing quantifications (==intensities) matrix to tsv file "
          f"({len(ordered_unique_mzs)} rows x"
          f" {matrices_dfs_list[0].shape[0]} columns)")
    write_united_matrix_to_file(ordered_unique_mzs,
                                matrices_dfs_list, out_file)
    # # NOTE: saving with pandas takes a long long long time.
    # # so keep the following lines commented (not active) :
    # cool = yield_matrix_dataframe_united(ordered_unique_mzs, matrices_dfs_list)
    # print(f"saving to tsv file ({cool.shape[0]} rows x {cool.shape[1]} columns)")
    # cool.to_csv(out_file, sep='\t', index=True, header=True )
    print("end")


if __name__ == "__main__":
    #  my vars
    metaspace_registry = "2022-08-16_10h14m23s"
    api_versions = {'HMDB': ("HMDB", "v4"),
                   'LipidMaps': ("LipidMaps", "2017-12-12")}
    api_fdr = 0.5
    out_dir = "~/spatial_thesis/spatial_metabolomics/data/whole-ellis/"
    filename = "whole-ellis_mx.tsv"

    # start !
    out_file = os.path.join(os.path.expanduser(out_dir),
                            filename)
    run_mx_from_metaspace_converter(metaspace_registry, api_versions, api_fdr,
                                    out_file)








