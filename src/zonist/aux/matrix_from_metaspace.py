# johaGL
# Get the aligned & recalibrated annotated data from Metaspace Server directly
#    which is preferred if the given raw imzML is unaligned & not-recalibrated
import os
from typing import List, Tuple

from tqdm import tqdm

import numpy as np
import pandas as pd
import anndata as ad

from metaspace_converter import metaspace_to_anndata
from scipy.sparse import csr_matrix
import zonist.helpers as helpers
from zonist.constants import mz_decimal_places
import gc
import logging

logger = logging.getLogger(__name__)
logger = helpers.setup_logger(logger)


def build_coord_names(df) -> List[str]:
    p1 = "x" + df['ion_image_pixel_x'].astype(str).to_numpy()
    p2 = "_y" + df['ion_image_pixel_y'].astype(str).to_numpy()
    pixels_names_l = p1 + p2
    return list(pixels_names_l)


def adatas_to_matrices_list(metaspace_registry, api_versions, api_fdr
)->(set, List[pd.DataFrame]):
    matrices_dfs_list = []
    annotations_dfs_dict = {}
    all_mz_with_multiplets = list()
    for key in api_versions.keys():
        adata = metaspace_to_anndata(
            dataset_id=metaspace_registry,
            fdr=api_fdr,
            database=api_versions[key]  # ("HMDB", "v4")
        )
        if "formula" not in adata.var.columns:
            adata.var["formula"] = adata.var.index.tolist()
        adata.var = adata.var.drop(columns=["isotopeImages"])
        annotations_dfs_dict[key] = adata.var.copy()
        Xdf = pd.DataFrame(adata.X)
        mz_arr = np.around(adata.var["mz"].astype(float).to_numpy(),
                           decimals=mz_decimal_places)
        Xdf.columns = mz_arr
        coordinates_names_list = build_coord_names(adata.obs)
        Xdf.index = coordinates_names_list
        matrices_dfs_list.append(Xdf)
        all_mz_with_multiplets += list(mz_arr)
        del adata

    return set(all_mz_with_multiplets), matrices_dfs_list, annotations_dfs_dict


def homogenize_rows_order_mx_dfs(mx_dfs) -> List[pd.DataFrame]:
    matrices_dfs_list = [mx_dfs[0]]  # initialize with the first df
    ref_rownames = mx_dfs[0].index.tolist()
    for i in range(1, len(mx_dfs)):
        tmp = mx_dfs[i].T  # transpose to select faster by cols
        tmp = tmp[ref_rownames]
        matrices_dfs_list.append(tmp.T)  # transpose to yield original shape

    return matrices_dfs_list


def yield_matrix_from_read_adatas(
        ordered_unique_mzs, matrices_dfs_list) -> pd.DataFrame:
    """yields matrix-dataframe with shape: coordinates x features"""
    united_mx_dict = {}
    spots_list = matrices_dfs_list[0].index.tolist()
    for mz in tqdm(ordered_unique_mzs):
        if not mz in united_mx_dict.keys():
            found_mz = None
            df_index = 0
            while (df_index < len(matrices_dfs_list)) and (found_mz is None):
                try:
                    found_mz = matrices_dfs_list[df_index][mz].to_numpy()
                    try:   # handle repeated rows if any (ok, checked)
                        if found_mz.shape[1] >= 2:
                            found_mz = found_mz[:,0]
                    except IndexError:
                        pass
                except KeyError:
                    pass
                df_index += 1
            # end while
            united_mx_dict[mz] = found_mz
        # end if
    # end for
    out_df = pd.DataFrame(united_mx_dict)
    out_df.index = spots_list
    out_df = out_df.T
    out_df.index.name = "mz"

    return out_df.T


def steps_get_united_df(metaspace_registry: str, api_versions: dict,
         api_fdr: float) -> Tuple[pd.DataFrame, dict]:
    """returns a df with shape: coordinates (rows) x features (columns)"""
    logger.info("fetching data with metaspace-converter (T Alexandrov lab package)")

    unique_mzs, matrices_dfs_list, annotations_dfs_dict = adatas_to_matrices_list(
        metaspace_registry, api_versions, api_fdr)
    # note: annotations from Metaspace python API: no offSampleProb floats (instead only offSample bool)
    matrices_dfs_list = homogenize_rows_order_mx_dfs(matrices_dfs_list)
    ordered_unique_mzs = sorted(list(unique_mzs))

    outdf = yield_matrix_from_read_adatas(ordered_unique_mzs,
                                          matrices_dfs_list)
    return outdf, annotations_dfs_dict


def save_quantifications_anndata(outdf, outfile) -> None:

    coordinates_arr = helpers.generate_coordinates_arr(list(outdf.index))
    out_obj = ad.AnnData(X=csr_matrix(outdf),
                         obs=pd.DataFrame({"obs":outdf.index.tolist()},
                                          index=outdf.index.tolist()),
                         var=pd.DataFrame(
                             {"mz": np.array(
                                 [float(x) for x in outdf.columns.tolist()])},
                             index=outdf.columns.tolist()
                         ),
                         obsm={"spatial": coordinates_arr})
    out_obj.uns['spatial'] = {
        'MY_SLICE7_NAME': {'scalefactors':
                               {'fiducial_diameter_fullres': 100,
                                'spot_diameter_fullres': 1,
                                'tissue_hires_scalef': 0.5,
                                'tissue_lowres_scalef': 0.5}}}
    del outdf
    gc.collect()
    outfile = outfile + ".h5ad"
    out_obj.write(outfile, compression=None, compression_opts=None)
    logger.info("done")


if __name__ == "__main__":
    metaspace_registry = "2023-10-24_17h13m08s"
    api_versions = {'HMDB': ("HMDB", "v4"),
                   'LipidMaps': ("LipidMaps", "2017-12-12")}
    api_fdr = 0.5
    out_dir = "../../../data"
    filename = "basictest" # no extension

    out_file = os.path.join(os.path.expanduser(out_dir),
                            filename)
    out_df, annotations_dfs_dict = steps_get_united_df(
        metaspace_registry, api_versions, api_fdr)
