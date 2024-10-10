# johaGL 2024
# usage:
#   python3 spectra2anndata_complex.py $DATA_DIR $FILEBASENAME $PREFIX $OUTDIR
# where:
#   $FILEBASENAME is the file name without extension
#       Note: both files .imzML and .ibd must have same FILEBASENAME
#   $OUTDIR is the output folder and must contain 'objects' subfolder
# example:
#  python3 spectra2anndata_complex.py $data_dir \
#     perfbrain1-5_04m_nedc_14x_50um_recal Pub $output_dir

import os
import sys
import numpy as np
import pandas as pd
import anndata as ad
import pyimzml.ImzMLParser as imzmlp
import pickle
from tqdm import tqdm

"""
must be used instead raw_anndata_converter.py only in the case that
mz 'indices' DO NOT match perfectly across all pixels
"""


def rescue_dict(coord_spectra_object):
    """
    from MSI coordinates object (read with pyimzml library),
    constructs a dictionary of dataframes, one by coordinate:
      each dataframe has columns mz and coordinate integer (first coordinate is 0)
    """
    print("Preparing dictionary of separated dataframes (one by coordinate)."
          " Preliminary step")
    dd = dict()
    for i in tqdm(range(len(coord_spectra_object.coordinates))):
        mzs, intensities = coord_spectra_object.getspectrum(i)
        features = np.around(np.array(mzs), 9)
        # decimals 9 as is the case in Metaspace outputs
        curr_df = pd.DataFrame({'mz': features, i: intensities})
        curr_df = curr_df.loc[curr_df[i] >= 50,:]
        dd[i] = curr_df
    return dd


def update_mz_dict(df, mz_dict):
    for u in df['rounded_mz'].unique():
        mzs = df.loc[df['rounded_mz'] == u, 'mz'].tolist()
        if u not in mz_dict.keys():
            mz_dict[u] = mzs
        else:
            for k in mzs:
                [ mz_dict[u].append(k) for k in mzs if k not in mz_dict[u] ]
    return mz_dict


def merge_inexact_matches_dfs(df1, df2, k):
    df2 = df2.drop(columns=["mz"])
    df2 = df2.sort_values(by=k, ascending=False)
    df2 = df2.drop_duplicates(subset='rounded_mz', keep='first')
    df1 = df1.merge(df2, on="rounded_mz",  how="outer")
    return df1


def run_merge_inexact_matches_dfs(df1, df2, k, mz_dict, MYFLOATPLACES):
    df2['rounded_mz'] = np.around(df2['mz'].to_numpy(), MYFLOATPLACES)
    mz_dict = update_mz_dict(df2, mz_dict)
    df3 = merge_inexact_matches_dfs(df1, df2, k)
    if df3.shape[0] > 200:
        df3 = df3.loc[df3.var(axis=1, ddof=0) >= 2, :]
        df3 = df3.loc[df3.sum(axis=1) >= 1000, :]
    return df3, mz_dict


def closest_match_heavy_step(dd, MYFLOATPLACES):
    firstcol = list(dd.keys())[0]
    # set first coordinate (mz and intensities) as dataframe:
    df1 = dd[firstcol].copy(deep=True)
    df1['rounded_mz'] = np.around(df1['mz'].to_numpy(), MYFLOATPLACES)
    mz_dict = dict()
    mz_dict = update_mz_dict(df1, mz_dict)
    df1 = df1.drop(columns=["mz"])
    df_final = pd.DataFrame(
        {'rounded_mz': df1.rounded_mz, firstcol: df1[firstcol]})
    df_final = df_final.sort_values(by=firstcol, ascending=False)
    df_final = df_final.drop_duplicates(subset='rounded_mz', keep='first')
    print("Assembling matrix, as dataframe. Be patient (up to hours)")
    # parse the rest of coordinates
    for i in tqdm(range(1, len(list(dd.keys())))):
        k = list(dd.keys())[i]
        df2 = dd[k]
        df_final, mz_dict = run_merge_inexact_matches_dfs(
            df_final, df2, k, mz_dict, MYFLOATPLACES)

    return df_final, mz_dict


def set_coords_2d_array(co_obj):
    coordinates_list = list()
    for i, item in enumerate(co_obj.coordinates):
        coordinates_list.append((item[0], item[1]))

    return np.array(coordinates_list)


def save__anndata_and_coord_dict(coords_2d_arr, final_df, file_basename,
                                 out_dir, barcodes_prefix):
    new_coord_names = list()
    for tup in coords_2d_arr:
        new_coord_names.append(f'{barcodes_prefix}_{tup[0]}_{tup[1]}')
    print(final_df.iloc[0:5, 0:3])
    final_df = final_df.set_index('rounded_mz')
    final_df = final_df.fillna(0)
    final_df.columns = new_coord_names
    anndatafile = os.path.join(out_dir, f"{file_basename}--anndata.pk")
    adata = ad.AnnData(final_df.T,
                   obsm={"spatial": np.array(coords_2d_arr)},
                   varm={"transcripts" : np.array(final_df.index.tolist())})
    with open(anndatafile, "wb") as f:
        pickle.dump(adata, f)
    print(f"saved anndata object (pickle) to: {anndatafile}")


def load_and_run(data_dir):
    data_dir = os.path.expanduser(sys.argv[1])
    out_dir = os.path.join(os.path.expanduser(sys.argv[4]), "objects")
    file_basename = sys.argv[2]
    MYFLOATPLACES = 3
    imzml_fn = os.path.join(data_dir, f"{file_basename}.imzML")
    ibd_fn = os.path.join(data_dir, f"{file_basename}.ibd")
    co_obj = imzmlp.ImzMLParser(filename=imzml_fn, parse_lib=None,
                                ibd_file=open(ibd_fn, "rb"),
                                include_spectra_metadata=None)
    dd = rescue_dict(co_obj)  # TODO: activate this !!
    final_df, mz_dict = closest_match_heavy_step(dd,
                                                 MYFLOATPLACES)  # TODO: activate this !!
    """
    launched with these options :
    datadir="~/spatial_thesis/spatial_metabolomics/tutorial_publicdata/metaspace_data"
    outputdir="~/spatial_thesis/spatial_metabolomics/local_output"
    python3 spectra2anndata_complex.py $data_dir perfbrain1-5_04m_nedc_14x_50um_recal Pub $output_dir    
    and added these lines:
    moo = os.path.expanduser("~/spatial_thesis/spatial_metabolomics/prototype/notebooks/hope.csv")  # TODO delete
    final_df = pd.read_csv(moo, sep='\t')  # TODO: delete
    to
    """
    coords_2d_arr = set_coords_2d_array(co_obj)
    # TODO: anndata save
    save__anndata_and_coord_dict(coords_2d_arr, final_df, file_basename,
                                 out_dir, barcodes_prefix=sys.argv[3])


if __name__ == "__main__":

    print("Converts imz+ibd files into anndata object,"
          "in the special case the mzs do not match, and saves it to file")

    print(" data_dir:", sys.argv[1], '\n',
          "base_filename:", sys.argv[2], '\n',
          "prefix:", sys.argv[3], '\n',
          "output_dir:", sys.argv[4])


    if 'objects' not in os.listdir(os.path.expanduser(sys.argv[4])):
        raise ValueError(f"objects/ subfolder in {sys.argv[4]} missing")



