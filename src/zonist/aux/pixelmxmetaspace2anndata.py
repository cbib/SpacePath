# johaGL 2024

# opening the pixel intensities .csv files (downloaded from metaspace result)
#
# .py $DATADIR $KEGGpixels $Lipidmaps_pixels $PREFIX $OUTDIR

# TODO: reuse chunks as needed and DELETE ALL this entire script !!!!!!!



import numpy as np
import pandas as pd
import anndata as ad
from typing import List
from tqdm import tqdm
import os
import sys
import pickle


def generate_coordinates(pixels_names):
    coord_l = list()
    for k in pixels_names:   # eg 'x0_y0'
        items = k.split("_")
        coord_l.append(np.array([int(items[0].replace("x", "")),
                       int(items[1].replace("y", ""))]))
    return np.array(coord_l)


def get_new_names_pixels(coordinates_arr, prefix):
    ol = list()
    for tup in coordinates_arr:
        ol.append(f'{prefix}_{tup[0]}_{tup[1]}')
    return ol


def best_1D_arr(the_2D_arr):
    sums = the_2D_arr.sum(axis=1)
    best_arg = np.argmax(sums)
    return the_2D_arr[best_arg]


def merge_pixels_to_df(
        kegg_px: pd.DataFrame, lipi_px: pd.DataFrame,
        all_mz_values: List[float], pixels_names: List[str]) -> pd.DataFrame:
    array_intensities = np.zeros(kegg_px.iloc[:, 5:].shape[1])
    for mz in tqdm(all_mz_values):
        if mz in kegg_px['mz'].tolist():
            tmp = kegg_px.loc[kegg_px['mz'] == mz, pixels_names].to_numpy()
        else:
            tmp = lipi_px.loc[lipi_px['mz'] == mz, pixels_names].to_numpy()
        if tmp.shape[0] > 1:
            tmp = best_1D_arr(tmp)
        array_intensities = np.vstack((array_intensities, tmp))
    intensities_df = pd.DataFrame(
        array_intensities[1:])  # exclude the 0, was np.zeros
    return intensities_df


def load_run(data_dir, kegg_pixels, lipidmaps_pixels,
             prefix, output_dir):
    kegg_pxf = os.path.join(data_dir, kegg_pixels)
    lipi_pxf = os.path.join(data_dir, lipidmaps_pixels )
    kegg_px = pd.read_csv(kegg_pxf, comment='#', header=0)
    lipi_px = pd.read_csv(lipi_pxf, comment='#', header=0)

    all_mz_values = list(set(kegg_px['mz'].astype(float).to_numpy()).union(
        set(lipi_px['mz'].astype(float).to_numpy())))
    # print(np.all(kegg_px.iloc[:, 5:].columns.to_numpy() == lipi_px.iloc[
    # :,5:].columns.to_numpy()))
    pixels_names = kegg_px.iloc[:, 5:].columns.tolist()
    intensities_df = merge_pixels_to_df(kegg_px, lipi_px,
                                        all_mz_values, pixels_names)
    coordinates_arr = generate_coordinates(pixels_names)
    #  Length mismatch: Expected axis has 3857 elements, new values have 3827 elements
    intensities_df.index = [str(i) for i in all_mz_values]
    intensities_df.columns = get_new_names_pixels(coordinates_arr, prefix)
    adata = ad.AnnData(intensities_df.T,
                       obsm={"spatial": np.array(coordinates_arr)},
                       varm={"transcripts": np.array(list(all_mz_values))})
    outfile = os.path.join(output_dir, f"{prefix}--anndata.pk")
    with open(outfile, "wb") as fopen:
        pickle.dump(adata, fopen)
    print(f"saved anndata object (pickle) to file:  {outfile}")


if __name__ == "__main__":
    print("Accepts Pixels intensities (that you obtained from metaspace result). Saves anndata (pickle)")

    print(" data_dir:", sys.argv[1], '\n',
          "kegg_pixel_df:", sys.argv[2], '\n',
          "lipidmaps_pixel_df:", sys.argv[3], '\n',
          "prefix:", sys.argv[4], '\n',
          "output_dir:", sys.argv[5])

    load_run(data_dir=sys.argv[1], kegg_pixels=sys.argv[2],
             lipidmaps_pixels=sys.argv[3],
             prefix=sys.argv[4], output_dir=sys.argv[5])

    # data_dir = "data/"
    # kegg_pixels = "pixel_mx_metaspace/perfbrain-public_pixel_intensities-KEGG.csv"  #uno.csv"
    # lipidmaps_pixels = "pixel_mx_metaspace/perfbrain-public_pixel_intensities-Lipidmaps.csv" #dos.csv"
    # prefix = "tpub"
    # output_dir = 'output/objects'
    #load_run(data_dir, kegg_pixels, lipidmaps_pixels, prefix, output_dir)


