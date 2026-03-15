# johaGL 2024
# Functions to transform imzML+ibd to anndata object
# NOTE: requires imzML+ibd spectra already aligned & recalibrated
#       if data is unaligned use 'matrix_from_metaspace' script
# Usage: within in python scripting
import gc
import logging
import os

import anndata as ad

import numpy as np

import pandas as pd

import pyimzml.ImzMLParser as imzmlp

from scipy.sparse import csr_matrix

from smpath.constants import mz_decimal_places

from tqdm import tqdm

logging.basicConfig(level=logging.INFO)


def open_spectra(imzml_fn, ibd_fn):
    """"credits: pyimzml by Alexandrov team"""
    logging.info("Reading imzML + ibd files with pyimzML package.")
    spectra_obj = imzmlp.ImzMLParser(
        filename=imzml_fn,
        parse_lib=None,
        ibd_file=open(ibd_fn, "rb"),
        include_spectra_metadata=None)
    return spectra_obj


def files2blocks(spectra_obj):
    """
    returns  features (np.array) , coordinates (np.array) and intensities_df.
    Note: intensities_df row names match perfectly
    -in content and order- to features list
    """
    mzs0, intensities = spectra_obj.getspectrum(0)
    coordinates_list = list()
    coord_str_list = list()
    intensities_list = list()

    unwarped_data_message = ("ERROR: unaligned features (please perform data "
                             "correction and warping before using this tool)")

    logging.info(
        "Assembling intensities arrays into 2D array. Parsing %s pixels",
        len(spectra_obj.coordinates))
    for i in tqdm(range(0, len(spectra_obj.coordinates))):  #
        mzs, intensities = spectra_obj.getspectrum(i)

        assert len(np.array(mzs)) == len(np.array(mzs0)), (
            unwarped_data_message)
        assert np.all(np.array(mzs) == np.array(mzs0)), (
            unwarped_data_message)

        try:
            if (intensities.shape[0] > 1) and (
                    intensities.shape[1] == len(mzs)):
                intensities = intensities[0]
        except IndexError:
            pass
        # only 2 dimensions:
        coord_str_list.append(
            f'x{spectra_obj.coordinates[i][0]}_y{spectra_obj.coordinates[i][1]}')
        coordinates_list.append([spectra_obj.coordinates[i][0],
                                 spectra_obj.coordinates[i][1]])
        intensities_list.append(intensities)
    # end for
    intensities_mx = np.array(intensities_list)
    intensities_df = pd.DataFrame(intensities_mx.T)
    del intensities_mx
    features = np.around(np.array(mzs, dtype=float),
                         decimals=mz_decimal_places)
    intensities_df.index = list(features)
    intensities_df.columns = coord_str_list

    return np.array(coordinates_list), features, intensities_df


def load_run(data_dir, base_filename, output_dir):
    data_dir = os.path.expanduser(data_dir)
    # file names
    imzml_fn = os.path.join(data_dir, f"{base_filename}.imzML")
    ibd_fn = os.path.join(data_dir, f"{base_filename}.ibd")

    spectra = open_spectra(imzml_fn, ibd_fn)
    coordinates_arr, features, intensities_df = files2blocks(spectra)

    del spectra
    gc.collect()
    adata = ad.AnnData(X=csr_matrix(intensities_df.T),
                       obs=pd.DataFrame(
                           {"obs": intensities_df.columns.tolist()},
                           index=intensities_df.columns.tolist()),
                       var=pd.DataFrame(
                           {"mz": features},
                           index=features
                       ),
                       obsm={"spatial": coordinates_arr})

    del intensities_df
    gc.collect()
    adata.uns['spatial'] = {
        'MY_SLICE7_NAME': {'scalefactors':
                               {'fiducial_diameter_fullres': 100,
                                'spot_diameter_fullres': 1,
                                'tissue_hires_scalef': 0.5,
                                'tissue_lowres_scalef': 0.5}}}

    outfile = os.path.join(output_dir, f'{base_filename}.h5ad')
    adata.write(outfile, compression=None, compression_opts=None)
    logging.info("done")


if __name__ == "__main__":
    # only for debugging with checkpoints in IDE
    load_run("rawdata/MY-DEMO-SAMPLE",
             "MY-DEMO-SAMPLE",
             "data/MY-DEMO-SAMPLE")
