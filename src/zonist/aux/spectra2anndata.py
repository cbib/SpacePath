# johaGL 2024
# NOTE: requires imzML+ibd spectra already aligned & recalibrated
#       if data is unaligned use 'matrix_from_metaspace' script
# usage:
#   spectra2anndata-v2.py $DATADIR $FILE_BASENAME $PREFIX $OUTDIR
# where:
#   $FILEBASENAME is the file name without extension
#       Note: both files .imzML and .ibd must have same FILEBASENAME
#   $OUTDIR is the relative path to output folder and must contain 'objects/'
#       Note: outdir must exist and be one level below working dir

# example (example data):
#  cd $my_working_dir
#  python spectra2anndata.py \
#  ~/spatial_thesis/spatial_metabolomics/tutorial_publicdata/example_data \
#  example_data ED output/objects

# TODO: transfomr this in ----> spectra2matrix -> and save a .h5ad file !!!!!

import numpy as np
import pandas as pd
import anndata as ad
import pyimzml.ImzMLParser as imzmlp
from tqdm import tqdm
import os
import sys
from scipy.sparse import csr_matrix
from zonist.constants import mz_decimal_places


def open_spectra(imzml_fn, ibd_fn):
    """"credits: pyimzml by Alexandrov team"""
    print("Reading imzML + ibd files with pyimzML package.")
    spectra_obj = imzmlp.ImzMLParser(filename=imzml_fn,
            parse_lib=None,
            ibd_file=open(ibd_fn, "rb"),
            include_spectra_metadata=None)
    return spectra_obj


def files2blocks(spectra_obj):
    """returns  features (np.array) , coordinates (np.array) and intensities_df.
    Note: intensities_df row names match perfectly -in content and order- to features list"""
    mzs0, intensities = spectra_obj.getspectrum(0)
    coordinates_list = list()
    coord_str_list = list()
    intensities_list = list()

    unwarped_data_message = ("ERROR: unaligned features (please perform data "
                             "correction and warping before using this tool)")

    print("Assembling intensities arrays into 2D array. Parsing",
          len(spectra_obj.coordinates), "pixels")
    for i in tqdm(range(0, len(spectra_obj.coordinates) )):  #
        mzs, intensities = spectra_obj.getspectrum(i)

        assert len(np.array(mzs)) == len(np.array(mzs0)), (
            unwarped_data_message)
        assert np.all(np.array(mzs) == np.array(mzs0)), (
            unwarped_data_message)

        try:
            if (intensities.shape[0] > 1) and (intensities.shape[1] == len(mzs)):
                intensities = intensities[0]
        except IndexError:
            pass
        # only 2 dimensions:
        coord_str_list.append(
            f'x{spectra_obj.coordinates[i][0]}_y{spectra_obj.coordinates[i][1]}')
        coordinates_list.append([spectra_obj.coordinates[i][0],
                                 spectra_obj.coordinates[i][1]])
        intensities_list.append(intensities) 
        # intensities_mx = np.vstack((intensities_mx, intensities))  # bad because makes multiple recomputation+realocation in memory of size-changing matrix
    # end for
    intensities_mx = np.array(intensities_list)
    intensities_df = pd.DataFrame(intensities_mx.T)
    intensities_df.index = list(np.around(np.array(mzs, dtype=float),
                                          decimals=mz_decimal_places))
    intensities_df.columns = coord_str_list
           
    return np.array(coordinates_list), np.array(mzs), intensities_df


def load_run(data_dir, base_filename, output_dir):
    data_dir = os.path.expanduser(data_dir)
    # file names
    imzml_fn = os.path.join(data_dir, f"{base_filename}.imzML")
    ibd_fn = os.path.join(data_dir, f"{base_filename}.ibd")

    spectra = open_spectra(imzml_fn, ibd_fn)
    coordinates_arr, features, intensities_df = files2blocks(
        spectra
    )
    print(intensities_df.head(3))

    adata = ad.AnnData(X = csr_matrix(intensities_df.T),
                       obs=pd.DataFrame(
                           {"obs": intensities_df.columns.tolist()},
                           index=intensities_df.columns.tolist()),
                       var=pd.DataFrame(
                           {"mz": np.array(features, dtype=np.float32)},
                           index=list(features)
                       ),
                       obsm={"spatial": coordinates_arr} )
    adata.uns['spatial'] = {
        'MY_SLICE7_NAME': {'scalefactors':
                               {'fiducial_diameter_fullres': 100,
                                'spot_diameter_fullres': 1,
                                'tissue_hires_scalef': 0.5,
                                'tissue_lowres_scalef': 0.5}}}

    outfile = base_filename + ".h5ad"
    adata.write(outfile, compression=None, compression_opts=None)
    print(f"End: file saved to '{outfile}'")


if __name__ == "__main__":
    print("Converts imz+ibd files into anndata object, and saves it to file")
    print("   important: requires aligned and recalibrated spectra (otherwise, use 'matrix_from_metaspace')")
    print(" data_dir:", sys.argv[1], '\n',
          "base_filename:", sys.argv[2], '\n',
          "output_dir:", sys.argv[3])

    load_run(data_dir=sys.argv[1], base_filename=sys.argv[2],
              output_dir=sys.argv[3])

    #load_run("data/imzML","Ctrl_NEDC_NEG", "output/objects")
