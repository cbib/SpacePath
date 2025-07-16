import os
import numpy as np
import pandas as pd
import anndata as ad
import zonist.data as zd
from scipy.sparse import csr_matrix

"""
This step is completely custom as the input consists of .tsv files, 
see explanation.txt file.

 Tasks:
 1) Reassign mz indexes to MALDI data

   The data in .tsv files were produces by R scripts of Heiland team (Freiburg):
    a.corrected batch effect
    b.repeated the rows in the quantifications data to accomodate to the fact that several hmdb correspond to one mz
    c.assigned HMDB indexes to data
    So this script re-assigns mz to quantifications data and deduplicates the rows

 2) Save the resulting adata objects

"""


def format_properly_df(df):
    df = df.sort_values(by='mz', ascending=True)
    df['mz'] = np.around(df['mz'].astype(float).to_numpy(), decimals=6)
    df.index.name = "hmdb"
    df = df.reset_index()
    df = df.drop(columns=["hmdb"])
    df = df.drop_duplicates()
    mz_arr = df['mz'].astype(float).to_numpy()
    df = df.set_index('mz')
    return df, mz_arr


def open_coords_file(sample):
    mycoords1 = pd.read_csv(
        os.path.join(data_centralized_folder, "heiland-msi", "ancient",
                     f"coords_{sample}.tsv"),
        sep="\t", comment="#", header=0, index_col=0)

    mycoords1 = mycoords1.set_index("barcodes")

    return mycoords1


def as_adata_custom(df, mz_arr, coords_df):
    df = df.T  # transpose is compulsory (AnnData requires)
    adata = ad.AnnData(
        X=csr_matrix(df),
        obs=pd.DataFrame({"obs": df.index.tolist()},
                         index=df.index.tolist()),
        var=pd.DataFrame({"mz": mz_arr},
                         index=df.columns.tolist()),
        obsm={"spatial": np.array(coords_df[["x", "y"]])}
    )
    # set pixels plot sizes
    adata.uns['spatial'] = {'my_slice_name': {
        'scalefactors': {'fiducial_diameter_fullres': 100,
                         'spot_diameter_fullres': 1,
                         'tissue_hires_scalef': 0.5,
                         'tissue_lowres_scalef': 0.5}}}
    return adata


def shift_all_matrices_to_right_and_save(dict_of_dfs, dict_of_mz_arrs,
                                         min_value, out_dir):
    """shifts to avoid negative floats
    and saves adata objects files"""
    print("shifting all matrices to the right")
    for k in dict_of_dfs.keys():
        df = dict_of_dfs[k]
        mz_arr = dict_of_mz_arrs[k]
        df = df + abs(min_value)
        coords_df = open_coords_file(k)
        adata = as_adata_custom(df, mz_arr, coords_df)
        subfolder = f"MALDI_{k}"
        file_name = f'MALDI_{k}.h5ad'
        out_path = os.path.join(out_dir, subfolder, file_name)
        adata.write_h5ad(filename=out_path)


## start

samples = [
    "248_T", "259_T", "260_T", "262_T", "269_T", "275_T"
]

data_centralized_folder = "/home/johanna/spatial_thesis/spatial_metabolomics/data"
# "/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"
annotations_dict = {'HMDB': "Metaspace_annotations_943peaks.csv"}

annots_obj = zd.MetabolomeAnnotations(name="heiland",
                                     data_dir=data_centralized_folder,
                                     subfolder="heiland-msi",
                                     annotation_files=annotations_dict)
annots_obj.load()

cool = annots_obj.no_translated[["mz",
                             "moleculeIds"]].drop_duplicates()  # in HMDB, mz and moleculeIds are 1:1
hmdb_mz_dict = dict()
for i, r in cool.iterrows():
    hmdb_mz_dict[r["moleculeIds"]] = r["mz"]

dict_of_dfs = dict()
dict_of_mz_arrs = dict()
print("finding the min value of all the matrices")
min_value = 999  # arbitrary initial value 999
for k in samples:
    df = pd.read_csv(
        os.path.join(data_centralized_folder, "heiland-msi",
                     "ancient",
                     f'MALDI_{k}.tsv'),
        sep="\t", comment="#", header=0, index_col=0)
    df.insert(0, 'mz', df.index.map(hmdb_mz_dict.get))
    df, mz_arr = format_properly_df(df)
    if df.min().min() < min_value:
        min_value = df.min().min()
    dict_of_mz_arrs[k] = mz_arr
    dict_of_dfs[k] = df



shift_all_matrices_to_right_and_save(dict_of_dfs, dict_of_mz_arrs, min_value,
                                    out_dir=data_centralized_folder)

print("done")

