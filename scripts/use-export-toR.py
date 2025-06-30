import os
import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import squidpy as sq
from local_helpers.quantifs_methods import coordinates_as_df


""""
"""
if __name__ == "__main__":

   # define here my files and dirs

    mydata_dir = "../../data"
    plot_parent_folder = "../../figures/"

    r_test_folder = "../../cardinal_microtest/data_for_cardinal"

    # start
    dict_opts = {0: 'TD-Lipido-635',
                 1: 'LP-brain-mouse-5'}

    for kopt in dict_opts.keys():
        SAMPLE = dict_opts[kopt]
        print("\nSample: ", SAMPLE)
        subfolder = SAMPLE
        plot_folder = os.path.join(plot_parent_folder, SAMPLE)


        adata_mz = ad.read_h5ad(filename=os.path.join(
                    mydata_dir, subfolder, f'{SAMPLE}.h5ad'))

        print(adata_mz.uns["history"])
        print("MANDATORY: take away bad quality m/z") # todo: how to make an alert ?
        adata_mz = adata_mz[:,
                   adata_mz.var.loc[adata_mz.var['final_quality'] == "good",
                   "mz"].astype(str).tolist()]
        try:
            os.mkdir(os.path.join(r_test_folder, SAMPLE))
        except:
            pass

        out_coords = os.path.join(r_test_folder,
                              SAMPLE, f"{SAMPLE}-coords.tsv")
        out_vars = os.path.join(r_test_folder,
                                  SAMPLE, f"{SAMPLE}-vars.tsv")
        out_matrix = os.path.join(r_test_folder,
                              SAMPLE, f"{SAMPLE}-matrix.tsv")

        # seeing if sincerely grateful is a good phrase ?

        coords_df = coordinates_as_df(adata_mz)
        coords_df.to_csv(out_coords, sep='\t', index=True, header=True)

        adata_mz.var.to_csv(out_vars,
                            sep='\t', index=True, header=True
                            )
        adata_mz.to_df().to_csv(out_matrix,  sep='\t',  index=True, header=True )


