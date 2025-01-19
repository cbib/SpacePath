import os
import numpy as np
import anndata as ad
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix

"""
get the matrix of KEGG given the annotations, and the adata_mz
"""

def build_kegg_resolved_adata(adata_mz, annots_df):
    """
    Given the adata_mz and the annots_df, returns a new adata object
    """
    mz_list = annots_df["mz"].unique()
    intensities_arrs_list = []
    kegg_list = list()
    for one_mz in mz_list:
        local_kegg = annots_df.loc[
            annots_df["mz"] == one_mz, "KEGG"].tolist()
        local_1D_adata_mz = adata_mz[ :, str(one_mz)].copy()
        for a_kegg in local_kegg:
            intensities_arrs_list.append(
                np.array(local_1D_adata_mz.X.todense().T).flatten() / float(
                    len(local_kegg)))
            kegg_list.append(a_kegg)
    adata = ad.AnnData(X=csr_matrix(np.array(intensities_arrs_list).T),
                       obs=adata_mz.obs,
                       var=pd.DataFrame({"ID": kegg_list},
                                        index=kegg_list),
                       obsm=adata_mz.obsm)
    adata.uns['spatial'] = {'my_slice_name': {
        'scalefactors': {'fiducial_diameter_fullres': 100,
                         'spot_diameter_fullres': 1,
                         'tissue_hires_scalef': 0.5,
                         'tissue_lowres_scalef': 0.5}}}

    return adata

""""""