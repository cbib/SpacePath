import os
import numpy as np
import anndata as ad
import pandas as pd
import scanpy as sc
from typing import List
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from sklearn.utils import sparsefuncs
import squidpy as sq

"""
get the matrix of KEGG given the annotations, and the adata_mz
"""



def yield_specific_typical_pathways() -> dict:
    biologically_chosen_pathways = {
        "hsa00600": "Sphingolipid metabolism",
        "hsa00030": "Pentose phosphate pathway",
        "hsa00010": "Glycolysis Gluconeogenesis", # attention, name here has been edited to avoid "/" in the middle
        "hsa00620" : "Pyruvate metabolism",
        "hsa00020"  : "TCA cycle",
        "hsa00900": "Terpenoid backbone biosynthesis",
        "hsa04727": "GABAergic synapse"}
    return biologically_chosen_pathways


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
                np.array(local_1D_adata_mz.X.todense().T).flatten() / len(local_kegg))
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


def visualizing_spatial_domains_before_decision(adata, leiden_resol = 0.4,
                                                alpha = 0.6  # weight for the spatial (coords) nngraph
):
    # TODO: see Warning /home/johanna/spatial_thesis/spatial_metabolomics/user-pipeline/scripts/heiland-ravi/datasetcomplement.py:58: FutureWarning: In the future, the default backend for leiden will be igraph instead of leidenalg.
    # TODO To achieve the future defaults please pass: flavor="igraph" and n_iterations=2.  directed must also be False to work with igraph's implementation.
    print("visualizing the spatial domains before decision about out-of-tissue pixels...")
    #  use Spatial domains as in https://www.sc-best-practices.org/spatial/domains.html
    means_vars, variance_computed = sparsefuncs.mean_variance_axis(adata.X, axis=0)
    adata.var['variance'] = variance_computed
    tmp = adata.var.sort_values(by='variance', ascending=False)
    sampled_features = tmp.head(500)["mz"].astype(str).tolist()
    miniadata = adata[:, sampled_features].copy()
    sc.pp.pca(miniadata)
    sc.pp.neighbors(miniadata)  # nearest neighbor graph
    sc.tl.leiden(miniadata, resolution=leiden_resol, key_added=f"leiden{leiden_resol}")
    nn_graph_genes = miniadata.obsp["connectivities"]
    sq.gr.spatial_neighbors(miniadata)    # spatial proximity graph
    nn_graph_space = miniadata.obsp["spatial_connectivities"]
    print("computing leiden clustering on joint graph...")

    joint_graph = (1 - alpha) * nn_graph_genes + alpha * nn_graph_space
    sc.tl.leiden(miniadata, resolution=leiden_resol, adjacency=joint_graph,
                 key_added=f"squidpy_domains(leiden{leiden_resol})")
    sq.pl.spatial_scatter(
        miniadata, color=[f"squidpy_domains(leiden{leiden_resol})"],
        img=False, title=f"before filtering out-of-tissue pixels. (leiden{leiden_resol})")
    plt.show()
    return miniadata.obs


def take_decision(adata, computed_domains_df,
                  not_tissue_spatial_domains: List[int],
                  used_leiden_resol = 0.4):
    """e.g.  not_tissue_spatial_domains = [2,7,8]   # DECISION (requires viz first)
    """
    print("\nGet rid of pixels that are not part of the tissue (this must be after visualizing):")
    not_tissue_spatial_domains = [str(x) for x in not_tissue_spatial_domains]
    pixels_to_keep = computed_domains_df.loc[
        ~computed_domains_df[f"squidpy_domains(leiden{used_leiden_resol})"].isin(
            not_tissue_spatial_domains), "obs"].tolist()
    adata4viz = adata.copy()
    adata4viz.obs["spatial_domain"] = computed_domains_df[f"squidpy_domains(leiden{used_leiden_resol})"].tolist()
    adata4viz = adata4viz[pixels_to_keep, :]
    print("visualizing the spatial plot of domains after decision")
    tmp = f"AFTER filtering out-of-tissue pixels (leiden{used_leiden_resol})"
    sq.pl.spatial_scatter(adata4viz, color=["spatial_domain"],
                          img=False, title=tmp)
    plt.show()

    return adata[pixels_to_keep, :]


