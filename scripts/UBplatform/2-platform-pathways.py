import sys
import os
import numpy as np
import pandas as pd
import anndata as ad
import zonist.data as zd
import zonist.method as zm
from scipy.sparse import csr_matrix
from sklearn.utils import sparsefuncs
import squidpy as sq


import sspa as sspa


sys.path.append(os.path.join(os.path.dirname(__file__),
                             '../bank4zonist'))
from quantifs_methods import *
from datasetcomplement import (yield_specific_typical_pathways)


"""
pathways scores on *bordeaux proteome example data* 

"""

SAMPLE = "UBcalibr2"#sys.argv[1]
#"UBcalibr2"
print(SAMPLE)
mydata_dir = "/home/johanna/spatial_thesis/spatial_metabolomics/data"
#"/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"
subfolder = "UBplatform"
#"/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"
anno_dict = {'HMDB': f"{SAMPLE}_hmdb.csv",
             "LipidMaps": f"{SAMPLE}_lipidmaps.csv"}
plot_folder = "../../../figures/platform"


assert os.path.exists(plot_folder), "Cannot continue, the plot_folder does not exist"

me_annots = zd.MetabolomeAnnotations(name=SAMPLE,
                            data_dir=mydata_dir,
                            subfolder=subfolder,
                            annotation_files=anno_dict)
me_annots.load()


adata_mz = ad.read_h5ad(filename=os.path.join(
    mydata_dir, subfolder, f'{SAMPLE}_mx.h5ad'))

verification(adata_mz, me_annots)

## reference mappings to get KEGG ids from HMDB and LipidMaps identifiers
me_refmap = zd.MetaboliteIdsMappings(
    data_dir=mydata_dir,
    subfolder="metabolites-mappings",
    mapping_files={'HMDB': "HMDB-to-KEGG_05-07-24.tsv",
                   'LipidMaps': "LipidMaps-to-KEGG_10-18-24.tsv" }
)
me_refmap.load()
me_annots.convert_ids_to_kegg(me_refmap.table)
me_annots.treat_mz_multiplets(adata_mz)

## pathways
pathwaysObj = zd.MetabolicPathways(
    data_dir=mydata_dir,
    subfolder="metabolites-mappings",
    official_file="KEGG_hsa_pathways_compounds_R110.gmt")

pathwaysObj.load()

pathwaysObj, me_annots = zm.SyncPathwaysAndAnnots.run(pathwaysObj,
                                                      me_annots,
                                                      min_pathway_size=2)  # default max is 500

adata_mz = adata_mz[:, me_annots.df["mz"].astype(str).unique()].copy()

# the KEGG id adata
adata = build_kegg_resolved_adata(adata_mz, annots_df=me_annots.df)
print("KEGG resolved adata: ", adata)


# compute kpca #
kpca_scores = zm.ComputePathwaysScores.run(pathwaysObj, adata)
print(kpca_scores.shape)
print(kpca_scores.iloc[:3,:5])


# #

do_plots = True  # # several plots :
if do_plots:
    biologically_chosen_pathways = yield_specific_typical_pathways()

    # TODO: is it always needed to shift the distr to the right ???
    kpca_scores = (kpca_scores + abs(kpca_scores.min().min())) * 1e3 # shift to right and unshrink

    kpca_score_adata = ad.AnnData(
        X=csr_matrix(kpca_scores),
        var=pd.DataFrame(
            {'pathway_id' : kpca_scores.columns.tolist()},
            index=kpca_scores.columns.tolist()),
        obs=adata.obs,
        obsm=adata.obsm,
        uns=adata.uns
    )

    # plot_selected_pathways(kpca_score_adata, biologically_chosen_pathways,
    #                  f"{SAMPLE}--Picked", plot_folder)

    kpca_score_adata = spatial_standard_computations(kpca_score_adata, pca_ncomps=12,
                                      plot_folder=plot_folder,
                                  prefix=f'{SAMPLE}PathwaysScoresV2')


    kpca_score_adata = compute_moran_i(kpca_score_adata, 5)

    spatially_top_pathways = kpca_score_adata.uns['moranI'].iloc[:6, :].index.tolist()

    df = kpca_score_adata.uns['moranI'].copy()
    tmp2 = pathwaysObj.pathways_ids[["pathway_id", "pathway_name"]].drop_duplicates()
    these_names = []
    for theid in df.index.tolist():
        these_names.append(tmp2.loc[
                               tmp2["pathway_id"] == theid, "pathway_name"].item())
    df.insert(0, 'Pathway_name', these_names)
    df.to_csv(os.path.join(plot_folder, f'{SAMPLE}-PathwaysMoranI.tsv'),
              header=True, index=True, sep="\t")
    del these_names, tmp2, df

    for m_id in spatially_top_pathways:
        pw_name = pathwaysObj.pathways_ids.loc[
            pathwaysObj.pathways_ids['pathway_id'] == m_id, "pathway_name"].item()

        tmpstr = os.path.join(plot_folder, f'{SAMPLE}-{pw_name}.png')  # TODO: improve destination folder
        try:
            sq.pl.spatial_scatter(kpca_score_adata, color=[m_id],
                                  img=False, title=f"{m_id} {pw_name}",
                                  save=os.path.abspath(tmpstr))
        except Exception as e:
            print(e)
            print(f"impossible to plot {m_id} {pw_name}")

print("the end")







