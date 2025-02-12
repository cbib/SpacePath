import sys
import os
import numpy as np
import pandas as pd
import anndata as ad
import zonist.data as zd
from scipy.sparse import csr_matrix
from sklearn.utils import sparsefuncs
import squidpy as sq

from scripts.bank4zonist.pathways_methods import compute_pathways_scores

sys.path.append(os.path.join(os.path.dirname(__file__),
                             '../bank4zonist'))
from quantifs_methods import *
from datasetcomplement import (build_kegg_resolved_adata,
                               yield_specific_typical_pathways
                               )
from pathways_methods import compute_kpca, reshape_pathways_for_kpca


"""
Ravi & Heiland data
use:
   python $MYPY $SAMPLE
where :
  $MYPY is the script file with .py extension (extension included)
  sample one  of: 
    samples = [   "248_T"   , "259_T", "260_T", "262_T", "269_T", "275_T"]

Or use Bash oneliner:
mylist=(248_T 259_T 260_T 262_T 269_T 275_T);for k in ${mylist[*]}; do echo $k; python $MYPY $k; done

"""

SAMPLE = "248_T" #sys.argv[1] #  "248_T"

mydata_dir = "/home/johanna/spatial_thesis/spatial_metabolomics/data"
#"/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"
anno_dict = {'HMDB': "Metaspace_annotations_943peaks.csv"}
plot_folder = "../../../figures/heiland"

assert os.path.exists(plot_folder), "Cannot continue, the plot_folder does not exist"

gathered_annots = zd.MetabolomeAnnotations(name="heiland",
                            data_dir=mydata_dir,
                            subfolder="heiland-msi",
                            annotation_files=anno_dict)
gathered_annots.load()

adata_mz = ad.read_h5ad(filename=os.path.join(
    mydata_dir, "heiland-msi", f'MALDI_{SAMPLE}.h5ad'))

print(adata_mz)
adata_mz = get_rid_of_poor_var_mz(adata_mz)  # TODO : highlight this because is proper to this Ravi&Heiland data
print("now:", adata_mz)

me_annots = update_annotation_tomatch_adatamz(gathered_annots, adata_mz)  # TODO highlight: proper to Ravi&Heiland data as annotation is 1 file for all 6 slices
verification(adata_mz, me_annots)

## reference mappings to get KEGG ids from HMDB and LipidMaps identifiers
me_refmap = zd.MetaboliteIdsMappings(
    data_dir=mydata_dir,
    subfolder="metabolites-mappings",
    mapping_files={'HMDB': "HMDB-to-KEGG_05-07-24.tsv",
                   'LipidMaps': "LipidMaps-to-KEGG_10-18-24.tsv" })

me_refmap.load()
me_annots.convert_ids_to_kegg(me_refmap.table)
me_annots.treat_mz_multiplets(adata_mz)

# print(len(me_annots.df["KEGG"].to_list())); print(len(me_annots.df["KEGG"].unique()))

## pathways
pathwaysObj = zd.MetabolicPathways(
    data_dir=mydata_dir,
    subfolder="metabolites-mappings",
    official_file="KEGG_hsa_pathways_compounds_R110.gmt")

pathwaysObj.load()
pathwaysObj.intersect_with_annotations(me_annots)

print("subset adata_mz for only mz mapping to a KEGG metabolite identifiers that in turn map to pathways")
me_annots_OKdf = me_annots.df.loc[
              me_annots.df["KEGG"].isin(pathwaysObj.piled["KEGG"].unique()), :]

adata_mz = adata_mz[:, me_annots_OKdf["mz"].astype(str).unique()].copy()

# the KEGG id adata
adata = build_kegg_resolved_adata(adata_mz, annots_df=me_annots_OKdf)
print("KEGG resolved adata: ", adata)

###
#  select pathways to visualize:
#
##
biologically_chosen_pathways = yield_specific_typical_pathways()

###
# compute kpca
##
kpca_scores = compute_pathways_scores(pathwaysObj, adata,)

print(kpca_scores.iloc[:3,:5])

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

plot_selected_pathways(kpca_score_adata, biologically_chosen_pathways,
                  f"{SAMPLE}--Picked", plot_folder)

kpca_score_adata = spatial_standard_computations(kpca_score_adata, pca_ncomps=12,
                                  plot_folder=plot_folder,
                              prefix=f'{SAMPLE}PathwaysScores')

kpca_score_adata = compute_moran_i(kpca_score_adata, 5)

spatially_top_pathways = kpca_score_adata.uns['moranI'].iloc[:6, :].index.tolist()

df = kpca_score_adata.uns['moranI'].copy()
tmp2 = pathwaysObj.pathways_ids[["pathway_id", "pathway_name"]].drop_duplicates()
these_names = []
for theid in df.index.tolist():
    these_names.append(tmp2.loc[
                           tmp2["pathway_id"] == theid, "pathway_name"].item())
df.insert(0, 'Pathway_name', these_names)
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

print("ok ")
