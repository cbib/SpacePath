import os
import anndata as ad
import zonist.data as zd
import zonist.method as zm
import zonist.processing.measures_methods as zmm
from scripts.local_helpers.local_helpers import provide_annotation_dict


mydata_dir = "../../../data"
plot_parent_folder = "../../figures/"


pathwaysObj = zd.MetabolicPathways(
        data_dir=mydata_dir,
        subfolder="metabolites-mappings",
        file="WIKIPATHWAYS_human_multi_id_redundant.mmt"
)
pathwaysObj.molecule_id_type = "moleculeIds"
pathwaysObj.load()

SAMPLE = "UBcalibr4"
subfolder = SAMPLE

min_pathway_size = 5
print("\nSample: ", SAMPLE)
subfolder = SAMPLE
plot_folder = os.path.join(plot_parent_folder, SAMPLE)
annotations_to_use = ["LipidMaps", "HMDB", "CoreMetabolome"]
anno_dict = provide_annotation_dict(SAMPLE, annotations_to_use)

me_annots = zd.MetabolomeAnnotations(
    name=SAMPLE, data_dir=mydata_dir, subfolder=subfolder,
    annotation_files=anno_dict)

me_annots.usable_id = "moleculeIds"

me_annots.load()
adata_mz = ad.read_h5ad(filename=os.path.join(
    mydata_dir, subfolder, f'{SAMPLE}.h5ad'))

print(len(me_annots.no_translated['mz'].unique()))
print(len(adata_mz.var['mz'].unique()))


print(adata_mz.uns["history"])
print("MANDATORY: take away bad quality m/z")  # todo: how to make an alert ?
adata_mz = adata_mz[:,
           adata_mz.var.loc[adata_mz.var['final_quality'] == "good",
           "mz"].astype(str).tolist()]
print(len(adata_mz.var['mz'].unique()))

allgoodmz = adata_mz.var['mz'].astype(float).to_numpy()
me_annots, adata_mz = zmm.match_adata_mz_and_annots(me_annots, adata_mz)

me_annots.verify()
me_annots.keep_unique_mz_id_tuples(adata_mz.var.copy())  #
me_annots.remove_redundant_aliases()  # when using multitype 'moleculeIds' annotations

all_mz = me_annots.df['mz'].unique()

initial_pws_dict = pathwaysObj.to_dict()  # reserver pour plustard ? simplement avoir les tailles des sets du départ.

pws_obj_2, annots_2 = zm.SyncPathwaysAndAnnots.run(
    pathwaysObj.model_copy(),
    me_annots.model_copy(),
    min_pathway_size=min_pathway_size,
    max_pathway_size=500,   # todo: highlight that changed !!!
    mode="mz_aware")