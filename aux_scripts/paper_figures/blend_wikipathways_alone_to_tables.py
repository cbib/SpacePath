#########################################################################
# Very important requirement:
# first run the quality tagging over ions, which saves adequate .h5ad object:
#      SpacePath $input_file $data_dir $figures_dir --mode 1
# --
# johaGL 2025
#########################################################################

import os
import pandas as pd
import anndata as ad
import SpacePath.core.data as spd
import SpacePath.core.method as spm
import SpacePath.core.processing.measures_methods as measures_methods
from aux_scripts.data_download_prep_spots.spots_purge import \
    annotation_types_to_file_names_dict

from blend_keggids_keggpathways_to_tables import (
counts_per_orig_db, load_datasets_parameters
)


def initialize_WIKIPathways(data_dir):
    """
    Initialize and load metabolic pathways of wikipathways.
    """
    pathways_obj = spd.MetabolicPathways(
        data_dir=data_dir,
        subfolder="metabolites-mappings",
        file="WIKIPATHWAYS_human_multi_id_redundant.mmt"
    )
    pathways_obj.molecule_id_type = 'moleculeIds'  # define before loading
    pathways_obj.load()

    return pathways_obj


def process_sample_wiki(SAMPLE, my_new_dict, pathwaysObj, annots_obj,
                          min_pathway_size,
                          name_sample, mz_mapped_to_pathways,
                          compound_ids_to_mz_ratio_after_pathway_mapped,
                          how_many_pathways):

    pw1, ann1_df = spm.SyncPathwaysAndAnnots.run(
        pathwaysObj.model_copy(),
        annots_obj.df.copy(),
        usable_id="moleculeIds",
        min_set_size=1,
        max_set_size=200)

    my_new_dict[SAMPLE] += counts_per_orig_db(ann1_df)

    pws_obj_2, annots_2_df = spm.SyncPathwaysAndAnnots.run(
        pathwaysObj.model_copy(),
        annots_obj.df.copy(),
        usable_id="moleculeIds",
        min_set_size=min_pathway_size,
        max_set_size=200,
        mode="mz_aware"
    )

    name_sample.append(SAMPLE)
    mz_ok = annots_2_df['mz'].unique()
    theIDok = annots_2_df[annots_obj.usable_id].unique()
    mz_mapped_to_pathways.append(len(mz_ok))
    try:
        compound_ids_to_mz_ratio_after_pathway_mapped.append(
            len(set(theIDok)) / len(set(mz_ok))
        )
    except ZeroDivisionError:
        compound_ids_to_mz_ratio_after_pathway_mapped.append(0)

    how_many_pathways.append(len(pws_obj_2.piled['pathway_id'].unique()))

    return (my_new_dict, name_sample, mz_mapped_to_pathways,
                   compound_ids_to_mz_ratio_after_pathway_mapped,
                   how_many_pathways)


def save_wiki_mappings_results(my_new_dict, selection, name_sample,
                               mz_mapped_to_pathways,
                               compound_ids_to_mz_ratio_after_pathway_mapped,
                               how_many_pathways):
    out_df1 = pd.DataFrame.from_dict(
        my_new_dict, orient='index')

    out_df1.reset_index(inplace=True)
    out_df1.columns = ["sample_name", "HMDB identifiers",
                       "LipidMaps identifiers", "KEGG identifiers",
                       "HMDB ids in pathways", "LipidMaps ids in pathways",
                       "KEGG ids in pathways"]

    out_df1.to_csv(
        f"../../figs/mapp-IDSmix-pathways-{selection}.csv")

    out_df2 = pd.DataFrame({
        "sample_name": name_sample,
        "mz_mapped_to_pathways": mz_mapped_to_pathways,
        "ratio_id_mz": compound_ids_to_mz_ratio_after_pathway_mapped,
        "pathways_number_final": how_many_pathways})
    out_df2.to_csv(
        f"../../figs/crit-IDSmix-counts-{selection}.csv")


def main(selection):
    print("\nRunning a general quantitative description of the WIKIPATHWAYS"
          " mappings across specified datasets. Output: .csv files with counts")

    data_centralized_folder = "../../../data"

    datasets_parameters_df = load_datasets_parameters(selection)

    pathwaysObj = initialize_WIKIPathways(data_centralized_folder)

    name_sample = []
    mz_mapped_to_pathways = []
    compound_ids_to_mz_ratio_after_pathway_mapped = []
    how_many_pathways = []
    my_new_dict = dict()

    assert os.path.exists("../../figs"), "Error, ../../figs must exist"

    for i, row in datasets_parameters_df.iterrows():
        SAMPLE = row["dataset_name"]
        min_pathway_size = row["min_pathway_size"]
        print("\nSample: ", SAMPLE)
        subfolder = SAMPLE
        annotations_to_use = row["annotations_to_use"].split(" ")
        annotations_dict = annotation_types_to_file_names_dict(SAMPLE,
                                                               annotations_to_use)

        annots_obj = spd.MetabolomeAnnotations(
            name=SAMPLE, data_dir=data_centralized_folder,
            subfolder=subfolder,
            annotation_files=annotations_dict)

        annots_obj.usable_id = "moleculeIds"
        annots_obj.load()

        adata_mz = ad.read_h5ad(filename=os.path.join(
            data_centralized_folder, subfolder, f'{SAMPLE}.h5ad'))

        adata_mz = adata_mz[:,
        adata_mz.var.loc[adata_mz.var['final_quality'] == "good",
        "mz"].astype(str).tolist()]

        annots_obj, adata_mz = measures_methods.match_adata_mz_and_annots(
            annots_obj, adata_mz)

        my_new_dict[SAMPLE] = counts_per_orig_db(annots_obj.no_translated)

        annots_obj.verify()
        annots_obj.keep_best_mz_id_pairs(adata_mz.var.copy())
        annots_obj.remove_redundant_aliases()  # because using multitype 'moleculeIds' annotations

        (my_new_dict, name_sample, mz_mapped_to_pathways,
                   compound_ids_to_mz_ratio_after_pathway_mapped,
                   how_many_pathways) = process_sample_wiki(
            SAMPLE, my_new_dict, pathwaysObj, annots_obj,
            min_pathway_size,
            name_sample, mz_mapped_to_pathways,
            compound_ids_to_mz_ratio_after_pathway_mapped,
            how_many_pathways
        )
    # end for

    save_wiki_mappings_results(my_new_dict, selection, name_sample,
                           mz_mapped_to_pathways,
                           compound_ids_to_mz_ratio_after_pathway_mapped,
                           how_many_pathways)


if __name__ == "__main__":
    main(selection="public")  # Change to "private" as needed






