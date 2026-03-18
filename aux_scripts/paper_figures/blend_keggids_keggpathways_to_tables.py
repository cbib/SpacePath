#########################################################################
# Very important requirement:
# first run the quality tagging over ions, which saves adequate .h5ad object:
#     SpacePath $input_file $data_dir $figures_dir --mode 1
# --
# johaGL 2025
#########################################################################

import os
import pandas as pd
import anndata as ad
import smpath.data as spd
import smpath.method as spm
import smpath.processing.measures_methods as measures_methods
from aux_scripts.data_download_prep_spots.spots_purge import \
    annotation_types_to_file_names_dict


def counts_per_orig_db(annots_df):
    """
    Counts the total of UNIQUE compound identifiers for each type of
    nomenclature (HMDB, LipidMaps, or KEGG identifiers).
    Returns a list of integers, each integer corresponds to the amount of
    UNIQUE identifiers respecting the order: HMDB, LipidMaps, KEGG.
    """
    out_l = [0, 0, 0]
    position_dict = {0: ["HMDB", "CoreMetabolome"], 1: ['LipidMaps'],
                     2: ["KEGG"]}

    for k in position_dict.keys():
        s_df = annots_df.loc[
            annots_df['database'].isin(position_dict[k]), :].copy()
        out_l[k] = len(s_df['moleculeIds'].unique())

    return out_l


def load_datasets_parameters(selection):
    """
    Load datasets parameters based on the selection (public/private).
    """
    file_path = f"../../datasets_parameters/{selection}-datasets-parameters.tsv"
    return pd.read_csv(file_path, header=0, index_col=None, sep='\t')


def initialize_metabolite_mappings(data_dir):
    """
    Initialize and load metabolite mappings.
    """
    me_refmap = spd.MetaboliteIdsMappings(
        data_dir=data_dir,
        subfolder="metabolites-mappings",
        mapping_files={
            'HMDB': "HMDB-to-KEGG_05-07-24.tsv",
            'LipidMaps': "LipidMaps-to-KEGG_10-18-24.tsv",
            'CoreMetabolome': "CoreMetabolome-to-KEGG_05-07-24.tsv"
        }
    )
    me_refmap.load()
    return me_refmap


def initialize_KEGG_pathways(data_dir):
    """
    Initialize and load metabolic pathways.
    """
    pathways_obj = spd.MetabolicPathways(
        data_dir=data_dir,
        subfolder="metabolites-mappings",
        file="KEGG_hsa_pathways_compounds_R110.gmt"
    )
    pathways_obj.load()
    return pathways_obj


def process_sample(row, me_refmap, pathways_obj, data_dir):
    """
    Process each sample and return the counts of identifiers.
    """
    SAMPLE = row["dataset_name"]
    min_pathway_size = row["min_pathway_size"]
    annotations_to_use = row["annotations_to_use"].split(" ")
    annotations_dict = annotation_types_to_file_names_dict(SAMPLE,
                                                           annotations_to_use)

    annots_obj = spd.MetabolomeAnnotations(
        name=SAMPLE, data_dir=data_dir, subfolder=SAMPLE,
        annotation_files=annotations_dict
    )
    annots_obj.load()

    adata_mz = ad.read_h5ad(
        filename=os.path.join(data_dir, SAMPLE, f'{SAMPLE}.h5ad'))

    adata_mz = adata_mz[:, adata_mz.var.loc[
        adata_mz.var['final_quality'] == "good", "mz"].astype(str).tolist()]

    annots_obj, adata_mz = measures_methods.match_adata_mz_and_annots(
        annots_obj, adata_mz)
    counts = counts_per_orig_db(annots_obj.no_translated)

    annots_obj.usable_id = "KEGG"
    annots_obj.convert_ids_to_kegg(me_refmap.table)
    annots_obj.keep_best_mz_id_pairs(adata_mz.var.copy())
    counts += counts_per_orig_db(annots_obj.df)

    pws1, ann1_df = spm.SyncPathwaysAndAnnots.run(
        pathways_obj.model_copy(),
        annots_obj.df.copy(),
        usable_id="KEGG",
        min_set_size=1,
        max_set_size=200,
        mode="mz_aware"
    )
    counts += counts_per_orig_db(ann1_df)

    pws_obj_2, annots_2_df = spm.SyncPathwaysAndAnnots.run(
        pathways_obj.model_copy(),
        annots_obj.df.copy(),
        usable_id="KEGG",
        min_set_size=min_pathway_size,
        max_set_size=200,
        mode="mz_aware"
    )

    return SAMPLE, counts, annots_2_df


def save_results(my_new_dict, name_sample, mz_mapped_to_pathways,
                 kegg_compound_ids_to_mz_ratio_after_pathway_mapped,
                 how_many_pathways, selection):
    """
    Save the results to CSV files.
    """
    out_df1 = pd.DataFrame.from_dict(my_new_dict, orient='index')
    out_df1.reset_index(inplace=True)
    out_df1.columns = ["sample_name", "HMDB identifiers all",
                       "LipidMaps identifiers all", "KEGG identifiers all",
                       "HMDB ids having KEGG ids",
                       "LipidMaps ids having KEGG ids", "KEGG ids KEGG ids",
                       "HMDB in pathways", "LipidMaps in pathways",
                       "KEGG in pathways"]
    out_df1.to_csv(f"../../figs/mapp-ids-pathways-KEGG-{selection}.csv")

    out_df2 = pd.DataFrame({
        "sample_name": name_sample,
        "mz_mapped_to_pathways": mz_mapped_to_pathways,
        "ratio_keggid_mz": kegg_compound_ids_to_mz_ratio_after_pathway_mapped,
        "pathways_number_final": how_many_pathways
    })
    out_df2.to_csv(f"../../figs/crit-KEGG-compounds-counts-{selection}.csv")


def main(selection: str):
    data_centralized_folder = "../../../data"

    assert os.path.exists("../../figs"), "Error, ../../figs must exist"

    datasets_parameters_df = load_datasets_parameters(selection)
    me_refmap = initialize_metabolite_mappings(data_centralized_folder)
    pathways_obj = initialize_KEGG_pathways(data_centralized_folder)

    my_new_dict = {}
    name_sample = []
    mz_mapped_to_pathways = []
    kegg_compound_ids_to_mz_ratio_after_pathway_mapped = []
    how_many_pathways = []

    for i, row in datasets_parameters_df.iterrows():
        SAMPLE, counts, annots_2_df = process_sample(row, me_refmap,
                                                            pathways_obj,
                                                            data_centralized_folder)
        my_new_dict[SAMPLE] = counts
        name_sample.append(SAMPLE)
        mz_ok = annots_2_df['mz'].unique()
        kg_ok = annots_2_df["KEGG"].unique()
        mz_mapped_to_pathways.append(len(mz_ok))
        kegg_compound_ids_to_mz_ratio_after_pathway_mapped.append(
            len(set(kg_ok)) / len(set(mz_ok)) if len(set(mz_ok)) > 0 else 0)
        how_many_pathways.append(len(
            spm.SyncPathwaysAndAnnots.run(pathways_obj.model_copy(),
                                          annots_2_df.copy(),
                                          usable_id="KEGG", min_set_size=row[
                    "min_pathway_size"], max_set_size=200, mode="mz_aware")[
                0].piled['pathway_id'].unique()))

    save_results(my_new_dict, name_sample, mz_mapped_to_pathways,
                 kegg_compound_ids_to_mz_ratio_after_pathway_mapped,
                 how_many_pathways, selection)


if __name__ == "__main__":
    main(selection="public")  # Change to "private" as needed


    # my_csv_table_path_1 = "../../datasets_parameters/public-datasets-parameters.tsv"  #
    # my_csv_table_path_2 = "../../datasets_parameters/private-datasets-parameters.tsv"  #

    #assert os.path.exists("../../figs"), "Error, ../../figs must exist"

    # df1 = pd.read_csv(my_csv_table_path_1, header=0,
    #                   index_col=None, sep='\t')
    # df2 = pd.read_csv(my_csv_table_path_2, header=0,
    #                   index_col=None, sep='\t')
    # datasets_parameters_df = pd.concat([df1, df2])
    # datasets_parameters_df.reset_index(inplace=True)

