import sys
import os

import anndata as ad

import zonist.data as zd
import zonist.method as zm

import zonist.processing.measures_methods as zmm
from local_helpers.local_helpers import provide_annotation_dict

from local_helpers.local_viz import *



def counts_per_orig_db(annotsdf):
    # the .df attr of annots object has always the database of origin of the metabolites
    out_l = [0, 0, 0]
    position_dict = {0: ["HMDB", "CoreMetabolome"], 1: ['LipidMaps'],
                     2: ["KEGG"]}
    for k in position_dict.keys():
        sdf = annotsdf.loc[
                        annotsdf['database'].isin(position_dict[k]), :].copy()
        out_l[k] = len(sdf['moleculeIds'].unique())
    return out_l


if __name__ == "__main__":
    name_sample = []
    mz_mapped_to_pathways = []
    kegg_compound_ids_to_mz_ratio_after_pathway_mapped = []
    how_many_pathways = []

    my_new_dict = dict()

    # define here my files and dirs
    mydata_dir = "../../data"
    plot_parent_folder = "../../figures/"

    my_csv_table_path_1 = "../public-datasets-parameters.tsv"  # "../../small_test.csv"
    my_csv_table_path_2 = "../private-datasets-parameters.tsv"  # "../../small_test.csv"

    df1 = pd.read_csv(my_csv_table_path_1, header=0,
                      index_col=None, sep='\t')
    df2 = pd.read_csv(my_csv_table_path_2, header=0,
                      index_col=None, sep='\t')
    datasets_parameters_df = pd.concat([df1, df2])

    # start

    me_refmap = zd.MetaboliteIdsMappings(
            data_dir=mydata_dir,
            subfolder="metabolites-mappings",
            mapping_files={'HMDB': "HMDB-to-KEGG_05-07-24.tsv",
                           'LipidMaps': "LipidMaps-to-KEGG_10-18-24.tsv",
                           'CoreMetabolome': "CoreMetabolome-to-KEGG_05-07-24.tsv"})
    me_refmap.load()

    pathwaysObj = zd.MetabolicPathways(
        data_dir=mydata_dir,
        subfolder="metabolites-mappings",
        file="KEGG_hsa_pathways_compounds_R110.gmt")

    pathwaysObj.load()

    print("\nRunning a general quantitative description of the KEGG mappings "
          "(both compound ids and pathways) across specified datasets."
          " Output: .csv files with counts")

    pwsizeONE = 1  # constant

    for i, row in datasets_parameters_df.iterrows():
        SAMPLE = row["dataset_name"]
        min_pathway_size = row["min_pathway_size"]
        print("\nSample: ", SAMPLE)
        subfolder = SAMPLE
        plot_folder = os.path.join(plot_parent_folder, SAMPLE)
        annotations_to_use = row["annotations_to_use"].split(" ")
        anno_dict = provide_annotation_dict(SAMPLE, annotations_to_use)

        me_annots = zd.MetabolomeAnnotations(
            name=SAMPLE, data_dir=mydata_dir, subfolder=subfolder,
            annotation_files=anno_dict)

        me_annots.usable_id = "KEGG"  # because will do translation to KEGG id

        me_annots.load()

        adata_mz = ad.read_h5ad(filename=os.path.join(
            mydata_dir, subfolder, f'{SAMPLE}.h5ad'))

        print(adata_mz.uns["history"])
        print("MANDATORY: take away bad quality m/z")  # todo: how to make an alert ?
        adata_mz = adata_mz[:,
                   adata_mz.var.loc[adata_mz.var['final_quality'] == "good",
                   "mz"].astype(str).tolist()]

        me_annots, adata_mz = zmm.match_adata_mz_and_annots(me_annots, adata_mz)

        my_new_dict[SAMPLE] = counts_per_orig_db(me_annots.no_translated)

        me_annots.convert_ids_to_kegg(me_refmap.table)
        me_annots.keep_best_mz_id_pairs(adata_mz.var.copy())  #

        my_new_dict[SAMPLE] += counts_per_orig_db(me_annots.df)

        initial_pws_dict = pathwaysObj.to_dict()  # reserver pour plustard ? simplement avoir les tailles des sets du départ.

        pws1, ann1_df = zm.SyncPathwaysAndAnnots.run(
            pathwaysObj.model_copy(),
            me_annots.df.copy(),
            usable_id="KEGG",
            min_set_size=pwsizeONE,
            max_set_size=200,
            mode="mz_aware"
        )
        my_new_dict[SAMPLE] += counts_per_orig_db(ann1_df)

        pws_obj_2, annots_2_df = zm.SyncPathwaysAndAnnots.run(
            pathwaysObj.model_copy(),
            me_annots.df.copy(),
            usable_id="KEGG",
            min_set_size=min_pathway_size,
            max_set_size=200,
            mode="mz_aware"
        )

        myratio = 0
        try:
            if len(pws_obj_2.piled["pathway_id"].unique()) != 0:
                myratio = len(annots_2_df['mz'].unique()) / len(pws_obj_2.piled['pathway_id'].unique())
        except:
            pass
        print("ratio mzs/pathwayNumber : ", myratio)

        name_sample.append(SAMPLE)
        mz_ok = annots_2_df['mz'].unique()
        kg_ok = annots_2_df["KEGG"].unique()
        mz_mapped_to_pathways.append(len(mz_ok))
        try:
            kegg_compound_ids_to_mz_ratio_after_pathway_mapped.append(
               len(set(kg_ok)) / len(set(mz_ok))
            )
        except ZeroDivisionError:
            kegg_compound_ids_to_mz_ratio_after_pathway_mapped.append(0)
        how_many_pathways.append(len(pws_obj_2.piled['pathway_id'].unique()))
    # end for

    out_df1 = pd.DataFrame.from_dict(
        my_new_dict, orient='index')

    out_df1.reset_index(inplace=True)
    out_df1.columns = ["sample_name", "HMDB identifiers all", "LipidMaps identifiers all",
                       "KEGG identifiers all",
                       "HMDB ids having KEGG ids",
                       "LipidMaps ids having KEGG ids",
                       "KEGG ids KEGG ids",
                       "HMDB in pathways",
                       "LipidMaps in pathways",
                       "KEGG in pathways"  ]

    out_df1.to_csv("../figs/initial-ids-pathways-KEGG.csv")
    print("saved", "../figs/initial-ids-pathways-KEGG.csv")

    out_df2 = pd.DataFrame({
        "sample_name": name_sample,
        "mz_mapped_to_pathways":mz_mapped_to_pathways,
    "ratio_keggid_mz": kegg_compound_ids_to_mz_ratio_after_pathway_mapped,
    "pathways_number_final" : how_many_pathways    })
    out_df2.to_csv("../figs/final-KEGG-compounds-counts.csv")
    print("saved", "../figs/final-KEGG-compounds-counts.csv")


