import sys
import os
import pandas as pd
import anndata as ad

import zonist.data as zd
import zonist.method as zm
import zonist.processing.measures_methods as zmm

from local_helpers.local_helpers import provide_annotation_dict

from blend_keggids_keggpathways_to_tables import counts_per_orig_db



if __name__ == "__main__":

    print("\nRunning a general quantitative description of the WIKIPATHWAYS"
          " mappings across specified datasets. Output: .csv files with counts")

    name_sample = []
    mz_mapped_to_pathways = []
    compound_ids_to_mz_ratio_after_pathway_mapped = []
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
    datasets_parameters_df.reset_index(inplace=True)
    # start
    print("doing with new version of 'aliases collapse' that is done before the pathways mapping")

    pathwaysObj = zd.MetabolicPathways(
        data_dir=mydata_dir,
        subfolder="metabolites-mappings",
        file="WIKIPATHWAYS_human_multi_id_redundant.mmt"
    )
    pathwaysObj.molecule_id_type = "moleculeIds"
    pathwaysObj.load()

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

        me_annots.usable_id = "moleculeIds"

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

        me_annots.verify()
        me_annots.keep_best_mz_id_pairs(adata_mz.var.copy())

        me_annots.remove_redundant_aliases()  # because using multitype 'moleculeIds' annotations

        initial_pws_dict = pathwaysObj.to_dict()  # reserver pour plustard ? simplement avoir les tailles des sets du départ.

        pw1, ann1_df = zm.SyncPathwaysAndAnnots.run(
            pathwaysObj.model_copy(),
            me_annots.df.copy(),
            usable_id="moleculeIds",
            min_set_size=pwsizeONE,
            max_set_size=200 )

        my_new_dict[SAMPLE] += counts_per_orig_db(ann1_df)

        pws_obj_2, annots_2_df = zm.SyncPathwaysAndAnnots.run(
            pathwaysObj.model_copy(),
            me_annots.df.copy(),
            usable_id="moleculeIds",
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
        theIDok = annots_2_df[me_annots.usable_id].unique()
        mz_mapped_to_pathways.append(len(mz_ok))
        try:
            compound_ids_to_mz_ratio_after_pathway_mapped.append(
                len(set(theIDok)) / len(set(mz_ok))
            )
        except ZeroDivisionError:
            compound_ids_to_mz_ratio_after_pathway_mapped.append(0)

        how_many_pathways.append(len(pws_obj_2.piled['pathway_id'].unique()))


    out_df1  = pd.DataFrame.from_dict(
        my_new_dict, orient='index')

    out_df1.reset_index(inplace=True)
    out_df1.columns = ["sample_name", "HMDB identifiers",
                      "LipidMaps identifiers", "KEGG identifiers",
        "HMDB ids in pathways", "LipidMaps ids in pathways",
        "KEGG ids in pathways"]

    out_df1.to_csv(f"../figs/initial-IDSmix-pathways-psize{pwsizeONE}.csv")

    out_df2 = pd.DataFrame({
        "sample_name": name_sample,
        "mz_mapped_to_pathways":mz_mapped_to_pathways,
    "ratio_id_mz": compound_ids_to_mz_ratio_after_pathway_mapped,
    "pathways_number_final" : how_many_pathways    })
    out_df2.to_csv(f"../figs/final-IDSmix-compounds-counts-psize{min_pathway_size}.csv")
