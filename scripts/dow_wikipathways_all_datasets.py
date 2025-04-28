import sys
import os
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
import seaborn as sns
import squidpy as sq
import zonist.data as zd
import zonist.method as zm

from bank4zonist.quantifs_methods import match_adata_mz_and_annots
from bank4zonist.constants import provide_annotation_dict

from do_kegg_translation__describe import (mz_and_compounds_none_pathway,
                                           give_list_compound_mapped_ok,
                                           give_list_one_single_compound_in_any_pathway,
                                           save_detailed_table)
from bank4zonist.visualization import *
from scripts.doplant_wikipathways import min_pathway_size


def thisfun1(ids_mixture_list):
    hmdb_l = [u for u in ids_mixture_list if u.startswith("HMDB")]
    lpmp_l = [u for u in ids_mixture_list if u.startswith("LM")]
    kegg_l = [u for u in ids_mixture_list if u.startswith("C")]

    return [len(hmdb_l), len(lpmp_l), len(kegg_l)]

"""
"""
if __name__ == "__main__":
    # also add to the anndata.var if a feature was yes Y or no N in pathways that will pass.
    # that will allow to "paste" the mz or kegg that were not in the pw trick (? if makes sense such "paste" ?)
    name_sample = []
    mz_mapped_to_pathways = []
    compound_ids_to_mz_ratio_after_pathway_mapped = []
    how_many_pathways = []

    my_new_dict = dict()

    # define here my files and dirs
    minpwsize = 1
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

    # me_refmap = zd.MetaboliteIdsMappings(
    #         data_dir=mydata_dir,
    #         subfolder="metabolites-mappings",
    #         mapping_files={'HMDB': "HMDB-to-KEGG_05-07-24.tsv",
    #                        'LipidMaps': "LipidMaps-to-KEGG_10-18-24.tsv",
    #                        'CoreMetabolome': "CoreMetabolome-to-KEGG_05-07-24.tsv"})
    # me_refmap.load()


    print("doing with wikipathways")

    pathwaysObj = zd.MetabolicPathways(
        data_dir=mydata_dir,
        subfolder="metabolites-mappings",
        official_file="WIKIPATHWAYS_human_multi_id_redundant.mmt"  # TODO: mark as new !!!
    )
    pathwaysObj.molecule_id_type = "moleculeIds"
    pathwaysObj.load()

    print("\nDo XXXX across specified datasets")

    for i, row in datasets_parameters_df.iterrows():
        SAMPLE = row["dataset_name"]
        #min_pathway_size = row["min_pathway_size"]
        min_pathway_size = minpwsize
        print("\nSample: ", SAMPLE)
        subfolder = SAMPLE
        plot_folder = os.path.join(plot_parent_folder, SAMPLE)
        annotations_to_use = row["annotations_to_use"].split(" ")
        anno_dict = provide_annotation_dict(SAMPLE, annotations_to_use)

        me_annots = zd.MetabolomeAnnotations(
            name=SAMPLE, data_dir=mydata_dir, subfolder=subfolder,
            annotation_files=anno_dict)

        me_annots.usable_id = "moleculeIds"  # todo: mark as new

        me_annots.load()

        adata_mz = ad.read_h5ad(filename=os.path.join(
            mydata_dir, subfolder, f'{SAMPLE}.h5ad'))

        print(adata_mz.uns["history"])
        print("MANDATORY: take away bad quality m/z")  # todo: how to make an alert ?
        adata_mz = adata_mz[:,
                   adata_mz.var.loc[adata_mz.var['final_quality'] == "good",
                   "mz"].astype(str).tolist()]

        me_annots, adata_mz = match_adata_mz_and_annots(me_annots, adata_mz)
        me_annots.verify()
        me_annots.keep_unique_mz_id_tuples(adata_mz.var.copy())  #

        my_new_dict[SAMPLE] = thisfun1(me_annots.df[me_annots.usable_id
                                      ].unique())

        initial_pws_dict = pathwaysObj.to_dict()  # reserver pour plustard ? simplement avoir les tailles des sets du départ.
        # none mapping
        mzno, theIDnono = mz_and_compounds_none_pathway(pathwaysObj,
                                                        me_annots,
                                                        id_col=me_annots.usable_id)
        save_detailed_table(adata_mz, mzno, theIDnono, me_annots, plot_folder,
                            SAMPLE, suffix="IDSmix-NotMapped", id_col=me_annots.usable_id)
        # correct mapping
        mz_ok, theIDok = give_list_compound_mapped_ok(
            pathwaysObj, me_annots,
            min_pathway_size=min_pathway_size,
        id_col=me_annots.usable_id)
        save_detailed_table(adata_mz, mz_ok, theIDok, me_annots, plot_folder,
                            SAMPLE, suffix="IDSmix-Mapped", id_col=me_annots.usable_id)
        # oneTOone mapping
        mz_one, theID_oneone = give_list_one_single_compound_in_any_pathway(
            pathwaysObj, me_annots, id_col=me_annots.usable_id)
        mz_one = list(set(mz_one) - set(mz_ok))
        theID_oneone = list(set(theID_oneone) - set(theIDok))
        save_detailed_table(adata_mz, mz_one, theID_oneone, me_annots, plot_folder,
                            SAMPLE, suffix="IDSmix-SinglesMapped", id_col=me_annots.usable_id)

        pws_obj_2, annots_2 = zm.SyncPathwaysAndAnnots.run(
            pathwaysObj.model_copy(),
            me_annots.model_copy(),
            min_pathway_size=min_pathway_size,
            max_pathway_size=200,
            mode="mz_aware"
        )

        my_new_dict[SAMPLE] += thisfun1(annots_2.df[me_annots.usable_id
                                      ].unique())

        print("one mz : one pw (not passed)")
        print("oneTOone ", len(mz_one))
        print("intersected oneTOone and those that passed:", len(set(
            mz_one
        ).intersection(set(mz_ok))))
        print("passsed")
        print("Mzs!!!!!", len(annots_2.df['mz'].unique()), len(mz_ok))
        print("ids?????", len(annots_2.df['moleculeIds'].unique()), len(theIDok))
        print("pathways %%%%%%", len(pws_obj_2.piled['pathway_id'].unique()))

        myratio = 0
        try:
            if len(pws_obj_2.piled["pathway_id"].unique()) != 0:
                myratio = len(annots_2.df['mz'].unique()) / len(pws_obj_2.piled['pathway_id'].unique())
        except:
            pass
        print("ratio mzs/pathwayNumber : ", myratio)
        # print(f"when mapped to KEGG, with min pathwaysize {min_pathway_size}, "
        #       f"remained {me_annots.df.shape[0]} different KEGG id")

        name_sample.append(SAMPLE)
        mz_mapped_to_pathways.append(len(mz_ok))
        try:
            compound_ids_to_mz_ratio_after_pathway_mapped.append(
                len(set(theIDok)) / len(set(mz_ok))
            )
        except ZeroDivisionError:
            compound_ids_to_mz_ratio_after_pathway_mapped.append(0)
        how_many_pathways.append(len(pws_obj_2.piled['pathway_id'].unique()))


    out_df0  = pd.DataFrame.from_dict(
        my_new_dict, orient='index')

    out_df0.reset_index(inplace=True)
    out_df0.columns = ["sample_name", "HMDB identifiers",
                      "LipidMaps identifiers", "KEGG identifiers",
        "HMDB ids in pathways", "LipidMaps ids in pathways",
        "KEGG ids in pathways"]

    print(out_df0)

    # out_df0 = pd.DataFrame({
    #     "sample_name": name_sample,
    #     "HMDB identifiers": hmdb_ids,
    #     "HMDB ids in pathways": hmdb_mapped,
    #     "LipidMaps identifiers": lipidmaps_ids,
    #     "LipidMaps ids in pathways": lipidmaps_mapped,
    #     "KEGG identifiers": kegg_ids,
    #     "KEGG ids in pathways" : kegg_mapped
    # })
    out_df0.to_csv(f"../initial-ids-compounds-IDSmix-psize{minpwsize}.csv")

    out_df = pd.DataFrame({
        "sample_name": name_sample,
        "mz_mapped_to_pathways":mz_mapped_to_pathways,
    "ratio_id_mz": compound_ids_to_mz_ratio_after_pathway_mapped,
    "pathways_number_final" : how_many_pathways    })
    out_df.to_csv(f"../final-IDSmix-compounds-counts-psize{minpwsize}.csv")


