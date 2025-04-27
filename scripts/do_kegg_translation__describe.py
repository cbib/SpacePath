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

from bank4zonist.visualization import *


def mz_and_compounds_none_pathway(pathwaysObj, me_annots, id_col="KEGG"):
    initial_kegg_all = me_annots.df[id_col].unique()
    initial_mz_all = me_annots.df['mz'].unique()
    pathways_o_2, annts2 = zm.SyncPathwaysAndAnnots.run(
        pathwaysObj.model_copy(),
        me_annots.model_copy(),
        min_pathway_size=1,
        max_pathway_size=200)

    not_mapping_kegg = set(initial_kegg_all) - set(annts2.df[id_col].unique())
    not_mapping_mz = set(initial_mz_all) - set(annts2.df['mz'].unique())

    return [str(u) for u in not_mapping_mz], list(not_mapping_kegg),


def give_list_one_single_compound_in_any_pathway(pathwaysObj, me_annots, id_col="KEGG"):
    p1, m1 = zm.SyncPathwaysAndAnnots.run(
        pathwaysObj.model_copy(),
        me_annots.model_copy(),
        min_pathway_size=1,
        max_pathway_size=1)

    return [str(u) for u in m1.df["mz"].unique()], m1.df[id_col].unique()


def give_list_compound_mapped_ok(pathwaysObj, me_annots, min_pathway_size, id_col="KEGG"):

    #zm.SyncPathwaysAndAnnotsConfig.mode = "mz_aware"  # TODO: use ?

    p2, m2 = zm.SyncPathwaysAndAnnots.run(
        pathwaysObj.model_copy(),
        me_annots.model_copy(),
        min_pathway_size=min_pathway_size,
        max_pathway_size=200,
        mode="mz_aware"
    )

    return [str(u) for u in m2.df["mz"].unique()], m2.df[id_col].unique()


def save_detailed_table(adata_mz, the_mz_list, the_kegg_list, me_annots, plot_folder,
                        SAMPLE, suffix: str, id_col="KEGG"):
    print(f"Saving things {suffix}")
    out_df = adata_mz.var.copy()
    out_df = out_df.loc[out_df['mz'].astype(str).isin(the_mz_list), :]
    out_df = out_df.drop_duplicates(subset='mz')
    tmp = me_annots.preloaded[
        ["mz", "moleculeIds", "moleculeNames", "database"]].copy()
    tmp['mzETdb'] = tmp['mz'].astype(str).str.cat(tmp['database'], sep="&")
    tmp = tmp.drop_duplicates(subset="mz")
    tmp['mz'] = tmp["mz"].astype(str).tolist()
    intersectionmzs = list(
        set(out_df["mz"].astype(str).tolist()).intersection(
            set(tmp["mz"].astype(str).tolist())))
    tmp = tmp.drop(columns='mzETdb')
    tmp = tmp.set_index("mz")
    out_df = out_df.loc[out_df["mz"].astype(str).isin(intersectionmzs)]
    tmp = tmp.loc[intersectionmzs, :]
    for col in tmp.columns:
        out_df[col] = tmp[col].tolist()
    out_df['the_ids'] = ''
    for i, r in out_df.iterrows():
        ki = me_annots.df.loc[me_annots.df['mz'] == r['mz'], id_col].tolist()
        curr_kegg_ids = list(set(the_kegg_list).intersection(set(
            ki )))
        out_df.loc[out_df['mz']==r['mz'], 'the_ids'] = " ".join(curr_kegg_ids)
    print(f"table dimensions: {out_df.shape[0]}")
    out_df.to_csv(os.path.join(
        plot_folder, f"{SAMPLE}-{suffix}.csv"), index=False)


"""
"""
if __name__ == "__main__":
    # also add to the anndata.var if a feature was yes Y or no N in pathways that will pass.
    # that will allow to "paste" the mz or kegg that were not in the pw trick (? if makes sense such "paste" ?)
    name_sample = []
    mz_mapped_to_pathways = []
    kegg_compound_ids_to_mz_ratio_after_pathway_mapped = []
    how_many_pathways = []

    # define here my files and dirs
    mydata_dir = "../../data"
    my_csv_table_path = "../datasets-parameters.tsv" #"../small-test-parameters.tsv"
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
        official_file="KEGG_hsa_pathways_compounds_R110.gmt")

    pathwaysObj.load()

    print("\nDo XXXX across specified datasets")

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
        me_annots.load()

        adata_mz = ad.read_h5ad(filename=os.path.join(
            mydata_dir, subfolder, f'{SAMPLE}.h5ad'))

        print(adata_mz.uns["history"])
        print("MANDATORY: take away bad quality m/z")  # todo: how to make an alert ?
        adata_mz = adata_mz[:,
                   adata_mz.var.loc[adata_mz.var['final_quality'] == "good",
                   "mz"].astype(str).tolist()]

        me_annots, adata_mz = match_adata_mz_and_annots(me_annots, adata_mz)

        me_annots.convert_ids_to_kegg(me_refmap.table)
        me_annots.keep_unique_mz_id_tuples(adata_mz.var.copy())  #

        initial_pws_dict = pathwaysObj.to_dict()  # reserver pour plustard ? simplement avoir les tailles des sets du départ.
        # none mapping
        mzno, kgno = mz_and_compounds_none_pathway(pathwaysObj,
                                                   me_annots)
        save_detailed_table(adata_mz, mzno, kgno, me_annots, plot_folder,
                            SAMPLE, suffix="KEGG-NotMapped")
        # correct mapping
        mz_ok, kg_ok = give_list_compound_mapped_ok(pathwaysObj, me_annots,
                                                    min_pathway_size=min_pathway_size)
        save_detailed_table(adata_mz, mz_ok, kg_ok, me_annots, plot_folder,
                            SAMPLE, suffix="KEGG-Mapped")
        # oneTOone mapping
        mz_one, kg_one = give_list_one_single_compound_in_any_pathway(pathwaysObj, me_annots)
        mz_one = list(set(mz_one) - set(mz_ok))
        kg_one = list(set(kg_one) - set(kg_ok))
        save_detailed_table(adata_mz, mz_one, kg_one, me_annots, plot_folder,
                            SAMPLE, suffix="KEGG-SinglesMapped")

        pws_obj_2, annots_2 = zm.SyncPathwaysAndAnnots.run(
            pathwaysObj.model_copy(),
            me_annots.model_copy(),
            min_pathway_size=min_pathway_size,
            max_pathway_size=200,
            mode="mz_aware"
        )

        print("one mz : one pw (not passed)")
        print("oneTOone ", len(mz_one))
        print("intersected oneTOone and those that passed:", len(set(
            mz_one
        ).intersection(set(mz_ok))))
        print("passsed")
        print("Mzs!!!!!", len(annots_2.df['mz'].unique()), len(mz_ok))
        print("Kegg?????", len(annots_2.df['KEGG'].unique()), len(kg_ok))
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
            kegg_compound_ids_to_mz_ratio_after_pathway_mapped.append(
               len(set(kg_ok)) / len(set(mz_ok))
            )
        except ZeroDivisionError:
            kegg_compound_ids_to_mz_ratio_after_pathway_mapped.append(0)
        how_many_pathways.append(len(pws_obj_2.piled['pathway_id'].unique()))


    out_df = pd.DataFrame({
        "sample_name": name_sample,
        "mz_mapped_to_pathways":mz_mapped_to_pathways,
    "ratio_keggid_mz": kegg_compound_ids_to_mz_ratio_after_pathway_mapped,
    "pathways_number_final" : how_many_pathways    })
    out_df.to_csv("../final-KEGG-compounds-counts.csv")


