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
import zonist.helpers as zh

from bank4zonist.quantifs_methods import match_adata_mz_and_annots
from bank4zonist.constants import provide_annotation_dict

def thisfun1(ids_mixture_list):
    hmdb_l = [u for u in ids_mixture_list if u.startswith("HMDB")]
    lpmp_l = [u for u in ids_mixture_list if u.startswith("LM")]
    kegg_l = [u for u in ids_mixture_list if u.startswith("C")]

    return [len(hmdb_l), len(lpmp_l), len(kegg_l)]


"""
"""
if __name__ == "__main__":
    name_sample = []
    hmdb_ids = []
    hmdb_mapped = []
    lipidmaps_ids = []
    lipidmaps_mapped= []

    hmdb_inpws = []
    lpmp_inpws = []

    # also add to the anndata.var if a feature was yes Y or no N in pathways that will pass.
    # that will allow to "paste" the mz or kegg that were not in the pw trick (? if makes sense such "paste" ?)

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
        print("\nSample: ", SAMPLE)
        subfolder = SAMPLE
        #plot_folder = os.path.join(plot_parent_folder, SAMPLE)
        annotations_to_use = row["annotations_to_use"].split(" ")
        anno_dict = provide_annotation_dict(SAMPLE, annotations_to_use)

        me_annots = zd.MetabolomeAnnotations(
            name=SAMPLE, data_dir=mydata_dir, subfolder=subfolder,
            annotation_files=anno_dict)
        me_annots.load()

        adata_mz = ad.read_h5ad(filename=os.path.join(
            mydata_dir, subfolder, f'{SAMPLE}.h5ad'))

        print(adata_mz.uns["history"])
        print("MANDATORY: take away bad quality m/z") # todo: how to make an alert ?
        adata_mz = adata_mz[:,
                   adata_mz.var.loc[adata_mz.var['final_quality'] == "good",
                   "mz"].astype(str).tolist()]

        me_annots, adata_mz = match_adata_mz_and_annots(me_annots, adata_mz)
        me_annots.convert_ids_to_kegg(me_refmap.table)
        me_annots.keep_unique_mz_id_tuples(adata_mz.var.copy())
        name_sample.append(SAMPLE)

        ##
        pw1, ann1 = zm.SyncPathwaysAndAnnots.run(
            pathwaysObj.model_copy(),
            me_annots.model_copy(),
            min_pathway_size=minpwsize,
            max_pathway_size=200)
        ##

        #print("!!", len(me_annots.no_translated['mz'].unique()))
        del adata_mz
        for annot_type in ['HMDB', 'LipidMaps']:
            sub_annots_df = me_annots.no_translated.loc[
                            me_annots.no_translated['database'] == annot_type, :]

            sub_annots_piled = zh.spreadout_annotation_df(
                sub_annots_df, mz_col="mz", formula_col="formula",
                molecule_ids_col="moleculeIds", molecule_names_col="moleculeNames")

            sub_refs = me_refmap.table.loc[
                       me_refmap.table['database'] == annot_type, :]

            poo = len(set(
                sub_refs['initial_id_unconverted'].unique()).intersection(
                set(sub_annots_piled['moleculeIds'].unique())))

            inpws = ann1.df.loc[ann1.df['database'] == annot_type, :]
            n_inpws = len(set(inpws['moleculeIds']))

            print("*", annot_type, ":  can be mapped to a KEGG compound ids:", poo, "from a total of", len(sub_annots_piled['moleculeIds'].unique()))
            if annot_type == "HMDB":
                hmdb_ids.append(len(sub_annots_piled['moleculeIds'].unique()))
                hmdb_mapped.append(poo)
                hmdb_inpws.append(n_inpws)
            elif annot_type == "LipidMaps":
                lipidmaps_ids.append(len(sub_annots_piled['moleculeIds'].unique()))
                lipidmaps_mapped.append(poo)
                lpmp_inpws.append(n_inpws)


     # end for
    out_df = pd.DataFrame({
        "sample_name": name_sample,
        "HMDB identifiers all": hmdb_ids,
    "HMDB ids having KEGG ids": hmdb_mapped,
    "LipidMaps identifiers all" : lipidmaps_ids,
    "LipidMaps ids having KEGG ids" : lipidmaps_mapped,
    })
    out_df.to_csv("../initial-ids-compounds-KEGG.csv")

    out_df = pd.DataFrame({
        "sample_name": name_sample,
        "HMDB identifiers all": hmdb_ids,
        "HMDB ids having KEGG ids": hmdb_mapped,
        "HMDB in pathways": hmdb_inpws,
        "LipidMaps identifiers all": lipidmaps_ids,
        "LipidMaps ids having KEGG ids": lipidmaps_mapped,
        "LipidMaps in pathways": lpmp_inpws
    })
    out_df.to_csv("../initial-pathways-KEGG.csv")


        # sub_annots_df = me_annots.no_translated.loc[me_annots.no_translated['database'] == "LipidMaps", :]
        #
        # sub_annots_piled = zh.spreadout_annotation_df(sub_annots_df, mz_col="mz", formula_col="formula",
        #                                               molecule_ids_col= "moleculeIds", molecule_names_col="moleculeNames")
        #
        # sub_refs = me_refmap.table.loc[me_refmap.table['database'] == "LipidMaps", :]
        #
        #
        # print(sub_refs)
        #
        # print("can be mapped to a KEGG compound ids:")
        # print(len(set(sub_refs['initial_id_unconverted'].unique()).intersection(
        #     set(sub_annots_piled['moleculeIds'].unique())
        # )))
        # print("from a total of")
        # print(len(sub_annots_piled['moleculeIds'].unique()))
