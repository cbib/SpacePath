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




"""
"""
if __name__ == "__main__":
    name_sample = []
    hmdb_ids = []
    hmdb_mapped = []
    lipidmaps_ids = []
    lipidmaps_mapped= []

    # also add to the anndata.var if a feature was yes Y or no N in pathways that will pass.
    # that will allow to "paste" the mz or kegg that were not in the pw trick (? if makes sense such "paste" ?)

    # define here my files and dirs
    mydata_dir = "../../data"
    my_csv_table_path = "../datasets-parameters.tsv" #"../small-test-parameters.tsv"
    plot_parent_folder = "../../figures/"

    # start
    datasets_parameters_df = pd.read_csv(my_csv_table_path, header=0,
                      index_col=None, sep='\t')

    me_refmap = zd.MetaboliteIdsMappings(
            data_dir=mydata_dir,
            subfolder="metabolites-mappings",
            mapping_files={'HMDB': "HMDB-to-KEGG_05-07-24.tsv",
                           'LipidMaps': "LipidMaps-to-KEGG_10-18-24.tsv",
                           'CoreMetabolome': "CoreMetabolome-to-KEGG_05-07-24.tsv"})
    me_refmap.load()

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
        print("MANDATORY: take away bad quality m/z")
        adata_mz = adata_mz[:,
                   adata_mz.var.loc[adata_mz.var['final_quality'] == "good",
                   "mz"].astype(str).tolist()]

        me_annots, adata_mz = match_adata_mz_and_annots(me_annots, adata_mz)

        name_sample.append(SAMPLE)

        #print("!!", len(me_annots.no_translated['mz'].unique()))
        del adata_mz
        for hoho in ['HMDB', 'LipidMaps']:
            #print(hoho)
            sub_annots_df = me_annots.no_translated.loc[
                            me_annots.no_translated['database'] == hoho, :]

            sub_annots_piled = zh.spreadout_annotation_df(
                sub_annots_df, mz_col="mz", formula_col="formula",
                molecule_ids_col="moleculeIds", molecule_names_col="moleculeNames")

            sub_refs = me_refmap.table.loc[
                       me_refmap.table['database'] == hoho, :]

            poo = len(set(
                sub_refs['initial_id_unconverted'].unique()).intersection(
                set(sub_annots_piled['moleculeIds'].unique())
            ))

            print("*", hoho, ":  can be mapped to a KEGG compound ids:", poo, "from a total of", len(sub_annots_piled['moleculeIds'].unique()))
            if hoho == "HMDB":
                hmdb_ids.append(len(sub_annots_piled['moleculeIds'].unique()))
                hmdb_mapped.append(poo)
            elif hoho == "LipidMaps":
                lipidmaps_ids.append(len(sub_annots_piled['moleculeIds'].unique()))
                lipidmaps_mapped.append(poo)

    out_df = pd.DataFrame({
        "sample_name": name_sample,
        "HMDB identifiers all": hmdb_ids,
    "HMDB ids having KEGG ids": hmdb_mapped,
    "LipidMaps identifiers all" : lipidmaps_ids,
    "LipidMaps ids having KEGG ids" : lipidmaps_mapped,
    })
    out_df.to_csv("../initial_ids_compounds.csv")


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
