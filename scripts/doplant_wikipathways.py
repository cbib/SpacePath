
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
from do_kegg_translation__describe import mz_and_compounds_none_pathway, save_detailed_table, give_list_one_single_compound_in_any_pathway, give_list_compound_mapped_ok
from bank4zonist.visualization import *

mydata_dir = "../../data"
my_csv_table_path = "../small-test-parameters.tsv"  # " ../small-test-parameters.tsv"
plot_parent_folder = "../../figures/"


# start
datasets_parameters_df = pd.read_csv(my_csv_table_path, header=0,
                                     index_col=None, sep='\t')

pathwaysObj = zd.MetabolicPathways(
        data_dir=mydata_dir,
        subfolder="metabolites-mappings",
        official_file="WIKIPATHWAYS_human_multi_id_redundant.mmt")


pathwaysObj.molecule_id_type = "moleculeIds"
pathwaysObj.load()

print("\nMarks bad quality m/z across specified datasets (.var property)")

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
    me_annots.usable_id = "moleculeIds"   # default is KEGG
    me_annots.load()

    adata_mz = ad.read_h5ad(filename=os.path.join(
        mydata_dir, subfolder, f'{SAMPLE}.h5ad'))

    print(adata_mz.uns["history"])
    print(
        "MANDATORY: take away bad quality m/z")  # todo: how to make an alert ?
    adata_mz = adata_mz[:,
               adata_mz.var.loc[adata_mz.var['final_quality'] == "good",
               "mz"].astype(str).tolist()]

    adata_mz

    me_annots, adata_mz = match_adata_mz_and_annots(me_annots, adata_mz)
    print(me_annots.df.head())

    me_annots.verify()
    me_annots.keep_unique_mz_id_tuples(adata_mz.var.copy())
    print("loading pathways dict")
    initial_pws_dict = pathwaysObj.to_dict()  # reserver pour plustard ? simplement avoir les tailles des sets du départ.

    p1, m1 = zm.SyncPathwaysAndAnnots.run(
        pathwaysObj.model_copy(),
        me_annots.model_copy(),
        min_pathway_size=1,
        max_pathway_size=1)


    print(me_annots.df.shape)

    print(m1.df.shape)
    haha = zh.compute_compounds_n_per_mz(adata_mz.var.copy(),  m1.df, 'mz', 'moleculeIds')  #var_df, ann_df, mz_col, id_col)
    haha = haha[['mz', 'fdr', 'I', 'compound_ids_n__per_mz']]
    haha = haha.loc[haha['compound_ids_n__per_mz'] != 0, :]
    hoho = 0



    # none mapping
    mzno, kgno = mz_and_compounds_none_pathway(pathwaysObj,
                                               me_annots, id_col=me_annots.usable_id)
    save_detailed_table(adata_mz, mzno, kgno, me_annots, plot_folder,
                        SAMPLE, suffix="XX-NotMapped", id_col=me_annots.usable_id)

    # correct mapping
    mz_ok, kg_ok = give_list_compound_mapped_ok(pathwaysObj, me_annots,
                                                min_pathway_size=min_pathway_size, id_col=me_annots.usable_id)
    save_detailed_table(adata_mz, mz_ok, kg_ok, me_annots, plot_folder,
                        SAMPLE, suffix="XX-Mapped", id_col=me_annots.usable_id)

    # oneTOone mapping
    mz_one, kg_one = give_list_one_single_compound_in_any_pathway(pathwaysObj,
                                                                  me_annots, id_col=me_annots.usable_id)
    mz_one = list(set(mz_one) - set(mz_ok))
    kg_one = list(set(kg_one) - set(kg_ok))
    save_detailed_table(adata_mz, mz_one, kg_one, me_annots, plot_folder,
                        SAMPLE, suffix="XX-SinglesMapped", id_col=me_annots.usable_id)

    pws_obj_2, annots_2 = zm.SyncPathwaysAndAnnots.run(
        pathwaysObj.model_copy(),
        me_annots.model_copy(),
        min_pathway_size=min_pathway_size,
        max_pathway_size=200,
        mode="mz_aware"
    )

    for u in pws_obj_2.to_dict().keys():
        print(pws_obj_2.pathways_ids.loc[pws_obj_2.pathways_ids['pathway_id']==u, "pathway_name"])



    print("one mz : one pw (not passed)")
    print("oneTOone ", len(mz_one))
    print("intersected oneTOone and those that passed:", len(set(
        mz_one
    ).intersection(set(mz_ok))))
    print("passsed")
    print("Mzs!!!!!", len(annots_2.df['mz'].unique()), len(mz_ok))
    print("IDS?????", len(annots_2.df[me_annots.usable_id].unique()), len(kg_ok))
    print("pathways %%%%%%", len(pws_obj_2.piled['pathway_id'].unique()))

    myratio = 0
    try:
        if len(pws_obj_2.piled["pathway_id"].unique()) != 0:
            myratio = len(annots_2.df['mz'].unique()) / len(
                pws_obj_2.piled['pathway_id'].unique())
    except:
        pass
    print("ratio mzs/pathwayNumber : ", myratio)
    # print(f"when mapped to KEGG, with min pathwaysize {min_pathway_size}, "
    #       f"remained {me_annots.df.shape[0]} different KEGG id")

    # name_sample.append(SAMPLE)
    # mz_mapped_to_pathways.append(len(mz_ok))
    # try:
    #     kegg_compound_ids_to_mz_ratio_after_pathway_mapped.append(
    #         len(set(kg_ok)) / len(set(mz_ok))
    #     )
    # except ZeroDivisionError:
    #     kegg_compound_ids_to_mz_ratio_after_pathway_mapped.append(0)
    # how_many_pathways.append(len(pws_obj_2.piled['pathway_id'].unique()))
    #
    # out_df = pd.DataFrame({
    #     "sample_name": name_sample,
    #     "mz_mapped_to_pathways": mz_mapped_to_pathways,
    #     "ratio_keggid_mz": kegg_compound_ids_to_mz_ratio_after_pathway_mapped,
    #     "pathways_number_final": how_many_pathways})
    # out_df.to_csv("../final-KEGG-compounds-counts.csv")