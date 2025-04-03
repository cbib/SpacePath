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

from bank4zonist.quantifs_methods import update_annotation_tomatch_adatamz
from bank4zonist.constants import provide_annotation_dict

from bank4zonist.visualization import *
from user_routines import prepare_harmonized_objects


def mz_and_kegg_compounds_not_mapping_at_all(pathwaysObj, me_annots):
    initial_kegg_all = me_annots.df['KEGG'].unique()
    initial_mz_all = me_annots.df['mz'].unique()
    pathways_o_2, annts2 = zm.SyncPathwaysAndAnnots.run(
        pathwaysObj.model_copy(),
        me_annots.model_copy(),
        min_pathway_size=1,
        max_pathway_size=200)

    not_mapping_kegg = set(initial_kegg_all) - set(annts2.df['KEGG'].unique())
    not_mapping_mz = set(initial_mz_all) - set(annts2.df['mz'].unique())

    return list(not_mapping_kegg), list([str(z) for z in not_mapping_mz])


def save_detailed_table(adata_mz, the_mz_list, the_kegg_list, me_annots, plot_folder,
                        suffix: str):
    print(f"Saving things {suffix}")
    unmaps = adata_mz.var.copy()
    unmaps = unmaps.loc[unmaps['mz'].astype(str).isin(the_mz_list), :]
    unmaps = unmaps.drop_duplicates(subset='mz')
    tmp = me_annots.preloaded[
        ["mz", "moleculeIds", "moleculeNames", "database"]].copy()
    tmp['mzETdb'] = tmp['mz'].astype(str).str.cat(tmp['database'], sep="&")
    me_annots.pull_mz_quality_df()
    tmp = tmp.drop_duplicates(subset="mz")
    tmp['mz'] = tmp["mz"].astype(str).tolist()
    intersectionmzs = list(
        set(unmaps["mz"].astype(str).tolist()).intersection(
            set(tmp["mz"].astype(str).tolist())))
    tmp = tmp.drop(columns='mzETdb')
    tmp = tmp.set_index("mz")
    unmaps = unmaps.loc[unmaps["mz"].astype(str).isin(intersectionmzs)]
    tmp = tmp.loc[intersectionmzs, :]
    for col in tmp.columns:
        unmaps[col] = tmp[col].tolist()

    unmaps['Kegg_not_mapped'] = ''
    for i, r in unmaps.iterrows():
        ki = me_annots.df.loc[me_annots.df['mz'] == r['mz'], "KEGG"].tolist()
        curr_nofound_kegg_list = list(set(the_kegg_list).intersection(set(
            ki  )))
        unmaps.loc[unmaps['mz']==r['mz'], 'Kegg_not_mapped'] = ", ".join(curr_nofound_kegg_list)
    print(f"table dimensions: {unmaps.shape[0]}")
    unmaps.to_csv(os.path.join(
        plot_folder, f"{SAMPLE}-{suffix}.csv"), index=False)


def give_list_one_single_kegg_in_any_pathway(pathwaysObj, me_annots):
    p2, m2 = zm.SyncPathwaysAndAnnots.run(
        pathwaysObj.model_copy(),
        me_annots.model_copy(),
        min_pathway_size=1,
        max_pathway_size=1)

    return m2.df["KEGG"].unique(), m2.df["mz"].unique()


def give_list_kegg_mapped_ok(pathwaysObj, me_annots, min_pathway_size):
    pathwaysObj, me_annots = zm.SyncPathwaysAndAnnots.run(
        pathwaysObj.model_copy(),
        me_annots.model_copy(),
        min_pathway_size=min_pathway_size,
        max_pathway_size=200)

    return me_annots.df["KEGG"].unique(), me_annots.df["mz"].unique()

"""
"""
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
        mydata_dir, subfolder, f'{SAMPLE}-v2.h5ad'))  # ==== todo (?) what to do ?

    print(adata_mz.uns["history"])
    print("MANDATORY: take away bad quality m/z")
    adata_mz = adata_mz[:,
               adata_mz.var.loc[adata_mz.var['final_quality'] == "good",
               "mz"].astype(str).tolist()]

    me_annots = update_annotation_tomatch_adatamz(me_annots, adata_mz)
    me_annots.convert_ids_to_kegg(me_refmap.table)
    me_annots.treat_mz_multiplets(adata_mz)  # cleans .df and also creates property .mz_redundant

    print("===>", set(me_annots.df['mz'].tolist()).intersection(me_annots.mz_redundant))

    print("Mzs?????", len(me_annots.df['mz'].unique()))

    print("Kegg?????", len(me_annots.df['KEGG'].unique()))

    initial_pws_dict = pathwaysObj.to_dict()  # reserver pour plustard ? simplement avoir les tailles des sets du départ.

    # kgno, mzno = mz_and_kegg_compounds_not_mapping_at_all(pathwaysObj,
    #                                                       me_annots)
    # save_detailed_table(adata_mz, mzno, kgno, me_annots, plot_folder,
    #                     suffix="NotMapped")
    #
    # kg_one, mz_one = give_list_one_single_kegg_in_any_pathway(pathwaysObj, me_annots)
    # save_detailed_table(adata_mz, kg_one, mz_one, me_annots, plot_folder,
    #                     suffix="SinglesMapped")

    kg_ok, mz_ok = give_list_kegg_mapped_ok(pathwaysObj, me_annots, min_pathway_size)
    save_detailed_table(adata_mz, kg_ok, mz_ok, me_annots, plot_folder,
                        suffix="Mapped")

    pws_obj_2, annots_2 = zm.SyncPathwaysAndAnnots.run(
        pathwaysObj.model_copy(),
        me_annots.model_copy(),
        min_pathway_size=min_pathway_size,
        max_pathway_size=200)

    #adata_mz = adata_mz[:, ll].copy()  # todo this is NOT sync, but like so
    myratio = 0
    try:
        if len(pws_obj_2.piled["pathway_id"].unique()) != 0:
            myratio = len(annots_2.df['mz'].unique()) / len(pws_obj_2.piled['pathway_id'].unique())
    except:
        pass
    print("ratio mzs/pathwayNumber : ", myratio)
    # print(f"when mapped to KEGG, with min pathwaysize {min_pathway_size}, "
    #       f"remained {me_annots.df.shape[0]} different KEGG id")






#   ll = [str(z) for z in sorted(me_annots.df["mz"].astype(float).unique())]
# adata_mz.var.index = adata_mz.var['mz'].astype(str).tolist()
# adata_mz.var['mz'] = adata_mz.var['mz'].astype(float).to_numpy()
# adata_mz.write_h5ad(os.path.join(mydata_dir, subfolder, f'{SAMPLE}.h5ad'))


#tmp2 = me_annots.mz_quality.copy()
#tmp2['mzETdb'] = tmp2['mz'].astype(str).str.cat(tmp2['database'], sep="&")
#tmp = tmp.loc[tmp['mzETdb'].isin(tmp2['mzETdb'].unique()), :]