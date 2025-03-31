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
    pathwaysObj, me_annots = zm.SyncPathwaysAndAnnots.run(
        pathwaysObj.model_copy(),
        me_annots.model_copy(),
        min_pathway_size=1,
        max_pathway_size=200)

    not_mapping_kegg = set(initial_kegg_all) - set(me_annots.df['KEGG'].unique())
    not_mapping_mz = set(initial_mz_all) - set(me_annots.df['mz'].unique())

    return list(not_mapping_kegg), list(not_mapping_mz)


def save_table_unmapped(adata_mz, notmapped_mz, me_annots):
    print("Saving things not mapping to any pathways")
    unmaps = adata_mz.var.copy()
    unmaps = unmaps.loc[unmaps['mz'].astype(str).isin(notmapped_mz), :]
    unmaps = unmaps.loc[
             ~unmaps["mz"].astype(str).isin(me_annots.mz_redundant), :]
    unmaps = unmaps.drop_duplicates(subset='mz')
    print(unmaps, ":::::")
    tmp = me_annots.preloaded[
        ["mz", "moleculeIds", "moleculeNames", "database"]].copy()
    tmp['mzETdb'] = tmp['mz'].astype(str).str.cat(tmp['database'], sep="&")
    me_annots.pull_mz_quality_df()
    tmp2 = me_annots.mz_quality.copy()
    tmp2['mzETdb'] = tmp2['mz'].astype(str).str.cat(tmp2['database'], sep="&")
    tmp = tmp.loc[tmp['mzETdb'].isin(tmp2['mzETdb'].unique()), :]
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
    print(f"Total of unique mz not mapped at all: {unmaps.shape[0]}")
    unmaps.to_csv(os.path.join(
        mydata_dir, subfolder, f"{SAMPLE}-notMappableAll.csv"),
        index=False)
    del tmp


def give_list_one_single_kegg_in_any_pathway(pathwaysObj, me_annots):
    pathwaysObj, me_annots = zm.SyncPathwaysAndAnnots.run(
        pathwaysObj.model_copy(),
        me_annots.model_copy(),
        min_pathway_size=1,
        max_pathway_size=1)

    return me_annots.df["KEGG"].unique(), me_annots.df["mz"].unique()


"""
"""
# also add to the anndata.var if a feature was yes Y or no N in pathways that will pass.
# that will allow to "paste" the mz or kegg that were not in the pw trick (? if makes sense such "paste" ?)


# define here my files and dirs
mydata_dir = "../../data"
my_csv_table_path = "../small-test-parameters.tsv"
plot_parent_folder = "../../figures/"

initial_min_pathway_size = 1

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

print("\nClear offSample and bad quality m/z across specified datasets")

for i, row in datasets_parameters_df.iterrows():
    SAMPLE = row["dataset_name"]

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
    me_annots = update_annotation_tomatch_adatamz(me_annots, adata_mz)  # todo: mark as new
    me_annots.convert_ids_to_kegg(me_refmap.table)


    me_annots.treat_mz_multiplets(adata_mz)  # cleans .df and also creates property .mz_redundant

    initial_pws_dict = pathwaysObj.to_dict()  # reserver pour plustard ? simplement avoir les tailles des sets du départ.

    kgno, mzno = mz_and_kegg_compounds_not_mapping_at_all(pathwaysObj,
                                                          me_annots)
    print(kgno, mzno)
    save_table_unmapped(adata_mz, mzno, me_annots)
    #adata_mz = adata_mz[:, ll].copy()  # todo this is NOT sync, but like so

    # print(f"when mapped to KEGG, with min pathwaysize {min_pathway_size}, "
    #       f"remained {me_annots.df.shape[0]} different KEGG id")






#   ll = [str(z) for z in sorted(me_annots.df["mz"].astype(float).unique())]
 # adata_mz.var.index = adata_mz.var['mz'].astype(str).tolist()
    # adata_mz.var['mz'] = adata_mz.var['mz'].astype(float).to_numpy()
    # adata_mz.write_h5ad(os.path.join(mydata_dir, subfolder, f'{SAMPLE}.h5ad'))