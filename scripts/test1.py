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
import zonist.processing.measures_methods as zmm
import zonist.helpers as zh

from local_helpers.local_helpers import provide_annotation_dict
from local_helpers.quantifs_methods import compute_metabolites_probs_df, score_as_adata



def build_ids_adata(adata_mz, X_metabolites_probs_df, annots):
    """
    Given the adata_mz and the annots_df, returns a new adata object
    """
    tmp_X = dict()
    for a_mz in adata_mz.var['mz'].to_numpy():
        the_kegg_here = annots.df.loc[
            annots.df['mz'] == a_mz, annots.usable_id].tolist()
        sub_adata = adata_mz[:, str(a_mz)].copy()
        i_mz_arr = np.array(sub_adata.X.todense().T).flatten()
        for a_kg in the_kegg_here:
            curr_prod_arr = \
                X_metabolites_probs_df[a_kg].to_numpy() * i_mz_arr
            tmp_X[a_kg] = np.array(curr_prod_arr, dtype=np.float32)
    # end for
    identifiers_list = list(tmp_X.keys())
    tmp_X = np.array(pd.DataFrame(tmp_X), dtype=np.float32)
    adataIds = ad.AnnData(X=csr_matrix(tmp_X),
                       obs=adata_mz.obs,
                       var=pd.DataFrame({"ID": identifiers_list},
                                        index=identifiers_list),
                       obsm=adata_mz.obsm)

    #adata.layers["X_probabilities"] = X_metabolites_probs_df
    return adataIds



def run_scores_chapter(adata_mz, pws_obj, annots):
    annots, adata_mz = zmm.match_adata_mz_and_annots(annots, adata_mz)
    X_df = compute_metabolites_probs_df(adata_mz, annots)
    adataIds = build_ids_adata(adata_mz, X_df, annots)
    scores_df = zm.ComputePathwaysScores.run(
        pws_obj.to_dict(), adataIds
    )
    scores_adata = score_as_adata(scores_df, adataIds,
                   pws_obj.pathways_ids,
                   pws_obj.piled)

    return scores_adata


""""
"""
if __name__ == "__main__":

   # define here my files and dirs
    pwtype = "KEGG"

    mydata_dir = "../../data"
    plot_parent_folder = "../../figures/"

    my_csv_table_path_1 = "../subset2.tsv"  # "../../small_test.csv"

    df2 = pd.read_csv(my_csv_table_path_1, header=0,
                      index_col=None, sep='\t')
    datasets_parameters_df = pd.read_csv(my_csv_table_path_1, header=0,
                      index_col=None, sep='\t')

    # start
    me_refmap = zd.MetaboliteIdsMappings(
        data_dir=mydata_dir,
        subfolder="metabolites-mappings",
        mapping_files={'HMDB': "HMDB-to-KEGG_05-07-24.tsv",
                       'LipidMaps': "LipidMaps-to-KEGG_10-18-24.tsv",
                       'CoreMetabolome': "CoreMetabolome-to-KEGG_05-07-24.tsv"})
    me_refmap.load()

    if pwtype == "WIKI":
        pathwaysObj = zd.MetabolicPathways(
            data_dir=mydata_dir,
            subfolder="metabolites-mappings",
            official_file="WIKIPATHWAYS_human_multi_id_redundant.mmt"
        )
        pathwaysObj.molecule_id_type = "moleculeIds"

    elif pwtype == "KEGG":
        pathwaysObj = zd.MetabolicPathways(
            data_dir=mydata_dir,
            subfolder="metabolites-mappings",
            official_file="KEGG_hsa_pathways_compounds_R110.gmt"
        )
        pathwaysObj.molecule_id_type = "KEGG"

    pathwaysObj.load()

    print("\nDo .....")

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

        if pwtype == "WIKI":
            me_annots.usable_id = "moleculeIds"
        elif pwtype == "KEGG":
            me_annots.usable_id = "KEGG"

        me_annots.load()
        adata_mz = ad.read_h5ad(filename=os.path.join(
            mydata_dir, subfolder, f'{SAMPLE}.h5ad'))

        print(adata_mz.uns["history"])
        print("MANDATORY: take away bad quality m/z")  # todo: how to make an alert ?
        adata_mz = adata_mz[:,
                   adata_mz.var.loc[adata_mz.var['final_quality'] == "good",
                   "mz"].astype(str).tolist()]
        allgoodmz = adata_mz.var['mz'].astype(float).to_numpy()
        me_annots, adata_mz = zmm.match_adata_mz_and_annots(me_annots, adata_mz)

        if pwtype == "WIKI":
            me_annots.verify()
        elif pwtype == "KEGG":
            me_annots.convert_ids_to_kegg(me_refmap.table)

        me_annots.keep_unique_mz_id_tuples(adata_mz.var.copy())  #
        if pwtype == "WIKI":
            me_annots.remove_redundant_aliases()  # when using multitype 'moleculeIds' annotations

        all_mz = me_annots.df['mz'].unique()

        initial_pws_dict = pathwaysObj.to_dict()  # reserver pour plustard ? simplement avoir les tailles des sets du départ.

        pws_obj_2, annots_2 = zm.SyncPathwaysAndAnnots.run(
            pathwaysObj.model_copy(),
            me_annots.model_copy(),
            min_pathway_size=min_pathway_size,
            max_pathway_size=200,
            mode="mz_aware")

        mz_not_pws = set(allgoodmz) - set(annots_2.df['mz'].unique())
        mz_in_pws = set(annots_2.df['mz'].unique())
        print("YES:", len(mz_in_pws), "    NO:", len(mz_not_pws))

        pd.DataFrame({'mz': np.array(list(mz_in_pws), dtype=float)}).to_csv(
            os.path.join(mydata_dir, subfolder, f'{SAMPLE}-{pwtype}-YES.txt'),
            index=False
        )

        pd.DataFrame({'mz': np.array(list(mz_not_pws), dtype=float)}).to_csv(
            os.path.join(mydata_dir, subfolder, f'{SAMPLE}-{pwtype}-NO.txt'),
            index=False
        )

        #
        # scores_adata = run_scores_chapter(adata_mz.copy(), pws_obj_2, annots_2)
        #
        # scores_adata.write_h5ad(
        #     os.path.join(mydata_dir, subfolder, f'{SAMPLE}-scores-{pwtype}.h5ad'))
        #
        # print("printing extra info")
        # for a_pw in scores_adata.var['pathway_id'].unique():
        #     print(a_pw, pws_obj_2.get_name(a_pw))
        #     cpw = scores_adata.uns['current_pathways_df']
        #     cpw = cpw.loc[cpw['pathway_id'] == a_pw].drop(
        #         columns = ['pathway_id']
        #     ).sort_values(by="mz")
        #     print(cpw)











