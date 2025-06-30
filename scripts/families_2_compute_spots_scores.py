import os
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix

import zonist.data as zd
import zonist.method as zm

from compute_pathways_scores import run_pathway_scores_steps
from local_helpers.local_helpers import provide_annotation_dict
from local_helpers.quantifs_methods import compute_metabolites_probs_df, score_as_adata

if __name__ == "__main__":

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

    too_heavy_do_later= ["UB-breast-cancer-1", "UB-breast-cancer-2", "LP-brain-mouse-3"]
    datasets_parameters_df = datasets_parameters_df.loc[~datasets_parameters_df['dataset_name'].isin(too_heavy_do_later),:]

    min_set_size = 3   # already did with 5

    print("\nDo .....")

    for i, row in datasets_parameters_df.iterrows():
        SAMPLE = row["dataset_name"]
        min_pathway_size = row["min_pathway_size"]
        k = row["k_for_kmeans"]
        print("\nSample: ", SAMPLE)
        subfolder = SAMPLE
        plot_folder = os.path.join(plot_parent_folder, SAMPLE)

        familiesObj = zd.MetabolicPathways(
            data_dir=mydata_dir,
            subfolder=subfolder,
            official_file=f"{SAMPLE}-molecules_families.mmt"
        )
        familiesObj.molecule_id_type = "moleculeIds"

        familiesObj.load()
        o_fi = os.path.join(mydata_dir,
                            subfolder, f'{SAMPLE}-WIKI-CleanAnnots.tsv')
        annots_df = pd.read_csv(o_fi,
                                sep='\t', header=0, index_col=None)

        mydata = ad.read_h5ad(filename=os.path.join(
            mydata_dir, subfolder, f'{SAMPLE}.h5ad'))

        print(mydata.uns["history"])
        print(
            "MANDATORY: take away bad quality m/z")  # todo: how to make an alert ?
        mydata = mydata[:,
                 mydata.var.loc[mydata.var['final_quality'] == "good",
                 "mz"].astype(str).tolist()].copy()

        adata_mz = mydata[:, annots_df['mz'].astype(str).unique()].copy()

        pws_objW, annots_Wdf = zm.SyncPathwaysAndAnnots.run(
            familiesObj.model_copy(),
            annots_df,
            usable_id="moleculeIds",
            min_set_size=min_set_size,
            max_set_size=700,
            mode="mz_aware")

        print("new name to outputs: now include the suffix -sz<MINSETSIZE>.h5ad")
        ofipw = os.path.join(mydata_dir, subfolder,
                             f'{SAMPLE}-scores-Families-sz{min_set_size}.h5ad')

        if not os.path.exists(ofipw):
            scores_adata = run_pathway_scores_steps(adata_mz,
                                                    pws_objW, annots_Wdf, "moleculeIds")

            scores_adata.write_h5ad(ofipw)
        else:
            print("nothing to do")
