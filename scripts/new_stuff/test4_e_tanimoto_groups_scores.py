
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

from  test1_pws_scores import build_ids_adata, run_scores_chapter



if __name__ == "__main__":

   # define here my files and dirs

    pwtype = "TANIMOTO"

    mydata_dir = "../../data"
    plot_parent_folder = "../../figures/"

    my_csv_table_path_1 =  "../public-datasets-parameters.tsv" # "../subset2.tsv"  # "../../small_test.csv"

    # df2 = pd.read_csv(my_csv_table_path_1, header=0,
    #                   index_col=None, sep='\t')
    # datasets_parameters_df = pd.read_csv(my_csv_table_path_1, header=0,
    #                   index_col=None, sep='\t')

    # start

    file_pseudopathways = "tanimoto_groups.mmt"

    pathwaysObj = zd.MetabolicPathways(
            data_dir="../../user-pipeline",
            subfolder="scripts",
            file=file_pseudopathways
    )
    pathwaysObj.molecule_id_type = "moleculeIds"
    pathwaysObj.load()

   # for ... # todo : explain justrunning one sampel for test

    SAMPLE = "TD-Lipido-635"
    subfolder = SAMPLE

    min_pathway_size = 5  # changed !!!
    print("\nSample: ", SAMPLE)
    subfolder = SAMPLE
    plot_folder = os.path.join(plot_parent_folder, SAMPLE)
    annotations_to_use = ["LipidMaps", "HMDB", "CoreMetabolome"]
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
    allgoodmz = adata_mz.var['mz'].astype(float).to_numpy()
    me_annots, adata_mz = zmm.match_adata_mz_and_annots(me_annots, adata_mz)

    me_annots.verify()
    me_annots.keep_unique_mz_id_tuples(adata_mz.var.copy())  #
    me_annots.remove_redundant_aliases()  # when using multitype 'moleculeIds' annotations

    all_mz = me_annots.df['mz'].unique()

    initial_pws_dict = pathwaysObj.to_dict()  # reserver pour plustard ? simplement avoir les tailles des sets du départ.

    pws_obj_2, annots_2 = zm.SyncPathwaysAndAnnots.run(
        pathwaysObj.model_copy(),
        me_annots.model_copy(),
        min_pathway_size=min_pathway_size,  # changed !!!
        max_pathway_size=500,   # changed !!!
        mode="mz_aware")

    mz_in_pws = set(annots_2.df['mz'].unique())
    print("YES:", len(mz_in_pws))

    pd.DataFrame({'mz': np.array(list(mz_in_pws), dtype=float)}).to_csv(
        os.path.join(mydata_dir, subfolder, f'{SAMPLE}-{pwtype}-YES.txt'),
        index=False
    )
   # remember that I deleted the NO part because no sense in tanimoto experiment

    ofipw = os.path.join(mydata_dir, subfolder, f'{SAMPLE}-scores-{pwtype}.h5ad')

    if not os.path.exists(ofipw):
        scores_adata = run_scores_chapter(adata_mz.copy(), pws_obj_2, annots_2)

        scores_adata.write_h5ad(ofipw)

        print("printing extra info")
        for a_pw in scores_adata.var['pathway_id'].unique():
            print(a_pw, pws_obj_2.get_name(a_pw))
            cpw = scores_adata.uns['current_pathways_df']
            cpw = cpw.loc[cpw['pathway_id'] == a_pw].drop(
                columns = ['pathway_id']
            ).sort_values(by="mz")
            print(cpw)
    else:
        print(f"{SAMPLE} : the scores file '{ofipw}' already exists. Nothing to do.")


