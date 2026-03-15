import os
import sys
from pathlib import Path
import anndata as ad
import pandas as pd
import squidpy as sq
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append(str(Path(__file__).resolve().parents[2]))

from SpacePath.src import join_scores_adatas
from SpacePath.src import (join_clean_annotations_from_files
                           )



def print_map(diff_sets_all, data_centralized_folder, subfolder, SAMPLE, pathways_types,
              mapping_info, mymz):
    """check map"""
    """ 
             pathway_id   moleculeIds          mz
        0     WP1403   HMDB0003540  461.165659

        """
    pathways_types = set(pathways_types) - {'FAMILIES'}


    annots = join_clean_annotations_from_files(pathways_types,
                    SAMPLE, subfolder, data_centralized_folder)

    diff_mz = pd.read_csv(os.path.join(data_centralized_folder, subfolder,
                                       f'{SAMPLE}-top-differential-mzs.tsv'),
                          sep='\t', header=0)

    diff_mz = diff_mz.loc[diff_mz['mz'].astype(float) == float(mymz),:]
    print(diff_mz)
    if diff_mz.shape[0] >= 1:
        print("yes, MZ is significant in at least one tissue zone")

    annots = annots.loc[annots['mz'].astype(float) == float(mymz),:]
    print(annots)

    mapping_info = mapping_info.loc[mapping_info['mz'].astype(float) == float(mymz),:]

    sets_names = diff_sets_all.loc[diff_sets_all['pathway_id'].isin(
        mapping_info['pathway_id'].tolist()), "pathway_name"].tolist()

    if len(sets_names) != 0:
        print(sets_names)
    else:
        print(mapping_info['pathway_id'].tolist())

    return 0

def join_clean_annotations_from_files(pathways_types: List[str],
                    SAMPLE: str, subfolder: str,
                   data_centralized_folder: str)-> pd.DataFrame:
    """
    concatenates annotations for each element in pathways_types
     (but not FAMILIES, as the annotations in this case are the same that served to 'WIKI' type of pathway)
     in other words, concatenates KEGG and no-KEGG (that we coined 'WIKI') annotations
    """
    assert "FAMILIES" not in pathways_types, "Error, FAMILIES not allowed"
    annots_dict = read_clean_annotations(SAMPLE, pathways_types,
                                         subfolder, data_centralized_folder)
    dfs_to_gather = list()
    for k in annots_dict.keys():
        if k == 'KEGG':
            df = annots_dict[k][["KEGG", "moleculeNames", "mz"]]
            df.columns = ["moleculeIds", "moleculeNames", "mz"]
        elif k == 'WIKI':
            df = annots_dict[k][["moleculeIds", "moleculeNames", "mz"]]
        dfs_to_gather.append(df)

    gathered_annotations = pd.concat(dfs_to_gather, axis=0)
    return gathered_annotations


if __name__ == "__main__":

    mysample = "Bladder-cancer-hum"     # "Brain-healthy-mouse"
    mymz = "690.506832" #"369.351578"

    # define here my files and dirs
    data_centralized_folder = "../../../data"

    my_csv_table_path_1 = "../../datasets_parameters/public-datasets-parameters.tsv"  #
    my_csv_table_path_2 = "../../datasets_parameters/private-datasets-parameters.tsv"  #

    df1 = pd.read_csv(my_csv_table_path_1, header=0,
                      index_col=None, sep='\t')
    df2 = pd.read_csv(my_csv_table_path_2, header=0,
                      index_col=None, sep='\t')
    datasets_parameters_df = pd.concat([df1, df2])

    datasets_parameters_df = datasets_parameters_df.loc[
        datasets_parameters_df["dataset_name"] == mysample,
    :].copy()

    for i, row in datasets_parameters_df.iterrows():
        SAMPLE = row["dataset_name"]
        min_pathway_size = row["min_pathway_size"]
        #print("\nSample: ", SAMPLE)
        subfolder = SAMPLE
        pathways_types = ["WIKI", "KEGG"]
        #try:
        try:
            adata_mz = ad.read_h5ad(filename=os.path.join(
                data_centralized_folder, subfolder, f'{SAMPLE}.h5ad'))
            # adata_mz = adata_mz[:,
            # adata_mz.var.loc[adata_mz.var['final_quality'] == "good",
            # "mz"].astype(str).tolist()]
            print("attention, bad and good quality mzs are present in the object !")

            print(adata_mz.uns['norm'])
            sns.set_context("paper", font_scale=1.5)

            if adata_mz.uns['norm'] == "TotalIonCount":
                sq.pl.spatial_scatter(
                    adata_mz, color=[str(mymz)], img=False,
                    title=f'$m/z$: {mymz}',
                    linewidth=0, edgecolors=(0, 0, 0, 0),
                    shape="circle", size=1.2,
                    palette='viridis')
                plt.show()
                print("ok")
            else:
                print("Error: only normalized data accepted ! ")
        except FileNotFoundError:
            print("no file (h5ad)")

        try:
            diff_sets_all = pd.read_csv(
                os.path.join(data_centralized_folder, subfolder,
                             f"{SAMPLE}-top-differential-SETS.tsv"),
                sep='\t', header=0
            )
            print("Moran's I : ", adata_mz.var.loc[adata_mz.var['mz'].astype(float) == float(mymz), 'I'].item)
            scores_concatenated = join_scores_adatas(pathways_types, SAMPLE,
                                                     subfolder,
                                                     data_centralized_folder)

            mapping_info = scores_concatenated.uns['current_pathways_df']

            print_map(diff_sets_all, data_centralized_folder, subfolder, SAMPLE,
                      pathways_types,
                      mapping_info, mymz)
        except FileNotFoundError:
            print(SAMPLE,
                 ": could not plot top differential SETS, insufficient data quality or annotations")




