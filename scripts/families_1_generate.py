import os
import pandas as pd
import numpy as np
import zonist.method as zm
import matplotlib.pyplot as plt
import seaborn as sns

import pickle


def differential_results_filter(differential_results, lfc=0.8,
                                pvals_adj=0.05, ntop_per_cluster=70):
    def get_top_n(group, *args):
        n = args[0]
        # todo: fix this warning:  /home/johanna/spatial_thesis/spatial_metabolomics/user-pipeline/scripts/families_1_generate.py:28: DeprecationWarning: DataFrameGroupBy.apply operated on the grouping columns. This behavior is deprecated, and in a future version of pandas the grouping columns will be excluded from the operation. Either pass `include_groups=False` to exclude the groupings or explicitly select the grouping columns after groupby to silence this warning.
        sorted_group = group.sort_values(
            by=['pvals_adj', 'abslfc'], ascending=[True, False] )

        return sorted_group.head(n)

    tmp = differential_results.loc[(abs(differential_results[
             'logfoldchanges']) >= lfc) & (
            differential_results[
                'pvals_adj'] <= pvals_adj), :].copy()

    tmp = tmp.assign(abslfc=abs(differential_results['logfoldchanges']))

    tmp = tmp.groupby('group', group_keys=False, as_index=False).apply(
        get_top_n, ntop_per_cluster).reset_index(drop=True)

    return tmp


if __name__ == "__main__":

    mydata_dir = "../../data"

    my_csv_table_path_1 = "../public-datasets-parameters.tsv"  # "../../small_test.csv"
    my_csv_table_path_2 = "../private-datasets-parameters.tsv"  # "../../small_test.csv"

    df1 = pd.read_csv(my_csv_table_path_1, header=0,
                      index_col=None, sep='\t')
    df2 = pd.read_csv(my_csv_table_path_2, header=0,
                      index_col=None, sep='\t')
    datasets_parameters_df = pd.concat([df1, df2])

    plot_parent_folder = "../../figures/"
    for i, row in datasets_parameters_df.iterrows():
        SAMPLE = row["dataset_name"]
        cutoff_tree = float(row["cut_tree_mol_families"])
        print("\nSample: ", SAMPLE)
        subfolder = SAMPLE
        plot_folder = os.path.join(plot_parent_folder, SAMPLE)

        differential_results = pd.read_csv(os.path.join(
                mydata_dir, subfolder, f"{SAMPLE}-differential-results.csv"),
                                           index_col=None)

        with open(os.path.join(mydata_dir,
                             subfolder,f"{SAMPLE}-WIKI-NO.txt"), "r") as f:
            WIKI_NO = f.readlines()
        WIKI_NO.pop(0)
        WIKI_NO = [float(t.strip()) for t in WIKI_NO]

        with open(os.path.join(mydata_dir,
                             subfolder,f"{SAMPLE}-KEGG-NO.txt"), "r") as f:
            KEGG_NO = f.readlines()
        KEGG_NO.pop(0)
        KEGG_NO = [float(t.strip()) for t in KEGG_NO]

        ALL_MZ_NOT_MAPPED = set(WIKI_NO).union(set(KEGG_NO))

        print("selecting the molecules based on the differential expression of their associated m/z")

        top_differential = differential_results_filter(differential_results)

        del differential_results

        ##
        o_fi1 = os.path.join(mydata_dir,
                             subfolder, f'{SAMPLE}-KEGG-CleanAnnots.tsv')
        o_fi2 = os.path.join(mydata_dir,
                             subfolder, f'{SAMPLE}-WIKI-CleanAnnots.tsv')

        annots_df1 = pd.read_csv(o_fi1,
                                 sep='\t', header=0, index_col=None)
        annots_df2 = pd.read_csv(o_fi2,
                                 sep='\t', header=0, index_col=None)

        print(annots_df1.shape, "KEGG")
        print(annots_df2.shape, "WIKI")
        print("I take so called 'wiki' annotations, because richer in lipids")

        print(
            "mz for further analyses: those not mapping to pathways AND differentially abundant between clusters")
        molecules_ids_df = pd.merge(top_differential, annots_df2, how = 'left',
        left_on = 'names', right_on = 'mz')
        molecules_ids_df.dropna(how='any', inplace=True)
        molecules_ids_df = molecules_ids_df.loc[
                           molecules_ids_df['mz'].astype(float).isin(list(ALL_MZ_NOT_MAPPED)),:]
        molecules_ids_df = molecules_ids_df[['moleculeIds', 'mz', 'moleculeNames']].copy()

        molecules_ids_df.drop_duplicates(inplace=True)

        if not os.path.exists(os.path.join(mydata_dir, subfolder,
                f"{SAMPLE}-mol-similarities-matrix.pkl")):
            mol_similarities_result = zm.ComputeMoleculesSimilarityMatrix.run(
                molecules_ids_list=molecules_ids_df["moleculeIds"]
            )

            with open(os.path.join(mydata_dir, subfolder,
                    f"{SAMPLE}-mol-similarities-matrix.pkl"), "wb") as f:
                pickle.dump(mol_similarities_result, f)
        else:
            with open(os.path.join(mydata_dir, subfolder,
                                   f"{SAMPLE}-mol-similarities-matrix.pkl"),
                      "rb") as f:
                mol_similarities_result = pickle.load(f)

        if not os.path.exists(os.path.join(
            mydata_dir, subfolder, f"{SAMPLE}-mol-fam-dict.pkl")):

            tmp = molecules_ids_df.loc[molecules_ids_df['moleculeIds'].isin(
                mol_similarities_result['the_ids']
            ), ["moleculeIds", "moleculeNames"]].copy()
            tmp.index = tmp['moleculeIds']
            tmp = tmp.T[mol_similarities_result['the_ids']]

            the_molecules_names = tmp.T['moleculeNames'].tolist()

            trunc_names = [i[:16] for i in the_molecules_names]

            arr_df = pd.DataFrame(data=mol_similarities_result['array'],
                                  columns=trunc_names,
                                  index=trunc_names)

            # values = np.arange(start=500, stop=5000, step=10)
            # scale_vals = [0.75, 2]
            # scaled_values = min(scale_vals) + (
            #         values - min(values)) * (max(scale_vals) - min(scale_vals)
            #                                  ) / (max(values) - min(values))
            # diffs = abs(values - len(mol_similarities_result['the_ids']))
            # position = diffs.argmin()
            # factor = np.sum(scale_vals)
            # cutoff_tree = np.around(factor - scaled_values[position], 1)

            print("number of compound ids that were used to generate families:", len(mol_similarities_result['the_ids']))
            print("using cutoff tree :", cutoff_tree)
            ##
            computed_dict = zm.ComputeMoleculesFamilies.run(mol_similarities_result,
                                            molecules_ids_df, mydata_dir, subfolder, SAMPLE,
                                            cutoff_tree=cutoff_tree)

            with open(os.path.join(
                mydata_dir, subfolder, f"{SAMPLE}-mol-fam-dict.pkl"), "wb") as f:
                pickle.dump(computed_dict, f)

            computed_dendrogram = computed_dict['computed_dendrogram']
            # try:  # todo: the plot 'clustermap' slows a lot the computations, better move to a different script !
            #     sns.clustermap(arr_df, figsize=(13, 13),
            #                    #row_linkage=computed_dendrogram,
            #                    cmap="Blues",
            #                    linewidths=0)
            #     plt.title(f"the clustering of molecules")
            #     plt.tight_layout()
            #     plt.savefig(
            #         os.path.join(plot_folder, f"{SAMPLE}-hierarchy-molecules.png"),
            #         dpi=600
            #     )
            # except:
            #     print("skipping clustermap (too heavy)")
        else:
            print('nothing to do')
        print("end")

 # mz_for_further_analyses = list(
        #     set(top_differential['names'].astype(
        #         float).tolist()).intersection(
        #         set(ALL_MZ_NOT_MAPPED))
        # )