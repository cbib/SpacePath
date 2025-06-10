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
from sklearn.cluster import KMeans

from local_helpers.local_helpers import provide_annotation_dict

# TODO: replace all png by pdf files here **** important!


def feed_dictionary(out_dict, df, obs_col, clusters_col, suffix_a, suffix_b):
    """
    input df:
        obs             kmeans_5
        xO_y0           1
        x6_y2           1
        x0_y12          2
        ...             ...
        x4000_y3333     4
        ...             ...
    output:
       a dictionary:
         {1:  [xO_y0 , x6_y2], ... , 4 : [ x4000_y3333, ...]}
    """
    df['tmp_num'] = df[clusters_col].astype(int).tolist()
    df = df.sort_values(by='tmp_num', ascending=True)
    unique_clusters = df[clusters_col].unique()
    for clus in unique_clusters:
        new_key = f'{clus}_{suffix_a}_{suffix_b}'
        out_dict[new_key] = df.loc[df[clusters_col] == clus, obs_col].unique()

    return out_dict


def no_control_groups(df_option, n, plot_parent_folder, SAMPLE, suffix, ww):
    """finding supposedly best matches, without controlling the type of groups"""
    mx0 = np.zeros(shape=(n, n))
    df_fill = pd.DataFrame(data=np.array(mx0, dtype=int),
                            index=list(foo_dict.keys()),
                            columns=list(foo_dict.keys()))
    for i, row in df_fill.iterrows():
        substuff = df_option.loc[i, :]
        names = substuff.index.tolist()
        quants = substuff.to_numpy()
        j_max = np.argmax(quants)
        df_fill.loc[i, names[j_max]] = 1

    sns.heatmap(df_fill, cmap="Blues")
    plt.title(f"BAD candidates by {ww} (not symmetrical)")
    plt.tight_layout()
    plt.savefig(os.path.join(plot_parent_folder, SAMPLE,
                             f"{SAMPLE}--candidates{suffix}.png"))
    plt.close()


""""
"""
if __name__ == "__main__":

   # define here my files and dirs

    mydata_dir = "../../data"
    plot_parent_folder = "../../figures/"

    my_csv_table_path_1 = "../private-datasets-parameters.tsv"  # "../../small_test.csv"

    datasets_parameters_df = pd.read_csv(my_csv_table_path_1, header=0,
                      index_col=None, sep='\t')

    # start
    dict_of_k = { 'NO': 5,
                  'Scores':5}

    k = 5

    print("\nDo .....")

    for i, row in datasets_parameters_df.iterrows():
        SAMPLE = row["dataset_name"]
        min_pathway_size = row["min_pathway_size"]
        print("\nSample: ", SAMPLE)
        subfolder = SAMPLE
        plot_folder = os.path.join(plot_parent_folder, SAMPLE)
        annotations_to_use = row["annotations_to_use"].split(" ")
        anno_dict = provide_annotation_dict(SAMPLE, annotations_to_use)

        mydata = ad.read_h5ad(filename=os.path.join(
            mydata_dir, subfolder, f'{SAMPLE}.h5ad'))

        print(mydata.uns["history"])

        foo_dict = dict()

        # add initial
        df = pd.read_csv(os.path.join(
            mydata_dir, subfolder,
            f'{SAMPLE}-clusters-initial.txt'))
        foo_dict = feed_dictionary(
            foo_dict,
            df, obs_col='obs', clusters_col=f'kmeans_{k}',
            suffix_a="", suffix_b="initial")
        del df


        for pwtype in ['WIKI', 'KEGG']:
            for SUFFIX in ['NO', 'Scores']:
                the_clus_file = os.path.join(
                    mydata_dir, subfolder,
                    f'{SAMPLE}-{pwtype}-clusters-{SUFFIX}-k-{k}.txt')
                df = pd.read_csv(the_clus_file, index_col=None)
                foo_dict = feed_dictionary(
                    foo_dict,
                    df, obs_col='obs', clusters_col=f'kmeans_{k}',
                    suffix_a=pwtype, suffix_b=SUFFIX)

        print(len(foo_dict))
        print(foo_dict.keys())

        matrix_intersections = np.zeros(shape=(len(foo_dict), len(foo_dict)))
        matrix_foo = np.zeros(shape=(len(foo_dict), len(foo_dict)))
        jaccardmx = np.zeros(shape=(len(foo_dict), len(foo_dict)))

        momo = list(foo_dict.keys())
        for i, i_name in enumerate(momo):
            for j, j_name in enumerate(momo):
                size_intersection =  len(
                    set(foo_dict[i_name]).intersection(set(foo_dict[j_name]))
                )
                matrix_intersections[i, j] = size_intersection
                size_union = len(set(foo_dict[i_name]).union(set(foo_dict[j_name])))
                jaccardmx[i, j] = size_intersection / size_union

        threshold_jaccard = 0.1   # todo: make it an internal parameter

        jaccardmx[jaccardmx < threshold_jaccard] = 0  # lets ignore all distances under a threshold
        np.fill_diagonal(jaccardmx, 0)

        print("jaccard df ok")

        jacc_df = pd.DataFrame(data=jaccardmx,
                            index=list(foo_dict.keys()),
                            columns=list(foo_dict.keys()))
        sns.heatmap(jacc_df, cmap="BuPu")
        plt.title(f"jaccard distances. (cutoff={threshold_jaccard})")
        plt.tight_layout()
        plt.savefig(os.path.join(plot_parent_folder, SAMPLE,
                                 f"{SAMPLE}--jaccard.png"))
        plt.close()

        mx_intersections_modified = matrix_intersections.copy()
        np.fill_diagonal(mx_intersections_modified, 0)
        intersects_df = pd.DataFrame(data=np.array(mx_intersections_modified, dtype=int),
                                     index=list(foo_dict.keys()),
                                     columns=list(foo_dict.keys()))
        sns.heatmap(intersects_df, cmap="BuPu")
        plt.title("intersections")
        plt.tight_layout()
        plt.savefig(os.path.join(plot_parent_folder, SAMPLE,
                                 f"{SAMPLE}--intersections.png"))
        plt.close()

        no_control_groups(intersects_df, len(foo_dict),
                          plot_parent_folder, SAMPLE, suffix="BADV1", ww="intersections")
        candidates1 = np.zeros(shape=(len(foo_dict), len(foo_dict)))
        candids1 = pd.DataFrame(data=np.array(candidates1, dtype=int),
                                index=list(foo_dict.keys()),
                                columns=list(foo_dict.keys()))
        for i, row in candids1.iterrows():
            substuff = intersects_df.loc[i, :]
            if str(i).endswith("Scores"):
                selnames = [n for n in candids1.columns.tolist() if n.endswith("NO")]
                sub_sorted = substuff[selnames].sort_values(ascending=False)
                tada = sub_sorted.head(2)
                for loo in tada.index.tolist():
                    candids1.loc[i, loo] = 1
            elif (str(i).endswith("NO")) or (str(i).endswith("initial")):
                selnames = [n for n in candids1.columns.tolist() if n.endswith("Scores")]
                sub_sorted = substuff[selnames].sort_values(ascending=False)
                tada = sub_sorted.head(2)
                for loo in tada.index.tolist():
                    candids1.loc[i, loo] = 1

        sns.heatmap(candids1, cmap="Blues")
        plt.title("candidates by intersections (not symmetrical)")
        plt.tight_layout()
        plt.savefig(os.path.join(plot_parent_folder, SAMPLE,
                                 f"{SAMPLE}--candidatesV1.png"))
        plt.close()

        # TODO: improve this hell
        no_control_groups(jacc_df, len(foo_dict),
                          plot_parent_folder, SAMPLE, suffix="BADV2", ww="jaccard")
        candidates2 = np.zeros(shape=(len(foo_dict), len(foo_dict)))
        candids2 = pd.DataFrame(data=np.array(candidates2, dtype=int),
                                index=list(foo_dict.keys()),
                                columns=list(foo_dict.keys()))
        for i, row in candids2.iterrows():
            substuff = jacc_df.loc[i, :]
            if str(i).endswith("Scores"):
                selnames = [n for n in candids2.columns.tolist() if
                            n.endswith("NO")]
                sub_sorted = substuff[selnames].sort_values(
                    ascending=False)
                tada = sub_sorted.head(2)
                for loo in tada.index.tolist():
                    candids2.loc[i, loo] = 1
            elif (str(i).endswith("NO")) or (str(i).endswith("initial")):
                selnames = [n for n in candids2.columns.tolist() if
                            n.endswith("Scores")]
                sub_sorted = substuff[selnames].sort_values(
                    ascending=False)
                tada = sub_sorted.head(2)
                for loo in tada.index.tolist():
                    candids2.loc[i, loo] = 1

        sns.heatmap(candids2, cmap="Blues")
        plt.title("candidates by jaccard distances (not symmetrical)")
        plt.tight_layout()
        plt.savefig(os.path.join(plot_parent_folder, SAMPLE,
                                 f"{SAMPLE}--candidatesV2.png"))
        plt.close()


###
        # candidates = np.zeros(shape=(len(foo_dict), len(foo_dict)))
        # for i in range(len(momo)):
        #     if mx_intersections_modified[i, :].sum() != 0:
        #         whichmax = np.argmax(np.array(mx_intersections_modified[i, :]))
        #         candidates[i, whichmax] = 1
        #
        # huhu = pd.DataFrame(data=np.array(candidates, dtype=int),
        #                     index=list(foo_dict.keys()),
        #                     columns=list(foo_dict.keys()))
        #
        # sns.heatmap(huhu)
        # plt.title("candidates")
        # plt.show()




        j = 99
        u = 66


