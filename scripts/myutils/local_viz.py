import seaborn as sns
import matplotlib.pyplot as plt


def customized_clustermap(arr_df, computed_dendrogram):
    cg = sns.clustermap(arr_df, figsize=(12.5, 11),
                        row_linkage=computed_dendrogram,
                        col_linkage=computed_dendrogram,
                        cmap="Blues",
                        dendrogram_ratio=0.15,
                        cbar_pos=(0.07, 0.84, .015, .14),
                        xticklabels=False)
    plt.text(-1.7, 0.23, "similarity index", rotation=90)

    cg.ax_row_dendrogram.set_visible(False)
    cg.ax_heatmap.set_ylabel("")

    return cg