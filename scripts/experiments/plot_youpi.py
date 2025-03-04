import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def custom_scater_two_featuresV0(df, feature1, feature2, width, height):
    """
    Creates two scatter plots for gene1 and gene2 with a shared color scale.

    Parameters:
    df (pd.DataFrame): DataFrame containing columns [x, y, gene1, gene2].
    """
    # Check if the required columns are in the DataFrame
    required_columns = ['x', 'y', feature1, feature2]
    if not all(column in df.columns for column in required_columns):
        raise ValueError(
            f"DataFrame must contain the following columns: {required_columns}")

    # Set up the figure and axes
    fig, axs = plt.subplots(1, 2,
                            width_ratios=(0.47, 0.51),
                             #             - fract_keybar_width,   0.5 + fract_keybar_width),
                            figsize=(width, height), sharey=False)

    # Create a color map for the scatter plots
    norm = plt.Normalize(df[[feature1, feature2]].min().min(),
                         df[[feature1, feature2]].max().max())
    my_min = df[[feature1, feature2]].min().min()
    my_max = df[[feature1, feature2]].max().max()

    #sm = plt.cm.ScalarMappable(cmap='viridis',
    #                           norm=norm)
    #sm.set_array([])

    scatter1 = axs[0].scatter(df['x'], df['y'], c=df[feature1], cmap='viridis',
                              #norm=norm  # instead norm, vmin and vmax
                              vmin=my_min,
                              vmax=my_max)
    axs[0].invert_yaxis()
    axs[0].set_title(f'{feature1}')

    scatter2 = axs[1].scatter(df['x'], df['y'], c=df[feature2], cmap='viridis',
                              #norm=norm  # instead norm, vmin and vmax
                              vmin=my_min,
                              vmax=my_max)
    axs[1].invert_yaxis()
    axs[1].set_title(f'{feature2}')

    # Add colorbar
    cbar = plt.colorbar(
        # sm: not anymore used because uses normalized between 0 and 1
        scatter2,  # not matter 1 or 2 gives same result
        orientation='vertical',
                        fraction=0.07,
                        aspect=30,
                        pad=0.02, anchor=(0, 0))
    for i_ax in axs:
        i_ax.xaxis.set_tick_params(labelbottom=False)
        i_ax.yaxis.set_tick_params(labelleft=False)
        i_ax.set_xticks([])
        i_ax.set_yticks([])
    # Show the plots
    plt.tight_layout()
    plt.show()


def tsne_for_clustering():
    # tested with small dataset 2, slightly better but still needs spatial refinement

    from sklearn.manifold import TSNE
    tsne = TSNE(n_components=2, perplexity=5, n_iter=1000, random_state=42)
    tsne_arr = tsne.fit_transform(np.asarray(Xnew))

    kmeans = KMeans(n_clusters=3, n_init=10, random_state=0
                    ).fit(
        # scores_adata.X
        tsne_arr  # np.asarray(Xnew)
    )
    scores_adata.obs[f'kmeans'] = kmeans.labels_.astype(str)
    closest, _ = pairwise_distances_argmin_min(kmeans.cluster_centers_,
                                               # scores_adata.X
                                               tsne_arr
                                               )
    tsne_df = pd.DataFrame(data=tsne_arr, columns=["TSNE1", "TSNE2"])
    tsne_df['Cluster'] = kmeans.labels_.astype(str)
    for inum_close in closest:
        tsne_df.iloc[inum_close, -1] = "CENTERS"
    # print(pc_df.head(2))

    sns.lmplot(x="TSNE1", y="TSNE2",
               data=tsne_df,
               fit_reg=False,
               hue='Cluster',  # color by cluster
               legend=True,
               scatter_kws={"s": 30})  # specify the point size
    plt.show()

    sq.pl.spatial_scatter(scores_adata, img=False, color="kmeans")
    plt.show()




# Example usage
if __name__ == "__main__":
    # Sample DataFrame
    data = {
        'x': [1, 2, 3, 4, 5],
        'y': [5, 4, 3, 2, 1],
        'm1': [400, 30, 45, 2, 20],
        'm2' : [100, 30, 1437854, 60, 250]
    }
    df = pd.DataFrame(data)

    # Call the function to create scatter plots
    custom_scater_two_featuresV0(df, "m1", "m2",
                                 10, 3)



    # fake = axs[2].scatter(df['x'], df['y'], c=df[feature1], cmap='viridis',
    #                           vmin=my_min,
    #                           vmax=my_max)
