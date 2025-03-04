import os

from typing import List, Dict, Union
import numpy as np
import matplotlib.pyplot as plt

import scanpy as sc
import squidpy as sq


def plot_standard_computs(adata,
                                  plot_folder=None, prefix='sample')-> None:
    k_cols = [u for u in adata.obs.columns.tolist() if u.startswith('kmeans')]
    k = adata.obs[k_cols[0]].astype(int).max() + 1
    if plot_folder is not None:
        tmpstr = os.path.join(plot_folder,
                              f'{prefix}--Kmeans-{k}.png')
        sq.pl.spatial_scatter(adata, color=[f'kmeans{k}'],
                              img=False, title=f"Kmeans (k={k})",
                              save=os.path.abspath(tmpstr))
    else:
        sq.pl.spatial_scatter(adata, color=[f'kmeans{k}'],
                              img=False, title=f"Kmeans (k={k})"
                              )
        plt.show()


def plot_selected_kegg_compound_ids(adata, me_annots,
                                    metabolites_kegg_ids: List[str], sample, plot_folder)-> None:

    for item in metabolites_kegg_ids:
        item_name = me_annots.df.loc[
            me_annots.df["KEGG"] == item, "moleculeNames"].item()
        mz_here = me_annots.df.loc[me_annots.df["KEGG"] == item, "mz"].item()
        tmpstr = os.path.join(plot_folder,
                              f'{sample}-{item_name}.png')
        sq.pl.spatial_scatter(adata, color=[item],
                              img=False,
                              title=f"({mz_here}) {item}: {item_name}",
                              save=os.path.abspath(tmpstr)
                              )


def plot_cosine_distances_random_spots(spots_cos_sim_df):
    import seaborn as sns
    sns.displot(data=spots_cos_sim_df,
                x="x111_y0",
                bins=60, ec="white",
                kde=True)
    plt.show()
    random_spots = np.random.choice(a=spots_cos_sim_df.columns.tolist(),
                                    size=5, replace=False)
    for rspot in random_spots:
        sns.displot(data=spots_cos_sim_df,
                    x=rspot,
                    bins=60, ec="white",
                    kde=True)
        plt.show()


def plot_selected_pathways(adata_scores, biologically_chosen_pathways,
                           SAMPLEANDDETALS, plot_folder) -> None:
    for m_id in biologically_chosen_pathways.keys():
        try:
            pw_name = biologically_chosen_pathways[m_id]
            tmpstr = os.path.join(plot_folder,
                                  f'{SAMPLEANDDETALS}-{pw_name}.png')  # TODO: improve destination folder
            sq.pl.spatial_scatter(adata_scores, color=[m_id],
                                  img=False, title=f"{m_id} {pw_name}",
                                  save=os.path.abspath(tmpstr))
            plt.close()
        except Exception as e:
            print("[plot_selected_pathways]Could NOT plot: {} {}".format(
                m_id, pw_name))



def custom_scatter_one_feat(df, feature: str, the_ax, spot_size=2):

    my_scatter = the_ax.scatter(df['x'], df['y'],
                                c=df[feature], cmap='viridis',
                                vmin=df[feature].min(),
                                vmax=df[feature].max(),
                                s=spot_size)
    the_ax.invert_yaxis()
    the_ax.set_title(f'{feature}')
    cbar = plt.colorbar(
        my_scatter,
        orientation='vertical',
        fraction=0.07,
        aspect=30,
        pad=0.02, anchor=(0, 0))
    the_ax.xaxis.set_tick_params(labelbottom=False)
    the_ax.yaxis.set_tick_params(labelleft=False)
    the_ax.set_xticks([])
    the_ax.set_yticks([])
    plt.tight_layout()

    return the_ax


def plot_with_table(tmp, i_mz, sub_df, out_file_path, spot_size=5, dpi=600) -> None:
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2)
    fig.set_figheight(3)
    fig.set_figwidth(10)

    ax1 = custom_scatter_one_feat(tmp, i_mz, the_ax=ax1,
                                  spot_size=spot_size)  # this replaces  sq.pl.spatial_scatter(adata_mz, color=[str(i_mz)], img=False,  ax=ax1)

    ax2.set_axis_off()
    auto_height = sub_df.shape[0] * 0.12
    if sub_df.shape[0] == 1:
        auto_height = 0.17
    if auto_height >= 10:
        auto_height = 10
    auto_width = 0.7
    table = ax2.table(  # https://www.statology.org/seaborn-table/
        cellText=sub_df.values, rowLabels=None,
        colLabels=sub_df.columns,
        bbox=(-0.08, 0.5, auto_width, auto_height), # [xmin, ymin, width, height]
        colWidths=(0.2, 0.5), cellLoc='center'
    )
    table.set_fontsize(10)

    plt.savefig(out_file_path, dpi=dpi, transparent=True)


def custom_scatter_two_features(df, feature1, feature2, width, height,
                                spot_size=20, save_path:Union[str, None]=None):
    """
    Creates two scatter plots for feature1 and feature with a shared color scale.
    - EXPERIMENTAL -
    Parameters:
    df (pd.DataFrame): DataFrame containing columns [x, y, feature1, feature2].
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

    my_min = df[[feature1, feature2]].min().min()
    my_max = df[[feature1, feature2]].max().max()

    # Create a color map for the scatter plots
    norm = plt.Normalize(df[[feature1, feature2]].min().min(),
                         df[[feature1, feature2]].max().max())

    #sm = plt.cm.ScalarMappable(cmap='viridis',
    #                           norm=norm)
    #sm.set_array([])

    scatter1 = axs[0].scatter(df['x'], df['y'], c=df[feature1], cmap='viridis',
                              #norm=norm  # instead norm, vmin and vmax
                              vmin=my_min, vmax=my_max, s=spot_size)
    axs[0].invert_yaxis()
    axs[0].set_title(f'{feature1}')

    scatter2 = axs[1].scatter(df['x'], df['y'], c=df[feature2], cmap='viridis',
                              #norm=norm  # instead norm, vmin and vmax
                              vmin=my_min, vmax=my_max, s=spot_size)
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
    if save_path is not None:
        plt.savefig(save_path)
    else:
        plt.show()
    plt.close()
