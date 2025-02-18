import os

from typing import List, Dict, Union
import numpy as np
import matplotlib.pyplot as plt

import scanpy as sc
import squidpy as sq


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



def custom_scatter_two_features(df, feature1, feature2, width, height,
                                spot_size=20, save_path:Union[str, None]=None):
    """
    Creates two scatter plots for feature1 and feature with a shared color scale.

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
