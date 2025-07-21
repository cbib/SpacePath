
import numpy as np
import pandas as pd

import squidpy as sq
from matplotlib.axes import Axes
from typing import Union

import matplotlib.pyplot as plt
from matplotlib.figure import Figure


np.random.seed(42)  # important to make the computations be reproducible



def filter_sets_per_cluster_for_plot(diff_copy_df, cluster_to_check,
                                     #cutoff_moran_for_plot=0.3,
                                     #cutoff_abslfc_for_plot=1.5,
                                     top_for_plotting=2):

    diff_copy_df = diff_copy_df.loc[
                 diff_copy_df['group'].astype(str) == str(cluster_to_check), :].copy()

    # diff_copy_df = diff_copy_df.loc[(diff_copy_df['I'].astype(
    #     float) >=cutoff_moran_for_plot) & (
    #         diff_copy_df['abslfc'] >= cutoff_abslfc_for_plot), :].copy()
    # # keep deactivated: moran I is high for artifacts, do not use it anymore

    diff_copy_df = diff_copy_df.sort_values(
        by=['abslfc', 'pvals_adj'], ascending=[False, True]).head(top_for_plotting)

    return diff_copy_df



def custom_scatter_one_feat(df, feature: str, title:str,  the_ax, spot_size=2,
                            marker="o", rasterized=False):

    my_scatter = the_ax.scatter(df['x'], df['y'],
                                c=df[feature], cmap='viridis',
                                vmin=df[feature].min(),
                                vmax=df[feature].max(),
                                s=spot_size,
                                marker=marker,
                                linewidth=0, rasterized=rasterized
                                )
    my_title = f'{feature}'
    if title is not None:
        my_title = title
    the_ax.invert_yaxis()
    the_ax.set_title(my_title)
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


def complex_plot(adata_mz, map_df, set_id, set_full_name , i_mz,
                 save_path=None, dpi=400) -> Union[Figure, None]:

    sup_title = f"m/z = {i_mz}"
    caption = f"ion in set: {set_id}, '{set_full_name[:25]}'   \n"

    curr_mz_df = map_df.loc[map_df['mz'] == i_mz, ['moleculeIds', 'moleculeNames' ]]
    curr_mz_df['moleculeNames'] = curr_mz_df['moleculeNames'].str.slice(0,25)
    mz_in_pw_df = curr_mz_df.loc[curr_mz_df['moleculeIds'].isin(
        map_df.loc[(map_df['mz'] == i_mz) & (map_df['pathway_id'] == set_id), 'moleculeIds'].tolist())  ,:]
    if curr_mz_df.shape[0] > 15:
        n_more = curr_mz_df.shape[0] - 15
        curr_mz_df = curr_mz_df.head(15)  # limit the lenght of the total annots of the current mz
        curr_mz_df.loc[15] = ["...", f"...(truncated {n_more} more)..."]
    tmp = adata_mz[:,str(i_mz)].to_df()

    tmp['x'] = np.array(adata_mz.obsm['spatial'])[:,0]
    tmp['y'] = np.array(adata_mz.obsm['spatial'])[:,1]
    spot_size=10
    marker="s"
    fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3)
    fig.set_figheight(5)
    fig.set_figwidth(15)
    fig.suptitle(sup_title, size=14)

    # ax1 = custom_scatter_one_feat(tmp, feature=str(i_mz),
    #                               title=f"m/z = {i_mz}\n", the_ax=ax1,
    #                               spot_size=spot_size,
    #                               marker=marker, rasterized=True)
    sq.pl.spatial_scatter(
        adata_mz, color=[str(i_mz)], img=False,
        title=f"m/z = {i_mz}\n", size=1.4,
        ax=ax1    #, rasterized=True
    )
    ax2.set_axis_off()
    ax2.set_rasterized(True)
    start_y = 0.09

    auto_height = mz_in_pw_df.shape[0] * 0.12
    if mz_in_pw_df.shape[0] == 1:
        auto_height = 0.17
    if auto_height >= 1.2:
        auto_height = 1.2
    auto_width = 0.9

    ax2.text(auto_width, 0.43, "Annotations mapped to set",
             rotation=270, va='center')

    table = ax2.table(  # https://www.statology.org/seaborn-table/
        cellText=mz_in_pw_df.values, rowLabels=None,
        colLabels=mz_in_pw_df.columns,
        bbox=(-0.08, start_y, auto_width, auto_height),
        # [xmin, ymin, width, height]
        colWidths=(0.43, 0.73), cellLoc='center'
    )
    table.auto_set_font_size(False)
    table.set_fontsize(10.5)
    #table.scale(2, 2)  # may help
    ##
    ax3.set_axis_off()
    ax3.text(auto_width, 0.5, "Total annotations of current m/z",
             rotation=270, va='center')
    auto_height = curr_mz_df.shape[0] * 0.12
    if curr_mz_df.shape[0] == 1:
        auto_height = 0.17
    if auto_height >= 1.2:
        auto_height = 1.2

    table = ax3.table(  # https://www.statology.org/seaborn-table/
        cellText=curr_mz_df.values, rowLabels=None,
        colLabels=curr_mz_df.columns,
        bbox=(-0.08, start_y, auto_width, auto_height),
        colWidths=(0.41, 0.75), cellLoc='center'
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9.5)
    ##
    fig.text(0.5, 0.01, caption, ha='center', fontsize=10)
    plt.tight_layout()
    #
    if save_path is not None:
        plt.savefig(save_path, dpi=dpi)
        plt.close()
    else:
        plt.show()






