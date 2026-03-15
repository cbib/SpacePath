import os
import pandas as pd
import numpy as np
import holoviews as hv
from holoviews import opts, dim
import matplotlib.pyplot as plt
import seaborn as sns
import squidpy as sq
import matplotlib as mpl

###### these  lines are mandatory to comply with Nature journals ######
mpl.rcParams['font.family'] = "sans-serif"
mpl.rcParams['font.sans-serif'] = "Arial"
plt.rcParams.update({'font.size': 18})

hv.extension('matplotlib')   # for saving as .pdf, the alluvial plot


def prepare_mzs_for_plotting(mapping_info, diff_sets, diff_mz_relevant ):
    prep_mzs = mapping_info[['pathway_id', 'mz']].copy()
    prep_mzs = prep_mzs.loc[prep_mzs['mz'].astype(float).isin(
        diff_mz_relevant['mz'].astype(float).tolist()
    ), :]
    prep_mzs = prep_mzs.drop_duplicates().reset_index(drop=True)
    prep_mzs = pd.merge(prep_mzs, diff_sets, how='left', on='pathway_id',
                        left_index=False, right_index=False)
    prep_mzs = prep_mzs[['pathway_id', 'pathway_name', 'mz']].copy()

    return prep_mzs


def run_alluvial_plot(SAMPLE, diff_sets, mapping_info, diff_mz_relevant,
    figures_centralized_folder, subfolder,
                      width_edges_layer1=1, width_edges_layer2=0.1) -> None:
    """
    diff_sets:
        group	names	scores	logfoldchanges	pvals	pvals_adj	abslfc	pathway_id	pathway_name	I
        0	WP4072	36.483738	2.57	2.07e-291	1.92e-29	2.137	WP4072	Glycogen metabolism	0.23
        1	WP4050	37.605724	1.78	1.78e-309	1.19e-307	1.3878	WP4050	Pentose phosphate pathway	0.61
        1	hsa00230	35.557163	-1.32	6.41e-277	1.46e-275	1.322	hsa00230	Purine metabolism	0.56

    mapping_info:
                pathway_id   moleculeIds          mz
           0     WP1403   HMDB0003540  461.165659
           1     WP1403   HMDB0000122  294.165962
           2      WP167  LMFA03010021  430.306419
           3      WP167  LMFA03020023  432.322069
    diff_mz_relevant:
       group  names  scores  logfoldchanges pvals pvals_adj abslfc	  mz
        1	294.165962	 4.59  3.1433194  4.3e-06  6.7e-06	3.14	294.165962
        2	760.548697	 40.9  -8.013631	     0.0	 0.0   8.011   760.548697
        2	834.621861	 39.2      7.   	 0.0    0.0	   7.6245   834.621861
        2	214.191389	 40.2  -2.8529496	0.0   	0.0	 2.85	 214.191389

    """
    # do a 'source' 'target' dataframe with all the types of information
    # for cluster,Pathway edges, depending on logfc: up 1, down 0
    np.random.seed(0)
    truncation_pw_name = 50

    delim_merged_pw_name = "  "

    prep_paths = diff_sets[
        ['group', 'logfoldchanges', 'pathway_id', 'pathway_name']].copy()
    prep_paths['clust_nb'] = prep_paths.group.astype(int).tolist()
    prep_paths = prep_paths.sort_values(by='clust_nb', ascending=False)  # sankey starts from lower-left corner
    prep_paths['source'] = prep_paths.group.astype(str).tolist()
    prep_paths['target'] = prep_paths.pathway_id.str.cat(
        prep_paths.pathway_name.str.slice(0, truncation_pw_name),
                                           sep=delim_merged_pw_name)
    prep_paths['theoretical_direction'] = [ (lambda x: 1.0 if x > 0.0
        else 0.0) (x) for x in prep_paths.logfoldchanges.tolist() ]

    prep_paths = prep_paths[['source', 'target','theoretical_direction', 'pathway_id' ]].copy()
    prep_paths['value'] = width_edges_layer1  # the width of the edges layer 1
    prep_paths = prep_paths.reset_index(drop=True)
    prep_mzs = prepare_mzs_for_plotting(mapping_info, diff_sets, diff_mz_relevant)
    prep_mzs['source'] = prep_mzs.pathway_id.str.cat(
        prep_mzs.pathway_name.str.slice(0, truncation_pw_name),
        sep=delim_merged_pw_name)
    prep_mzs = prep_mzs.dropna()
    prep_mzs['target'] = prep_mzs['mz'].astype(str).tolist()
    prep_mzs = prep_mzs[['source', 'target', 'pathway_id']].copy()
    prep_mzs['theoretical_direction'] = np.nan  # it varies across clusters
    prep_mzs['value'] = width_edges_layer2  # the width of the edges layer 2

    ##
    intersec_pws = list(set(
        prep_paths['pathway_id'].astype(str).tolist()
    ).intersection(set(
        prep_mzs['pathway_id'].astype(str).tolist()
    )))

    prep_paths = prep_paths.loc[prep_paths.pathway_id.astype(str).isin(intersec_pws),:].copy()
    prep_mzs = prep_mzs.loc[
        prep_mzs.pathway_id.astype(str).isin(intersec_pws), :].copy()

    df = pd.concat([prep_paths, prep_mzs], axis=0, ignore_index=False)
    df.reset_index(inplace=True)
    df = df.drop_duplicates()

    df['edge_colors'] = [(lambda x: 'lightsteelblue' if x == 0.0
        else ('sandybrown' if x == 1.0 else 'silver'))(x)
        for x in df['theoretical_direction'].astype(float).tolist()]

    #pastel_palette = sns.color_palette("pastel").as_hex()

    sankey = hv.Sankey(df,  ['source', 'target'], ['value', 'edge_colors'])

    sankey.opts(
        edge_color='edge_colors',
        node_width=20,  # 4x the default (default is usually 10)
        label_position='right',
        # width=800,  Option 'width' for Sankey type not valid for selected backend ('matplotlib')
        # height=600,  Option 'height' for Sankey type not valid for selected backend ('matplotlib')
        node_sort=False,
        node_padding=5,
        show_values=False,
        #cmap=pastel_palette  # this dis not work:   TypeError: unhashable type: '_ColorPalette'
    )
    plt_fig = hv.render(sankey, backend='matplotlib')
    path_sankey_plot = os.path.join(figures_centralized_folder,
                                    subfolder, f"{SAMPLE}-result-holistic.pdf")
    plt_fig.savefig(path_sankey_plot)


def run_individual_mz_scatters(SAMPLE, diff_sets, mapping_info, diff_mz_relevant,
                               adata_mz,
                      figures_centralized_folder, subfolder)->None:
    mzs_df = prepare_mzs_for_plotting(mapping_info, diff_sets, diff_mz_relevant)
    mzs_df = mzs_df.drop_duplicates()
    plot_folder = os.path.join(figures_centralized_folder,
                         subfolder, "prominent")
    if mzs_df.shape[0] != 0:
        for i_mz in mzs_df['mz'].astype(str).unique().tolist():
            inner_plot_path = os.path.join(
                plot_folder, f'{SAMPLE}-top-signif-mz{i_mz}.pdf')

            fig, ax = plt.subplots(1,1)
            sq.pl.spatial_scatter(
                adata_mz, color=[str(i_mz)], img=False,
                title=f'$m/z$: {i_mz}',
                linewidth=0, edgecolors=(0, 0, 0, 0),
                shape="circle", size=1.2,
                palette='viridis',
                ax=ax
            )
            ax.collections[0].set_rasterized(True)
            figure_colorbar = fig.axes[-1]  # take cbar, for scient. notation
            figure_colorbar.ticklabel_format(style='sci', scilimits=(0, 0))
            plt.savefig(inner_plot_path,
                        format="pdf", dpi=400,
                        bbox_inches='tight')
            plt.close()
    else:
        print("the selection of top ions yielded zero results")


def truncate_df_with_text_last_cell(df, n):
    df = df.reset_index(drop=True)
    df = df.drop_duplicates(ignore_index=True)
    if df.shape[0] > n:
        n_more = df.shape[0] - n
        df = df.head(n)  # limit the lenght of the total annots of the current mz
        df.loc[n] = ["...", f"...(truncated {n_more} more)...", "..."]

    return df


def do_for_i_mz(i_mz_df, i_mz, out_plot ):
    auto_height = i_mz_df.shape[0] * 0.2
    if i_mz_df.shape[0] <= 6:
        auto_height = i_mz_df.shape[0] * 0.4
    if i_mz_df.shape[0] == 1:
        auto_height = i_mz_df.shape[0] * 0.7

    auto_width = 3.5
    fig, ax = plt.subplots(nrows=1, ncols=1,
                              figsize=(auto_width+1, auto_height+1))
    ax.set_axis_off()
    table = ax.table(  # https://www.statology.org/seaborn-table/
        cellText=i_mz_df.values, rowLabels=None,
        colLabels=i_mz_df.columns,
        bbox=(-0.08, 0.3, auto_width, auto_height),   # [xmin, ymin, w, h]
        colWidths=(0.2, 0.43, 0.37),
        cellLoc='center')

    table.auto_set_font_size(False)
    table.set_fontsize(20)
    # Left-align the third column
    for (row, col), cell in table.get_celld().items():
        if col == 2:  # Third column (0-based index)
            cell.set_text_props(ha='left')
    ##
    fig.text(0.1, 0.03, f'$m/z$:{i_mz}', ha='left', fontsize=23)
    try:
        plt.tight_layout()
    except:
        pass
    plt.savefig(out_plot, pad_inches=0.5, bbox_inches='tight')
    plt.close()


def run_individual_tabular_plots(
        SAMPLE, diff_sets, mapping_info, diff_mz_relevant, concatenated_clean_annotations,
        figures_centralized_folder, subfolder
    ):
    mzs_df = prepare_mzs_for_plotting(mapping_info, diff_sets,
                                      diff_mz_relevant)
    mzs_df = mzs_df.drop_duplicates()
    # print(mzs_df.head(1))
    # print(mapping_info.head(1))
    # print(diff_mz_relevant.head(1))
    # print(concatenated_clean_annotations.head(1))
    lili = mzs_df['mz'].astype(float).unique().tolist()
    for i_mz in lili:
        curr_mz_df = concatenated_clean_annotations.loc[
            concatenated_clean_annotations['mz'].astype(float) == float(i_mz), ['moleculeIds', 'moleculeNames']].copy()
        curr_mz_df['moleculeNames'] = curr_mz_df['moleculeNames'].str.slice(0, 35)
        my_list_awsome = list()
        for mole_id in curr_mz_df['moleculeIds'].unique().tolist():
            the_pws = mapping_info.loc[mapping_info['moleculeIds'] == mole_id, 'pathway_id'].unique().tolist()
            the_pws_trunc = the_pws.copy()
            if len(the_pws) > 3:
                the_pws_trunc = the_pws[:3] + ["..."]
            my_list_awsome.append(', '.join(the_pws_trunc))

        curr_mz_df['pathway_id'] = my_list_awsome
        curr_mz_df = curr_mz_df.loc[curr_mz_df["pathway_id"] != "",:].copy()
        curr_mz_df = curr_mz_df.sort_values(by='moleculeIds', ascending=True)
        #print(curr_mz_df.head())
        print(curr_mz_df["pathway_id"].unique())
        if curr_mz_df.shape[0] > 0:
            LIMIT_TAB = 15
            if curr_mz_df.shape[0] > LIMIT_TAB:
                curr_mz_df = truncate_df_with_text_last_cell(curr_mz_df, LIMIT_TAB)

            out_plot= os.path.join(
                figures_centralized_folder, subfolder, "prominent",
                f"{SAMPLE}-tab-plot-{i_mz}.pdf"
            )
            do_for_i_mz(curr_mz_df, i_mz, out_plot)
