#########################################################################
# requires running first the two 'blend_..._to_tables.py' scripts
# --
# johaGL 2025
#########################################################################
import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


###### these three lines are mandatory to comply with Nature journals ######
import matplotlib as mpl
mpl.rcParams['font.family'] = "sans-serif"
mpl.rcParams['font.sans-serif'] = "Arial"

# todo: add here the global plot params: Arial etc
def update_names_TD_samples(df):
    names_dict = {'TD-Lipido-635': 'GB-x-lipids-1',
                  'TD-Lipido-636': 'GB-x-lipids-2',
    'TD-Lipido-827': 'GB-x-lipids-3',
    'TD-Lipido-836': 'GB-x-lipids-4'
    }
    new_names_samples = list()
    for u in df['sample_name'].tolist():
        if u in list(names_dict.keys()):
            v = names_dict[u]
        else:
            v = u
        new_names_samples.append(v)
    df['sample_name'] = new_names_samples

    return df


def transform2plot_v2(df, mode="all"):
    cols = [col for col in df.columns.tolist() if col.endswith(mode)]
    col_hmdb = [k for k in cols if k.startswith("HMDB")]
    df1 = pd.melt(df, id_vars="sample_name",  # ,
                  value_vars=col_hmdb[0], #"HMDB identifiers all",
                  value_name="amount",
                  var_name="to_drop")
    df1['database'] = "HMDB"
    df1 = df1.drop(columns=["to_drop"])
    col_lipids = [k for k in cols if k.startswith("Lipid")]
    df2 = pd.melt(df, id_vars="sample_name",
                  value_vars=col_lipids[0], #"LipidMaps identifiers all",
                  value_name="amount",
                  var_name="to_drop")
    df2['database'] = "LipidMaps"
    df2 = df2.drop(columns=["to_drop"])

    col_kegg = [k for k in cols if k.startswith("KEGG")]
    df3 = pd.melt(df, id_vars="sample_name",  # ,
                  value_vars=col_kegg[0],  # "HMDB identifiers all",
                  value_name="amount",
                  var_name="to_drop")
    df3['database'] = "KEGG"

    if not np.isnan(df3['amount'].to_numpy()).all():
        df_pl = pd.concat([df1, df2, df3], axis=0)
    else:
        df_pl = pd.concat([df1, df2], axis=0)

    return df_pl


def bubble_plot_kegg_metabolite_id_conversion_success(
        df_pl, df_pl2, size_range_markers1, size_range_markers2,
        fig_width,  fig_height
) -> None:
    sns.set_context("paper", font_scale=1.5)
    sns.set_style("whitegrid")
    plt.figure(figsize=(fig_width, fig_height))
    ax = sns.scatterplot(df_pl, x="database", y="sample_name", size="amount",
                         lw=0.4, sizes=size_range_markers1,
                         # sizes=(200, 1000),
                         color="silver", alpha=0.5,
                         edgecolors="black")
    ax.grid(color='grey', linestyle='--', alpha=0.7, zorder=1)
    sns.scatterplot(df_pl2, ax=ax,
                    x="database", y="sample_name", size="amount",
                    lw=0.4, sizes=size_range_markers2,  # sizes=(15, 60),
                    color="sandybrown", alpha=1,
                    edgecolors="black").legend(loc="upper left", title="count",
                                               bbox_to_anchor=(1,
                                                               1.03))  # .02))
    plt.margins(x=1)
    sns.despine(fig=None, ax=None,
                top=True, right=True, left=True, bottom=True, offset=None,
                trim=False)
    plt.tight_layout()
    plt.savefig(f"../../figs/initial-ids-compounds-KEGG-{my_selection}.pdf")
    plt.close()


def bubble_plot_for_pathways(df_pl, df_pl3, size_range_markers1,
                             size_range_markers3, out_path,
                             fig_width, fig_height) -> None:
    sns.set_context("paper", font_scale=1.5)
    sns.set_style("whitegrid")
    plt.figure(figsize=(fig_width, fig_height))
    ax = sns.scatterplot(df_pl, x="database", y="sample_name", size="amount",
                         lw=0.4, sizes=size_range_markers1,
                         # sizes=(200, 1000),
                         color="silver", alpha=0.5,
                         edgecolors="black")
    ax.grid(color='grey', linestyle='--', alpha=0.7, zorder=1)
    sns.scatterplot(df_pl3, ax=ax,
                    x="database", y="sample_name", size="amount",
                    lw=0.4, sizes=size_range_markers3,  # sizes=(15, 60),
                    color="darkcyan", alpha=1,
                    edgecolors="black").legend(loc="upper left", title="count",
                                               bbox_to_anchor=(1,
                                                               1.03))  # .02))
    plt.margins(x=1)
    sns.despine(fig=None, ax=None,
                top=True, right=True, left=True, bottom=True, offset=None,
                trim=False)
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()


def good_samples_by_criteria(my_selection: str, df):
    tmp = df.loc[(df['ratio_id_mz'] <= 2.2) & (
            df['pathways_number_final'] >= 10) & (
            df['mz_mapped_to_pathways'] >= 10),:].copy()

    if my_selection == 'private':   # avoid in the public cohort to keep Cerebellum tissue
        tmp = tmp.loc[tmp['ratio_mz_pathway_count'] >= 1,:].copy()

    good_samples = tmp["sample_name"].unique()
    print("good tissues:\n", good_samples)

    return good_samples


## start
##  note: this script is a fast draft, more files loading at lines

###  dot plots (as PairGrid) comparing KEGG and WIKIPATHWAYS

my_selection = "public"   #  private
min_pathway_size = 2
print("dot plots (as PairGrid) comparing KEGG and WIKIPATHWAYS")

if my_selection == "public":
    my_height = 9
    my_width = 5.7
else:
    my_height = 3.3
    my_width = 5.3


df = pd.read_csv(f"../../figs/crit-KEGG-compounds-counts-{my_selection}.csv", index_col=None)
df2 = pd.read_csv(f"../../figs/crit-IDSmix-counts-psize{min_pathway_size}-{my_selection}.csv", index_col=None)

df = df[['sample_name', 'ratio_keggid_mz', 'mz_mapped_to_pathways', 'pathways_number_final']]
df2 = df2[['sample_name', 'ratio_id_mz', 'mz_mapped_to_pathways', 'pathways_number_final']]
df = df.replace(0, np.nan)
df.columns = ['sample_name', 'ratio_id_mz', 'mz_mapped_to_pathways', 'pathways_number_final']

df['PType'] = "KEGG_pathways"
df2['PType'] = "Wikipathways"

df = pd.concat([df, df2])
df = update_names_TD_samples(df)

df["ratio_mz_pathway_count"] = (df['mz_mapped_to_pathways'].to_numpy(
) + 0.001) / (df['pathways_number_final'].to_numpy() + 0.001)

df = df[['sample_name', 'ratio_id_mz', 'mz_mapped_to_pathways',
         'pathways_number_final', 'ratio_mz_pathway_count', 'PType']]

print("excluding from the plots the low quality datasets")

good_samples = good_samples_by_criteria(my_selection, df)

df = df.loc[df["sample_name"].isin(good_samples), :]

print(df)

kk = df['sample_name'] + '&' +df['PType']
df.index = kk

#df = df.loc[order_samples,: ]
df['data_type'] = "Lipidomics"
df.loc[df['sample_name'].str.startswith("MALDI"), "data_type"] = "SmallMolecules"
datasmalls = ["TD-635-1-A", "TD-hemi-635-1-B", "LP-heart-mouse-2"]
df.loc[df['sample_name'].isin(datasmalls), "data_type"] = "SmallMolecules"


sns.set_theme(style="whitegrid",font_scale=1.4)
#plt.figure(figsize=(60, 8))   # very difficult to set the dimensions, instead saving manually from display window

g = sns.PairGrid(df,
                 x_vars=df.columns[1:5], y_vars="sample_name",
                 height=10,
                 aspect=.25)

# Draw a dot plot using the stripplot function
g.map(sns.stripplot, data=df,  size=10, orient="h", jitter=False,
hue="PType",
      palette="cividis", linewidth=0.7, edgecolor="black"
      )
g.add_legend()

# Use semantically meaningful titles for the columns
titles = ["IDs per m/z", "m/z mapped to pathways", "pathways count", "m/z count per pathway"]

for ax, title in zip(g.axes.flat, titles):
    # Set a different title for each axes
    ax.set(title=title)
    # Make the grid horizontal instead of vertical
    ax.xaxis.grid(True)
    ax.yaxis.grid(True)

sns.despine(left=True, bottom=True)
#plt.savefig("../plot_wiki_vs_kegg.pdf")
plt.show()


## === bubble plots

df = pd.read_csv(f"../../figs/initial-ids-pathways-KEGG-{my_selection}.csv")  # open file
df = df.replace(0, np.nan)

df = update_names_TD_samples(df)
df = df.loc[df["sample_name"].isin(good_samples), :]


### bubble plots of the translation of all types of identifiers to KEGG compounds ids

print("bubble plots of the translation of all types of identifiers to KEGG compounds ids")
df_pl = transform2plot_v2(df, mode="all")
size_range_markers1 = (df_pl['amount'].min() / 10,
                       df_pl['amount'].max() / 10)

df_pl2 = transform2plot_v2(df, mode="ids")
size_range_markers2 = (df_pl2['amount'].min() / 10,
                       df_pl2['amount'].max() / 10)

bubble_plot_kegg_metabolite_id_conversion_success(df_pl, df_pl2,
                                                  size_range_markers1,
                                                  size_range_markers2,
                                                  fig_width=my_width,
                                                  fig_height=my_height)

###  bubble plots of the final mapping of all types of identifiers to KEGG pathways

print("bubble plots of the final mapping of all types of identifiers to KEGG pathways")
df_pl3 = transform2plot_v2(df, mode="pathways")
size_range_markers3 = (df_pl3['amount'].min() / 10, df_pl3['amount'].max() / 10)
bubble_plot_for_pathways(df_pl, df_pl3, size_range_markers1,
                         size_range_markers3,
                         f"../../figs/initial-ids-pathways-KEGG-{my_selection}.pdf",
                         fig_width=my_width,
                         fig_height=my_height)


### WIKIPATHWAYS bubble plots
df = pd.read_csv(f"../../figs/initial-IDSmix-pathways-psize1-{my_selection}.csv")  # open file
df = df.replace(0, np.nan)
print(df.columns)
df = update_names_TD_samples(df)

df = df.loc[df["sample_name"].isin(good_samples), :]   # exclude bad_samples
df = update_names_TD_samples(df)

df_pl = transform2plot_v2(df, mode="identifiers")
size_range_markers1 = (df_pl['amount'].min() / 10, df_pl['amount'].max() / 10)
print(df_pl.head(10))
print(df_pl['amount'].min())

df_pl2 = transform2plot_v2(df, mode="pathways")
size_range_markers2 = (df_pl2['amount'].min() / 10, df_pl2['amount'].max() / 10)

bubble_plot_for_pathways(df_pl, df_pl2, size_range_markers1,
                         size_range_markers2,
                         f"../../figs/initial-ids-pws-WIKIPATHWAYS-{my_selection}.pdf",
                         fig_width=my_width,
                         fig_height=my_height)



print("done")


# ancient chunk
#df = df.loc[df['ratio_id_mz'] <= 2.2, : ].copy()
#df = df.loc[df['pathways_number_final'] >= 5, :].copy()
#df = df.loc[df['ratio_mz_pathway_count'] >=0.7,: ].copy()
#df = df.loc[df['mz_mapped_to_pathways'] >= 5,: ].copy()

# ancient chunk latest
# bad_samples = df.loc[(df['ratio_id_mz'] > 2.2) | (
#    df['pathways_number_final'] < 10) | (df['mz_mapped_to_pathways'] < 10) , "sample_name" ].unique()