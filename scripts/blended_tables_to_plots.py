import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


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
        df_pl, df_pl2, size_range_markers1, size_range_markers2) -> None:

    sns.set_style("whitegrid")
    plt.figure(figsize=(6, 11))
    ax = sns.scatterplot(df_pl, x="database", y="sample_name", size="amount",
                         lw=0.4, sizes=size_range_markers1,
                         # sizes=(200, 1000),
                         color="gray", alpha=0.1,
                         edgecolors="black")
    ax.grid(color='grey', linestyle='--', alpha=0.4, zorder=1)
    sns.scatterplot(df_pl2, ax=ax,
                    x="database", y="sample_name", size="amount",
                    lw=0.4, sizes=size_range_markers2,  # sizes=(15, 60),
                    color="sandybrown", alpha=0.7,
                    edgecolors="black").legend(loc="upper left", title="size",
                                               bbox_to_anchor=(1,
                                                               1.03))  # .02))
    plt.margins(x=1)
    sns.despine(fig=None, ax=None,
                top=True, right=True, left=True, bottom=True, offset=None,
                trim=False)
    plt.tight_layout()
    plt.savefig("../figs/initial-ids-compounds-KEGG.pdf")
    plt.close()


def bubble_plot_for_pathways(df_pl, df_pl3, size_range_markers1,
                             size_range_markers3, out_path) -> None:
    sns.set_style("whitegrid")
    plt.figure(figsize=(6, 11))
    ax = sns.scatterplot(df_pl, x="database", y="sample_name", size="amount",
                         lw=0.4, sizes=size_range_markers1,
                         # sizes=(200, 1000),
                         color="gray", alpha=0.1,
                         edgecolors="black")
    ax.grid(color='grey', linestyle='--', alpha=0.4, zorder=1)
    sns.scatterplot(df_pl3, ax=ax,
                    x="database", y="sample_name", size="amount",
                    lw=0.4, sizes=size_range_markers3,  # sizes=(15, 60),
                    color="darkcyan", alpha=0.7,
                    edgecolors="black").legend(loc="upper left", title="size",
                                               bbox_to_anchor=(1,
                                                               1.03))  # .02))
    plt.margins(x=1)
    sns.despine(fig=None, ax=None,
                top=True, right=True, left=True, bottom=True, offset=None,
                trim=False)
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()



## start

df = pd.read_csv("../figs/initial-ids-pathways-KEGG.csv")
df = df.replace(0, np.nan)



### bubble plots of the translation of all types of identifiers to KEGG compounds ids

print("bubble plots of the translation of all types of identifiers to KEGG compounds ids")
df_pl = transform2plot_v2(df, mode="all")
size_range_markers1 = (df_pl['amount'].min() / 10,
                       df_pl['amount'].max() / 10)

df_pl2 = transform2plot_v2(df, mode="ids")
size_range_markers2 = (df_pl2['amount'].min() / 10,
                       df_pl2['amount'].max() / 10)

bubble_plot_kegg_metabolite_id_conversion_success(df_pl, df_pl2, size_range_markers1, size_range_markers2)



###  bubble plots of the final mapping of all types of identifiers to KEGG pathways

print("bubble plots of the final mapping of all types of identifiers to KEGG pathways")
df_pl3 = transform2plot_v2(df, mode="pathways")
size_range_markers3 = (df_pl3['amount'].min() / 10, df_pl3['amount'].max() / 10)
bubble_plot_for_pathways(df_pl, df_pl3, size_range_markers1,
                         size_range_markers3, "../figs/initial-ids-pathways-KEGG.pdf")



### WIKIPATHWAYS bubble plots
df = pd.read_csv("../figs/initial-IDSmix-pathways-psize1.csv")
df = df.replace(0, np.nan)
print(df.columns)

df_pl = transform2plot_v2(df, mode="tifiers")
size_range_markers1 = (df_pl['amount'].min() / 10, df_pl['amount'].max() / 10)
print(df_pl.head(10))
print(df_pl['amount'].min())

df_pl2 = transform2plot_v2(df, mode="pathways")
size_range_markers2 = (df_pl2['amount'].min() / 10, df_pl2['amount'].max() / 10)

bubble_plot_for_pathways(df_pl, df_pl2, size_range_markers1,
                         size_range_markers2, "../figs/initial-ids-pws-WIKIPATHWAYS.pdf")



###  dot plots (as PairGrid) comparing KEGG and WIKIPATHWAYS

print("dot plots (as PairGrid) comparing KEGG and WIKIPATHWAYS")

df = pd.read_csv("../figs/final-KEGG-compounds-counts.csv", index_col=None)
df2 = pd.read_csv("../figs/final-IDSmix-compounds-counts-psize2.csv", index_col=None)

df = df[['sample_name', 'ratio_keggid_mz', 'mz_mapped_to_pathways', 'pathways_number_final']]
df2 = df2[['sample_name', 'ratio_id_mz', 'mz_mapped_to_pathways', 'pathways_number_final']]
df = df.replace(0, np.nan)
df.columns = ['sample_name', 'ratio_id_mz', 'mz_mapped_to_pathways', 'pathways_number_final']

df['PType'] = "KEGG_pathways"
df2['PType'] = "Wikipathways"

df = pd.concat([df, df2])

df["ratio_mz_pathway_count"] = (df['mz_mapped_to_pathways'].to_numpy(
) + 0.001) / (df['pathways_number_final'].to_numpy() + 0.001)

df = df[['sample_name', 'ratio_id_mz', 'mz_mapped_to_pathways', 'pathways_number_final', 'ratio_mz_pathway_count', 'PType']]

bad_samples = df.loc[(df['ratio_id_mz'] > 2.25) | (
    df['pathways_number_final'] < 10) | (df['mz_mapped_to_pathways'] < 10) , "sample_name" ].unique()

print("bad tissues:\n", bad_samples)

df = df.loc[~df["sample_name"].isin(bad_samples), :]

kk = df['sample_name'] + '&' +df['PType']
df.index = kk

#df = df.loc[order_samples,: ]
df['data_type'] = "Lipidomics"
df.loc[df['sample_name'].str.startswith("MALDI"), "data_type"] = "SmallMolecules"
datasmalls = ["TD-635-1-A", "TD-hemi-635-1-B", "LP-heart-mouse-2"]
df.loc[df['sample_name'].isin(datasmalls), "data_type"] = "SmallMolecules"


sns.set_theme(style="whitegrid")
#plt.figure(figsize=(60, 8))   # very difficult to set the dimensions, instead saving manually from display window

g = sns.PairGrid(df,
                 x_vars=df.columns[1:5], y_vars="sample_name",
                 height=10, aspect=.25)

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

print("done")


# ancient chunk
#df = df.loc[df['ratio_id_mz'] <= 2.2, : ].copy()
#df = df.loc[df['pathways_number_final'] >= 5, :].copy()
#df = df.loc[df['ratio_mz_pathway_count'] >=0.7,: ].copy()
#df = df.loc[df['mz_mapped_to_pathways'] >= 5,: ].copy()