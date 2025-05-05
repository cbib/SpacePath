#
# select only the datasets that have at least 10 pathways,  less than 2.2
# of ambiguities, and a ratio superior to 1.2 ?



import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

df = pd.read_csv("../final-KEGG-compounds-counts.csv", index_col=None)
df2 = pd.read_csv("../final-IDSmix-compounds-counts-psize2.csv", index_col=None)

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


df = df.loc[df['ratio_id_mz'] <= 2.2, : ].copy()
df = df.loc[df['pathways_number_final'] >= 10, :].copy()
#df = df.loc[df['ratio_mz_pathway_count'] >=0.7,: ].copy()
df = df.loc[df['mz_mapped_to_pathways'] >= 20,: ].copy()


# df['quality_metric'] = (df["ratio_mz_pathway_count"] / df['ratio_keggid_mz']) * (
#     df['pathways_number_final'] * 1/10)
# df = df.sort_values("quality_metric", ascending=False)

kk = df['sample_name'] + '&' +df['PType']
df.index = kk

#df = df.loc[order_samples,: ]
df['data_type'] = "Lipidomics"
df.loc[df['sample_name'].str.startswith("MALDI"), "data_type"] = "SmallMolecules"
datasmalls = ["TD-635-1-A", "TD-hemi-635-1-B", "LP-heart-mouse-2"]
df.loc[df['sample_name'].isin(datasmalls), "data_type"] = "SmallMolecules"

print(df)

sns.set_theme(style="whitegrid")

# Make the PairGrid
g = sns.PairGrid(df,
                 x_vars=df.columns[1:5], y_vars="sample_name",
                 height=10, aspect=.25)

# Draw a dot plot using the stripplot function
g.map(sns.stripplot, data=df,  size=10, orient="h", jitter=False,
hue="PType",
      palette="cividis", linewidth=0.7, edgecolor="black"
      )
g.add_legend()

# Use the same x axis limits on all columns and add better labels
#g.set(xlim=(0, 25), xlabel="hohoho", ylabel="")

# Use semantically meaningful titles for the columns
titles = ["IDs per m/z", "m/z mapped to pathways", "pathways count", "m/z count per pathway", "new metric"]

for ax, title in zip(g.axes.flat, titles):

    # Set a different title for each axes
    ax.set(title=title)

    # Make the grid horizontal instead of vertical
    ax.xaxis.grid(True)
    ax.yaxis.grid(True)

sns.despine(left=True, bottom=True)
plt.show()