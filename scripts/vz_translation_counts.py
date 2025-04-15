import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


df = pd.read_csv("../final_compounds_counts.csv", index_col=None)
df = df.replace(0, np.nan)
# out_df = pd.DataFrame({
#     "sample_name": sample_name,
#     "mz_mapped_to_pathways": mz_mapped_to_pathways,
#     "ratio_keggid_mz": kegg_compound_ids_to_mz_ratio_after_pathway_mapped,
#     "": how_many_pathways})
# reorder
df = df[['sample_name', 'ratio_keggid_mz', 'mz_mapped_to_pathways', 'pathways_number_final']]

df["ratio_mz_pathway_count"] = df['mz_mapped_to_pathways'].to_numpy(
) / df['pathways_number_final'].to_numpy()

# df['quality_metric'] = (df["ratio_mz_pathway_count"] / df['ratio_keggid_mz']) * (
#     df['pathways_number_final'] * 1/10)
# df = df.sort_values("quality_metric", ascending=False)

order_samples = ["MET-breast-cancer-1", "MET-skin-cancer-1",
                 "MET-lung", "MET-stomach", "MET-gut", "MET-skin-cancer-2",
                 "MET-brain-mouse", "MET-brain-rat", "MET-eye",
                  "MET-heart", "UBcalibr2",
                  "MET-kidney2","MET-kidney1", "MET-liver-cancer",
                 "UBcalibr1", "UBcalibr4",  "MALDI_259_T",
                 "MALDI_269_T", "MALDI_275_T", "MALDI_262_T", "MALDI_260_T",  "MALDI_248_T"]
df.index = df['sample_name'].tolist()
df = df.loc[order_samples,: ]
df['data_type'] = "Lipidomics"
df.loc[df['sample_name'].str.startswith("MALDI"), "data_type"] = "SmallMolecules"
df.loc[df['sample_name'] == "UBcalibr1", "data_type"] = "SmallMolecules"
df.loc[df['sample_name'] == "UBcalibr4", "data_type"] = "SmallMolecules"


sns.set_theme(style="whitegrid")

# Make the PairGrid
g = sns.PairGrid(df,
                 x_vars=df.columns[1:5], y_vars="sample_name",
                 height=10, aspect=.25)

# Draw a dot plot using the stripplot function
g.map(sns.stripplot, data=df,  size=10, orient="h", jitter=False,
hue="data_type",
      palette="muted", linewidth=0.7, edgecolor="black"
      )
g.add_legend()

# Use the same x axis limits on all columns and add better labels
#g.set(xlim=(0, 25), xlabel="hohoho", ylabel="")

# Use semantically meaningful titles for the columns
titles = ["KEGGids per m/z", "m/z mapped to pathways", "pathways count", "m/z count per pathway", "new metric"]

for ax, title in zip(g.axes.flat, titles):

    # Set a different title for each axes
    ax.set(title=title)

    # Make the grid horizontal instead of vertical
    ax.xaxis.grid(True)
    ax.yaxis.grid(True)

sns.despine(left=True, bottom=True)
plt.show()