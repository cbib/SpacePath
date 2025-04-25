import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def transform2plot(df, mode="all"):
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
    df_pl = pd.concat([df1, df2], axis=0)
    return df_pl


df = pd.read_csv("../initial-ids-compounds-KEGG.csv")
df = df.replace(0, np.nan)

df_pl = transform2plot(df, mode="all")
size_range_markers1 = (df_pl['amount'].min() / 10, df_pl['amount'].max() / 10)
print(df_pl.head(10))
print(df_pl['amount'].min())

df_pl2 = transform2plot(df, mode="ids")
size_range_markers2 = (df_pl2['amount'].min() / 10, df_pl2['amount'].max() / 10)


sns.set_style("whitegrid")
plt.figure(figsize=(4.8 , 11))
ax = sns.scatterplot(df_pl, x="database", y="sample_name", size="amount",
                lw=0.4, sizes = size_range_markers1, #sizes=(200, 1000),
                color="gray", alpha=0.1,
                edgecolors="black")
ax.grid(color='grey', linestyle='--', alpha=0.4, zorder=1)
sns.scatterplot(df_pl2, ax= ax,
                x="database", y="sample_name", size="amount",
                lw=0.4, sizes=size_range_markers2, # sizes=(15, 60),
                color="darkcyan", alpha=0.7,
                edgecolors="black").legend(loc="upper left", title="size",
                             bbox_to_anchor=(1,1.03))  #.02))
plt.margins(x=1)
sns.despine(fig=None, ax=None,
            top=True, right=True, left=True, bottom=True, offset=None, trim=False)
plt.tight_layout()
plt.savefig("../initial-ids-compounds-KEGG.pdf")