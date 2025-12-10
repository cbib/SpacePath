#import smpath.data as spd
import pandas as pd
import numpy as np

my_selection = "public"
min_pathway_size = 1


df = pd.read_csv(f"../../figs/initial-ids-pathways-KEGG-{my_selection}.csv",
                 header=0, index_col=0)
df = df.replace(0, np.nan)

df2 = pd.read_csv(
    f"../../figs/initial-IDSmix-pathways-psize{min_pathway_size}-{my_selection}.csv",
header=0, index_col=0)
df2 = df2.replace(0, np.nan)


uu = 99
print(df.head())
print(df2.head())

## compute the global fraction mapped for KEGG pathways
global_KEGG_mapped =  (df['HMDB in pathways'].mean(
) + df['LipidMaps in pathways'].mean()
                       ) / (
        df['HMDB identifiers all'].mean() + df['LipidMaps identifiers all'].mean())

print("ids mapped KEGG : ", global_KEGG_mapped)

## compute the global fraction mapped for WIKIPATHWAYS
global_WIKI_mapped =  (df2['HMDB ids in pathways'].mean(
) + df2['LipidMaps ids in pathways'].mean()
                       ) / (
        df2['HMDB identifiers'].mean() + df2['LipidMaps identifiers'].mean())

print("ids mapped WIKI : ", global_WIKI_mapped)


global_both_databases = (global_KEGG_mapped + global_WIKI_mapped) / 2

unmapped_fraction = 1 - global_both_databases
print("===============")
print("unmapped fraction overall: ", unmapped_fraction)
print("===============")



