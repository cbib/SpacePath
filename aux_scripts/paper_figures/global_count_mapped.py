import smpath.data as spd
import pandas as pd
import numpy as np
import os

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


####### now test how many were now associated to families
total_nb = df2['HMDB identifiers'].mean() + df2['LipidMaps identifiers'].mean()

print("total number of hmdb and lm IDs", total_nb)
set_ids_all_families = set()
data_dir = "../../../data"

my_selection = "public"

file_n = f"../../datasets_parameters/{my_selection}-datasets-parameters.tsv"
datasets_parameters_df = pd.read_csv(file_n, header=0,
                                     index_col=None, sep='\t')
for u in datasets_parameters_df['dataset_name'].tolist():
    print(os.path.exists(os.path.join(data_dir, u, f'{u}-molecules_families.mmt')))
    pathwaysObj = spd.MetabolicPathways(
        data_dir=data_dir,
        subfolder=u,
        file=f"{u}-molecules_families.mmt"
    )
    pathwaysObj.molecule_id_type = "moleculeIds"
    pathwaysObj.load()
    tmp_dict = pathwaysObj.to_dict()
    for innerkey in tmp_dict:
        set_ids_all_families = set_ids_all_families.union(set(tmp_dict[innerkey]))

print("overall total unique IDs in families:", len(set_ids_all_families))

print("fraction families relative to all ids: ", total_nb / len(set_ids_all_families))

