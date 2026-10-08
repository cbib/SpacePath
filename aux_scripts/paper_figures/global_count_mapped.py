import SpacePath.core.data as spd
from typing import List
import pandas as pd
import numpy as np
import os
import anndata as ad

data_types = ["public", "private"]
sets_coverage_dict= {}
data_dir = "../../../apollo-data"


def kegg_id_retrieve_origin_list(kegg_ids: List[str], good_annots_kegg: pd.DataFrame):
    foo = good_annots_kegg.loc[
        good_annots_kegg['KEGG'].isin(kegg_ids),
        "moleculeIds"].tolist()

    return foo


for data_type in data_types:
    print("opening the", data_type, "collection and secondary associated tables")
    # reuse the tables produced for the figures of coverage of the paper
    df = pd.read_csv(f"../../figs/mapp-ids-pathways-KEGG-{data_type}.csv",
                     header=0, index_col=0)
    df = df.replace(0, np.nan)

    df2 = pd.read_csv(
        f"../../figs/mapp-IDSmix-pathways-{data_type}.csv",
    header=0, index_col=0)
    df2 = df2.replace(0, np.nan)

    ## compute the global fraction mapped for KEGG pathways
    global_KEGG_mapped =  (df['HMDB in pathways'].mean(
    ) + df['LipidMaps in pathways'].mean()
                           ) / (
            df['HMDB identifiers all'].mean() + df['LipidMaps identifiers all'].mean())

    print("ids mapped KEGG (fraction): ", global_KEGG_mapped)

    ## compute the global fraction mapped for WIKIPATHWAYS
    global_WIKI_mapped =  (df2['HMDB ids in pathways'].mean(
    ) + df2['LipidMaps ids in pathways'].mean()
                           ) / (
            df2['HMDB identifiers'].mean() + df2['LipidMaps identifiers'].mean())

    print("ids mapped WIKI (fraction): ", global_WIKI_mapped)

    canonical_coverage_avg = (global_KEGG_mapped + global_WIKI_mapped) / 2

    sets_coverage_dict[data_type] = {'cov canonical' : canonical_coverage_avg,
                                     'cov canonical KEGG': global_KEGG_mapped,
                                     'cov canonical WIKI': global_WIKI_mapped}

    ####### now test how many were now associated to families or pathways or both
    fractions_coverage_families_list = list()
    fractions_cov_pathways_list = list()
    fractions_coverage_all_sets = list()
    coverage_pathways_list = list()
    set_ids_all_families = set()

    file_n = f"../../datasets_parameters/{data_type}-datasets-parameters.tsv"
    datasets_parameters_df = pd.read_csv(file_n, header=0,
                                         index_col=None, sep='\t')

    for u in datasets_parameters_df['dataset_name'].tolist():
        families_effective = set()
        good_annots_wiki = pd.read_csv(
            os.path.join(data_dir, u, f'{u}-WIKI-CleanAnnotations.tsv'), sep='\t', index_col=0)

        good_annots_kegg = pd.read_csv(
            os.path.join(data_dir, u, f'{u}-KEGG-CleanAnnotations.tsv'), sep='\t', index_col=0)

        local_total_list = good_annots_wiki['moleculeIds'].unique().tolist()

        try:
            fams_scores = ad.read_h5ad(os.path.join(data_dir, u, f'{u}-scores-FAMILIES.h5ad'))
            ids_in_fam = fams_scores.uns['current_pathways_df']['moleculeIds'].unique().tolist() #scores.var.pathway_id.unique().tolist()
        except FileNotFoundError:
            ids_in_fam = list()

        try:
            kegg_scores = ad.read_h5ad(os.path.join(data_dir, u, f'{u}-scores-KEGG.h5ad'))
            ids_in_kegg = kegg_scores.uns['current_pathways_df']['KEGG'].unique().tolist()
            ids_in_kegg = kegg_id_retrieve_origin_list(ids_in_kegg, good_annots_kegg)
        except FileNotFoundError:
            ids_in_kegg = list()

        try:
            wiki_scores = ad.read_h5ad(os.path.join(data_dir, u, f'{u}-scores-WIKI.h5ad'))
            ids_in_wiki = wiki_scores.uns['current_pathways_df']['moleculeIds'].unique().tolist()
        except FileNotFoundError:
            ids_in_wiki = list()

        fractions_coverage_families_list.append(
            len(ids_in_fam) / len(local_total_list)
        )

        fractions_cov_pathways_list.append(
            ((len(set(ids_in_kegg)) / len(local_total_list)) + (
                len(set(ids_in_wiki)) / len(local_total_list)
            )) / 2  # because two mapping systems, kegg and wiki
        )

        fractions_coverage_all_sets.append(
            len(set(ids_in_fam).union(set(ids_in_kegg)).union(set(ids_in_wiki))) / len(local_total_list)
        )

    #print("fractions of families coverage, one per dataset:", fractions_coverage_families_list)
    print("avg coverage by families alone: ", np.mean(np.array(fractions_coverage_families_list)))
    #print("fractions of pathways coverage, one per dataset:", fractions_cov_pathways_list)
    print("avg coverage by pathways: ", np.mean(np.array(fractions_cov_pathways_list)))
    #print("fractions of all sets coverage, one per dataset:", fractions_coverage_all_sets)
    print("avg coverage by all sets : ", np.mean(np.array(fractions_coverage_all_sets)))
    sets_coverage_dict[data_type]['cov all sets avg'] = np.mean(np.array(fractions_coverage_all_sets))
    sets_coverage_dict[data_type]['cov poo'] = np.mean(np.array(fractions_cov_pathways_list))
    sets_coverage_dict[data_type]['cov families avg'] = np.mean(np.array(fractions_coverage_families_list))

    print(" ==> end of ", data_type, "collection")
    # end of a very long for

# Final report over both types of data collections
print("\n=========== Final report ==============")
print("coverage by canonical pathways alone, max:",
      np.max([sets_coverage_dict["private"]['cov canonical'],
       sets_coverage_dict["public"]["cov canonical"]]))   #  max 0.0136

print("coverage by canonical pathways alone, avg over public and private data:",
      (sets_coverage_dict["private"]['cov canonical'] +
       sets_coverage_dict["public"]["cov canonical"]) / 2)

print("coverage by KEGG:",
      (sets_coverage_dict["private"]['cov canonical KEGG'] +
       sets_coverage_dict["public"]["cov canonical KEGG"]) / 2)

print("coverage by WIKI:",
      (sets_coverage_dict["private"]['cov canonical WIKI'] +
       sets_coverage_dict["public"]["cov canonical WIKI"]) / 2)


print("coverage by all sets: pathways and SpacePath families:",
      (sets_coverage_dict["private"]['cov all sets avg'] +
       sets_coverage_dict["public"]["cov all sets avg"]) / 2)
print("=======================================")


# =========== Final report ==============

# =======================================