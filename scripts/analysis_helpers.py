import os
import pandas as pd
import anndata as ad
import scanpy as sc
from myutils.quantifs_methods import differential_results_filter
from myutils.local_viz import customized_clustermap
import matplotlib.pyplot as plt


def differential_routine(adata, group_str, ntop_per_cluster):
    sc.tl.rank_genes_groups(adata, groupby=group_str,
                            groups='all',
                            reference='rest',
                            method='wilcoxon',
                            key_added="top_feats_by_cluster",
                            corr_method='benjamini-hochberg')

    differential_results_by_cluster = sc.get.rank_genes_groups_df(
        adata, key="top_feats_by_cluster",
        group=None)  # returns all groups

    top_differential = differential_results_filter(
        differential_results_by_cluster, ntop_per_cluster=ntop_per_cluster)

    return top_differential


def safe_load_processed_data(SAMPLE, subfolder, data_centralized_folder):

    adata_mz = ad.read_h5ad(filename=os.path.join(
        data_centralized_folder, subfolder, f'{SAMPLE}.h5ad'))
    adata_mz = adata_mz[:,
               adata_mz.var.loc[adata_mz.var['final_quality'] == "good",
               "mz"].astype(str).tolist()]


    return adata_mz


def read_clean_annotations(SAMPLE, pathways_types, subfolder, data_centralized_folder):
    annots_dict = {}
    for pwtype in pathways_types:
        tmpfi = os.path.join(data_centralized_folder,
                             subfolder,
                             f'{SAMPLE}-{pwtype}-CleanAnnotations.tsv')
        annots_dict[pwtype] = pd.read_csv(tmpfi,
                                          sep='\t', header=0, index_col=None)

    return annots_dict


def safe_load_data_and_clean_annotations(SAMPLE,
                                         pathways_types,
                                         subfolder, data_centralized_folder):

    adata_mz = safe_load_processed_data(SAMPLE, subfolder, data_centralized_folder)
    annots_dict = read_clean_annotations(SAMPLE, pathways_types, subfolder, data_centralized_folder)

    adequate_mz = set()
    for pwtype in pathways_types:
        adequate_mz = adequate_mz.union(set(annots_dict[pwtype]['mz'].astype(str).unique()))

    adata_mz = adata_mz[:, list(adequate_mz)].copy()
    # if other_adequate_mzs is not None:
    #     adata_mz = adata_mz[: other_adequate_mzs].copy()

    return adata_mz, annots_dict


def clean_annotations(adata_mz, annots_obj, ref_KEGG_ids_obj, pwtype) -> None:

    if pwtype == "WIKI":
        annots_obj.usable_id = "moleculeIds"
    elif pwtype == "KEGG":
        annots_obj.usable_id = "KEGG"

    if pwtype == "WIKI":
        annots_obj.verify()
    elif pwtype == "KEGG":
        annots_obj.convert_ids_to_kegg(ref_KEGG_ids_obj.table)

    annots_obj.keep_best_mz_id_pairs(adata_mz.var.copy())

    if pwtype == "WIKI":
        annots_obj.remove_redundant_aliases()  # when using multitype 'moleculeIds' annotations

    return annots_obj


def prepare_mols_df(top_differential, clean_ann_df):
    top_differential['names'] = top_differential['names'].astype(float).tolist()
    molecules_ids_df = pd.merge(top_differential, clean_ann_df,
                                how='left', left_on='names', right_on='mz')
    molecules_ids_df.dropna(how='any', inplace=True)
    molecules_ids_df = molecules_ids_df[
        #['moleculeIds', 'moleculeNames']].copy()
        ['moleculeIds', 'mz', 'moleculeNames']].copy()
    molecules_ids_df.drop_duplicates(inplace=True)
    molecules_ids_df.reset_index(inplace=True)
    molecules_ids_df['moleculeNames'] = molecules_ids_df['moleculeNames'].astype(str).tolist()

    return molecules_ids_df


def message_ids_how_many_mz(the_ids, molecules_ids_df):
    molecules_ids_df = molecules_ids_df.loc[
        molecules_ids_df['moleculeIds'].isin(the_ids)
    ]
    return f"representing a total of mzs: {len(molecules_ids_df['mz'].unique())}"


def plot_families(molecules_ids_df, mol_similarities_result, computed_dendrogram,
                  plot_folder, SAMPLE):
    tmp = molecules_ids_df.loc[molecules_ids_df['moleculeIds'].isin(
        mol_similarities_result['the_ids']
    ), ["moleculeIds", "moleculeNames"]].copy()
    tmp.reset_index(inplace=True)
    tmp = tmp.drop_duplicates(subset="moleculeIds")
    tmp.index = tmp['moleculeIds']
    tmp = tmp.loc[mol_similarities_result['the_ids'],:].copy()
    the_molecules_names = tmp['moleculeNames'].tolist()
    trunc_names = [i[:30] for i in the_molecules_names]

    arr_df = pd.DataFrame(data=mol_similarities_result['array'],
                          columns=trunc_names,
                          index=trunc_names)
    gc = customized_clustermap(arr_df, computed_dendrogram)
    plt.savefig(
        os.path.join(plot_folder, f"{SAMPLE}-hierarchy-molecules.png"),
        dpi=400)
    plt.close()


def files_of_pathways_defaults():
    # stored in the same centralized data folder (same level of the tissues' folders)
    pathways_files_defaults = {
        "KEGG" : {'file': "KEGG_hsa_pathways_compounds_R110.gmt",
                      'subfolder': "metabolites-mappings",
                      'id_type': 'KEGG'},
        "WIKI" : {'file': "WIKIPATHWAYS_human_multi_id_redundant.mmt",
                          'subfolder': "metabolites-mappings",
                          'id_type': 'moleculeIds'}}
    return pathways_files_defaults


# functions to read the table of parameters

def provide_annotation_dict(SAMPLE, annotation_types_list):
    out_d = {}
    for annotation_type in annotation_types_list:
        out_d[annotation_type] = f"{SAMPLE}-{annotation_type}.csv"

    return out_d


def extract_params_per_dataset(row):
    SAMPLE = row["dataset_name"]
    annotations_to_use = row["annotations_to_use"].split(" ")
    annotations_dict = provide_annotation_dict(SAMPLE, annotations_to_use)
    offSample_cutoff = row["offSample_cutoff"]
    msm_cutoff = row["msm_cutoff"]
    fdr_cutoff = row["fdr_cutoff"]
    moran_I_cutoff = row["Moran_I_cutoff"]
    moran_padj_cutoff = row["Moran_padj_cutoff"]
    clustering_params = row['cluster_params'].split(":")
    clustering_parameters = {'choice': clustering_params[0].lower(),
                             'numeric_parameter': clustering_params[1]}
    pathways_types = row["pathways_types"].split(" ")
    min_pathway_size = row["min_pathway_size"]
    cutoff_tree = row['cut_tree_mol_families']
    min_family_size = row["min_family_size"]

    return (SAMPLE, annotations_dict, offSample_cutoff,
            msm_cutoff, fdr_cutoff, moran_I_cutoff,
            moran_padj_cutoff, clustering_parameters,
            pathways_types, min_pathway_size, cutoff_tree,
            min_family_size)



# def reference_noKEGG_KEGG_identifiers():
#     # stored in the same centralized data folder (same level of the tissues' folders)
#     subfolder="metabolites-mappings",
#     mapping_files={'HMDB': "HMDB-to-KEGG_05-07-24.tsv",
#                    'LipidMaps': "LipidMaps-to-KEGG_10-18-24.tsv",
#                    'CoreMetabolome': "CoreMetabolome-to-KEGG_05-07-24.tsv"})
#
#     return subfolder, mapping_files