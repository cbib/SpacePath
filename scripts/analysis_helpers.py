import os
import pandas as pd
import anndata as ad
import scanpy as sc
import zonist.data as zd
import zonist.method as zm
import zonist.processing.measures_methods as zmm
from local_helpers.local_helpers import provide_annotation_dict
from local_helpers.quantifs_methods import differential_results_filter
from local_helpers.quantifs_methods import run_pathway_scores_steps
from local_helpers.local_viz import customized_clustermap
import matplotlib.pyplot as plt
import seaborn as sns


def differential_routine(adata, group_str):
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
        differential_results_by_cluster)

    return top_differential


def safe_load_processed_data(SAMPLE, subfolder, data_centralized_folder):

    adata_mz = ad.read_h5ad(filename=os.path.join(
        data_centralized_folder, subfolder, f'{SAMPLE}.h5ad'))
    adata_mz = adata_mz[:,
               adata_mz.var.loc[adata_mz.var['final_quality'] == "good",
               "mz"].astype(str).tolist()]
    print(adata_mz.uns["history"])


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
    print("shape of computed simils:::::", tmp.shape)
    the_molecules_names = tmp['moleculeNames'].tolist()
    print(len(list(mol_similarities_result['the_ids'])))
    print(len(the_molecules_names))
    trunc_names = [i[:30] for i in the_molecules_names]

    arr_df = pd.DataFrame(data=mol_similarities_result['array'],
                          columns=trunc_names,
                          index=trunc_names)
    gc = customized_clustermap(arr_df, computed_dendrogram)
    plt.savefig(
        os.path.join(plot_folder, f"{SAMPLE}-hierarchy-molecules.png"),
        dpi=400)
    plt.close()


# todo: deprecated ?
def save_scores(adata_mz, annots_dict,
                pwtype, sets_args_dict, SAMPLE, data_dir,
                not_for_families: list=None, only_mz_good_moran=True):
    assert pwtype in ["KEGG", "WIKI", "FAMILIES"], "Error, type of sets not recognized"
    annots_df = annots_dict["WIKI"] # valid also for families

    if pwtype == "KEGG":
        annots_df = annots_dict["KEGG"]

    if only_mz_good_moran:
        print("taking into account only the annotations of mz having I >= 0 !!")
        okmoran = adata_mz.var.loc[adata_mz.var['I'] >= 0, 'mz'].astype(str).tolist()
        annots_df = annots_df.loc[annots_df['mz'].astype(str).isin(okmoran),:].copy()

    sets_obj = zd.MetabolicPathways(
        data_dir=data_dir,
        subfolder=sets_args_dict[pwtype]['subfolder'],
        file=sets_args_dict[pwtype]['file']
    )
    sets_obj.molecule_id_type = sets_args_dict[pwtype]['id_type']

    sets_obj.load()

    if pwtype == "FAMILIES":
        sets_obj.intersect_with_annotations(annots_df,
                                            usable_id=sets_args_dict[pwtype]['id_type'],
                                            mz_col='mz')
        if not_for_families is not None:
            sets_obj.eject_molecules(not_for_families, column="mz")

    sets_obj, annots_df = zm.SyncPathwaysAndAnnots.run(
        sets_obj.model_copy(),
        annots_df,
        usable_id=sets_args_dict[pwtype]['id_type'],
        min_set_size=sets_args_dict[pwtype]['min_set_size'],
        mode="mz_aware")

    ofipw = os.path.join(data_centralized_folder, subfolder,
                         f'{SAMPLE}-scores-{pwtype}.h5ad')

    scores_adata = run_pathway_scores_steps(adata_mz.copy(),
                                            sets_obj,
                                            annots_df,
                                            usable_id=sets_args_dict[pwtype]['id_type'])
    scores_adata.write_h5ad(ofipw)
    print("youpi, saved scores")
    mapped_elements = scores_adata.uns['current_pathways_df']['mz'].unique(
        ).tolist()

    return mapped_elements

# def reference_noKEGG_KEGG_identifiers():
#     # stored in the same centralized data folder (same level of the tissues' folders)
#     subfolder="metabolites-mappings",
#     mapping_files={'HMDB': "HMDB-to-KEGG_05-07-24.tsv",
#                    'LipidMaps': "LipidMaps-to-KEGG_10-18-24.tsv",
#                    'CoreMetabolome': "CoreMetabolome-to-KEGG_05-07-24.tsv"})
#
#     return subfolder, mapping_files

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
