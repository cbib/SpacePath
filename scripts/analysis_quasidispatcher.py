import os
import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc

from typing import Dict, List, Union, Tuple

from scipy.sparse import csr_matrix
from smpath.data import MetabolicPathways, MetabolomeAnnotations
import smpath.method as spm
from smpath.method import ComputePathwaysScoresConfig

import smpath.data as spd
from sklearn.cluster import KMeans

import pickle
from analysis_helpers import *
import smpath.processing.measures_methods as measures_methods
import matplotlib.pyplot as plt

import warnings
warnings.simplefilter('ignore',category=FutureWarning)
import squidpy as sq


np.random.seed(42)  # important to make reproducible the computations

def display_minireport(top_differential, SAMPLE):
    all_l = top_differential.names.unique()
    all_fam = [k for k in all_l if k.startswith("MolFam")]
    wiki_p = [k for k in all_l if k.startswith("WP")]
    kegg_p = [k for k in all_l if k.startswith("hsa")]
    print("SAMPLE:", SAMPLE)
    print("   differential families:", len(all_fam))
    print("   differential KEGG pathways:", len(kegg_p))
    print("   differential WIKIPATHWAYS:", len(wiki_p))


def compute_differential_families_and_pathways(pathways_types,
        clustering_parameters, SAMPLE, subfolder, data_centralized_folder):
    pathways_types += ["FAMILIES"]
    str_out = f"{clustering_parameters['choice']}-{clustering_parameters['numeric_parameter']}"

    list_of_scores_objects =  list()
    for u in pathways_types:
        adata_scores = ad.read_h5ad(filename=os.path.join(data_centralized_folder, subfolder,
                         f'{SAMPLE}-scores-{u}.h5ad'))
        adata_scores.var['set_type'] = u
        list_of_scores_objects.append(adata_scores)
        del adata_scores
    joined = ad.concat(adatas=list_of_scores_objects, join='outer', axis=1)

    clusters_df = pd.read_csv(os.path.join(data_centralized_folder, subfolder,
        f'{SAMPLE}-clusters-{str_out}.txt'), header=0)
    joined.obs[str_out] = pd.Categorical(clusters_df[str_out].tolist(),
                                   categories=clusters_df[
                                       str_out].unique())
    joined.X = csr_matrix(joined.X.todense() + abs(joined.X.min().min()) ) #  shift to avoid negative values
    top_differential = extract_top_differential_features(joined, group_str=str_out,
                                                         ntop_per_cluster=10)
    top_differential.to_csv(os.path.join(data_centralized_folder, subfolder,
                                         f'{SAMPLE}-top-differential-SETS.tsv'),
                            sep='\t', header=True, index=False)
    display_minireport(top_differential, SAMPLE)


def compute_families_scores( min_family_size,
          SAMPLE, subfolder, data_centralized_folder):

    families_obj = spd.MetabolicPathways(
        data_dir=data_centralized_folder,
        subfolder=subfolder,
        file=f"{SAMPLE}-molecules_families.mmt"
    )
    families_obj.molecule_id_type = 'moleculeIds'
    families_obj.load()

    adata_mz = safe_load_processed_data(SAMPLE, subfolder,
                                        data_centralized_folder)

    clean_annots_wikipathways = pd.read_csv(
        os.path.join(data_centralized_folder, subfolder,
                     f'{SAMPLE}-WIKI-CleanAnnotations.tsv'),
        sep='\t', header=0, index_col=None)

    families_obj, annots_df = spm.SyncPathwaysAndAnnots.run(
        families_obj.model_copy(),
        clean_annots_wikipathways,
        usable_id='moleculeIds',
        min_set_size=min_family_size)
    ##
    metabolites_computed_df, coords = spm.ComputeAssignmentsWeights.run(
        adata_mz, annots_df, usable_id='moleculeIds' )
    del adata_mz  # free some space, we can reload it later

    scores_adata = run_pathway_scores_steps_produce_anndata(metabolites_computed_df, coords,
                                             families_obj,
                                             min_set_size=min_family_size,
                                             type_abbr="FAMILIES"
         )
    out_file = os.path.join(data_centralized_folder, subfolder,
                         f'{SAMPLE}-scores-FAMILIES.h5ad')
    scores_adata.write_h5ad(out_file)



def compute_molecules_families( cutoff_tree,
                               SAMPLE, subfolder, data_centralized_folder,  plot_folder):
    top_differential = pd.read_csv(
        os.path.join(data_centralized_folder, subfolder,
                     f'{SAMPLE}-top-differential-mzs.tsv'), sep='\t', header=0)

    clean_annots_nokegg_ids = pd.read_csv(
        os.path.join(data_centralized_folder, subfolder,
                     f'{SAMPLE}-WIKI-CleanAnnotations.tsv'),
        sep='\t', header=0, index_col=None)

    molecules_ids_df = prepare_mols_df(top_differential,
                                       clean_ann_df=clean_annots_nokegg_ids)

    mol_similarities_result = spm.ComputeMoleculesSimilarityMatrix.run(
        molecules_ids_list=molecules_ids_df["moleculeIds"])

    with open(os.path.join(data_centralized_folder, subfolder,
                           f"{SAMPLE}-mol-similarities-matrix.pkl"),
              "wb") as f:
        pickle.dump(mol_similarities_result, f)

    updated_config = spm.ComputeMoleculesFamiliesConfig(cutoff_tree=cutoff_tree)
    computed_fam_results = spm.ComputeMoleculesFamilies.run(
        mol_similarities_result,
        molecules_ids_df,
        SAMPLE,
        subfolder,
        data_centralized_folder,
        plot_folder,
        config=updated_config)  # all results: maybe of relevance for the user

    with open(os.path.join(
            data_centralized_folder, subfolder, f"{SAMPLE}-mol-fam-dict.pkl"),
            "wb") as f:
        pickle.dump(computed_fam_results, f)



def compute_pathway_scores( pathways_types, min_pathway_size,
         SAMPLE, subfolder, data_centralized_folder ):
    """doctstring here?"""

    sets_args_dict = files_of_pathways_defaults()
    #mapped_elements = list()
    for pwtype in pathways_types:
        adata_mz, clean_annots_dict = safe_load_data_and_clean_annotations(
            SAMPLE, pathways_types, subfolder, data_centralized_folder
        )
        annots_df = clean_annots_dict[pwtype]
        pathways_obj = spd.MetabolicPathways(
            data_dir=data_centralized_folder,
            subfolder=sets_args_dict[pwtype]['subfolder'],
            file=sets_args_dict[pwtype]['file']
        )
        pathways_obj.molecule_id_type = sets_args_dict[pwtype]['id_type']

        pathways_obj.load()
        pathways_obj, annots_df = spm.SyncPathwaysAndAnnots.run(
            pathways_obj.model_copy(),
            annots_df,
            usable_id=sets_args_dict[pwtype]['id_type'],
            min_set_size=min_pathway_size)

        out_file = os.path.join(data_centralized_folder, subfolder,
                             f'{SAMPLE}-scores-{pwtype}.h5ad')

        metabolites_computed_df, coords = spm.ComputeAssignmentsWeights.run(adata_mz, annots_df,
                                   usable_id=sets_args_dict[pwtype]['id_type'])
        del adata_mz  # free some space, we  reload it

        scores_adata = run_pathway_scores_steps_produce_anndata(
            metabolites_computed_df, coords,  pathways_obj,
            min_set_size=min_pathway_size,
             type_abbr=pwtype  )

        scores_adata.write_h5ad(out_file); print("success, saved scores\n")
        #mapped_elements += scores_adata.uns['current_pathways_df']['mz'].unique(
        #).tolist()


def compute_differential_analysis_of_mzs(pathways_types,
    clustering_parameters, SAMPLE,  subfolder, data_centralized_folder
):
    adata_mz, clean_annots_dict = safe_load_data_and_clean_annotations(
        SAMPLE, pathways_types, subfolder, data_centralized_folder
    )
    str_out = f"{clustering_parameters['choice']}-{clustering_parameters['numeric_parameter']}"
    clusters_df = pd.read_csv(os.path.join(data_centralized_folder, subfolder,
        f'{SAMPLE}-clusters-{str_out}.txt'), header=0)
    adata_mz.obs[str_out] = pd.Categorical(clusters_df[str_out].tolist(),
                                           categories=clusters_df[str_out].unique())
    sc.pp.log1p(adata_mz)
    print(
        "computing the differential abundances of all the m/z features, cluster to cluster, for all clusters")
    top_differential = extract_top_differential_features(adata_mz, group_str=str_out,
                                                         ntop_per_cluster=50)
    top_differential.to_csv(os.path.join(data_centralized_folder, subfolder,
            f'{SAMPLE}-top-differential-mzs.tsv'),
                              sep='\t', header=True, index=False)


def plot_clustering_results(adata_mz, clustering_parameters,
                             SAMPLE, plot_folder):
    str_out = f"{clustering_parameters['choice']}-{clustering_parameters['numeric_parameter']}"
    plot_out_path = os.path.join(plot_folder,
                                 f'{SAMPLE}-clusters-{str_out}.png')
    sq.pl.spatial_scatter(
        adata_mz, color=[str_out], img=False,
        title="Clustering of ions of good quality",
        linewidth=0, edgecolors=(0,0,0,0),
        shape="circle",  #  'circle', 'square', 'hex'
        save=os.path.abspath(plot_out_path),
        dpi=600)
    plt.close()
    print("saved figure")


def compute_leiden__weighted_by_spatial_graph(
    datacopy, resolution_param, alpha_gr=0.25
) -> Tuple[ad.AnnData, List[str]]:
    sc.pp.neighbors(datacopy)
    nn_graph_genes = datacopy.obsp["connectivities"]
    sq.gr.spatial_neighbors(datacopy)
    nn_graph_space = datacopy.obsp["spatial_connectivities"]

    alpha = alpha_gr  # the weight of the spatial graph
    joint_graph = (1 - alpha) * nn_graph_genes + alpha * nn_graph_space
    sc.tl.leiden(datacopy,
                 adjacency=joint_graph,
                 n_iterations=2,
                 flavor="igraph",
                 directed=False, key_added="Leiden_sg",
                 resolution=resolution_param)
    the_labels_str_l = datacopy.obs['Leiden_sg'].astype(str).tolist()
    return datacopy, the_labels_str_l


def compute_clustering_of_spots(pathways_types, clustering_parameters,
                                SAMPLE, subfolder, data_centralized_folder,  plot_folder ):
    # subset to agree with the clean annotations
    adata_mz, clean_annots_dict = safe_load_data_and_clean_annotations(
          SAMPLE, pathways_types, subfolder, data_centralized_folder
    )
    top_moran_i = adata_mz.var.loc[
         adata_mz.var['I'] >= 0.1, "mz"].astype(str).tolist()
    datacopy = adata_mz[:,top_moran_i].copy()  # only do the clustering with this selection ? yes
    #datacopy = adata_mz.copy()
    sc.pp.log1p(datacopy)
    sc.tl.pca(datacopy, svd_solver="arpack")
    # adata_mz.obsm['X_pca'] = datacopy.obsm['X_pca']
    if clustering_parameters['choice'] == "kmeans":
        user_k = int(clustering_parameters['numeric_parameter'])
        kmeans = KMeans(n_clusters=user_k, random_state=0).fit(datacopy.obsm['X_pca'])
        the_labels_str_l = kmeans.labels_.astype(str)
    elif clustering_parameters['choice'] == "leiden_sg":
        resolution = float(clustering_parameters['numeric_parameter'])
        datacopy, the_labels_str_l = compute_leiden__weighted_by_spatial_graph(
            datacopy, resolution)

    str_out = f"{clustering_parameters['choice']}-{clustering_parameters['numeric_parameter']}"
    adata_mz.obs[str_out] = the_labels_str_l
    adata_mz.obs[["obs", str_out]].to_csv(
        os.path.join(data_centralized_folder, subfolder,
              f'{SAMPLE}-clusters-{str_out}.txt'), index=False)
    plot_clustering_results(adata_mz, clustering_parameters, SAMPLE,
                             plot_folder)


def clean_annotations_and_save_tsv(annots_obj, pathways_types,
                                   SAMPLE, subfolder, data_centralized_folder
                                   )->None:
    """
    cleans the loaded annotations using smpath library,
      then saves to tab delimited .tsv files.
      one file for each type of pathway (KEGG or WIKIPATHWAYS) is saved
    """
    adata_mz = safe_load_processed_data(SAMPLE, subfolder, data_centralized_folder)
    annots_obj, adata_mz = measures_methods.match_adata_mz_and_annots(annots_obj, adata_mz)

    # object required to translate annotations to KEGG ids when required
    ref_KEGG_ids_obj = spd.MetaboliteIdsMappings(
        data_dir=data_centralized_folder,
        subfolder="metabolites-mappings",
        mapping_files={'HMDB': "HMDB-to-KEGG_05-07-24.tsv",
                       'LipidMaps': "LipidMaps-to-KEGG_10-18-24.tsv",
                       'CoreMetabolome': "CoreMetabolome-to-KEGG_05-07-24.tsv"})
    ref_KEGG_ids_obj.load()

    # for each type of pathway, save the respective clean annotations
    for pwtype in pathways_types:
        annots_current = annots_obj.model_copy() # copy is compulsory
        # cleaning annotations using `smpath` library:
        if pwtype == "WIKI":
            annots_current.usable_id = "moleculeIds"
            annots_current.verify()
        elif pwtype == "KEGG":
            annots_current.usable_id = "KEGG"
            annots_current.convert_ids_to_kegg(ref_KEGG_ids_obj.table)

        annots_current.keep_best_mz_id_pairs(adata_mz.var.copy())
        annots_current.remove_redundant_aliases()

        # saving to tsv file
        out_file = os.path.join(data_centralized_folder, subfolder,
            f'{SAMPLE}-{pwtype}-CleanAnnotations.tsv')
        annots_current.df.to_csv(out_file, sep='\t', header=True, index=False)
        del annots_current
    # end for


def normalize_classify_and_save_anndata_mz(
        annots_obj, offSample_cutoff,
        msm_cutoff, fdr_cutoff, moran_I_cutoff,
        moran_padj_cutoff, SAMPLE, subfolder,
        data_centralized_folder, plot_folder) -> None:
    """
    performs the first step consisting of intensities' normalisation
    and the classification of the ions (m/z) by their quality metrics,
    using the library smpath
    """
    adata_mz = ad.read_h5ad(filename=os.path.join(
        data_centralized_folder, subfolder, f'{SAMPLE}.h5ad'))

    try:  # resets adata.var columns if user re-runs analysis  # todo: how to avoid using this in the future ?
        columns_to_reset = set(adata_mz.var.columns.tolist()) - set(["mz"])
        adata_mz.var = adata_mz.var.drop(columns=list(columns_to_reset))
    except KeyError:
        pass

    annots_obj, adata_mz = measures_methods.match_adata_mz_and_annots(annots_obj, adata_mz)
    adata_mz = spm.NormalizeIntensities.run(adata_mz)
    adata_mz = spm.ClassifyMZProtocol.run(
        adata_mz, annots_obj, msm_cutoff, offSample_cutoff, fdr_cutoff,
        moran_I_cutoff, moran_padj_cutoff, SAMPLE, plot_folder)

    adata_mz.write_h5ad(filename=os.path.join(
        data_centralized_folder, subfolder, f"{SAMPLE}.h5ad"))