import logging
import warnings

warnings.simplefilter('ignore', category=FutureWarning)

from Spacepath.core.helpers import setup_logger  # to add extra info to the log

import SpacePath.core.processing.measures_methods as measures_methods
from sklearn.cluster import KMeans

import pickle

from SpacePath.src.constants import (pathways_files_location_info,
                                     options_palette_supported)

from SpacePath.src.helpers import *

logger = logging.getLogger(__name__)
setup_logger(logger)

np.random.seed(42)  # important to make the computations be reproducible

""""Wrappers for a complete workflow as performed in our benchmark"""


def display_report_sets(top_differential, SAMPLE) -> None:
    all_l = top_differential.names.unique()
    all_fam = [k for k in all_l if k.startswith("MolFam")]
    wiki_p = [k for k in all_l if k.startswith("WP")]
    kegg_p = [k for k in all_l if k.startswith("hsa")]

    logger.info("SAMPLE: %s", SAMPLE)
    logger.info("   differential KEGG pathways: %s", len(kegg_p))
    logger.info("   differential WIKIPATHWAYS: %s", len(wiki_p))
    logger.info("   differential Families: %s", len(all_fam))


def compute_differential_families_and_pathways(
        pathways_types: List[str],
        clustering_parameters: Dict[str, str | int | float],
        SAMPLE: str, subfolder: str, data_centralized_folder: str,
        n_top_per_cluster: int = 30) -> None:  # before: n_top_per_cluster 25
    if not "FAMILIES" in pathways_types:
        pathways_types.append("FAMILIES")
    str_out = f"{clustering_parameters['choice']}-{clustering_parameters['numeric_parameter']}"

    scores_merged = join_scores_adatas(pathways_types, SAMPLE, subfolder,
                                       data_centralized_folder)

    clusters_df = pd.read_csv(os.path.join(data_centralized_folder, subfolder,
                                           f'{SAMPLE}-clusters-{str_out}.txt'),
                              header=0)

    scores_merged = safe_insert_clusters_column(scores_merged,
                                                clusters_df,
                                                str_out)
    #  shift to avoid negative values
    scores_merged.X = csr_matrix(
        scores_merged.X.todense() + abs(scores_merged.X.min().min()))

    differential_results_by_cluster = squidpy_differential_analysis(
        scores_merged, group_str=str_out)

    top_differential = differential_results_filter(
        differential_results_by_cluster, n_top_per_cluster=n_top_per_cluster)

    # new : moran filter for the sets
    scores_merged = measures_methods.compute_moran_i(scores_merged)
    scores_merged.var['I'] = scores_merged.uns['moranI'].loc[
        scores_merged.var.index.tolist(), 'I'].astype(float).tolist()
    top_sets_by_moran_i = scores_merged.uns['moranI'].loc[
        scores_merged.uns['moranI']['I'] >= 0.2, :].copy().index.tolist()
    top_differential = top_differential.loc[
        top_differential['names'].isin(top_sets_by_moran_i), :]

    top_differential = pd.merge(
        top_differential,
        scores_merged.var[['pathway_id', 'pathway_name', 'I']].copy(),
        how='left', left_on="names", right_on="pathway_id")

    top_differential.to_csv(os.path.join(data_centralized_folder, subfolder,
                                         f'{SAMPLE}-top-differential-SETS.tsv'),
                            sep='\t', header=True, index=False)
    display_report_sets(top_differential, SAMPLE)


def compute_families_scores(min_family_size,
                            SAMPLE, subfolder,
                            data_centralized_folder) -> None:
    """Compute families scores exclusively, and save to .h5ad file"""
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

    metabolites_computed_df, coords = spm.ComputeAssignmentsWeights.run(
        adata_mz, annots_df, usable_id='moleculeIds')
    del adata_mz  # free some space, we can reload it later

    scores_adata, time_elapsed = run_pathway_scores_steps_produce_anndata(
        metabolites_computed_df, coords, families_obj,
        min_set_size=min_family_size, type_abbr="FAMILIES"
    )
    out_file = os.path.join(data_centralized_folder, subfolder,
                            f'{SAMPLE}-scores-FAMILIES.h5ad')
    scores_adata.write_h5ad(out_file)


def compute_molecules_families(cutoff_tree, SAMPLE, subfolder,
                               data_centralized_folder, plot_folder) -> None:
    np.random.seed(42)
    top_differential_mzs = pd.read_csv(
        os.path.join(data_centralized_folder, subfolder,
                     f'{SAMPLE}-top-differential-mzs.tsv'), sep='\t',
        header=0)

    # to ensure consistency of families,  refilter for 25 top mz per cluster
    top_differential = differential_results_filter(top_differential_mzs,
                                                   n_top_per_cluster=25)

    clean_annots_nokegg_ids = pd.read_csv(
        os.path.join(data_centralized_folder, subfolder,
                     f'{SAMPLE}-WIKI-CleanAnnotations.tsv'),
        sep='\t', header=0, index_col=None)

    molecules_ids_df = prepare_molecular_ids_df(
        top_differential, clean_ann_df=clean_annots_nokegg_ids)

    mol_similarities_result = spm.ComputeMoleculesSimilarityMatrix.run(
        molecules_ids_list=molecules_ids_df["moleculeIds"])

    with open(os.path.join(data_centralized_folder, subfolder,
                           f"{SAMPLE}-mol-similarities-matrix.pkl"),
              "wb") as f:
        pickle.dump(mol_similarities_result, f)

    updated_config = spm.ComputeMoleculesFamiliesConfig(
        cutoff_tree=cutoff_tree)
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


def compute_pathway_scores(pathways_types, min_pathway_size,
                           SAMPLE, subfolder,
                           data_centralized_folder) -> None:
    """
    Compute pathways scores exclusively,
    saves each type of scores (e.g. KEGG or WIKI) to distinct .h5ad files
    """
    for pw_type in pathways_types:
        adata_mz, clean_annots_dict = safe_load_data_and_clean_annotations(
            SAMPLE, pathways_types, subfolder, data_centralized_folder
        )
        annots_df = clean_annots_dict[pw_type]
        pathways_obj = spd.MetabolicPathways(
            data_dir=data_centralized_folder,
            subfolder=pathways_files_location_info[pw_type]['subfolder'],
            file=pathways_files_location_info[pw_type]['file']
        )
        pathways_obj.molecule_id_type = pathways_files_location_info[pw_type][
            'id_type']

        pathways_obj.load()

        pathways_obj, annots_df = spm.SyncPathwaysAndAnnots.run(
            pathways_obj.model_copy(),
            annots_df,
            usable_id=pathways_files_location_info[pw_type]['id_type'],
            min_set_size=min_pathway_size)

        metabolites_computed_df, coords = spm.ComputeAssignmentsWeights.run(
            adata_mz, annots_df,
            usable_id=pathways_files_location_info[pw_type]['id_type'])
        del adata_mz  # free some space

        scores_adata, time_elapsed = run_pathway_scores_steps_produce_anndata(
            metabolites_computed_df, coords, pathways_obj,
            min_set_size=min_pathway_size,
            type_abbr=pw_type)
        logger.info(
            "success! scores for sample %s (set type %s) computed in seconds: %s",
            SAMPLE, pw_type, time_elapsed)

        out_file = os.path.join(data_centralized_folder, subfolder,
                                f'{SAMPLE}-scores-{pw_type}.h5ad')
        scores_adata.write_h5ad(out_file)
        logger.info("scores saved to file %s", out_file)


def compute_differential_analysis_of_mzs(pathways_types,
                                         clustering_parameters, SAMPLE,
                                         subfolder, data_centralized_folder,
                                         n_top_per_cluster=100  # new
                                         ) -> None:
    logger.info("computing differential ions between clusters")
    adata_mz, clean_annots_dict = safe_load_data_and_clean_annotations(
        SAMPLE, pathways_types, subfolder, data_centralized_folder
    )
    str_out = f"{clustering_parameters['choice']}-{clustering_parameters['numeric_parameter']}"
    clusters_df = pd.read_csv(os.path.join(data_centralized_folder, subfolder,
                                           f'{SAMPLE}-clusters-{str_out}.txt'),
                              header=0)

    adata_mz = safe_insert_clusters_column(adata_mz, clusters_df, str_out)

    sc.pp.log1p(adata_mz)
    logger.info("computing the differential abundances of all the m/z features,"
          " cluster to cluster, for all clusters")
    differential_results_by_cluster = squidpy_differential_analysis(
        adata_mz, group_str=str_out)

    top_differential = differential_results_filter(
        differential_results_by_cluster, n_top_per_cluster=n_top_per_cluster)

    top_differential['mz'] = top_differential['names'].tolist()
    top_differential.to_csv(os.path.join(data_centralized_folder, subfolder,
                                         f'{SAMPLE}-top-differential-mzs.tsv'),
                            sep='\t', header=True, index=False)


def compute_clustering_of_spots(
        pathways_types: List[str],
        clustering_parameters: Dict[str, int | float],
        SAMPLE: str, subfolder: str,
        data_centralized_folder: str, plot_folder: str) -> None:
    """Performs clustering with the selected method parameters"""
    logger.info("computing clustering of spots with method: %s ; "
                "user parameter: %s", clustering_parameters['choice'],
                clustering_parameters['numeric_parameter'])
    # subset to agree with the clean annotations
    adata_mz, clean_annots_dict = safe_load_data_and_clean_annotations(
        SAMPLE, pathways_types, subfolder, data_centralized_folder
    )
    top_mz_by_moran_i = adata_mz.var.loc[
        adata_mz.var['I'] >= 0.1, "mz"].astype(str).tolist()
    datacopy = adata_mz[:, top_mz_by_moran_i].copy()
    sc.pp.log1p(datacopy)
    sc.tl.pca(datacopy, svd_solver="arpack")
    if clustering_parameters['choice'] == "kmeans":
        user_k = int(clustering_parameters['numeric_parameter'])
        kmeans = KMeans(n_clusters=user_k, random_state=42).fit(
            datacopy.obsm['X_pca'])
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
    plot_clustering_results(
        adata_mz, clustering_parameters, SAMPLE, plot_folder,
        colors_spots=options_palette_supported["default"])


def clean_annotations_and_save_tsv(annots_obj: spd.MetabolomeAnnotations,
                                   pathways_types: List[str],
                                   SAMPLE: str, subfolder: str,
                                   data_centralized_folder: str
                                   ) -> None:
    """
    cleans the loaded annotations using core library,
      then saves to tab delimited .tsv files.
      one file for each type of pathway (KEGG or WIKIPATHWAYS) is saved
    """
    adata_mz = safe_load_processed_data(SAMPLE, subfolder,
                                        data_centralized_folder)
    annots_obj, adata_mz = measures_methods.match_adata_mz_and_annots(
        annots_obj, adata_mz)

    # object required to translate annotations to KEGG ids when required
    ref_KEGG_ids_obj = spd.MetaboliteIdsMappings(
        data_dir=data_centralized_folder,
        subfolder="metabolites-mappings",
        mapping_files={'HMDB': "HMDB-to-KEGG_05-07-24.tsv",
                       'LipidMaps': "LipidMaps-to-KEGG_10-18-24.tsv",
                       'CoreMetabolome': "CoreMetabolome-to-KEGG_05-07-24.tsv"})
    ref_KEGG_ids_obj.load()

    # for each type of pathway, save the respective clean annotations
    for pw_type in pathways_types:
        annots_current = annots_obj.model_copy()  # copy is compulsory
        # cleaning annotations using core library:
        if pw_type == "WIKI":
            annots_current.usable_id = "moleculeIds"
            annots_current.verify()
        elif pw_type == "KEGG":
            annots_current.usable_id = "KEGG"
            annots_current.convert_ids_to_kegg(ref_KEGG_ids_obj.table)

        annots_current.keep_best_mz_id_pairs(adata_mz.var.copy())
        annots_current.remove_redundant_aliases()

        # saving to tsv file
        out_file = os.path.join(data_centralized_folder, subfolder,
                                f'{SAMPLE}-{pw_type}-CleanAnnotations.tsv')
        annots_current.df.to_csv(out_file, sep='\t', header=True, index=False)
        del annots_current
    # end for
    logger.info("saved annotations tables after processing")


def normalize_classify_and_save_anndata_mz(
        annots_obj: spd.MetabolomeAnnotations, offSample_cutoff: float,
        msm_cutoff: float, fdr_cutoff: float, moran_I_cutoff: float,
        moran_padj_cutoff: float, SAMPLE: str, subfolder: str,
        data_centralized_folder: str, plot_folder: str) -> None:
    """
    performs the first step consisting of intensities' normalisation
    and the classification of the ions (m/z) by their quality metrics,
    using the core library
    """
    adata_mz = ad.read_h5ad(filename=os.path.join(
        data_centralized_folder, subfolder, f'{SAMPLE}.h5ad'))

    annots_obj, adata_mz = measures_methods.match_adata_mz_and_annots(
        annots_obj, adata_mz)
    adata_mz = spm.NormalizeIntensities.run(adata_mz)

    updated_config = spm.ClassifyMZProtocolConfig(
        msm_cutoff=msm_cutoff, offSample_cutoff=offSample_cutoff,
        fdr_cutoff=fdr_cutoff, moran_I_cutoff=moran_I_cutoff,
        moran_padj_cutoff=moran_padj_cutoff,
    )
    adata_mz = spm.ClassifyMZProtocol.run(
        adata_mz, annots_obj, SAMPLE, plot_folder, config=updated_config)

    adata_mz.write_h5ad(filename=os.path.join(
        data_centralized_folder, subfolder, f"{SAMPLE}.h5ad"))
