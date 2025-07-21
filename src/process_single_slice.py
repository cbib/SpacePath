from typing import  Tuple
import logging
import warnings
warnings.simplefilter('ignore',category=FutureWarning)

from smpath.data import MetabolomeAnnotations
from smpath.processing.measures_methods import compute_moran_i
from smpath.helpers import setup_logger  # to add extra info to the smpath log
import smpath.data as spd
from sklearn.cluster import KMeans

import pickle

import smpath.processing.measures_methods as measures_methods

from .single_slice_helpers import *
from .single_slice_visualization import *

logger = logging.getLogger(__name__)
setup_logger(logger)


np.random.seed(42)  # important to make the computations be reproducible


pathways_files_location_info = {
        "KEGG" : {'file': "KEGG_hsa_pathways_compounds_R110.gmt",
                      'subfolder': "metabolites-mappings",
                      'id_type': 'KEGG'},
        "WIKI" : {'file': "WIKIPATHWAYS_human_multi_id_redundant.mmt",
                          'subfolder': "metabolites-mappings",
                          'id_type': 'moleculeIds'}}


def join_clean_annotations(pathways_types: List[str],
                    SAMPLE: str, subfolder: str,
                   data_centralized_folder: str)-> pd.DataFrame:
    """
    concatenates annotations for each element in pathways_types
     (but not FAMILIES, as the annotations in this case are the same that served to 'WIKI' type of pathway)
     in other words, concatenates KEGG and no-KEGG (that we coined 'WIKI') annotations
    """
    pathways_types = list(set(pathways_types) - {"FAMILIES"})
    annots_dict = read_clean_annotations(SAMPLE, pathways_types, subfolder, data_centralized_folder)
    dfs_to_gather = list()
    for k in annots_dict.keys():
        if k == 'KEGG':
            df = annots_dict[k][["KEGG", "moleculeNames"]]
            df.columns = ["moleculeIds", "moleculeNames"]
        elif k == 'WIKI':
            df = annots_dict[k][["moleculeIds", "moleculeNames"]]
        dfs_to_gather.append(df)

    gathered_annotations = pd.concat(dfs_to_gather, axis=0)
    return gathered_annotations


def save_top_plots(pathways_types,
                   clustering_parameters, SAMPLE, subfolder,
                   data_centralized_folder, plot_folder)-> None:
    logger.info("saving top pathways and families to plots")
    if not "FAMILIES" in pathways_types:
        pathways_types.append("FAMILIES")
    str_out = f"{clustering_parameters['choice']}-{clustering_parameters['numeric_parameter']}"

    joined = join_scores_adatas(pathways_types, SAMPLE, subfolder, data_centralized_folder)
    map_df = joined.uns['current_pathways_df'].copy()

    gathered_annotations = join_clean_annotations(pathways_types,
                    SAMPLE, subfolder,
                   data_centralized_folder)
    map_df = pd.merge(map_df, gathered_annotations, how='left', left_on='moleculeIds', right_on='moleculeIds')
    diff_sets = pd.read_csv(os.path.join(data_centralized_folder, subfolder,
                                         f'{SAMPLE}-top-differential-SETS.tsv'),
                            sep='\t', header=0)

    clusters_df = pd.read_csv(os.path.join(data_centralized_folder, subfolder,
                                           f'{SAMPLE}-clusters-{str_out}.txt'),
                              header=0)

    # important: add cluster to joined adata
    joined.obs[str_out] = pd.Categorical(clusters_df[str_out].tolist(),
                                         categories=clusters_df[
                                             str_out].unique())

    # mzs will be required in inner loop
    adata_mz = safe_load_processed_data(SAMPLE, subfolder,
                                        data_centralized_folder)

    for cluster in diff_sets['group'].unique().tolist():
        diff_sets_current_cluster = filter_sets_per_cluster_for_plot(diff_sets.copy(),
                                                                     cluster)
        pid = diff_sets_current_cluster['pathway_id'].tolist()
        pnm = diff_sets_current_cluster['pathway_name'].tolist()
        moran_vals = diff_sets_current_cluster['I'].tolist()

        for set_id, set_full_name, mor_v in zip(pid, pnm, moran_vals):
            plot_pathway_path = os.path.join(
                 plot_folder, str(cluster), f'{set_id}-cluster-{cluster}--{SAMPLE}-png')
            try:
                os.makedirs(os.path.join(plot_folder, str(cluster)))
            except:
                pass # skips if already exists

            sq.pl.spatial_scatter(
                joined, color=[set_id], img=False,
                title=f'{set_full_name[:45]}\n{set_id}',
                size=1.2,
                save=os.path.abspath(plot_pathway_path)
            )
            plt.close()

            # plot also the mz belonging to the current pathway
            mzs_in_pw = map_df.loc[
                map_df['pathway_id'] == set_id, 'mz'].unique().tolist()

            if (set_id.startswith('MolFam')) or (set_id.startswith('WP')):
                map_copy = map_df.loc[
                         (map_df['moleculeIds'].str.startswith("HMDB")
                          ) | (map_df['moleculeIds'].str.startswith("LM")), :].copy()
            elif set_id.startswith('hsa'):  # KEGG prefix for H sapiens
                map_copy = map_df.loc[
                    map_df['moleculeIds'].str.startswith("C"), :].copy()
            for i_mz in mzs_in_pw:
                # inner loop
                inner_plot_path = os.path.join(plot_folder,
                                               str(cluster), f'{set_id}-mz{i_mz}--{SAMPLE}.pdf' )
                complex_plot(
                    adata_mz, map_copy, set_id, set_full_name , i_mz,
                    save_path=inner_plot_path
                )


def display_minireport(top_differential, SAMPLE):
    all_l = top_differential.names.unique()
    all_fam = [k for k in all_l if k.startswith("MolFam")]
    wiki_p = [k for k in all_l if k.startswith("WP")]
    kegg_p = [k for k in all_l if k.startswith("hsa")]

    logger.info("SAMPLE: %s", SAMPLE)
    logger.info("   differential KEGG pathways: %s", len(kegg_p))
    logger.info("   differential WIKIPATHWAYS: %s", len(wiki_p))
    logger.info("   differential Families: %s", len(all_fam))


def compute_differential_families_and_pathways(
        pathways_types: List[str], clustering_parameters: Dict[str, str|int|float],
        SAMPLE: str, subfolder: str, data_centralized_folder: str) -> None:
    if not "FAMILIES" in pathways_types:
        pathways_types.append("FAMILIES")
    str_out = f"{clustering_parameters['choice']}-{clustering_parameters['numeric_parameter']}"

    joined = join_scores_adatas(pathways_types, SAMPLE, subfolder,
                                data_centralized_folder)

    clusters_df = pd.read_csv(os.path.join(data_centralized_folder, subfolder,
        f'{SAMPLE}-clusters-{str_out}.txt'), header=0)
    joined.obs[str_out] = pd.Categorical(clusters_df[str_out].tolist(),
                                   categories=clusters_df[
                                       str_out].unique())

    #  shift to avoid negative values
    joined.X = csr_matrix(joined.X.todense() + abs(joined.X.min().min()) )

    differential_results_by_cluster = squidpy_differential_analysis(
        joined, group_str=str_out)

    top_differential = differential_results_filter(
        differential_results_by_cluster, ntop_per_cluster=20)

    # new : moran filter for the sets
    joined = compute_moran_i(joined)
    joined.var['I'] = joined.uns['moranI'].loc[joined.var.index.tolist(),
        'I'].astype(float).tolist()
    top_sets_by_moran_i = joined.uns['moranI'].loc[
                          joined.uns['moranI']['I'] >= 0.2,:].copy().index.tolist()
    top_differential = top_differential.loc[
        top_differential['names'].isin(top_sets_by_moran_i), :]

    top_differential = pd.merge(top_differential,
                                joined.var[['pathway_id', 'pathway_name', 'I']].copy(),
                                how='left',
                                left_on="names", right_on="pathway_id")

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

    scores_adata = run_pathway_scores_steps_produce_anndata(
        metabolites_computed_df, coords, families_obj,
        min_set_size=min_family_size,  type_abbr="FAMILIES"
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

    molecules_ids_df = prepare_molecular_ids_df(top_differential,
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

    #mapped_elements = list()
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
        pathways_obj.molecule_id_type = pathways_files_location_info[pw_type]['id_type']

        pathways_obj.load()
        pathways_obj, annots_df = spm.SyncPathwaysAndAnnots.run(
            pathways_obj.model_copy(),
            annots_df,
            usable_id=pathways_files_location_info[pw_type]['id_type'],
            min_set_size=min_pathway_size)

        metabolites_computed_df, coords = spm.ComputeAssignmentsWeights.run(
            adata_mz, annots_df,
            usable_id=pathways_files_location_info[pw_type]['id_type'])
        del adata_mz  # free some space, we  reload it

        scores_adata = run_pathway_scores_steps_produce_anndata(
            metabolites_computed_df, coords,  pathways_obj,
            min_set_size=min_pathway_size,
             type_abbr=pw_type  )

        out_file = os.path.join(data_centralized_folder, subfolder,
                                f'{SAMPLE}-scores-{pw_type}.h5ad')
        scores_adata.write_h5ad(out_file)
        logging.info("success, saved scores\n")
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
    differential_results_by_cluster = squidpy_differential_analysis(
        adata_mz, group_str=str_out)

    top_differential = differential_results_filter(
        differential_results_by_cluster, ntop_per_cluster=50)

    top_differential['mz'] = top_differential['names'].tolist()
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
        title="Clustering of spots",
        linewidth=0, edgecolors=(0,0,0,0),
        shape="circle",  size=1,#  'circle', 'square', 'hex'
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


def compute_clustering_of_spots(
        pathways_types: List[str],
        clustering_parameters: Dict[str, str|int|float],
        SAMPLE: str, subfolder: str,
        data_centralized_folder: str, plot_folder: str ):
    # subset to agree with the clean annotations
    adata_mz, clean_annots_dict = safe_load_data_and_clean_annotations(
          SAMPLE, pathways_types, subfolder, data_centralized_folder
    )
    top_mz_by_moran_i = adata_mz.var.loc[
         adata_mz.var['I'] >= 0.1, "mz"].astype(str).tolist()
    datacopy = adata_mz[:,top_mz_by_moran_i].copy()
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
    plot_clustering_results(adata_mz, clustering_parameters, SAMPLE, plot_folder)


def clean_annotations_and_save_tsv(annots_obj: MetabolomeAnnotations,
                                   pathways_types: List[str],
                                   SAMPLE: str, subfolder: str, data_centralized_folder: str
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
    for pw_type in pathways_types:
        annots_current = annots_obj.model_copy() # copy is compulsory
        # cleaning annotations using `smpath` library:
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


def normalize_classify_and_save_anndata_mz(
        annots_obj: MetabolomeAnnotations, offSample_cutoff: float,
        msm_cutoff: float, fdr_cutoff: float, moran_I_cutoff: float,
        moran_padj_cutoff: float, SAMPLE: str, subfolder: str,
        data_centralized_folder: str, plot_folder: str) -> None:
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

    updated_config = spm.ClassifyMZProtocolConfig(
        msm_cutoff=msm_cutoff, offSample_cutoff=offSample_cutoff,
        fdr_cutoff=fdr_cutoff, moran_I_cutoff=moran_I_cutoff,
        moran_padj_cutoff=moran_padj_cutoff,
    )
    adata_mz = spm.ClassifyMZProtocol.run(
        adata_mz, annots_obj, SAMPLE, plot_folder, config=updated_config)

    adata_mz.write_h5ad(filename=os.path.join(
        data_centralized_folder, subfolder, f"{SAMPLE}.h5ad"))


# big call

def run_analysis_one_tissue(annotations_dict: Dict[str, str],
                            offSample_cutoff: float ,
                            msm_cutoff: float, fdr_cutoff: float,
                            moran_I_cutoff: float, moran_padj_cutoff: float,
                            clustering_parameters: Dict[str, str|int|float],
                            pathways_types: List[str], min_pathway_size: int,
                            cutoff_tree: float,
                            min_family_size: int,  SAMPLE: str,
                            data_centralized_folder: str,
                            figures_centralized_folder: str)-> None:
    """
    runs step by step the analysis for a single given tissue slice dataset.
    """

    logger.info('\nRunning sample: %s', SAMPLE)

    subfolder = SAMPLE
    plot_folder = os.path.join(figures_centralized_folder, SAMPLE)
    if not os.path.exists(plot_folder):
            os.mkdir(plot_folder)

    annots_obj = spd.MetabolomeAnnotations(
        name=SAMPLE, data_dir=data_centralized_folder,
        subfolder=subfolder, annotation_files=annotations_dict)
    annots_obj.load()  # load annotations with smpath library

    normalize_classify_and_save_anndata_mz(
        annots_obj, offSample_cutoff, msm_cutoff, fdr_cutoff, moran_I_cutoff,
        moran_padj_cutoff, SAMPLE, subfolder,
        data_centralized_folder, plot_folder)

    clean_annotations_and_save_tsv(annots_obj, pathways_types, SAMPLE, subfolder,
                                   data_centralized_folder); del annots_obj  # obj no longer used: use clean .tsv files from now

    compute_clustering_of_spots(pathways_types, clustering_parameters,
        SAMPLE, subfolder, data_centralized_folder, plot_folder)

    compute_differential_analysis_of_mzs(pathways_types,
         clustering_parameters, SAMPLE, subfolder, data_centralized_folder)

    try:
        compute_pathway_scores(pathways_types, min_pathway_size, SAMPLE,
             subfolder, data_centralized_folder)
    except Exception as e:
        logger.warning("failed compute_pathway_scores, %s", e)

    try:
        compute_molecules_families(cutoff_tree,
               SAMPLE, subfolder, data_centralized_folder, plot_folder)
    except Exception as e:
        logger.warning("failed compute_molecules_families , %s", e)

    try:
        compute_families_scores( min_family_size,
            SAMPLE, subfolder, data_centralized_folder)
    except Exception as e:
        logger.warning("failed compute_families_scores , %s", e)

    try:
        compute_differential_families_and_pathways(pathways_types,
            clustering_parameters, SAMPLE, subfolder, data_centralized_folder)
    except Exception as e:
        logger.warning("failed differential families %s", e)

    try:
        save_top_plots(pathways_types,
                       clustering_parameters, SAMPLE, subfolder,
                       data_centralized_folder, plot_folder)
    except Exception as e:
        logger.warning("failed save_top_plots: %s", e)

    logger.info("processed sample %s \n\n", SAMPLE)