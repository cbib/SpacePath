import os
import pandas as pd
import numpy as np
import anndata as ad
import zonist.data as zd
import zonist.method as zm
import zonist.processing.measures_methods as zmm
from scipy.sparse import csr_matrix

from local_helpers.local_helpers import provide_annotation_dict
from local_helpers.quantifs_methods import differential_results_filter
from local_helpers.quantifs_methods import run_pathway_scores_steps
from local_helpers.local_viz import customized_clustermap

from sklearn.cluster import KMeans
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
from analysis_helpers import *
"""
All this workflow relies on the library XXXX
@github     # here the github TODO
"""

"""

to use when doing differential of mzs : 
clusters_df = pd.read_csv(
        
"""

np.random.seed(42)  # important to make reproducible the computations

def do_the_next():
    print("ok")
    return 0

def compute_differential_families_and_pathways(SAMPLE,   pathways_types,
        clustering_parameters, data_centralized_folder,subfolder):
    str_out = f"{clustering_parameters['choice']}-{clustering_parameters['numeric_parameter']}"
    clusters_df = pd.read_csv(os.path.join(data_centralized_folder, subfolder,
                                           f'{SAMPLE}-clusters-{str_out}.txt'),
                              header=0)
    pathways_types += ["FAMILIES"]
    list_of_scores_objects =  list()
    for u in pathways_types:
        adata_scores = ad.read_h5ad(filename=os.path.join(data_centralized_folder, subfolder,
                         f'{SAMPLE}-scores-{u}.h5ad'))
        adata_scores.var['set_type'] = u
        list_of_scores_objects.append(adata_scores)
        del adata_scores
    joined = ad.concat(adatas=list_of_scores_objects,
              join='outer', axis=1)
    joined.obs[str_out] = pd.Categorical(clusters_df[str_out].tolist(),
                                   categories=clusters_df[
                                       str_out].unique())
    joined.X = csr_matrix(joined.X.todense() + abs(joined.X.min().min()) ) #  shift to avoid negative values
    print("\n========\n", joined, "\n========\n")
    top_differential = differential_routine(joined, group_str=str_out)
    all_l = top_differential.names.unique()
    all_fam = [k for k in all_l if k.startswith("MolFam")]
    wiki_p = [k for k in all_l if k.startswith("WP")]
    kegg_p = [k for k in all_l if k.startswith("hsa")]


    print("differential_stuff:", SAMPLE  )
    print("   families:", len(all_fam)  )
    print("   kegg:", len(kegg_p)  )
    print("   wikipw:", len(wiki_p)  )






def compute_families_scores(SAMPLE,  min_family_size,
          data_centralized_folder,  subfolder):

    families_obj = zd.MetabolicPathways(
        data_dir=data_centralized_folder,
        subfolder=subfolder,
        file=f"{SAMPLE}-molecules_families.mmt"
    )
    families_obj.molecule_id_type = 'moleculeIds'
    families_obj.load()

    adata_mz = safe_load_processed_data(SAMPLE, subfolder,
                                        data_centralized_folder)

    clean_annots_nokegg_ids = pd.read_csv(
        os.path.join(data_centralized_folder, subfolder,
                     f'{SAMPLE}-WIKI-CleanAnnotations.tsv'),
        sep='\t', header=0, index_col=None)

    # families_obj.intersect_with_annotations(clean_annots_nokegg_ids,
    #                                     usable_id='moleculeIds',
    #                                     mz_col='mz')
    families_obj, annots_df = zm.SyncPathwaysAndAnnots.run(
        families_obj.model_copy(),
        clean_annots_nokegg_ids,
        usable_id='moleculeIds',
        min_set_size=min_family_size)

    scores_adata = run_pathway_scores_steps(adata_mz.copy(),
                                            families_obj,
                                            annots_df,
                                            usable_id='moleculeIds')
    
    out_file = os.path.join(data_centralized_folder, subfolder,
                         f'{SAMPLE}-scores-FAMILIES.h5ad')
    scores_adata.write_h5ad(out_file)


def compute_molecules_families(SAMPLE,  cutoff_tree,
          data_centralized_folder,  subfolder, plot_folder):
    top_differential = pd.read_csv(
        os.path.join(data_centralized_folder, subfolder,
                     f'{SAMPLE}-top-differential-mzs.tsv'), sep='\t', header=0)

    clean_annots_nokegg_ids = pd.read_csv(
        os.path.join(data_centralized_folder, subfolder,
                     f'{SAMPLE}-WIKI-CleanAnnotations.tsv'),
        sep='\t', header=0, index_col=None)

    molecules_ids_df = prepare_mols_df(top_differential,
                                       clean_ann_df=clean_annots_nokegg_ids)

    mol_similarities_result = zm.ComputeMoleculesSimilarityMatrix.run(
        molecules_ids_list=molecules_ids_df["moleculeIds"])

    with open(os.path.join(data_centralized_folder, subfolder,
                           f"{SAMPLE}-mol-similarities-matrix.pkl"),
              "wb") as f:
        pickle.dump(mol_similarities_result, f)

    print("number of compound ids that were used to generate families:",
        len(mol_similarities_result['the_ids']))
    print(message_ids_how_many_mz(mol_similarities_result['the_ids'],
                                  molecules_ids_df))

    computed_dict = zm.ComputeMoleculesFamilies.run(
        mol_similarities_result,
        molecules_ids_df, data_centralized_folder, subfolder, SAMPLE,
        cutoff_tree=cutoff_tree)

    with open(os.path.join(
            data_centralized_folder, subfolder, f"{SAMPLE}-mol-fam-dict.pkl"),
            "wb") as f:
        pickle.dump(computed_dict, f)

    computed_dendrogram = computed_dict['computed_dendrogram']

    plot_families(molecules_ids_df, mol_similarities_result,
                  computed_dendrogram, plot_folder, SAMPLE)


def compute_pathway_scores(SAMPLE, pathways_types, min_pathway_size,
         data_centralized_folder,  subfolder):
    adata_mz, clean_annots_dict = safe_load_data_and_clean_annotations(
        SAMPLE, pathways_types, subfolder, data_centralized_folder
    )
    sets_args_dict = files_of_pathways_defaults()
    #mapped_elements = list()
    for pwtype in pathways_types:
        annots_df = clean_annots_dict[pwtype]
        pathways_obj = zd.MetabolicPathways(
            data_dir=data_centralized_folder,
            subfolder=sets_args_dict[pwtype]['subfolder'],
            file=sets_args_dict[pwtype]['file']
        )
        pathways_obj.molecule_id_type = sets_args_dict[pwtype]['id_type']

        pathways_obj.load()
        pathways_obj, annots_df = zm.SyncPathwaysAndAnnots.run(
            pathways_obj.model_copy(),
            annots_df,
            usable_id=sets_args_dict[pwtype]['id_type'],
            min_set_size=min_pathway_size)

        out_file = os.path.join(data_centralized_folder, subfolder,
                             f'{SAMPLE}-scores-{pwtype}.h5ad')

        scores_adata = run_pathway_scores_steps(adata_mz.copy(), pathways_obj,
            annots_df,  usable_id=sets_args_dict[pwtype]['id_type'])
        print(scores_adata.var.shape)
        scores_adata.write_h5ad(out_file)
        print("youpi, saved scores")
        print("")
        #mapped_elements += scores_adata.uns['current_pathways_df']['mz'].unique(
        #).tolist()


def compute_differential_analysis_of_mzs(SAMPLE, pathways_types,
    clustering_parameters, data_centralized_folder,  subfolder
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
    top_differential = differential_routine(adata_mz, group_str=str_out)
    top_differential.to_csv(os.path.join(data_centralized_folder, subfolder,
            f'{SAMPLE}-top-differential-mzs.tsv'),
                              sep='\t', header=True, index=False)


def plot_clustering_results(SAMPLE, adata_mz, clustering_parameters,
                             plot_folder):
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
    datacopy, resolution_param
):
    sc.pp.neighbors(datacopy)
    nn_graph_genes = datacopy.obsp["connectivities"]
    sq.gr.spatial_neighbors(datacopy)
    nn_graph_space = datacopy.obsp["spatial_connectivities"]

    alpha = 0.25  # the weight of the spatial graph
    joint_graph = (1 - alpha) * nn_graph_genes + alpha * nn_graph_space
    sc.tl.leiden(datacopy,
                 adjacency=joint_graph,
                 n_iterations=2,
                 flavor="igraph",
                 directed=False, key_added="Leiden_sg",
                 resolution=resolution_param)
    the_labels_str_l = datacopy.obs['Leiden_sg'].astype(str).tolist()
    return datacopy, the_labels_str_l


def compute_clustering_of_spots(SAMPLE, pathways_types, clustering_parameters,
       data_centralized_folder,  subfolder, plot_folder ):
    # subset to agree with the clean annotations
    adata_mz, clean_annots_dict = safe_load_data_and_clean_annotations(
          SAMPLE, pathways_types, subfolder, data_centralized_folder
    )
    # top_moran_i = adata_mz.var.loc[
    #     adata_mz.var['I'] >= 0.1, "mz"].astype(str).tolist()
    # datacopy = adata_mz[:,top_moran_i].copy()  # only do the clustering with this selection
    datacopy = adata_mz.copy()
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
    plot_clustering_results(SAMPLE, adata_mz, clustering_parameters,
                             plot_folder)


def clean_annotations_and_save(SAMPLE, annots_obj, pathways_types,
    data_centralized_folder, subfolder
):
    adata_mz = safe_load_processed_data(SAMPLE, subfolder, data_centralized_folder)
    annots_obj, adata_mz = zmm.match_adata_mz_and_annots(annots_obj, adata_mz)

    ref_KEGG_ids_obj = zd.MetaboliteIdsMappings(
        data_dir=data_centralized_folder,
        subfolder="metabolites-mappings",
        mapping_files={'HMDB': "HMDB-to-KEGG_05-07-24.tsv",
                       'LipidMaps': "LipidMaps-to-KEGG_10-18-24.tsv",
                       'CoreMetabolome': "CoreMetabolome-to-KEGG_05-07-24.tsv"})
    ref_KEGG_ids_obj.load()  # required to translate to KEGG ids

    for pwtype in pathways_types:
        annots_current = annots_obj.model_copy() # copy required
        clean_annotations(adata_mz, annots_current, ref_KEGG_ids_obj, pwtype)
        out_file = os.path.join(data_centralized_folder, subfolder,
            f'{SAMPLE}-{pwtype}-CleanAnnotations.tsv')
        annots_current.df.to_csv(out_file, sep='\t', header=True, index=False)
        del annots_current
    # end for


def normalize_classify_and_save(SAMPLE, annots_obj, offSample_cutoff,
    msm_cutoff, fdr_cutoff,  moran_I_cutoff, moran_padj_cutoff,
    data_centralized_folder, subfolder, plot_folder) -> None:
    """
    performs the first step consisting of intensities' normalisation
    and the classification of the ions (m/z) by their quality metrics.
    """

    adata_mz = ad.read_h5ad(filename=os.path.join(
        data_centralized_folder, subfolder, f'{SAMPLE}.h5ad'))

    try:  # resets adata.var columns if user re-runs analysis
        columns_to_reset = set(adata_mz.var.columns.tolist()) - set(["mz"])
        adata_mz.var = adata_mz.var.drop(columns=list(columns_to_reset))
    except KeyError:
        pass

    annots_obj, adata_mz = zmm.match_adata_mz_and_annots(annots_obj, adata_mz)
    adata_mz = zm.NormalizeIntensities.run(adata_mz)
    adata_mz = zm.ClassifyMZProtocol.run(
        adata_mz, annots_obj, msm_cutoff, offSample_cutoff, fdr_cutoff,
        moran_I_cutoff, moran_padj_cutoff, SAMPLE, plot_folder)

    adata_mz.write_h5ad(filename=os.path.join(
        data_centralized_folder, subfolder, f"{SAMPLE}.h5ad"))



# the function that calls all of the above

def run_analysis_one_tissue(SAMPLE, annotations_dict, offSample_cutoff,
                            msm_cutoff, fdr_cutoff,
                            moran_I_cutoff, moran_padj_cutoff,
                            clustering_parameters,
                            pathways_types, min_pathway_size, cutoff_tree,
                            min_family_size,
                            data_centralized_folder, figures_centralized_folder)-> None:
    """runs every step for a given dataset and its parameters"""
    subfolder = SAMPLE
    plot_folder = os.path.join(figures_centralized_folder, SAMPLE)
    if not os.path.exists(plot_folder):
            os.mkdir(plot_folder)

    annots_obj = zd.MetabolomeAnnotations(
        name=SAMPLE, data_dir=data_centralized_folder,
        subfolder=subfolder, annotation_files=annotations_dict)
    annots_obj.load()  # load annotations

    normalize_classify_and_save(SAMPLE, annots_obj, offSample_cutoff,
        msm_cutoff, fdr_cutoff,  moran_I_cutoff, moran_padj_cutoff,
        data_centralized_folder, subfolder, plot_folder)

    clean_annotations_and_save(SAMPLE, annots_obj, pathways_types,
        data_centralized_folder, subfolder)
    del annots_obj   # no longer needed as the clean versions are now as .tsv

    compute_clustering_of_spots(SAMPLE, pathways_types, clustering_parameters,
        data_centralized_folder,  subfolder, plot_folder)

    compute_differential_analysis_of_mzs(SAMPLE, pathways_types,
         clustering_parameters,  data_centralized_folder,  subfolder)

    compute_pathway_scores(SAMPLE, pathways_types, min_pathway_size,
         data_centralized_folder,  subfolder)

    compute_molecules_families(SAMPLE,  cutoff_tree,
           data_centralized_folder,  subfolder, plot_folder)

    compute_families_scores(SAMPLE,  min_family_size,
          data_centralized_folder,  subfolder)

    compute_differential_families_and_pathways(SAMPLE,   pathways_types,
        clustering_parameters, data_centralized_folder,subfolder)

    # todo for documentation: explain the structure of data/ folder,
    #  give the data/ folder with minimal data and with the folder metabolites-mappings !!!!

    # EL FIN !!! 