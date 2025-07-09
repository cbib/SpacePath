import os
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import scanpy as sc
import squidpy as sq

import zonist.data as zd
import zonist.method as zm
import zonist.processing.measures_methods as zmm
from local_helpers.local_helpers import provide_annotation_dict
from local_helpers.quantifs_methods import differential_results_filter
from local_helpers.quantifs_methods import run_pathway_scores_steps
from local_helpers.local_viz import customized_clustermap

import pickle


def save_clean_annots(adata_mz, me_annots, me_refmap, pwtype) -> None:
    if pwtype == "WIKI":
        me_annots.usable_id = "moleculeIds"
    elif pwtype == "KEGG":
        me_annots.usable_id = "KEGG"

    if pwtype == "WIKI":
        me_annots.verify()
    elif pwtype == "KEGG":
        me_annots.convert_ids_to_kegg(me_refmap.table)

    #me_annots.keep_best_mz_id_pairs(adata_mz.var.copy())  #

    if pwtype == "WIKI":
        me_annots.remove_redundant_aliases()  # when using multitype 'moleculeIds' annotations

    o_fi = os.path.join(mydata_dir,
                        subfolder, f'{SAMPLE}-{pwtype}-CleanAnnots.tsv')

    me_annots.df.to_csv(o_fi, sep='\t', header=True, index=False)


def wrapper_clustering(adata_mz, k,
                       moran_cutoff_features=0,   # recommended and set by default
                       col_features="mz") -> pd.DataFrame:
    top_moran_i = adata_mz.var.loc[
        adata_mz.var['I'] >= moran_cutoff_features, col_features].astype(
        str).tolist()
    mydata = adata_mz[:,
             top_moran_i].copy()  # only do the clustering with this selection
    sc.pp.log1p(mydata)
    sc.tl.pca(mydata)
    kmeans = KMeans(n_clusters=k, random_state=0).fit(mydata.obsm['X_pca'])
    the_labels_str_l = kmeans.labels_.astype(str)
    adata_mz.obs[f'kmeans{k}'] = the_labels_str_l
    adata_mz.obsm['X_pca'] = mydata.obsm['X_pca']
    del mydata

    return adata_mz


def wrapper_differential(adata_mz, k):
    sc.pp.log1p(adata_mz)
    print(
        "computing the differential abundances of all the m/z features, cluster to cluster, for all clusters")
    sc.tl.rank_genes_groups(adata_mz, groupby=f'kmeans{k}', groups='all',
                            reference='rest',
                            method='wilcoxon',  # use logreg if ugly
                            key_added="top_feats_by_cluster",
                            corr_method='benjamini-hochberg')

    differential_results_by_cluster = sc.get.rank_genes_groups_df(
        adata_mz, key="top_feats_by_cluster",
        group=None)  # returns all groups

    return differential_results_by_cluster


def open_clean_annots(SAMPLE, subfolder, pwtypes=["KEGG", "WIKI"]):
    annots_dict = {}
    for pwtype in pwtypes:
        tmpfi = os.path.join(mydata_dir,
                         subfolder, f'{SAMPLE}-{pwtype}-CleanAnnots.tsv')
        annots_dict[pwtype] = pd.read_csv(tmpfi,
                             sep='\t', header=0, index_col=None)
    return annots_dict


def prepare_mols_df(top_differential, clean_ann_df):
    top_differential['names'] = top_differential['names'].astype(float).tolist()
    molecules_ids_df = pd.merge(top_differential, clean_ann_df,
                                how='left', left_on='names', right_on='mz')
    molecules_ids_df.dropna(how='any', inplace=True)
    molecules_ids_df = molecules_ids_df[
        ['moleculeIds', 'mz', 'moleculeNames']].copy()

    molecules_ids_df.drop_duplicates(inplace=True)

    return molecules_ids_df


def booboo(the_ids, molecules_ids_df):
    molecules_ids_df = molecules_ids_df.loc[
        molecules_ids_df['moleculeIds'].isin(the_ids)
    ]
    print("representing a total of mzs:", len(molecules_ids_df['mz'].unique()))


def plot_families(molecules_ids_df, mol_similarities_result,
                  plot_folder, SAMPLE):
    tmp = molecules_ids_df.loc[molecules_ids_df['moleculeIds'].isin(
        mol_similarities_result['the_ids']
    ), ["moleculeIds", "moleculeNames"]].copy()
    tmp.index = tmp['moleculeIds']
    tmp = tmp.T[mol_similarities_result['the_ids']]

    the_molecules_names = tmp.T['moleculeNames'].tolist()

    trunc_names = [i[:30] for i in the_molecules_names]

    arr_df = pd.DataFrame(data=mol_similarities_result['array'],
                          columns=trunc_names,
                          index=trunc_names)
    gc = customized_clustermap(arr_df, computed_dendrogram)
    plt.savefig(
        os.path.join(plot_folder, f"{SAMPLE}-hierarchy-molecules.png"),
        dpi=400)
    plt.close()





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
        max_set_size=sets_args_dict[pwtype]['max_set_size'],
        mode="mz_aware")

    ofipw = os.path.join(mydata_dir, subfolder,
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
#


if __name__ == "__main__":
    """
    bad or good quality of m/z features
    """
    # define here my files and dirs
    mydata_dir = "../../data"
    my_csv_table_path = "../small-test-parameters.tsv"
    plot_parent_folder = "../../figures/"

    force_for_these = []   # fill if need of re-computing Moran's I  # e.g. "LP-example0"

    # start

    datasets_parameters_df = pd.read_csv(my_csv_table_path, header=0,
                      index_col=None, sep='\t')
    datasets_parameters_df = datasets_parameters_df.loc[
        datasets_parameters_df['dataset_name'] == "LP-example0"
    ]

    me_refmap = zd.MetaboliteIdsMappings(
        data_dir=mydata_dir,
        subfolder="metabolites-mappings",
        mapping_files={'HMDB': "HMDB-to-KEGG_05-07-24.tsv",
                       'LipidMaps': "LipidMaps-to-KEGG_10-18-24.tsv",
                       'CoreMetabolome': "CoreMetabolome-to-KEGG_05-07-24.tsv"})
    me_refmap.load()

    print("\nfull pipeline in minimal example data")


    for i, row in datasets_parameters_df.iterrows():
        SAMPLE = row["dataset_name"]
        offSample_cutoff = row["offSample_cutoff"]
        msm_cutoff = row["msm_cutoff"]
        fdr_cutoff = row["fdr_cutoff"]
        moran_I_cutoff = row["Moran_I_cutoff"]
        moran_padj_cutoff = row["Moran_padj_cutoff"]
        k = row['k_for_kmeans']
        cutoff_tree = row['cut_tree_mol_families']

        print("\nSample: ", SAMPLE)

        subfolder = SAMPLE

        plot_folder = os.path.join(plot_parent_folder, SAMPLE)
        if not os.path.exists(plot_folder):
            os.mkdir(plot_folder)

        annotations_to_use = row["annotations_to_use"].split(" ")
        anno_dict = provide_annotation_dict(SAMPLE, annotations_to_use)
        me_annots = zd.MetabolomeAnnotations(
            name=SAMPLE, data_dir=mydata_dir,
            subfolder=subfolder, annotation_files=anno_dict)
        me_annots.load()  # load annotations

        adata_mz = ad.read_h5ad(filename=os.path.join(
            mydata_dir, subfolder, f'{SAMPLE}.h5ad'))
        print("Loaded object", adata_mz.uns["history"])


        ##  1. classify mzs as good or bad quality

        if (SAMPLE not in force_for_these) and ('features_n_after_mz_clean_2' in list(adata_mz.uns['history'].keys())):
            print("Nothing to do for this sample.")
        else:
            try:
                columns_to_reset = set(adata_mz.var.columns.tolist()) - set(["mz"])
                print(columns_to_reset)
                adata_mz.var = adata_mz.var.drop(
                    columns=list(columns_to_reset))
            except KeyError:
                pass

            me_annots, adata_mz = zmm.match_adata_mz_and_annots(me_annots,
                                                            adata_mz)

            adata_mz = zm.NormalizeIntensities.run(adata_mz)

            adata_mz = zm.ClassifyMZProtocol.run(adata_mz,
                me_annots, msm_cutoff, offSample_cutoff, fdr_cutoff,
                moran_I_cutoff, moran_padj_cutoff,
                SAMPLE, plot_folder
            )
            print(adata_mz)
            print(adata_mz.uns['history'])

            adata_mz.write_h5ad(filename=os.path.join(
                mydata_dir, subfolder, f"{SAMPLE}.h5ad"))

        ## 2. Clean annotations that will be used:

        # for each patwhay type: kegg or wiki
        print(adata_mz.uns["history"])
        print(
            "MANDATORY: take away bad quality m/z")  # todo: how to make an alert ?
        adata_mz = adata_mz[:,
                   adata_mz.var.loc[adata_mz.var['final_quality'] == "good",
                   "mz"].astype(str).tolist()]

        me_annots, adata_mz = zmm.match_adata_mz_and_annots(me_annots,
                                                            adata_mz)

        for pwtype in ["KEGG", "WIKI"]:
            save_clean_annots(adata_mz, me_annots, me_refmap, pwtype)
        del me_annots # will no longer be used because clean ones are in tsv now

        ##  3. Clustering step with the good quality mzs

        adata_mz = wrapper_clustering(adata_mz, k=k)  # modified to make moran cutoff as zero :)

        adata_mz.obs[["obs", f"kmeans{k}"]].to_csv(os.path.join(mydata_dir, subfolder,
                             f'{SAMPLE}-clusters-k-{k}.txt'), index=False)

        plot_out_path = os.path.join(plot_folder, f'{SAMPLE}-.png')
        # sq.pl.spatial_scatter(
        #     adata_mz, color=[f'kmeans{k}'], img=False,
        #     title='All m/z (good quality)',
        #     save=os.path.abspath(plot_out_path),
        #     dpi=600)
        # plt.close()


        ##  4. Differential

        sq.pl.spatial_scatter(adata_mz, color=[f'kmeans{k}'],
                              img=False,
                              title="pop up (verification)")
        plt.show()

        ada4diff = adata_mz.copy()
        k_for_diff = 12  # internally defined !

        ada4diff = wrapper_clustering(ada4diff, k=k_for_diff)

        differential_results_by_cluster = wrapper_differential(ada4diff, k_for_diff)

        ## 5. Families of molecules, dark (largest fraction) and no dark included

        annots_dict = open_clean_annots(SAMPLE, subfolder, pwtypes=["KEGG", "WIKI"])

        top_differential = differential_results_filter(differential_results_by_cluster)

        molecules_ids_df = prepare_mols_df(top_differential,
                                           clean_ann_df=annots_dict["WIKI"])

        mol_similarities_result = zm.ComputeMoleculesSimilarityMatrix.run(
            molecules_ids_list=molecules_ids_df["moleculeIds"] )

        with open(os.path.join(mydata_dir, subfolder,
                               f"{SAMPLE}-mol-similarities-matrix.pkl"),"wb") as f:
            pickle.dump(mol_similarities_result, f)

        print(
            "number of compound ids that were used to generate families:",
            len(mol_similarities_result['the_ids']))

        print(booboo(mol_similarities_result['the_ids'], molecules_ids_df))

        print("using cutoff tree :", cutoff_tree)

        computed_dict = zm.ComputeMoleculesFamilies.run(
            mol_similarities_result,
            molecules_ids_df, mydata_dir, subfolder, SAMPLE,
            cutoff_tree=cutoff_tree)

        with open(os.path.join(
                mydata_dir, subfolder, f"{SAMPLE}-mol-fam-dict.pkl"),
                "wb") as f:
            pickle.dump(computed_dict, f)

        computed_dendrogram = computed_dict['computed_dendrogram']

        doplotfam = True
        if doplotfam:
            plot_families(molecules_ids_df, mol_similarities_result,
                          plot_folder, SAMPLE)

        ## scores of pathways and families
        sets_args_dict = {
            "KEGG" : {'file': "KEGG_hsa_pathways_compounds_R110.gmt",
                      'subfolder': "metabolites-mappings",
                      'id_type': 'KEGG',
                      'min_set_size': 2, 'max_set_size':200},
            "WIKI" : {'file': "WIKIPATHWAYS_human_multi_id_redundant.mmt",
                      'subfolder': "metabolites-mappings",
                      'id_type': 'moleculeIds',
                      'min_set_size': 2, 'max_set_size': 200
                      },
            "FAMILIES" : {'file': f"{SAMPLE}-molecules_families.mmt",
                      'subfolder': SAMPLE,
                      'id_type': 'moleculeIds',
                      'min_set_size': 2, 'max_set_size': 200
                      }
        }
        mapped_elems = {}
        for pwtype in ["WIKI","KEGG"]:   # todo activate with the kegg
            mapped_elems[pwtype] = save_scores(adata_mz, annots_dict, pwtype,
                                       sets_args_dict, SAMPLE, mydata_dir)
        # now families *without* the things that mapped in the wiki and kegg
        set_mapped = set(mapped_elems['WIKI']).union(set(mapped_elems["KEGG"]))  # todo activate with the union
        mapped_in_families = save_scores(adata_mz, annots_dict, "FAMILIES",
                    sets_args_dict, SAMPLE, mydata_dir,
                    not_for_families=list(set_mapped))

        ##   Annotation step (new) (what to do???)
        # open the scores objects, take WIKI-scores and families-scores
        # combine them
        # on the selected cluster, do a nested clustering k=2 , make a test
        # plot
        scores_families = ad.read_h5ad(filename=os.path.join(mydata_dir, subfolder,
                         f'{SAMPLE}-scores-FAMILIES.h5ad'))
        chosen_type_pathway = "KEGG" # or WIKI if desired
        scores_pathways = ad.read_h5ad(filename=os.path.join(mydata_dir, subfolder,
                         f'{SAMPLE}-scores-{chosen_type_pathway}.h5ad'))

        scores_families.var['set_type'] = "families"
        scores_pathways.var['set_type'] = "pathways"

        scores_full  = ad.concat((scores_pathways, scores_families),
                                 join='outer', axis=1)

        # # todo : open kegg and families and run:
        # ## for selected clusters, do something with the scores: "nested clster" within the selected cluster
        # my_selected_clusters =  [2]
        # i_clus = 2
        # the_cluster_spots = adata_mz.obs.loc[adata_mz.obs[f'kmeans{k}'] == i_clus,'obs']
        # i_data = adata_mz[the_cluster_spots,:].copy()
        # i_scores = scores_full[the_cluster_spots,:].copy()
        #
        # innerk = 2
        #
        # #  sc.pp.log1p()    sc.tl.pca()   # no log  nor pca to be done as scores are already PCs
        #
        # kmeans = KMeans(n_clusters=innerk, random_state=0).fit(
        #     i_scores.X.todense())
        # the_labels_str_l = kmeans.labels_.astype(str)
        # i_scores.obs[f'kmeans{innerk}'] = the_labels_str_l
        #
        # sq.pl.spatial_scatter(
        #     i_scores, color=[f'kmeans{innerk}'], img=False,
        #     title='???????',
        #     frameon=False,
        #     #save=os.path.abspath(plot_out_path), dpi=600
        #     )
        # #plt.close()
        # plt.show()


        # print to a table the subclustering : say which pathways or families are "on" in which of the two clusters => PCA?