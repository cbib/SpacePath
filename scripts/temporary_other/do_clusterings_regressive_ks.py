import os
import pandas as pd
import anndata as ad
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import scanpy as sc
import squidpy as sq

import zonist.data as zd
import zonist.method as zm
from local_helpers.local_helpers import provide_annotation_dict
from local_helpers.quantifs_methods import run_pathway_scores_steps
from local_helpers.local_viz import customized_clustermap


def save_clean_annots(adata_mz, me_annots, me_refmap, pwtype) -> None:
    if pwtype == "WIKI":
        me_annots.usable_id = "moleculeIds"
    elif pwtype == "KEGG":
        me_annots.usable_id = "KEGG"

    if pwtype == "WIKI":
        me_annots.verify()
    elif pwtype == "KEGG":
        me_annots.convert_ids_to_kegg(me_refmap.table)

    me_annots.keep_unique_mz_id_tuples(adata_mz.var.copy())  #

    if pwtype == "WIKI":
        me_annots.remove_redundant_aliases()  # when using multitype 'moleculeIds' annotations

    o_fi = os.path.join(mydata_dir,
                        subfolder, f'{SAMPLE}-{pwtype}-CleanAnnots.tsv')

    me_annots.df.to_csv(o_fi, sep='\t', header=True, index=False)


def wrapper_clustering(adata_mz, k,
                       moran_cutoff_features=0,   # recommended and set by default
                       col_features="mz") -> ad.AnnData:
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
    plot_parent_folder = "../../figures/"


    picked_sample = "TD-Lipido-635" # ,  and, LP-gut  #  LP-brain-rat-2   LP-brain-mouse-5

    my_csv_table_path_1 = "../public-datasets-parameters.tsv"  # "../../small_test.csv"
    my_csv_table_path_2 = "../private-datasets-parameters.tsv"  # "../../small_test.csv"

    df1 = pd.read_csv(my_csv_table_path_1, header=0,
                      index_col=None, sep='\t')
    df2 = pd.read_csv(my_csv_table_path_2, header=0,
                      index_col=None, sep='\t')

    datasets_parameters_df = pd.concat([df1, df2])
    force_for_these = []
    # start

    datasets_parameters_df = datasets_parameters_df.loc[
        datasets_parameters_df['dataset_name'] == picked_sample
    ]

    me_refmap = zd.MetaboliteIdsMappings(
        data_dir=mydata_dir,
        subfolder="metabolites-mappings",
        mapping_files={'HMDB': "HMDB-to-KEGG_05-07-24.tsv",
                       'LipidMaps': "LipidMaps-to-KEGG_10-18-24.tsv",
                       'CoreMetabolome': "CoreMetabolome-to-KEGG_05-07-24.tsv"})
    me_refmap.load()

    print("\nregressive kas for the ")


    for i, row in datasets_parameters_df.iterrows():
        SAMPLE = row["dataset_name"]
        normalize_yes_no = row["TIC_normalization"]
        offSample_cutoff = row["offSample_cutoff"]
        msm_cutoff = row["msm_cutoff"]
        fdr_cutoff = row["fdr_cutoff"]
        moran_I_cutoff = row["Moran_I_cutoff"]
        moran_padj_cutoff = row["Moran_padj_cutoff"]
        user_k = row['k_for_kmeans']
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


        ##  1. classify mzs as good or bad quality: already done


        ## 2. Clean annotations that will be used:  already done by the clean annots to tables

        # for each patwhay type: kegg or wiki
        print(adata_mz.uns["history"])
        print(
            "MANDATORY: take away bad quality m/z")  # todo: how to make an alert ?

        adata_mz = adata_mz[:,
                   adata_mz.var.loc[adata_mz.var['final_quality'] == "good",
                   "mz"].astype(str).tolist()]
        annots_dict = open_clean_annots(SAMPLE, subfolder, pwtypes=["KEGG", "WIKI"])
        anns_df = annots_dict['WIKI']
        print(adata_mz)
        #adata_mz = adata_mz[:, anns_df['mz'].astype(str).unique().tolist()].copy() # todo: mark this as important !
        #print(adata_mz)
        ##  3. Clustering step with the good quality mzs
        moran_cut = 0.2
        if SAMPLE == "TD-gut":   # this tissue is special as has tiny true structures sparsed
            moran_cut = 0.0

        for user_k in [8, 7, 6, 5]:  #   , 4, 3
            adata_mz = wrapper_clustering(adata_mz, k=user_k,
                                          moran_cutoff_features=moran_cut)  # modified to make moran cutoff as zero :)
            # save uer diff result, will be used in another script
            adata_mz.obs[["obs", f"kmeans{user_k}"]].to_csv(os.path.join(mydata_dir, subfolder,
                                 f'{SAMPLE}-REGRESSIVES-K--{user_k}.txt'), index=False)

            plot_out_path = os.path.join(plot_folder, f'{SAMPLE}-MZ-kmeans{user_k}.png')
            sq.pl.spatial_scatter(
                adata_mz, color=[f'kmeans{user_k}'], img=False,
                title=f'k={user_k}',
                save=os.path.abspath(plot_out_path),
                dpi=600)
            plt.close()