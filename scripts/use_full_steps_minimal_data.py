import sys
import os
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import squidpy as sq

import zonist.data as zd
import zonist.method as zm
import zonist.processing.measures_methods as zmm
from local_helpers.local_helpers import provide_annotation_dict
from families_1_generate import differential_results_filter
import pickle



def hahaha(adata_mz, me_annots, me_refmap, pwtype) -> None:
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
                       moran_cutoff_features=0.2  # recommended and set by default
                       ) -> pd.DataFrame:
    top_moran_i = adata_mz.var.loc[
        adata_mz.var['I'] >= moran_cutoff_features, "mz"].astype(
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

    trunc_names = [i[:16] for i in the_molecules_names]

    arr_df = pd.DataFrame(data=mol_similarities_result['array'],
                          columns=trunc_names,
                          index=trunc_names)
    try:  # todo: the plot 'clustermap' slows a lot the computations, make it optional or seek if file to avoid rerun
        sns.clustermap(arr_df, figsize=(13, 13),
                       # row_linkage=computed_dendrogram,
                       cmap="Blues",
                       linewidths=0)
        plt.title(f"the clustering of molecules")
        plt.tight_layout()
        plt.savefig(
            os.path.join(plot_folder, f"{SAMPLE}-hierarchy-molecules.png"),
            dpi=600
        )
    except:
        print("skipping clustermap (too heavy)")

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
        normalize_yes_no = row["TIC_normalization"]
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
            if (normalize_yes_no.lower() in ['y', 'yes']) and (
                    not "is_X_TIC_normalized" in adata_mz.uns.keys()):
                adata_mz.X = csr_matrix(zmm.tic_normalization(
                    adata_mz.X.copy(), adata_mz.obs['obs'].tolist(),
                    adata_mz.var['mz'].tolist()))

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


        for pwtype in ["KEGG"]:  # TODO: add wiki latter  :::: , "WIKI"]:
            hahaha(adata_mz, me_annots, me_refmap, pwtype)
        del me_annots # will no longer be used because clean ones are in tsv now


        ##  3. Clustering step with the good quality mzs

        adata_mz = wrapper_clustering(adata_mz, k=k)

        adata_mz.obs[["obs", f"kmeans{k}"]].to_csv(os.path.join(mydata_dir, subfolder,
                             f'{SAMPLE}-clusters-k-{k}.txt'), index=False)
        plot_out_path = os.path.join(plot_folder,
                                     f'{SAMPLE}-.png')
        sq.pl.spatial_scatter(
            adata_mz, color=[f'kmeans{k}'], img=False,
            title='All m/z (good quality)',
            save=os.path.abspath(plot_out_path),
            dpi=600)
        plt.close()


        ##  4. Differential

        sq.pl.spatial_scatter(adata_mz, color=[f'kmeans{k}'],
                              img=False,
                              title="pop up (verification)")
        plt.show()

        differential_results_by_cluster = wrapper_differential(adata_mz, k)

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

        doplotfam = False
        if doplotfam:
            plot_families(molecules_ids_df, mol_similarities_result,
                          plot_folder, SAMPLE)


        ## scores of pathways and families
        for pwtype in ["KEGG", "WIKI", "FAMILIES"]:
            pass
            #do_scores(adata_mz, pwtype=pwtype)


        ## for selected clusters, do something with the scores: "nested clster" within the selected cluster
        my_selected_clusters =  [2]

        # open the scores objects, take WIKI-scores and families-scores
        # combine them
        # on the selected cluster, do a nested clustering k=2 , make a test
        # plot
        # print to a table the subclustering : say which pathways or families are "on" in which of the two clusters => PCA?




















