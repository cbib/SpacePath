import os
import numpy as np
import anndata as ad
import squidpy as sq
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
from sklearn.utils import sparsefuncs



def verification(adata_mz, me_annots):
    print("verifying mz match")
    try:
        bool_equal = np.all(np.array(
            sorted(me_annots.unfiltered["mz"].unique())) == np.array(
            sorted(adata_mz.var["mz"].tolist())))
        print(f"perfect match between mz in me_a.df "
              f"and adata.X: => {bool_equal} <=")
    except ValueError:
        print("mz do not match between adata_mz and annotation table")
        print("unique annotated mz: ", len(me_annots.unfiltered["mz"].unique(
        )))
        print("unique mz in adata_mz: ", len(adata_mz.var["mz"].unique()))


def fix_annotation_to_adata_mz(gathered_annots, adata_mz):
    # fix the annotations to match with the adata_mz mz features
    me_annots = gathered_annots.model_copy()  # copy is deprecated
    okmz = adata_mz.var['mz']
    me_annots.preloaded = me_annots.preloaded.loc[
        me_annots.preloaded['mz'].isin(okmz)]
    me_annots.unfiltered = me_annots.unfiltered.loc[
        me_annots.unfiltered['mz'].isin(okmz)]

    return me_annots


def get_rid_of_poor_var_mz(adata_mz):
    print("Excluding mz features with extremely low variance")
    means_vars, variance_computed = sparsefuncs.mean_variance_axis(
        adata_mz.X, axis=0)
    #print(pd.DataFrame(variance_computed).describe())
    min_accept_tile = np.percentile(variance_computed, 5)
    #low_variance = variance_computed[variance_computed >= min_accept_centile]
    adata_mz.var["computed_var"] = variance_computed
    good_mz = adata_mz.var.loc[adata_mz.var["computed_var"
                     ] > min_accept_tile, "mz"].to_numpy()
    adata_mz = adata_mz[:, [str(i) for i in good_mz]]

    return adata_mz


def spatial_standard_computations(adata, pca_ncomps=20,
                                  plot_folder=None, prefix='sample'):
    #Xorigin = adata_mz.X.copy()
    # compute metrics
    print("computing variable metrics with squidy")
    try:
        sc.pp.highly_variable_genes(adata,
                                    flavor="seurat_v3",
                                    check_values=False)
    except ValueError:
        pass  # TODO :: add by hand variances and means ?
    print("computing pca")
    sc.pp.pca(adata, n_comps=pca_ncomps)

    kmeans = KMeans(n_clusters=pca_ncomps, n_init=10, random_state=0).fit(
        adata.obsm['X_pca'])
    adata.obs[f'kmeans{pca_ncomps}'] = kmeans.labels_.astype(str)
    if plot_folder is not None:
        tmpstr = os.path.join(plot_folder, f'{prefix}--Kmeans-{pca_ncomps}.png')
        sq.pl.spatial_scatter(adata, color=[f'kmeans{pca_ncomps}'],
                              img=False, title=f"Kmeans (k={pca_ncomps})",
                              save=os.path.abspath(tmpstr))

    return adata


def compute_moran_i(adata_mz, n_neighs=10):
    print("computing spatial neighbors")
    print('with n_neighs: ', n_neighs)
    sq.gr.spatial_neighbors(adata_mz, n_neighs=n_neighs)
    print("Computing Moran's I")
    sq.gr.spatial_autocorr(adata_mz, mode="moran", n_perms=100, n_jobs=1)
    print(adata_mz.uns['moranI'].head())

    return adata_mz


def get_target_mz_pathway(subdf):
    """
    The  subdf represents only *one mz* dataframe (cols: pathway_id, KEGG).
    KEGG is the metabolite identifier.
    For the given distinct KEGG ids and pathway ids in the dataframe,
    first constructs a dict of sets of KEGG ids for each pathway.
    Then, compares each possible pair (k, h) of pathways from the dict,
    to check if there are KEGG ids that are exclusive to k (same for h).
    Returns a dict where the keys are the compared pairs of pathways, e.g.
    {('hsa00230', 'hsa04976'): {'hsa00230': ['C00212', 'C00330'], 'hsa04976': ['C07210']}
    Note: if every (k, h) pathways had same metabolites, output is empty set.
    """
    adict = dict()
    outdict = dict()
    for i, row in subdf.iterrows():
        if row["pathway_id"] not in adict.keys():
            adict[row["pathway_id"]] = {row["KEGG"]}
        adict[row["pathway_id"]].update([row["KEGG"]])
    for k in adict.keys():
        for h in adict.keys():
            if k != h:
                a = adict[k] - adict[h]
                b = adict[h] - adict[k]
                if len(a) > 0 and len(b) > 0:
                    outdict[(k,h)] = {k: list(a), h: list(b)}
    return outdict

def ugly_long_function(tdf) -> None:
    one_kegg_one_pathway = 0
    one_kegg_multiple_pathways = 0
    multiple_kegg_one_pathway = 0
    susceptible_of_refinement = 0
    susceptible_of_refinement_list = []

    for mz in tdf["mz"].unique():
        subdf = tdf.loc[tdf["mz"] == mz, :]
        subdf = subdf.drop(columns=["mz"])
        subdf = subdf.drop_duplicates()
        if subdf.shape[0] == 1:
            one_kegg_one_pathway += 1
        elif (len(subdf["KEGG"].unique()) == 1) and (
                len(subdf["pathway_id"].unique()) > 1):
           one_kegg_multiple_pathways += 1
        elif (len(subdf["KEGG"].unique()) > 1) and (
                len(subdf["pathway_id"].unique()) == 1):
            multiple_kegg_one_pathway += 1
        elif (len(subdf["KEGG"].unique()) > 1) and (
                len(subdf["pathway_id"].unique()) > 1):
            subdf = subdf.drop_duplicates("KEGG")
            out = get_target_mz_pathway(subdf)
            if len(out) > 0:
                susceptible_of_refinement += 1
                susceptible_of_refinement_list.append(mz)
            else:  # not real multi-kegg multi-pathway, see target_mz_pathway
                multiple_kegg_one_pathway += 1

    print("Total mz: ", len(tdf["mz"].unique()))
    ## ok
    print("one_kegg_one_pathway: ", one_kegg_one_pathway)
    print("one_kegg_multiple_pathways: ", one_kegg_multiple_pathways)
    print("multiple_kegg_one_pathway: ", multiple_kegg_one_pathway)
    print("susceptible_of_refinement: ", susceptible_of_refinement)

    return susceptible_of_refinement_list


def assemble_adata_from_scores_result(adatakegg, scores_df):
    adata_scores = ad.AnnData(
            X=csr_matrix(scores_df),
            var=pd.DataFrame(
                {'pathway_id': scores_df.columns.tolist()},
                index=scores_df.columns.tolist()),
            obs=adatakegg.obs,
            obsm=adatakegg.obsm,
            uns=adatakegg.uns
        )
    return adata_scores


def pathways_dict_coherce_min_len(pathways: dict, min_len_admitted: int= 2):
    good_keys = [mykey for mykey, myval in pathways.items()
                     if len(pathways[mykey]) >= min_len_admitted]
    subset_dict = dict(filter(lambda item:item[0] in good_keys, pathways.items()))
    return subset_dict


def modules_scores_by_scanpy(adata, pathways: dict, do_log=False):
    adatacopy = adata.copy()
    if do_log:
        if adatacopy.X.min().min() < 0:
            print("Warning: negative numbers before log")
            print("applying shift of the log values (all to be positives)")
            adatacopy.X = adatacopy.X + abs(adatacopy.X.min().min())
        sc.pp.log1p(adatacopy)

    pathway_names = list()
    scores_lists = list()
    for pw in pathways.keys():
        pathway_names.append(pw)
        mygenelist = pathways[pw]
        miniadata = sc.tl.score_genes(adatacopy, gene_list=mygenelist,
                                      ctrl_size=len(mygenelist), copy=True)
        scores_lists.append(miniadata.obs['score'].tolist())

    df = pd.DataFrame(np.array(scores_lists),
                      index=pathway_names,
                      columns=adatacopy.obs['obs'].tolist(),
                      )
    return df.T