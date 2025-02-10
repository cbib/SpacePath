import os
from typing import Dict

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
from scipy.sparse import csr_matrix
from sklearn.cluster import KMeans
from sklearn.utils import sparsefuncs


# todo: a new module named assignments_refine in :  zonist.src.zonist.processing
# def xxxx
def build_kegg_resolved_adata(adata_mz, annots_df):
    """
    Given the adata_mz and the annots_df, returns a new adata object
    """
    mz_list = annots_df["mz"].unique()
    intensities_arrs_list = []
    kegg_list = list()
    for one_mz in mz_list:
        local_kegg = annots_df.loc[
            annots_df["mz"] == one_mz, "KEGG"].tolist()
        local_1D_adata_mz = adata_mz[ :, str(one_mz)].copy()
        for a_kegg in local_kegg:
            intensities_arrs_list.append(
                np.array(local_1D_adata_mz.X.todense().T).flatten() / len(local_kegg))
            kegg_list.append(a_kegg)
    adata = ad.AnnData(X=csr_matrix(np.array(intensities_arrs_list).T),
                       obs=adata_mz.obs,
                       var=pd.DataFrame({"ID": kegg_list},
                                        index=kegg_list),
                       obsm=adata_mz.obsm)
    adata.uns['spatial'] = {'my_slice_name': {
        'scalefactors': {'fiducial_diameter_fullres': 100,
                         'spot_diameter_fullres': 1,
                         'tissue_hires_scalef': 0.5,
                         'tissue_lowres_scalef': 0.5}}}

    return adata


def adata_mz_filter_negligible_abundance(adata_mz, spots_non_zero=5):   #TODO change to 10?
    # TODO: evaluate the nb of zeroes in entire matrix, if too few , fix a minval for converting to zero all values <= minval ??
    non_zero_counts = adata_mz.X.getnnz(axis=0)
    mz_df = adata_mz.var.copy()
    mz_df['bool_arr'] = [1 if x >= spots_non_zero
                         else 0 for x in non_zero_counts]
    adata_mz = adata_mz[:,
               mz_df.loc[mz_df['bool_arr'] == 1, "mz"].astype(str).tolist()]
    return adata_mz


def update_annotation_as_adata_mz(annotationsObj, adata_mz):  # zd.MetabolomeAnnotations
    # update the annotations to match with the adata_mz mz features
    me_annots = annotationsObj.model_copy()  # copy is deprecated
    okmz = adata_mz.var['mz']
    me_annots.preloaded = me_annots.preloaded.loc[
        me_annots.preloaded['mz'].isin(okmz)]
    me_annots.unfiltered = me_annots.unfiltered.loc[
        me_annots.unfiltered['mz'].isin(okmz)]

    return me_annots


def verification(adata_mz, me_annots):  # zd.MetabolomeAnnotations
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


def get_rid_of_poor_var_mz(adata_mz):
    print("Excluding mz features with extremely low variance")
    means_vars, variance_computed = sparsefuncs.mean_variance_axis(
        adata_mz.X, axis=0)
    #logger.info("variance info: %s", pd.DataFrame(variance_computed).describe())
    min_accept_tile = np.percentile(variance_computed, 5)
    #low_variance = variance_computed[variance_computed >= min_accept_centile]
    adata_mz.var["computed_var"] = variance_computed
    good_mz = adata_mz.var.loc[adata_mz.var["computed_var"
                     ] > min_accept_tile, "mz"].to_numpy()
    adata_mz = adata_mz[:, [str(i) for i in good_mz]]

    return adata_mz


def compute_mz_binary_df(adata_mz: ad.AnnData):
    # TODO: evaluate the nb of zeroes in entire matrix, if not enough, fix a minval for converting to zero all values <= minval
    X_presence_absence_mask = adata_mz.X.todense() > 0
    X_mz_binary_df = pd.DataFrame(
        np.array(X_presence_absence_mask, dtype=int),
        columns= adata_mz.var['mz'].to_numpy(),
        index=adata_mz.obs.index.tolist())

    return X_mz_binary_df


def compute_metabolites_probs_df(X_mz_binary_df, annots_OKdf):
    preX_probs = []
    col_names_kegg = []
    for a_mz in X_mz_binary_df.columns:
        i_arr = X_mz_binary_df[a_mz].to_numpy()
        the_kegg_here = annots_OKdf.loc[
            annots_OKdf['mz'] == a_mz, "KEGG"].tolist()
        n = len(the_kegg_here)
        for kid in the_kegg_here:
            col_names_kegg.append(kid)
            initial_prob = 1 / n
            # res = i_arr * np.array([initial_prob for v in range(len(i_arr))])
            res = i_arr * np.repeat(initial_prob, len(i_arr))
            preX_probs.append(np.array(res, dtype=np.float32))

    X_metabolites_probs_df = pd.DataFrame(np.array(preX_probs).T,
                                          columns=col_names_kegg,
                                          index=X_mz_binary_df.index)
    return X_metabolites_probs_df


def trunc_custom(values, decs=0):
    # never used (yet)
    # credits https://stackoverflow.com/questions/42021972/truncating-decimal-digits-numpy-array-of-floats
    return np.trunc(values*10**decs)/(10**decs)


def is_it_well_done(X_metabolites_probs_df, adata_mz, adata, annots_OKdf) -> None:
    newly_comp = dict()
    for a_mz in adata_mz.var['mz'].to_numpy():
        the_kegg_here = annots_OKdf.loc[
            annots_OKdf['mz'] == a_mz, "KEGG"].tolist()
        sub_adata = adata_mz[:, str(a_mz)].copy()
        i_mz_arr = np.array(sub_adata.X.todense().T).flatten()
        for a_kg in the_kegg_here:
            curr_prod_arr = \
                X_metabolites_probs_df[a_kg].to_numpy() * i_mz_arr
            newly_comp[a_kg] = np.array(curr_prod_arr, dtype=np.float32)

    newly_comp_matrix = np.apply_over_axes(trunc_custom, np.array(pd.DataFrame(newly_comp)), axes=0)
    tmp = np.apply_over_axes(trunc_custom, np.array(adata.X.todense()), axes=0)
    if np.all(tmp == newly_comp_matrix):
        print("OK! computed metabolite abundances"
              " by two distinct methods are strictly identical")
    else:
        print("computed metabolite abundances by two distinct"
              " methods are not identical, but please verify decimals")


def closeness_by_z_scores(value1: float, value2: float,
                     values_arr: np.array, z_score_threshold: float=1.0):
    """
    Checks if two values are close to each other with respect to their positions in a distribution.
    :param value1: First value
    :param value2: Second value
    :param values_arr: np.array-like of values representing the distribution
    :param z_score_threshold: Z-score difference threshold for closeness
    :return: True if the values are close, False otherwise
    """
    mean = np.mean(values_arr)
    std_dev = np.std(values_arr)
    z_score_1 = (value1 - mean) / std_dev
    z_score_2 = (value2 - mean) / std_dev

    return abs(z_score_1 - z_score_2) <= z_score_threshold


def threshold_of_presence_by_distribution(values_array: np.array,
                                          z_score_threshold: float=0):
    sorted_arr = np.sort(values_array)
    print(sorted_arr)
    z_scores = (sorted_arr - np.mean(sorted_arr)
                    ) / np.std(sorted_arr)
    tp = [u for u, v in zip(sorted_arr, z_scores) if v >= z_score_threshold]
    print(tp)

    return tp[0]


def closeness_by_fraction(value1: float, value2: float, fraction: float=0.8) -> bool:
    """
    Checks if two values are close to each other, by their max/min quotient

    :param value1: First value
    :param value2: Second value
    :param fraction: The x in [0.5, 1] expected for max_value/min_value ratio
    :return: True if the values are close (ratio >= fraction), False otherwise
    """
    assert 0.5 <= fraction <= 1, "'fraction' must be in the interval [0.5, 1]"
    max_value = max(value1, value2)
    min_value = min(value1, value2)

    return (max_value / min_value) >= fraction

# todo: change name of this fun:
def spatial_standard_computations(adata, pca_ncomps=20,
                                  plot_folder=None, prefix='sample'):
    # TODO: do not mix 'processing' and 'visualization', instead, do a separate funciton for figures
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
    if plot_folder is not None:  # TODO:  move to a visualization fun
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


def get_target_pathways_of_one_mz(subdf, pathways_dict):
    """
    :param subdf: represents only *one mz* dataframe (cols: pathway_id, KEGG),
       being KEGG the metabolite identifier.
    :param: pathways_dict: keys are pathway ids, values are all the kegg ids
    :return: a dict where the keys are the compared pairs of pathways, e.g.
       {('hsa0230', 'hsa4976'): {'hsa0230': ['C00212'], 'hsa4976': ['C07210']}

    Explanation:

    For the given distinct KEGG ids and pathway ids in the dataframe,
    first constructs a dict (of the current mz) of sets of KEGG ids for each pathway.
    Then, compares each possible pair (k, h) of pathways from the 'mz' dict,
    to check if there are KEGG ids that are exclusive to k (same for h).
    Also takes into account the total KEGG ids of k and h pathways (not ony the current mz).

    If (k, h) pathways had an empty intersection related to the kegg of the current mz,
    and if there are, in each pathway, a total of kegg ids that is superior or equal to the
    n+1 kegg ids of the current mz, then both (k, h) pathways go into target dict.
    """
    adict_curr_mz: Dict[str, set] = dict()
    pathways_to_target = dict()
    for i, row in subdf.iterrows():
        if row["pathway_id"] not in adict_curr_mz.keys():
            adict_curr_mz[row["pathway_id"]] = {row["KEGG"]}
        adict_curr_mz[row["pathway_id"]].update([row["KEGG"]])
    for k in adict_curr_mz.keys():
        for h in adict_curr_mz.keys():
            if k != h:
                pw1_allkegg = set(pathways_dict[k])  # total k
                pw2_allkegg = set(pathways_dict[h])   # total h
                a: set = adict_curr_mz[k] - adict_curr_mz[h]
                b: set = adict_curr_mz[h] - adict_curr_mz[k]
                intersect_k_h :set = adict_curr_mz[h].intersection(adict_curr_mz[k])

                if (len(a) >= 1 ) and (len(b) >= 1):   ##  (len(a) == 1 ) and (len(b) == 1): #
                    if (len(intersect_k_h) == 0) and (
                            len(pw1_allkegg) >= (len(a) + 1)) and (
                            len(pw2_allkegg) >= (len(b) + 1)):
                        pathways_to_target[(k,h)] = {k: list(a), h: list(b)}

    return pathways_to_target


def message_factory_1(total_nb, kegg_oneormore__one_pathway,
                       one_kegg_multiple_pathways,  multiple_kegg_shared,
                      susceptible_of_refinement) -> None:
    print("Total mz: ", total_nb)
    print("kegg_oneormore__one_pathway: ", kegg_oneormore__one_pathway)
    print("one_kegg_multiple_pathways: ", one_kegg_multiple_pathways)
    print("multiple_kegg_shared: ", multiple_kegg_shared)
    print("susceptible_of_refinement: ", susceptible_of_refinement)


def ugly_long_function(mz_met_pathway_piled: pd.DataFrame, pathways_dict) :
    """
    :param mz_met_pathway_piled: The df with the correspondences of mz,
        metabolites KEGG ids, and pathways KEGG ids.
    :return: a dictionary with the mz and pairs of pathways that
        are candidates for metabolites assignments refinement, e.g.
        '555.5555999' : {('hsa0230', 'hsa4976'): {
                          'hsa0230': ['C00212'], 'hsa4976': ['C07210']}
    """
    kegg_oneormore__one_pathway = 0
    one_kegg_multiple_pathways = 0
    multiple_kegg_shared = 0
    susceptible_of_refinement = 0
    final_dict: Dict[str, Dict[tuple, Dict[str, list]]] = dict()
    for mz in mz_met_pathway_piled["mz"].unique():
        subdf = mz_met_pathway_piled.loc[mz_met_pathway_piled["mz"] == mz, :]
        subdf = subdf.drop(columns=["mz"])
        subdf = subdf.drop_duplicates()
        if subdf.shape[0] == 1:
            kegg_oneormore__one_pathway += 1
        elif (len(subdf["KEGG"].unique()) == 1) and (
                len(subdf["pathway_id"].unique()) > 1):
           one_kegg_multiple_pathways += 1
        elif (len(subdf["KEGG"].unique()) > 1) and (
                len(subdf["pathway_id"].unique()) == 1):
            kegg_oneormore__one_pathway += 1
        elif (len(subdf["KEGG"].unique()) > 1) and (
                len(subdf["pathway_id"].unique()) > 1):
            subdf = subdf.drop_duplicates("KEGG")
            pws_ok = get_target_pathways_of_one_mz(subdf, pathways_dict)
            if len(pws_ok) > 0:
                susceptible_of_refinement += 1
                final_dict[str(mz)] = pws_ok
            else:  # not real multi-kegg multi-pathway, see target_mz_pathway
                multiple_kegg_shared += 1

    message_factory_1(len(mz_met_pathway_piled["mz"].unique()), kegg_oneormore__one_pathway,
                      one_kegg_multiple_pathways, multiple_kegg_shared,
                      susceptible_of_refinement)

    return final_dict


def my_switch(x: int):
    match x:
        case x if x == 0:
            return "unchanged"
        case x if (x >= 1) and (x <= 5):
            return "1to5Kegg"
        case x if (x >= 6) and (x <= 10):
            return "6to10Kegg"
        case x if x > 10:
            return "11andmoreKegg"
        case _:
            return "FAILED"


def evaluate_spot1_change(refinement_mz_pws_dict, Xbinary_pathways,
                             X_metabolites_probs_df, spot1, spot2):

    present_pathway = 1
    absent_pathway = 0
    the_mz_for_refinement_list = list(refinement_mz_pws_dict.keys())
    modifiable_keggs = set()
    for one_mz in the_mz_for_refinement_list:
        pws_tuples_l = list(refinement_mz_pws_dict[str(one_mz)].keys())
        for pws_tuple in pws_tuples_l:
            pw1 = pws_tuple[0];   pw2 = pws_tuple[1]
            metabolites_pw1 = refinement_mz_pws_dict[str(one_mz)][pws_tuple][pw1]
            metabolites_pw2 = refinement_mz_pws_dict[str(one_mz)][pws_tuple][pw2]
            pw1andpw2_metabolites = list(set(
                metabolites_pw1).union(set(metabolites_pw2)))
            pw_score = Xbinary_pathways.loc[[spot1, spot2], [pw1, pw2]]
            metabo_probs_local = X_metabolites_probs_df.loc[
                [spot1, spot2], pw1andpw2_metabolites]
            if (metabo_probs_local.sum().sum()) == 0:  # impossible to proceed, all values are zero
                continue
            else:
                if pw_score.loc[spot1, pw1] == present_pathway:
                    if pw_score.loc[spot2, pw1] == present_pathway:
                        modifiable_keggs = modifiable_keggs.union(set(metabolites_pw1))

    return modifiable_keggs


def complicated_stuff_1(pw_score, met_prob, spot1, spot2, pw1, pw2, metabolite1, metabolite2):
    green_light_pw1 = False
    green_light_pw2 = False
    happy1 = 0
    happy2 = 0
    distr_vals = pw_score.loc[[spot1, spot2], [pw1, pw2]].to_numpy().flatten()
    presence_thr = np.percentile(distr_vals, 50)
    barre_pw1 = 10000 # todo: replace, use a percentile 50 of the local distribution of the pw1.&.pw2 scores or something relative such: min(pw_score.loc[spot1, pw1] / 1.2,  pw_score.loc[spot2, pw1] / 1.2)
    dummy_ammount = 10  # bad threshold TODO fix this  because values have wide range ! 0 - 60K
    if max(pw_score.loc[spot1, pw1],
           pw_score.loc[
               spot2, pw1]) > presence_thr:#dummy_ammount:  # pw1 in at least one spot is strong
        if abs(pw_score.loc[spot1, pw1] - pw_score.loc[
            spot2, pw1]) <= barre_pw1:  # and both have pw1 similarly strong
            green_light_pw1 = True
    # we want to 'steal' weight from pw2, but first we have to see if pw2 is weaker to pw1 to be able to do so, in spot 1
    if pw_score.loc[spot1, pw2] < pw_score.loc[spot1, pw1]:
        green_light_pw2 = True

    if green_light_pw1 or green_light_pw2:
        gain_or_pain = 0.1
        if met_prob.loc[spot1, metabolite2] > 0:
            happy2 = float(0)
            if met_prob.loc[spot1, metabolite2] >= gain_or_pain:
                happy2 = met_prob.loc[spot1, metabolite2] - gain_or_pain
            else:
                gain_or_pain = met_prob.loc[spot1, metabolite2]
            happy1 = met_prob.loc[spot1, metabolite1] + gain_or_pain

    return happy1, happy2


def complicated_stuff_2(refinement_mz_pws_dict, pathways_scores_df, X_metabolites_probs_df, spot1,spot2):
    the_mz_for_refinement_list = list(refinement_mz_pws_dict.keys())
    updated_probs_dict = dict()
    for one_mz in the_mz_for_refinement_list:
        pws_tuples_l = list(refinement_mz_pws_dict[str(one_mz)].keys())
        for pws_tuple in pws_tuples_l:
            pw1 = pws_tuple[0];   pw2 = pws_tuple[1]
            metabolite1 = refinement_mz_pws_dict[str(one_mz)][pws_tuple][pw1][
                0]
            metabolite2 = refinement_mz_pws_dict[str(one_mz)][pws_tuple][pw2][
                0]
            pw_score = pathways_scores_df.loc[[spot1, spot2], [pw1, pw2]]
            met_prob = X_metabolites_probs_df.loc[
                [spot1, spot2], [metabolite1, metabolite2]]
            # green_light_pw1 = False
            # green_light_pw2 = False
            if (met_prob.loc[spot1, metabolite1] + met_prob.loc[
                spot1, metabolite2]) == 0:  # impossible to proceed, both are zero
                pass
            elif (met_prob.loc[spot1, metabolite1] + met_prob.loc[
                spot1, metabolite2] > 0):
                happy1, happy2 = complicated_stuff_1(pw_score, met_prob,
                                                     spot1, spot2, pw1, pw2,
                                                     metabolite1, metabolite2)
                if happy1 + happy2 > 0:
                    updated_probs_dict = {
                            metabolite1: happy1,
                            metabolite2: happy2}
                    break   #  it is ok if able to update with compared pathways , we can stop ??? # Todo check this !! not sure :/

    return updated_probs_dict


def subset_m_probs_curr_pw(metabolites_probs_here, pathwaysObj, pw):  # TODO is this used ?

    theor_compounds_in_pw = pathwaysObj.to_dict()[pw]
    real_compounds_in_pw = set(theor_compounds_in_pw).intersection(
        set(metabolites_probs_here.columns.tolist()))

    return metabolites_probs_here[list(real_compounds_in_pw)].T

if __name__ == "__main__":
    print(my_switch(0))