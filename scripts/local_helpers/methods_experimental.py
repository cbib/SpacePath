



def inner_fun_probs_updater(pw_score: pd.DataFrame, spot1: str, spot2: str, pw1: str, pw2: str,
                            metabolites1: List[str], metabolites2: List[str],
                            met_prob: pd.DataFrame,
                            stealing_fraction: float,
                            evaluated_keggs_curr_spot_curr_mz: List[str],
                            updated_probs_spot1: Dict[str,float]):
    """
    :param pw_score: e.g. :
             hsa00100  hsa00120
        x44_y0  0.576450 -0.056398
        x43_y0 -0.287701  0.670142

    :param spot1:
    :param spot2:
    :param pw1:
    :param pw2:
    :param metabolites1:
    :param metabolites2:
    :param met_prob: e.g.:
                    C01673  C17339  C17331  C17336  C17332
        x44_y0     0.2     0.2     0.2     0.2     0.2
        x43_y0     0.0     0.0     0.0     0.0     0.0


    :param stealing_fraction:
    :param evaluated_keggs_curr_spot_curr_mz:
    :param updated_probs_spot1:
    :return:

    note: as metabolites1 and metabolits2 are often of different lenghts,
    the stealing is done in the pooled values

    """

    if pw_score.loc[spot2, pw1] >= pw_score.loc[spot1, pw1]:
        if pw_score.loc[spot1, pw1] > pw_score.loc[spot1, pw2]:
            steal_prop_of_mt2 = met_prob.loc[
                                    spot1, metabolites2].to_numpy() * stealing_fraction  # * 0.5 or 1
            steal_pool = steal_prop_of_mt2.sum() / len(met_prob.loc[spot1, metabolites1].to_numpy())
            new_pr_2s = met_prob.loc[spot1, metabolites2].to_numpy() - steal_prop_of_mt2
            new_pr_1s = met_prob.loc[spot1, metabolites1].to_numpy() + steal_pool

            # print("the old was:", met_prob.loc[spot1, metabolites2].to_numpy())
            for a in range(len(metabolites1)):
                if metabolites1[a] not in evaluated_keggs_curr_spot_curr_mz:
                    updated_probs_spot1[metabolites1[a]] = new_pr_1s[a]
                    evaluated_keggs_curr_spot_curr_mz.append(metabolites1[a])
            for b in range(len(metabolites2)):
                if metabolites2[b] not in evaluated_keggs_curr_spot_curr_mz:
                    updated_probs_spot1[metabolites2[b]] = new_pr_2s[b]
                    evaluated_keggs_curr_spot_curr_mz.append(metabolites2[b])
            # print("the new group1 is: ", new_pr_1s, metabolites1)
            # print("the new group2 is: ", new_pr_2s, metabolites2)
            # print("")

    return updated_probs_spot1, evaluated_keggs_curr_spot_curr_mz


def probabilities_updater(refinement_mz_pws_dict, pathways_scores_df,
                          X_metabolites_probs_df, spot1, spot2, stealing_fraction: float):
    the_mz_for_refinement_list = list(refinement_mz_pws_dict.keys())
    updated_probs_spot1 = dict()
    for i_mz in the_mz_for_refinement_list:
        pws_tuples_l = list(refinement_mz_pws_dict[str(i_mz)].keys())
        x = 0
        evaluated_keggs_curr_spot_curr_mz = list()
        while (x < len(pws_tuples_l)):
            pws_tuple = pws_tuples_l[x]
            pw1 = pws_tuple[0]
            pw2 = pws_tuple[1]
            metabolites1 = refinement_mz_pws_dict[str(i_mz)][pws_tuple][pw1]
            metabolites2 = refinement_mz_pws_dict[str(i_mz)][pws_tuple][pw2]
            pw_score = pathways_scores_df.loc[[spot1, spot2], [pw1, pw2]]
            met_prob = X_metabolites_probs_df.loc[
                [spot1, spot2], metabolites1 + metabolites2]

            # in spot2: if pw1 is largely superior to pw2 (at least twice)
            if (pw_score.loc[spot2, pw1] >= (pw_score.loc[spot2, pw2] * 2)) and (
                    (met_prob.loc[spot1, metabolites1].sum() + met_prob.loc[
                spot1, metabolites2].sum()) > 0):

                updated_probs_spot1, evaluated_keggs_curr_spot_curr_mz = inner_fun_probs_updater(
                    pw_score, spot1, spot2, pw1, pw2,
                    metabolites1, metabolites2, met_prob, stealing_fraction,
                    evaluated_keggs_curr_spot_curr_mz, updated_probs_spot1)

            x += 1
        # end while

    return updated_probs_spot1


def get_target_pathways_of_one_mz(my_mz: str, subdf, pathways_dict, pws_mz_dict):
    """
    :param: my_mz: current mz
    :param subdf: represents only *one mz* dataframe (cols: pathway_id, KEGG),
       being KEGG the metabolite identifier.
    :param: pathways_dict: keys are pathway ids, values are all the kegg ids
    :param: pws_mz_dict: keys are pathway ids, values are all the mzs
    :return: a dict where the keys are the compared pairs of pathways, e.g.
       {('hsa0230', 'hsa4976'): {'hsa0230': ['C00212'], 'hsa4976': ['C07210']}

    Explanation:

    For the given distinct KEGG ids and pathway ids in the dataframe,
    first constructs a dict (of the current mz) of sets of KEGG ids for each pathway.
    Then, compares each possible pair (k, h) of pathways from the 'mz' dict,
    to check if there are KEGG ids that are exclusive to k (same for h).
    Also takes into account the total KEGG ids of k and h pathways (not ony the current mz).

    If (k, h) pathways had an empty intersection related to the kegg of the current mz,
               !?                 if there are, in each pathway, a total of kegg ids that is superior or equal to the
                 ?               n+1 kegg ids of the current mz,
    and if mzs other than current mz are in each pathway,
    then both (k, h) pathways go into target dict.
    """
    adict_curr_mz: Dict[str, set] = dict()
    pathways_to_target = dict()
    for i, row in subdf.iterrows():
        if row["pathway_id"] not in adict_curr_mz.keys():
            adict_curr_mz[row["pathway_id"]] = {row["KEGG"]}
        adict_curr_mz[row["pathway_id"]].update([row["KEGG"]])
    for k in adict_curr_mz.keys():  # k is a pathway id
        for h in adict_curr_mz.keys():  # h is also a pathway id
            if k != h:
                other_mz_pw1 = set(pws_mz_dict[k]) - set([my_mz])
                other_mz_pw2 = set(pws_mz_dict[h]) - set([my_mz])
                a: set = adict_curr_mz[k] - adict_curr_mz[h]  # a is a set of kegg ids
                b: set = adict_curr_mz[h] - adict_curr_mz[k]  # b is also a set of kegg ids
                intersect_k_h :set = adict_curr_mz[h].intersection(adict_curr_mz[k])
                if (len(a) >= 1) and (len(b) >= 1) and (len(intersect_k_h) == 0):
                    if ((len(other_mz_pw1) >= 1) and (len(other_mz_pw2) >= 2) or (
                            len(other_mz_pw1) >= 2) and (len(other_mz_pw2) >= 1)
                    )and (len(other_mz_pw1.intersection(other_mz_pw2)) == 0):
                        pathways_to_target[(k,h)] = {k: list(a), h: list(b)}

    return pathways_to_target



def message_factory_1(total_nb, kegg_oneormore__one_pathway,
                       one_kegg_multiple_pathways,  multiple_kegg_shared,
                      susceptible_of_refinement) -> None:
    print("Total mz: ", total_nb)
    print("kegg_oneormore__one_pathway: ", kegg_oneormore__one_pathway)
    print("one_kegg_multiple_pathways: ", one_kegg_multiple_pathways)
    print("multiple_kegg_shared: ", multiple_kegg_shared)
    # print("susceptible_of_refinement: ", susceptible_of_refinement)


def find_features_for_refinement(mz_met_pathway_piled: pd.DataFrame):
    """
    :param mz_met_pathway_piled: The df with the correspondences of mz,
        metabolites KEGG ids, and pathways KEGG ids.
    :return: a dictionary with the mz and pairs of pathways that
        are candidates for metabolites assignments refinement, e.g.
        '555.5555999' : {('hsa0230', 'hsa4976'): {
                          'hsa0230': ['C00212'], 'hsa4976': ['C07210']}
    """
    pathways_dict = zh.pws_piled_to_dict(mz_met_pathway_piled,
                                         'pathway_id', 'KEGG')
    pws_mz_dict = zh.pws_piled_to_dict(mz_met_pathway_piled,
                                         'pathway_id', 'mz')
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
            pws_ok = get_target_pathways_of_one_mz(str(mz), subdf, pathways_dict, pws_mz_dict)
            if len(pws_ok) > 0:
                susceptible_of_refinement += 1
                final_dict[str(mz)] = pws_ok
            else:  # not real multi-kegg multi-pathway, see target_mz_pathway
                multiple_kegg_shared += 1

    message_factory_1(len(mz_met_pathway_piled["mz"].unique()), kegg_oneormore__one_pathway,
                      one_kegg_multiple_pathways, multiple_kegg_shared,
                      susceptible_of_refinement)

    return final_dict


def closest_match_mz(mz_list1, mz_list2, tolerance=0.002):
    """ mz_list1 are the mzs as they appear in annotations
        mz_list2 are the mzs as they appear in intensities matrix"""
    df = pd.DataFrame(
        {"v1_f": np.array(mz_list1, dtype=float),
         "v2_f": [np.nan for i in mz_list1],
         "v1_str": mz_list1,
         "v2_str": [np.nan for i in mz_list1]}
    )
    v2_floats = np.array(mz_list2, dtype=float)
    for i, row in df.iterrows():
        tmp = abs(v2_floats - row["v1_f"])
        which_min = np.argmin(tmp)
        if tmp[which_min] <= tolerance:
            df.iloc[i, 1] = v2_floats[which_min]
            df.iloc[i, 3] = str(v2_floats[which_min])
    print("pre", df.shape)
    df = df.dropna()
    print("after", df.shape)
    return df


def spatial_standard_computations(
        adata, k: int=8,  pca_ncomps: int=12, as_log: bool=True,
        top_variable_features_fraction: float=0.5):
    if (adata.var.shape[0] > pca_ncomps) and (adata.var.shape[0] > 20):
        if as_log:
            adata.layers['X_raw'] = adata.X.copy()
            sc.pp.log1p(adata)
        # end if as_log
        topn = int(adata.var.shape[0] * top_variable_features_fraction)
        try:
            sc.pp.highly_variable_genes(
                adata, flavor="seurat_v3",
                check_values=False, n_top_genes=topn)
        except ValueError:
            print("Could not compute highly variable feats with scanpy. Computing with sparsefuncs")
            try:
                means_vars, variance_computed = sparsefuncs.mean_variance_axis(
                    adata.X, axis=0)
                adata.var['variance'] = variance_computed
                adata.var['tmp_index'] = range(adata.var.shape[0])
                tmpdf = adata.var.copy()
                tmpdf = tmpdf.sort_values(by='variance', ascending=False)
                tmpdf = tmpdf.head(topn).copy()
                adata.var['highly_variable'] = False
                adata.var.loc[adata.var['tmp_index'].isin(tmpdf['tmp_index'].tolist()), 'highly_variable'] = True
                adata.var = adata.var.drop(columns=['tmp_index'])
            except ValueError:
                print("impossible to compute highly variable features, marked all as True")
                adata.var['highly_variable'] = True
        print("computing pca")
        sc.pp.pca(adata, n_comps=pca_ncomps)  # use_highly_variable is true by default, but .var['highly_variable'] col was True in all values
        print("computing kmeans")
        kmeans = KMeans(n_clusters=k, n_init=10, random_state=0).fit(
            adata.obsm['X_pca'])
    else:
        print("skipping PCA: too few features")
        X = np.array(adata.X.todense())
        if as_log:
            X = np.log1p(X)
        kmeans = KMeans(n_clusters=k, n_init=10, random_state=0
                        ).fit(X)
    adata.obs[f'kmeans{k}'] = kmeans.labels_.astype(str)

    return adata