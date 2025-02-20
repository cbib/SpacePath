


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


def probabilities_updaterV0(refinement_mz_pws_dict, pathways_scores_df,
                          X_metabolites_probs_df, spot1, spot2):
    # todo: all this is bad because only works with 1 metbolite on each side, and not several
    stealing_fraction = 0.5  # or 1
    the_mz_for_refinement_list = list(refinement_mz_pws_dict.keys())
    updated_probs_spot1 = dict()
    for i_mz in the_mz_for_refinement_list:
        pws_tuples_l = list(refinement_mz_pws_dict[str(i_mz)].keys())
        x = 0
        while (x < len(pws_tuples_l)):
            pws_tuple = pws_tuples_l[x]
            pw1 = pws_tuple[0]
            pw2 = pws_tuple[1]
            # todo: stop taking only one metabolite for each side, because this is not handling cases where more than one metabolite is present
            metabolite1 = refinement_mz_pws_dict[str(i_mz)][pws_tuple][pw1][
                0]
            metabolite2 = refinement_mz_pws_dict[str(i_mz)][pws_tuple][pw2][
                0]
            pw_score = pathways_scores_df.loc[[spot1, spot2], [pw1, pw2]]
            met_prob = X_metabolites_probs_df.loc[
                [spot1, spot2], [metabolite1, metabolite2]]

            if (met_prob.loc[spot1, metabolite1] + met_prob.loc[
                spot1, metabolite2]) == 0:  # impossible to proceed, all zero
                continue
            if (met_prob.loc[spot1, metabolite1] + met_prob.loc[
                spot1, metabolite2]) > 0:
                if  pw_score.loc[spot2, pw1] >= pw_score.loc[spot1, pw1]:
                    if pw_score.loc[spot1, pw1] > pw_score.loc[spot1, pw2]:
                        steal_prop_of_mt2 = met_prob.loc[
                spot1, metabolite2] * stealing_fraction   #   * 0.5 or 1
                        new_pr_1 = met_prob.loc[spot1, metabolite1] + steal_prop_of_mt2
                        new_pr_2 = met_prob.loc[spot1, metabolite2] - steal_prop_of_mt2

                        updated_probs_spot1[metabolite1] = new_pr_1
                        updated_probs_spot1[metabolite2] = new_pr_2
            x += 1
        # end while

    return updated_probs_spot1




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

def trunc_custom(values, decs=0):
    # never used (yet)
    # credits https://stackoverflow.com/questions/42021972/truncating-decimal-digits-numpy-array-of-floats
    return np.trunc(values*10**decs)/(10**decs)

