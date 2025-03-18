


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

def exclude_bad_isobars(df, mz_col: str, offsample_col: str,
                         total_intens_col: str, fdr_column: str,
                         msm_column, isobars_column):
    """

    """
    df = df.sort_values(by="mz")
    tmpdf = df.loc[df[isobars_column] != 0,:].copy()
    for

    df = df.sort_values(by="mz")

    return df

###############################

# wrong_mz_manually_defined_dict = {
#     'test635_80-400': [135.029896921, 138.032729825, 140.048379889, 153.043628857, 163.027978793,
#                              167.059278921, 166.038877825, 180.054527889, 202.050966713, 204.054527889,
#                              219.978325273, 259.022442561,
#                              274.033341593, 267.072155657, 277.059672845, 218.992655905,  292.076136075,
#                              299.116793797, 330.125943161 , 341.108935089,
#                              377.085612821, # bad border
#                               351.119774545, 266.070381999,
#                               368.092667239, 374.965451185, 274.033341593,
#                              349.037828277, 333.042977395,   # small ball
#                              295.110645797, 297.126295861, # too much noise
#                               267.081538569
#                              ],   ##  note: this is nice clustered but with artifact :  346.055808353, 300.048991657,  kept for now
#     'test635_200-1K': [] # nothing yet for this other slice
#                              }
###
wrong_mz_manually_defined = [135.029896921, 138.032729825, 140.048379889, 153.043628857, 163.027978793,
                             167.059278921, 166.038877825, 180.054527889, 202.050966713, 204.054527889,
                             219.978325273, 259.022442561,
                             274.033341593, 267.072155657, 277.059672845, 218.992655905,  292.076136075,
                             299.116793797, 330.125943161 , 341.108935089,
                             377.085612821, # bad border
                              351.119774545, 266.070381999,
                              368.092667239, 374.965451185, 274.033341593,
                             349.037828277, 333.042977395,   # small ball
                             295.110645797, 297.126295861, # too much noise
                              267.081538569
                             ]

import plotly.express as px
sns.set_style("whitegrid")
fig = px.scatter(df, x="msm", y="rawOffSampleProb", color='fdr')
fig.show()  # plt.savefig(tmpstr)  does not work !
##################
def haha_from_tic_and_mzqual():
    tmpdf = me_annots.mz_quality.copy()
    tmpdf.index = [str(x) for x in tmpdf['mz'].tolist()]
    tmpdf = tmpdf.loc[adata_mz.var.index.tolist(), :]
    adata_mz.var["fdr"] = tmpdf[me_annots.fdr_column].tolist()
    del tmpdf
    print(adata_mz.var.head())
    ## 3: "bad" by manual definition
    wrong_mz_manually_defined = wrong_mz_manually_defined_dict[SAMPLE]
    wrong_mz_manually_defined = [str(z) for z in wrong_mz_manually_defined]
    good_mz_spatial = list(
        set(good_mz_spatial) - set(wrong_mz_manually_defined))

    print(adata_mz.var.columns)

    adata_mz.var.loc[adata_mz.var['mz'].isin(
        [float(i) for i in good_mz_spatial]), "mz_initial_quality"] = "good"

def failed():
    print("v1::", type(koko['mz']))
    print("q1::", type(me_annots.mz_quality['mz']))
    print("m1::", type(moran_mz_df['mz']))
    print("v::", type(koko['mz'].to_numpy()[10]))
    print("q::", type(me_annots.mz_quality['mz'].to_numpy()[10]))
    print("m::", type(moran_mz_df['mz'].to_numpy()[10]))
    koko2 = koko.join(me_annots.mz_quality[['mz', 'msm',
                                            'rawOffSampleProb', 'fdr',
                                            'quality_msm_outprob']],
                      how="left",
                      on="mz",
                      lsuffix="_x", rsuffix="_y")
    print(koko2.shape)
    koko3 = koko2.join(moran_mz_df[['mz', 'I', 'pval_sim_fdr_bh']],
                       how="left",
                       on="mz",
                       lsuffix="_a", rsuffix="_b")
    print(koko3.shape)
    """nor"""
    foo = me_annots.mz_quality[['mz', 'msm',
                                'rawOffSampleProb', 'fdr',
                                'quality_msm_outprob']].drop_duplicates(
        subset='mz'
    )
    koko = pd.concat([koko, foo
                      ], join="left",
                     keys="mz", ignore_index=True)
    print(koko.shape)
    koko = pd.concat([koko, moran_mz_df[['mz', 'I', 'pval_sim_fdr_bh']]
                      ], join="left",
                     keys="mz", ignore_index=True
                     )
    print(koko.shape)


def ok_yesok():
    print(all_quality_df.shape)
    print(len(set(all_quality_df['mz'].tolist()) - set(
        me_annots.mz_quality['mz'].tolist())))
    print("hahaha", len(set(all_quality_df['mz'].tolist()).intersection(
        set(me_annots.mz_quality['mz'].tolist()))))
    print(len(me_annots.mz_quality['mz'].tolist()))

def candidates_oldies():
    print("mz post sync", adata_mz)
    print("number of pathways:",
          len(pathwaysObj.piled['pathway_id'].unique()))
    # print("plotting log(values)")
    # sc.pp.log1p(adata_mz)
    refinement_mz_pws_dict = find_features_for_refinement(pathwaysObj.piled)

    import squidpy as sq
    j = 0
    for imz in refinement_mz_pws_dict.keys():
        tmpstr = os.path.join(plot_folder, f"Candidate-{j}.png")
        sq.pl.spatial_scatter(adata_mz,
                              color=str(imz), img=False,
                              title=f'{float(imz)}',
                              # (I={imorani}), padj={padj_i}',

                              save=os.path.abspath(tmpstr), dpi=200)
        plt.close()
        j += 1

def simpledisplots():
    sns.displot(x=tmpmo['I'].to_numpy(), bins=30)
    plt.show()
    sns.displot(np.log(tmpmo['var_sim'].to_numpy()), bins=30)
    plt.show()
    sns.displot(np.log(tmpmo['pval_z_sim'].to_numpy()), bins=30)
    plt.show()