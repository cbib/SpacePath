

test_pathways_dict = {'p1': [7, 4, 8], 'p2': [8, 5], 'p3': [99]}
res = pathways_dict_coherce_min_len(test_pathways_dict, min_len_admitted=3)   # pathways_dict_coherce_min_len
res


##

fake = ad.AnnData(X=csr_matrix(np.array([[1,0,11,0], [55,88,6,7], [0,2,4,0], [78,0,0,0]], )),
                  var=pd.DataFrame({"mz": [9.1, 8.1, 7.1, 6.1]},
                                   index=[str(x) for x in [9.1, 8.1, 7.1, 6.1]]),
                   obs = pd.DataFrame({"obs": ["x1y1", "x2y2", "x3y3", "x4y4"]}))
fake = adata_mz_filter_negligible_abundance(fake)
print(fake.to_df())
#    9.1  8.1  7.1
# 0    1    0   11
# 1   55   88    6
# 2    0    2    4
# 3   78    0    0


##

U = list(pathwaysObj.to_dict().keys())
V = pathwaysObj.piled["pathway_id"].unique()
empty_a = set(U) - set(V)
empty_b = set(V) - set(U)
if len(empty_a) == 0 and len(empty_b) == 0:
    print("Yeeeeeeeeeeeeeeeeeeeeeeeeah")

##

floatsofme = np.array([8.98354354, 5.31687654])
print(np.trunc(floatsofme, 4))


data = np.array([0,0,0, 3000, 50, 90, 5700, 9700, 25000, 50000])
value1 = 2200
value2 = 3000
estupido = np.percentile(data, 50)
kak = threshold_of_presence_by_distribution(data)

result = closeness_by_interval(value1, value2)
print(f"Are {value1} and {value2} close? {result}")




##
if happy1 + happy2 > 0:
    inner_dict = {
        (one_mz, pws_tuple): {
            metabolite1: (
                met_prob.loc[spot1, metabolite1], happy1),
            metabolite2: (met_prob.loc[spot1, metabolite1],
                          happy2)}}
# ...
jojo = "C00877"
print(len(new_probabs))
print(new_probabs)
increased = 0
decreased = 0
nothingnew = 0
for one_spot in new_probabs.keys():
    for a_mz__pws_tup in new_probabs[one_spot].keys():
        if new_probabs[one_spot][a_mz__pws_tup][jojo][1] == 0.6:
            increased += 1
        elif new_probabs[one_spot][a_mz__pws_tup][jojo][1] == 0.4:
            decreased += 1
        else:
            nothingnew +=1

        # 'C00877'
print("increased ", increased)
print("decreased ", decreased)
print("shit ", nothingnew)


##
better = []
for k in range(4, 8):
    for h in range(3, 7):
        better.append(f"x{k}_y{h}")

print(better)



df1 = pd.DataFrame(data=np.array([[55,88,6,7], [0,2,4,0],
                                 [78,0,0,0], [78,0,0,0]]),
                  columns=["c05", "c04", "c11", "c12"],
                  index=['a', 'u', 'v', 'o'])

df2 = pd.DataFrame(data=np.array([[1,0,0],  [0,2,0], [78,0,0]]),
                  columns=["c10", "c11", "c12"], index=['e', 'i', 'w'])

dfs_dicts = {0: df1, 1: df2}
ki = verticalassembly(dfs_dicts, ['a', 'e', 'i', 'o', 'u', 'v', 'w'])
print(ki)

## ## ##
metabo_sums = np.array(X_metabolites_probs_df).sum(axis=0) # sum of pixels intensities by compound, final dimension will be nb of compounds

##
new_probabs = dict()

spot1 = "x72_y45"
spot2 = "x71_y45"
tmp_mx1 = adata_mz[spot1, :].X.todense()  # already 2D array as required by cosine_similarity
tmp_mx2 = adata_mz[spot2, :].X.todense()
local_cosine_sim = cosine_similarity(np.array(tmp_mx1),
                                             np.array(tmp_mx2))
# new_probabs[spot1] = dict()
#
# #X_metabolites_probs_df
# for one_mz in the_mz_for_refinement_list:
#     pws_tuples_l = list(refinement_mz_pws_dict[str(one_mz)].keys())
#     for pws_tuple in pws_tuples_l:
#         pw1 = pws_tuple[0] ;      pw2 = pws_tuple[1]   # todo: what if I say that I choose the second element as pw1 ? (inverse sense?)
#         metabolite1 = refinement_mz_pws_dict[str(one_mz)][pws_tuple][pw1][0]
#         metabolite2 = refinement_mz_pws_dict[str(one_mz)][pws_tuple][pw2][0]
#         pw_score = pathways_pseudoscores.loc[[spot1, spot2],  [pw1, pw2]]
#         met_prob = X_metabolites_probs_df.loc[[spot1, spot2],  [metabolite1, metabolite2]]
#         # green_light_pw1 = False
#         # green_light_pw2 = False
#
#         if (met_prob.loc[spot1, metabolite1] + met_prob.loc[
#             spot1, metabolite2] > 0):  # first check if probs are tunable: must be distinct than zero both
#                 happy1, happy2 = complicated_stuff_1(pw_score, met_prob, spot1, spot2, pw1, pw2, metabolite1, metabolite2)
#                 if happy1 + happy2 > 0:
#                     new_probabs[spot1] = { (one_mz,pws_tuple): {metabolite1: (met_prob.loc[spot1, metabolite1], happy1),
#                                       metabolite2: (met_prob.loc[spot1, metabolite1], happy2)}}
#
#                 # we do not touch spot 2
##

pws01, anns01 = zm.SyncPathwaysAndAnnots.run(pathwaysObj.model_copy(),
                                             me_annots.model_copy(),
                                             min_pathway_size=1,
                                             max_pathway_size=1)

# print(len(pws01.piled['pathway_id'].tolist()))
# print(len(pws01.piled['pathway_id'].unique()))
lost_forever = adata_mz[:, anns01.df["mz"].astype(str).unique()].copy()
lost_forever = spatial_standard_computations(lost_forever, k=12)
plot_standard_computs(lost_forever, plot_folder=plot_folder, f"{SAMPLE}RejectedMz")

# sns.scatterplot(x=np.linspace(3, 10),  y=np.linspace(5, 10),
#               ax=ax1)