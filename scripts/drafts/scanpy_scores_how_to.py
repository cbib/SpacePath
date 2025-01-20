
#   this chunk was before compute kpca in 1-platform-pixel-select.py
###
# compute scanpy module score
##
def use_scanpy_scores(adata, pathways: dict):
    coherced_pathways = pathways_dict_coherce_min_len(
        pathways, min_len_admitted=2
    )
    df = modules_scores_by_scanpy(adata, coherced_pathways, do_log=True)

    return df

scores_by_scanpy_df = use_scanpy_scores(adata, pathwaysObj.pathways)
# TODO: is it always needed to shift the distr to the right ???
scores_by_scanpy_df = scores_by_scanpy_df + abs(scores_by_scanpy_df.min().min())  # shift to right

adata_scores_scanpy = assemble_adata_from_scores_result(adata, scores_by_scanpy_df)

plot_selected_pathways(adata_scores_scanpy, biologically_chosen_pathways,
                  f"{SAMPLE}-ScoresScanpy", plot_folder)