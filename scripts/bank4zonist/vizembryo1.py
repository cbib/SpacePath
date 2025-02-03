import os
import squidpy as sq


def plot_selected_kegg_compound_ids(adata, me_annots, metabolites_to_throw_away):
    for item in metabolites_to_throw_away:
        item_name = me_annots.df.loc[
            me_annots.df["KEGG"] == item, "moleculeNames"].item()
        mz_here = me_annots.df.loc[me_annots.df["KEGG"] == item, "mz"].item()
        sq.pl.spatial_scatter(adata, color=[item],
                              img=False,
                              title=f"({mz_here}) {item}: {item_name}",
                              # save=os.path.abspath(tmpstr)
                              )
        plt.show()

def plot_cosine_distances_random_spots(spots_cos_sim_df):
    import seaborn as sns
    sns.displot(data=spots_cos_sim_df,
                x="x111_y0",
                bins=60, ec="white",
                kde=True)
    plt.show()
    random_spots = np.random.choice(a=spots_cos_sim_df.columns.tolist(),
                                    size=5, replace=False)
    for rspot in random_spots:
        sns.displot(data=spots_cos_sim_df,
                    x=rspot,
                    bins=60, ec="white",
                    kde=True)
        plt.show()

def plot_selected_pathways(adata_scores, biologically_chosen_pathways,
                           SAMPLEANDDETALS, plot_folder) -> None:
    for id in biologically_chosen_pathways.keys():
        try:
            pw_name = biologically_chosen_pathways[id]
            tmpstr = os.path.join(plot_folder,
                                  f'{SAMPLEANDDETALS}-{pw_name}.png')  # TODO: improve destination folder
            sq.pl.spatial_scatter(adata_scores, color=[id],
                                  img=False, title=f"{id} {pw_name}",
                                  save=os.path.abspath(tmpstr))
        except Exception as e:
            print("[plot_selected_pathways]Could NOT plot: {} {}".format(
                id, pw_name))