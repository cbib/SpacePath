
import numpy as np

import squidpy as sq

from analysis_helpers import *




np.random.seed(42)  # important to make reproducible the computations


SAMPLE = "LP-cerebellum-rat"
data_centralized_folder = "../../data"
subfolder = SAMPLE
pathways_types = ['FAMILIES', 'KEGG', 'WIKI']
#str_out = "leiden_sg-0.4"
str_out = "kmeans-4"
diff_sets = pd.read_csv(os.path.join(data_centralized_folder, subfolder,
                                         f'{SAMPLE}-top-differential-SETS.tsv'),
                            sep='\t', header=0)
list_of_scores_objects = list()
for u in pathways_types:
    adata_scores = ad.read_h5ad(
        filename=os.path.join(data_centralized_folder, subfolder,
                              f'{SAMPLE}-scores-{u}.h5ad'))
    adata_scores.var['set_type'] = u
    list_of_scores_objects.append(adata_scores)
joined = ad.concat(adatas=list_of_scores_objects, join='outer', axis=1)

clusters_df = pd.read_csv(os.path.join(data_centralized_folder, subfolder,
                                       f'{SAMPLE}-clusters-{str_out}.txt'),
                          header=0)
joined.obs[str_out] = pd.Categorical(clusters_df[str_out].tolist(),
                                     categories=clusters_df[
                                         str_out].unique())
print(diff_sets.head())
diff_sets = diff_sets.loc[diff_sets['group'].astype(str) == '2',:].copy()
print(diff_sets)
print(diff_sets.shape)
print(joined.var[['pathway_id', 'pathway_name']].head(20))
joined.obsm = adata_scores.obsm
joined.uns['spatial'] = {'my_slice_name': {
        'scalefactors': {'fiducial_diameter_fullres': 100,
                         'spot_diameter_fullres': 1,
                         'tissue_hires_scalef': 0.5,
                         'tissue_lowres_scalef': 0.5}}}

print(diff_sets.names.unique().tolist())
pid = joined.var['pathway_id'].tolist()
print(joined.var['pathway_name'])
pnm = [t for t in joined.var['pathway_name'].astype(str).to_numpy()]
for set_id, set_full_name in zip(pid, pnm) :
    #set_full_name =
    print(set_full_name)
    sq.pl.spatial_scatter(
        joined, color=[set_id], img=False,
        title=f'{set_full_name[:45]}\n{set_id}',
        # save=os.path.abspath(plot_out_path),  dpi=600
        )
    plt.show()


