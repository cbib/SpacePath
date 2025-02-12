import sys
import os
import numpy as np
import pandas as pd
import anndata as ad

import matplotlib.pyplot as plt
import zonist.data as zd
import zonist.method as zm

import squidpy as sq
import networkx as nx
from tqdm import tqdm
from sklearn.metrics.pairwise import cosine_similarity

from scripts.user_routines import prepare_harmonized_objects
from scripts.bank4zonist.quantifs_methods import find_category_by_nb_kegg

sys.path.append(os.path.join(os.path.dirname(__file__),
                             '../bank4zonist'))
from quantifs_methods import *

"""
using presence-absence of pathways: deprecated
"""

##  copied chunk
SAMPLE = "UBcalibr2"#sys.argv[1]
print(SAMPLE)
mydata_dir = "/home/johanna/spatial_thesis/spatial_metabolomics/data"
subfolder = "UBplatform"
anno_dict = {'HMDB': f"{SAMPLE}_hmdb.csv",
             "LipidMaps": f"{SAMPLE}_lipidmaps.csv"}

pathways_scores_df = pd.read_csv(os.path.join(mydata_dir,subfolder,
                                f"{SAMPLE}-PathwaysScores.tsv"), sep='\t',
                                 index_col=0)


plot_folder = "../../../figures/platform"


assert os.path.exists(plot_folder), "Cannot continue, the plot_folder does not exist"

me_annots = zd.MetabolomeAnnotations(name=SAMPLE,
                            data_dir=mydata_dir,
                            subfolder=subfolder,
                            annotation_files=anno_dict)
me_annots.load()

adata_mz = ad.read_h5ad(filename=os.path.join(
    mydata_dir, subfolder, f'{SAMPLE}_mx.h5ad'))

verification(adata_mz, me_annots)


pathwaysObj, me_annots, adata_mz, orig_pws_dict = prepare_harmonized_objects(
    mydata_dir, me_annots, adata_mz)
print("mz post sync", adata_mz)


adata = build_kegg_resolved_adata(adata_mz, annots_df=me_annots.df)
print("KEGG resolved adata: ", adata)

## end copied chunk

X_mz_binary_df = compute_mz_binary_df(adata_mz)
X_metabolites_probs_df = compute_metabolites_probs_df(X_mz_binary_df,
                                                      me_annots.df)
X_presence_absence_mask = pathways_scores_df >= 0  # score of 0 or more is presence (negative scores are absence)
Xbinary_pathways = pd.DataFrame(
    np.array(X_presence_absence_mask, dtype=int),
    columns=pathways_scores_df.columns.tolist(),
    index=pathways_scores_df.index.tolist())


#is_it_well_done(X_metabolites_probs_df, adata_mz, adata, me_annots.df)

refinement_mz_pws_dict = ugly_long_function(pathwaysObj.piled.copy(), pathwaysObj.to_dict())

for t in refinement_mz_pws_dict.keys():
    print(t)
    print("   ", refinement_mz_pws_dict[t])
# Note: the refinement_mz_pws_dict has all sense and anti-sense pathways comparisons (pw1 vs pw2) and (pw2 vs pw1)
the_mz_for_refinement_list = [float(i) for i in list(refinement_mz_pws_dict.keys())]

########
# spatial info
RADIUS = 1
sq.gr.spatial_neighbors(adata, radius=RADIUS, coord_type="generic")
connectivities_sparse_mx = adata.obsp["spatial_connectivities"]

coordinates_names_str = adata.obs["obs"].tolist()
coords_matcher_dict = {k: v for k, v in zip(
    range(0, len(coordinates_names_str)), coordinates_names_str)}
coords_matcher_rev = {v: k for k, v in coords_matcher_dict.items()}

print(f"parsing spots : for each spots, neighbors at d = {RADIUS} are evaluated")
COSINE_SIMILARITY_THRESHOLD = 0.6
G = nx.from_scipy_sparse_array(connectivities_sparse_mx)
nodes_list = G.nodes()

categ_will_modif = dict()

for v1 in tqdm(nodes_list):
    total_kegg_will_change = set()
    neighs_of_v1 = G.neighbors(v1)
    tmp_mx1 = adata[coords_matcher_dict[v1], :].X.todense()
    spot1 = coords_matcher_dict[v1]
    for v2 in neighs_of_v1:
        tmp_mx2 = adata[coords_matcher_dict[v2], :].X.todense()
        local_cosine_sim = cosine_similarity(np.array(tmp_mx1),
                                             np.array(tmp_mx2)) # 1 x 1 mx
        local_cosine_sim = float(local_cosine_sim.reshape(
            local_cosine_sim.shape[1])[0])
        if local_cosine_sim >= COSINE_SIMILARITY_THRESHOLD:
            spot2 = coords_matcher_dict[v2]
            modifiablein_s1_whencomparing_s2 = evaluate_spot1_change(refinement_mz_pws_dict,
                                       Xbinary_pathways,
                                       X_metabolites_probs_df,
                                       spot1, spot2)
            total_kegg_will_change = total_kegg_will_change.union(
                modifiablein_s1_whencomparing_s2
            )

    categ_will_modif[spot1] = str(find_category_by_nb_kegg(len(total_kegg_will_change))) # unchanged, change1to5Kegg, change6to10Kegg, change11andmoreKegg



oh = pd.DataFrame({
    'obs': [tup[0] for tup in categ_will_modif.items()],
    'will_modif': [str(tup[1]) for tup in categ_will_modif.items()]
})


if np.all(adata.obs['obs'].to_numpy() == oh['obs'].to_numpy()):
    adata.obs = adata.obs.assign(will_modif=oh["will_modif"].tolist())
    print("ok !")
else:
    raise ValueError("Abort ! ")

sq.pl.spatial_scatter(adata, color='will_modif',
                              img=False, title="weird analysis",
                             # save=os.path.abspath(tmpstr)
                      )
plt.show()




