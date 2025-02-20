import sys
import os
import numpy as np
import pandas as pd
import anndata as ad
import squidpy as sq
import networkx as nx
from tqdm import tqdm
import matplotlib.pyplot as plt
import zonist.method as zm


SAMPLE = "UBcalibr2"
#"UBcalibr2"
print(SAMPLE)
mydata_dir = "/home/johanna/spatial_thesis/spatial_metabolomics/data"
#"/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"
subfolder = "UBplatform"
#"/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"



adata = ad.read_h5ad(os.path.join(mydata_dir, subfolder,
                                         f"{SAMPLE}-KEGG_adata.h5ad"))


scores = ad.read_h5ad(os.path.join(mydata_dir, subfolder,
                                         f"{SAMPLE}-PathwayScores.h5ad"))

# look metabolites of
my_pathway= "hsa00900"
tmp = scores.var.copy()
my_name = tmp.loc[tmp['pathway_id'] == my_pathway, "pathway_name"]
print(my_name)
joo = scores.uns["current_pathways_df"].loc[
    scores.uns["current_pathways_df"]["KEGG"] == "C00332", :
]
print(joo)
for metaboliteid in scores.uns["former_pathways_dict"][my_pathway]:
    try:
        # mymz = scores.uns["current_pathways_df"].loc[
        #     scores.uns["current_pathways_df"]["KEGG"] == metaboliteid, "mz"
        # ].item()
        # print(metaboliteid, "===", mymz)
        sq.pl.spatial_scatter(adata, color=[metaboliteid], img=False)
        plt.show()
        print("in data: ", metaboliteid)
        #plt.close()
    except:
        plt.close()

#for feature in scores.uns['']


