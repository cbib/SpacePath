import os
import anndata as ad
import zonist.data as zd
import pandas as pd
import numpy as np
import scanpy as sc
import squidpy as sq
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt

mydata_dir = "../../../data"
plot_parent_folder = "../../figures/"

SAMPLE = "LP-brain-rat-2" #" #"TD-Lipido-635"    LP-brain-mouse-5


mydata_dir = "../../../data"
my_csv_table_path_1 = "../../public-datasets-parameters.tsv"  # "../../small_test.csv"
my_csv_table_path_2 = "../../private-datasets-parameters.tsv"  # "../../small_test.csv"

df1 = pd.read_csv(my_csv_table_path_1, header=0,
                  index_col=None, sep='\t')
df2 = pd.read_csv(my_csv_table_path_2, header=0,
                  index_col=None, sep='\t')
datasets_parameters_df = pd.concat([df1, df2])

me_refmap = zd.MetaboliteIdsMappings(
    data_dir=mydata_dir,
    subfolder="metabolites-mappings",
    mapping_files={'HMDB': "HMDB-to-KEGG_05-07-24.tsv",
                   'LipidMaps': "LipidMaps-to-KEGG_10-18-24.tsv",
                   'CoreMetabolome': "CoreMetabolome-to-KEGG_05-07-24.tsv"})
me_refmap.load()

import sys
sys.exit()


for i, r in datasets_parameters_df.iterrows():
    SAMPLE = r['dataset_name']
    subfolder = SAMPLE
    adata_mz = ad.read_h5ad(filename=os.path.join(
            mydata_dir, subfolder, f'{SAMPLE}.h5ad'))
    try:
        print(adata_mz.uns['norm'])
        print(adata_mz.uns['moranI'].head())
    except Exception as e:
        print(SAMPLE, "oh no here no raw, continue")
        continue
    print("")
    #print(SAMPLE, adata_mz.uns['history']['spot_number_retained'])

import sys
sys.exit()

if SAMPLE == "UBcalibr4":
    sq.pl.spatial_scatter(
        adata_mz, color=[str(715.575949)], img=False,
    )
    plt.show()
    sq.pl.spatial_scatter(
        adata_mz, color=[str(715.575949)], img=False,
        layer='raw_X'
    )
    plt.show()

# plot_folder = os.path.join(plot_parent_folder, SAMPLE)
#
# familiesObj = zd.MetabolicPathways(
#     data_dir=mydata_dir,
#     subfolder=subfolder,
#     file=f"{SAMPLE}-molecules_families.mmt"
# )
# familiesObj.molecule_id_type = "moleculeIds"
#
# familiesObj.load()
#
# #print(familiesObj)
#
# #print(len(familiesObj.to_dict()))

