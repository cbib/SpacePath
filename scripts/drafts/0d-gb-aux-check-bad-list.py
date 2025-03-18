import sys
import os
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import zonist.data as zd
import zonist.method as zm

sys.path.append(os.path.join(os.path.dirname(__file__),
                             '../bank4zonist'))

sys.path.append(os.path.join(os.path.dirname(__file__),
                             '../'))

from user_routines import prepare_harmonized_objects

from quantifs_methods import *
from visualization import *


"""only for file with mz range 80-400"""

SAMPLE = "test635_80-400"

print(SAMPLE)
mydata_dir = "/home/johanna/spatial_thesis/spatial_metabolomics/data"
#"/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"
subfolder = "gb-brain"

plot_folder = "../../../figures/gb-brain"


assert os.path.exists(plot_folder), "Cannot continue, the plot_folder does not exist"



adata_mz = ad.read_h5ad(filename=os.path.join(
    mydata_dir, subfolder, f'{SAMPLE}.h5ad'))
  ##  note: this is nice clusered but with artifact :  346.055808353, 300.048991657,  kept for now
wrong_mz_manually_defined = [str(z) for z in wrong_mz_manually_defined]



print(adata_mz)
import scanpy as sc
sc.pp.log1p(adata_mz)
import squidpy as sq
for i in wrong_mz_manually_defined:
    sq.pl.spatial_scatter(adata_mz, color=str(i), img=False, title=str(i),
                          save=os.path.abspath(os.path.join(plot_folder, f'{i}--BAD.png')))
    plt.show()


## aux: visual picking todo: think if doing appart
if do_step_viz_pick_wrong:  # step 1: manually complete wrong list:
    def plotem(adata_mz, me_annots, good_mz_spatial) -> None:
        import scanpy as sc
        import squidpy as sq
        adata_tmp = adata_mz.copy()

        adata_tmp = adata_tmp[:, good_mz_spatial]
        pathwaysObj, me_annots, adata_tmp, initial_pws_dic = prepare_harmonized_objects(
            mydata_dir, me_annots, adata_tmp, prob_out=0.9,
            min_pw_size=2)
        foo = pathwaysObj.piled["mz"].unique()
        foo = np.array(foo)
        print(
            f"visualizing {len(foo)} the features  (only those mapping to any pathway)"
            "after visualizing, please manually complete the list 'wrong_mz_manually_defined'")
        sc.pp.log1p(adata_tmp)
        for i in foo:
            sq.pl.spatial_scatter(adata_tmp, color=str(i), img=False,
                                  title=str(i))
            plt.show()
        # returns none


    plotem(adata_mz, me_annots, good_mz_spatial)

    answer = input("finished step1, continue?")
    if not answer.lower() in ["y", "yes"]:
        print("Aborted")
    # end if
# endif