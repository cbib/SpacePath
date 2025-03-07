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
                             ]   ##  note: this is nice clusered but with artifact :  346.055808353, 300.048991657,  kept for now
wrong_mz_manually_defined = [str(z) for z in wrong_mz_manually_defined]



print(adata_mz)
import scanpy as sc
sc.pp.log1p(adata_mz)
import squidpy as sq
for i in wrong_mz_manually_defined:
    sq.pl.spatial_scatter(adata_mz, color=str(i), img=False, title=str(i),
                          save=os.path.abspath(os.path.join(plot_folder, f'{i}--BAD.png')))
    plt.show()