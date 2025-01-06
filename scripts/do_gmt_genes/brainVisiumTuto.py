#  https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_visium_hne.html

import numpy as np
import pandas as pd

import anndata as ad
import scanpy as sc
import squidpy as sq



# load the pre-processed dataset
img = sq.datasets.visium_hne_image()
adata = sq.datasets.visium_hne_adata()   # NOTE: too slow, and need to test sevearl times the scatter
# went to notebooks/appendix to do it there

sq.pl.spatial_scatter(adata,
                      color="cluster")