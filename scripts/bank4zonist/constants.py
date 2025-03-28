import os
import numpy as np
import anndata as ad
import pandas as pd
import scanpy as sc
from typing import List
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from sklearn.utils import sparsefuncs
import squidpy as sq

"""
get the matrix of KEGG given the annotations, and the adata_mz
"""




def yield_specific_typical_pathways() -> dict:
    biologically_chosen_pathways = {
        "hsa00600": "Sphingolipid metabolism",
        "hsa00030": "Pentose phosphate pathway",
        "hsa00010": "Glycolysis Gluconeogenesis", # attention, name here has been edited to avoid "/" in the middle
        "hsa00620" : "Pyruvate metabolism",
        "hsa00020"  : "TCA cycle",
        "hsa00900": "Terpenoid backbone biosynthesis",
        "hsa04727": "GABAergic synapse"}
    return biologically_chosen_pathways


