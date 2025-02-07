import os
import numpy as np
import anndata as ad
import squidpy as sq
import scanpy as sc
import pandas as pd
from typing import Dict, List
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
from sklearn.utils import sparsefuncs


def is_log_transformed(sparse_matrix):

    if sparse_matrix.min().min() < 0:
        return True

    return False


