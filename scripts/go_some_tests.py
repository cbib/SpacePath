
import os
import pandas as pd
import numpy as np
import zonist.method as zm
import matplotlib.pyplot as plt
import seaborn as sns

from families_1_generate import differential_results_filter
import pickle



mydata_dir = "../../data"
SAMPLE = "LP-gut"

subfolder = SAMPLE

differential_results = pd.read_csv(os.path.join(
        mydata_dir, subfolder, f"{SAMPLE}-differential-results.csv"),
                                   index_col=None)

print(differential_results.shape)

differential_results = differential_results_filter(
    differential_results, ntop_per_cluster=5)

print(differential_results.shape)
print(differential_results.head())
print(differential_results[['group', 'logfoldchanges']])