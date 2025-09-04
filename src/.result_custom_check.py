
import numpy as np

import squidpy as sq

from single_slice_helpers import *

import matplotlib.pyplot as plt


np.random.seed(42)  # important to make the computations be reproducible


SAMPLE = "LP-cerebellum"#"TD-Lipido-635"
data_centralized_folder = "../../data"
subfolder = SAMPLE
pathways_types = ['FAMILIES', 'KEGG', 'WIKI']
parent_plot_folder = "../../figures"

clustering_parameters = {'choice': "kmeans", 'numeric_parameter': 4}
plot_folder = os.path.join(parent_plot_folder, SAMPLE)
from src.process_single_slice import save_top_plots, save_top_as_heterogeneous_files

# save_top_plots(pathways_types,
#                    clustering_parameters, SAMPLE, subfolder,
#                    data_centralized_folder, plot_folder)


# new:
save_top_as_heterogeneous_files(pathways_types,
                   clustering_parameters, SAMPLE, subfolder,
                   data_centralized_folder, plot_folder)


