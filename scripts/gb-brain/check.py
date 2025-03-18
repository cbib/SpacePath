

# check if mz in one annotation are not in the other
import sys
import os
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
import squidpy as sq
import zonist.data as zd
import zonist.method as zm

from scripts.bank4zonist.quantifs_methods import tic_normalization

sys.path.append(os.path.join(os.path.dirname(__file__),
                             '../bank4zonist'))

sys.path.append(os.path.join(os.path.dirname(__file__),
                             '../'))

from user_routines import prepare_harmonized_objects

from quantifs_methods import *
from visualization import *

"""
clear off not-spatially relevant, or misleading, mz
"""

SAMPLE = "test635_80-400"   # "test635_200-1K"  sys.argv[1]## todo: pending to sweep the larger mz range file


print(SAMPLE)
mydata_dir = "/home/johanna/spatial_thesis/spatial_metabolomics/data"
#"/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data"
subfolder = "gb-brain"
anno_dict = {'HMDB': f"{SAMPLE}_hmdb.csv",
             "KEGG": f"{SAMPLE}_kegg.csv",
             "CoreMetabolome": f"{SAMPLE}_CoreMetabolome.csv"}

me_annots = zd.MetabolomeAnnotations(name=SAMPLE,
                            data_dir=mydata_dir,
                            subfolder=subfolder,
                            annotation_files=anno_dict)
me_annots.load()
print(me_annots.preloaded['database'].unique())
core = me_annots.preloaded.loc[me_annots.preloaded["database"]=="CoreMetabolome",:]
hmdb = me_annots.preloaded.loc[me_annots.preloaded["database"]=="HMDB",:]

print(set(core["mz"].unique())- set(hmdb["mz"].unique()))   # ok zero