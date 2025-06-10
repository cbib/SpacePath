
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


import pickle
import json

df = pd.read_csv("tanimoto_groups.mmt",
                 sep='\t', header=None)
df = df.iloc[:, 2:].T
df = df.fillna(0)


myannots = pd.read_csv("identifs_for_tanimoto.csv", index_col=None)

hmdbs = []
lipids = []

for i in range(df.shape[1]):
    elems = df.iloc[:, i].tolist()
    elems = [el for el in elems if el != 0]
    hmdbs += [h for h in elems if h.startswith('HMDB')]
    lipids += [h for h in elems if not h.startswith('HMDB')]

fake_drugs = ['HMDB0015486', 'LMGP01011742']

print(len(hmdbs), "total of hmdbs")
print(len(lipids), "total of lipids in lipidmaps")

for g in hmdbs:
    name = myannots.loc[myannots['moleculeIds']==g,'moleculeNames'].item()
    if (name.startswith("SM(")) or name.startswith("PE(") or (
        name.startswith("PE-")) or name.startswith("PC(") or (
            name.startswith("PC-")) or name.startswith("Lyso") or (
            name.startswith("PS(")) or name.startswith("PEP") or (
            name.startswith("PA")) or name.startswith("TG(")  or (
            name.startswith("PG")) or name.startswith("PIP") or name.startswith("PI("):
        pass
    elif name.startswith("LMPK12"):
        print('flavonoid!', g, name)
    else:
        print(g, name)



