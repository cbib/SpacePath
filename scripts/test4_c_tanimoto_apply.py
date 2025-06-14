import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from rdkit import DataStructs
from rdkit import Chem
from rdkit.Chem import AllChem
import pickle

def give_ok_order(smiles_df, the_identifiers_may):
    mycopy = smiles_df.copy()
    mycopy = mycopy.drop_duplicates(subset=['moleculeIds', 'SMILES'], keep='first')
    mycopy = mycopy[['moleculeIds', 'SMILES']].copy()
    mycopy.index = mycopy['moleculeIds'].tolist()
    print("hourray ! ")
    mycopy = mycopy.T
    mycopy = mycopy[list(set(the_identifiers_may))]
    mycopy = mycopy.T

    return mycopy['moleculeIds'].tolist(), mycopy['SMILES'].tolist()


smiles_lipidmaps = pd.read_csv(
    "../resources/LipidMaps-SMILES_.06-08-25tsv",
     sep='\t', encoding='latin1',
    header=0)

smiles_hmdb =  pd.read_csv(
    "../resources/HMDB-SMILES_05-07-24.tsv",
     sep='\t', encoding='latin1',
    header=0)

smiles_df = pd.concat([smiles_lipidmaps, smiles_hmdb])


print("selecting reference with not none smiles")
smiles_df = smiles_df.loc[~smiles_df['SMILES'].isnull(),:].copy()

print("====<<<<<", smiles_df.loc[smiles_df['moleculeIds'] == 'HMDB0113179', "SMILES"])

##
myannots = pd.read_csv("identifs_for_tanimoto.csv", index_col=None)

the_identifiers_may = list(set(
    smiles_df['moleculeIds'].unique().tolist()
).intersection(set(
    myannots['moleculeIds'].unique().tolist()
)))
print(len(the_identifiers_may), ":) ")

the_identifiers_almost, the_smiles_good_order = give_ok_order(smiles_df, the_identifiers_may)
print(len(the_identifiers_almost), ":D")

print("*****", the_smiles_good_order[-7:])

yy = 99
arr = np.zeros(shape=(len(the_identifiers_almost),
                      len(the_identifiers_almost)))

print(arr.shape)

print(the_smiles_good_order[-7:])
print(the_identifiers_almost[-7:])

the_identifiers_that_enter = the_identifiers_almost#[]
the_smiles_enter = the_smiles_good_order#[]
#fps = []

fpgen = AllChem.GetRDKitFPGenerator()

the_ms = [Chem.MolFromSmiles(k) for k in the_smiles_enter]

fps = [fpgen.GetFingerprint(x) for x in the_ms]

# carajo = [fpgen.GetFingerprint(x) for x in ['[C@H](O)([C@H]1OC(=O)C=C1)CCCCCCCCCCCCCCC']]
# print(carajo, "fingerprint §§§§§§§")
# all_exeptions = 0
# for i in range(len(the_identifiers_almost)):
#     i_smile = the_smiles_good_order[i]
#     i_name = the_identifiers_almost[i]
#     try :
#         fps.append(fpgen.GetFingerprint(i_smile))
#         the_identifiers_that_enter.append(i_name)
#         the_smiles_enter(i_smile)
#     except:
#         all_exeptions += 1
#         continue

#print("oh no", all_exeptions)

print("filling the array")
for u in range(len(the_identifiers_that_enter)):
    for v in range(u + 1, len(the_identifiers_that_enter)):
        tanimoto_similarity = DataStructs.TanimotoSimilarity(fps[u], fps[v])
        arr[u,v] = tanimoto_similarity
        arr[v,u] = tanimoto_similarity

np.fill_diagonal(arr, 1)
dsave = {"arr": arr, "the_ids": the_identifiers_that_enter, "the_smiles": the_smiles_enter}
with open("tanimoto_similarities.pkl", "wb") as f:
    pickle.dump(dsave, f)

print("end")


