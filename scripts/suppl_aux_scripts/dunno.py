import os
import pandas as pd
from local_helpers.quantifs_methods import differential_results_filter

mydata_dir = "../../../data"

# define here my input
my_csv_table_path_1 = "../../public-datasets-parameters.tsv"  #"../../small_test.csv"
my_csv_table_path_2 = "../../private-datasets-parameters.tsv"  #"../../small_test.csv"


df1 = pd.read_csv(my_csv_table_path_1, header=0,
                  index_col=None, sep='\t')
df2 = pd.read_csv(my_csv_table_path_2, header=0,
                  index_col=None, sep='\t')
datasets_parameters_df = pd.concat([df1,df2])

name_l = []
spot_n = []
features_n = []
for i, r in datasets_parameters_df.iterrows():
    SAMPLE = r["dataset_name"]

    subfolder = SAMPLE

    differential_results = pd.read_csv(os.path.join(
        mydata_dir, subfolder, f"{SAMPLE}-differential-results.csv"),
        index_col=None)

    with open(os.path.join(mydata_dir,
                           subfolder, f"{SAMPLE}-WIKI-NO.txt"), "r") as f:
        WIKI_NO = f.readlines()
    WIKI_NO.pop(0)
    WIKI_NO = [float(t.strip()) for t in WIKI_NO]

    with open(os.path.join(mydata_dir,
                           subfolder, f"{SAMPLE}-KEGG-NO.txt"), "r") as f:
        KEGG_NO = f.readlines()
    KEGG_NO.pop(0)
    KEGG_NO = [float(t.strip()) for t in KEGG_NO]

    ALL_MZ_NOT_MAPPED = set(WIKI_NO).union(set(KEGG_NO))

    print(
        "selecting the molecules based on the differential expression of their associated m/z")

    top_differential = differential_results_filter(differential_results)

    del differential_results

    ##
    o_fi1 = os.path.join(mydata_dir,
                         subfolder, f'{SAMPLE}-KEGG-CleanAnnots.tsv')
    o_fi2 = os.path.join(mydata_dir,
                         subfolder, f'{SAMPLE}-WIKI-CleanAnnots.tsv')

    annots_df1 = pd.read_csv(o_fi1,
                             sep='\t', header=0, index_col=None)
    annots_df2 = pd.read_csv(o_fi2,
                             sep='\t', header=0, index_col=None)

    print(annots_df1.shape, "KEGG")
    print(annots_df2.shape, "WIKI")
    print("I take so called 'wiki' annotations, because richer in lipids")

    print(
        "mz for further analyses: those not mapping to pathways AND differentially abundant between clusters")
    molecules_ids_df = pd.merge(top_differential, annots_df2, how='left',
                                left_on='names', right_on='mz')
    molecules_ids_df.dropna(how='any', inplace=True)
    molecules_ids_df = molecules_ids_df.loc[
                       molecules_ids_df['mz'].astype(float).isin(
                           list(ALL_MZ_NOT_MAPPED)), :]
    molecules_ids_df = molecules_ids_df[
        ['moleculeIds', 'mz', 'moleculeNames']].copy()

    molecules_ids_df.drop_duplicates(inplace=True)

    print("\nsample",SAMPLE, ":",  molecules_ids_df.shape[0], "\n")