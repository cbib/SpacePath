
import os
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
import squidpy as sq
import zonist.data as zd
import zonist.method as zm
import zonist.processing.measures_methods as zmm
from local_helpers.local_helpers import provide_annotation_dict

from zonist.processing.measures_methods import part_a__mz_quality_mark

if __name__ == "__main__":
    mydata_dir = "../../data"
    plot_parent_folder = "../../figures/"
    my_csv_table_path_1 = "../public-datasets-parameters.tsv"  # "../../small_test.csv"
    my_csv_table_path_2 = "../private-datasets-parameters.tsv"  # "../../small_test.csv"

    df1 = pd.read_csv(my_csv_table_path_1, header=0,
                      index_col=None, sep='\t')
    df2 = pd.read_csv(my_csv_table_path_2, header=0,
                      index_col=None, sep='\t')
    datasets_parameters_df = pd.concat([df1, df2])

    # start
    print("\nUUUUUUU")

    h_tuples = set()

    for i, row in datasets_parameters_df.iterrows():
        SAMPLE = row["dataset_name"]
        offSample_cutoff = row["offSample_cutoff"]
        msm_cutoff = row["msm_cutoff"]
        fdr_cutoff = row["fdr_cutoff"]
        print("\nSample: ", SAMPLE)
        subfolder = SAMPLE
        #plot_folder = os.path.join(plot_parent_folder, SAMPLE)
        annotations_to_use = row["annotations_to_use"].split(" ")
        anno_dict = provide_annotation_dict(SAMPLE, annotations_to_use)

        me_annots = zd.MetabolomeAnnotations(
            name=SAMPLE, data_dir=mydata_dir, subfolder=subfolder,
            annotation_files=anno_dict)
        me_annots.usable_id = 'moleculeIds'
        me_annots.load()
        me_annots.verify()
        me_annots.remove_redundant_aliases()

        me_annots, mz_to_keep_part1 = part_a__mz_quality_mark(
            me_annots, msm_cutoff, offSample_cutoff, fdr_cutoff
        )

        my_df = me_annots.df.loc[me_annots.df['mz'].astype(str).isin(
            mz_to_keep_part1
        ), :].copy()

        print(my_df.shape)
        print(my_df.head(2))

        my_df = my_df.loc[my_df['database'] == 'HMDB', ['moleculeIds', 'moleculeNames']]
        my_df = my_df.drop_duplicates()
        myly = my_df.values.tolist()
        myly = set([tuple(tup) for tup in myly])
        h_tuples = h_tuples.union(myly)

    mytext = ""
    for u in list(h_tuples):
        mytext += f'{u[0]};{u[1]}\n'

    with open("zoomyhey.txt", "w") as f:
        f.write(mytext)




