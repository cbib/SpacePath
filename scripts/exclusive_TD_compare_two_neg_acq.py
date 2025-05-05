import sys
import os
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
import seaborn as sns
import squidpy as sq
import zonist.data as zd
import zonist.method as zm
import zonist.helpers as zh

from local_helpers.quantifs_methods import match_adata_mz_and_annots
from local_helpers.local_helpers import provide_annotation_dict


def only_good_adatamz(adata_mz):
    adata_mz = adata_mz[:,
    adata_mz.var.loc[adata_mz.var['final_quality'] == "good",
    "mz"].astype(str).tolist()].copy()

    return adata_mz


"""
"""
if __name__ == "__main__":

    print(
        "make sure already ran 'analyze_classify_spots' and 'classify_normalize_features'")

    # also add to the anndata.var if a feature was yes Y or no N in pathways that will pass.
    # that will allow to "paste" the mz or kegg that were not in the pw trick (? if makes sense such "paste" ?)

    # define here my files and dirs
    mydata_dir = "../../data"

    plot_parent_folder = "../../figures/"

    samples = ["TD-hemi-635-1-A", "TD-hemi-635-1-B"]   #  ["TD-635-1-A", "TD-635-1-B"]

    annotations_to_use_specifics = {"TD-hemi-635-1-A": ["HMDB", "LipidMaps", "CoreMetabolome"],
                                    "TD-hemi-635-1-B": ["HMDB", "LipidMaps", "CoreMetabolome"],
                                    "TD-635-1-A": ["HMDB", "KEGG", "CoreMetabolome"],
                                    "TD-635-1-B": ["HMDB", "KEGG", "LipidMaps", "CoreMetabolome"]}

    data1 = ad.read_h5ad(filename=os.path.join(
            mydata_dir, samples[0] , f'{samples[0]}.h5ad'))
    data2 = ad.read_h5ad(filename=os.path.join(
        mydata_dir, samples[1], f'{samples[1]}.h5ad'))


    print("MANDATORY: take away bad quality m/z")

    data1 = only_good_adatamz(data1)
    data2 = only_good_adatamz(data2)

    # annots 1 and 2
    print("processing annotations for sample:", samples[0])
    annots1 = zd.MetabolomeAnnotations(
            name=samples[0], data_dir=mydata_dir, subfolder=samples[0],
            annotation_files=provide_annotation_dict(samples[0], annotations_to_use_specifics[samples[0]]))

    annots1.load()

    print("processing annotations for sample:", samples[1])
    annots2 = zd.MetabolomeAnnotations(
        name=samples[1], data_dir=mydata_dir, subfolder=samples[1],
        annotation_files=provide_annotation_dict(samples[1],
                                                 annotations_to_use_specifics[
                                                     samples[1]]))
    annots2.load()

    annots1, data1 = match_adata_mz_and_annots(annots1, data1)
    annots2, data2 = match_adata_mz_and_annots(annots2, data2)

    print("first mz of both annotations:")
    print("   in ", samples[0], ":", annots1.no_translated["mz"].min())
    print("   in ", samples[1], ":", annots2.no_translated["mz"].min())


    print("subsetting both annotations to a common m/z range (between 100 and 400 m/z)")
    df1 = annots1.no_translated.loc[annots1.no_translated['mz'] >= 100,: ]
    df2 = annots2.no_translated.loc[annots2.no_translated['mz'] <= 400,: ]

    print("annotations in", samples[0], ", that are not in ", samples[1])
    print(len(set(df1["moleculeIds"].tolist()) - set(df2["moleculeIds"].tolist())))
    print(set(df1["moleculeIds"].tolist()) - set(
        df2["moleculeIds"].tolist()))
    thelist = list(set(df1["moleculeIds"].tolist()) - set(
        df2["moleculeIds"].tolist()))
    out_df = df1.loc[df1["moleculeIds"].isin(thelist), ["mz", "moleculeIds",
                                               "moleculeNames"]]
    out_df = out_df.drop_duplicates()
    print(out_df)

    # mz moleculeNames  moleculeIds

    #print(len(set(df1[""]) - set([""]))
    #print(len(set(annots1.no_translated[""]) - set(annots2.no_translated[""]))

    # me_refmap = zd.MetaboliteIdsMappings(
    #         data_dir=mydata_dir,
    #         subfolder="metabolites-mappings",
    #         mapping_files={'HMDB': "HMDB-to-KEGG_05-07-24.tsv",
    #                        'LipidMaps': "LipidMaps-to-KEGG_10-18-24.tsv",
    #                        'CoreMetabolome': "CoreMetabolome-to-KEGG_05-07-24.tsv"})
    # me_refmap.load()

