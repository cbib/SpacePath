import os
import pandas as pd
import anndata as ad

import zonist.data as zd
import zonist.method as zm
import zonist.processing.measures_methods as zmm
from local_helpers.local_helpers import provide_annotation_dict

"""
Saves updated adata_mz object with classified m/z by the quality of the assigned formulas
Specifically, adata_mz.var acquires new columns, with a helpful column 'final_quality' 
"""

if __name__ == "__main__":

    # define here my files and dirs
    mydata_dir = "../../data"
    my_csv_table_path_1 = "../public-datasets-parameters.tsv"  #
    my_csv_table_path_2 = "../small-test-parameters.tsv"  # ../private-datasets-parameters.tsv"  # "

    df1 = pd.read_csv(my_csv_table_path_1, header=0,
                      index_col=None, sep='\t')
    df2 = pd.read_csv(my_csv_table_path_2, header=0,
                      index_col=None, sep='\t')

    datasets_parameters_df = pd.concat([df1, df2])


    plot_parent_folder = "../../figures/"

    force_for_these = ["LP-brain-rat-2"]    # keep empty if no sample to force   # e.g.  "TD-Lipido-635", "TD-Lipido-636", "TD-Lipido-836", "TD-Lipido-827", "WW-breast-cancer-1"
    # start

    print("\nMarks bad quality m/z across specified datasets (.var property)")

    for i, row in datasets_parameters_df.iterrows():
        SAMPLE = row["dataset_name"]
        offSample_cutoff = row["offSample_cutoff"]
        msm_cutoff = row["msm_cutoff"]
        fdr_cutoff = row["fdr_cutoff"]
        moran_I_cutoff = row["Moran_I_cutoff"]
        moran_padj_cutoff = row["Moran_padj_cutoff"]
        print("\nSample: ", SAMPLE)

        subfolder = SAMPLE
        plot_folder = os.path.join(plot_parent_folder, SAMPLE)
        if not os.path.exists(plot_folder):
            os.mkdir(plot_folder)

        annotations_to_use = row["annotations_to_use"].split(" ")
        anno_dict = provide_annotation_dict(SAMPLE, annotations_to_use)
        me_annots = zd.MetabolomeAnnotations(
            name=SAMPLE, data_dir=mydata_dir,
            subfolder=subfolder, annotation_files=anno_dict)

        adata_mz = ad.read_h5ad(filename=os.path.join(
            mydata_dir, subfolder, f'{SAMPLE}.h5ad'))


        print("Loaded my intensities data:", adata_mz.uns["history"])

        if (SAMPLE not in force_for_these) and ('features_n_after_mz_clean_2' in list(adata_mz.uns['history'].keys())):
            print("Nothing to do for this sample.")
            print(adata_mz.layers.keys())
        else:
            try:
                columns_to_reset = set(adata_mz.var.columns.tolist()) - set(["mz"])
                adata_mz.var = adata_mz.var.drop(
                    columns=list(columns_to_reset))
            except KeyError:
                pass

            me_annots.load()  # load annotations
            me_annots, adata_mz = zmm.match_adata_mz_and_annots(me_annots,
                                                            adata_mz)
            adata_mz = zm.NormalizeIntensities.run(adata_mz)

            adata_mz = zm.ClassifyMZProtocol.run(adata_mz,
                me_annots, msm_cutoff, offSample_cutoff, fdr_cutoff,
                moran_I_cutoff, moran_padj_cutoff,
                SAMPLE, plot_folder
            )
            print(adata_mz)
            print(adata_mz.uns['history'])
            adata_mz.write_h5ad(filename=os.path.join(
                mydata_dir, subfolder, f"{SAMPLE}.h5ad"))
            print("saved object")


##

    # my_csv_table_path = "../private-datasets-parameters.tsv"   #  small-test-parameters.tsv"

    # datasets_parameters_df = pd.read_csv(my_csv_table_path, header=0,
    #                   index_col=None, sep='\t')
