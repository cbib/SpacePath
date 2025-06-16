import os
import pandas as pd
import anndata as ad

import zonist.data as zd
import zonist.processing.measures_methods as zmm

from local_helpers.local_helpers import provide_annotation_dict


if __name__ == "__main__":

   # define here my files and dirs

    pwtype = "WIKI"  # set KEGG or WIKI  # TODO: think of a best alternative to do this

    mydata_dir = "../../data"


    my_csv_table_path_1 = "../public-datasets-parameters.tsv"  # "../../small_test.csv"
    my_csv_table_path_2 = "../private-datasets-parameters.tsv"  # "../../small_test.csv"

    df1 = pd.read_csv(my_csv_table_path_1, header=0,
                      index_col=None, sep='\t')
    df2 = pd.read_csv(my_csv_table_path_2, header=0,
                      index_col=None, sep='\t')
    datasets_parameters_df = pd.concat([df1, df2])

# start
    me_refmap = zd.MetaboliteIdsMappings(
        data_dir=mydata_dir,
        subfolder="metabolites-mappings",
        mapping_files={'HMDB': "HMDB-to-KEGG_05-07-24.tsv",
                       'LipidMaps': "LipidMaps-to-KEGG_10-18-24.tsv",
                       'CoreMetabolome': "CoreMetabolome-to-KEGG_05-07-24.tsv"})
    me_refmap.load()

    if pwtype == "WIKI":
        pathwaysObj = zd.MetabolicPathways(
            data_dir=mydata_dir,
            subfolder="metabolites-mappings",
            official_file="WIKIPATHWAYS_human_multi_id_redundant.mmt"
        )
        pathwaysObj.molecule_id_type = "moleculeIds"

    elif pwtype == "KEGG":
        pathwaysObj = zd.MetabolicPathways(
            data_dir=mydata_dir,
            subfolder="metabolites-mappings",
            official_file="KEGG_hsa_pathways_compounds_R110.gmt"
        )
        pathwaysObj.molecule_id_type = "KEGG"

    pathwaysObj.load()

    print("\nDo .....")

    for i, row in datasets_parameters_df.iterrows():
        SAMPLE = row["dataset_name"]
        min_pathway_size = row["min_pathway_size"]
        print("\nSample: ", SAMPLE)
        subfolder = SAMPLE

        annotations_to_use = row["annotations_to_use"].split(" ")
        anno_dict = provide_annotation_dict(SAMPLE, annotations_to_use)

        me_annots = zd.MetabolomeAnnotations(
            name=SAMPLE, data_dir=mydata_dir, subfolder=subfolder,
            annotation_files=anno_dict)

        o_fi = os.path.join(mydata_dir,
                             subfolder, f'{SAMPLE}-{pwtype}-CleanAnnots.tsv')
        #ok= True
        if not os.path.exists(o_fi):

            if pwtype == "WIKI":
                me_annots.usable_id = "moleculeIds"
            elif pwtype == "KEGG":
                me_annots.usable_id = "KEGG"

            me_annots.load()

            adata_mz = ad.read_h5ad(filename=os.path.join(
                mydata_dir, subfolder, f'{SAMPLE}.h5ad'))

            print(adata_mz.uns["history"])
            print("MANDATORY: take away bad quality m/z")  # todo: how to make an alert ?
            adata_mz = adata_mz[:,
                       adata_mz.var.loc[adata_mz.var['final_quality'] == "good",
                       "mz"].astype(str).tolist()]
            allgoodmz = adata_mz.var['mz'].astype(float).to_numpy()
            me_annots, adata_mz = zmm.match_adata_mz_and_annots(me_annots, adata_mz)

            if pwtype == "WIKI":
                me_annots.verify()
            elif pwtype == "KEGG":
                me_annots.convert_ids_to_kegg(me_refmap.table)

            me_annots.keep_unique_mz_id_tuples(adata_mz.var.copy())  #

            if pwtype == "WIKI":
                me_annots.remove_redundant_aliases()  # when using multitype 'moleculeIds' annotations

            me_annots.df.to_csv(o_fi, sep='\t', header=True, index=False)
            print("saved ", o_fi)

        else:
            print(f"{SAMPLE} : the clean annotations file already exists. Nothing to do.")
    print("done")


