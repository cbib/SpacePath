import os
import squidpy as sq
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix
import zonist.data as zd
import zonist.method as zm
import matplotlib.pyplot as plt

from bank4zonist.quantifs_methods import match_adata_mz_and_annots, verification

def prepare_harmonized_objects(mydata_dir, me_annots, adata_mz, prob_out=0.87,
                               min_pw_size=1, max_pw_size=100):  # todo: deprecate this!
    # reference mappings are the same for all my projects
    me_refmap = zd.MetaboliteIdsMappings(
        data_dir=mydata_dir,
        subfolder="metabolites-mappings",
        mapping_files={'HMDB': "HMDB-to-KEGG_05-07-24.tsv",
                       'LipidMaps': "LipidMaps-to-KEGG_10-18-24.tsv",
                       'CoreMetabolome': "CoreMetabolome-to-KEGG_05-07-24.tsv"})

    me_refmap.load()
    ##
    me_annots, adata_mz = match_adata_mz_and_annots(me_annots, adata_mz)
    verification(adata_mz, me_annots)  # todo del
    print("Excludes m/z with out-of-tissue probability >=", prob_out)
    me_annots.discard_out_of_tissue_mz(prob_out=prob_out)
    me_annots.convert_ids_to_kegg(me_refmap.table)
    me_annots.keep_unique_mz_id_tuples(adata_mz.var.copy())

    ## pathways are also the same for all my projects
    pathwaysObj = zd.MetabolicPathways(
        data_dir=mydata_dir,
        subfolder="metabolites-mappings",
        official_file="KEGG_hsa_pathways_compounds_R110.gmt")

    pathwaysObj.load()
    initial_pws_dict = pathwaysObj.to_dict()  # reserver pour plustard ? simplement avoir les tailles des sets du départ.

    pathwaysObj, me_annots = zm.SyncPathwaysAndAnnots.run(pathwaysObj,
                                                          me_annots,
                                                          min_pathway_size=min_pw_size,
                                                          max_pathway_size=max_pw_size,
                                                          mode=TOIMPLEMENT)
    ll = [str(z) for z in sorted(me_annots.df["mz"].astype(float).unique())]

    adata_mz = adata_mz[:, ll].copy()  # todo this is NOT sync, but like so

    return pathwaysObj, me_annots, adata_mz, initial_pws_dict


def plot_susceptible_of_refinement(refinement_mz_pws_dict, adata_mz, SAMPLE, plot_folder):
    the_mz_for_refinement_list = list(refinement_mz_pws_dict.keys())
    adata_mz = compute_moran_i(adata_mz, 12)
    sub_moran = adata_mz.uns['moranI'].loc[
        adata_mz.uns['moranI'].index.isin(
            [str(x) for x in the_mz_for_refinement_list])]

    sub_moran.to_csv(
        os.path.join(plot_folder, f'{SAMPLE}--moran-subset.tsv'),
        header=True, index=True, sep='\t')

    tmpstr = os.path.join(plot_folder,
                          f'{SAMPLE}--susceptible_of_refinement.png')
    sq.pl.spatial_scatter(adata_mz,
                          color=[str(x) for x in the_mz_for_refinement_list],
                          img=False,
                          save=os.path.abspath(tmpstr))


