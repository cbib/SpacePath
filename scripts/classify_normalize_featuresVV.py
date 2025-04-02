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

from bank4zonist.quantifs_methods import (tic_normalization,
                                          custom_join_left, compute_moran_i,
                                          update_annotation_tomatch_adatamz)
from bank4zonist.constants import provide_annotation_dict

from bank4zonist.visualization import *


def categorize_intervals(values):
    """
    Categorize an array of values into specified intervals.
    """
    labels = ['0.05', '0.1', '0.2', '0.5']

    # Create an empty array to store the categorical labels
    categorical_labels = np.empty(len(values), dtype='<U10')

    # Categorize each value into the appropriate interval
    for i, value in enumerate(values):
        if 0 <= value <= 0.05:
            categorical_labels[i] = labels[0]
        elif 0.05 < value <= 0.1:
            categorical_labels[i] = labels[1]
        elif 0.1 < value <= 0.2:
            categorical_labels[i] = labels[2]
        elif 0.2 < value <= 0.501:
            categorical_labels[i] = labels[3]

    return np.array(categorical_labels)


def plot_mz_given_quality_metrics(df, SAMPLE, plot_folder,
                                  msm_cutoff, offSample_cutoff,
                                  text_notice):

    fileoutplotstr = os.path.join(plot_folder, f"{SAMPLE}-msm_&_out.pdf")
    #print(df['fdr'].astype(float).to_numpy())
    df['FDR'] = categorize_intervals(df['fdr'].astype(float).to_numpy())
    df['FDR'] = pd.Categorical(df['FDR'].tolist(),
                               categories= ['0.05', '0.1', '0.2', '0.5'])
    #print(df.head())

    sns.set_style("whitegrid")

    mypalette = ["#580F41", "#7A2F68", "#907359", "palegoldenrod"]
    #plt.figure(figsize=(5.5, 5))
    sns.scatterplot(df, x="msm", y="rawOffSampleProb",
                    hue="FDR",  s=16,
                    palette=sns.color_palette(mypalette, len(mypalette)),
                    ).legend(loc="upper left",title="FDR",
                             bbox_to_anchor=(1,1.02))
    plt.axvline(x=msm_cutoff, color='teal', linestyle='--')
    plt.axhline(y=offSample_cutoff, color='teal', linestyle='--')
    plt.title(SAMPLE + '\n')
    plt.tight_layout()
    plt.text(x=0.93, y=-0.12, s=text_notice,
             fontsize=9)
             #bbox=dict(facecolor='white', alpha=0.5))
    plt.savefig(fileoutplotstr)
    sns.set_style() # reset style
    plt.close()
    print(f"check the plot in {fileoutplotstr} (cutoffs as dotted lines)")


def plot__moran_info(df, plot_folder, SAMPLE):
    df["log_var_sim"] = np.log(df["var_sim"].to_numpy())
    df['padj class'] = f" > {moran_padj_cutoff}"
    df.loc[df[
                  "pval_sim_fdr_bh"] <= moran_padj_cutoff, 'padj class'] = f" <= {moran_padj_cutoff}"
    df['padj class'] = pd.Categorical(df['padj class'],
                                      categories=[f" <= {moran_padj_cutoff}",
                                                  f" > {moran_padj_cutoff}"] )
    sns.set_style("whitegrid")
    fig = sns.scatterplot(df, x="log_var_sim", y="I",
                          hue="padj class",
                          style="padj class",
                          palette=["teal", "slategray"],
                          markers=["o", "s"],
                          s=20
                          )
    fig.axhline(moran_I_cutoff, ls="--", alpha=0.5)
    fig.legend(loc="upper left", title="$p_{adj}$ class"
               )
    fileoutplotstr = os.path.join(plot_folder, f"{SAMPLE}-moran_scores.pdf")
    plt.title(SAMPLE)
    plt.savefig(fileoutplotstr)
    plt.close()
    sns.set_style()


def part_b_clean_with_moran_params(adata_mz, me_annots, moran_I_cutoff, moran_padj_cutoff):
    print(f"Part b) defining spatial quality of m/z features, by :\n"
        f"Moran's I value >= {moran_I_cutoff}\n"
          f"Adjusted p-value (BH) of Moran's I <= {moran_padj_cutoff}")

    if not "moranI" in list(adata_mz.uns.keys()):
        adata_mz = compute_moran_i(adata_mz, n_neighs=20, n_perms=50)

    plot__moran_info(adata_mz.uns['moranI'].copy(), plot_folder, SAMPLE)
    all_quality_df = adata_mz.var.copy()
    all_quality_df["mz"] = all_quality_df["mz"].astype(str).tolist()
    me_annots.mz_quality['mz'] = me_annots.mz_quality['mz'].astype(
        str).tolist()
    all_quality_df = custom_join_left(
        all_quality_df, me_annots.mz_quality[
            ['mz', 'msm',  'rawOffSampleProb', 'fdr',
            'quality_msm_off_fdr']],  keycol='mz')
    moran_mz_df = adata_mz.uns['moranI'].copy()
    moran_mz_df['mz'] = moran_mz_df.index.astype(str).tolist()
    all_quality_df = custom_join_left(all_quality_df, moran_mz_df[
        ['mz', 'I', 'pval_sim_fdr_bh']],  keycol='mz')

    keep_part2 = all_quality_df.loc[
        (all_quality_df['I'] >= moran_I_cutoff) & (
                all_quality_df['pval_sim_fdr_bh'] <= moran_padj_cutoff),
       "mz"].astype(str).tolist()

    if np.all(np.array(adata_mz.var.index.tolist()) == np.array(
            all_quality_df['mz'].tolist())):
            all_quality_df.index = all_quality_df['mz'].astype(str).tolist()
            all_quality_df['mz'] = all_quality_df['mz'].astype(float).tolist()
            adata_mz.var = all_quality_df
            adata_mz = adata_mz[:, keep_part2].copy()
            return adata_mz
    else:
        raise AssertionError


"""
clear off not-spatially relevant, or misleading, mz
"""
# define here my files and dirs
mydata_dir = "../../data"
my_csv_table_path = "../small-test-parameters.tsv" # "../datasets-parameters.tsv"
plot_parent_folder = "../../figures/"

overwrite_file = False  # todo : better define what to do about it


# start
datasets_parameters_df = pd.read_csv(my_csv_table_path, header=0,
                  index_col=None, sep='\t')

print("\nClear offSample and bad quality m/z across specified datasets")

for i, row in datasets_parameters_df.iterrows():
    SAMPLE = row["dataset_name"]
    normalize_yes_no = row["TIC_normalization"]
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

    adata_mz = ad.read_h5ad(filename=os.path.join(
        mydata_dir, subfolder, f'{SAMPLE}.h5ad'))
    print("loaded object", adata_mz.uns["history"])

    if 'features_n_after_mz_clean_2' in list(adata_mz.uns['history'].keys()):
        print("Nothing to do for this sample.")
    else:
        print(f"Part a) select m/z features with parameters: "
              f"MSM >= {msm_cutoff}, offSample <= {offSample_cutoff}")
        annotations_to_use = row["annotations_to_use"].split(" ")
        anno_dict = provide_annotation_dict(SAMPLE, annotations_to_use)
        me_annots = zd.MetabolomeAnnotations(
            name=SAMPLE, data_dir=mydata_dir,
            subfolder=subfolder, annotation_files=anno_dict)
        me_annots.load()  # load annotations
        if SAMPLE.startswith("MALDI_"):
            me_annots = update_annotation_tomatch_adatamz(me_annots, adata_mz)
        me_annots.pull_mz_quality_df()
        me_annots.mz_quality['quality_msm_off_fdr'] = "bad"
        me_annots.mz_quality.loc[
            (me_annots.mz_quality['msm'] >= msm_cutoff) & (
                    me_annots.mz_quality['rawOffSampleProb'] <= offSample_cutoff
            ) & (
                    me_annots.mz_quality['fdr'] <= fdr_cutoff), 'quality_msm_off_fdr'] = "good"

        ok_mz_keep_part1 = me_annots.mz_quality.loc[
            me_annots.mz_quality['quality_msm_off_fdr'] == "good", "mz"].astype(str).tolist()
        text_notice = f"admitted\n n = {len(ok_mz_keep_part1)}"
        plot_mz_given_quality_metrics(me_annots.mz_quality, SAMPLE,
                                      plot_folder,
                                      msm_cutoff, offSample_cutoff, text_notice)

        answer = input(f"A total of {len(set(ok_mz_keep_part1))} m/z remain after msm, offSample and fdr cutoffs."
                       f"\nDo you agree to continue?")
        if not answer.lower() in ["y", "yes"]:
            raise AssertionError("Aborted by user. Bye.")
        adata_mz = adata_mz[:, ok_mz_keep_part1].copy()
        adata_mz.uns["history"]['features_n_after_mz_clean_1'] = adata_mz.var.shape[0]
        if (normalize_yes_no.lower() in ['y', 'yes']) and (
                not "is_X_TIC_normalized" in adata_mz.uns.keys()):
            adata_mz.layers['raw_X'] = adata_mz.X.copy()
            adata_mz.X = csr_matrix(tic_normalization(
                adata_mz.X.copy(), adata_mz.obs['obs'].tolist(),
                adata_mz.var['mz'].tolist()))
            adata_mz.uns["is_X_TIC_normalized"] = True
            print("Updated .X with TIC normalized matrix (raw stored in .layer['raw_X']).")

        adata_mz = part_b_clean_with_moran_params(adata_mz, me_annots, moran_I_cutoff, moran_padj_cutoff)
        adata_mz.uns["history"]['features_n_after_mz_clean_2'] = adata_mz.var.shape[0]

        if overwrite_file:
            adata_mz.write_h5ad(filename=os.path.join(
            mydata_dir, subfolder, f"{SAMPLE}.h5ad"))
        else:
            adata_mz.write_h5ad(filename=os.path.join(
                mydata_dir, subfolder, f"{SAMPLE}-v2.h5ad"))
        print("saved object")



# adata_mz.write_h5ad(filename=os.path.join(
#    mydata_dir, subfolder, f"{SAMPLE}.h5ad"))

#################

# print(set(me_annots.mz_quality['mz'].astype(str).tolist()).intersection(
#             set(
#                 adata_mz.var.index.tolist()
#             )
#         ))
#         print(
#             len(set(me_annots.mz_quality['mz'].astype(str).tolist()).intersection(
#                 set(
#                     adata_mz.var.index.tolist()
#                 )
#             )))
#         df = closest_match_mz(me_annots.mz_quality['mz'].astype(str).tolist(),
#                               adata_mz.var.index.tolist())
