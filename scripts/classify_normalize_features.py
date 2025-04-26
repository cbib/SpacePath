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
                                          match_adata_mz_and_annots)
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
                                  text_notice) -> None:
    fdr_categories = ['0.05', '0.1', '0.2', '0.5']
    df['FDR'] = categorize_intervals(df['fdr'].astype(float).to_numpy())
    df['FDR'] = pd.Categorical(df['FDR'].tolist(),
                               categories= fdr_categories)
    sns.set_theme(style="whitegrid", font_scale=1.2)
    mypalette = ["#580F41", "#74414D",  "#907359", "#cfcb94"]
    ax = sns.scatterplot(df.loc[df['FDR'] == '0.5', :],
                         x="msm", y="rawOffSampleProb",
                         hue="FDR", alpha=.9, zorder=1,
                         s=30,
                         palette=sns.color_palette(mypalette, len(mypalette)),
                         legend=False)
    sns.scatterplot(df.loc[df['FDR'] != '0.5', :],
                    ax=ax, x="msm", y="rawOffSampleProb",
                    hue="FDR",  alpha=.9, zorder=2,
                    s=30,
                    palette=sns.color_palette(mypalette, len(mypalette)),
                    legend=True)
    ax.legend(loc="upper left",title="FDR",
                             bbox_to_anchor=(1,1.03))  #.02))
    plt.xlim(-0.04, 1.04)
    plt.ylim(-0.04, 1.04)
    plt.axvline(x=msm_cutoff, color='teal', linestyle='--', lw=0.8)
    plt.axhline(y=offSample_cutoff, color='teal', linestyle='--', lw=0.8)
    plt.title(SAMPLE + '\n', fontsize=18)
    plt.text(x=1.1, y=0.1, s=text_notice,
             fontsize=15)  #bbox=dict(facecolor='white', alpha=0.5))
    plt.tight_layout()
    fileoutplotstr = os.path.join(plot_folder, f"{SAMPLE}-msm_&_out.pdf")
    plt.savefig(fileoutplotstr)
    plt.close()
    print(f"check the plot in {fileoutplotstr} (cutoffs as dotted lines)")


def part_a_filters(me_annots, msm_cutoff, offSample_cutoff, fdr_cutoff,
                   SAMPLE, plot_folder):
    print(f">> Part a): Select m/z features with parameters: "
          f"MSM >= {msm_cutoff}, offSample <= {offSample_cutoff}, FDR <= {fdr_cutoff}")

    me_annots.pull_mz_quality_df()
    me_annots.mz_quality['quality_msm_off_fdr'] = "bad"
    me_annots.mz_quality.loc[
        (me_annots.mz_quality['msm'] >= msm_cutoff
         ) & (me_annots.mz_quality['rawOffSampleProb'] <= offSample_cutoff
         ) & (me_annots.mz_quality['fdr'] <= fdr_cutoff),
        'quality_msm_off_fdr'] = "good"

    ok_mz_keep_part1 = me_annots.mz_quality.loc[
        me_annots.mz_quality['quality_msm_off_fdr'] == "good", "mz"].astype(
        str).tolist()
    text_notice = f"valid $m/z$\n n = {len(ok_mz_keep_part1)}"
    plot_mz_given_quality_metrics(me_annots.mz_quality, SAMPLE,  # todo: take the plot part outside, and call correctly in all usages
                                  plot_folder,
                                  msm_cutoff, offSample_cutoff, text_notice)
    print(f"A total of {len(set(ok_mz_keep_part1))} m/z pass MSM, offSample and FDR cutoffs.")

    return me_annots, ok_mz_keep_part1


def plot__moran_info(df, plot_folder, SAMPLE, text_notice) -> None:
    df["log_var_sim"] = np.log(df["var_sim"].to_numpy())
    df['padj class'] = f" > {moran_padj_cutoff}"
    df.loc[df[
                  "pval_sim_fdr_bh"] <= moran_padj_cutoff, 'padj class'] = f" <= {moran_padj_cutoff}"
    df['padj class'] = pd.Categorical(df['padj class'],
                                      categories=[f" <= {moran_padj_cutoff}",
                                                  f" > {moran_padj_cutoff}"] )
    sns.set_theme(style="whitegrid", font_scale=1.2)
    plt.figure(figsize=(5.2, 5.05))
    fig = sns.scatterplot(df, x="log_var_sim", y="I",
                          hue="padj class",
                          style="padj class", lw=0,
                          palette=["teal", "slategray"],
                          markers=["o", "s"], alpha=0.4,
                          s=20)
    fig.axhline(moran_I_cutoff, ls="--", lw=0.7, color='teal')
    fig.legend(loc="upper left", title="$p_{adj}$ class")
    plt.title(SAMPLE + "\n" + text_notice + "\n" )
    plt.tight_layout()
    fileoutplotstr = os.path.join(plot_folder, f"{SAMPLE}-moran_scores.pdf")
    plt.savefig(fileoutplotstr)
    plt.close()


def part_b_filters_with_autocorr(adata_mz, ok_part1, me_annots, moran_I_cutoff, moran_padj_cutoff):
    print(f">> Part b): Classifying m/z features, by :\n"
          f"Moran's I value >= {moran_I_cutoff}, "
          f"(BH) adjusted p-value <= {moran_padj_cutoff}")
    # prepare quality merged df
    all_quality_df = adata_mz.var.copy()
    all_quality_df["mz"] = all_quality_df["mz"].astype(str).tolist()
    me_annots.mz_quality['mz'] = me_annots.mz_quality['mz'].astype(
        str).tolist()
    all_quality_df = custom_join_left(
        all_quality_df, me_annots.mz_quality[
            ['mz', 'msm', 'rawOffSampleProb', 'fdr',
             'quality_msm_off_fdr' # todo: add formula? this will go into .var
             ]], keycol='mz')

    print("(Moran's metrics are only computed for the m/z that passed the a) filters)")
    copy_mz = adata_mz[:, ok_part1].copy()   # only compute moran for the ok_part1 features
    try:
        if copy_mz.uns['moranI'].shape[0] != len(ok_part1):
            del copy_mz.uns['moranI']
    except Exception as e:
        pass

    if not "moranI" in list(copy_mz.uns.keys()):
        n_perms = 50
        if copy_mz.obs.shape[0] > 100000:
            n_perms = 3
        copy_mz = compute_moran_i(copy_mz, n_neighs=20, n_perms=n_perms)
    else:
        print(f"will use the existing .uns['moranI'] with dimension: {copy_mz.uns['moranI'].shape}")
        assert copy_mz.uns['moranI'].shape[0] == len(ok_part1), "Fatal: existing moranI dimensions do not match with length of m/z list passing a) filters. Killed! (solution: delete .uns['moranI'] dataframe and try again)"
    moran_mz_df = copy_mz.uns['moranI'].copy()
    moran_mz_df['mz'] = moran_mz_df.index.astype(str).tolist()
    all_quality_df = pd.merge(all_quality_df,
                              moran_mz_df[['mz', 'I', 'pval_sim_fdr_bh']],
                              how="left", left_on="mz", right_on="mz")
    keep_part2 = all_quality_df.loc[
        (all_quality_df['I'] >= moran_I_cutoff) & (
                all_quality_df['pval_sim_fdr_bh'] <= moran_padj_cutoff),
        "mz"].astype(str).tolist()

    text_notice = f"valid $m/z$: n = {len(keep_part2)}"
    plot__moran_info(copy_mz.uns['moranI'].copy(), plot_folder, SAMPLE,
                     text_notice)

    all_quality_df['quality_autocorr'] = 'bad'
    all_quality_df.loc[all_quality_df['mz'].astype(str).isin(keep_part2),
     'quality_autocorr'] = "good"

    return all_quality_df, keep_part2, moran_mz_df


if __name__ == "__main__":
    """
    bad or good quality of m/z features
    """
    # define here my files and dirs
    mydata_dir = "../../data"
    my_csv_table_path = "../small-test-parameters.tsv"  # " ../small-test-parameters.tsv"
    plot_parent_folder = "../../figures/"

    force_for_these = []    # keep empty if no sample to force   # e.g.  , "TD-635-1-A", "TD-635-1-B"   "TD-hemi-635-1-A", "TD-hemi-635-1-B"

    # start
    datasets_parameters_df = pd.read_csv(my_csv_table_path, header=0,
                      index_col=None, sep='\t')

    print("\nMarks bad quality m/z across specified datasets (.var property)")

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

        annotations_to_use = row["annotations_to_use"].split(" ")
        anno_dict = provide_annotation_dict(SAMPLE, annotations_to_use)
        me_annots = zd.MetabolomeAnnotations(
            name=SAMPLE, data_dir=mydata_dir,
            subfolder=subfolder, annotation_files=anno_dict)

        adata_mz = ad.read_h5ad(filename=os.path.join(
            mydata_dir, subfolder, f'{SAMPLE}.h5ad'))
        try:
            columns_to_reset = set(adata_mz.var.columns.tolist()) - set(["mz"])
            adata_mz.var = adata_mz.var.drop(columns=list(columns_to_reset))
        except KeyError:
            pass
        print("Loaded object", adata_mz.uns["history"])

        if (SAMPLE not in force_for_these) and ('features_n_after_mz_clean_2' in list(adata_mz.uns['history'].keys())):
            print("Nothing to do for this sample.")
        else:
            me_annots.load()  # load annotations
            me_annots, adata_mz = match_adata_mz_and_annots(me_annots,
                                                            adata_mz)
            if (normalize_yes_no.lower() in ['y', 'yes']) and (
                    not "is_X_TIC_normalized" in adata_mz.uns.keys()):
                adata_mz.layers['raw_X'] = adata_mz.X.copy()
                adata_mz.X = csr_matrix(tic_normalization(
                    adata_mz.X.copy(), adata_mz.obs['obs'].tolist(),
                    adata_mz.var['mz'].tolist()))
                adata_mz.uns["is_X_TIC_normalized"] = True
                print("Updated .X with TIC normalized matrix (raw stored in .layer['raw_X']).")

            me_annots, ok_mz_keep_part1 = part_a_filters(
                me_annots, msm_cutoff, offSample_cutoff, fdr_cutoff , SAMPLE, plot_folder)
            adata_mz.uns["history"]['features_n_after_mz_clean_1'] = len(ok_mz_keep_part1)

            all_quality_df, keep_part2, moran_mz_df = part_b_filters_with_autocorr(adata_mz, ok_mz_keep_part1,
                                                                                   me_annots,
                                                                                   moran_I_cutoff,
                                                                                   moran_padj_cutoff)
            adata_mz.uns['moranI'] = moran_mz_df
            adata_mz.uns["history"]['features_n_after_mz_clean_2'] = len(keep_part2)
            all_quality_df['final_quality'] = "bad"
            all_quality_df.loc[(all_quality_df['quality_msm_off_fdr'] == "good") & (
                all_quality_df['quality_autocorr'] == "good"), 'final_quality'
            ] = "good"
            okok = all_quality_df.loc[all_quality_df['final_quality']=="good",:].shape[0]
            if np.all(np.array(adata_mz.var.index.tolist()) == np.array(
                    all_quality_df['mz'].tolist())):
                all_quality_df.index = all_quality_df['mz'].astype(
                    str).tolist()
                all_quality_df['mz'] = all_quality_df['mz'].astype(
                    float).tolist()
                adata_mz.var = all_quality_df
                print(adata_mz)
                print(adata_mz.uns['history'])
            else:
                raise AssertionError("Fatal: could not assembly final object. Aborted.")
            answer = input(
                f"A total of {okok} m/z have good final_quality (pass both a) and b) filters)."
                f"\nDo you agree to continue?")
            if not answer.lower() in ["y", "yes"]:
                print("going to the next, nothing saved")
                continue
            else:
                print("Usable m/z elements (n of features):", okok)
                adata_mz.write_h5ad(filename=os.path.join(
                   mydata_dir, subfolder, f"{SAMPLE}.h5ad"))
                print("saved object")

