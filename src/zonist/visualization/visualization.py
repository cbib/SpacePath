import os

from typing import List, Dict, Union
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import squidpy as sq



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



def plot__moran_info(df, plot_folder, SAMPLE, text_notice,
                      moran_I_cutoff, moran_padj_cutoff) -> None:
    df["log_var_sim"] = np.log(df["var_sim"].to_numpy())
    df['padj class'] = f" > {moran_padj_cutoff}"
    df.loc[df["pval_sim_fdr_bh"] <= moran_padj_cutoff,
        'padj class'] = f" <= {moran_padj_cutoff}"
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
