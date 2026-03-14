#########################################################################
# requires running first the two 'blend_..._to_tables.py' scripts
# --
# johaGL 2025
#########################################################################
import sys
import os
from pathlib import Path
from typing import Tuple
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl


###### these three lines are mandatory to comply with Nature journals ######
mpl.rcParams['font.family'] = "sans-serif"
mpl.rcParams['font.sans-serif'] = "Arial"
plt.rcParams.update({'font.size': 22})

FIG_DIR = Path("../../figs")

# ---------------------------------------------------------------------
# utilities
# ---------------------------------------------------------------------

def compute_marker_size_range(df: pd.DataFrame) -> Tuple[float, float]:
    return df["amount"].min() / 10, df["amount"].max() / 10


def load_csv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    return df.replace(0, np.nan)


# ---------------------------------------------------------------------
# dataframe transformations
# ---------------------------------------------------------------------

def transform_to_plot(df: pd.DataFrame, mode: str = "all") -> pd.DataFrame:
    cols = [c for c in df.columns if c.endswith(mode)]

    def melt_database(prefix: str, label: str) -> pd.DataFrame:
        col = next(c for c in cols if c.startswith(prefix))
        tmp = pd.melt(
            df,
            id_vars="sample_name",
            value_vars=col,
            value_name="amount",
            var_name="tmp",
        )
        tmp["database"] = label
        return tmp.drop(columns="tmp")

    frames = [
        melt_database("HMDB", "HMDB"),
        melt_database("Lipid", "LipidMaps"),
    ]

    kegg_cols = [c for c in cols if c.startswith("KEGG")]
    if kegg_cols:
        df3 = melt_database("KEGG", "KEGG")
        if not np.isnan(df3["amount"].to_numpy()).all():
            frames.append(df3)

    return pd.concat(frames)


def compute_good_samples(selection: str, df: pd.DataFrame) -> np.ndarray:
    tmp = df.loc[
        (df["ratio_id_mz"] <= 2.2)
        & (df["pathways_number_final"] >= 10)
        & (df["mz_mapped_to_pathways"] >= 10)
    ].copy()

    if selection == "private":
        tmp = tmp.loc[tmp["ratio_mz_pathway_count"] >= 1].copy()

    good_samples = tmp["sample_name"].unique()
    print("good tissues:\n", good_samples)
    return good_samples

# ---------------------------------------------------------------------
# plotting
# ---------------------------------------------------------------------

def bubble_plot(
    df_bg: pd.DataFrame,
    df_fg: pd.DataFrame,
    size_bg: Tuple[float, float],
    size_fg: Tuple[float, float],
    fig_size: Tuple[float, float],
    color_fg: str,
    out_path: Path,
) -> None:

    sns.set_context("paper", font_scale=1.5)
    sns.set_style("whitegrid")

    plt.figure(figsize=fig_size)

    ax = sns.scatterplot(
        df_bg,
        x="database",
        y="sample_name",
        size="amount",
        sizes=size_bg,
        color="silver",
        alpha=0.5,
        lw=0.4,
        edgecolors="black",
    )

    ax.grid(color="grey", linestyle="--", alpha=0.7, zorder=1)

    sns.scatterplot(
        df_fg,
        ax=ax,
        x="database",
        y="sample_name",
        size="amount",
        sizes=size_fg,
        color=color_fg,
        lw=0.4,
        edgecolors="black",
    ).legend(loc="upper left", title="count", bbox_to_anchor=(1, 1.03))

    plt.margins(x=1)
    sns.despine(top=True, right=True, left=True, bottom=True)

    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()


def pairgrid_plot(df: pd.DataFrame) -> None:

    sns.set_theme(style="whitegrid", font_scale=1.4)

    g = sns.PairGrid(
        df,
        x_vars=df.columns[1:5],
        y_vars="sample_name",
        height=10,
        aspect=0.25,
    )

    g.map(
        sns.stripplot,
        data=df,
        size=10,
        orient="h",
        jitter=False,
        hue="PType",
        palette="cividis",
        linewidth=0.7,
        edgecolor="black",
    )

    g.add_legend()

    titles = [
        "IDs per m/z",
        "m/z mapped to pathways",
        "pathways count",
        "m/z count per pathway",
    ]

    for ax, title in zip(g.axes.flat, titles):
        ax.set(title=title)
        ax.xaxis.grid(True)
        ax.yaxis.grid(True)

    sns.despine(left=True, bottom=True)
    plt.show()


# ---------------------------------------------------------------------
# main workflow
# ---------------------------------------------------------------------

def load_criteria_tables(selection: str) -> pd.DataFrame:

    df_kegg = load_csv(FIG_DIR / f"crit-KEGG-compounds-counts-{selection}.csv")
    df_mix = load_csv(FIG_DIR / f"crit-IDSmix-counts-{selection}.csv")

    df_kegg = df_kegg[
        ["sample_name", "ratio_keggid_mz", "mz_mapped_to_pathways", "pathways_number_final"]
    ]

    df_mix = df_mix[
        ["sample_name", "ratio_id_mz", "mz_mapped_to_pathways", "pathways_number_final"]
    ]

    df_kegg.columns = [
        "sample_name",
        "ratio_id_mz",
        "mz_mapped_to_pathways",
        "pathways_number_final",
    ]

    df_kegg["PType"] = "KEGG_pathways"
    df_mix["PType"] = "Wikipathways"

    df = pd.concat([df_kegg, df_mix])

    df["ratio_mz_pathway_count"] = (
        df["mz_mapped_to_pathways"] + 0.001
    ) / (df["pathways_number_final"] + 0.001)

    return df[["sample_name", "ratio_id_mz",
            "mz_mapped_to_pathways",  "pathways_number_final",
            "ratio_mz_pathway_count", "PType"]]


def main(selection: str, only_passing_qc: bool) -> None:

    fig_size = (5.7, 9) if selection == "public" else (5.3, 3.3)

    df = load_criteria_tables(selection)

    if only_passing_qc:
        good_samples = compute_good_samples(selection, df)
    else:
        good_samples = df['sample_name'].tolist()   # all if not QC filter

    df = df.loc[df["sample_name"].isin(good_samples)]

    df.index = df["sample_name"] + "&" + df["PType"]

    pairgrid_plot(df)

    # KEGG bubble plots
    df_kegg = load_csv(FIG_DIR / f"mapp-ids-pathways-KEGG-{selection}.csv")
    df_kegg = df_kegg.loc[df_kegg["sample_name"].isin(good_samples)]

    df_pl = transform_to_plot(df_kegg, "all")
    df_pl2 = transform_to_plot(df_kegg, "ids")
    df_pl3 = transform_to_plot(df_kegg, "pathways")

    bubble_plot(
        df_pl,
        df_pl2,
        compute_marker_size_range(df_pl),
        compute_marker_size_range(df_pl2),
        fig_size,
        "sandybrown",
        FIG_DIR / f"mapp-ids-compounds-KEGG-{selection}.pdf",
    )

    bubble_plot(
        df_pl,
        df_pl3,
        compute_marker_size_range(df_pl),
        compute_marker_size_range(df_pl3),
        fig_size,
        "darkcyan",
        FIG_DIR / f"mapp-ids-pathways-KEGG-{selection}.pdf",
    )

    # Wikipathways
    df_wp = load_csv(FIG_DIR / f"mapp-IDSmix-pathways-{selection}.csv")
    df_wp = df_wp.loc[df_wp["sample_name"].isin(good_samples)]

    df_pl = transform_to_plot(df_wp, "identifiers")
    df_pl2 = transform_to_plot(df_wp, "pathways")

    bubble_plot(
        df_pl,
        df_pl2,
        compute_marker_size_range(df_pl),
        compute_marker_size_range(df_pl2),
        fig_size,
        "darkcyan",
        FIG_DIR / f"mapp-ids-pws-WIKIPATHWAYS-{selection}.pdf",
    )

    print("done")


if __name__ == "__main__":
    main(selection="private",  # set as "private" as needed
         only_passing_qc=True)