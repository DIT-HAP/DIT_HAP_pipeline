"""
This script is to compare the PBL and PBR reads and summarize the insertion number.
"""

import sys
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

from scipy.stats import pearsonr
from sklearn.metrics import mean_squared_error

import matplotlib as mpl
import matplotlib.pyplot as plt


def main(argv):
    # add arguments
    parser = argparse.ArgumentParser(description="Compare PBL and PBR reads.")
    parser.add_argument("-pbl", "--pbl_reads", help="PBL reads file", required=True)
    parser.add_argument("-pbr", "--pbr_reads", help="PBR reads file", required=True)
    parser.add_argument("-reads", help="Reads file", required=True)
    parser.add_argument("-pdf", help="Output pdf", required=True)

    # get arguments
    args = parser.parse_args()
    pbl_reads = Path(args.pbl_reads)
    pbr_reads = Path(args.pbr_reads)
    reads = Path(args.reads)
    pdf = Path(args.pdf)

    PBL_reads = pd.read_csv(pbl_reads, header=0, index_col=[0, 1, 2, 3, 4])
    PBR_reads = pd.read_csv(pbr_reads, header=0, index_col=[0, 1, 2, 3, 4])
    reads = pd.read_csv(reads, header=0, index_col=[0, 1, 2, 3, 4])
    timepoints = reads.columns.tolist()

    concat_reads = pd.concat(
        [PBL_reads, PBR_reads, reads],
        axis=1,
        keys=["PBL", "PBR", "Reads"],
        join="outer",
    )

    log2_concat_reads = concat_reads.applymap(lambda x: 0 if x == 0 else np.log2(x))

    plot_comparison(log2_concat_reads, timepoints, pdf)


def plot_comparison(log2_concat_reads, timepoints, pdf):
    rows = len(timepoints)
    max_PBL_PBR = (
        max(
            log2_concat_reads["PBL"].max(axis=1).max(),
            log2_concat_reads["PBR"].max(axis=1).max(),
        )
        + 1
    )

    fig, ax = plt.subplots(rows, 2, figsize=(14, 7 * rows))
    fig.tight_layout(h_pad=7, w_pad=7, pad=5)
    for row, tp in enumerate(timepoints):
        X = log2_concat_reads["PBL"][tp]
        Y = log2_concat_reads["PBR"][tp]
        ax[row, 0].scatter(
            X, Y, s=10, facecolor="none", edgecolor="black", alpha=0.5, rasterized=True
        )
        ax[row, 0].plot(
            [0, max_PBL_PBR],
            [0, max_PBL_PBR],
            color="red",
            linestyle="--",
            lw=2,
            alpha=0.5,
        )
        ax[row, 0].set_xlabel("PBL reads (log2)", fontsize=16)
        ax[row, 0].set_ylabel("PBR reads (log2)", fontsize=16)
        ax[row, 0].set_title(f"PBL vs PBR reads ({tp})", fontsize=16)
        ax[row, 0].tick_params(
            axis="both", which="major", labelsize=14, labelleft=True, labelbottom=True
        )

        # calculate the Pearson correlation coefficient
        PCC = pearsonr(X, Y)[0]
        RMSE = mean_squared_error(X, Y, squared=False)
        ax[row, 0].text(
            0.05,
            0.95,
            f"PCC = {PCC:.2f}\nRMSE = {RMSE:.2f}\nn = {log2_concat_reads.shape[0]}",
            transform=ax[row, 0].transAxes,
            fontsize=12,
            verticalalignment="top",
        )

        # plot the histogram
        max_reads = log2_concat_reads.max(axis=1).max()
        read_bins = np.linspace(0, max_reads, 50)
        ax[row, 1].hist(
            log2_concat_reads["Reads"][tp],
            bins=read_bins,
            edgecolor="gray",
            facecolor="lightgray",
        )
        ax[row, 1].set_xlabel("Reads (log2)", fontsize=16)
        ax[row, 1].set_ylabel("Frequency", fontsize=16)
        ax[row, 1].set_title(f"Reads ({tp})", fontsize=16)
        ax[row, 1].tick_params(
            axis="both", which="major", labelsize=14, labelleft=True, labelbottom=True
        )

    # rasterize the scatter plot
    fig.savefig(pdf, dpi=300, bbox_inches="tight")
    plt.close()


if __name__ == "__main__":
    main(sys.argv[1:])
