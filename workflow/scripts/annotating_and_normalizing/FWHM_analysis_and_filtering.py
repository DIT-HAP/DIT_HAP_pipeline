"""
This script is to do hard filtering based on the FWHM analysis.
"""

import sys
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from scipy import interpolate
from matplotlib.backends.backend_pdf import PdfPages as PdfPage
from utils.normalization_functions import calculte_MA, MA_plot


def main(args):
    init_timepoint = args.init_timepoint
    sample_name = args.input.name.split(".")[0]

    # read data
    normalized_reads = pd.read_csv(args.input, header=0, index_col=[0, 1, 2, 3, 4])
    log2_reads = np.log2(normalized_reads[init_timepoint] + 1)

    # FWHM plot
    FWHM, cutoff_x, peak_x = FWHM_calculate(log2_reads, sample_name, args.FWHM)
    # FWHM statistic
    FWHM_statistic = pd.DataFrame()
    FWHM_statistic.loc[sample_name, "FWHM"] = FWHM
    FWHM_statistic.loc[sample_name, "cutoff_x"] = cutoff_x
    FWHM_statistic.loc[sample_name, "peak_x"] = peak_x
    FWHM_statistic.rename_axis("Sample").to_csv(args.report, header=True, index=True, float_format="%.3f")

    # FWHM filtering
    filtered_reads = normalized_reads[log2_reads > cutoff_x].copy()

    # MA value calculation
    Ms, As = calculte_MA(filtered_reads, init_timepoint)

    # save files
    filtered_reads.to_csv(
        args.output_reads, header=True, index=True, float_format="%.3f"
    )
    Ms.to_csv(args.output_M, header=True, index=True, float_format="%.3f")
    As.to_csv(args.output_A, header=True, index=True, float_format="%.3f")

    # MA plot
    MA_plot(Ms, As, args.MA)

def FWHM_calculate(values, sample_name, pdf_file):
    kde_x, kde_y = sns.histplot(values, bins=50, kde=True).get_lines()[0].get_data()
    density = stats.gaussian_kde(values)
    density_y = density(kde_x)
    factor = np.max(kde_y) / np.max(density_y)

    x_smooth = np.linspace(kde_x.min(), kde_x.max(), 10000)
    y_smooth = density(x_smooth) * factor

    y_max = y_smooth.max()
    y_half = y_max / 2

    max_idx = np.argmax(y_smooth)
    left_idx = np.argmin(np.abs(y_smooth[:max_idx] - y_half))
    right_idx = np.argmin(np.abs(y_smooth[max_idx:] - y_half)) + max_idx

    peak_x = x_smooth[max_idx]
    left_x = x_smooth[left_idx]
    right_x = x_smooth[right_idx]

    FWHM = right_x - left_x

    cutoff_x = peak_x * 5 / 8

    plot_FWHM_distribution(
        values,
        x_smooth,
        y_smooth,
        peak_x,
        y_max,
        left_x,
        right_x,
        cutoff_x,
        sample_name,
        pdf_file,
    )

    return FWHM, cutoff_x, peak_x


def plot_FWHM_distribution(
    values, x, y, peak_x, peak_y, left_x, right_x, cutoff_x, sample_name, pdf_file
):
    half_y = peak_y / 2
    peak_reads = 2**peak_x - 1
    cutoff_reads = 2**cutoff_x - 1
    FWHM = right_x - left_x

    fig, ax = plt.subplots(figsize=(8, 6))

    ax.hist(values, bins=50, facecolor="lightgray", edgecolor="gray")
    ax.plot(x, y, color="black", linewidth=3)
    ax.scatter(
        [left_x, right_x],
        [half_y, half_y],
        facecolor="darkblue",
        edgecolor="darkred",
        s=50,
    )
    ax.plot(
        [left_x, right_x],
        [half_y, half_y],
        color="darkblue",
        linewidth=3,
        ls="--",
        label=f"FWHM: {FWHM:.2f}",
    )
    ax.axvline(
        x=peak_x, color="darkred", linewidth=3, ls="--", label=f"Peak: {peak_reads:.2f}"
    )
    ax.axvline(
        x=cutoff_x,
        color="darkgreen",
        linewidth=3,
        ls="--",
        label=f"Cutoff: {cutoff_reads:.2f}",
    )

    ax.set_xlabel("log2(reads+1)", fontsize=14)
    ax.set_ylabel("Frequency", fontsize=14)
    ax.tick_params(axis="both", which="major", labelsize=12)
    ax.set_title(sample_name, fontsize=16)
    ax.legend(loc="upper right", fontsize=12)

    fig.savefig(pdf_file, dpi=300, bbox_inches="tight")
    plt.close()


if __name__ == "__main__":
    # add arguments
    parser = argparse.ArgumentParser(
        description="This script is to do hard filtering based on the FWHM analysis."
    )
    parser.add_argument(
        "-i",
        "--input",
        dest="input",
        type=Path,
        help="input insertion read file",
    )
    parser.add_argument(
        "-itp",
        "--init-timepoint",
        dest="init_timepoint",
        type=str,
        help="initial timepoint",
    )
    parser.add_argument(
        "-f",
        "--FWHM",
        dest="FWHM",
        type=Path,
        help="FWHM plot file",
    )
    parser.add_argument(
        "-r",
        "--report",
        dest="report",
        type=Path,
        help="report file",
    )
    parser.add_argument(
        "-ma",
        "--MA",
        dest="MA",
        type=Path,
        help="MA plot after FWHM filtering",
    )
    parser.add_argument(
        "-or",
        "--output-Reads",
        dest="output_reads",
        type=Path,
        help="output read file",
    )
    parser.add_argument(
        "-om",
        "--output-M",
        dest="output_M",
        type=Path,
        help="output M value",
    )
    parser.add_argument(
        "-oa",
        "--output-A",
        dest="output_A",
        type=Path,
        help="output A value",
    )

    args = parser.parse_args()

    main(args)
