import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def normalizaiton_factor_report(
    reads_before_normalization, reads_after_normalization, sample
):
    median_before_normalization = reads_before_normalization.median()
    median_after_normalization = reads_after_normalization.median()
    normalization_report = median_before_normalization.to_frame(
        name="Median before normalization"
    ).astype(int)
    normalization_report[
        "Median after normalization"
    ] = median_after_normalization.round(2)
    normalization_report[
        "Read number before normalization"
    ] = reads_before_normalization.sum().round(2)
    normalization_report[
        "Read number after normalization"
    ] = reads_after_normalization.sum().round(2)
    normalization_report["Normalization factor"] = (
        median_before_normalization.min() / (median_before_normalization)
    ).round(2)
    normalization_report = normalization_report.rename_axis(
        ["Timepoint"], axis=0
    ).reset_index()
    normalization_report["Sample"] = sample
    normalization_report = normalization_report.set_index(
        ["Sample", "Timepoint"])

    return normalization_report


def normalizaiton_factor_report_multisample(
    reads_before_normalization, reads_after_normalization
):
    median_before_normalization = reads_before_normalization.median()
    median_after_normalization = reads_after_normalization.median()
    normalization_report = median_before_normalization.to_frame(
        name="Median before normalization"
    ).astype(int)
    normalization_report[
        "Median after normalization"
    ] = median_after_normalization.round(2)
    normalization_report[
        "Read number before normalization"
    ] = reads_before_normalization.sum().round(2)
    normalization_report[
        "Read number after normalization"
    ] = reads_after_normalization.sum().round(2)
    normalization_report["Normalization factor"] = (
        median_before_normalization.min() / (median_before_normalization)
    ).round(2)
    normalization_report = normalization_report.rename_axis(
        ["Sample", "Timepoint"], axis=0
    )

    return normalization_report


def calculte_MA(reads, init_timepoint):

    Ms = (
        -(reads + 1)
        .div((reads[init_timepoint] + 1), axis=0)
        .map(np.log2)
    )

    As = (reads + 1).mul(
        (reads[init_timepoint] + 1), axis=0
    ).map(np.log2) * 0.5

    return Ms, As


def MA_plot(M, A, MA_pdf):
    timepoint = M.columns.tolist()
    n_row = len(timepoint)

    fig, ax = plt.subplots(n_row, 1, figsize=(
        8, 8 * n_row), sharex=True, sharey=True)
    fig.tight_layout(h_pad=6, pad=5)
    for row, row_tp in enumerate(timepoint):
        M_values = M[row_tp]
        A_values = A[row_tp]

        ax[row].scatter(
            M_values,
            A_values,
            s=10,
            facecolor="none",
            edgecolor="black",
            alpha=0.5,
            rasterized=True,
        )
        ax[row].axvline(0, c="r", ls="--", lw=2, alpha=0.5)
        ax[row].set_xlabel("M value", fontsize=16)
        ax[row].set_ylabel("A value", fontsize=16)
        ax[row].set_title(f"MA plot - {row_tp}", fontsize=16)

    fig.savefig(MA_pdf, dpi=300, bbox_inches="tight")
    plt.close()
