"""
Insertion-level depletion analysis for non-replicates data
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from loguru import logger
from pydantic import BaseModel, Field, field_validator


# ======================== Configuration & Models ========================


class InsertionLevelDepletionAnalysisConfig(BaseModel):
    """Configuration for insertion-level depletion analysis."""

    counts_file: Path = Field(..., description="Path to the counts file")
    control_insertions_file: Path = Field(
        ..., description="Path to the control insertions file"
    )
    output_file: Path = Field(..., description="Path to the output file")

    @field_validator("counts_file", "control_insertions_file")
    def validate_input_exists(cls, v, field):
        if not v.exists():
            raise ValueError(f"Input file {v} does not exist")
        return v

    @field_validator("output_file")
    def validate_output_dir(cls, v):
        if not v.parent.exists():
            raise ValueError(f"Output directory {v.parent} does not exist")
        return v


# ======================== Setup Logging ========================


def setup_logging(log_level: str = "INFO") -> None:
    """Configure loguru for genomic annotation."""
    logger.remove()
    logger.add(
        sys.stdout,
        format="{time:HH:mm:ss} | {level: <8} | {message}",
        level=log_level,
        colorize=False,
    )


# ======================== Core Functions ========================


@logger.catch
def load_and_preprocess_data(
    counts_file: Path, control_insertions_file: Path
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    counts_df = pd.read_csv(
        counts_file, sep="\t", index_col=[0, 1, 2, 3], header=[0, 1]
    )
    control_insertions_df = pd.read_csv(
        control_insertions_file, sep="\t", index_col=[0, 1, 2, 3]
    )

    return counts_df, control_insertions_df


@logger.catch
def median_based_normalization(
    counts_df: pd.DataFrame, control_insertions_df: pd.DataFrame
) -> pd.DataFrame:
    median_values = counts_df.loc[control_insertions_df.index].median()
    min_median_values = median_values.min()
    counts_df_after_normalization = counts_df.mul(min_median_values).div(median_values)

    return counts_df_after_normalization


@logger.catch
def calculte_MA(
    counts_df_after_normalization: pd.DataFrame, init_timepoint: str
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    Ms = (
        -(counts_df_after_normalization + 1)
        .div((counts_df_after_normalization[init_timepoint] + 1), axis=0)
        .map(np.log2)
    )

    As = (counts_df_after_normalization + 1).mul(
        (counts_df_after_normalization[init_timepoint] + 1), axis=0
    ).map(np.log2) * 0.5

    return Ms, As


@logger.catch
def MA_plot(M: pd.DataFrame, A: pd.DataFrame, MA_pdf: Path):
    timepoint = M.columns.tolist()
    n_row = len(timepoint)

    fig, ax = plt.subplots(n_row, 1, figsize=(8, 8 * n_row), sharex=True, sharey=True)
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


# ======================== Main Function ========================


def main():
    """Main entry point for the script."""
    setup_logging()

    parser = argparse.ArgumentParser(description="Insertion-level depletion analysis")
    parser.add_argument(
        "-i", "--counts_file", type=Path, required=True, help="Path to the counts file"
    )
    parser.add_argument(
        "-c",
        "--control_insertions_file",
        type=Path,
        required=True,
        help="Path to the control insertions file",
    )
    parser.add_argument(
        "-o", "--output_file", type=Path, required=True, help="Path to the output file"
    )
    args = parser.parse_args()

    counts_df, control_insertions_df = load_and_preprocess_data(
        args.counts_file, args.control_insertions_file
    )

    counts_df_after_normalization = median_based_normalization(
        counts_df, control_insertions_df
    )

    Ms, As = calculte_MA(counts_df_after_normalization, "0")

    MA_plot(Ms, As, args.output_file)


if __name__ == "__main__":
    main()
