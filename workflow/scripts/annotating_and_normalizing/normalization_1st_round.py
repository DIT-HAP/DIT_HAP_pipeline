"""
This script is to normalize the insertion reads by median based normalization.
"""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from utils.normalization_functions import normalizaiton_factor_report, calculte_MA, MA_plot


def main(args):
    
    sample = args.input.name.split(".")[0]
    init_timepoint = args.init_timepoint

    # read all insertion reads files
    reads_before_normalization = pd.read_csv(
        args.input, header=0, index_col=[0, 1, 2, 3, 4]
    )

    # normalize the insertion reads by median based normalization
    reads_after_normalization = normalization(reads_before_normalization)
    normalization_report = normalizaiton_factor_report(
        reads_before_normalization, reads_after_normalization, sample
    )
    normalization_report.to_csv(args.normalization_factor_report)


    # MA values
    normalized_reads = reads_after_normalization.dropna(axis=0).round(3)
    Ms, As = calculte_MA(normalized_reads, init_timepoint)

    Ms.to_csv(args.output, header=True, index=True, float_format="%.3f")

    # MA plot
    MA_plot(Ms, As, args.MA_plot)

def normalization(reads_before_normalization):

    median_values = reads_before_normalization.median()
    min_median_values = median_values.min()

    reads_after_normalization = reads_before_normalization.mul(min_median_values).div(
        median_values
    )
    return reads_after_normalization

if __name__ == "__main__":

    # add arguments
    parser = argparse.ArgumentParser(
        description="Normalize the insertion reads by median based normalization."
    )
    # input files
    parser.add_argument(
        "-i",
        "--input",
        dest="input",
        type=Path,
        help="input insertion reads file",
    )
    parser.add_argument(
        "-t",
        "--init-timepoint",
        dest="init_timepoint",
        type=str,
        help="initial timepoint",
    )
    parser.add_argument(
        "-nf",
        "--normalization-factor-report",
        dest="normalization_factor_report",
        type=Path,
        help="normalization factor report",
    )
    parser.add_argument(
        "-MA", "--MA-plot", dest="MA_plot", type=Path, help="MA plot pdf file"
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        type=Path,
        help="output M values for normalized insertion reads file",
    )

    # parse arguments
    args = parser.parse_args()

    main(args)
