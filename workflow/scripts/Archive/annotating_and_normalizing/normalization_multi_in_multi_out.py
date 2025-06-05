"""
This script is to normalize the insertion reads by median based normalization.
"""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from utils.normalization_functions import normalizaiton_factor_report_multisample, calculte_MA, MA_plot


def parse_args():
    # add arguments
    parser = argparse.ArgumentParser(
        description="Normalize the insertion reads by median based normalization."
    )
    # input files
    parser.add_argument(
        "-i",
        "--input",
        dest="input",
        nargs="+",
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
        "-a",
        "--annotation",
        dest="annotation",
        type=Path,
        help="Insertion annotation",
    )
    # parser.add_argument(
    #     "-nf",
    #     "--normalization-factor-report",
    #     dest="normalization_factor_report",
    #     type=Path,
    #     help="normalization factor report",
    # )

    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        type=Path,
        help="output M values for normalized insertion reads file",
    )

    # parse arguments
    return parser.parse_args()


def main():

    args = parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "normalized_values").mkdir(parents=True, exist_ok=True)
    (args.output / "M_values").mkdir(parents=True, exist_ok=True)
    (args.output / "A_values").mkdir(parents=True, exist_ok=True)

    init_timepoint = args.init_timepoint

    # read all insertion reads files
    samples = {
        f.name.split(".")[0]: pd.read_csv(
            f,
            header=0,
            index_col=[0, 1, 2, 3]
        ) for f in args.input
    }

    reads_before_normalization = pd.concat(samples, axis=1)

    # normalize the insertion reads by median based normalization

    annotations = pd.read_csv(
        args.annotation, header=0, index_col=[0, 1, 2, 3])
    # normalization_report.to_csv(args.normalization_factor_report)
    intergenic_insertions_filtered = annotations[(annotations["Type"]
                                                 == "Intergenic region") & (annotations["Distance_to_region_start"] > 500) & (annotations["Distance_to_region_end"] > 500)].index

    median_values = reads_before_normalization[reads_before_normalization.index.isin(intergenic_insertions_filtered
                                                                                     )].median()
    min_median_values = median_values.min()

    reads_after_normalization = reads_before_normalization.mul(min_median_values).div(
        median_values
    )
    # MA values
    for sample in reads_after_normalization.columns.levels[0]:
        normalized_reads = reads_after_normalization[sample].dropna(
            axis=0).round(3)
        normalized_reads.to_csv(args.output / "normalized_values" / f"{sample}.normalized.csv",
                                header=True, index=True, float_format="%.3f")
        Ms, As = calculte_MA(normalized_reads, init_timepoint)
        Ms.to_csv(args.output / "M_values" / f"{sample}.M.csv",
                  header=True, index=True, float_format="%.3f")
        As.to_csv(args.output / "A_values" / f"{sample}.A.csv",
                  header=True, index=True, float_format="%.3f")


if __name__ == "__main__":

    main()
