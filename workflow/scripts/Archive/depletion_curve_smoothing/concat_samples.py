"""
This script is to concatenate the samples from different datasets.
"""

import sys
import argparse
import pandas as pd
import numpy as np
from pathlib import Path


def parse_args():
    # Parse arguments
    parser = argparse.ArgumentParser(
        description="Concatenate samples from different datasets."
    )
    parser.add_argument(
        "-M",
        "--M-files",
        dest="M_files",
        required=True,
        nargs="*",
        type=Path,
        help="Files of M values",
    )
    parser.add_argument(
        "-A",
        "--A-files",
        dest="A_files",
        required=True,
        nargs="*",
        type=Path,
        help="Files of A values",
    )
    parser.add_argument(
        "-N",
        "--normalized-reads",
        dest="normalized_reads",
        required=True,
        nargs="*",
        type=Path,
        help="Files of normalized reads",
    )
    parser.add_argument(
        "-uTP",
        "--unused-timepoints",
        dest="unused_timepoints",
        nargs="*",
        default=[],
        help="Unused timepoints",
    )
    parser.add_argument(
        "-itp",
        "--initial-timepoint",
        dest="initial_timepoint",
        help="Initial timepoint",
    )

    parser.add_argument(
        "-oM",
        "--output-M",
        dest="output_M",
        required=True,
        type=Path,
        help="Output M file with concated values",
    )
    parser.add_argument(
        "-oA",
        "--output-A",
        dest="output_A",
        required=True,
        type=Path,
        help="Output A file with concated values",
    )
    parser.add_argument(
        "-oN",
        "--output-normalized",
        dest="output_normalized",
        required=True,
        type=Path,
        help="Output normalized file with concated values",
    )
    parser.add_argument(
        "-oCS",
        "--output-confidence-score",
        dest="output_CS",
        required=True,
        type=Path,
        help="Output confience score file with concated values",
    )

    return parser.parse_args()


def main():

    args = parse_args()

    unused_tps = args.unused_timepoints
    initial_tp = args.initial_timepoint

    M_dict = {
        Path(f)
        .name.split(".")[0]: pd.read_csv(f, index_col=[0, 1, 2, 3], header=0)
        .drop(columns=unused_tps, errors="ignore")
        for f in args.M_files
    }
    A_dict = {
        Path(f)
        .name.split(".")[0]: pd.read_csv(f, index_col=[0, 1, 2, 3], header=0)
        .drop(columns=unused_tps, errors="ignore")
        for f in args.A_files
    }
    N_dict = {
        Path(f)
        .name.split(".")[0]: pd.read_csv(f, index_col=[0, 1, 2, 3], header=0)
        .drop(columns=unused_tps, errors="ignore")
        for f in args.normalized_reads
    }

    concated_M = pd.concat(M_dict, axis=1).rename_axis(
        ["Sample", "Timepoint"], axis=1)
    concated_A = pd.concat(A_dict, axis=1).rename_axis(
        ["Sample", "Timepoint"], axis=1)
    concated_N = pd.concat(N_dict, axis=1).rename_axis(
        ["Sample", "Timepoint"], axis=1)
    Confidence_score = concated_N.xs(
        initial_tp, axis=1, level=1).apply(np.sqrt)

    concated_M.to_csv(args.output_M, index=True, float_format="%.3f")
    concated_A.to_csv(args.output_A, index=True, float_format="%.3f")
    concated_N.to_csv(args.output_normalized, index=True, float_format="%.3f")
    Confidence_score.to_csv(args.output_CS, index=True, float_format="%.3f")


if __name__ == "__main__":

    main()
