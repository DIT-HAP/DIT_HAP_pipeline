"""

"""

import sys
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
from util.calculate_weight_averaged_M import calculate_weight_averaged_M


def main(args):
    initial_timepoint = args.initial_timepoint

    M_values = pd.read_csv(args.M_value_file, index_col=[0, 1, 2, 3], header=[0, 1])
    Confidence_score = pd.read_csv(
        args.Confidence_score_file, index_col=[0, 1, 2, 3], header=[0]
    )

    transformed_confidence_score = (
        Confidence_score.rename_axis("Sample", axis=1)
        .stack(level="Sample")
        .to_frame("Confidence_score")
    )

    site_level_weighted_Ms = M_values.apply(
        lambda row: calculate_weight_averaged_M(
            row.name, M_values, transformed_confidence_score, initial_timepoint
        ),
        axis=1,
    )

    site_level_weighted_Ms.xs("Weighted M", axis=1).to_csv(
        args.output_M, index=True, float_format="%.3f"
    )
    site_level_weighted_Ms.xs("Statistic", axis=1).to_csv(
        args.output_statistic, index=True, float_format="%.3f"
    )


if __name__ == "__main__":
    # Parse arguments
    parser = argparse.ArgumentParser(description="Calculate SDR.")
    parser.add_argument(
        "-M",
        "--M-values",
        dest="M_value_file",
        required=True,
        type=Path,
        help="File of M values",
    )
    parser.add_argument(
        "-CS",
        "--Confidence-score",
        dest="Confidence_score_file",
        required=True,
        type=Path,
        help="File of Confidence score",
    )
    parser.add_argument(
        "-itp",
        "--initial-timepoint",
        dest="initial_timepoint",
        required=True,
        type=str,
        help="Initial timepoint",
    )
    parser.add_argument(
        "-o",
        "--output-M",
        dest="output_M",
        required=True,
        type=Path,
        help="Output file for M",
    )
    parser.add_argument(
        "-os",
        "--output-statistic",
        dest="output_statistic",
        required=True,
        type=Path,
        help="Output file for statistic",
    )

    args = parser.parse_args()

    main(args)
