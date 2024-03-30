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
    annotated_insertions = pd.read_csv(
        args.annotated_insertion_file, index_col=[0, 1, 2, 3], header=[0]
    )
    Confidence_score = pd.read_csv(
        args.Confidence_score_file, index_col=[0, 1, 2, 3], header=[0]
    )

    transformed_confidence_score = (
        Confidence_score.rename_axis("Sample", axis=1)
        .stack(level="Sample")
        .to_frame("Confidence_score")
    )

    in_gene_index = annotated_insertions[
        annotated_insertions["Type"] != "Intergenic region"
    ].index

    in_gene_with_M_index = M_values.index.intersection(in_gene_index)

    bin_level_weighted_M = calculate_bin_level_weighted_M(
        in_gene_with_M_index,
        annotated_insertions,
        M_values,
        Confidence_score,
        transformed_confidence_score,
        initial_timepoint,
    )
    bin_level_weighted_M.rename_axis(
        ["#Chr", "Coordinate", "Strand", "Target", "Systematic ID"],
        axis=0,
        inplace=True,
    )
    bin_level_weighted_M.xs("Weighted M", axis=1).to_csv(
        args.output_M, header=True, index=True, float_format="%.3f"
    )
    bin_level_weighted_M.xs("Statistic", axis=1).to_csv(
        args.output_statistic, header=True, index=True, float_format="%.3f"
    )


def calculate_bin_level_weighted_M(
    in_gene_with_M_index,
    annotated_insertions,
    M_values,
    Confidence_score,
    transformed_confidence_score,
    initial_timepoint,
):
    bin_level_weighted_M = pd.DataFrame()

    for gene_ID, gene_df in annotated_insertions.loc[in_gene_with_M_index].groupby(
        "Systematic ID"
    ):
        deasceding_confidence_score = (
            Confidence_score.loc[gene_df.index].sum(axis=1).sort_values(ascending=False)
        )

        while deasceding_confidence_score.shape[0] > 0:
            coordinates = deasceding_confidence_score.index.get_level_values(
                "Coordinate"
            )
            most_confident_coordinate = coordinates[0]
            most_confident_index = deasceding_confidence_score.index[0]
            most_confident_index_with_geneID = (*most_confident_index, gene_ID)

            bin = deasceding_confidence_score.loc[
                (coordinates >= (most_confident_coordinate - 15))
                & (coordinates <= (most_confident_coordinate + 15))
            ]

            BWM = calculate_weight_averaged_M(
                bin.index, M_values, transformed_confidence_score, initial_timepoint
            )

            BWM_series = pd.Series(BWM, name=most_confident_index_with_geneID)

            bin_level_weighted_M = pd.concat(
                [bin_level_weighted_M, BWM_series.to_frame().T]
            )

            deasceding_confidence_score = deasceding_confidence_score.drop(
                index=bin.index
            )

    return bin_level_weighted_M


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
        "-anno",
        "--annotated-insertion",
        dest="annotated_insertion_file",
        required=True,
        type=Path,
        help="File of annotated insertion",
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
