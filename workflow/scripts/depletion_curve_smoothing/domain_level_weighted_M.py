"""

"""

import sys
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
from util.calculate_weight_averaged_M import calculate_weight_averaged_M
from util.protein_domain_functions import assign_protein_domain


def main(args):
    initial_timepoint = args.initial_timepoint

    M_values = pd.read_csv(args.M_value_file, index_col=[0, 1, 2, 3], header=[0, 1])
    annotated_insertions = pd.read_csv(
        args.annotated_insertion_file, index_col=[0, 1, 2, 3], header=[0]
    )
    Confidence_score = pd.read_csv(
        args.Confidence_score_file, index_col=[0, 1, 2, 3], header=[0]
    )
    domain = pd.read_csv(args.domain_file, header=[0], sep="\t")

    annotated_insertions[["domain_id", "domain_residues"]] = annotated_insertions.apply(
        lambda row: assign_protein_domain(row, domain), axis=1, result_type="expand"
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

    domain_level_weighted_M = calculate_domain_level_weighted_M(
        in_gene_with_M_index,
        annotated_insertions,
        M_values,
        Confidence_score,
        domain,
        transformed_confidence_score,
        initial_timepoint,
    )
    domain_level_weighted_M.rename_axis(
        ["Systematic ID", "domain_id", "domain_residues"], axis=0, inplace=True
    )
    domain_level_weighted_M.xs("Weighted M", axis=1).to_csv(
        args.output_M, header=True, index=True, float_format="%.3f"
    )
    domain_level_weighted_M.xs("Statistic", axis=1).to_csv(
        args.output_statistic, header=True, index=True, float_format="%.3f"
    )


def calculate_domain_level_weighted_M(
    in_gene_with_M_index,
    annotated_insertions,
    M_values,
    Confidence_score,
    domain,
    transformed_confidence_score,
    initial_timepoint,
):
    domain_level_weighted_M = pd.DataFrame()

    for domain_index, domain_df in annotated_insertions.loc[
        in_gene_with_M_index
    ].groupby(["Systematic ID", "domain_id", "domain_residues"]):
        DWM = calculate_weight_averaged_M(
            domain_df.index, M_values, transformed_confidence_score, initial_timepoint
        )

        DWM_series = pd.Series(DWM, name=domain_index)

        domain_level_weighted_M = pd.concat(
            [domain_level_weighted_M, DWM_series.to_frame().T]
        )

    return domain_level_weighted_M


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
        "-d",
        "--domain-file",
        dest="domain_file",
        required=True,
        type=Path,
        help="File of domain",
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
