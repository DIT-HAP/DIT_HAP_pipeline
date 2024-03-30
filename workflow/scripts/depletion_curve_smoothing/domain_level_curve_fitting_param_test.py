"""

"""

import sys
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
import math
from scipy.optimize import curve_fit
from util.calculate_weight_averaged_M import calculate_weight_averaged_M
from util.curve_fitting_functions_test import curve_fitting
from util.protein_domain_functions import assign_protein_domain


def main(args):
    initial_timepoint = args.initial_timepoint

    generation = (
        pd.read_csv(args.generation_file, index_col=[
                    0], header=[0]).mean(axis=1).round(3)
    )

    M_values = pd.read_csv(args.M_value_file, index_col=[
                           0, 1, 2, 3], header=[0, 1])
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

    domain_level_weighted_M = curve_fitting_for_DL_DR(
        in_gene_with_M_index,
        annotated_insertions,
        M_values,
        Confidence_score,
        domain,
        transformed_confidence_score,
        initial_timepoint,
        generation,
        args.xatol,
        args.using_weighted_M,
    )
    domain_level_weighted_M.rename_axis(
        ["Systematic ID", "domain_id", "domain_residues"], axis=0, inplace=True
    )
    domain_level_weighted_M.to_csv(
        args.output_M, header=True, index=True, float_format="%.3f"
    )


def curve_fitting_for_DL_DR(
    in_gene_with_M_index,
    annotated_insertions,
    M_values,
    Confidence_score,
    domain,
    transformed_confidence_score,
    initial_timepoint,
    generation,
    xatol,
    using_weighted_M
):
    domain_level_weighted_M = pd.DataFrame()

    for domain_index, domain_df in annotated_insertions.loc[
        in_gene_with_M_index
    ].groupby(["Systematic ID", "domain_id", "domain_residues"]):

        DWM = curve_fitting(
            domain_df.index, M_values, generation, transformed_confidence_score, initial_timepoint, xatol, useWeightedM=using_weighted_M
        )

        DWM_series = pd.Series(DWM, name=domain_index)

        domain_level_weighted_M = pd.concat(
            [domain_level_weighted_M, DWM_series.to_frame().T]
        )

    return domain_level_weighted_M


if __name__ == "__main__":
    # Parse arguments
    parser = argparse.ArgumentParser(
        description="Calculate domain level curve fitting.")
    parser.add_argument(
        "-G",
        "--generation",
        dest="generation_file",
        required=True,
        type=Path,
        help="File of generation",
    )
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
        "-WM",
        "--using-weighted-M",
        dest="using_weighted_M",
        required=True,
        type=bool,
        help="Using weighted M or not",
    )
    parser.add_argument(
        "-xatol",
        "--xatol",
        dest="xatol",
        required=True,
        type=float,
        help="xatol",
    )
    parser.add_argument(
        "-o",
        "--output-M",
        dest="output_M",
        required=True,
        type=Path,
        help="Output file for M",
    )

    args = parser.parse_args()

    main(args)
