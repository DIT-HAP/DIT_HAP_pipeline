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
    DDR = pd.read_csv(args.DDR, header=[0])
    DDR = DDR.groupby("Systematic ID").apply(calcualte_DDR_ratio).reset_index(drop=True)
    domains_ratio_gt_point6 = pd.Index(
        DDR[(DDR["DR_ratio"] >= 0.6) & (DDR["Confidence score"] >= 50)][
            ["Systematic ID", "domain_id", "domain_residues"]
        ]
    )

    annotated_insertions[["domain_id", "domain_residues"]] = annotated_insertions.apply(
        lambda row: assign_protein_domain(row, domain), axis=1, result_type="expand"
    )

    transformed_confidence_score = (
        Confidence_score.rename_axis("Sample", axis=1)
        .stack(level="Sample")
        .to_frame("Confidence_score")
    )

    in_ratio_gt_point6_domains_index = annotated_insertions[
        pd.Index(
            annotated_insertions[["Systematic ID", "domain_id", "domain_residues"]]
        ).isin(domains_ratio_gt_point6)
    ].index

    in_gene_with_M_index = M_values.index.intersection(in_ratio_gt_point6_domains_index)

    gene_level_weighted_M = calculate_gene_level_weighted_M(
        in_gene_with_M_index,
        annotated_insertions,
        M_values,
        transformed_confidence_score,
        initial_timepoint,
    )
    gene_level_weighted_M.rename_axis("Systematic ID", axis=0, inplace=True)
    gene_level_weighted_M.xs("Weighted M", axis=1).to_csv(
        args.output_M, header=True, index=True, float_format="%.3f"
    )
    gene_level_weighted_M.xs("Statistic", axis=1).to_csv(
        args.output_statistic, header=True, index=True, float_format="%.3f"
    )


def calcualte_DDR_ratio(DDR_sub_df):
    max_DDR = DDR_sub_df["DR"].max()
    DDR_sub_df["DR_ratio"] = round(DDR_sub_df["DR"] / max_DDR, 3)

    return DDR_sub_df


def calculate_gene_level_weighted_M(
    in_gene_with_M_index,
    annotated_insertions,
    M_values,
    transformed_confidence_score,
    initial_timepoint,
):
    gene_level_weighted_M = pd.DataFrame()

    for gene, gene_df in annotated_insertions.loc[in_gene_with_M_index].groupby(
        "Systematic ID"
    ):

        not_terminal_df = gene_df[gene_df["Distance_to_stop_codon"]>=4].copy()

        GWM = calculate_weight_averaged_M(
            not_terminal_df.index, M_values, transformed_confidence_score, initial_timepoint
        )

        GWM_series = pd.Series(GWM, name=gene)

        gene_level_weighted_M = pd.concat(
            [gene_level_weighted_M, GWM_series.to_frame().T]
        )

    return gene_level_weighted_M


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
        "-ddr",
        "--DDR",
        dest="DDR",
        required=True,
        type=Path,
        help="File of domain-level DR",
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
