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
from util.curve_fitting_functions import curve_fitting
from util.protein_domain_functions import assign_protein_domain


def main(args):

    GMs = pd.read_csv(args.GM, index_col=[0, 1, 2, 3], header=0)

    annotated_insertions = pd.read_csv(
        args.annotated_insertion_file, index_col=[0, 1, 2, 3], header=[0]
    )
    domain = pd.read_csv(args.domain_file, header=[0], sep="\t")

    DDR = pd.read_csv(args.DDR, header=[0])
    DDR = DDR.groupby("Systematic ID").apply(
        calcualte_DDR_ratio).reset_index(drop=True)

    domains_ratio_gt_point6 = pd.Index(
        DDR[(DDR["DR_ratio"] >= 0.6) & (DDR["Confidence score"] >= 50)][
            ["Systematic ID", "domain_id", "domain_residues"]
        ]
    )

    annotated_insertions[["domain_id", "domain_residues"]] = annotated_insertions.apply(
        lambda row: assign_protein_domain(row, domain), axis=1, result_type="expand"
    )

    in_ratio_gt_point6_domains_index = annotated_insertions[
        pd.Index(
            annotated_insertions[["Systematic ID",
                                  "domain_id", "domain_residues"]]
        ).isin(domains_ratio_gt_point6)
    ].index

    in_gene_with_M_index = GMs.index.intersection(
        in_ratio_gt_point6_domains_index)

    gene_level_weighted_M = curve_fitting_for_DL_DR(
        in_gene_with_M_index,
        annotated_insertions,
        GMs,
        args.fatol,
        using_weighted_M=True
    )
    gene_level_weighted_M.rename_axis("Systematic ID", axis=0, inplace=True)
    gene_level_weighted_M.to_csv(
        args.output, header=True, index=True, float_format="%.3f"
    )


def calcualte_DDR_ratio(DDR_sub_df):
    max_DDR = DDR_sub_df["DR"].max()
    DDR_sub_df["DR_ratio"] = round(DDR_sub_df["DR"] / max_DDR, 3)

    return DDR_sub_df


def curve_fitting_for_DL_DR(
    in_gene_with_M_index,
    annotated_insertions,
    GM_df,
    fatol,
    using_weighted_M
):
    gene_level_weighted_M = pd.DataFrame()

    for gene, gene_df in annotated_insertions.loc[in_gene_with_M_index].groupby(
        "Systematic ID"
    ):

        not_terminal_df = gene_df[gene_df["Distance_to_stop_codon"] >= 4].copy(
        )

        GWM = curve_fitting(
            not_terminal_df.index, GM_df, fatol, useWeightedM=using_weighted_M
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
        "-anno",
        "--annotated-insertion",
        dest="annotated_insertion_file",
        required=True,
        type=Path,
        help="File of annotated insertion",
    )
    parser.add_argument(
        "-GM",
        "--GM",
        dest="GM",
        required=True,
        type=Path,
        help="File of M values",
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
        "-fatol",
        "--fatol",
        dest="fatol",
        required=True,
        type=float,
        help="fatol",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        required=True,
        type=Path,
        help="Output file for M",
    )

    args = parser.parse_args()

    main(args)
