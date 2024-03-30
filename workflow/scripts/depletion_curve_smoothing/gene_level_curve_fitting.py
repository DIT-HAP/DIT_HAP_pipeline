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
from util.curve_fitting_functions import MG_curve, curve_fitting
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

    transformed_confidence_score = (
        Confidence_score.rename_axis("Sample", axis=1)
        .stack(level="Sample")
        .to_frame("Confidence_score")
    )

    in_ratio_gt_point6_domains_index = annotated_insertions[
        pd.Index(
            annotated_insertions[["Systematic ID",
                                  "domain_id", "domain_residues"]]
        ).isin(domains_ratio_gt_point6)
    ].index

    in_gene_with_M_index = M_values.index.intersection(
        in_ratio_gt_point6_domains_index)

    gene_level_weighted_M = curve_fitting_for_DL_DR(
        in_gene_with_M_index,
        annotated_insertions,
        M_values,
        Confidence_score,
        domain,
        transformed_confidence_score,
        initial_timepoint,
        generation,
        args.using_weighted_M,
        args.cf_method,
        args.xtol
    )
    gene_level_weighted_M.rename_axis("Systematic ID", axis=0, inplace=True)
    gene_level_weighted_M.to_csv(
        args.output_M, header=True, index=True, float_format="%.3f"
    )


def calcualte_DDR_ratio(DDR_sub_df):
    max_DDR = DDR_sub_df["DR"].max()
    DDR_sub_df["DR_ratio"] = round(DDR_sub_df["DR"] / max_DDR, 3)

    return DDR_sub_df


def curve_fitting_for_DL_DR(
    in_gene_with_M_index,
    annotated_insertions,
    M_values,
    Confidence_score,
    domain,
    transformed_confidence_score,
    initial_timepoint,
    generation,
    using_weighted_M,
    cf_method,
    xtol
):
    gene_level_weighted_M = pd.DataFrame()

    for gene, gene_df in annotated_insertions.loc[in_gene_with_M_index].groupby(
        "Systematic ID"
    ):

        not_terminal_df = gene_df[gene_df["Distance_to_stop_codon"] >= 4].copy(
        )

        GWM = curve_fitting(
            not_terminal_df.index, M_values, generation, transformed_confidence_score, initial_timepoint, useWeightedM=using_weighted_M, cf_method=cf_method, xtol=xtol
        )

        GWM_series = pd.Series(GWM, name=gene)

        gene_level_weighted_M = pd.concat(
            [gene_level_weighted_M, GWM_series.to_frame().T]
        )

    return gene_level_weighted_M


# def curve_fitting(
#     index, M_values, generation, transformaed_CS_values, initial_timepoint
# ):
#     GWM = calculate_weight_averaged_M(
#             index, M_values, transformaed_CS_values, initial_timepoint
#         )
#     Weighted_M = GWM["Weighted M"].values[1:]
#     generations = generation.values[1:]

#     if isinstance(index, tuple):
#         index = [index]
#     else:
#         index = index.tolist()

#     sub_Ms = M_values.loc[index]
#     sub_Ms = sub_Ms.stack(level="Sample").stack(level="Timepoint").to_frame("M")

#     generation_index = sub_Ms.index.get_level_values("Timepoint")
#     sub_generation = generation.loc[generation_index].values

#     CS_index = sub_Ms.index.droplevel(level="Timepoint")
#     sub_CS = transformaed_CS_values.loc[CS_index, "Confidence_score"].tolist()

#     sub_Ms["G"] = sub_generation
#     sub_Ms["Confidence_score"] = sub_CS

#     confidence_weights = 1/sub_Ms["Confidence_score"].values

#     init_guess = [0,0]

#     # popt, pcov = curve_fit(MG_curve, sub_Ms["G"].values, sub_Ms["M"].values, method='trf', p0=init_guess, maxfev=100000)
#     try:
#         popt, pcov, infodict, mesg, ier = curve_fit(MG_curve, generations, Weighted_M, method='dogbox', p0=init_guess, maxfev=100000, full_output=True)

#         DL, DR = popt
#     except RuntimeError:
#         DL, DR = np.nan, np.nan

#     # DL_pcov, DR_pcov = np.diag(pcov)
#     # linalg_cond_pcov = np.linalg.cond(pcov)


#     confidence_score = sub_Ms.xs(initial_timepoint, level="Timepoint")[
#         "Confidence_score"
#     ]

#     statistics = pd.Series()
#     individual_insertions = confidence_score.shape[0]
#     sum_confidence_score = confidence_score.sum()
#     statistics.loc["Individual insertions"] = individual_insertions
#     statistics.loc["Confidence score"] = sum_confidence_score
#     statistics.loc["DL"] = round(DL, 3)
#     statistics.loc["DR"] = round(DR, 3)

#     # statistics.loc["DL_pcov"] = round(DL_pcov, 3)
#     # statistics.loc["DR_pcov"] = round(DR_pcov, 3)
#     # statistics.loc["linalg_cond_pcov"] = round(linalg_cond_pcov, 3)

#     return statistics

if __name__ == "__main__":
    # Parse arguments
    parser = argparse.ArgumentParser(description="Calculate SDR.")
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
        "-WM",
        "--using-weighted-M",
        dest="using_weighted_M",
        required=True,
        type=bool,
        help="Using weighted M or not",
    )
    parser.add_argument(
        "-cfm",
        "--curve-fitting-method",
        dest="cf_method",
        required=True,
        type=str,
        help="Curve fitting method",
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
        "-xtol",
        "--xtol",
        dest="xtol",
        required=False,
        type=float,
        help="xtol for curve fitting",
    )

    args = parser.parse_args()

    main(args)
