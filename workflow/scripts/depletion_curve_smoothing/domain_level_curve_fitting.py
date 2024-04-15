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

    annotated_insertions[["domain_id", "domain_residues"]] = annotated_insertions.apply(
        lambda row: assign_protein_domain(row, domain), axis=1, result_type="expand"
    )

    in_gene_index = annotated_insertions[
        annotated_insertions["Type"] != "Intergenic region"
    ].index

    in_gene_with_M_index = GMs.index.intersection(in_gene_index)

    domain_level_weighted_M = curve_fitting_for_DL_DR(
        in_gene_with_M_index,
        annotated_insertions,
        GMs,
        args.fatol,
        using_weighted_M=True,
    )
    domain_level_weighted_M.rename_axis(
        ["Systematic ID", "domain_id", "domain_residues"], axis=0, inplace=True
    )
    domain_level_weighted_M.to_csv(
        args.output, header=True, index=True, float_format="%.3f"
    )


def curve_fitting_for_DL_DR(
    in_gene_with_M_index,
    annotated_insertions,
    GM_df,
    fatol,
    using_weighted_M
):
    domain_level_weighted_M = pd.DataFrame()

    for domain_index, domain_df in annotated_insertions.loc[
        in_gene_with_M_index
    ].groupby(["Systematic ID", "domain_id", "domain_residues"]):

        DWM = curve_fitting(
            domain_df.index, GM_df, fatol, useWeightedM=using_weighted_M
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
        "-fatol",
        "--fatol",
        dest="fatol",
        required=True,
        type=float,
        help="fatol",
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
        "-o",
        "--output",
        dest="output",
        required=True,
        type=Path,
        help="Output file",
    )

    args = parser.parse_args()

    main(args)
