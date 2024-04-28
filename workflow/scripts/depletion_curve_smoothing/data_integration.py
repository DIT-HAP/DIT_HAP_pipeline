"""

"""

import sys
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
# from util.calculate_weight_averaged_M import calculate_weight_averaged_M
import matplotlib.pyplot as plt
from multiprocessing import Pool
import statsmodels.api as sm

from util.protein_domain_functions import assign_protein_domain
lowess = sm.nonparametric.lowess


def parse_args():
    # Parse arguments
    parser = argparse.ArgumentParser(description="Depletion curve smoothing")
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
        "-n",
        "--norm",
        dest="norm",
        required=True,
        type=Path,
        help="File of normalized reads",
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
        "-g",
        "--generation",
        dest="generation",
        required=True,
        type=Path,
        help="Generation file",
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
        "-o",
        "--output-GMs",
        dest="output_GMs",
        required=True,
        type=Path,
        help="Output GM file",
    )
    parser.add_argument(
        "-ao",
        "--annotation-output",
        dest="annotation_output",
        required=True,
        type=Path,
        help="Output annotation file with domains",
    )

    return parser.parse_args()


def main():

    args = parse_args()

    initial_timepoint = args.initial_timepoint

    M_values = pd.read_csv(args.M_value_file, index_col=[
        0, 1, 2, 3], header=[0, 1])
    Confidence_score = pd.read_csv(
        args.Confidence_score_file, index_col=[0, 1, 2, 3], header=[0]
    )

    transformed_confidence_score = (
        Confidence_score.rename_axis("Sample", axis=1)
        .stack(level="Sample")
        .to_frame("Confidence_score")
    )
    generation = pd.read_csv(
        args.generation, index_col=[0], header=[0]
    ).T.stack().rename("Generation").rename_axis(["Sample", "Timepoint"])

    normalized_reads = pd.read_csv(
        args.norm, index_col=[0, 1, 2, 3], header=[0, 1]
    ).stack(level="Sample").stack(level="Timepoint").to_frame("Normalized_reads")

    GMs = M_values.stack(level="Sample").stack(level="Timepoint").to_frame("M")
    insertion_sample_index = GMs.index.droplevel(level="Timepoint")
    tp_sample_index = GMs.index.droplevel(
        level=["#Chr", "Coordinate", "Strand", "Target"])

    formated_CS = transformed_confidence_score.loc[insertion_sample_index].values
    formated_generation = generation.loc[tp_sample_index].values
    formated_normalized_reads = normalized_reads.loc[GMs.index,
                                                     "Normalized_reads"]

    GMs["G"] = formated_generation
    GMs["CS"] = formated_CS
    GMs["normalized_reads"] = formated_normalized_reads
    GMs = GMs.reset_index(drop=False).set_index(
        ["#Chr", "Coordinate", "Strand", "Target"])

    GMs.to_csv(args.output_GMs, index=True,
               header=True, float_format="%.3f")

    annotated_insertions = pd.read_csv(
        args.annotated_insertion_file, index_col=[0, 1, 2, 3], header=[0]
    )
    domain = pd.read_csv(args.domain_file, header=[0], sep="\t")

    annotated_insertions[["domain_id", "domain_residues"]] = annotated_insertions.apply(
        lambda row: assign_protein_domain(row, domain), axis=1, result_type="expand"
    )
    annotated_insertions.to_csv(
        args.annotation_output, header=True, index=True)


if __name__ == "__main__":

    main()
