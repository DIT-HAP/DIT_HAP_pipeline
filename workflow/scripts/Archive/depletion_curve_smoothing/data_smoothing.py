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
        "-C",
        "--cores",
        dest="cores",
        required=False,
        type=int,
        default=4,
        help="Number of cores",
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
        "-ob",
        "--output-before",
        dest="output_before",
        required=True,
        type=Path,
        help="Output file",
    )
    parser.add_argument(
        "-oa",
        "--output-after",
        dest="output_after",
        required=True,
        type=Path,
        help="Output file",
    )

    return parser.parse_args()


def LOWESS_smoothing(sub_GMs_df):

    # use LOESS for curve smoothing
    for insertion, sub_df in sub_GMs_df.groupby(["#Chr", "Coordinate", "Strand", "Target"]):

        sub_df_no_0h = sub_df[sub_df["G"] > 0].copy()
        smoothed = lowess(sub_df_no_0h["M"], sub_df_no_0h["G"],
                          frac=0.4, return_sorted=False, it=10, is_sorted=False)
        sub_GMs_df.loc[sub_df_no_0h.index, "M"] = smoothed

    return sub_GMs_df


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

    GMs = M_values.stack(level="Sample").stack(level="Timepoint").to_frame("M")
    insertion_sample_index = GMs.index.droplevel(level="Timepoint")
    tp_sample_index = GMs.index.droplevel(
        level=["#Chr", "Coordinate", "Strand", "Target"])

    formated_CS = transformed_confidence_score.loc[insertion_sample_index].values
    formated_generation = generation.loc[tp_sample_index].values

    GMs["G"] = formated_generation
    GMs["CS"] = formated_CS
    GMs = GMs.reset_index(drop=False).set_index(
        ["#Chr", "Coordinate", "Strand", "Target"])

    GMs.to_csv(args.output_before, index=True,
               header=True, float_format="%.3f")

    with Pool(args.cores) as pool:
        GMs_index_unique = GMs.index.unique()
        GMs_index_unique_split = np.array_split(GMs_index_unique, args.cores)
        sub_GMs = [GMs.loc[index].copy().set_index(["Sample", "Timepoint"], append=True)
                   for index in GMs_index_unique_split]

        results = pool.map(LOWESS_smoothing, sub_GMs)
        GMs = pd.concat(results, axis=0)

    GMs.to_csv(args.output_after, index=True, header=True, float_format="%.3f")


if __name__ == "__main__":

    main()
