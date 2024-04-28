"""

"""

import argparse
from pathlib import Path
import pandas as pd
import numpy as np
from multiprocessing import Pool
from util.curve_fitting_functions import curve_fitting


def parse_args():
    # Parse arguments
    parser = argparse.ArgumentParser(
        description="Calculate domain level curve fitting.")
    parser.add_argument(
        "-C",
        "--cores",
        dest="cores",
        type=int,
        default=8,
        help="Number of cores to use",
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
        "--use-weighted-M",
        dest="use_weighted_M",
        type=int,
        default=1,
        help="Use weighted M or not",
    )
    parser.add_argument(
        "--use-LOWESS",
        dest="use_LOWESS",
        type=int,
        default=1,
        help="Use LOWESS or not",
    )
    parser.add_argument(
        "-f",
        "--fitting",
        dest="fitting",
        required=True,
        type=str,
        default="points5",
        help="fitting",
    )
    parser.add_argument(
        "-p",
        "--prediction",
        dest="prediction",
        required=True,
        type=str,
        default="points5",
        help="prediction",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        required=True,
        type=Path,
        help="Output file",
    )

    return parser.parse_args()


def main():

    args = parse_args()

    GMs = pd.read_csv(args.GM, index_col=[0, 1, 2, 3], header=0)

    annotated_insertions = pd.read_csv(
        args.annotated_insertion_file, index_col=[0, 1, 2, 3], header=[0]
    )

    in_gene_index = annotated_insertions[
        annotated_insertions["Type"] != "Intergenic region"
    ].index

    in_gene_with_M_index = GMs.index.intersection(in_gene_index)

    chunked_domain_insertions = chunk_domain_level_insertions(
        args.cores, in_gene_with_M_index, annotated_insertions)

    generations = GMs.reset_index(drop=True)[["Sample", "Timepoint", "G"]].drop_duplicates().groupby("Timepoint").apply(
        lambda x: np.average(x["G"])).sort_values()

    use_LOWESS = bool(int(args.use_LOWESS))
    use_weighted_M = bool(args.use_weighted_M)
    domain_level_fitting_params = [
        (generations, domain_df, GMs, args.fatol, use_weighted_M, use_LOWESS, args.fitting, args.prediction) for domain_df in chunked_domain_insertions]

    with Pool(args.cores) as pool:
        results = pool.starmap(
            domain_level_fitting, domain_level_fitting_params)

    domain_level_weighted_M = pd.concat(results, axis=0)
    domain_level_weighted_M.rename_axis(
        ["Systematic ID", "domain_id", "domain_residues"], axis=0, inplace=True
    )
    domain_level_weighted_M.to_csv(
        args.output, header=True, index=True, float_format="%.3f"
    )


def chunk_domain_level_insertions(cores, in_gene_with_M_index, annotated_insertions):

    annotations = annotated_insertions.loc[in_gene_with_M_index].copy()
    annotations = annotations.reset_index(drop=False).set_index(
        ["Systematic ID", "domain_id", "domain_residues"])
    uniqued_domains = annotations.index.unique()

    chunked_domains = np.array_split(uniqued_domains, cores)

    chunked_domain_insertions = [
        annotations.loc[domain].copy().reset_index(
            drop=False).set_index(["#Chr", "Coordinate", "Strand", "Target"])
        for domain in chunked_domains
    ]
    return chunked_domain_insertions


def domain_level_fitting(
    generation,
    sub_domain_insertions,
    GMs,
    fatol,
    using_weighted_M,
    LOWESS_smoothing,
    fitting,
    prediction
):
    domain_level_weighted_M = pd.DataFrame()

    for domain_index, domain_df in sub_domain_insertions.groupby(["Systematic ID", "domain_id", "domain_residues"]):

        insertions_in_current_level = GMs.loc[GMs.index.isin(
            domain_df.index)].copy()

        DWM = curve_fitting(
            domain_index, generation, insertions_in_current_level, fatol, useWeightedM=using_weighted_M, useLOWESS=LOWESS_smoothing, fitting=fitting, prediction=prediction)

        DWM_series = pd.Series(DWM, name=domain_index)

        domain_level_weighted_M = pd.concat(
            [domain_level_weighted_M, DWM_series.to_frame().T]
        )

    return domain_level_weighted_M


if __name__ == "__main__":

    main()
