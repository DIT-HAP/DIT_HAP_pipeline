import argparse
from pathlib import Path
import pandas as pd
import numpy as np
from util.curve_fitting_functions import curve_fitting
from multiprocessing import Pool


def parse_args():
    # Parse arguments
    parser = argparse.ArgumentParser(description="Depletion curve smoothing")
    parser.add_argument(
        "-GM",
        "--GM",
        dest="GM",
        required=True,
        type=Path,
        help="File of M values",
    )
    parser.add_argument(
        "-C",
        "--cores",
        dest="cores",
        type=int,
        default=8,
        help="Number of cores to use",
    )
    parser.add_argument(
        "-fatol",
        "--fatol",
        dest="fatol",
        type=float,
        default=5e-3,
        help="fatol for curve fitting",
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


def site_level_curve_fitting_for_sub_df(sub_GMs_df, fatol):

    statistics = pd.DataFrame()
    for insertion, sub_df in sub_GMs_df.groupby(["#Chr", "Coordinate", "Strand", "Target"]):
        insertion_stat = curve_fitting(
            pd.MultiIndex.from_tuples([insertion]), sub_GMs_df, fatol, useWeightedM=True)
        insertion_stat.rename(insertion, inplace=True)
        insertion_stat = insertion_stat.to_frame().T
        statistics = pd.concat([statistics, insertion_stat], axis=0)

    return statistics


def main():

    args = parse_args()

    GMs = pd.read_csv(args.GM, index_col=[0, 1, 2, 3], header=0)

    pool = Pool(args.cores)

    GMs_index_unique = GMs.index.unique()
    GMs_index_unique_split = np.array_split(GMs_index_unique, args.cores)
    sub_GMs_with_fatol = [(GMs.loc[index].copy(), args.fatol)
                          for index in GMs_index_unique_split]

    results = pool.starmap(site_level_curve_fitting_for_sub_df,
                           sub_GMs_with_fatol)
    pool.close()

    statistics = pd.concat(results, axis=0)

    statistics.to_csv(args.output, header=True,
                      index=True, float_format="%.3f")


if __name__ == "__main__":
    main()
