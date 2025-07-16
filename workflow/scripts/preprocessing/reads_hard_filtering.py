"""
This script is to filter the insertion reads by hard filtering.
"""
import sys
import argparse
from pathlib import Path
import pandas as pd


def parse_args():
    # add arguments
    parser = argparse.ArgumentParser(
        description="Filter the insertion reads by hard filtering."
    )
    parser.add_argument(
        "-i",
        "--input",
        dest="input",
        type=Path,
        help="input insertion reads file",
    )
    parser.add_argument(
        "-itp",
        "--init-timepoint",
        dest="init_timepoint",
        type=str,
        help="initial timepoint",
    )
    parser.add_argument(
        "-c",
        "--cutoff",
        dest="cutoff",
        type=int,
        help="cutoff value",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        type=Path,
        help="output insertion reads file",
    )

    # parse arguments
    return parser.parse_args()


def main():

    args = parse_args()
    # read insertion reads file
    raw_reads = pd.read_csv(args.input, index_col=[0, 1, 2, 3], sep="\t", header=[0,1])
    shape_before_filtering = raw_reads.shape[0]
    print("*** Before filtering")
    print(raw_reads.head())
    print(raw_reads.columns)
    print("*** Init timepoint: ", args.init_timepoint)
    print("*** Cutoff: ", args.cutoff)
    # filtering
    filtered_reads = {}
    for sample, sample_reads in raw_reads.groupby(level="Sample", axis=1):
        filtered_reads[sample] = sample_reads[sample_reads[sample][args.init_timepoint] >= args.cutoff].copy()
    
    filtered_reads = pd.concat(list(filtered_reads.values()), axis=1)
    shape_after_filtering = filtered_reads.shape[0]
    print("*** After filtering")
    print(filtered_reads.head())
    print(filtered_reads.columns)
    # write to file
    filtered_reads.to_csv(args.output, sep="\t", header=True, index=True)

    print("### Hard filtering completed ###")
    print("*** Shape before filtering: ", shape_before_filtering)
    print("*** Shape after filtering: ", shape_after_filtering)
    print("*** Retention rate: ", round(shape_after_filtering / shape_before_filtering * 100, 4), "%")


if __name__ == "__main__":
    main()
