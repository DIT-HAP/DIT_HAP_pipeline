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
    raw_reads = pd.read_csv(args.input, header=0, index_col=[0, 1, 2, 3, 4])
    # filtering
    filtered_reads = raw_reads[raw_reads[args.init_timepoint]
                               > args.cutoff].copy()
    # write to file
    filtered_reads.to_csv(args.output, sep=",", header=True, index=True)


if __name__ == "__main__":
    main()
