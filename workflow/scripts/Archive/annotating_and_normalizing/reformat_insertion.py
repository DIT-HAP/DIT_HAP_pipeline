"""

"""

import argparse
from pathlib import Path
import pandas as pd


def main(args):

    input_format = pd.read_csv(args.input, header=0, index_col=[0, 1, 2, 3, 4])

    output_format = input_format.droplevel(level="End", axis=0).rename_axis(
        ["#Chr", "Coordinate", "Strand", "Target"], axis=0
    )
    output_format.to_csv(args.output, header=True, index=True, float_format="%.3f")


if __name__ == "__main__":
    
    # add arguments
    parser = argparse.ArgumentParser(
        description="This script is to reformat the insertions."
    )

    parser.add_argument(
        "-i",
        "--input",
        dest="input",
        type=Path,
        help="input file",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        type=Path,
        help="output file",
    )

    args = parser.parse_args()

    main(args)
