"""
This script is to merge the read summary files for all samples.
"""

import argparse
import sys
from pathlib import Path
import pandas as pd


def main(args):
    # add arguments
    parser = argparse.ArgumentParser(
        description="Merge read summary files for all samples."
    )
    parser.add_argument(
        "-i", "--input", dest="input", help="Input directory of read summary files."
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        help="Output directory of merged read summary files.",
    )
    args = parser.parse_args()

    # get input and output directory
    input_dir = Path(args.input)
    output_dir = Path(args.output)

    # get all read summary files
    read_summary_files = input_dir.glob("*/*reads_summary.csv")

    df = pd.concat([pd.read_csv(f) for f in read_summary_files], axis=0)
    df.sort_values(by=df.columns[0], inplace=True)

    df.to_csv(output_dir / "reads_summary.csv", index=False)


if __name__ == "__main__":
    main(sys.argv[1:])
