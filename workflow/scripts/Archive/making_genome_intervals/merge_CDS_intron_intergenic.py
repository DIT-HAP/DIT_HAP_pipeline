"""
This script is to merge the CDS, intron, and intergenic regions into a single dataframe.
"""

import sys
import argparse
import pandas as pd


def parse_args():
    # add arguments
    parser = argparse.ArgumentParser(
        description="Merge the CDS, intron, and intergenic regions into a single dataframe."
    )
    parser.add_argument(
        "-ci",
        "--CDS-intron",
        dest="input_CDS_intron",
        metavar="FILE",
        help="Input CDS and intron regions.",
        required=True,
    )
    parser.add_argument(
        "-i",
        "--intergenic",
        dest="input_intergenic",
        metavar="FILE",
        help="Input intergenic regions.",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        metavar="FILE",
        help="Output merged regions.",
        required=True,
    )
    return parser.parse_args()


def main():

    args = parse_args()

    # read in
    CDS_intron_df = pd.read_csv(args.input_CDS_intron, sep="\t", header=0)
    intergenic_regions_df = pd.read_csv(
        args.input_intergenic, sep="\t", header=0
    )

    # Merge the CDS, intron, and intergenic regions into a single dataframe
    merged_df = pd.concat(
        [CDS_intron_df, intergenic_regions_df], ignore_index=True
    ).sort_values(by=["#Chr", "Start", "End"])

    merged_df.loc[(merged_df["Transcript"].isna()), "Start_Region"] = merged_df.loc[(
        merged_df["Transcript"].isna()), "Start"]
    merged_df.loc[(merged_df["Transcript"].isna()), "End_Region"] = merged_df.loc[(
        merged_df["Transcript"].isna()), "End"]
    merged_df.loc[(merged_df["Transcript"].isna()), "Length_Region"] = merged_df.loc[(
        merged_df["Transcript"].isna()), "Length"]

    merged_df = merged_df.astype(
        {"Start_Region": int, "End_Region": int, "Length_Region": int})

    # save
    merged_df.to_csv(args.output, sep="\t", header=True, index=False)


if __name__ == "__main__":
    main()
