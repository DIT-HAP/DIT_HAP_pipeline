import argparse
import pybedtools
from pybedtools import BedTool
import numpy as np
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description="Annotate intergenic regions")
    parser.add_argument(
        "-ic", "--input_coding_gene", required=True, help="Input coding gene BED file"
    )
    parser.add_argument(
        "-iol",
        "--input_overlapped_gene_regions",
        required=True,
        help="Input overlapped gene regions BED file",
    )
    parser.add_argument(
        "-f", "--fai", required=True, help="FASTA .fai file for genome size"
    )
    parser.add_argument(
        "-o", "--output", required=True, help="Output intergenic regions BED file"
    )
    return parser.parse_args()


def main():

    args = parse_args()
    # Load the coding gene BED file
    genes = BedTool(args.input_coding_gene)
    overlapping_regions = BedTool(args.input_overlapped_gene_regions)

    # Create a genome file from the FASTA .fai
    genome = {}
    with open(args.fai, "r") as fai_file:
        for line in fai_file:
            chrom, length, _, _, _ = line.strip().split("\t")
            genome[chrom] = (0, int(length))

    # Find the complement (intergenic regions)
    intergenic = genes.sort().complement(g=genome)

    # non-overlapping coding regions
    non_overlapping_coding_regions = genes.subtract(overlapping_regions)

    non_overlapping_coding_regions_df = non_overlapping_coding_regions.to_dataframe(
        names=["#Chr", "Start", "End", "Systematic ID",
               "Type", "Strand", "Length"]
    )
    intergenic_df = intergenic.to_dataframe(
        names=["#Chr", "Start", "End", "Systematic ID",
               "Type", "Strand", "Length"]
    )
    intergenic_df["Type"] = "Intergenic region"

    # Merge the non-overlapping coding regions and intergenic regions
    merged_df = (
        pd.concat([non_overlapping_coding_regions_df,
                  intergenic_df], ignore_index=True)
        .sort_values(by=["#Chr", "Start", "End"])
        .reset_index(drop=True)
    )

    merged_df[["-1Strand", "-1Region"]] = merged_df.shift(1)[
        ["Strand", "Systematic ID"]
    ]
    merged_df[["+1Strand", "+1Region"]] = merged_df.shift(-1)[
        ["Strand", "Systematic ID"]
    ]

    merged_df.fillna(".", inplace=True)

    merged_df["Systematic ID"] = merged_df["-1Region"] + \
        "|" + merged_df["+1Region"]
    merged_df["Strand"] = merged_df["-1Strand"] + "|" + merged_df["+1Strand"]

    merged_df = merged_df.drop(
        columns=["-1Strand", "-1Region", "+1Strand", "+1Region"])

    intergenic_regions = merged_df[merged_df["Type"]
                                   == "Intergenic region"].copy()

    # modify the intergenic regions at the ends of the chromosomes
    for chr, length in genome.items():
        left_boundary = (intergenic_regions["#Chr"] == chr) & (
            intergenic_regions["Start"] == 0
        )
        right_boundary = (intergenic_regions["#Chr"] == chr) & (
            intergenic_regions["End"] == length[1]
        )

        if left_boundary.sum() != 0:
            left_boundary_intergenic_ID = intergenic_regions.loc[
                left_boundary, "Systematic ID"
            ].values[0]
            left_boundary_intergenic_strand = intergenic_regions.loc[
                left_boundary, "Strand"
            ].values[0]
            intergenic_regions.loc[left_boundary, "Systematic ID"] = (
                "." + "|" + left_boundary_intergenic_ID.split("|")[1]
            )
            intergenic_regions.loc[left_boundary, "Strand"] = (
                "." + "|" + left_boundary_intergenic_strand.split("|")[1]
            )

        if right_boundary.sum() != 0:
            right_boundary_intergenic_ID = intergenic_regions.loc[
                right_boundary, "Systematic ID"
            ].values[0]
            right_boundary_intergenic_strand = intergenic_regions.loc[
                right_boundary, "Strand"
            ].values[0]
            intergenic_regions.loc[right_boundary, "Systematic ID"] = (
                right_boundary_intergenic_ID.split("|")[0] + "|" + "."
            )
            intergenic_regions.loc[right_boundary, "Strand"] = (
                right_boundary_intergenic_strand.split("|")[0] + "|" + "."
            )

    intergenic_regions["Length"] = (
        intergenic_regions["End"] - intergenic_regions["Start"]
    )
    # save the intergenic regions
    intergenic_regions.to_csv(args.output, sep="\t", index=False, header=True)


if __name__ == "__main__":
    main()
