"""
This script is used to merge the genome regions (coding regions and intergenic regions) into a single dataframe.
"""

import argparse
import sys
import pandas as pd
from pybedtools import BedTool


def parse_args():
    # Parse the command line arguments
    parser = argparse.ArgumentParser(
        description="Merge genome regions (coding regions and intergenic regions) into a single dataframe"
    )
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
        "-iIGR",
        "--input_intergenic_regions",
        required=True,
        help="Input intergenic regions BED file",
    )
    parser.add_argument(
        "-o1",
        "--output_one",
        required=True,
        help="Output merged genome regions BED file (coding regions and intergenic regions)",
    )
    parser.add_argument(
        "-o2",
        "--output_two",
        required=True,
        help="Output merged genome regions BED file (non-overlapped coding regions, overlapped coding regions, intergenic regions)",
    )
    return parser.parse_args()


def main():

    args = parse_args()
    # Load the coding gene BED file
    genes = BedTool(args.input_coding_gene)
    overlapping_regions = BedTool(args.input_overlapped_gene_regions)
    intergenic_regions = BedTool(args.input_intergenic_regions)

    genes_df = pd.read_csv(args.input_coding_gene, sep="\t", header=0)
    overlapping_regions_df = pd.read_csv(
        args.input_overlapped_gene_regions, sep="\t", header=0
    )
    intergenic_regions_df = pd.read_csv(
        args.input_intergenic_regions, sep="\t", header=0
    )

    # non-overlapping coding regions
    non_overlapping_coding_regions = genes.subtract(overlapping_regions)
    non_overlapping_coding_regions_df = non_overlapping_coding_regions.to_dataframe(
        names=["#Chr", "Start", "End", "Systematic ID",
               "Type", "Strand", "Length"]
    )

    # Merge the coding regions and intergenic regions into a single dataframe
    merge_coding_regions_and_intergenic_regions(
        genes_df, intergenic_regions_df, args.output_one
    )

    # Merge the non-overlapping coding regions, overlapping coding regions, and intergenic regions into a single dataframe
    merge_non_overlapping_overlapping_intergenic_regions(
        non_overlapping_coding_regions_df,
        overlapping_regions_df,
        intergenic_regions_df,
        args.output_two,
    )


def merge_coding_regions_and_intergenic_regions(
    genes_df, intergenic_regions_df, output_one
):
    """
    Merge the coding regions and intergenic regions into a single dataframe.
    """
    # Merge the coding regions and intergenic regions into a single dataframe
    genome_regions_df = pd.concat(
        [genes_df, intergenic_regions_df], ignore_index=True
    ).sort_values(by=["#Chr", "Start", "End"])

    # save
    genome_regions_df.to_csv(output_one, sep="\t", header=True, index=False)


def merge_non_overlapping_overlapping_intergenic_regions(
    non_overlapping_coding_regions_df,
    overlapping_regions_df,
    intergenic_regions_df,
    output_two,
):
    """
    Merge the non-overlapping coding regions, overlapping coding regions, and intergenic regions into a single dataframe.
    """
    # Merge the non-overlapping coding regions, overlapping coding regions, and intergenic regions into a single dataframe
    genome_regions_df = pd.concat(
        [
            non_overlapping_coding_regions_df,
            overlapping_regions_df,
            intergenic_regions_df,
        ],
        ignore_index=True,
    ).sort_values(by=["#Chr", "Start", "End"])

    # save
    genome_regions_df.to_csv(output_two, sep="\t", header=True, index=False)


if __name__ == "__main__":
    main()
