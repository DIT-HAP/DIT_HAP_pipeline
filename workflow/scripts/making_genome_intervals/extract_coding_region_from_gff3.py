"""Extract coding region from gff3 file
This script is used to extract coding region from gff3 file.
Attention:
GFF is 1-based, but bed is 0-based.
Usage: python extract_coding_region_from_gff3.py -g <gff3 file> -p <peptide stats file> -o <output file>
"""
import sys
import argparse
import re
from pathlib import Path
import pandas as pd
import numpy as np
from utils.genome_utils import parse_gff_file, extract_coding_gene, extract_representative_coding_transcript


def parse_args():
    """
    Parse command line arguments
    """
    parser = argparse.ArgumentParser(
        description="Extract coding region from gff3 file"
    )
    parser.add_argument(
        "-g", "--gff", type=Path, required=True, help="Path to gff3 file"
    )
    parser.add_argument(
        "-p", "--peptide", type=Path, required=True, help="Path to peptide stats file")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Path to output file (Coding region bed)",
    )
    return parser.parse_args()


def main():

    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    # Parse gff file
    GFF = parse_gff_file(args.gff)

    # Extract coding gene
    region_of_representative_coding_transcript = extract_representative_coding_transcript(
        GFF, args.peptide)

    # Output
    region_of_representative_coding_transcript.to_csv(
        args.output, sep="\t", index=False, header=True)


if __name__ == "__main__":
    main()
