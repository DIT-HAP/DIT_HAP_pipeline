"""Extract CDS and intron from gff3 file
This script is used to extract CDS and intron from gff3 file.
Attention:
GFF is 1-based, but bed is 0-based.
"""
import sys
import argparse
import re
from pathlib import Path
import pandas as pd
import numpy as np
from utils.genome_utils import parse_gff_file, extract_coding_gene, extract_CDS_intron, cal_accumlated_CDS_bases


def parse_args():
    # Parse arguments
    parser = argparse.ArgumentParser(
        description="Extract CDS and intron from gff3 file"
    )
    parser.add_argument(
        "-g", "--gff", type=Path, required=True, help="Path to gff3 file")
    parser.add_argument(
        "-p", "--peptide", type=Path, required=True, help="Path to peptide stats file")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Path to output file (CDS and intron bed)",
    )
    return parser.parse_args()


def main():

    args = parse_args()

    # switch to path obj
    GFF = Path(args.gff)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    # Parse gff file
    GFF = parse_gff_file(GFF)

    # Extract coding gene
    representative_CDS_and_intron = extract_CDS_intron(GFF, args.peptide)

    # Output
    representative_CDS_and_intron.to_csv(
        output, sep="\t", index=False, header=True)

    accumulated_CDS_base = representative_CDS_and_intron.groupby(
        "Systematic ID").apply(cal_accumlated_CDS_bases).droplevel(0).sort_index()

    representative_CDS_and_intron_with_accumulated_CDS_base = representative_CDS_and_intron.join(
        accumulated_CDS_base["Accumulated_CDS_bases"])

    representative_CDS_and_intron_with_accumulated_CDS_base.to_csv(
        output, sep="\t", index=False, header=True)


if __name__ == "__main__":
    main()
