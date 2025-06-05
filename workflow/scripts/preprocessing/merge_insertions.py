"""
This script is to merge the insertions from the PBL and PBR reads.
Usage: python merge_insertions.py -il <inputPBL> -ir <inputPBR> -o <output>

Updated to handle the new format: Chr, Coordinate, +, - (single coordinate with strand counts)
"""
from pathlib import Path
import argparse
import pandas as pd
import numpy as np


def parse_args():
    """
    Parse the arguments.
    """
    parser = argparse.ArgumentParser(
        description="This script is to merge the insertions from the PBL and PBR reads. Input format: Chr, Coordinate, +, -"
    )
    parser.add_argument("-il", "--inputPBL",
                        help="The input file of PBL reads.", type=Path)
    parser.add_argument("-ir", "--inputPBR",
                        help="The input file of PBR reads.", type=Path)
    parser.add_argument(
        "-o", "--output", help="The output file of merged insertions.", type=Path)
    args = parser.parse_args()

    return args


def main():
    args = parse_args()

    # Merge the insertions
    PBLs, PBRs = read_PBL_and_PBR_insertion(args.inputPBL, args.inputPBR)
    PBL_PBRs = merge_PBL_and_PBR_insertion(PBLs, PBRs)
    PBL_PBRs.to_csv(args.output, sep='\t', index=True, header=True)
    
    print(f"Merged insertions written to: {args.output}")
    print(f"Total insertion sites: {len(PBL_PBRs)}")


def read_PBL_and_PBR_insertion(PBL_file, PBR_file):
    """
    Read the PBL and PBR insertions with the new format: Chr, Coordinate, +, -
    """
    if PBL_file.exists():
        print(f"Reading PBL file: {PBL_file}")
        PBLs = pd.read_csv(PBL_file, sep="\t", header=0)
        print(f"PBL insertions loaded: {len(PBLs)} sites")
        
        # Validate expected columns
        expected_cols = ["Chr", "Coordinate", "+", "-"]
        if not all(col in PBLs.columns for col in expected_cols):
            print(f"Warning: PBL file missing expected columns. Found: {list(PBLs.columns)}")
            print(f"Expected: {expected_cols}")
    else:
        print("The PBL file does not exist.")
        PBLs = pd.DataFrame(columns=["Chr", "Coordinate", "+", "-"])

    if PBR_file.exists():
        print(f"Reading PBR file: {PBR_file}")
        PBRs = pd.read_csv(PBR_file, sep="\t", header=0)
        print(f"PBR insertions loaded: {len(PBRs)} sites")
        
        # Validate expected columns
        expected_cols = ["Chr", "Coordinate", "+", "-"]
        if not all(col in PBRs.columns for col in expected_cols):
            print(f"Warning: PBR file missing expected columns. Found: {list(PBRs.columns)}")
            print(f"Expected: {expected_cols}")
    else:
        print("The PBR file does not exist.")
        PBRs = pd.DataFrame(columns=["Chr", "Coordinate", "+", "-"])

    return PBLs, PBRs


def merge_PBL_and_PBR_insertion(PBLs, PBRs):
    """
    Merge the PBL and PBR insertions with the new format: Chr, Coordinate, +, -
    """
    print("Merging PBL and PBR insertions...")
    
    # Merge on Chr and Coordinate (instead of Chr, Start, End)
    PBL_PBRs = pd.merge(
        PBLs,
        PBRs,
        how="outer",
        on=["Chr", "Coordinate"],
        suffixes=("_PBL", "_PBR"),
    )

    # fill na with 0
    PBL_PBRs.fillna(0, inplace=True)

    plusInsertion = PBL_PBRs[["Chr", "Coordinate", "-_PBL", "+_PBR"]].copy()
    plusInsertion["Strand"] = "+"
    plusInsertion.rename(
        columns={"-_PBL": "PBL", "+_PBR": "PBR"}, inplace=True)

    minusInsertion = PBL_PBRs[["Chr", "Coordinate", "+_PBL", "-_PBR"]].copy()
    minusInsertion["Strand"] = "-"
    minusInsertion.rename(
        columns={"+_PBL": "PBL", "-_PBR": "PBR"}, inplace=True)

    PBL_PBRs = (
        pd.concat([plusInsertion, minusInsertion], axis=0)
        .set_index(["Chr", "Coordinate", "Strand"])
        .astype(int)
        .sort_index()
    )
    PBL_PBRs["Reads"] = PBL_PBRs["PBL"] + PBL_PBRs["PBR"]

    return PBL_PBRs


if __name__ == "__main__":
    main()
