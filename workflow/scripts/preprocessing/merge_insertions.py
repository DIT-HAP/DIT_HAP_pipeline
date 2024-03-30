"""
This script is to merge the insertions from the PBL and PBR reads.
Usage: python merge_insertions.py -il <inputPBL> -ir <inputPBR> -o <output>
"""
from pathlib import Path
import argparse
import sys
from matplotlib import path
import numpy as np
import pandas as pd


def parse_args():
    """
    Parse the arguments.
    """
    parser = argparse.ArgumentParser(
        description="This script is to merge the insertions from the PBL and PBR reads."
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
    PBL_PBRs.to_csv(args.output, index=True, header=True)


def read_PBL_and_PBR_insertion(PBL_file, PBR_file):
    """
    Read the PBL and PBR insertions.
    """
    if PBL_file.exists():
        PBLs = pd.read_csv(PBL_file, sep="\t", header=0)
    else:
        print("The PBL file does not exist.")
        PBLs = pd.DataFrame(columns=["#Chr", "Start", "End", "+", "-"])

    if PBR_file.exists():
        PBRs = pd.read_csv(PBR_file, sep="\t", header=0)
    else:
        print("The PBR file does not exist.")
        PBRs = pd.DataFrame(columns=["#Chr", "Start", "End", "+", "-"])

    return PBLs, PBRs


def merge_PBL_and_PBR_insertion(PBLs, PBRs):
    """
    Merge the PBL and PBR insertions.
    """
    PBL_PBRs = pd.merge(
        PBLs,
        PBRs,
        how="outer",
        on=["#Chr", "Start", "End"],
        suffixes=("_PBL", "_PBR"),
    ).replace(0, np.nan)

    plusInsertion = PBL_PBRs[["#Chr", "Start", "End", "-_PBL", "+_PBR"]].copy()
    plusInsertion.dropna(axis=0, how="all", subset=[
                         "-_PBL", "+_PBR"], inplace=True)
    plusInsertion["Strand"] = "+"
    plusInsertion.rename(
        columns={"-_PBL": "PBL", "+_PBR": "PBR"}, inplace=True)

    minusInsertion = PBL_PBRs[[
        "#Chr", "Start", "End", "+_PBL", "-_PBR"]].copy()
    minusInsertion.dropna(axis=0, how="all", subset=[
                          "+_PBL", "-_PBR"], inplace=True)
    minusInsertion["Strand"] = "-"
    minusInsertion.rename(
        columns={"+_PBL": "PBL", "-_PBR": "PBR"}, inplace=True)

    PBL_PBRs = (
        pd.concat([plusInsertion, minusInsertion], axis=0)
        .set_index(["#Chr", "Start", "End", "Strand"])
        .fillna(0)
        .astype(int)
        .sort_index()
    )
    PBL_PBRs["Reads"] = PBL_PBRs["PBL"] + PBL_PBRs["PBR"]

    return PBL_PBRs


if __name__ == "__main__":
    main()
