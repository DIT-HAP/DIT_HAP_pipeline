# %%
"""
This script concatenates the timepoints of the same subject into one file.
"""
from email.parser import Parser
import sys
import re
import argparse
import pandas as pd
from pathlib import Path
from Bio import SeqIO


def parse_args():
    # add arguments
    parser = argparse.ArgumentParser(
        description="Concatenates the timepoints of the same subject into one file."
    )
    parser.add_argument("-s", "--sample", help="Sample name.")
    parser.add_argument(
        "-in",
        "--input-files",
        dest="input_files",
        nargs="+",
        type=Path,
        help="input insertion reads file",
    )
    parser.add_argument(
        "-tp",
        "--timepoints",
        dest="timepoints",
        nargs="+",
        type=str,
        help="timepoints",
    )
    parser.add_argument(
        "-g", "--genome", type=Path, help="fasta file for the genome for adding target sequence"
    )
    parser.add_argument("-ol", "--outputPBL", type=Path,
                        help="Output PBL file.")
    parser.add_argument("-or", "--outputPBR", type=Path,
                        help="Output PBR file.")
    parser.add_argument("-o", "--outputReads", type=Path,
                        help="Output Reads file.")
    return parser.parse_args()


def main():

    args = parse_args()

    # find timepoints
    print(args.timepoints)
    print(args.input_files)

    ref = SeqIO.to_dict(SeqIO.parse(args.genome, "fasta"))
    concat_timepoints(args.input_files, args.timepoints, args.outputPBL,
                      args.outputPBR, args.outputReads, ref)


def concat_timepoints(input, timepoints, outputPBL, outputPBR, outputReads, ref):
    """
    Concatenates the timepoints of a subject into one file.
    """

    # concatenate timepoints
    concated_files = (
        pd.concat(
            [
                pd.read_csv(file, header=0, index_col=[0, 1, 2, 3])
                for file in input
            ],
            axis=1,
            keys=timepoints,
            join="outer",
        )
        .sort_index(level=0, axis=1, key=lambda x: x.str.lower())
        .sort_index(axis=0)
    )

    def get_target(chr, coordinate, ref): return str(
        ref[chr].seq[coordinate - 4: coordinate]
    )
    target = concated_files.index.to_frame().apply(
        lambda row: get_target(row["#Chr"], row["Start"], ref), axis=1
    )

    # insert target sequence
    concated_files = concated_files.set_index(
        target.rename("Target"), append=True)

    # save concatenated files
    concated_files.xs("PBL", level=1, axis=1).fillna(0).astype(int).to_csv(
        outputPBL, index=True
    )
    concated_files.xs("PBR", level=1, axis=1).fillna(0).astype(int).to_csv(
        outputPBR, index=True
    )
    concated_files.xs("Reads", level=1, axis=1).fillna(0).astype(int).to_csv(
        outputReads, index=True
    )


if __name__ == "__main__":
    main()
