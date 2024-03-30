"""
This script is to dissect the PBL and PBR reads from the trimmed reads.
"""
import gzip
import sys
import argparse
import re
from pathlib import Path
from Bio import SeqIO
import pandas as pd


def main(argv):
    # add argument parser
    parser = argparse.ArgumentParser(
        description="This script is to dissect the PBL and PBR reads from the trimmed reads."
    )
    parser.add_argument("-i", "--input", help="The input file of trimmed reads.")
    parser.add_argument(
        "-ol", "--outputPBL", help="The output file of dissected PBL reads."
    )
    parser.add_argument(
        "-or", "--outputPBR", help="The output file of dissected PBR reads."
    )
    parser.add_argument("-s", "--summary", help="The summary file of dissected reads.")
    parser.add_argument("-n", "--cores", help="The number of cores to use.")
    args = parser.parse_args()

    # Get the input files
    inputFile = Path(args.input)
    outputPBLFile = Path(args.outputPBL)
    outputPBRFile = Path(args.outputPBR)
    summaryFile = Path(args.summary)

    if not summaryFile.exists():
        # Add headers to the summary file
        summaryHeader = [
            "Fastq File",
            "Total Reads",
            "Reads with PBL",
            "Reads with PBR",
            "Fraction of reads with PBL",
            "Fraction of reads with PBR",
            "Reads with PBL/PBR",
            "Fraction of reads with PBL/PBR",
            "Reads with TTAA",
            "Fraction of reads with TTAA",
        ]
        with open(summaryFile, "w") as summary_handle:
            summary_handle.write("\t".join(summaryHeader) + "\n")

    # Dissect the reads
    dissect_reads(inputFile, outputPBLFile, outputPBRFile, summaryFile)


def dissect_reads(fqIn, PBLout, PBRout, summaryInfo):
    # Define the output files
    PBLout_handle = gzip.open(PBLout, "wt")
    PBRout_handle = gzip.open(PBRout, "wt")

    # Define PBL and PBR seq
    PBL = "CATGCGTCAATTTTACGCAGACTATCTTTCTAGGG"
    PBR = "ACGCATGATTATCTTTAACGTACGTCACAATATGATTATCTTTCTAGGG"
    Insertion_pattern = re.compile(f"({PBL}|{PBR})([ATGC]{{4}})(\S+)")

    # Prepare the dictionary for storing reads information
    readsInfo = {
        "Fastq File": fqIn.name,
        "Total Reads": 0,
        "Reads with PBL": 0,
        "Reads with PBR": 0,
        "Fraction of reads with PBL": 0,
        "Fraction of reads with PBR": 0,
        "Reads with PBL/PBR": 0,
        "Fraction of reads with PBL/PBR": 0,
        "Reads with TTAA": 0,
        "Fraction of reads with TTAA": 0,
    }

    # analyze the input file
    with gzip.open(fqIn, "rt") as fqIn_handle:
        for record in SeqIO.parse(fqIn_handle, "fastq"):
            readsInfo["Total Reads"] += 1
            read_seq = str(record.seq)
            Insertion = re.search(Insertion_pattern, read_seq)

            if Insertion:
                PB_seq = Insertion.group(1)
                PB_target = Insertion.group(2)
                start, end = Insertion.span(2)[0], Insertion.span(3)[1]
                genome_seq_record = record[start:end]

                if (PB_seq == PBL) and (len(genome_seq_record) >= 30):
                    readsInfo["Reads with PBL"] += 1
                    SeqIO.write(genome_seq_record, PBLout_handle, "fastq")
                elif (PB_seq == PBR) and (len(genome_seq_record) >= 30):
                    readsInfo["Reads with PBR"] += 1
                    SeqIO.write(genome_seq_record, PBRout_handle, "fastq")

                if PB_target == "TTAA":
                    readsInfo["Reads with TTAA"] += 1

    # Close the output files
    PBLout_handle.close()
    PBRout_handle.close()

    # Calculate the fraction of reads with PBL and PBR
    readsInfo["Fraction of reads with PBL"] = (
        readsInfo["Reads with PBL"] / readsInfo["Total Reads"]
    )
    readsInfo["Fraction of reads with PBR"] = (
        readsInfo["Reads with PBR"] / readsInfo["Total Reads"]
    )
    readsInfo["Reads with PBL/PBR"] = (
        readsInfo["Reads with PBL"] + readsInfo["Reads with PBR"]
    )
    readsInfo["Fraction of reads with PBL/PBR"] = (
        readsInfo["Reads with PBL/PBR"] / readsInfo["Total Reads"]
    )
    readsInfo["Fraction of reads with TTAA"] = (
        readsInfo["Reads with TTAA"] / readsInfo["Total Reads"]
    )

    # Write the reads information to the summary file and keep two decimal places
    with summaryInfo.open("a") as summary_handle:
        summary_handle.write(
            f"{readsInfo['Fastq File']}\t"
            + f"{readsInfo['Total Reads']}\t"
            + f"{readsInfo['Reads with PBL']}\t"
            + f"{readsInfo['Reads with PBR']}\t"
            + f"{readsInfo['Fraction of reads with PBL']:.2%}\t"
            + f"{readsInfo['Fraction of reads with PBR']:.2%}\t"
            + f"{readsInfo['Reads with PBL/PBR']}\t"
            + f"{readsInfo['Fraction of reads with PBL/PBR']:.2%}\t"
            + f"{readsInfo['Reads with TTAA']}\t"
            + f"{readsInfo['Fraction of reads with TTAA']:.2%}\n"
        )


if __name__ == "__main__":
    main(sys.argv[1:])
