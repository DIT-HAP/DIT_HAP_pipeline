# %%
"""
This script extracts the transposon insertions from the filtered reads
"""

# Importing modules
import sys
import argparse
from multiprocessing import Pool
from pathlib import Path
import pandas as pd
import re


def main(argv):
    # add argument parser
    parser = argparse.ArgumentParser(
        description="This script is to extract the transposon insertions from the filtered reads."
    )
    parser.add_argument("-i", "--input", help="The input sam file of filtered reads.")
    parser.add_argument(
        "-o", "--output", help="The output bed file of transposon insertions."
    )
    parser.add_argument("-n", "--cores", help="The number of cores to use.")
    args = parser.parse_args()

    # Get the input files
    inputSAM = Path(args.input)
    outputBed = Path(args.output)
    cores = int(args.cores)

    # Extract the transposon insertions
    insertions = chunk_the_sam_file_concat_insertions(inputSAM, cores)

    insertions.to_csv(outputBed, sep="\t", index=True, header=True)


def chunk_the_sam_file_concat_insertions(sam_input, cores, chunk_size=10000000):
    """
    This function chunks the sam file
    """
    chunked_sams = pd.read_csv(
        sam_input,
        sep="\t",
        comment="@",
        names=["QNAME", "FLAG", "RNAME", "POS", "CIGAR"],
        usecols=(0, 1, 2, 3, 5),
        engine="c",
        chunksize=chunk_size,
    )

    # multi-processing
    pool = Pool(cores)
    # map the function to the chunked sam file
    chunked_insertions = pool.map(mapping_reads_to_insertions, chunked_sams)
    # close the pool
    pool.close()
    # concatenate the chunked insertion dataframes
    insertions = pd.concat(chunked_insertions)
    insertion_count = (
        insertions.value_counts()
        .to_frame("Count")
        .sort_index()
        .unstack("Strand")
        .droplevel(0, axis=1)
        .fillna(0)
        .astype(int)
    )

    return insertion_count


def mapping_reads_to_insertions(df_input):
    """
    This function maps the reads to the insertions
    """
    chunk_size = 100000000

    Strand = lambda flag: "+" if flag == 0 else "-"
    Target = lambda flag, Seq: Seq[:4] if flag == 0 else Seq[-4:]
    # For uniquely & perfectly mapped reads
    Coordinate1 = (
        lambda flag, posi, cigar: int(posi) + 3
        if flag == 0
        else int(posi) + int(cigar[:-1]) - 1
    )
    # For uniquely & perfectly/imperfectly mapped reads
    Coordinate2 = (
        lambda flag, posi, cigar: int(posi) + 3
        if flag == 0
        else int(posi)
        + sum(map(int, re.findall("(\d+)M", cigar)))
        + sum(map(int, re.findall("(\d+)D", cigar)))
        - 1
    )

    insertions = pd.DataFrame(columns=["#Chr", "Start", "End", "Strand"])
    insertions["#Chr"] = df_input["RNAME"]
    insertions["Strand"] = df_input["FLAG"].apply(Strand)
    # insertionBED["Target"] = sam_DF.apply(lambda row: Target(row["FLAG"],row["SEQ"]),axis=1)
    # set the end of TTAA as the coordinate of this TTAA site
    try:
        insertions["Start"] = df_input.apply(
            lambda row: Coordinate1(row["FLAG"], row["POS"], row["CIGAR"]), axis=1
        )
    # For indel
    except ValueError:
        insertions["Start"] = df_input.apply(
            lambda row: Coordinate2(row["FLAG"], row["POS"], row["CIGAR"]), axis=1
        )

    insertions["End"] = insertions["Start"]
    return insertions


if __name__ == "__main__":
    main(sys.argv[1:])
# %%
