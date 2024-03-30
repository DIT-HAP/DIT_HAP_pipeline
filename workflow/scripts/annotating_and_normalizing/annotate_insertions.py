"""
This script is to annnotate the insertions with genome regions.
"""

import sys
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
from pybedtools import BedTool


def main(args):
    # args.output.parent.mkdir(parents=True, exist_ok=True)
    # args.vv_sites.parent.mkdir(parents=True, exist_ok=True)

    # read file
    insertions = pd.read_csv(args.input, header=0, usecols=[0, 1, 2, 3, 4])
    genome_region = pd.read_csv(args.genome_region, sep="\t", header=0)

    # transform the dataframe to bedtools
    insertions_bed = BedTool.from_dataframe(insertions)
    genome_region_bed = BedTool.from_dataframe(genome_region)

    # intersect the insertions with genome regions
    insertion_names = insertions.columns.tolist()
    genome_region_names = add_suffix(
        genome_region.columns.tolist(), insertion_names, "_Interval"
    )
    insertions_with_annotation = insertions_bed.intersect(
        genome_region_bed, wa=True, wb=True
    ).to_dataframe(names=insertion_names + genome_region_names)

    insertions_with_annotation.rename(
        columns={"#Chr_Interval": "Chr_Interval"}, inplace=True
    )

    # replace the cell of "." with np.nan
    insertions_with_annotation.replace(
        r"^\.$", np.nan, inplace=True, regex=True)

    insertions_with_annotation["Distance_to_region_start"] = (
        insertions_with_annotation["Start"] -
        insertions_with_annotation["Start_Region"]
    )
    insertions_with_annotation["Distance_to_region_end"] = (
        insertions_with_annotation["End_Region"] -
        insertions_with_annotation["End"]
    )
    insertions_with_annotation["Fraction_to_region_start"] = (
        insertions_with_annotation["Distance_to_region_start"]
        / insertions_with_annotation["Length_Region"]
    )
    insertions_with_annotation["Fraction_to_region_end"] = (
        insertions_with_annotation["Distance_to_region_end"]
        / insertions_with_annotation["Length_Region"]
    )

    name_distance = [
        "Distance_to_start_codon",
        "Distance_to_stop_codon",
        "Fraction_to_start_codon",
        "Fraction_to_stop_codon",
    ]

    insertions_with_annotation[name_distance] = insertions_with_annotation.apply(
        calculate_distance_to_start_stop_codon, name_distance=name_distance, axis=1
    )

    residue_stat = ["Residue_affected", "Residue_frame"]
    insertions_with_annotation[residue_stat] = insertions_with_annotation.apply(
        cal_residues_affected, residue_stat=residue_stat, axis=1
    )
    insertions_with_annotation[
        "Insertion_direction"
    ] = insertions_with_annotation.apply(assign_insertion_direction, axis=1)

    insertions_with_annotation_no_duplicates = (
        insertions_with_annotation.groupby(["#Chr", "Start", "Strand"])
        .apply(drop_duplicated_insertions_around_boundary)
        .reset_index(drop=True)
    )

    insertions_with_annotation_no_duplicates.to_csv(
        args.output, index=False, header=True, float_format="%.3f"
    )

    # VV_sites = insertions_with_annotation_no_duplicates[
    #     insertions_with_annotation_no_duplicates["Essentiality"] == "V|V"
    # ][["#Chr", "Start", "End", "Strand", "Target"]]
    # VV_sites.to_csv(args.vv_sites, index=False, header=True)


def calculate_distance_to_start_stop_codon(row, name_distance):
    if row["Type"] != "Intergenic region":
        if row["Strand_Interval"] == "+":
            distance_values = [
                row["Distance_to_region_start"],
                row["Distance_to_region_end"],
                row["Fraction_to_region_start"],
                row["Fraction_to_region_end"],
            ]
        else:
            distance_values = [
                row["Distance_to_region_end"],
                row["Distance_to_region_start"],
                row["Fraction_to_region_end"],
                row["Fraction_to_region_start"],
            ]
    else:
        distance_values = [np.nan] * 4

    return pd.Series(distance_values, index=name_distance)


def cal_residues_affected(row, residue_stat):
    if row["Type"] == "Intergenic region":
        Residue_affected = np.nan
        Residue_frame = np.nan
    else:
        CDS_base = float(row["Accumulated_CDS_bases"])
        if (row["Type"] == "CDS") and (row["Strand_Interval"] == "+"):
            CDS_base = CDS_base + \
                int(row["Start"]) - int(row["Start_Interval"])
        elif (row["Type"] == "CDS") and (row["Strand_Interval"] == "-"):
            CDS_base = CDS_base + int(row["End_Interval"] - row["End"])
        Residue_affected = CDS_base // 3 + 1
        Residue_frame = CDS_base % 3
    return pd.Series([Residue_affected, Residue_frame], index=residue_stat)


def assign_insertion_direction(row):
    if row["Type"] == "Intergenic region":
        Insertion_direction = np.nan
    else:
        if row["Strand"] == row["Strand_Interval"]:
            Insertion_direction = "Forward"
        else:
            Insertion_direction = "Reverse"
    return Insertion_direction


def add_suffix(a_name_list, b_name_list, suffix):
    new_name_list = []
    for a_name in a_name_list:
        if a_name in b_name_list:
            new_name_list.append(a_name + suffix)
        else:
            new_name_list.append(a_name)
    return new_name_list


def drop_duplicated_insertions_around_boundary(sub_df):
    """
    This function is to drop duplicated insertions around boundary.
    """
    coding_region_index = sub_df[
        (sub_df["Distance_to_start_codon"] == 0)
        | (sub_df["Distance_to_stop_codon"] == 0)
    ].index
    sub_df = sub_df.drop(index=coding_region_index)
    if sub_df["Type"].unique().shape[0] == 1:
        return sub_df
    else:
        coding_region_index = sub_df[sub_df["Type"]
                                     != "Intergenic region"].index
        return sub_df.drop(index=coding_region_index)


if __name__ == "__main__":
    # add arguments
    parser = argparse.ArgumentParser(
        description="Annotate the insertions with genome regions."
    )
    parser.add_argument(
        "-i",
        "--input",
        dest="input",
        type=Path,
        help="input insertion reads files",
    )
    parser.add_argument(
        "-g",
        "--genome-region",
        dest="genome_region",
        type=Path,
        help="genome region file",
    )
    parser.add_argument("-o", "--output", dest="output",
                        type=Path, help="output file")
    # parser.add_argument(
    #     "-v",
    #     "--vv-sites",
    #     dest="vv_sites",
    #     type=Path,
    #     help="vv sites file",
    # )
    args = parser.parse_args()

    main(args)
