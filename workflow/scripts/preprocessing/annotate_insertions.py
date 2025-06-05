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

    # read file
    insertions = pd.read_csv(args.input, header=0, usecols=[0, 1, 2, 3], sep="\t").rename(columns={"Coordinate": "End"})
    insertions.insert(1, "Start", insertions["End"])
    print("*** For the insertions, the start and end are the same for bed processing")
    print(insertions.head())
    genome_region = pd.read_csv(args.genome_region, sep="\t", header=0).rename(columns={"#Chr": "Chr"})
    print("*** Loaded genome region file")
    print(genome_region.head())

    # transform the dataframe to bedtools
    insertions_bed = BedTool.from_dataframe(insertions)
    genome_region_bed = BedTool.from_dataframe(genome_region)

    # intersect the insertions with genome regions
    insertion_names = insertions.columns.tolist()
    print("*** Insertion names:")
    print(insertion_names)
    genome_region_names = add_suffix(
        genome_region.columns.tolist(), insertion_names, "_Interval"
    )
    print("*** Genome region names:")
    print(genome_region_names)
    insertions_with_annotation = insertions_bed.intersect(
        genome_region_bed, wa=True, wb=True)

    insertions_with_annotation = insertions_with_annotation.to_dataframe(
        names=insertion_names + genome_region_names)
    
    insertions_with_annotation.drop(columns=["Start"], inplace=True)
    insertions_with_annotation.rename(columns={"End": "Coordinate"}, inplace=True)
    print("*** Insertions with annotation:")

    # replace the cell of "." with np.nan
    print("After pybedtools, the NaN values are '.'")
    print(insertions_with_annotation.head())
    insertions_with_annotation.replace(
        r"^\.$", np.nan, inplace=True, regex=True)
    print("After replacing the '.' with np.nan, the NaN values are np.nan")
    print(insertions_with_annotation.head())
    print("*** Calculating the distance to the region start and end")
    insertions_with_annotation["Distance_to_region_start"] = (
        insertions_with_annotation["Coordinate"] -
        insertions_with_annotation["ParentalRegion_start"]
    )
    insertions_with_annotation["Distance_to_region_end"] = (
        insertions_with_annotation["ParentalRegion_end"] -
        insertions_with_annotation["Coordinate"]
    )
    insertions_with_annotation["Fraction_to_region_start"] = (
        insertions_with_annotation["Distance_to_region_start"]
        / insertions_with_annotation["ParentalRegion_length"]
    )
    insertions_with_annotation["Fraction_to_region_end"] = (
        insertions_with_annotation["Distance_to_region_end"]
        / insertions_with_annotation["ParentalRegion_length"]
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
    print("*** Calculated the distance to the start and stop codon")
    print(insertions_with_annotation.head())
    print("*** Calculating the residues affected")
    residue_stat = ["Residue_affected", "Residue_frame"]
    insertions_with_annotation[residue_stat] = insertions_with_annotation.apply(
        cal_residues_affected, residue_stat=residue_stat, axis=1
    )

    print("*** Assigning the insertion direction for insertions in the coding region")
    insertions_with_annotation[
        "Insertion_direction"
    ] = insertions_with_annotation.apply(assign_insertion_direction, axis=1)

    print("*** Dropping the duplicated insertions around boundary")
    insertions_with_annotation_no_duplicates = (
        insertions_with_annotation.groupby(["Chr", "Coordinate", "Strand"])
        .apply(drop_duplicated_insertions_around_boundary)
        .reset_index(drop=True)
    )

    print("*** Saving the annotated insertions")
    insertions_with_annotation_no_duplicates.to_csv(
        args.output, index=False, header=True, float_format="%.3f", sep="\t"
    )

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
        if (row["Feature"] == "CDS") and (row["Strand_Interval"] == "+"):
            CDS_base = CDS_base + \
                int(row["Coordinate"]) - int(row["Start_Interval"])
        elif (row["Feature"] == "CDS") and (row["Strand_Interval"] == "-"):
            CDS_base = CDS_base + int(row["End_Interval"] - row["Coordinate"])
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
    args = parser.parse_args()

    main(args)
