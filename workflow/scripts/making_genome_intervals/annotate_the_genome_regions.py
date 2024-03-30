"""
This sript is to annotate the genome regions for insertion annotations and normalization with VV insertions.
"""
import sys
import argparse
from pathlib import Path
import pandas as pd


def parse_args():
    # add arguments
    parser = argparse.ArgumentParser(
        description="Annotate the genome regions for insertion annotations and normalization with VV insertions."
    )
    parser.add_argument(
        "-i",
        "--input",
        dest="input",
        metavar="input",
        type=Path,
        help="input genome region bed file",
    )
    parser.add_argument(
        "-n",
        "--name",
        dest="name",
        metavar="name",
        type=Path,
        help="File storing the name of the coding genes",
    )
    parser.add_argument(
        "-e",
        "--essentiality",
        dest="essentiality",
        metavar="essentiality",
        type=Path,
        help="File storing the essentiality of the coding genes (excel file)",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        metavar="output",
        type=Path,
        help="output of the annotated genome region bed file",
    )
    return parser.parse_args()


def main():

    args = parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)

    genome_regions = pd.read_csv(args.input, sep="\t", header=0)
    gene_names = pd.read_csv(args.name, sep="\t", header=0)
    essentiality = pd.read_excel(args.essentiality, header=0)[
        ["Systematic ID", "Gene dispensability. This study"]
    ]

    gene_names_with_essentiality = update_systemaic_id_and_add_info(
        gene_names, essentiality
    )

    genome_regions = add_essentialiyt_and_name_for_the_regions(
        genome_regions, gene_names_with_essentiality
    )

    genome_regions.to_csv(
        args.output, sep="\t", index=False
    )


def update_systemaic_id_and_add_info(gene_names, essentiality):
    """
    Update the systematic id in the essentialiyt file.
    """

    # get the protein coding gene names and their synonyms
    coding_gene_names = gene_names[gene_names["gene_type"] == "protein coding gene"][
        ["gene_systematic_id", "gene_name"]
    ].copy()
    coding_gene_synonyms = gene_names[gene_names["gene_type"] == "protein coding gene"][
        ["gene_systematic_id", "synonyms"]
    ].copy()

    # set the synonyms as the index and their systematic id as the value
    coding_gene_synonyms = (
        coding_gene_synonyms.set_index("gene_systematic_id")["synonyms"]
        .str.split(",", expand=True)
        .stack()
        .droplevel(level=1)
        .rename("synonyms")
        .reset_index()
        .set_index("synonyms")["gene_systematic_id"]
        .to_dict()
    )

    # print the intersection of synonyms and gene_systematic_id
    print(
        "The intersection of synonyms and gene_systematic_id is (should be empty): ",
        set(coding_gene_synonyms.keys()).intersection(
            set(coding_gene_names["gene_systematic_id"])
        ),
    )

    # fill the empty "gene_name" with "gene_systematic_id"
    coding_gene_names["gene_name"] = coding_gene_names["gene_name"].fillna(
        coding_gene_names["gene_systematic_id"]
    )

    # replace the "Systematic ID" in the "synonyms" with "gene_systematic_id"
    essentiality["Systematic ID"] = essentiality["Systematic ID"].replace(
        coding_gene_synonyms
    )

    gene_names_with_essentiality = (
        pd.merge(
            coding_gene_names,
            essentiality,
            how="left",
            left_on="gene_systematic_id",
            right_on="Systematic ID",
        )
        .drop(columns=["Systematic ID"])
        .rename(
            columns={
                "gene_systematic_id": "Systematic ID",
                "gene_name": "Name",
                "Gene dispensability. This study": "Essentiality",
            }
        )
        .set_index("Systematic ID")
    )

    # fill the NA in the Essentiality column with "Unknown"
    gene_names_with_essentiality["Essentiality"] = gene_names_with_essentiality[
        "Essentiality"
    ].fillna("Unknown")

    # genes in the essentiality file but not in the name file
    print(
        "Genes in the essentiality file but not in the name file: ",
        set(essentiality["Systematic ID"]).difference(
            set(coding_gene_names["gene_systematic_id"])
        ),
    )

    return gene_names_with_essentiality


def add_info_based_on_systematic_id(SysID, info_column, gene_names_with_essentiality):
    """
    Add the info based on the systematic id.
    """
    if "|" in SysID:
        left_id, right_id = SysID.split("|")
        if (left_id != ".") and (right_id != "."):
            left_info = gene_names_with_essentiality.loc[left_id, info_column]
            right_info = gene_names_with_essentiality.loc[right_id, info_column]
        elif (left_id != ".") and (right_id == "."):
            left_info = gene_names_with_essentiality.loc[left_id, info_column]
            right_info = "."
        elif (left_id == ".") and (right_id != "."):
            left_info = "."
            right_info = gene_names_with_essentiality.loc[right_id, info_column]
        return left_info + "|" + right_info
    else:
        return gene_names_with_essentiality.loc[SysID, info_column]


def add_essentialiyt_and_name_for_the_regions(
    genome_regions, gene_names_with_essentiality
):
    genome_regions["Name"] = genome_regions["Systematic ID"].apply(
        add_info_based_on_systematic_id,
        info_column="Name",
        gene_names_with_essentiality=gene_names_with_essentiality,
    )
    genome_regions["Essentiality"] = genome_regions["Systematic ID"].apply(
        add_info_based_on_systematic_id,
        info_column="Essentiality",
        gene_names_with_essentiality=gene_names_with_essentiality,
    )
    genome_regions["Length"] = genome_regions["End"] - genome_regions["Start"]

    return genome_regions


if __name__ == "__main__":
    main()
