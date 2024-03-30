import numpy as np
import pandas as pd
import re


def parse_gff_file(gff_file):
    GFF = pd.read_csv(
        gff_file,
        sep="\t",
        comment="#",
        names=[
            "Chr",
            "Source",
            "Feature",
            "Start",
            "End",
            "Score",
            "Strand",
            "Frame",
            "Attribute",
        ],
    )

    extract_systematic_ID_pattern = re.compile(
        r"ID=(\S+?)(?:$|(?:|\.\d(?::\S+|));)")
    GFF["Systematic ID"] = GFF["Attribute"].str.extract(
        extract_systematic_ID_pattern, expand=False
    )
    GFF["Transcript"] = GFF["Attribute"].str.extract(
        "Parent=(\S+?)$", expand=False)

    return GFF


def extract_coding_gene(gff):
    coding_genes = gff[(gff["Feature"] == "mRNA")
                       ]["Systematic ID"].unique().tolist()

    # switch 1-based to 0-based
    region_of_coding_genes = (
        gff[gff["Systematic ID"].isin(coding_genes) & (gff["Feature"] == "CDS")][
            ["Chr", "Start", "End", "Strand", "Systematic ID"]
        ]
        .groupby(["Chr", "Strand", "Systematic ID"])
        .aggregate({"Start": (lambda x: np.min(x) - 1), "End": "max"})
        .reset_index()[["Chr", "Start", "End", "Systematic ID", "Strand"]]
        .rename(columns={"Chr": "#Chr"})
        .sort_values(["#Chr", "Start", "End"])
    )
    region_of_coding_genes.insert(4, "Type", "Coding gene")
    region_of_coding_genes["Length"] = (
        region_of_coding_genes["End"] - region_of_coding_genes["Start"]
    )

    return region_of_coding_genes


def obtain_representative_coding_transcript(gff, peptide_stats):
    peptide_stats = pd.read_csv(peptide_stats, sep="\t")

    coding_genes = gff[(gff["Feature"] == "mRNA")
                       ]["Systematic ID"].unique().tolist()

    coding_transcripts = gff[(gff["Feature"] == "CDS") & gff["Systematic ID"].isin(coding_genes)
                             ]["Transcript"].unique().tolist()

    # switch 1-based to 0-based
    protein_length_of_coding_transcripts = (
        gff[gff["Transcript"].isin(coding_transcripts) & (gff["Feature"] == "CDS")][
            ["Chr", "Start", "End", "Strand", "Systematic ID", "Transcript"]
        ]
        .groupby(["Chr", "Strand", "Systematic ID", "Transcript"])
        # remove the stop codon
        .apply(lambda sub_df: (((sub_df["End"]-sub_df["Start"]).abs()+1).sum()-3)/3)
        .rename("Protein_length")
    )

    transcript_and_protein_residues = pd.merge(protein_length_of_coding_transcripts.reset_index(), peptide_stats,
                                               left_on="Systematic ID", right_on="Systematic_ID", how="outer")
    transcript_and_protein_residues["Delta"] = transcript_and_protein_residues["Residues"] - \
        transcript_and_protein_residues["Protein_length"]
    representative_transcript = transcript_and_protein_residues[transcript_and_protein_residues["Delta"] == 0].copy(
    )
    return representative_transcript


def extract_representative_coding_transcript(gff, peptide_stats):

    representative_transcript = obtain_representative_coding_transcript(
        gff, peptide_stats)

    # switch 1-based to 0-based
    region_of_representative_coding_transcript = (
        gff[gff["Transcript"].isin(representative_transcript["Transcript"]) & (gff["Feature"] == "CDS")][
            ["Chr", "Start", "End", "Strand", "Systematic ID"]
        ]
        .groupby(["Chr", "Strand", "Systematic ID"])
        .aggregate({"Start": (lambda x: np.min(x) - 1), "End": "max"})
        .reset_index()[["Chr", "Start", "End", "Systematic ID", "Strand"]]
        .rename(columns={"Chr": "#Chr"})
        .sort_values(["#Chr", "Start", "End"])
    )
    region_of_representative_coding_transcript.insert(4, "Type", "Coding gene")
    region_of_representative_coding_transcript["Length"] = (
        region_of_representative_coding_transcript["End"] -
        region_of_representative_coding_transcript["Start"]
    )

    return region_of_representative_coding_transcript


def extract_CDS_intron(gff, peptide_stats):

    representative_transcript = obtain_representative_coding_transcript(
        gff, peptide_stats)

    # CDS and intron of coding genes
    CDS_intron_of_representative_transcripts = gff[
        (gff["Transcript"].isin(representative_transcript["Transcript"]))
        & (gff["Feature"].isin(["CDS", "intron"]))
    ][
        ["Chr", "Start", "End", "Strand", "Feature", "Systematic ID", "Transcript"]
    ].rename(
        columns={"Feature": "Type"}
    ).reset_index()

    # switch 1-based to 0-based
    CDS_intron_of_representative_transcripts["Start"] = CDS_intron_of_representative_transcripts["Start"] - 1

    # calculate length
    CDS_intron_of_representative_transcripts["Length"] = (
        CDS_intron_of_representative_transcripts["End"] -
        CDS_intron_of_representative_transcripts["Start"]
    )

    CDS_intron_of_representative_transcripts = CDS_intron_of_representative_transcripts[[
        "Chr", "Start", "End", "Strand", "Type", "Systematic ID", "Transcript", "Length"]]

    # # use the longest transcript as representative transcript
    # representative_transcripts = (
    #     coding_gene_transcript.to_frame("Transcript_length")
    #     .reset_index(drop=False)
    #     .groupby("Systematic ID")
    #     .apply(
    #         lambda sub_df: sub_df.loc[
    #             sub_df["Transcript_length"].idxmax(), "Transcript"
    #         ]
    #     )
    # )

    # extarct the gene region
    region_of_coding_representative_transcripts = (
        gff[gff["Transcript"].isin(representative_transcript["Transcript"]) & (gff["Feature"] == "CDS")][
            ["Chr", "Start", "End", "Strand", "Systematic ID"]
        ]
        .groupby(["Chr", "Strand", "Systematic ID"])
        .aggregate({"Start": (lambda x: np.min(x) - 1), "End": "max"})
        .reset_index()[["Chr", "Start", "End", "Systematic ID", "Strand"]]
        .sort_values(["Chr", "Start", "End"])
    )

    region_of_coding_representative_transcripts["Length"] = (
        region_of_coding_representative_transcripts["End"] -
        region_of_coding_representative_transcripts["Start"]
    )

    # merge gene region with CDS and intron
    CDS_intron_of_coding_representative_transcripts_with_gene_info = pd.merge(
        CDS_intron_of_representative_transcripts,
        region_of_coding_representative_transcripts[[
            "Systematic ID", "Start", "End", "Length"]],
        on=["Systematic ID"],
        how="left",
        suffixes=["", "_Region"],
    ).rename(
        columns={"Chr": "#Chr"}
    )

    # filter out intron outside of gene region
    CDS_intron_of_coding_representative_transcripts_with_gene_info_filtered = CDS_intron_of_coding_representative_transcripts_with_gene_info[
        (CDS_intron_of_coding_representative_transcripts_with_gene_info["Start"] >= CDS_intron_of_coding_representative_transcripts_with_gene_info["Start_Region"]) & (
            CDS_intron_of_coding_representative_transcripts_with_gene_info[
                "End"] <= CDS_intron_of_coding_representative_transcripts_with_gene_info["End_Region"]
        )
    ][["#Chr", "Start", "End", "Systematic ID", "Type", "Strand", "Length", "Transcript", "Start_Region", "End_Region", "Length_Region"]].sort_values(
        ["#Chr", "Start", "End"]
    ).reset_index(drop=True)

    return CDS_intron_of_coding_representative_transcripts_with_gene_info_filtered


def cal_accumlated_CDS_bases(sub_df):

    sub_df_sorted = sub_df.sort_values(["Start"]).copy()
    accumulated_CDS_bases = 0
    strand = sub_df_sorted["Strand"].iloc[0]
    if strand == "+":
        index_order = sub_df_sorted.index
    else:
        index_order = sub_df_sorted.index[::-1]
    for idx in index_order:
        if sub_df_sorted.loc[idx, "Type"] == "CDS":
            sub_df_sorted.loc[idx,
                              "Accumulated_CDS_bases"] = accumulated_CDS_bases
            accumulated_CDS_bases += sub_df_sorted.loc[idx, "Length"]
        else:
            sub_df_sorted.loc[idx,
                              "Accumulated_CDS_bases"] = accumulated_CDS_bases

    return (sub_df_sorted)
