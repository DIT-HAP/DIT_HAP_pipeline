"""

"""
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
from multiprocessing import Pool
from util.curve_fitting_functions import curve_fitting


def parse_args():
    # Parse arguments
    parser = argparse.ArgumentParser(description="Calculate SDR.")
    parser.add_argument(
        "-C",
        "--cores",
        dest="cores",
        type=int,
        default=8,
        help="Number of cores to use",
    )
    parser.add_argument(
        "-anno",
        "--annotated-insertion",
        dest="annotated_insertion_file",
        required=True,
        type=Path,
        help="File of annotated insertion",
    )
    parser.add_argument(
        "-GM",
        "--GM",
        dest="GM",
        required=True,
        type=Path,
        help="File of M values",
    )
    parser.add_argument(
        "--use-weighted-M",
        dest="use_weighted_M",
        type=int,
        default=1,
        help="Use weighted M or not",
    )
    parser.add_argument(
        "--use-LOWESS",
        dest="use_LOWESS",
        type=int,
        default=1,
        help="Use LOWESS or not",
    )
    parser.add_argument(
        "-ddr",
        "--DDR",
        dest="DDR",
        required=True,
        type=Path,
        help="File of domain-level DR",
    )
    parser.add_argument(
        "-fatol",
        "--fatol",
        dest="fatol",
        required=True,
        type=float,
        help="fatol",
    )
    parser.add_argument(
        "-f",
        "--fitting",
        dest="fitting",
        required=True,
        type=str,
        default="points5",
        help="fitting",
    )
    parser.add_argument(
        "-p",
        "--prediction",
        dest="prediction",
        required=True,
        type=str,
        default="points5",
        help="prediction",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        required=True,
        type=Path,
        help="Output file for M",
    )

    return parser.parse_args()


def chunk_gene_level_insertions(cores, in_gene_with_M_index, annotated_insertions):

    annotations = annotated_insertions.loc[in_gene_with_M_index].copy()
    annotations = annotations.reset_index(drop=False).set_index(
        "Systematic ID")
    uniqued_genes = annotations.index.unique()

    chunked_genes = np.array_split(uniqued_genes, cores)

    chunked_gene_insertions = [
        annotations.loc[genes].copy().reset_index(
            drop=False).set_index(["#Chr", "Coordinate", "Strand", "Target"])
        for genes in chunked_genes
    ]
    return chunked_gene_insertions


def main():

    args = parse_args()

    GMs = pd.read_csv(args.GM, index_col=[0, 1, 2, 3], header=0)

    annotated_insertions = pd.read_csv(
        args.annotated_insertion_file, index_col=[0, 1, 2, 3], header=[0]
    )

    DDR = pd.read_csv(args.DDR, header=[0])
    DDR = DDR.groupby("Systematic ID").apply(
        calcualte_DDR_ratio).reset_index(drop=True)

    domains_ratio_gt_point6 = pd.Index(
        DDR[(DDR["DR_ratio"] >= 0.6) & (DDR["Confidence score"] >= 50)][
            ["Systematic ID", "domain_id", "domain_residues"]
        ]
    )

    in_ratio_gt_point6_domains_index = annotated_insertions[
        pd.Index(
            annotated_insertions[["Systematic ID",
                                  "domain_id", "domain_residues"]]
        ).isin(domains_ratio_gt_point6) & (annotated_insertions["Distance_to_stop_codon"] >= 4)
    ].index

    in_gene_with_M_index = GMs.index.intersection(
        in_ratio_gt_point6_domains_index)

    chunked_gene_insertions = chunk_gene_level_insertions(
        args.cores, in_gene_with_M_index, annotated_insertions)

    generations = GMs.reset_index(drop=True)[["Sample", "Timepoint", "G"]].drop_duplicates().groupby("Timepoint").apply(
        lambda x: np.average(x["G"])).sort_values()

    use_LOWESS = bool(int(args.use_LOWESS))
    use_weighted_M = bool(args.use_weighted_M)
    gene_level_fitting_params = [
        (generations, gene_df, GMs, args.fatol, use_weighted_M, use_LOWESS, args.fitting, args.prediction) for gene_df in chunked_gene_insertions]
    with Pool(args.cores) as pool:
        results = pool.starmap(
            gene_level_fitting, gene_level_fitting_params)
    gene_level_weighted_M = pd.concat(results, axis=0)
    gene_level_weighted_M.rename_axis("Systematic ID", axis=0, inplace=True)
    gene_level_weighted_M.to_csv(
        args.output, header=True, index=True, float_format="%.3f"
    )


def calcualte_DDR_ratio(DDR_sub_df):
    max_DDR = DDR_sub_df["DR"].max()
    DDR_sub_df["DR_ratio"] = round(DDR_sub_df["DR"] / max_DDR, 3)

    return DDR_sub_df


def gene_level_fitting(
    generation,
    sub_gene_insertions,
    GMs,
    fatol,
    using_weighted_M,
    LOWESS_smoothing,
    fitting,
    prediction
):
    gene_level_weighted_M = pd.DataFrame()

    for gene, gene_df in sub_gene_insertions.groupby("Systematic ID"):

        insertions_in_current_level = GMs.loc[GMs.index.isin(
            gene_df.index)].copy()

        GWM = curve_fitting(
            gene, generation, insertions_in_current_level, fatol, useWeightedM=using_weighted_M, useLOWESS=LOWESS_smoothing, fitting=fitting, prediction=prediction)

        GWM_series = pd.Series(GWM, name=gene)

        gene_level_weighted_M = pd.concat(
            [gene_level_weighted_M, GWM_series.to_frame().T]
        )

    return gene_level_weighted_M


if __name__ == "__main__":

    main()
