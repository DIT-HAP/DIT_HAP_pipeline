# %%
from pathlib import Path
import argparse
import numpy as np
import pandas as pd

def filter_insertions(insertion_annotations):
    intergenic_insertions_filtered = insertion_annotations[
        (insertion_annotations["Type"] == "Intergenic region") & 
        (insertion_annotations["Distance_to_region_start"] > 500) & 
        (insertion_annotations["Distance_to_region_end"] > 500)
    ].index

    in_gene_insertions = insertion_annotations.query(
        "Type != 'Intergenic region' & Distance_to_stop_codon > 3"
    ).index

    return intergenic_insertions_filtered, in_gene_insertions

def transfer_FR_index(idxs):
    idxs = list(idxs)
    idxs[2] = "+" if idxs[2] == "-" else "-"
    return tuple(idxs)

def impute_missing_values(in_gene_counts_df):
    stacked_df = in_gene_counts_df.stack(level=0, dropna=False)
    stacked_dropna_df = in_gene_counts_df.stack(level=0, dropna=True)
    
    in_gene_isna_idx = stacked_df[stacked_df.isna().all(axis=1)].index
    in_gene_complementary_idx = [transfer_FR_index(idx) for idx in in_gene_isna_idx]
    in_gene_index_for_imputation = list(set(in_gene_complementary_idx) & set(stacked_dropna_df.index))
    in_gene_has_complementary_idxs = [transfer_FR_index(idx) for idx in in_gene_index_for_imputation]

    stacked_df.loc[in_gene_has_complementary_idxs, :] = stacked_df.loc[in_gene_index_for_imputation, :].values

    return stacked_df.unstack().reorder_levels([1, 0], axis=1), in_gene_has_complementary_idxs

def parse_args():
    parser = argparse.ArgumentParser(description="Impute missing values using FR.")
    parser.add_argument("-i", "--input", type=Path, required=True, help="Path to the input file.")
    parser.add_argument("-a", "--annotation", type=Path, required=True, help="Path to the annotation file.")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Path to the output file.")
    return parser.parse_args()

def main():

    args = parse_args()

    counts_df = pd.read_csv(args.input, index_col=[0, 1, 2, 3], sep="\t", header=[0,1])
    insertion_annotations = pd.read_csv(args.annotation, index_col=[0, 1, 2, 3], sep="\t")
    
    _, in_gene_insertions = filter_insertions(insertion_annotations)

    in_gene_counts_df = counts_df[counts_df.index.isin(in_gene_insertions)].copy()
    intergenic_counts_df = counts_df[~counts_df.index.isin(in_gene_insertions)].copy().dropna(axis=0, how="any")
    print("Insertions with all replicates available in intergenic regions:", intergenic_counts_df.shape[0])
    print("Insertions with at least one replicate available in coding genes:", in_gene_counts_df.shape[0])
    print("Insertions with all replicates available in coding genes:", in_gene_counts_df.dropna(axis=0, how="any").shape[0])

    imputed_in_gene_counts_df, in_gene_has_complementary_idxs = impute_missing_values(in_gene_counts_df)

    imputed_in_gene_counts_df_noNA = imputed_in_gene_counts_df.dropna(axis=0, how="any")
    print("Insertions with all replicates available in coding genes after imputation:", imputed_in_gene_counts_df_noNA.shape[0])

    imputed_counts = pd.concat([intergenic_counts_df, imputed_in_gene_counts_df_noNA], axis=0)

    imputed_counts.to_csv(args.output, index=True, sep="\t")

    print("### Impute missing values using FR completed ###")
    insertion_num = counts_df.shape[0]
    ingene_num = in_gene_counts_df.shape[0]
    intergenic_num = counts_df[~counts_df.index.isin(in_gene_insertions)].shape[0]

    noNA_insertion_num = counts_df.dropna(axis=0, how="any").shape[0]
    noNA_ingene_num = in_gene_counts_df.dropna(axis=0, how="any").shape[0]
    noNA_intergenic_num = intergenic_counts_df.shape[0]

    noNA_imputed_ingene_num = imputed_in_gene_counts_df_noNA.shape[0]
    increased_ingene_num = noNA_imputed_ingene_num - noNA_ingene_num

    print("*** Total insertions:", insertion_num)
    print("*** Insertions in coding genes", ingene_num, "({:.2f}%)".format(ingene_num / insertion_num * 100))
    print("*** Insertions in intergenic regions:", intergenic_num, "({:.2f}%)".format(intergenic_num / insertion_num * 100))

    print("*** Insertions with all replicates available:", noNA_insertion_num, "({:.2f}%)".format(noNA_insertion_num / insertion_num * 100))
    print("*** Insertions with all replicates available in coding genes:", noNA_ingene_num, "({:.2f}%)".format(noNA_ingene_num / noNA_insertion_num * 100), "Compared with insertions with all replicates available")
    print("*** Insertions with all replicates available in intergenic regions:", noNA_intergenic_num, "({:.2f}%)".format(noNA_intergenic_num / noNA_insertion_num * 100), "Compared with insertions with all replicates available")

    print("*** Insertions with all replicates available in coding genes after imputation:", noNA_imputed_ingene_num, "({:.2f}%)".format(noNA_imputed_ingene_num / noNA_insertion_num * 100))
    print("*** Inscrease in insertions with all replicates available in coding genes:", increased_ingene_num, "({:.2f}%)".format(increased_ingene_num / noNA_insertion_num * 100))

    print("*** Insertions with all replicates available in coding genes after imputation in all in-gene insertions:", noNA_imputed_ingene_num, "({:.2f}%)".format(noNA_imputed_ingene_num / ingene_num * 100))

if __name__ == "__main__":
    main()