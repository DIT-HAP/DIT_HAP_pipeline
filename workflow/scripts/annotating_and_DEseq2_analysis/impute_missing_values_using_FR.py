# %%
from pathlib import Path
import argparse
import numpy as np
import pandas as pd

def load_raw_reads(input_files):
    return {
        file.stem: pd.read_csv(file, index_col=[0, 1, 2, 3])
        for file in input_files
    }

def filter_counts_by_replicates(counts_df):
    count_df_3rep = counts_df.loc[counts_df.notna().sum(axis=1) == 18].copy()
    count_df_2rep = counts_df.loc[counts_df.notna().sum(axis=1) == 12].copy()
    count_df_1rep = counts_df.loc[counts_df.notna().sum(axis=1) == 6].copy()

    count_df_2rep = count_df_2rep[(count_df_2rep.filter(like="0h") > 60).sum(axis=1) > 1]
    count_df_1rep = count_df_1rep[(count_df_1rep.filter(like="0h") > 100).sum(axis=1) > 0]

    return pd.concat([count_df_3rep, count_df_2rep, count_df_1rep], axis=0)

def load_insertion_annotations(file_path):
    return pd.read_csv(file_path, index_col=[0, 1, 2, 3])

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

    return stacked_df.unstack().reorder_levels([1, 0], axis=1)

def parse_args():
    parser = argparse.ArgumentParser(description="Impute missing values using FR.")
    parser.add_argument("-i", "--input", type=Path, nargs='+', required=True, help="Path to the input file.")
    parser.add_argument("-a", "--annotation", type=Path, required=True, help="Path to the annotation file.")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Path to the output file.")
    return parser.parse_args()

def main():

    args = parse_args()

    raw_reads = load_raw_reads(args.input)
    counts_df = pd.concat(raw_reads, axis=1)
    counts_df = filter_counts_by_replicates(counts_df)

    insertion_annotations = load_insertion_annotations(args.annotation)

    _, in_gene_insertions = filter_insertions(insertion_annotations)

    in_gene_counts_df = counts_df[counts_df.index.isin(in_gene_insertions)].copy()
    imputed_in_gene_counts_df = impute_missing_values(in_gene_counts_df)

    imputed_in_gene_count_df_3rep = imputed_in_gene_counts_df.loc[imputed_in_gene_counts_df.notna().sum(axis=1) == 18].copy()
    imputed_in_gene_count_df_2rep = imputed_in_gene_counts_df.loc[imputed_in_gene_counts_df.notna().sum(axis=1) == 12].copy()
    imputed_in_gene_count_df_1rep = imputed_in_gene_counts_df.loc[imputed_in_gene_counts_df.notna().sum(axis=1) == 6].copy()

    new_counts_3p_df = pd.concat([counts_df.loc[counts_df.notna().sum(axis=1) == 18], imputed_in_gene_count_df_3rep],
                                 axis=0).drop_duplicates(keep="first")

    new_counts_3p_df.to_csv(args.output, index=True)

if __name__ == "__main__":
    main()