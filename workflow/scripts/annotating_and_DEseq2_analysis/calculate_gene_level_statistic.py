# %%
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
from typing import List, Tuple
from scipy.stats import permutation_test
import scipy.stats as stats
import multiprocessing as mp
from functools import partial
from sklearn.utils import resample
import warnings
from scipy.stats import ConstantInputWarning
warnings.filterwarnings("ignore", category=ConstantInputWarning)

def load_data(lfc_path: str, annotations_path: str) -> Tuple[pd.DataFrame, pd.DataFrame, pd.Index]:
    """Load and preprocess the data."""
    results = pd.read_csv(lfc_path, index_col=[0, 1, 2, 3], header=[0, 1])
    insertion_annotations = pd.read_csv(annotations_path, index_col=[0, 1, 2, 3])
    in_gene_insertions = insertion_annotations.query(
        "Type != 'Intergenic region' & Distance_to_stop_codon > 4").index
    return results, insertion_annotations, in_gene_insertions

def process_lfc_data(results: pd.DataFrame, in_gene_insertions: pd.Index) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Process LFC data."""
    LFCs = results.xs("log2FoldChange", level=1, axis=1)
    pvalues = results.xs("padj", level=1, axis=1)
    in_gene_LFCs = LFCs[LFCs.index.isin(in_gene_insertions)].copy()
    in_gene_pvalues = pvalues[pvalues.index.isin(in_gene_insertions)].copy()
    return in_gene_LFCs, in_gene_pvalues

def calculate_weights(GMs: pd.DataFrame) -> pd.DataFrame:
    """Calculate weights for LFC values."""
    GMs["Weights"] = -np.log10(GMs["Padj"].apply(lambda x: 1e-10 if x <= 1e-10 else 1 if x > 1-1e-10 else x))
    return GMs

def merge_and_annotate(GMs: pd.DataFrame, insertion_annotations: pd.DataFrame) -> pd.DataFrame:
    """Merge GMs with annotations and normalize weights."""
    GMs_annotated = pd.merge(GMs, insertion_annotations, how="left", left_index=True, right_index=True)
    GMs_annotated = GMs_annotated[GMs_annotated["Weights"].notna()].copy()
    GMs_annotated["Normalized_weights"] = GMs_annotated.groupby(["Systematic ID", "Timepoint"])["Weights"].transform(
        lambda x: x/x.sum())
    return GMs_annotated

def test_statistic(x, weights):
    return np.sum(x * weights)

def fisher_method(pvalues):
    chi_square = -2 * np.sum(np.log(pvalues))
    combined_pvalue = 1 - stats.chi2.cdf(chi_square, df=2*len(pvalues))
    return combined_pvalue

def stouffer_method(pvalues):
    z_scores = stats.norm.ppf(1 - np.array(pvalues))
    z_combined = np.sum(z_scores) / np.sqrt(len(pvalues))
    return 1 - stats.norm.cdf(z_combined)

def brown_method(pvalues, correlations):
    k = len(pvalues)
    fisher_stat = -2 * np.sum(np.log(pvalues))
    
    # Calculate S (covariance term)
    S = np.sum(correlations) - k
    
    # Calculate c (scaling factor)
    c = (k + S) / (2 * k)
    
    # Adjusted Fisher statistic
    fisher_stat_adj = fisher_stat / c
    
    # Adjusted degrees of freedom
    df_adj = (2 * k) / c
    
    # Calculate combined p-value
    p_combined = 1 - stats.chi2.cdf(fisher_stat_adj, df_adj)
    
    return p_combined

def bootstrap_correlation(x, y, n_bootstrap=1000):
    correlations = []
    for _ in range(n_bootstrap):
        indices = np.random.choice(len(x), len(x), replace=True)
        corr, _ = stats.pearsonr(x[indices], y[indices])
        correlations.append(corr)
    # remove nan
    correlations = [corr for corr in correlations if not np.isnan(corr)]
    return np.mean(correlations)

def calculate_correlation(i, j, gene_fitness, n_bootstrap):
        corr = bootstrap_correlation(gene_fitness.iloc[i].values, gene_fitness.iloc[j].values, n_bootstrap)
        return i, j, corr

def estimate_correlation_matrix(gene_fitness, n_bootstrap=1000, cores=16):

    genes = gene_fitness.index
    correlation_matrix = pd.DataFrame(np.eye(len(genes)), index=genes, columns=genes)
    
    with mp.Pool(cores) as pool:
        partial_calc = partial(calculate_correlation, gene_fitness=gene_fitness, n_bootstrap=n_bootstrap)
        results = pool.starmap(partial_calc, [(i, j) for i in range(gene_fitness.shape[0]) for j in range(i+1, gene_fitness.shape[0])])

    for i, j, corr in results:
        correlation_matrix.iloc[i, j] = corr
        correlation_matrix.iloc[j, i] = corr

    return correlation_matrix

def gene_level_lfc_and_pvalue(gene_df: pd.DataFrame, n_permutations: int = 1000) -> Tuple[List[float], List[float]]:
    """Perform permutation test for LFC values."""
    # drop the index duplicated rows
    gene_df = gene_df.reset_index().groupby(
        ["#Chr", "Coordinate", "Strand", "Target", "Timepoint", "Systematic ID", "Name", "Essentiality", "Transcript"]
        ).filter(lambda x: x["Transcript"].unique()[0].endswith(".1")).set_index(
            ["#Chr", "Coordinate", "Strand", "Target", "Timepoint"])
    unstacked_gene_df = gene_df.unstack(level="Timepoint")
    # correlation_matrix = estimate_correlation_matrix(unstacked_gene_df["LFC"]).values
    gene_level_lfc_pvalues = pd.DataFrame()
    for tp, tp_df in gene_df.groupby("Timepoint"):
        lfcs = tp_df["LFC"].values
        padjs = tp_df["Padj"].values
        normalized_weights = tp_df["Normalized_weights"].values
        observed_gene_level_lfc = np.average(lfcs, weights=normalized_weights)
        
        try:
            # res = permutation_test(
            #     (lfcs, normalized_weights),
            #     test_statistic,
            #     n_resamples=n_permutations,
            #     alternative='two-sided',
            #     permutation_type='pairings',
            # )
            gene_level_lfc_pvalue = fisher_method(padjs)
            # gene_level_lfc_pvalue = stouffer_method(padjs)
            # gene_level_lfc_pvalue = brown_method(padjs, correlation_matrix)
        except ValueError as e:
            gene_level_lfc_pvalue = padjs[0]
        gene_level_lfc_pvalues.loc[tp, "LFC"] = observed_gene_level_lfc
        gene_level_lfc_pvalues.loc[tp, "pvalue"] = gene_level_lfc_pvalue
        
    return gene_level_lfc_pvalues.sort_index()["LFC"].tolist(), gene_level_lfc_pvalues.sort_index()["pvalue"].tolist()

def calculate_gene_level_statistics(GMs_annotated: pd.DataFrame) -> pd.DataFrame:
    """Calculate gene-level statistics including weighted average LFC and p-value."""
    GWMs = pd.DataFrame()
    n = 0
    total_genes = len(GMs_annotated["Systematic ID"].unique())
    for (sysID, gene, ess), gene_df in GMs_annotated.groupby(["Systematic ID", "Name", "Essentiality"]):
        lfc, pvalue = gene_level_lfc_and_pvalue(gene_df)
        
        GWMs.loc[sysID, "Name"] = gene
        GWMs.loc[sysID, "Essentiality"] = ess
        GWMs.loc[sysID, ("YES0", "YES1", "YES2", "YES3", "YES4")] = lfc
        GWMs.loc[sysID, ("YES0_pvalue", "YES1_pvalue", "YES2_pvalue", "YES3_pvalue", "YES4_pvalue")] = pvalue
        n += 1
        print(f"################# {n}/{total_genes} ################", end="\n", flush=True)
    
    return GWMs[GWMs.notna().all(axis=1)].copy()

def parse_args():
    parser = argparse.ArgumentParser(description="Calculate gene-level statistics")
    parser.add_argument("-l", "--lfc_path", type=Path, required=True, help="Path to the LFC file")
    parser.add_argument("-a", "--annotations_path", type=Path, required=True, help="Path to the annotations file")
    parser.add_argument("-o", "--output_path", type=Path, required=True, help="Path to the output file")
    return parser.parse_args()

# Main execution
def main():
    args = parse_args()

    results, insertion_annotations, in_gene_insertions = load_data(args.lfc_path, args.annotations_path)
    in_gene_LFCs, in_gene_pvalues = process_lfc_data(results, in_gene_insertions)

    GMs = pd.merge(in_gene_LFCs.stack().to_frame("LFC"), in_gene_pvalues.stack().to_frame("Padj"), 
                how="left", left_index=True, right_index=True).rename_axis(["#Chr", "Coordinate", "Strand", "Target", "Timepoint"])
    GMs = calculate_weights(GMs)

    GMs_annotated = merge_and_annotate(GMs, insertion_annotations)
    GWMs = calculate_gene_level_statistics(GMs_annotated)

    GWMs.to_csv(args.output_path, index=True)

if __name__ == "__main__":
    main()
# %%
