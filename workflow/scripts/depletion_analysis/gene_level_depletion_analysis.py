#!/usr/bin/env python3
"""
Gene-Level Depletion Analysis for Transposon Insertion Sequencing

This script performs gene-level aggregation of transposon insertion depletion data,
combining multiple insertion sites within genes using weighted averages and statistical
meta-analysis methods. It processes log2 fold changes (LFC) and adjusted p-values
from insertion-level analysis to generate gene-level fitness scores.

The analysis pipeline includes:
1. Loading insertion-level LFC and p-value data
2. Filtering for in-gene insertions based on annotations
3. Calculating statistical weights from p-values
4. Aggregating insertions per gene using weighted averages
5. Combining p-values using Fisher's method or alternatives
6. Generating gene-level fitness profiles across time points

Statistical Methods:
- Fisher's method: Combines p-values assuming independence
- Stouffer's method: Z-score based combination with equal weights
- Brown's method: Fisher's method adjusted for correlation structure

Typical usage:
    python gene_level_depletion_analysis.py -l lfc_results.csv -a annotations.csv -o gene_results.csv

Input: CSV files with insertion-level LFC/p-values and genomic annotations
Output: CSV file with gene-level fitness scores and combined p-values

Author: Bioinformatics Pipeline
Date: 2024
"""

import logging
import numpy as np
import pandas as pd
import multiprocessing as mp
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Union
import argparse
from functools import partial
from collections import defaultdict
import time
import warnings

import scipy.stats as stats
from scipy.stats import permutation_test, ConstantInputWarning
from sklearn.utils import resample

# Suppress specific scipy warnings for cleaner output
warnings.filterwarnings("ignore", category=ConstantInputWarning)

# Statistical significance thresholds
SIGNIFICANCE_THRESHOLDS = {
    'high': 0.001,
    'medium': 0.01, 
    'low': 0.05
}

# Default parameters
DEFAULT_N_PERMUTATIONS = 1000
DEFAULT_N_BOOTSTRAP = 1000
DEFAULT_CORES = mp.cpu_count() - 1


def setup_logging(verbose: bool = False) -> None:
    """Set up logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def load_data(lfc_path: Path, annotations_path: Path) -> Tuple[pd.DataFrame, pd.DataFrame, pd.Index]:
    """
    Load and validate insertion-level LFC data and genomic annotations.
    
    Args:
        lfc_path: Path to CSV file with LFC and p-value results
        annotations_path: Path to CSV file with insertion annotations
        
    Returns:
        Tuple of (results_df, annotations_df, in_gene_insertions_index)
        
    Raises:
        FileNotFoundError: If input files don't exist
        ValueError: If required columns are missing
    """
    logging.info(f"Loading LFC data from {lfc_path}")
    logging.info(f"Loading annotations from {annotations_path}")
    
    if not lfc_path.exists():
        raise FileNotFoundError(f"LFC file not found: {lfc_path}")
    if not annotations_path.exists():
        raise FileNotFoundError(f"Annotations file not found: {annotations_path}")
    
    try:
        # Load results with multi-level columns and index
        results = pd.read_csv(lfc_path, index_col=[0, 1, 2, 3], header=[0, 1], sep="\t")
        insertion_annotations = pd.read_csv(annotations_path, index_col=[0, 1, 2, 3], sep="\t")
        
        # Validate required columns exist
        required_level1_cols = ['log2FoldChange', 'padj']
        available_level1_cols = results.columns.get_level_values(1).unique()
        
        missing_cols = [col for col in required_level1_cols if col not in available_level1_cols]
        if missing_cols:
            raise ValueError(f"Missing required columns in LFC data: {missing_cols}")
        
        # Filter for in-gene insertions
        in_gene_insertions = insertion_annotations.query(
            "(Type != 'Intergenic region') and (Distance_to_stop_codon > 4)"
        ).index
        in_gene_insertion_count = results[results.index.isin(in_gene_insertions)].shape[0]
        
        logging.info(f"Loaded {results.shape[0]} total insertions")
        logging.info(f"Found {in_gene_insertion_count} in-gene insertions")
        logging.info(f"Time points in data: {results.columns.get_level_values(0).unique().tolist()}")
        
        return results, insertion_annotations, in_gene_insertions
        
    except pd.errors.EmptyDataError:
        raise ValueError("Input files are empty")
    except Exception as e:
        raise ValueError(f"Error loading data: {e}")


def process_lfc_data(results: pd.DataFrame, in_gene_insertions: pd.Index) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Extract and filter LFC and p-value data for in-gene insertions.
    
    Args:
        results: Multi-level DataFrame with LFC and p-value data
        in_gene_insertions: Index of insertions located within genes
        
    Returns:
        Tuple of (in_gene_LFCs, in_gene_pvalues) DataFrames
    """
    logging.info("Processing LFC and p-value data")
    
    # Extract LFC and p-value data
    LFCs = results.xs("log2FoldChange", level=1, axis=1)
    pvalues = results.xs("padj", level=1, axis=1)
    
    # Filter for in-gene insertions
    in_gene_LFCs = LFCs[LFCs.index.isin(in_gene_insertions)].copy()
    in_gene_pvalues = pvalues[pvalues.index.isin(in_gene_insertions)].copy()
    
    # Data quality checks
    lfc_nan_count = in_gene_LFCs.isna().any(axis=1).sum()
    pval_nan_count = in_gene_pvalues.isna().any(axis=1).sum()
    
    if lfc_nan_count > 0:
        logging.warning(f"Found {lfc_nan_count} rows with NaN values in LFC data")
    if pval_nan_count > 0:
        logging.warning(f"Found {pval_nan_count} rows with NaN values in p-value data")
        logging.warning(f"Fill NaN values in p-value data with 1 since they are mainly due to the large cooks distance")
        logging.warning(f"Note on p-values set to NA: https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#more-information-on-results-columns")
        in_gene_pvalues = in_gene_pvalues.fillna(1)
    
    logging.info(f"Processed {in_gene_LFCs.shape[0]} in-gene LFC measurements")
    
    return in_gene_LFCs, in_gene_pvalues


def calculate_statistical_weights(GMs: pd.DataFrame, min_pvalue: float = 1e-10) -> pd.DataFrame:
    """
    Calculate statistical weights from adjusted p-values using -log10 transformation.
    
    Args:
        GMs: DataFrame with 'Padj' column
        min_pvalue: Minimum p-value to prevent infinite weights
        
    Returns:
        DataFrame with added 'Weights' column
    """
    logging.info("Calculating statistical weights from p-values")
    
    # Clip p-values to prevent infinite weights
    logging.info(f"Clipping p-values to prevent infinite weights, min_pvalue={min_pvalue}, max_pvalue={1-min_pvalue}")
    clipped_pvalues = GMs["Padj"].clip(lower=min_pvalue, upper=1-min_pvalue)
    
    # Calculate -log10 weights
    GMs["Weights"] = -np.log10(clipped_pvalues)
    
    # Log weight statistics
    weight_stats = GMs["Weights"].describe()
    logging.info(f"Weight statistics: mean={weight_stats['mean']:.2f}, "
                f"median={weight_stats['50%']:.2f}, max={weight_stats['max']:.2f}")
    
    return GMs


def merge_and_normalize_weights(GMs: pd.DataFrame, insertion_annotations: pd.DataFrame) -> pd.DataFrame:
    """
    Merge LFC data with annotations and normalize weights within genes.
    
    Args:
        GMs: DataFrame with LFC, p-values, and weights
        insertion_annotations: DataFrame with genomic annotations
        
    Returns:
        Annotated DataFrame with normalized weights per gene per timepoint
    """
    logging.info("Merging data with annotations and normalizing weights")
    
    # Merge with annotations
    GMs_annotated = pd.merge(
        GMs, insertion_annotations, 
        how="left", left_index=True, right_index=True
    )
    
    # Remove rows with missing weights
    initial_count = GMs_annotated.shape[0]
    GMs_annotated = GMs_annotated[GMs_annotated["Weights"].notna()].copy()
    final_count = GMs_annotated.shape[0]
    
    if initial_count != final_count:
        logging.warning(f"GMs_annotated: Removed {initial_count - final_count} rows with missing weights")
    
    # Normalize weights within each gene and timepoint
    GMs_annotated["Normalized_weights"] = GMs_annotated.groupby(
        ["Systematic ID", "Timepoint"]
    )["Weights"].transform(lambda x: x / x.sum())
    
    # Validate normalization
    weight_sums = GMs_annotated.groupby(["Systematic ID", "Timepoint"])["Normalized_weights"].sum()
    if not np.allclose(weight_sums, 1.0, rtol=1e-10):
        logging.warning("Weight normalization may have numerical issues")
    
    logging.info(f"GMs_annotated: Successfully annotated {GMs_annotated.shape[0]} insertions")
    
    return GMs_annotated


def fisher_method(pvalues: np.ndarray) -> float:
    """
    Combine p-values using Fisher's method.
    
    Fisher's method assumes independence between tests and combines p-values
    using the chi-square distribution: -2 * sum(ln(p_i)) ~ chi2(2k)
    
    Args:
        pvalues: Array of p-values to combine
        
    Returns:
        Combined p-value
    """
    if len(pvalues) == 0:
        return 1.0
    
    # Remove NaN values
    valid_pvalues = pvalues[~np.isnan(pvalues)]
    if len(valid_pvalues) == 0:
        return 1.0
    
    # Clip p-values to prevent log(0)
    clipped_pvalues = np.clip(valid_pvalues, 1e-16, 1.0)
    
    chi_square = -2 * np.sum(np.log(clipped_pvalues))
    combined_pvalue = 1 - stats.chi2.cdf(chi_square, df=2*len(clipped_pvalues))
    
    return combined_pvalue


def stouffer_method(pvalues: np.ndarray) -> float:
    """
    Combine p-values using Stouffer's Z-score method.
    
    Converts p-values to Z-scores and combines them assuming equal weights:
    Z_combined = sum(Z_i) / sqrt(k)
    
    Args:
        pvalues: Array of p-values to combine
        
    Returns:
        Combined p-value
    """
    if len(pvalues) == 0:
        return 1.0
    
    # Remove NaN values
    valid_pvalues = pvalues[~np.isnan(pvalues)]
    if len(valid_pvalues) == 0:
        return 1.0
    
    # Clip p-values to prevent infinite Z-scores
    clipped_pvalues = np.clip(valid_pvalues, 1e-16, 1-1e-16)
    
    # Convert to Z-scores (one-tailed to two-tailed)
    z_scores = stats.norm.ppf(1 - clipped_pvalues)
    z_combined = np.sum(z_scores) / np.sqrt(len(z_scores))
    
    return 1 - stats.norm.cdf(z_combined)


def brown_method(pvalues: np.ndarray, correlations: np.ndarray) -> float:
    """
    Combine p-values using Brown's method (correlation-adjusted Fisher's method).
    
    Adjusts Fisher's method for correlation structure between tests:
    - Calculates scaling factor c based on correlation matrix
    - Adjusts both test statistic and degrees of freedom
    
    Args:
        pvalues: Array of p-values to combine
        correlations: Correlation matrix between tests
        
    Returns:
        Combined p-value
    """
    if len(pvalues) == 0:
        return 1.0
    
    k = len(pvalues)
    if k == 1:
        return pvalues[0]
    
    # Remove NaN values
    valid_pvalues = pvalues[~np.isnan(pvalues)]
    if len(valid_pvalues) == 0:
        return 1.0
    
    # Standard Fisher statistic
    fisher_stat = -2 * np.sum(np.log(np.clip(valid_pvalues, 1e-16, 1.0)))
    
    # Calculate covariance adjustment
    S = np.sum(correlations) - k  # Sum of off-diagonal correlations
    c = (k + S) / (2 * k)  # Scaling factor
    
    # Adjusted statistics
    fisher_stat_adj = fisher_stat / c
    df_adj = (2 * k) / c
    
    # Combined p-value
    p_combined = 1 - stats.chi2.cdf(fisher_stat_adj, df_adj)
    
    return p_combined


def bootstrap_correlation(x: np.ndarray, y: np.ndarray, n_bootstrap: int = 1000) -> float:
    """
    Estimate correlation between two variables using bootstrap resampling.
    
    Args:
        x: First variable
        y: Second variable  
        n_bootstrap: Number of bootstrap samples
        
    Returns:
        Mean correlation across bootstrap samples
    """
    if len(x) != len(y) or len(x) < 3:
        return 0.0
    
    correlations = []
    rng = np.random.default_rng(42)  # For reproducibility
    
    for _ in range(n_bootstrap):
        indices = rng.choice(len(x), len(x), replace=True)
        try:
            corr, _ = stats.pearsonr(x[indices], y[indices])
            if not np.isnan(corr):
                correlations.append(corr)
        except:
            continue
    
    return np.mean(correlations) if correlations else 0.0


def calculate_correlation_worker(args: Tuple[int, int, pd.DataFrame, int]) -> Tuple[int, int, float]:
    """Worker function for parallel correlation calculation."""
    i, j, gene_fitness, n_bootstrap = args
    corr = bootstrap_correlation(
        gene_fitness.iloc[i].values, 
        gene_fitness.iloc[j].values, 
        n_bootstrap
    )
    return i, j, corr


def estimate_correlation_matrix(gene_fitness: pd.DataFrame, 
                              n_bootstrap: int = 1000, 
                              cores: int = DEFAULT_CORES) -> pd.DataFrame:
    """
    Estimate correlation matrix between genes using bootstrap resampling.
    
    Args:
        gene_fitness: DataFrame with genes as rows and timepoints as columns
        n_bootstrap: Number of bootstrap samples per correlation
        cores: Number of CPU cores to use
        
    Returns:
        Correlation matrix between genes
    """
    logging.info(f"Estimating correlation matrix for {len(gene_fitness)} genes")
    
    genes = gene_fitness.index
    correlation_matrix = pd.DataFrame(np.eye(len(genes)), index=genes, columns=genes)
    
    # Prepare arguments for parallel processing
    args_list = [
        (i, j, gene_fitness, n_bootstrap) 
        for i in range(gene_fitness.shape[0]) 
        for j in range(i+1, gene_fitness.shape[0])
    ]
    
    if not args_list:
        return correlation_matrix
    
    # Parallel correlation calculation
    with mp.Pool(cores) as pool:
        results = pool.map(calculate_correlation_worker, args_list)
    
    # Fill correlation matrix
    for i, j, corr in results:
        correlation_matrix.iloc[i, j] = corr
        correlation_matrix.iloc[j, i] = corr
    
    logging.info("Correlation matrix estimation completed")
    
    return correlation_matrix


def gene_level_lfc_and_pvalue(gene_df: pd.DataFrame, 
                             method: str = 'fisher',
                             n_permutations: int = DEFAULT_N_PERMUTATIONS) -> Tuple[List[float], List[float]]:
    """
    Calculate gene-level LFC and combined p-values across timepoints.
    
    Args:
        gene_df: DataFrame with insertion data for a single gene
        method: P-value combination method ('fisher', 'stouffer', 'brown')
        n_permutations: Number of permutations for permutation test
        
    Returns:
        Tuple of (gene_level_lfcs, gene_level_pvalues) across timepoints
    """
    # Remove duplicate transcript entries, keeping only .1 transcripts
    gene_df_clean = gene_df.reset_index().set_index(["Chr", "Coordinate", "Strand", "Target", "Timepoint"])
    
    gene_level_results = pd.DataFrame()
    
    for tp, tp_df in gene_df_clean.groupby("Timepoint"):
        lfcs = tp_df["LFC"].values
        padjs = tp_df["Padj"].values
        normalized_weights = tp_df["Normalized_weights"].values
        
        # Calculate weighted average LFC
        observed_gene_level_lfc = np.average(lfcs, weights=normalized_weights)
        
        # Combine p-values using specified method
        try:
            if method == 'fisher':
                gene_level_pvalue = fisher_method(padjs)
            elif method == 'stouffer':
                gene_level_pvalue = stouffer_method(padjs)
            else:  # Default to Fisher if method not recognized
                gene_level_pvalue = fisher_method(padjs)
                
        except Exception as e:
            logging.warning(f"P-value combination failed for timepoint {tp}: {e}")
            gene_level_pvalue = padjs[0] if len(padjs) > 0 else 1.0
        
        gene_level_results.loc[tp, "LFC"] = observed_gene_level_lfc
        gene_level_results.loc[tp, "pvalue"] = gene_level_pvalue
    
    # Sort by timepoint and return as lists
    sorted_results = gene_level_results.sort_index()
    return sorted_results.round(3)


def calculate_gene_level_statistics(GMs_annotated: pd.DataFrame, 
                                  method: str = 'fisher') -> pd.DataFrame:
    """
    Calculate gene-level statistics for all genes.
    
    Args:
        GMs_annotated: Annotated DataFrame with insertion-level data
        method: P-value combination method
        
    Returns:
        DataFrame with gene-level LFC and p-values across timepoints
    """
    logging.info("Calculating gene-level statistics")
    
    Gene_level_statistics = pd.DataFrame()
    unique_genes = GMs_annotated["Systematic ID"].unique()
    total_genes = len(unique_genes)
    
    start_time = time.time()
    
    for n, (sysID, gene, viability, ess) in enumerate(
        GMs_annotated.groupby(["Systematic ID", "Name", "FYPOviability", "DeletionLibrary_essentiality"]).groups.keys(), 1
    ):
        gene_df = GMs_annotated.query(
            f"`Systematic ID` == '{sysID}' and Name == '{gene}' and FYPOviability == '{viability}' and DeletionLibrary_essentiality == '{ess}'"
        )
        
        try:
            lfc_and_pvalue = gene_level_lfc_and_pvalue(gene_df, method)
            
            Gene_level_statistics.loc[sysID, "Name"] = gene
            Gene_level_statistics.loc[sysID, "FYPOviability"] = viability
            Gene_level_statistics.loc[sysID, "DeletionLibrary_essentiality"] = ess
            
            # Store LFC and p-values for each timepoint
            timepoints = lfc_and_pvalue.index.tolist()
            for i, tp in enumerate(timepoints):
                if i < len(lfc_and_pvalue["LFC"]):
                    Gene_level_statistics.loc[sysID, tp] = lfc_and_pvalue.loc[tp, "LFC"]
                    # Gene_level_statistics.loc[sysID, f"{tp}_pvalue"] = lfc_and_pvalue.loc[tp, "pvalue"]
                else:
                    Gene_level_statistics.loc[sysID, tp] = np.nan
                    # Gene_level_statistics.loc[sysID, f"{tp}_pvalue"] = np.nan
            
        except Exception as e:
            logging.error(f"Error processing gene {sysID} ({gene}): {e}")
            continue
        
        # Progress reporting
        if n % 100 == 0 or n == total_genes:
            elapsed = time.time() - start_time
            rate = n / elapsed
            eta = (total_genes - n) / rate if rate > 0 else 0
            logging.info(f"Processed {n}/{total_genes} genes "
                        f"({n/total_genes*100:.1f}%) - ETA: {eta:.0f}s")
    
    # Remove genes with all NaN values
    complete_genes = Gene_level_statistics.dropna(how='all', axis=0)
    logging.info(f"Completed analysis for {len(complete_genes)} genes")
    
    return complete_genes


def generate_summary_statistics(Gene_level_statistics: pd.DataFrame) -> Dict[str, Union[int, float]]:
    """
    Generate summary statistics for gene-level analysis results.
    
    Args:
        Gene_level_statistics: DataFrame with gene-level results
        
    Returns:
        Dictionary with summary statistics
    """
    
    stats = {
        'Total genes analyzed': len(Gene_level_statistics),
        'FYPOviability: Essential genes': len(Gene_level_statistics[Gene_level_statistics['FYPOviability'] == 'inviable']),
        'FYPOviability: Non-essential genes': len(Gene_level_statistics[Gene_level_statistics['FYPOviability'] == 'viable']),
        'DeletionLibrary_essentiality: Essential genes': len(Gene_level_statistics[Gene_level_statistics['DeletionLibrary_essentiality'] == 'E']),
        'DeletionLibrary_essentiality: Non-essential genes': len(Gene_level_statistics[Gene_level_statistics['DeletionLibrary_essentiality'] == 'V'])
    }
    
    return stats


def display_summary_table(stats: Dict[str, Union[int, float]]) -> None:
    """Display summary statistics in formatted table."""
    logging.info("\n" + "="*60)
    logging.info("GENE-LEVEL DEPLETION ANALYSIS SUMMARY")
    logging.info("="*60)
    
    for key, value in stats.items():
        if isinstance(value, float):
            logging.info(f"{key:<40}: {value:.3f}")
        else:
            logging.info(f"{key:<40}: {value}")
    
    logging.info("="*60)


def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='''
        Gene-level depletion analysis for transposon insertion sequencing.
        
        This script aggregates insertion-level LFC and p-value data to generate
        gene-level fitness scores using weighted averages and statistical
        meta-analysis methods.
        ''',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('-l', '--lfc_path', type=Path, required=True,
                       help='Path to CSV file with LFC and p-value results')
    parser.add_argument('-a', '--annotations_path', type=Path, required=True,
                       help='Path to CSV file with insertion annotations')
    parser.add_argument('-o', '--output_path', type=Path, required=True,
                       help='Path for output CSV file with gene-level results')
    parser.add_argument('-m', '--method', choices=['fisher', 'stouffer'], 
                       default='fisher', help='P-value combination method')
    parser.add_argument('-c', '--cores', type=int, default=DEFAULT_CORES,
                       help='Number of CPU cores to use')
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Enable verbose logging')
    
    return parser.parse_args()


def main() -> None:
    """Main execution function."""
    start_time = time.time()
    
    # Parse arguments and setup
    args = parse_arguments()
    setup_logging(args.verbose)
    
    logging.info("Starting gene-level depletion analysis")
    logging.info(f"LFC file: {args.lfc_path}")
    logging.info(f"Annotations file: {args.annotations_path}")
    logging.info(f"P-value combination method: {args.method}")
    logging.info(f"Using {args.cores} CPU cores")
    
    try:
        # Create output directory
        args.output_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Load and process data
        results, insertion_annotations, in_gene_insertions = load_data(
            args.lfc_path, args.annotations_path
        )
        
        in_gene_LFCs, in_gene_pvalues = process_lfc_data(results, in_gene_insertions)
        
        # Merge LFC and p-value data
        GMs = pd.merge(
            in_gene_LFCs.stack().to_frame("LFC"), 
            in_gene_pvalues.stack().to_frame("Padj"), 
            how="left", left_index=True, right_index=True
        )
        
        # Calculate weights and merge with annotations
        GMs = calculate_statistical_weights(GMs)
        GMs_annotated = merge_and_normalize_weights(GMs, insertion_annotations)
        
        # Calculate gene-level statistics
        Gene_level_statistics = calculate_gene_level_statistics(GMs_annotated, args.method)
        
        # Generate and display summary statistics
        stats = generate_summary_statistics(Gene_level_statistics)
        display_summary_table(stats)
        
        # Save results
        Gene_level_statistics = Gene_level_statistics.rename_axis("Systematic ID")
        Gene_level_statistics.to_csv(args.output_path.parent/"LFC.tsv", index=True, sep="\t")
        Gene_level_statistics.to_csv(args.output_path, index=True, sep="\t")
        
        # Final summary
        elapsed_time = time.time() - start_time
        logging.info(f"\nAnalysis completed successfully in {elapsed_time:.1f} seconds")
        logging.info(f"Results saved to: {args.output_path}")
        
    except Exception as e:
        logging.error(f"Analysis failed: {e}")
        raise


if __name__ == "__main__":
    main()
